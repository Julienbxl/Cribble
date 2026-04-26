/*
 * beal_compaction.cuh  —  Compaction des survivants (v3)
 *
 * Opti 3 (maintenu) : itérateur fantôme CUB → pas de d_pairs_in
 * Opti new-2 : compacter des indices 32-bit au lieu de paires 64-bit
 *              → trafic de sortie K1 divisé par 2
 *              → K2 reconstruit A,B depuis l'index (div/mod)
 * Opti new-3 : K2 append direct via atomicAdd warp-aggregated
 *              → supprime flags2 + CUB pour K2 (densité finale ~0.1%)
 */
#pragma once
#include <cstdint>
#include <iterator>
#include <cuda_runtime.h>
#include <cub/cub.cuh>

// ── Encodage paires (conservé pour compatibilité K2 → output final) ───────────

__host__ __device__ __forceinline__
uint64_t beal_pack_pair(uint32_t A, uint32_t B) {
    return ((uint64_t)A << 32) | (uint64_t)B;
}

__host__ __device__ __forceinline__
void beal_unpack_pair(uint64_t packed, uint32_t& A, uint32_t& B) {
    A = (uint32_t)(packed >> 32);
    B = (uint32_t)(packed & 0xFFFFFFFFULL);
}

// ── Reconstruction d'une paire depuis un index 32-bit ─────────────────────────

__host__ __device__ __forceinline__
void beal_idx_to_pair(uint32_t idx, uint32_t N_B,
                      uint64_t A_start, uint64_t B_start,
                      uint32_t& A, uint32_t& B) {
    A = (uint32_t)(A_start + idx / N_B);
    B = (uint32_t)(B_start + idx % N_B);
}

// ── Itérateur compteur 32-bit ─────────────────────────────────────────────────

struct BealCountingIter32 {
    typedef uint32_t                        value_type;
    typedef uint32_t                        reference;
    typedef uint32_t*                       pointer;
    typedef ptrdiff_t                       difference_type;
    typedef std::random_access_iterator_tag iterator_category;

    uint32_t val;
    __host__ __device__ explicit BealCountingIter32(uint32_t v = 0) : val(v) {}
    __host__ __device__ reference operator*()  const { return val; }
    __host__ __device__ reference operator[](difference_type n) const { return (uint32_t)(val + n); }
    __host__ __device__ BealCountingIter32  operator+(difference_type n)  const { return BealCountingIter32((uint32_t)(val + n)); }
    __host__ __device__ BealCountingIter32& operator+=(difference_type n) { val += (uint32_t)n; return *this; }
    __host__ __device__ BealCountingIter32& operator++()    { ++val; return *this; }
    __host__ __device__ BealCountingIter32  operator++(int) { BealCountingIter32 t(*this); ++val; return t; }
    __host__ __device__ difference_type     operator-(const BealCountingIter32& o) const { return (difference_type)(val - o.val); }
    __host__ __device__ bool operator==(const BealCountingIter32& o) const { return val == o.val; }
    __host__ __device__ bool operator!=(const BealCountingIter32& o) const { return val != o.val; }
};

// ── Itérateur de transformation minimal ──────────────────────────────────────

template<typename OutputT, typename FunctorT, typename InputIterT>
struct BealTransformIter {
    typedef OutputT                         value_type;
    typedef OutputT                         reference;
    typedef OutputT*                        pointer;
    typedef ptrdiff_t                       difference_type;
    typedef std::random_access_iterator_tag iterator_category;

    InputIterT it;
    FunctorT   fn;

    __host__ __device__
    BealTransformIter(InputIterT it, FunctorT fn) : it(it), fn(fn) {}

    __host__ __device__ reference operator*()  const { return fn(*it); }
    __host__ __device__ reference operator[](difference_type n) const { return fn(it[n]); }
    __host__ __device__ BealTransformIter operator+(difference_type n) const {
        return BealTransformIter(it + n, fn);
    }
    __host__ __device__ BealTransformIter& operator+=(difference_type n) { it += n; return *this; }
    __host__ __device__ BealTransformIter& operator++()    { ++it; return *this; }
    __host__ __device__ BealTransformIter  operator++(int) { BealTransformIter t(*this); ++it; return t; }
    __host__ __device__ difference_type operator-(const BealTransformIter& o) const { return it - o.it; }
    __host__ __device__ bool operator==(const BealTransformIter& o) const { return it == o.it; }
    __host__ __device__ bool operator!=(const BealTransformIter& o) const { return it != o.it; }
};

// Foncteur identité pour indices 32-bit (CUB a besoin d'un foncteur)
struct BealIdentity32 {
    __host__ __device__ __forceinline__
    uint32_t operator()(const uint32_t& x) const { return x; }
};

using BealIdxIter = BealTransformIter<uint32_t, BealIdentity32, BealCountingIter32>;

// ── BealCompactor ─────────────────────────────────────────────────────────────

struct BealCompactor {
    void*  d_cub_tmp;
    size_t cub_tmp_bytes;
    int*   d_num_selected;
};

inline void beal_compact_alloc(BealCompactor& comp, uint64_t N_total) {
    comp.cub_tmp_bytes = 0;
    comp.d_cub_tmp     = nullptr;

    // Calibrer avec les vrais types (indices 32-bit)
    BealCountingIter32 c(0);
    BealIdentity32     fn;
    BealIdxIter        it(c, fn);

    cub::DeviceSelect::Flagged(
        nullptr, comp.cub_tmp_bytes,
        it, (uint8_t*)nullptr, (uint32_t*)nullptr, (int*)nullptr,
        (int)N_total);

    cudaMalloc(&comp.d_cub_tmp,      comp.cub_tmp_bytes);
    cudaMalloc(&comp.d_num_selected, sizeof(int));
}

inline void beal_compact_free(BealCompactor& comp) {
    if (comp.d_cub_tmp)      { cudaFree(comp.d_cub_tmp);      comp.d_cub_tmp      = nullptr; }
    if (comp.d_num_selected) { cudaFree(comp.d_num_selected); comp.d_num_selected = nullptr; }
    comp.cub_tmp_bytes = 0;
}

// ── Compaction K1 → indices 32-bit ───────────────────────────────────────────
//
// Sortie : d_idx_out[i] = index linéaire 32-bit du i-ème survivant
//          K2 reconstruit A,B via beal_idx_to_pair()
//
inline int beal_compact_pairs(
    BealCompactor& comp,
    const uint8_t* d_flags,
    uint64_t A_start, uint64_t B_start,   // conservés pour compatibilité API
    uint64_t N_A, uint64_t N_B,
    uint32_t* d_idx_out,                  // CHANGE : uint32_t* au lieu de uint64_t*
    cudaStream_t stream = 0)
{
    const uint64_t N_total = N_A * N_B;

    BealCountingIter32 c(0);
    BealIdentity32     fn;
    BealIdxIter        it(c, fn);

    cub::DeviceSelect::Flagged(
        comp.d_cub_tmp, comp.cub_tmp_bytes,
        it, d_flags, d_idx_out,
        comp.d_num_selected, (int)N_total, stream);

    int h_count = 0;
    cudaMemcpyAsync(&h_count, comp.d_num_selected, sizeof(int),
                    cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
    return h_count;
}

// ── K2 : append direct warp-aggregated (remplace flags2 + CUB) ───────────────
//
// Kernel appelé sur les indices survivants K1.
// Les threads qui survivent K2 font un atomicAdd warp-aggregé
// et écrivent leur paire directement dans d_pairs_out.
// Pas de flags2, pas de deuxième CUB.
//
// Note : d_count doit être initialisé à 0 avant l'appel.
//
__global__ void k_beal_k2_append(
    const uint32_t* __restrict__ d_idx,    // indices survivants K1
    int n_idx,
    uint64_t A_start, uint64_t B_start,
    uint32_t N_B,
    int x, int y,
    const void* __restrict__ d_primes_raw, // BealK2Prime*
    int n_primes,
    uint64_t* __restrict__ d_pairs_out,    // sortie finale
    int* __restrict__ d_count)             // compteur atomique
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_idx) return;

    // Reconstruire A, B depuis l'index
    const uint32_t idx = d_idx[tid];
    uint32_t A, B;
    beal_idx_to_pair(idx, N_B, A_start, B_start, A, B);

    // Filtres K2 (identique à k_beal_filter_heavy)
    struct K2P { uint64_t p, r2, ni, euler_exp; };
    const K2P* primes = (const K2P*)d_primes_raw;

    uint8_t alive = 1;
    for (int pi = 0; pi < n_primes && alive; pi++) {
        const uint64_t p  = primes[pi].p;
        const uint64_t r2 = primes[pi].r2;
        const uint64_t ni = primes[pi].ni;
        const uint64_t eu = primes[pi].euler_exp;

        const uint64_t Ax = beal_pow_dispatch((uint64_t)A, x, p, r2, ni);
        const uint64_t By = beal_pow_dispatch((uint64_t)B, y, p, r2, ni);
        const uint64_t S  = beal_add_mod(Ax, By, p);
        if (S == 0) continue;

        uint64_t result = beal_to_mont(1ULL, p, r2, ni);
        uint64_t base   = beal_to_mont(S,    p, r2, ni);
        uint64_t e      = eu;
        while (e > 0) {
            if (e & 1ULL) result = beal_mont_mul(result, base, p, ni);
            base = beal_mont_mul(base, base, p, ni);
            e >>= 1;
        }
        if (beal_from_mont(result, p, ni) != 1ULL) alive = 0;
    }

    // Warp-aggregated append
    const unsigned mask    = __ballot_sync(0xFFFFFFFF, alive);
    const int      n_alive = __popc(mask);
    const int      lane    = threadIdx.x & 31;

    int base_slot = 0;
    if (lane == 0 && n_alive > 0)
        base_slot = atomicAdd(d_count, n_alive);
    base_slot = __shfl_sync(0xFFFFFFFF, base_slot, 0);

    if (alive) {
        // Position locale dans le warp parmi les survivants
        const int local_rank = __popc(mask & ((1u << lane) - 1));
        d_pairs_out[base_slot + local_rank] = beal_pack_pair(A, B);
    }
}
