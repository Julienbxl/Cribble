/*
 * beal_kernel1.cu  —  Tamis rapide stride + bitmaps (v3 : CUB Killer)
 *
 * Opti v3 : suppression de d_flags + CUB pour K1.
 * Les survivants sont directement appendés dans d_idx_out via warp-aggregated
 * atomics, sans écrire un tableau de flags de 100MB.
 *
 * Suppressions :
 *   - d_flags1 (100MB par batch)
 *   - BealCompactor comp1 + CUB DeviceSelect (2-3 kernels internes)
 *   Total économie : ~200MB de trafic VRAM par batch
 */
#pragma once
#include "beal_bitmaps.cuh"

#ifndef BEAL_CHUNK_SIZE
#define BEAL_CHUNK_SIZE 256
#endif
#ifndef BEAL_BLOCK_Y
#define BEAL_BLOCK_Y 128
#endif

struct BealBitmapPrime {
    uint64_t       p, r2, ni;
    const uint8_t* d_bitmap;
};

__global__ void k_beal_sieve(
    uint32_t* __restrict__ d_idx_out,
    int*      __restrict__ d_count,
    uint64_t A_start, uint64_t B_start,
    uint64_t N_A,     uint64_t N_B,
    int x, int y,
    BealBitmapPrime bp0,
    BealBitmapPrime bp1)
{
    const uint64_t chunk_idx = (uint64_t)blockIdx.x;
    const uint64_t b_idx     = (uint64_t)blockIdx.y * blockDim.x + threadIdx.x;

    const uint64_t A_base = A_start + chunk_idx * BEAL_CHUNK_SIZE;
    if (A_base >= A_start + N_A) return;

    const uint64_t chunk_end = min(
        (uint64_t)BEAL_CHUNK_SIZE,
        N_A - chunk_idx * BEAL_CHUNK_SIZE);

    // ── Shared : A^x mod P pour tout le chunk, calculé une fois par bloc ──────
    __shared__ uint64_t sAx0[BEAL_CHUNK_SIZE];
    __shared__ uint64_t sAx1[BEAL_CHUNK_SIZE];

    if (threadIdx.x == 0) {
        uint64_t d0[BEAL_MAX_EXP + 1];
        uint64_t d1[BEAL_MAX_EXP + 1];
        beal_stride_init_v2(A_base, x, bp0.p, bp0.r2, bp0.ni, d0);
        beal_stride_init_v2(A_base, x, bp1.p, bp1.r2, bp1.ni, d1);
        for (uint64_t i = 0; i < chunk_end; i++) {
            sAx0[i] = d0[0];
            sAx1[i] = d1[0];
            beal_stride_step(d0, x, bp0.p);
            beal_stride_step(d1, x, bp1.p);
        }
    }
    __syncthreads();

    // ── By calculé une fois par thread ───────────────────────────────────────
    // Les threads hors borne participent au ballot mais ne survivent jamais
    const bool b_valid = (b_idx < N_B);
    const uint64_t By0 = b_valid ? beal_pow_dispatch(B_start + b_idx, y, bp0.p, bp0.r2, bp0.ni) : 0;
    const uint64_t By1 = b_valid ? beal_pow_dispatch(B_start + b_idx, y, bp1.p, bp1.r2, bp1.ni) : 0;

    const int      lane      = (int)(threadIdx.x & 31);
    const unsigned full_mask = 0xFFFFFFFFu;

    // ── Boucle chunk : warp-aggregated append ─────────────────────────────────
    for (uint64_t i = 0; i < chunk_end; i++) {

        uint32_t flag = 0;
        if (b_valid) {
            const uint64_t S0 = beal_add_mod(sAx0[i], By0, bp0.p);
            if (beal_bitmap_query(bp0.d_bitmap, S0)) {
                const uint64_t S1 = beal_add_mod(sAx1[i], By1, bp1.p);
                if (beal_bitmap_query(bp1.d_bitmap, S1))
                    flag = 1;
            }
        }

        const unsigned alive_mask = __ballot_sync(full_mask, flag);
        if (alive_mask == 0) continue;

        const int n_alive    = __popc(alive_mask);
        const int first_lane = __ffs(alive_mask) - 1;

        // Un seul atomicAdd par warp
        int base_slot = 0;
        if (lane == first_lane)
            base_slot = atomicAdd(d_count, n_alive);
        base_slot = __shfl_sync(full_mask, base_slot, first_lane);

        if (flag) {
            const int local_rank = __popc(alive_mask & ((1u << lane) - 1));
            const uint64_t a_local = chunk_idx * BEAL_CHUNK_SIZE + i;
            d_idx_out[base_slot + local_rank] =
                (uint32_t)(a_local * N_B + b_idx);
        }
    }
}

struct BealSieveParams {
    uint64_t A_start, B_start, N_A, N_B;
    int x, y;
    BealBitmapPrime bp0;
    BealBitmapPrime bp1;
};

// Lance K1, retourne le nombre de survivants dans d_idx_out
// d_idx_out : buffer uint32_t préalloué (N_A*N_B * 15% max)
// d_count   : int* device, initialisé à 0 avant l'appel
inline int beal_sieve_launch(
    uint32_t* d_idx_out,
    int*      d_count,
    const BealSieveParams& p,
    cudaStream_t stream = 0)
{
    cudaMemsetAsync(d_count, 0, sizeof(int), stream);

    uint64_t nc = (p.N_A + BEAL_CHUNK_SIZE - 1) / BEAL_CHUNK_SIZE;
    uint64_t nb = (p.N_B + BEAL_BLOCK_Y   - 1) / BEAL_BLOCK_Y;
    dim3 grid((unsigned)nc, (unsigned)nb);
    size_t smem = 2 * BEAL_CHUNK_SIZE * sizeof(uint64_t);

    k_beal_sieve<<<grid, BEAL_BLOCK_Y, smem, stream>>>(
        d_idx_out, d_count,
        p.A_start, p.B_start, p.N_A, p.N_B,
        p.x, p.y, p.bp0, p.bp1);

    int h_count = 0;
    cudaMemcpyAsync(&h_count, d_count, sizeof(int),
                    cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
    return h_count;
}
