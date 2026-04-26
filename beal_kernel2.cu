/*
 * =============================================================================
 * beal_kernel2.cu  —  Tamis lourd : filtres 64-bit sur survivants compactés
 * =============================================================================
 *
 * Entrée  : tableau de paires (A,B) compactées par beal_compaction.cuh
 *           (~11% de l'espace initial après kernel1)
 * Sortie  : flags[i] = 1 si la paire survit à tous les primes 64-bit
 *           (~1.2% des entrées pour 4 primes, z=3)
 *
 * Différence vs Kernel 1 :
 *   K1 : stride sur A contigu + bitmap lookup (accès mémoire régulier)
 *   K2 : paires éparpillées + pow_mod direct + test Euler (pur ALU)
 *        Pas de bitmap (les primes 64-bit nécessiteraient 2^58 bytes !)
 *        Pas de stride (les A ne sont plus contigus après compaction)
 *
 * Test Euler :
 *   S est une puissance z-ième mod P  ⟺  S^euler_exp ≡ 1 (mod P)
 *   euler_exp = (P-1) / gcd(z, P-1)   précalculé sur CPU
 *
 * Primes 64-bit validés (z=3, gcd(3,P-1)=3) :
 *   P0 = 2305843009213693951  (M61 = 2^61-1, Mersenne !)
 *   P1 = 2305843009213694017
 *   P2 = 2305843009213694149
 *   P3 = 2305843009213694173
 *
 * Tests : voir beal_kernel2_test.cu
 * =============================================================================
 */
#pragma once
#include "beal_bitmaps.cuh"
#include "beal_kernel1.cu"
#include "beal_compaction.cuh"

// ── Primes 64-bit et constantes ───────────────────────────────────────────────

#define BEAL_K2_MAX_PRIMES 19

struct BealK2Prime {
    uint64_t p;           // le prime
    uint64_t r2;          // R² mod P (pour Montgomery)
    uint64_t ni;          // -P^{-1} mod 2^64 (pour Montgomery)
    uint64_t euler_exp;   // (P-1) / gcd(z, P-1)
};

// 19 primes 64-bit pour z=3, gcd(3,P-1)=3 pour chacun.
// Taux de survie : (1/3)^19 ≈ 8.6e-10 → ~860 candidats sur 10^12 paires.
// Générés et vérifiés par primes.py + Python (constantes Montgomery validées).
static const BealK2Prime BEAL_K2_PRIMES_Z3[19] = {
    // 4 primes originaux
    { 2305843009213693951ULL, 64ULL,
      0x2000000000000001ULL, 768614336404564650ULL },
    { 2305843009213694017ULL, 270400ULL,
      0x103F03F03F03F03FULL, 768614336404564672ULL },
    { 2305843009213694149ULL, 2483776ULL,
      0x56942462C2EC81F3ULL, 768614336404564716ULL },
    { 2305843009213694173ULL, 3125824ULL,
      0x488B01288B01288BULL, 768614336404564724ULL },
    // 15 primes supplémentaires
    { 2305843009213694257ULL, 5953600ULL,
      0x100D6DF43FCA482FULL, 768614336404564752ULL },
    { 2305843009213694317ULL, 8526400ULL,
      0xD38CF9B00B38CF9BULL, 768614336404564772ULL },
    { 2305843009213694323ULL, 8809024ULL,
      0x3FBDC1E5C76AF445ULL, 768614336404564774ULL },
    { 2305843009213694443ULL, 15429184ULL,
      0x1166B67C16F0E13DULL, 768614336404564814ULL },
    { 2305843009213694491ULL, 18593344ULL,
      0x1499E28D8942F7EDULL, 768614336404564830ULL },
    { 2305843009213694497ULL, 19009600ULL,
      0x983FE1F00783FE1FULL, 768614336404564832ULL },
    { 2305843009213694569ULL, 24364096ULL,
      0xD90CDCBB8A2A5227ULL, 768614336404564856ULL },
    { 2305843009213694683ULL, 34199104ULL,
      0xCD0059A70C41D6ADULL, 768614336404564894ULL },
    { 2305843009213694791ULL, 45050944ULL,
      0x270493AE39F94989ULL, 768614336404564930ULL },
    { 2305843009213694851ULL, 51724864ULL,
      0xFD42596A373E5CD5ULL, 768614336404564950ULL },
    { 2305843009213694887ULL, 55950400ULL,
      0x378BE90046178BE9ULL, 768614336404564962ULL },
    { 2305843009213694917ULL, 59598400ULL,
      0x6FD9CC88E3157CF3ULL, 768614336404564972ULL },
    { 2305843009213695001ULL, 70425664ULL,
      0xABD6065857DAE7D7ULL, 768614336404565000ULL },
    { 2305843009213695187ULL, 97614400ULL,
      0x5510CA5003510CA5ULL, 768614336404565062ULL },
    { 2305843009213695349ULL, 124902976ULL,
      0x2BBA5D7518BD1D23ULL, 768614336404565116ULL },
};

// ── Kernel principal ──────────────────────────────────────────────────────────
//
// 1 thread = 1 paire (A,B) survivante
// Pour chaque prime : calculer S = A^x + B^y mod P, tester S^euler_exp == 1
// Si un prime échoue → flag = 0, continuer (pas de return pour éviter divergence)
//
__global__ void k_beal_filter_heavy(
    uint8_t*        __restrict__ flags,       // output [n_pairs], 0 = éliminé
    const uint64_t* __restrict__ pairs,       // input : paires compactées
    int n_pairs,
    int x, int y,
    const BealK2Prime* __restrict__ primes,   // tableau de primes en device
    int n_primes)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_pairs) return;

    uint32_t A, B;
    beal_unpack_pair(pairs[idx], A, B);

    uint8_t alive = 1;

    for (int pi = 0; pi < n_primes && alive; pi++) {
        const uint64_t p  = primes[pi].p;
        const uint64_t r2 = primes[pi].r2;
        const uint64_t ni = primes[pi].ni;
        const uint64_t eu = primes[pi].euler_exp;

        // S = A^x + B^y mod P
        const uint64_t Ax = beal_pow_mod((uint64_t)A, x, p, r2, ni);
        const uint64_t By = beal_pow_mod((uint64_t)B, y, p, r2, ni);
        const uint64_t S  = beal_add_mod(Ax, By, p);

        // S == 0 → toujours puissance z-ième (0 = 0^z)
        if (S == 0) continue;

        // Test Euler : S^euler_exp mod P == 1 ?
        // euler_exp peut dépasser int → on passe par beal_pow_mod avec cast
        // Note : euler_exp est un uint64_t, mais beal_pow_mod prend int pour exp
        // On fait l'exponentiation manuelle pour les grands exposants
        uint64_t result = beal_to_mont(1ULL, p, r2, ni);
        uint64_t base   = beal_to_mont(S,    p, r2, ni);
        uint64_t e      = eu;
        while (e > 0) {
            if (e & 1ULL) result = beal_mont_mul(result, base, p, ni);
            base = beal_mont_mul(base, base, p, ni);
            e >>= 1;
        }
        const uint64_t S_pow = beal_from_mont(result, p, ni);

        if (S_pow != 1ULL) alive = 0;
    }

    flags[idx] = alive;
}

// ── Fonctions de lancement ────────────────────────────────────────────────────

struct BealK2State {
    BealK2Prime* d_primes;
    int n_primes;
    int* d_count;   // compteur atomique pour l'append direct
};

inline void beal_k2_init(BealK2State& state, const BealK2Prime* h_primes, int n) {
    state.n_primes = n;
    cudaMalloc(&state.d_primes, n * sizeof(BealK2Prime));
    cudaMemcpy(state.d_primes, h_primes, n * sizeof(BealK2Prime),
               cudaMemcpyHostToDevice);
    cudaMalloc(&state.d_count, sizeof(int));
}

inline void beal_k2_free(BealK2State& state) {
    if (state.d_primes) { cudaFree(state.d_primes); state.d_primes = nullptr; }
    if (state.d_count)  { cudaFree(state.d_count);  state.d_count  = nullptr; }
}

// Lance le filtre lourd sur la liste compactée
// flags doit être alloué et zéré (taille n_pairs)
inline void beal_k2_launch(
    uint8_t* d_flags,
    const uint64_t* d_pairs,
    int n_pairs,
    int x, int y,
    const BealK2State& state,
    cudaStream_t stream = 0)
{
    if (n_pairs == 0) return;
    const int threads = 256;
    const int blocks  = (n_pairs + threads - 1) / threads;
    k_beal_filter_heavy<<<blocks, threads, 0, stream>>>(
        d_flags, d_pairs, n_pairs, x, y,
        state.d_primes, state.n_primes);
}

// Note : beal_k2_compact supprimé — remplacé par beal_k2_launch_append()
// défini dans beal_compaction.cuh (append direct warp-aggregated).

// ── Wrapper de lancement K2 append direct ────────────────────────────────────
inline int beal_k2_launch_append(
    const uint32_t* d_idx,
    int n_idx,
    uint64_t A_start, uint64_t B_start,
    uint32_t N_B,
    int x, int y,
    const BealK2State& state,
    uint64_t* d_pairs_out,
    cudaStream_t stream = 0)
{
    if (n_idx == 0) return 0;

    cudaMemsetAsync(state.d_count, 0, sizeof(int), stream);

    const int threads = 256;
    const int blocks  = (n_idx + threads - 1) / threads;

    k_beal_k2_append<<<blocks, threads, 0, stream>>>(
        d_idx, n_idx,
        A_start, B_start, N_B,
        x, y,
        state.d_primes, state.n_primes,
        d_pairs_out, state.d_count);

    int h_count = 0;
    cudaMemcpyAsync(&h_count, state.d_count, sizeof(int),
                    cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
    return h_count;
}
