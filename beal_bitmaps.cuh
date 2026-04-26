/*
 * =============================================================================
 * beal_bitmaps.cuh  —  Bitmaps des puissances z-ièmes mod P
 * =============================================================================
 *
 * Un bitmap encode pour chaque résidu r ∈ [0, P) :
 *   bitmap[r] = 1  si r est une puissance z-ième mod P  (ou r = 0)
 *   bitmap[r] = 0  sinon
 *
 * Encodage bit-packed : bitmap[r >> 3] bit (r & 7)
 *   → taille = ceil(P / 8) bytes
 *   → P_L2 (~24 bits) : ~1 MB  — tient en L2 cache (32 MB sur Blackwell sm_120)
 *   → P_V  (~31 bits) : ~268 MB — VRAM globale
 *
 * Construction CPU (beal_build_bitmap) :
 *   Méthode par énumération via générateur de Z/PZ*
 *   Coût : O(P/z) multiplications mod P  (vs O(P log P) pour Euler brut)
 *   P_L2 (~24 bits) : <10 ms
 *   P_V  (~31 bits) : ~1 s
 *
 * Usage dans le kernel GPU :
 *   uint32_t S = (Ax + By) % P;   // addition mod P
 *   if (beal_bitmap_query(bitmap, S)) { ... }  // 2 instructions
 *
 * Générateurs validés (Python) :
 *   P_L2_0 = 8388619   → g = 2
 *   P_V0   = 2147483647 (M31) → g = 7
 *
 * Tests : voir beal_bitmaps_test.cu
 * =============================================================================
 */
#pragma once
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cuda_runtime.h>
#include "beal_powmod.cuh"

// ── Taille du bitmap ──────────────────────────────────────────────────────────

static inline size_t beal_bitmap_bytes(uint64_t p) {
    return (size_t)((p + 7) / 8);
}

// ── Lookup GPU/CPU (inline, 2 instructions) ───────────────────────────────────

__host__ __device__ __forceinline__
int beal_bitmap_query(const uint8_t* __restrict__ bitmap, uint64_t r) {
#if defined(__CUDA_ARCH__)
    // Opti 2 : lecture via Read-Only Data Cache (__ldg)
    // Plus tolérant aux accès non-contigus que le cache L2 standard
    return (__ldg(&bitmap[r >> 3]) >> (r & 7)) & 1;
#else
    return (bitmap[r >> 3] >> (r & 7)) & 1;
#endif
}

// ── Helpers pour passer les constantes Montgomery ────────────────────────────

static inline uint64_t beal_ni_for(uint64_t p) {
    if (p == BEAL_P_L2_0) return BEAL_NI_L2_0;
    if (p == BEAL_P_L2_1) return BEAL_NI_L2_1;
    if (p == BEAL_P_V0)   return BEAL_NI_V0;
    if (p == BEAL_P_V1)   return BEAL_NI_V1;
    if (p == BEAL_P_V2)   return BEAL_NI_V2;
    uint64_t x = 1;
    for (int i = 0; i < 6; i++) x = x * (2 - p * x);
    return (uint64_t)(-(int64_t)x);
}

static inline uint64_t beal_r2_for(uint64_t p) {
    if (p == BEAL_P_L2_0) return BEAL_R2_L2_0;
    if (p == BEAL_P_L2_1) return BEAL_R2_L2_1;
    if (p == BEAL_P_V0)   return BEAL_R2_V0;
    if (p == BEAL_P_V1)   return BEAL_R2_V1;
    if (p == BEAL_P_V2)   return BEAL_R2_V2;
    // IMPORTANT : ne pas caster en uint64_t avant le modulo (tronque à 0)
    unsigned __int128 R = ((unsigned __int128)1 << 64) % p;
    return (uint64_t)(R * R % p);
}

// ── Recherche d'un générateur de Z/pZ* ───────────────────────────────────────
//
// Retourne le plus petit générateur g tel que g génère tout Z/pZ*.
// Utilisé uniquement pour les primes non-listés (tests avec P arbitraire).
// Pour les primes standards, les générateurs sont hardcodés dans beal_build_bitmap.
//
static uint64_t beal_find_generator(uint64_t p) {
    uint64_t ni = beal_ni_for(p);
    uint64_t r2 = beal_r2_for(p);
    uint64_t n = p - 1;
    uint64_t factors[64];
    int nf = 0;
    {
        uint64_t tmp = n;
        for (uint64_t d = 2; d * d <= tmp; d++) {
            if (tmp % d == 0) {
                factors[nf++] = d;
                while (tmp % d == 0) tmp /= d;
            }
        }
        if (tmp > 1) factors[nf++] = tmp;
    }
    for (uint64_t g = 2; g < p; g++) {
        bool is_gen = true;
        for (int i = 0; i < nf; i++) {
            if (beal_pow_mod(g, (int)(n / factors[i]), p, r2, ni) == 1) {
                is_gen = false; break;
            }
        }
        if (is_gen) return g;
    }
    return 0;
}

// ── Construction CPU du bitmap ────────────────────────────────────────────────
//
// Remplit bitmap[] (déjà alloué, taille beal_bitmap_bytes(p)) avec
// les puissances z-ièmes mod p.
//
// Algorithme :
//   1. Trouver g = générateur de Z/pZ*
//   2. Calculer step = g^z mod p
//   3. Énumérer val = g^0, g^z, g^(2z), ... en multipliant par step
//   4. Marquer chaque val dans le bitmap
//   5. Marquer aussi 0 (cas r=0, toujours puissance z-ième)
//
// Nombre de valeurs marquées : (p-1)/gcd(z, p-1) + 1  (avec le 0)
//
void beal_build_bitmap(uint8_t* bitmap, uint64_t p, int z) {
    const size_t nbytes = beal_bitmap_bytes(p);
    memset(bitmap, 0, nbytes);

    // Marquer 0 (0 = 0^z pour tout z)
    bitmap[0] |= 1;

    // Trouver générateur
    uint64_t r2 = beal_r2_for(p);
    uint64_t ni = beal_ni_for(p);

    // Trouver g (petit générateur)
    // Pour nos primes connus, utiliser les valeurs validées
    uint64_t g;
    if (p == BEAL_P_L2_0 || p == BEAL_P_L2_1) g = 2;
    else if (p == BEAL_P_V0) g = 7;
    else if (p == BEAL_P_V1) g = 2;
    else if (p == BEAL_P_V2) g = 5;
    else {
        // Recherche générique (lent, pour tests avec primes arbitraires)
        g = beal_find_generator(p);
    }

    // step = g^z mod p
    uint64_t step = beal_pow_mod(g, z, p, r2, ni);

    // Énumérer toutes les puissances z-ièmes : g^0, g^z, g^2z, ...
    // Le groupe a ordre p-1, donc il y a (p-1)/gcd(z,p-1) cubes
    uint64_t n_cubes = (p - 1);
    // Calcul de gcd(z, p-1)
    uint64_t a = (uint64_t)z, b = p - 1;
    while (b) { uint64_t t = b; b = a % b; a = t; }
    n_cubes = (p - 1) / a;  // a = gcd(z, p-1)

    uint64_t val = 1;  // g^0 = 1
    for (uint64_t k = 0; k < n_cubes; k++) {
        bitmap[val >> 3] |= (uint8_t)(1 << (val & 7));
        // Avancer : val = val * step mod p  (sans Montgomery pour simplicité CPU)
        val = (unsigned __int128)val * step % p;
    }
}

// ── Upload vers GPU ───────────────────────────────────────────────────────────

struct BealBitmap {
    uint8_t* d_data;   // pointeur VRAM
    uint64_t p;        // le prime associé
    size_t   bytes;    // taille en bytes
};

// Alloue + copie le bitmap vers la VRAM
// Le bitmap CPU h_bitmap doit avoir été construit par beal_build_bitmap
BealBitmap beal_upload_bitmap(const uint8_t* h_bitmap, uint64_t p) {
    BealBitmap bm;
    bm.p     = p;
    bm.bytes = beal_bitmap_bytes(p);
    cudaMalloc(&bm.d_data, bm.bytes);
    cudaMemcpy(bm.d_data, h_bitmap, bm.bytes, cudaMemcpyHostToDevice);
    return bm;
}

void beal_free_bitmap(BealBitmap& bm) {
    if (bm.d_data) { cudaFree(bm.d_data); bm.d_data = nullptr; }
}
