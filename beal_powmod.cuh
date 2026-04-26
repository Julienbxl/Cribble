/*
 * =============================================================================
 * beal_powmod.cuh  —  Arithmétique modulaire 64-bit pour le crible Beal
 * =============================================================================
 *
 * Fournit :
 *   1. Multiplication de Montgomery 64-bit (PTX sm_120 + fallback __int128)
 *   2. Exponentiation modulaire par square-and-multiply
 *   3. Initialisation des deltas pour les différences finies (stride)
 *   4. Avancement d'un step dans la boucle stride
 *
 * Le stride remplace l'exponentiation dans la boucle principale du Kernel 1 :
 *   - Setup  : (x+1) appels pow_mod_64  →  coût O(x log x) par thread
 *   - Boucle : x additions mod P par step  →  coût O(x) par candidat
 *   - Gain   : ×8 pour x=3, ×4 pour x=15 vs square-and-multiply
 *
 * Formule des deltas (différences finies d'ordre k) :
 *   Δ^k f(A) = Σ_{j=0}^{k} (-1)^{k-j} * C(k,j) * f(A+j)
 *   avec f(A) = A^x mod P
 *
 * Le dernier delta est constant : Δ^x f = x!
 *
 * Avancement d'un step (ordre CROISSANT obligatoire) :
 *   for k in 0..x-1: delta[k] = (delta[k] + delta[k+1]) % P
 *   delta[x] reste inchangé (= x!)
 *
 * Primes validés (script primes.py, x=y=z=3) :
 *   Vitesse L2  (~24 bits) : P_L2_0=8388619,   P_L2_1=8388637
 *   Vitesse VRAM (~32 bits) : P_V0=2147483647,  P_V1=2147483659,  P_V2=2147483713
 *
 * Tests : voir beal_powmod_test.cu
 * =============================================================================
 */
#pragma once
#include <cstdint>
#include <cuda_runtime.h>

// ── Constantes des primes validés ────────────────────────────────────────────

// Vitesse 1 : Cache L2 (~24 bits, bitmap = 2 MB, permanent en L2)
static constexpr uint64_t BEAL_P_L2_0  = 8388619ULL;
static constexpr uint64_t BEAL_NI_L2_0 = 0x59A0FC4C2A51745DULL;  // -P^{-1} mod 2^64
static constexpr uint64_t BEAL_R2_L2_0 = 6072010ULL;               // R² mod P, R=2^64

static constexpr uint64_t BEAL_P_L2_1  = 8388637ULL;
static constexpr uint64_t BEAL_NI_L2_1 = 0x1A45D80C2D0D3DCBULL;
static constexpr uint64_t BEAL_R2_L2_1 = 5455139ULL;

// Vitesse 2 : VRAM (~32 bits, bitmap = 268 MB)
static constexpr uint64_t BEAL_P_V0  = 2147483647ULL;   // M31 (Mersenne !)
static constexpr uint64_t BEAL_NI_V0 = 0x4000000080000001ULL;
static constexpr uint64_t BEAL_R2_V0 = 16ULL;

static constexpr uint64_t BEAL_P_V1  = 2147483659ULL;
static constexpr uint64_t BEAL_NI_V1 = 0xBCD391FBC5D1745DULL;
static constexpr uint64_t BEAL_R2_V1 = 234256ULL;

static constexpr uint64_t BEAL_P_V2  = 2147483713ULL;
static constexpr uint64_t BEAL_NI_V2 = 0xF2B71BB0BF03F03FULL;
static constexpr uint64_t BEAL_R2_V2 = 285610000ULL;

// Exposant max supporté pour le stride (dimensionne le tableau de deltas)
static constexpr int BEAL_MAX_EXP = 20;

// ── Multiplication de Montgomery 64-bit ──────────────────────────────────────
//
// Calcule a * b * R^{-1} mod P  avec R = 2^64
// Entrées : a, b en représentation Montgomery (< P)
// Sortie  : résultat en représentation Montgomery (< P)
//
#pragma diag_suppress 550
__host__ __device__ __forceinline__
uint64_t beal_mont_mul(uint64_t a, uint64_t b, uint64_t p, uint64_t ni) {
#if defined(__CUDA_ARCH__)
    uint64_t lo, hi;
    asm("mul.lo.u64 %0,%2,%3; mul.hi.u64 %1,%2,%3;"
        : "=l"(lo), "=l"(hi) : "l"(a), "l"(b));
    uint64_t m = lo * ni, thi, dummy;
    asm("mad.lo.cc.u64 %0,%2,%3,%4; madc.hi.u64 %1,%2,%3,%5;"
        : "=l"(dummy), "=l"(thi) : "l"(m), "l"(p), "l"(lo), "l"(hi));
    (void)dummy;
    return (thi >= p) ? (thi - p) : thi;
#else
    unsigned __int128 ab = (unsigned __int128)a * b;
    uint64_t m  = (uint64_t)ab * ni;
    uint64_t r  = (uint64_t)((ab + (unsigned __int128)m * p) >> 64);
    return (r >= p) ? (r - p) : r;
#endif
}

__host__ __device__ __forceinline__
uint64_t beal_add_mod(uint64_t a, uint64_t b, uint64_t p) {
    uint64_t r = a + b;
    return (r >= p) ? (r - p) : r;
}

__host__ __device__ __forceinline__
uint64_t beal_sub_mod(uint64_t a, uint64_t b, uint64_t p) {
    return (a >= b) ? (a - b) : (a + p - b);
}

// ── Conversion vers/depuis Montgomery ────────────────────────────────────────

// Convertit x natif → représentation Montgomery (x * R mod P)
__host__ __device__ __forceinline__
uint64_t beal_to_mont(uint64_t x, uint64_t p, uint64_t r2, uint64_t ni) {
    return beal_mont_mul(x, r2, p, ni);
}

// Convertit représentation Montgomery → valeur native
__host__ __device__ __forceinline__
uint64_t beal_from_mont(uint64_t x, uint64_t p, uint64_t ni) {
    return beal_mont_mul(x, 1ULL, p, ni);
}

// ── Exponentiation modulaire générique (square-and-multiply) ─────────────────
//
// Calcule base^exp mod P — utilisé quand l'exposant n'est pas connu à la compile
//
__host__ __device__ __forceinline__
uint64_t beal_pow_mod(uint64_t base, int exp, uint64_t p, uint64_t r2, uint64_t ni) {
    if (exp == 0) return 1ULL;
    if (base == 0) return 0ULL;

    uint64_t result = beal_to_mont(1ULL, p, r2, ni);
    uint64_t b      = beal_to_mont(base, p, r2, ni);

    int e = exp;
    while (e > 0) {
        if (e & 1)
            result = beal_mont_mul(result, b, p, ni);
        b = beal_mont_mul(b, b, p, ni);
        e >>= 1;
    }

    return beal_from_mont(result, p, ni);
}

// ── Exponentiation modulaire template (exposant connu à la compile) ───────────
//
// Calcule base^EXP mod P avec EXP fixe — zéro boucle, zéro branchement.
// Séquence optimale de mulmods Montgomery pour EXP = 3..11 :
//
//   EXP=3  : b² → b³                          (2 mulmods)
//   EXP=4  : b² → b⁴                          (2 mulmods)
//   EXP=5  : b² → b⁴ → b⁵                    (3 mulmods)
//   EXP=6  : b² → b³ → b⁶                    (3 mulmods)
//   EXP=7  : b² → b³ → b⁶ → b⁷              (4 mulmods)
//   EXP=8  : b² → b⁴ → b⁸                    (3 mulmods)
//   EXP=9  : b² → b⁴ → b⁸ → b⁹              (4 mulmods)
//   EXP=11 : b² → b³ → b⁶ → b¹¹             (4 mulmods + 1 mul séparé)
//
// vs square-and-multiply générique : 4 mulmods pour EXP=3, 5 pour EXP=5, etc.
//
template<int EXP>
__device__ __forceinline__
uint64_t beal_pow_const_mod(uint64_t base, uint64_t p, uint64_t r2, uint64_t ni) {
    if (base == 0) return 0ULL;

    uint64_t b = beal_to_mont(base, p, r2, ni);
    uint64_t r;

    if constexpr (EXP == 3) {
        r = beal_mont_mul(b, b, p, ni);   // b²
        r = beal_mont_mul(r, b, p, ni);   // b³
    }
    else if constexpr (EXP == 4) {
        r = beal_mont_mul(b, b, p, ni);   // b²
        r = beal_mont_mul(r, r, p, ni);   // b⁴
    }
    else if constexpr (EXP == 5) {
        r = beal_mont_mul(b, b, p, ni);   // b²
        uint64_t b4 = beal_mont_mul(r, r, p, ni);  // b⁴
        r = beal_mont_mul(b4, b, p, ni);  // b⁵
    }
    else if constexpr (EXP == 6) {
        r = beal_mont_mul(b, b, p, ni);   // b²
        r = beal_mont_mul(r, b, p, ni);   // b³
        r = beal_mont_mul(r, r, p, ni);   // b⁶
    }
    else if constexpr (EXP == 7) {
        r = beal_mont_mul(b, b, p, ni);   // b²
        r = beal_mont_mul(r, b, p, ni);   // b³
        uint64_t b6 = beal_mont_mul(r, r, p, ni);  // b⁶
        r = beal_mont_mul(b6, b, p, ni);  // b⁷
    }
    else if constexpr (EXP == 8) {
        r = beal_mont_mul(b, b, p, ni);   // b²
        r = beal_mont_mul(r, r, p, ni);   // b⁴
        r = beal_mont_mul(r, r, p, ni);   // b⁸
    }
    else if constexpr (EXP == 9) {
        r = beal_mont_mul(b, b, p, ni);   // b²
        r = beal_mont_mul(r, r, p, ni);   // b⁴
        uint64_t b8 = beal_mont_mul(r, r, p, ni);  // b⁸
        r = beal_mont_mul(b8, b, p, ni);  // b⁹
    }
    else if constexpr (EXP == 11) {
        r = beal_mont_mul(b, b, p, ni);         // b²
        uint64_t b3 = beal_mont_mul(r, b, p, ni);  // b³
        uint64_t b6 = beal_mont_mul(b3, b3, p, ni); // b⁶
        uint64_t b9 = beal_mont_mul(b6, b3, p, ni); // b⁹
        r = beal_mont_mul(b9, r, p, ni);        // b¹¹
    }
    else {
        // Fallback générique pour les exposants non listés
        return beal_pow_mod(base, EXP, p, r2, ni);
    }

    return beal_from_mont(r, p, ni);
}

// ── Dispatch runtime → template ───────────────────────────────────────────────
//
// Permet d'appeler le template avec un exposant connu au runtime mais
// fixe pendant toute la campagne. À utiliser dans les kernels via switch.
//
__device__ __forceinline__
uint64_t beal_pow_dispatch(uint64_t base, int exp, uint64_t p, uint64_t r2, uint64_t ni) {
    switch (exp) {
        case  3: return beal_pow_const_mod<3> (base, p, r2, ni);
        case  4: return beal_pow_const_mod<4> (base, p, r2, ni);
        case  5: return beal_pow_const_mod<5> (base, p, r2, ni);
        case  6: return beal_pow_const_mod<6> (base, p, r2, ni);
        case  7: return beal_pow_const_mod<7> (base, p, r2, ni);
        case  8: return beal_pow_const_mod<8> (base, p, r2, ni);
        case  9: return beal_pow_const_mod<9> (base, p, r2, ni);
        case 11: return beal_pow_const_mod<11>(base, p, r2, ni);
        default: return beal_pow_mod(base, exp, p, r2, ni);
    }
}

// ── Stride : initialisation des deltas ───────────────────────────────────────
//
// Calcule les (x+1) deltas pour la boucle stride à partir de A_base.
//
// Algorithme : application successive de l'opérateur différence Δ sur
//   le tableau f[j] = (A_base+j)^x mod P
//
// À la fin : delta[d] = Δ^d f(A_base) pour d = 0..x
//   delta[0] = f(A_base) = A_base^x mod P
//   delta[x] = x!  (constant, ne change pas pendant la boucle)
//
// Coût setup : (x+1) appels pow_mod, x*(x+1)/2 soustractions mod P
//
// Validé Python : x ∈ {3,5,7,10,15}, multiples A_base, 20 steps chacun ✅
//
__host__ __device__ __forceinline__
void beal_stride_init_v2(uint64_t A_base, int x, uint64_t p, uint64_t r2, uint64_t ni,
                          uint64_t* delta) {
    // Table f[j] = (A_base + j)^x mod P  pour j = 0..x
    uint64_t f[BEAL_MAX_EXP + 1];
    for (int j = 0; j <= x; j++)
        f[j] = beal_pow_mod(A_base + (uint64_t)j, x, p, r2, ni);

    // Application successive de l'opérateur différence Δ
    // Après d applications : f[j] = Δ^d f(A_base + j)
    // On sauvegarde f[0] = Δ^d f(A_base) à chaque étape
    for (int d = 0; d <= x; d++) {
        delta[d] = f[0];
        if (d < x) {
            // Appliquer Δ une fois de plus : f[j] = f[j+1] - f[j]
            for (int j = 0; j < x - d; j++)
                f[j] = beal_sub_mod(f[j+1], f[j], p);
        }
    }
}

// ── Stride : avancement d'un step ────────────────────────────────────────────
//
// Avance tous les deltas d'une position (A → A+1).
// Règle : delta[k] += delta[k+1]  pour k = 0, 1, ..., x-1
// L'ordre croissant est obligatoire.
// delta[x] = x! ne change pas.
//
// Après l'appel :
//   delta[0] contient f(A+1) = (A+1)^x mod P
//
__host__ __device__ __forceinline__
void beal_stride_step(uint64_t* delta, int x, uint64_t p) {
    for (int k = 0; k < x; k++)
        delta[k] = beal_add_mod(delta[k], delta[k+1], p);
}
