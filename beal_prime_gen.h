/*
 * =============================================================================
 * beal_prime_gen.h  —  Génération des primes de filtrage au démarrage
 * =============================================================================
 *
 * Génère les primes adaptés à une campagne (x, y, z) quelconque.
 * Remplace les constantes hardcodées dans beal_kernel2.cu.
 *
 * Usage dans beal_cpu.cu :
 *   BealPrimeSet ps = beal_gen_prime_set(cfg.x, cfg.y, cfg.z);
 *   // ps.k2 contient les primes pour beal_k2_init()
 *   // ps.l2 et ps.vram contiennent les primes pour les bitmaps
 *
 * Temps de génération : < 1s pour toute combinaison (x, y, z)
 *
 * Algorithme :
 *   On cherche des premiers P de la forme P = k * lcm(x,y,z) + 1
 *   Cela garantit que lcm(x,y,z) | (P-1), donc gcd(z, P-1) >= z,
 *   ce qui maximise le taux de rejet (~66.7% par prime pour z=3).
 *
 *   Pour chaque candidat, on vérifie :
 *     1. Miller-Rabin (16 rounds, déterministe jusqu'à 3.3×10^24)
 *     2. gcd(z, P-1) == z  (condition de rejet maximal)
 *
 *   On calcule ensuite les constantes Montgomery inline.
 * =============================================================================
 */
#pragma once
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <cassert>

// ── Arithmétique 128-bit portable ────────────────────────────────────────────

static inline uint64_t mulhi64(uint64_t a, uint64_t b) {
    return (unsigned __int128)a * b >> 64;
}
static inline uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t m) {
    return (unsigned __int128)a * b % m;
}

// ── GCD 64-bit ───────────────────────────────────────────────────────────────

static inline uint64_t gcd64(uint64_t a, uint64_t b) {
    while (b) { uint64_t t = b; b = a % b; a = t; }
    return a;
}

static inline uint64_t lcm64(uint64_t a, uint64_t b) {
    return a / gcd64(a, b) * b;
}

// ── Miller-Rabin déterministe 64-bit ─────────────────────────────────────────
//
// Déterministe jusqu'à 3.3×10^24 avec les témoins {2,3,5,7,11,13,17,19,23,29,31,37}.
// Source : https://en.wikipedia.org/wiki/Miller-Rabin_primality_test#Testing_against_small_sets
//

static uint64_t powmod64(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = mulmod64(result, base, mod);
        base = mulmod64(base, base, mod);
        exp >>= 1;
    }
    return result;
}

static bool miller_rabin_witness(uint64_t n, uint64_t a, uint64_t d, int r) {
    uint64_t x = powmod64(a, d, n);
    if (x == 1 || x == n - 1) return true;
    for (int i = 0; i < r - 1; i++) {
        x = mulmod64(x, x, n);
        if (x == n - 1) return true;
    }
    return false;
}

static bool is_prime_64(uint64_t n) {
    if (n < 2)  return false;
    if (n == 2 || n == 3 || n == 5 || n == 7) return true;
    if (n % 2 == 0 || n % 3 == 0 || n % 5 == 0) return false;

    // Écrire n-1 = d * 2^r
    uint64_t d = n - 1;
    int r = 0;
    while (d % 2 == 0) { d >>= 1; r++; }

    // Témoins déterministes pour n < 3.3×10^24
    static const uint64_t witnesses[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (uint64_t a : witnesses) {
        if (a >= n) continue;
        if (!miller_rabin_witness(n, a, d, r)) return false;
    }
    return true;
}

// ── Constantes Montgomery ─────────────────────────────────────────────────────

static inline uint64_t montgomery_ni(uint64_t p) {
    // ni = -P^{-1} mod 2^64  via itérations de Newton
    uint64_t x = 1;
    for (int i = 0; i < 6; i++) x = x * (2 - p * x);
    return (uint64_t)(-(int64_t)x);
}

static inline uint64_t montgomery_r2(uint64_t p) {
    // R² mod P  avec R = 2^64
    // IMPORTANT : ne pas caster en uint64_t avant le modulo
    unsigned __int128 R = ((unsigned __int128)1 << 64) % p;
    return (uint64_t)(R * R % p);
}

// ── Structure résultat ────────────────────────────────────────────────────────

// Réutilise BealK2Prime depuis beal_kernel2.cu
// (ce header est inclus AVANT beal_kernel2.cu dans beal_cpu.cu)
struct BealPrimeInfo {
    uint64_t p;
    uint64_t r2;
    uint64_t ni;
    uint64_t euler_exp;
    int      bits;
    double   rejection_rate;  // estimé théorique : 1 - 1/z
};

struct BealPrimeSet {
    std::vector<BealPrimeInfo> l2;    // primes ~24 bits (bitmaps L2)
    std::vector<BealPrimeInfo> vram;  // primes ~31 bits (bitmaps VRAM)
    std::vector<BealPrimeInfo> k2;    // primes ~62 bits (filtres K2)
    int x, y, z;
    uint64_t lcm_xyz;
};

// ── Générateur principal ──────────────────────────────────────────────────────

static BealPrimeInfo make_prime_info(uint64_t p, int z) {
    BealPrimeInfo pi;
    pi.p          = p;
    pi.r2         = montgomery_r2(p);
    pi.ni         = montgomery_ni(p);
    pi.euler_exp  = (p - 1) / gcd64(z, p - 1);
    pi.bits       = 0;
    uint64_t tmp  = p;
    while (tmp) { pi.bits++; tmp >>= 1; }
    pi.rejection_rate = 1.0 - 1.0 / gcd64(z, p - 1);
    return pi;
}

// Génère n premiers de la forme k*lcm+1 dans la plage [target_bits-1, target_bits+1]
static std::vector<BealPrimeInfo> gen_primes_of_size(
    int target_bits, int n, uint64_t lcm, int z)
{
    std::vector<BealPrimeInfo> result;
    uint64_t lo = (uint64_t)1 << (target_bits - 1);
    uint64_t hi = ((uint64_t)1 << target_bits) - 1;

    // Point de départ
    uint64_t k = lo / lcm;
    if (k == 0) k = 1;

    while ((int)result.size() < n) {
        uint64_t p = k * lcm + 1;

        // Vérifier qu'on reste dans la plage de bits (tolérance ±1 bit)
        if (p > (hi << 1)) break;  // trop grand, abandonner

        if (p >= lo && is_prime_64(p)) {
            uint64_t g = gcd64(z, p - 1);
            if (g == (uint64_t)z) {  // condition de rejet maximal
                result.push_back(make_prime_info(p, z));
            }
        }
        k++;
    }
    return result;
}

// ── API publique ──────────────────────────────────────────────────────────────

BealPrimeSet beal_gen_prime_set(
    int x, int y, int z,
    int n_l2   = 2,   // primes pour bitmaps L2  (~24 bits)
    int n_vram = 3,   // primes pour bitmaps VRAM (~31 bits)
    int n_k2   = 19)  // primes pour filtres K2  (~62 bits)
{
    BealPrimeSet ps;
    ps.x = x; ps.y = y; ps.z = z;
    ps.lcm_xyz = lcm64(x, lcm64(y, z));

    printf("[primes] lcm(%d,%d,%d) = %llu\n", x, y, z,
           (unsigned long long)ps.lcm_xyz);

    // Génération des primes par taille cible
    ps.l2   = gen_primes_of_size(24, n_l2,   ps.lcm_xyz, z);
    ps.vram = gen_primes_of_size(31, n_vram,  ps.lcm_xyz, z);
    ps.k2   = gen_primes_of_size(62, n_k2,   ps.lcm_xyz, z);

    // Vérifications
    if ((int)ps.l2.size()   < n_l2)   fprintf(stderr, "⚠️  Seulement %d/%d primes L2 trouvés\n",   (int)ps.l2.size(),   n_l2);
    if ((int)ps.vram.size() < n_vram)  fprintf(stderr, "⚠️  Seulement %d/%d primes VRAM trouvés\n", (int)ps.vram.size(), n_vram);
    if ((int)ps.k2.size()   < n_k2)   fprintf(stderr, "⚠️  Seulement %d/%d primes K2 trouvés\n",   (int)ps.k2.size(),   n_k2);

    printf("[primes] L2=%d  VRAM=%d  K2=%d  (taux rejet ~%.1f%% chacun)\n",
           (int)ps.l2.size(), (int)ps.vram.size(), (int)ps.k2.size(),
           (1.0 - 1.0/z) * 100.0);

    // Survie théorique totale
    int total = (int)(ps.l2.size() + ps.vram.size() + ps.k2.size());
    double surv = 1.0;
    for (int i = 0; i < total; i++) surv /= z;
    printf("[primes] Survie théorique : (1/%d)^%d = %.2e\n", z, total, surv);

    return ps;
}

// Affiche tous les primes générés (pour debug/log)
static void beal_print_prime_set(const BealPrimeSet& ps) {
    printf("\n── Primes L2 (%d) ──\n", (int)ps.l2.size());
    for (auto& p : ps.l2)
        printf("  P=%llu (%d bits) R2=%llu NI=0x%016llX\n",
               (unsigned long long)p.p, p.bits,
               (unsigned long long)p.r2, (unsigned long long)p.ni);
    printf("── Primes VRAM (%d) ──\n", (int)ps.vram.size());
    for (auto& p : ps.vram)
        printf("  P=%llu (%d bits) R2=%llu NI=0x%016llX\n",
               (unsigned long long)p.p, p.bits,
               (unsigned long long)p.r2, (unsigned long long)p.ni);
    printf("── Primes K2 (%d) ──\n", (int)ps.k2.size());
    for (auto& p : ps.k2)
        printf("  P=%llu (%d bits)\n", (unsigned long long)p.p, p.bits);
    printf("\n");
}

// Convertit BealPrimeInfo → BealK2Prime (pour beal_k2_init)
// BealK2Prime doit être défini dans beal_kernel2.cu
static std::vector<BealK2Prime> beal_prime_set_to_k2(
    const std::vector<BealPrimeInfo>& primes)
{
    std::vector<BealK2Prime> result;
    for (auto& p : primes) {
        BealK2Prime k2p;
        k2p.p          = p.p;
        k2p.r2         = p.r2;
        k2p.ni         = p.ni;
        k2p.euler_exp  = p.euler_exp;
        result.push_back(k2p);
    }
    return result;
}
