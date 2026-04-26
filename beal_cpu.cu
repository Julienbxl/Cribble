/*
 * =============================================================================
 * beal_cpu.cpp  —  Orchestrateur de campagne Beal
 * =============================================================================
 *
 * Usage :
 *   ./beal_cpu --x 3 --y 3 --z 3 --amax 50000 --bmax 50000
 *   ./beal_cpu --x 3 --y 3 --z 3 --amax 50000 --bmax 50000 --resume beal.ckpt
 *
 * Architecture :
 *   L'espace (A,B) est découpé en batches 2D de taille N_BATCH × N_BATCH.
 *   Pour chaque batch :
 *     1. k_beal_sieve (K1)         — flags[N²] sur GPU
 *     2. beal_compact_pairs         — liste dense de survivants K1
 *     3. k_beal_filter_heavy (K2)  — flags2[n_k1] sur GPU
 *     4. beal_k2_compact            — liste finale de candidats
 *     5. drain_ring                 — écriture dans le fichier candidats
 *
 * Checkpointing : fichier .ckpt sauvé après chaque batch.
 * Reprise : le curseur (A_cursor, B_cursor) reprend où on s'était arrêté.
 *
 * Sortie :
 *   - Console : progression, stats, candidats trouvés
 *   - Fichier candidats.txt : une ligne "A B" par candidat final
 *     (à traiter par beal_verify.py)
 *
 * Paramètres de batch :
 *   N_BATCH = 10 000 → VRAM ~1.3 GB, ~100ms/batch sur RTX 5060
 *   Pour vitrine A,B < 50 000 : 25 batches × ~100ms ≈ 3 secondes
 *   Pour campagne A,B < 10^6  : 10 000 batches ≈ 11 jours
 *
 * Compile :
 *   nvcc -O2 -arch=sm_120 -std=c++17 beal_cpu.cpp -o beal_cpu
 * =============================================================================
 */

#include "beal_kernel2.cu"
#include "beal_prime_gen.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <ctime>
#include <cassert>
#include <chrono>
#include <string>
#include <thread>

// ── Paramètres ────────────────────────────────────────────────────────────────

#ifndef BEAL_N_BATCH
#define BEAL_N_BATCH 10000          // taille d'un côté du batch 2D
#endif

static constexpr uint32_t CKPT_MAGIC   = 0xBEA1BEA1u;
static constexpr uint32_t CKPT_VERSION = 1u;

// ── Structures ────────────────────────────────────────────────────────────────

struct BealConfig {
    int      x, y, z;
    uint64_t A_min, B_min;  // borne inférieure de l'espace de recherche
    uint64_t A_max, B_max;
    uint64_t N_batch;
    const char* ckpt_path;
    const char* out_path;
};

struct BealCheckpoint {
    uint32_t magic;
    uint32_t version;
    int      x, y, z;
    uint64_t A_min, B_min;
    uint64_t A_max, B_max;
    uint64_t A_cursor;
    uint64_t B_cursor;
    uint64_t candidates_found;
    uint64_t batches_done;
    uint64_t elapsed_ms;
};

// ── Checkpoint I/O ────────────────────────────────────────────────────────────

static bool ckpt_save(const BealCheckpoint& ck, const char* path) {
    FILE* f = fopen(path, "wb");
    if (!f) { fprintf(stderr, "Cannot write checkpoint %s\n", path); return false; }
    fwrite(&ck, sizeof(ck), 1, f);
    fclose(f);
    return true;
}

static bool ckpt_load(BealCheckpoint& ck, const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) return false;
    size_t n = fread(&ck, sizeof(ck), 1, f);
    fclose(f);
    if (n != 1 || ck.magic != CKPT_MAGIC || ck.version != CKPT_VERSION) {
        fprintf(stderr, "Checkpoint invalide ou corrompu : %s\n", path);
        return false;
    }
    return true;
}

// ── Runner ────────────────────────────────────────────────────────────────────

struct BealRunner {
    BealConfig    cfg;
    BealK2State   k2state;
    // comp1 supprimé — K1 fait maintenant un append direct (CUB Killer)

    struct BitmapEntry { uint8_t* h_data; BealBitmap gpu; };
    std::vector<BitmapEntry> bitmaps;
    BealBitmapPrime          bp_l2;
    BealBitmapPrime          bp_vram;

    uint64_t N_batch_sq;
    // d_flags1 supprimé — plus de tableau de flags 100MB
    int*      d_k1_count; // compteur atomique K1 (device)
    uint32_t* d_pairs1;   // indices 32-bit survivants K1
    uint64_t* d_pairs2;   // paires 64-bit finales (sortie K2)
    FILE*     out_file;
    uint64_t  candidates_found, batches_done;

    void init() {
        N_batch_sq = cfg.N_batch * cfg.N_batch;

        BealPrimeSet ps = beal_gen_prime_set(cfg.x, cfg.y, cfg.z,
                                              1, 1, 19);

        {
            auto& pi = ps.l2[0];
            size_t nb = beal_bitmap_bytes(pi.p);
            BitmapEntry e;
            e.h_data = (uint8_t*)malloc(nb);
            printf("[init] Bitmap L2   P=%llu (%zu MB)...\n",
                   (unsigned long long)pi.p, nb >> 20);
            beal_build_bitmap(e.h_data, pi.p, cfg.z);
            e.gpu = beal_upload_bitmap(e.h_data, pi.p);
            bp_l2 = BealBitmapPrime{ pi.p, pi.r2, pi.ni, e.gpu.d_data };
            bitmaps.push_back(std::move(e));
        }
        {
            auto& pi = ps.vram[0];
            size_t nb = beal_bitmap_bytes(pi.p);
            BitmapEntry e;
            e.h_data = (uint8_t*)malloc(nb);
            printf("[init] Bitmap VRAM P=%llu (%zu MB)...\n",
                   (unsigned long long)pi.p, nb >> 20);
            beal_build_bitmap(e.h_data, pi.p, cfg.z);
            e.gpu = beal_upload_bitmap(e.h_data, pi.p);
            bp_vram = BealBitmapPrime{ pi.p, pi.r2, pi.ni, e.gpu.d_data };
            bitmaps.push_back(std::move(e));
        }
        printf("[init] 2 bitmaps prêts ✅\n");

        auto k2v = beal_prime_set_to_k2(ps.k2);
        printf("[init] K2: %d primes, P[0]=%llu euler[0]=%llu\n",
               (int)k2v.size(),
               k2v.empty() ? 0ULL : (unsigned long long)k2v[0].p,
               k2v.empty() ? 0ULL : (unsigned long long)k2v[0].euler_exp);
        beal_k2_init(k2state, k2v.data(), (int)k2v.size());

        // Buffers GPU — d_flags1 supprimé (économie ~100MB)
        cudaMalloc(&d_k1_count, sizeof(int));
        cudaMalloc(&d_pairs1,   N_batch_sq * 4);  // uint32_t indices K1
        cudaMalloc(&d_pairs2,   N_batch_sq * 8);  // uint64_t paires K2

        out_file = fopen(cfg.out_path, "a");
        if (!out_file) { fprintf(stderr,"Cannot open %s\n",cfg.out_path); exit(1); }
        candidates_found = batches_done = 0;
        printf("[init] Buffers GPU alloués (%.0f MB)  [d_flags supprimé ✅]\n",
               (double)(N_batch_sq * 12) / 1e6);  // idx(4) + pairs(8)
    }

    void free_all() {
        for (auto& e : bitmaps) { beal_free_bitmap(e.gpu); ::free(e.h_data); }
        beal_k2_free(k2state);
        cudaFree(d_k1_count);
        cudaFree(d_pairs1);
        cudaFree(d_pairs2);
        if (out_file) { fclose(out_file); out_file = nullptr; }
    }

    int run_batch(uint64_t A_start, uint64_t N_A, uint64_t B_start, uint64_t N_B) {
        BealSieveParams sp;
        sp.A_start = A_start; sp.B_start = B_start;
        sp.N_A     = N_A;     sp.N_B     = N_B;
        sp.x       = cfg.x;  sp.y        = cfg.y;
        sp.bp0     = bp_l2;
        sp.bp1     = bp_vram;

        // K1 : append direct — plus de d_flags, plus de CUB
        int cnt1 = beal_sieve_launch(d_pairs1, d_k1_count, sp);
        if (cnt1 == 0) return 0;

        // K2 : append direct
        int cnt2 = beal_k2_launch_append(
            d_pairs1, cnt1,
            A_start, B_start, (uint32_t)N_B,
            cfg.x, cfg.y, k2state,
            d_pairs2);
        if (cnt2 == 0) return 0;

        std::vector<uint64_t> h_pairs(cnt2);
        cudaMemcpy(h_pairs.data(), d_pairs2, cnt2 * 8, cudaMemcpyDeviceToHost);

        for (int i = 0; i < cnt2; i++) {
            uint32_t A, B;
            beal_unpack_pair(h_pairs[i], A, B);
            fprintf(out_file, "%u %u\n", A, B);
            printf("  🎯 Candidat : A=%u B=%u\n", A, B);
        }
        fflush(out_file);
        return cnt2;
    }

    void run(uint64_t A_cursor_start = 1, uint64_t B_cursor_start = 1,
             uint64_t candidates_init = 0, uint64_t batches_init = 0) {
        candidates_found = candidates_init;
        batches_done     = batches_init;

        // Calcul exact du nombre de batches restants depuis le curseur
        uint64_t total_batches = 0;
        for (uint64_t b = B_cursor_start; b <= cfg.B_max; b += cfg.N_batch) {
            uint64_t a_start = (b == B_cursor_start) ? A_cursor_start : cfg.A_min;
            total_batches += (cfg.A_max - a_start) / cfg.N_batch + 1;
        }

        printf("[run] Espace : A ∈ [%llu,%llu], B ∈ [%llu,%llu], x=%d y=%d z=%d\n",
               (unsigned long long)cfg.A_min, (unsigned long long)cfg.A_max,
               (unsigned long long)cfg.B_min, (unsigned long long)cfg.B_max,
               cfg.x, cfg.y, cfg.z);
        printf("[run] Batch : %llu×%llu, ~%llu batches total\n",
               (unsigned long long)cfg.N_batch,
               (unsigned long long)cfg.N_batch,
               (unsigned long long)total_batches);
        printf("[run] Reprise depuis A=%llu B=%llu\n",
               (unsigned long long)A_cursor_start,
               (unsigned long long)B_cursor_start);

        auto t_start = std::chrono::steady_clock::now();

        for (uint64_t B_start = B_cursor_start;
             B_start <= cfg.B_max;
             B_start += cfg.N_batch)
        {
            uint64_t N_B = std::min(cfg.N_batch, cfg.B_max - B_start + 1);

            uint64_t A_start_loop = (B_start == B_cursor_start) ? A_cursor_start : cfg.A_min;

            for (uint64_t A_start = A_start_loop;
                 A_start <= cfg.A_max;
                 A_start += cfg.N_batch)
            {
                uint64_t N_A = std::min(cfg.N_batch, cfg.A_max - A_start + 1);

                auto t0 = std::chrono::steady_clock::now();
                int found = run_batch(A_start, N_A, B_start, N_B);
                auto t1 = std::chrono::steady_clock::now();
                double ms = std::chrono::duration<double,std::milli>(t1-t0).count();

                candidates_found += found;
                batches_done++;

                // Progression
                if (batches_done % 10 == 0 || found > 0) {
                    double pct = 100.0 * batches_done / total_batches; (void)pct;
                    auto elapsed = std::chrono::duration<double>(t1-t_start).count();
                    // batches_done compte depuis le début de ce run (pas depuis 0)
                    uint64_t done_this_run = batches_done - batches_init;
                    double eta = (done_this_run > 0)
                        ? (elapsed / done_this_run) * (total_batches - done_this_run)
                        : 0.0;
                    double pct_run = (total_batches > 0)
                        ? 100.0 * done_this_run / total_batches
                        : 100.0;
                    printf("[%6.2f%%] A=[%llu,%llu] B=[%llu,%llu] "
                           "| %.0fms | cand=%llu | ETA=%.0fs\n",
                           pct_run,
                           (unsigned long long)A_start,
                           (unsigned long long)(A_start+N_A-1),
                           (unsigned long long)B_start,
                           (unsigned long long)(B_start+N_B-1),
                           ms,
                           (unsigned long long)candidates_found,
                           eta);
                }

                // Checkpoint après chaque batch
                BealCheckpoint ck;
                ck.magic   = CKPT_MAGIC;
                ck.version = CKPT_VERSION;
                ck.x = cfg.x; ck.y = cfg.y; ck.z = cfg.z;
                ck.A_min = cfg.A_min; ck.B_min = cfg.B_min;
                ck.A_max = cfg.A_max; ck.B_max = cfg.B_max;
                // Le prochain batch commence à A_start + N_batch
                uint64_t next_A = A_start + cfg.N_batch;
                uint64_t next_B = B_start;
                if (next_A > cfg.A_max) { next_A = 1; next_B = B_start + cfg.N_batch; }
                ck.A_cursor = next_A;
                ck.B_cursor = next_B;
                ck.candidates_found = candidates_found;
                ck.batches_done     = batches_done;
                ck.elapsed_ms = (uint64_t)std::chrono::duration<double,std::milli>(
                    t1 - t_start).count();
                ckpt_save(ck, cfg.ckpt_path);
            }
        }

        auto t_end = std::chrono::steady_clock::now();
        double total_s = std::chrono::duration<double>(t_end - t_start).count();
        printf("\n[done] Campagne terminée en %.1fs\n", total_s);
        printf("[done] %llu batches, %llu candidats finaux\n",
               (unsigned long long)batches_done,
               (unsigned long long)candidates_found);
        printf("[done] Candidats écrits dans : %s\n", cfg.out_path);
    }
};

// ── main ──────────────────────────────────────────────────────────────────────

static void print_usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s --x X --y Y --z Z --amax N --bmax N\n"
        "          [--amin N] [--bmin N] [--batch B] [--resume FILE] [--out FILE]\n"
        "\n"
        "Options:\n"
        "  --x X       Exposant de A (défaut: 3)\n"
        "  --y Y       Exposant de B (défaut: 3)\n"
        "  --z Z       Exposant de C (défaut: 3)\n"
        "  --amin N    Valeur minimale de A (défaut: 1)\n"
        "  --bmin N    Valeur minimale de B (défaut: 1)\n"
        "  --amax N    Valeur maximale de A (défaut: 10000)\n"
        "  --bmax N    Valeur maximale de B (défaut: 10000)\n"
        "  --batch B   Taille de batch (défaut: %d)\n"
        "  --resume F  Reprendre depuis le checkpoint F\n"
        "  --out F     Fichier de sortie des candidats (défaut: candidates.txt)\n",
        prog, BEAL_N_BATCH);
}

int main(int argc, char** argv) {
    BealConfig cfg;
    cfg.x = 3; cfg.y = 3; cfg.z = 3;
    cfg.A_min = 1; cfg.B_min = 1;
    cfg.A_max = 10000; cfg.B_max = 10000;
    cfg.N_batch   = BEAL_N_BATCH;
    cfg.ckpt_path = "beal.ckpt";
    cfg.out_path  = "candidates.txt";
    bool resume = false;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"--x")      && i+1<argc) cfg.x      = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--y")      && i+1<argc) cfg.y      = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--z")      && i+1<argc) cfg.z      = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--amin")   && i+1<argc) cfg.A_min  = atoll(argv[++i]);
        else if (!strcmp(argv[i],"--bmin")   && i+1<argc) cfg.B_min  = atoll(argv[++i]);
        else if (!strcmp(argv[i],"--amax")   && i+1<argc) cfg.A_max  = atoll(argv[++i]);
        else if (!strcmp(argv[i],"--bmax")   && i+1<argc) cfg.B_max  = atoll(argv[++i]);
        else if (!strcmp(argv[i],"--batch")  && i+1<argc) cfg.N_batch= atoll(argv[++i]);
        else if (!strcmp(argv[i],"--resume") && i+1<argc) { cfg.ckpt_path = argv[++i]; resume = true; }
        else if (!strcmp(argv[i],"--out")    && i+1<argc) cfg.out_path = argv[++i];
        else if (!strcmp(argv[i],"--help"))  { print_usage(argv[0]); return 0; }
        else { fprintf(stderr, "Argument inconnu: %s\n", argv[i]); print_usage(argv[0]); return 1; }
    }

    // Validation
    if (cfg.x < 3 || cfg.y < 3 || cfg.z < 3) {
        fprintf(stderr, "x, y et z doivent être >= 3 (conjecture de Beal)\n");
        return 1;
    }
    if (cfg.x > BEAL_MAX_EXP || cfg.y > BEAL_MAX_EXP) {
        fprintf(stderr, "x et y doivent être <= %d (limite du stride GPU)\n", BEAL_MAX_EXP);
        return 1;
    }
    if (cfg.N_batch > 20000) {
        fprintf(stderr, "N_batch > 20000 peut dépasser 4GB VRAM. Continuer ? (y/n) ");
        char c = getchar();
        if (c != 'y') return 1;
    }

    // Le curseur de départ = A_min/B_min (sauf si checkpoint reprend)
    uint64_t A_cursor = cfg.A_min, B_cursor = cfg.B_min;
    uint64_t candidates_init = 0, batches_init = 0;
    if (resume) {
        BealCheckpoint ck;
        if (ckpt_load(ck, cfg.ckpt_path)) {
            if (ck.x != cfg.x || ck.y != cfg.y || ck.z != cfg.z) {
                fprintf(stderr, "Checkpoint incompatible : exposants différents\n");
                return 1;
            }
            cfg.A_min = ck.A_min; cfg.B_min = ck.B_min;
            cfg.A_max = ck.A_max; cfg.B_max = ck.B_max;
            A_cursor        = ck.A_cursor;
            B_cursor        = ck.B_cursor;
            candidates_init = ck.candidates_found;
            batches_init    = ck.batches_done;
            printf("[resume] Reprise : A=%llu B=%llu, %llu candidats déjà trouvés\n",
                   (unsigned long long)A_cursor,
                   (unsigned long long)B_cursor,
                   (unsigned long long)candidates_init);
        } else {
            fprintf(stderr, "Checkpoint non trouvé ou invalide, démarrage depuis le début\n");
        }
    }

    BealRunner runner;
    runner.cfg = cfg;
    runner.init();
    runner.run(A_cursor, B_cursor, candidates_init, batches_init);
    runner.free_all();
    return 0;
}
