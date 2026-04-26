#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════╗
║   BEAL COMMANDER v3 — Taxonomie Canonique (Primitive / Lift / Tags)║
╚══════════════════════════════════════════════════════════════════════╝

Philosophie de la v3
--------------------
La v2 mélangeait trois niveaux logiques :
  1) generating families (primitive forms),
  2) lift modes / homotheties,
  3) transverse arithmetic properties (p-adiques, symétries, etc.).

La v3 sépare explicitement ces couches :
  • primitive_family : mécanisme générateur fondamental
  • lift_type        : comment la primitive est relevée / homogénéisée
  • tags             : propriétés arithmétiques descriptives

Le but n'est pas seulement de “tout classer”, mais de classer proprement.
"""

import sys, os, math, time, subprocess, csv
from datetime import datetime
from collections import defaultdict, Counter

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BEAL_BIN   = "./beal_cpu"

X_RANGE    = [3, 4, 5, 6, 7, 8, 9, 10]
Y_RANGE    = [3, 4, 5, 6, 7, 8, 9, 10]
Z_RANGE    = [3, 4, 5, 6, 7, 8, 9, 10]

SCAN_AMAX  = 1_000_000
SCAN_BMAX  = 1_000_000

CANDIDATES_DIR = "candidates"
BIBLE_FILE     = "bible_v4.csv"
LOG_FILE       = "beal_commander_v4.log"

PADIC_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23]




LIFT_LABELS = {
    "LT0": "none",
    "LT1": "gcd-common",
    "LT2": "homogeneous-g",
    "LT3": "power-normalized",
    "LT4": "asym-residual",
}

MACRO_LABELS = {
    "M1": "symmetric",
    "M2": "residual",
    "M3": "binomial",
    "F0": "sporadic",
}

SUB_LABELS = {
    "M1.1": "sym-pure",
    "M1.2": "next-power",
    "M1.3": "scaled-sum",
    "M1.4": "weighted-lift",
    "M2.1": "direct-residual",
    "M2.2": "mixed-residual",
    "M3.G0": "binomial-G0",
    "M3.Gv": "binomial-Gv",
    "M3.G1": "binomial-G1",
    "M3.R":  "binomial-off-corridor",
    "M3.S":  "binomial-sym",
    "F0":    "sporadic",
}

SUB_REGIMES = {
    "M1.1": "prime", "M1.2": "prime", "M1.3": "mixed",
    "M1.4": "even",  "M2.1": "mixed", "M2.2": "mixed",
    "M3.G0":  "even",  "M3.G0d": "prime",
    "M3.Gv":  "even",  "M3.Gvd": "prime",
    "M3.G1":  "mixed", "M3.R":   "mixed", "M3.S": "even",
    "F0":   "both",
}

# ==============================================================================
# UTILITAIRES
# ==============================================================================

def log(msg):
    ts = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    line = f"{ts} {msg}"
    print(line)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")


def iroot(n: int, k: int):
    """Racine k-ième entière exacte via Newton."""
    if n < 0 or k < 1:
        return None
    if n == 0:
        return 0
    if n == 1:
        return 1
    if k == 1:
        return n
    try:
        x = int(round(n ** (1.0 / k)))
    except OverflowError:
        x = 2
    x = max(x, 2)
    while True:
        xk1 = x ** (k - 1)
        xnew = ((k - 1) * x + n // xk1) // k
        if xnew >= x:
            break
        x = xnew
    for c in [x - 1, x, x + 1]:
        if c > 0 and c ** k == n:
            return c
    return None


def is_perfect_power(n: int, z_max: int = 40):
    for z in range(z_max, 2, -1):
        c = iroot(n, z)
        if c is not None:
            return c, z
    return None


def normalize_base(val, exp):
    """Réduit val à sa base primitive, en absorbant les puissances parfaites."""
    res = is_perfect_power(val)
    if res is not None:
        base, k = res
        return base, exp * k
    return val, exp


def val_p(n: int, p: int) -> int:
    v = 0
    while n > 0 and n % p == 0:
        n //= p
        v += 1
    return v


# ==============================================================================
# CLASSIFICATION v4 — Macro-familles M1/M2/M3/F0
# ==============================================================================

"""
Classification v4 — Macro-families M1/M2/M3/F0
"""
import math

def iroot(n, k):
    if n < 0 or k < 1: return None
    if n == 0: return 0
    if k == 1: return n
    try: x = int(round(n ** (1.0/k)))
    except OverflowError: x = 2
    x = max(x, 2)
    while True:
        xk1 = x**(k-1)
        xnew = ((k-1)*x + n//xk1)//k
        if xnew >= x: break
        x = xnew
    for c in [x-1, x, x+1]:
        if c > 0 and c**k == n: return c
    return None

def normalize_base(val, exp):
    for z in range(40, 2, -1):
        c = iroot(val, z)
        if c is not None:
            return c, exp*z
    return val, exp

def val_p(n, p):
    v = 0
    while n > 0 and n % p == 0:
        n //= p; v += 1
    return v

PADIC_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23]


class Solution:
    """Précalculs centralisés pour la classification."""
    def __init__(self, A, B, x, y, C, z):
        self.A, self.B, self.x, self.y = A, B, x, y
        self.C, self.z = C, z

        # gcd et réduction
        self.g = math.gcd(math.gcd(A, B), C)
        self.a = A // self.g
        self.b = B // self.g
        self.c = C // self.g

        # Residual signature (α,β,γ)
        self.m = min(x, y, z)
        self.alpha = x - self.m
        self.beta  = y - self.m
        self.gamma = z - self.m

        # Bases normalisées (détecte A=p^k déguisés)
        self.An, self.xn = normalize_base(A, x)
        self.Bn, self.yn = normalize_base(B, y)
        self.Cn, self.zn = normalize_base(C, z)

        # Canonicalisation z (ex: 2048^5=32^11 → z_canonical=11)
        res = None
        for p in range(2, 40):
            c2 = iroot(C, p)
            if c2 is not None and c2**p == C and c2 < C:
                # garder le plus grand exposant
                if res is None or p > res[1]:
                    res = (c2, p)
        if res and res[1] > z:
            self.C_canonical = res[0]
            self.z_canonical = res[1]
        else:
            self.C_canonical = C
            self.z_canonical = z

    def ratio(self):
        """Retourne (r, small, large) si A|B ou B|A, sinon None."""
        if self.A > 0 and self.B > 0:
            if self.A % self.B == 0:
                return self.A // self.B, self.B, self.A, self.y, self.x
            if self.B % self.A == 0:
                return self.B // self.A, self.A, self.B, self.x, self.y
        return None

    def residual_eq(self):
        """Équation résiduelle : a^x·g^α + b^y·g^β = c^z·g^γ"""
        try:
            lhs = self.a**self.x * self.g**self.alpha + \
                  self.b**self.y * self.g**self.beta
            rhs = self.c**self.z * self.g**self.gamma
            return lhs == rhs
        except OverflowError:
            return False


def _geometry(s: Solution):
    """Detect the geometric corridor.
    G0d/Gvd = diagonal case r=1 (A==B on ratio side).
    G0/Gv   = rational/valuation case with r>1.
    G1      = Mordell corridor (z=3, y=x+2, 3|x).
    """
    r_info = s.ratio()
    if r_info is None:
        return "none"
    r, small, large, x_small, x_large = r_info
    diagonal = (r == 1)

    # G0 : grand terme a x_large==z, petit a x_small==x_large+1
    if x_large == s.z and x_small == x_large + 1:
        return "G0d" if diagonal else "G0"

    # Gv : exposants égaux, z==x+1
    if x_small == x_large and s.z == x_small + 1:
        return "Gvd" if diagonal else "Gv"

    # G1 : z=3, |y-x|=2, min(x,y)%3==0
    nx, ny = min(s.x, s.y), max(s.x, s.y)
    if s.z == 3 and ny - nx == 2 and nx % 3 == 0:
        return "G1"

    return "none"


def classify(A, B, x, y, C, z):
    """
    Returns (macro_family, sub_family, geometry, proof_status, certificate).

    Priority order:
      1. Mordell-like (5,6,3) → M3.G1 with elliptic note
      2. M1: symmetric (x==y), geometry computed but often none
      3. M3: binomials (integer ratio, x≠y)
      4. M2: residuals (all remaining non-symmetric non-binomial)
      5. F0 : sporadic
    """
    s = Solution(A, B, x, y, C, z)
    geo = _geometry(s)  # calculé une seule fois, transmis à tous les return

    # ── Mordell-like (5,6,3) → M3.G1 avec certificat elliptique ────────
    for (nx, ny) in [(s.x, s.y), (s.y, s.x)]:
        if nx == 5 and ny == 6 and s.z == 3:
            Bv = s.B if ny == 6 else s.A
            Av = s.A if nx == 5 else s.B
            if Bv > 0 and Av % Bv == 0:
                a_m = Av // Bv
                k3 = 3 * a_m**2 + 8
                k = iroot(k3, 3)
                if k is not None:
                    cert = (f"M3.G1·mordell-sporadic(a={a_m},k={k}) : "
                            f"3·{a_m}²+8={k}³ — isolated point on Y²=X³-r^5")
                    return "M3", "M3.G1", "G1", "identity_proven", cert

    # ── M1 : Symétriques (x==y) — géométrie transmise ──────────────────
    if s.x == s.y:

        # M1.1 sym-pure : A==B (bases normalisées)
        if s.An == s.Bn and s.xn == s.yn:
            try:
                if 2 * s.An**s.xn == s.Cn**s.zn:
                    cert = (f"M1.1·sym-pure(A={s.An},k={s.xn}) : "
                            f"2·{s.An}^{s.xn}={s.Cn}^{s.zn}")
                    return "M1", "M1.1", geo, "identity_proven", cert
            except OverflowError:
                pass

        # M1.2 next-power : (aS)^k+(bS)^k=S^(k+1), S=gcd=C
        if s.z == s.x + 1 and s.C > 0 and s.A % s.C == 0 and s.B % s.C == 0:
            a, b = s.A // s.C, s.B // s.C
            try:
                if a**s.x + b**s.x == s.C:
                    cert = (f"M1.2·next-power(a={a},b={b},k={s.x},S={s.C}) : "
                            f"({a}S)^{s.x}+({b}S)^{s.x}=S^{s.x+1}")
                    return "M1", "M1.2", geo, "identity_proven", cert
            except OverflowError:
                pass

        # M1.3 somme-scalee : (a^n+b^n)·g^n = C^z
        if s.g > 0:
            try:
                if (s.a**s.x + s.b**s.y) * s.g**s.x == s.C**s.z:
                    sc = s.a**s.x + s.b**s.y
                    cert = (f"M1.3·scaled-sum(a={s.a},b={s.b},"
                            f"n={s.x},g={s.g},s={sc})")
                    return "M1", "M1.3", geo, "identity_proven", cert
            except OverflowError:
                pass

        # M1.4 lift-pondéré : a^k+b^k = g·c^(k+1) ou résiduel
        if s.g > 0:
            if s.z == s.x + 1:
                try:
                    if s.a**s.x + s.b**s.y == s.g * s.c**s.z:
                        cert = (f"M1.4·weighted-lift(a={s.a},b={s.b},"
                                f"c={s.c},g={s.g},k={s.x})")
                        return "M1", "M1.4", geo, "identity_proven", cert
                except OverflowError:
                    pass
            exp = s.z - s.x
            if exp > 0:
                try:
                    if s.a**s.x + s.b**s.y == s.g**exp * s.c**s.z:
                        cert = (f"M1.4·weighted-lift-residual(a={s.a},b={s.b},"
                                f"c={s.c},g={s.g},coeff=g^{exp})")
                        return "M1", "M1.4", geo, "identity_proven", cert
                except OverflowError:
                    pass

    # ── M3 : Binômes (ratio A/B entier, x≠y) ────────────────────────────
    r_info = s.ratio()
    if r_info is not None and s.x != s.y:
        r, small, large, x_small, x_large = r_info
        am = s.a if s.A >= s.B else s.b
        aM = s.b if s.A >= s.B else s.a

        # M3.S binome-sym : A==B, x≠y (dx=1)
        if s.A == s.B:
            xm, xM = min(s.x, s.y), max(s.x, s.y)
            if xM - xm == 1:
                cert = f"M3.S·binom-sym(A={s.A},dx=1)"
                return "M3", "M3.S", geo, "identity_proven", cert

        # G0/Gv/G1 or M3.R (off-corridor)
        if geo in ("G0", "Gv", "G1"):
            sub = f"M3.{geo}"
            # G1 → sage_pending (à certifier par SageMath)
            proof = "sage_integral_points_pending" if geo == "G1" else "identity_proven"
            cert = (f"{sub}·binom(a={am},b={aM},r={r},"
                    f"g={s.g},x={s.x},y={s.y},z={s.z})")
            return "M3", sub, geo, proof, cert

        cert = (f"M3.R·binom(a={am},b={aM},r={r},"
                f"g={s.g},x={s.x},y={s.y},z={s.z})")
        return "M3", "M3.R", geo, "identity_proven", cert

    # ── M2 : Résiduels (α≠β, non symétrique, non binomial) ──────────────
    if s.g > 0:
        # M2.1 résiduel-direct : signature (0,β,0)
        if s.alpha == 0 and s.gamma == 0 and s.beta >= 1:
            try:
                ge = s.g**s.beta
                if s.a**s.x + s.b**s.y * ge == s.c**s.z:
                    cert = (f"M2.1·direct-residual(a={s.a},b={s.b},c={s.c},"
                            f"α=0,β={s.beta},γ=0,g={s.g})")
                    return "M2", "M2.1", geo, "identity_proven", cert
                if s.a**s.x * ge + s.b**s.y == s.c**s.z:
                    cert = (f"M2.1·direct-residual-mirror(a={s.a},b={s.b},c={s.c},"
                            f"g={s.g},β={s.beta})")
                    return "M2", "M2.1", geo, "identity_proven", cert
            except OverflowError:
                pass

        # M2.2 résiduel-mixte : équation résiduelle générale
        try:
            lhs = (s.a**s.x * s.g**s.alpha +
                   s.b**s.y * s.g**s.beta)
            rhs = s.c**s.z * s.g**s.gamma
            if lhs == rhs:
                cert = (f"M2.2·mixed-residual(a={s.a},b={s.b},c={s.c},"
                        f"α={s.alpha},β={s.beta},γ={s.gamma},g={s.g})")
                return "M2", "M2.2", geo, "identity_proven", cert
        except OverflowError:
            pass

    # ── F0 : Sporadique ──────────────────────────────────────────────────
    return "F0", "F0", geo, "unclassified", (
        f"F0·sporadic(A={A},B={B},x={x},y={y},C={C},z={z},g={s.g},"
        f"α={s.alpha},β={s.beta},γ={s.gamma})"
    )


# ==============================================================================
# SCAN
# ==============================================================================

def scan(resume: bool):
    log("=" * 60)
    log(f"SCAN — grille x={X_RANGE} y={Y_RANGE} z={Z_RANGE}")
    log(f"       espace [1,{SCAN_AMAX:,}]²  resume={resume}")
    log("=" * 60)

    if not resume and os.path.exists(CANDIDATES_DIR):
        import shutil
        shutil.rmtree(CANDIDATES_DIR)
    os.makedirs(CANDIDATES_DIR, exist_ok=True)

    combos = [(x, y, z) for x in X_RANGE for y in Y_RANGE for z in Z_RANGE if x <= y]
    t0 = time.time()

    for i, (x, y, z) in enumerate(combos):
        out_file  = os.path.join(CANDIDATES_DIR, f"cand_{x}_{y}_{z}.txt")
        ckpt_file = os.path.join(CANDIDATES_DIR, f"ckpt_{x}_{y}_{z}.bin")
        if resume and os.path.exists(out_file) and os.path.exists(ckpt_file):
            continue

        cmd = [
            BEAL_BIN,
            "--x", str(x),
            "--y", str(y),
            "--z", str(z),
            "--amax", str(SCAN_AMAX),
            "--bmax", str(SCAN_BMAX),
            "--out", out_file,
            "--resume", ckpt_file,
        ]
        log(f"  [{i+1}/{len(combos)}] x={x} y={y} z={z}")
        subprocess.run(cmd, capture_output=True, text=True)

    log(f"Scan terminé en {(time.time() - t0) / 60:.1f} min")


# ==============================================================================
# VERIFY + CLASSIFY → BIBLE
# ==============================================================================

def verify_and_classify():
    log("=" * 60)
    log(f"VERIFY + CLASSIFY → {BIBLE_FILE}")
    log("=" * 60)

    if not os.path.exists(CANDIDATES_DIR):
        log(f"  Dossier {CANDIDATES_DIR}/ introuvable.")
        return

    solutions = []
    seen = set()

    for fname in sorted(os.listdir(CANDIDATES_DIR)):
        if not fname.startswith("cand_") or not fname.endswith(".txt"):
            continue

        parts = fname.replace("cand_", "").replace(".txt", "").split("_")
        if len(parts) != 3:
            continue
        x, y, z = map(int, parts)

        fpath = os.path.join(CANDIDATES_DIR, fname)
        if os.path.getsize(fpath) == 0:
            continue

        pairs = []
        with open(fpath) as f:
            for line in f:
                ab = line.strip().split()
                if len(ab) >= 2:
                    try:
                        pairs.append((int(ab[0]), int(ab[1])))
                    except ValueError:
                        continue

        if not pairs:
            continue

        log(f"  x={x} y={y} z={z} : {len(pairs)} candidats")

        for A, B in pairs:
            for cur_x, cur_y in [(x, y), (y, x)]:
                try:
                    S = A ** cur_x + B ** cur_y
                except OverflowError:
                    continue

                # Chercher z_found : tester d'abord le z du scan (forme géométrique),
                # puis fallback sur is_perfect_power (forme canonique max exposant).
                # Cela préserve la géométrie : 21952^3 reste z=3 même si 21952=28^3.
                C_raw = iroot(S, z)
                if C_raw is not None and C_raw > 1:
                    # z du scan est valide → garder z brut (géométrie préservée)
                    C, z_found = C_raw, z
                else:
                    res = is_perfect_power(S)
                    if res is None:
                        continue
                    C, z_found = res
                g = math.gcd(math.gcd(A, B), C)

                # Déduplication : canonicalise (A^x, B^y) par ordre lexicographique
                # pour éviter de confondre A^3+B^5 avec A^5+B^3
                left1 = (A, cur_x)
                left2 = (B, cur_y)
                if left2 < left1:
                    left1, left2 = left2, left1
                key = (left1, left2, C, z_found)
                if key in seen:
                    continue
                seen.add(key)

                macro_family, sub_family, geometry, proof_status, certificate = \
                    classify(A, B, cur_x, cur_y, C, z_found)

                # Calcul des colonnes dérivées
                sol_obj = Solution(A, B, cur_x, cur_y, C, z_found)

                solutions.append({
                    "A": A, "B": B, "x": cur_x, "y": cur_y,
                    "C": C, "z": z_found,
                    "C_canonical": sol_obj.C_canonical,
                    "z_canonical": sol_obj.z_canonical,
                    "gcd": g,
                    "a": sol_obj.a, "b": sol_obj.b, "c": sol_obj.c,
                    "m": sol_obj.m,
                    "alpha": sol_obj.alpha,
                    "beta":  sol_obj.beta,
                    "gamma": sol_obj.gamma,
                    "macro_family": macro_family,
                    "sub_family":   sub_family,
                    "regime":       SUB_REGIMES.get(sub_family, "?"),
                    "geometry":     geometry,
                    "proof_status": proof_status,
                    "certificate":  certificate,
                })

    solutions.sort(key=lambda s: (s["macro_family"], s["sub_family"], s["A"] * s["B"]))

    with open(BIBLE_FILE, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "A", "B", "x", "y", "C", "z",
                "C_canonical", "z_canonical",
                "gcd", "a", "b", "c",
                "m", "alpha", "beta", "gamma",
                "macro_family", "sub_family", "regime", "geometry",
                "proof_status", "certificate",
            ],
        )
        writer.writeheader()
        for s in solutions:
            writer.writerow(s)

    log(f"  {len(solutions)} solutions uniques → {BIBLE_FILE}")
    _print_report(solutions)


# ==============================================================================
# RAPPORT
# ==============================================================================

def _print_report(solutions):
    by_macro = defaultdict(lambda: defaultdict(list))
    for s in solutions:
        by_macro[s["macro_family"]][s["sub_family"]].append(s)

    print()
    print("=" * 72)
    print("BIBLE v4 — SUMMARY BY MACRO-FAMILY")
    print("=" * 72)
    print()
    print("  Residual signature : a^x·g^α + b^y·g^β = c^z·g^γ")
    print("  Regimes : ● prime  ◑ mixed  ○ even  ◈ both")

    total = len(solutions)
    MACRO_ORDER = ["M1", "M2", "M3", "F0"]

    for macro in MACRO_ORDER:
        if macro not in by_macro:
            continue
        subs = by_macro[macro]
        macro_total = sum(len(v) for v in subs.values())
        print(f"\n  {'─'*68}")
        print(f"  {macro} — {MACRO_LABELS.get(macro,'?').upper()} : {macro_total} solutions")

        for sub in sorted(subs.keys()):
            members = subs[sub]
            regime = SUB_REGIMES.get(sub, "?")
            rmark = {"prime":"●","even":"○","mixed":"◑","both":"◈"}.get(regime,"?")
            sub_label = SUB_LABELS.get(sub, sub)
            print(f"\n    {sub} — {sub_label} : {len(members)} solutions  {rmark} {regime}")

            # Afficher la signature (α,β,γ) commune si uniforme
            sigs = set((s["alpha"], s["beta"], s["gamma"]) for s in members)
            if len(sigs) == 1:
                a,b,g = list(sigs)[0]
                print(f"      signature : (α={a}, β={b}, γ={g})")

            top = sorted(members, key=lambda s: s["A"] * s["B"], reverse=True)[:3]
            for s in top:
                geo = f"  [{s['geometry']}]" if s["geometry"] != "none" else ""
                print(
                    f"      {s['A']:>12,}^{s['x']} + {s['B']:>12,}^{s['y']} "
                    f"= {s['C']:>12,}^{s['z']}{geo}"
                )
                print(f"      → {s['certificate'][:70]}")

    print()
    print(f"  Total: {total} solutions")
    print()
    _print_elliptic_report(solutions)



def _print_elliptic_report(solutions):
    """
    Rapport géométrique — utilise le champ 'geometry' de la v4.
    G0 / Gv / G1 sont déjà calculés par classify().
    """
    g0  = {}   # r > 1
    g0d = {}   # r == 1 (diagonal)
    gv  = {}
    gvd = {}
    g1  = {}

    for s in solutions:
        x, y, z = s["x"], s["y"], s["z"]
        geo = s.get("geometry", "none")
        if geo == "none":
            continue
        A, B, C = s["A"], s["B"], s["C"]
        if A % B == 0:
            big, small, nx, ny = A, B, x, y
        elif B % A == 0:
            big, small, nx, ny = B, A, y, x
        else:
            continue
        r = big // small
        key = (nx, ny, z, r)
        if geo == "G0":   g0.setdefault(key, []).append((A, B, C))
        elif geo == "G0d": g0d.setdefault(key, []).append((A, B, C))
        elif geo == "Gv":  gv.setdefault(key, []).append((A, B, C))
        elif geo == "Gvd": gvd.setdefault(key, []).append((A, B, C))
        elif geo == "G1":  g1.setdefault(key, []).append((A, B, C))

    if not g0 and not g0d and not gv and not gvd and not g1:
        return

    print()
    print("=" * 72)
    print("GEOMETRY OF SOLUTIONS WITH RATIO A=rB")
    print("=" * 72)

    def _print_corridor(d, label, formula, note=""):
        if not d:
            return
        print()
        print(f"  {label}")
        print("  ─" * 35)
        if note:
            print(f"  {note}")
        by_trip = {}
        for (nx, ny, nz, r), exs in sorted(d.items()):
            by_trip.setdefault((nx, ny, nz), []).append((r, len(exs)))
        for (nx, ny, nz), rs in sorted(by_trip.items()):
            total = sum(n for _, n in rs)
            rlist = ", ".join(f"r={r}({n})" for r, n in sorted(rs)[:5])
            more = f" +{len(rs)-5} more" if len(rs) > 5 else ""
            print(f"    ({nx},{ny},{nz})  {total} solutions — {rlist}{more}")

    _print_corridor(g0d, "G0d — RATIONAL DIAGONAL CORRIDOR  r=1  →  B = X^x - 1",
                    "", "Case A=B: B=X^x-1, C=XB (e.g. 2B^x=C^x)")
    _print_corridor(g0,  "G0  — RATIONAL CORRIDOR  (x, x+1, x)  →  B = X^x - r^x",
                    "", "Universal parametrization : B=X^x-r^x, A=rB, C=X·B")
    _print_corridor(gvd, "Gvd — VALUATION DIAGONAL CORRIDOR  r=1  →  2·B^x = C^(x+1)",
                    "", "Case A=B: 2B^x=C^(x+1)")
    _print_corridor(gv,  "Gv  — VALUATION CORRIDOR  (x, x, x+1)  →  (r^x+1)·B^x = C^(x+1)",
                    "", "Structure: D=r^x+1, B=B₀·t^(x+1), C=C₀·t^x")

    # ── G1 ────────────────────────────────────────────────────────────────
    if g1:
        print()
        print("  G1 — ELLIPTIC CORRIDOR  (3m, 3m+2, 3)  →  Y² = X³ - r^x  (Mordell)")
        print("  ─" * 35)
        print("  Lift : A=rY, B=Y, C=X·Y^(x/3)")
        print()
        sage_lines = []
        by_trip = {}
        for (nx, ny, nz, r), exs in sorted(g1.items()):
            by_trip.setdefault((nx, ny, nz), []).append((r, exs))
        for (nx, ny, nz), rs in sorted(by_trip.items()):
            print(f"    ({nx},{ny},{nz})  →  Y²=X³-r^{nx}")
            for r, exs in sorted(rs):
                n_mordell = -(r**nx)
                print(f"      r={r:6}  n={n_mordell}  ({len(exs)} sol. GPU)")
                A0, B0, C0 = exs[0]
                print(f"        ex: {A0}^{nx}+{B0}^{ny}={C0}^{nz}")
                sage_lines.append(
                    f"EllipticCurve([0, {n_mordell}]).integral_points()"
                    f"  # ({nx},{ny},{nz}) r={r}"
                )
        if sage_lines:
            print()
            print("  ┌─ SAGEMATH (sagecell.sagemath.org) " + "─" * 34)
            print("  │  # Commands to certify the complete list of integer points")
            print("  │  # (sage_integral_points_pending → run to certify)")
            print("  │")
            for line in sage_lines:
                print(f"  │  {line}")
            print("  └" + "─" * 70)




def main():
    resume      = "--resume"      in sys.argv
    verify_only = "--verify-only" in sys.argv
    inject      = "--inject"      in sys.argv  # lit candidates.txt racine

    print(r"""
╔══════════════════════════════════════════════════════════════════════╗
║  BEAL COMMANDER v4 — Scan → Verify → Classify → Bible               ║
╚══════════════════════════════════════════════════════════════════════╝""")
    print(f"  x: {X_RANGE}  y: {Y_RANGE}  z: {Z_RANGE}")
    print(f"  Window : [1,{SCAN_AMAX:,}]²")

    if inject:
        # Mode injection : lit candidates.txt à la racine + x y z depuis args
        # Usage : python3 beal_commander_v4.py --inject --x 3 --y 5 --z 3
        xi = yi = zi = None
        for i, arg in enumerate(sys.argv):
            if arg == "--x" and i+1 < len(sys.argv): xi = int(sys.argv[i+1])
            if arg == "--y" and i+1 < len(sys.argv): yi = int(sys.argv[i+1])
            if arg == "--z" and i+1 < len(sys.argv): zi = int(sys.argv[i+1])
        if xi is None or yi is None or zi is None:
            print("Usage : --inject --x X --y Y --z Z")
            sys.exit(1)
        # Copier candidates.txt dans le dossier candidates/ avec le bon nom
        import shutil
        os.makedirs(CANDIDATES_DIR, exist_ok=True)
        dst = os.path.join(CANDIDATES_DIR, f"cand_{xi}_{yi}_{zi}.txt")
        shutil.copy("candidates.txt", dst)
        print(f"  Mode   : INJECT candidates.txt → {dst}")
        print(f"  Output : {BIBLE_FILE} (append)")
        print()
        verify_and_classify()

    elif verify_only:
        print(f"  Mode   : VERIFY ONLY (pas de scan)")
        print(f"  Output : {BIBLE_FILE}")
        print()
        verify_and_classify()

    else:
        print(f"  Mode   : {'reprise (--resume)' if resume else 'nouveau départ'}")
        print(f"  Output : {BIBLE_FILE}")
        print()
        scan(resume)
        verify_and_classify()


if __name__ == "__main__":
    main()
