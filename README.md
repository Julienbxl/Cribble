# Cribble — GPU Explorer for A^x + B^y = C^z

A GPU-accelerated tool to enumerate and classify integer solutions of the generalized Fermat equation:

```
A^x + B^y = C^z      with  A, B, C, x, y, z positive integers,  x, y, z ≥ 3
```

This equation is at the heart of the **Beal conjecture** (1993), which states that every solution has gcd(A, B, C) > 1. Cribble does not attempt to prove the conjecture — it systematically explores the solution space and builds a structured catalog of all found solutions.

---

## What it does

**Scan (GPU — `beal_cpu`)** — A CUDA kernel sweeps A, B over a configurable range [1, N]² for each exponent triple (x, y, z). A multi-level sieve (L2 cache + VRAM bitmap + K2 primorial filter) rejects ~99.9999999% of pairs before any expensive power computation, making large scans practical on a consumer GPU.

**Classify (Python — `beal_commander_v4.py`)** — Each GPU candidate is verified arithmetically, then placed into a three-level taxonomy:

| Level | What it captures |
|-------|-----------------|
| Macro-family (M1/M2/M3) | The algebraic generating mechanism |
| Sub-family | Specific parametric form |
| Geometry (G0/Gv/G1) | Corridor structure when A = r·B |

Results are written to `bible_v4.csv` with full certificates and residual signatures.

---

## Taxonomy

After extracting g = gcd(A, B, C) and writing a = A/g, b = B/g, c = C/g, the equation becomes:

```
a^x · g^α  +  b^y · g^β  =  c^z · g^γ      where  α=x-m, β=y-m, γ=z-m, m=min(x,y,z)
```

The triple **(α, β, γ)** is the residual signature — the fundamental invariant.

### Macro-families

**M1 — Symmetric** (x == y): both left terms share the same power of g.
- M1.1 `sym-pure` — A = B, x = y → 2A^k = C^z
- M1.2 `next-power` — (aS)^k + (bS)^k = S^(k+1) where S = a^k + b^k, infinite by construction
- M1.3 `somme-scalee` — (a^n + b^n)·g^n = C^z

**M2 — Residual** (α ≠ β): asymmetric powers of g.
- M2.1 `residuel-direct` — signature (0, e, 0): a^x + b^y·g^e = c^x
- M2.2 `residuel-mixte` — other (α, β, γ)

**M3 — Binomial** (A/B integer ratio r):
- M3.G0 — rational corridor (x, x+1, x): B = X^x − r^x, universal parametrization
- M3.Gv — valuation corridor (x, x, x+1): (r^x + 1)·B^x = C^(x+1)
- M3.G1 — **Mordell corridor** (3m, 3m+2, 3): reduces to Y² = X³ − r^x
- M3.R — binomial, no identified corridor
- M3.S — symmetric binomial, A = B, x ≠ y

**F0 — Sporadic**: catch-all for genuinely unclassified solutions (observed count: 0).

### Regime labels

Each sub-family is tagged with its activation regime:
- `● prime` — appears with prime exponents only
- `◑ mixed` — appears with both prime and composite exponents
- `○ even` — requires at least one composite (even) exponent

---

## The Mordell corridor (G1) — highlight result

For any exponent triple (x, y, 3) with A = r·B and 3 | x, substituting A = rY, B = Y and setting C = X·Y^(x/3) yields the universal reduction:

```
Y² = X³ − rˣ
```

a Mordell elliptic curve. Every integer point (X, Y) lifts to a Beal solution:

```
A = rY,   B = Y,   C = X·Y^(x/3)
```

This holds for all G1 families observed:

| Triple | Curve | Example |
|--------|-------|---------|
| (3, 5, 3) | Y² = X³ − r³ | 91³ + 13⁵ = 104³  (r=7) |
| (5, 6, 3) | Y² = X³ − r⁵ | 518616⁵ + 24696⁶ = 3354408288³  (r=21) |

For each ratio r discovered by the GPU scan, SageMath certifies the *complete* list of integer points — transforming empirical search into provable enumeration for each Mordell family:

```python
EllipticCurve([0, -r^3]).integral_points()   # e.g. r = 6, 7, 11, 23, ...
```

### Taxonomy completeness

On the full search space (A, B ≤ 1M, exponents 3–10), **F0 = 0**: every solution with gcd(A, B, C) > 1 falls into M1, M2, or M3. No solution required the sporadic catch-all family. This empirically validates the completeness of the three-family taxonomy over this search window.

### Disclaimer

This project searches exclusively for solutions where x, y, z ≥ 3. All found solutions satisfy gcd(A, B, C) > 1, consistent with the Beal conjecture. It does not contradict Fermat's Last Theorem: the case x = y = z = n with gcd = 1 is excluded by construction, and Pythagorean triples (z = 2) are outside the search space.

---

## Results (scan A, B ≤ 1M, exponents 3–10)

| Family | Count | Regime |
|--------|-------|--------|
| M1.1 sym-pure | 530 | ● prime |
| M1.2 next-power | 602 | ● prime |
| M1.3 somme-scalee | 445 | ◑ mixed |
| M2.1 residuel-direct | 491 | ◑ mixed |
| M2.2 residuel-mixte | 33 | ◑ mixed |
| M3.G0 binome-G0 | 1201 | ○ even |
| M3.G1 binome-G1 | 57 | ◑ mixed |
| M3.R binome-hors-corridor | 743 | ◑ mixed |
| M3.S binome-sym | 323 | ○ even |
| **Total** | **4425** | |

Geometric breakdown (solutions with integer ratio A = r·B):

| Corridor | Description | Solutions |
|----------|-------------|-----------|
| G0d | Rational diagonal r=1 | 165 |
| G0 | Rational r>1 : B = X^x − r^x | 1201 |
| Gvd | Valuation diagonal r=1 | 65 |
| Gv | Valuation r>1 : (r^x+1)·B^x = C^(x+1) | 162 |
| **G1** | **Mordell: Y² = X³ − r^x** | **57** |

The G1 corridor covers 35 distinct Mordell curves. F0 (sporadic / unclassified) = 0.

## Validation

The scan was cross-validated against SageMath for all 35 Mordell curves in the G1 corridor:

- **47/47** in-window solutions confirmed — scan is complete
- **0 false positives** — every GPU candidate satisfies A^x + B^y = C^z exactly
- **21 additional solutions** certified by SageMath beyond the GPU window (A > 1M), available in `g1_out_of_window.csv`

---

## Usage

```bash
# 1. Build the GPU kernel
make

# 2. Run a scan
./beal_cpu --x 3 --y 5 --z 3 --amin 1 --amax 1000000 --bmin 1 --bmax 1000000

# 3. Classify and build the bible
python3 beal_commander_v4.py

# 4. Resume an interrupted scan
python3 beal_commander_v4.py --resume

# 5. Verify only (no scan)
python3 beal_commander_v4.py --verify-only

# 6. Inject a targeted scan result
./beal_cpu --x 3 --y 5 --z 3 --amin 556808 --amax 568056 --bmin 5354 --bmax 5462
python3 beal_commander_v4.py --inject --x 3 --y 5 --z 3
```

### Key configuration (beal_commander_v4.py)

```python
X_RANGE    = list(range(3, 11))   # exponent ranges
Y_RANGE    = list(range(3, 11))
Z_RANGE    = list(range(3, 11))
SCAN_AMAX  = 1_000_000            # search window
SCAN_BMAX  = 1_000_000
BIBLE_FILE = "bible_v4.csv"
```

---

## Output: bible_v4.csv

Each row contains:

| Column | Description |
|--------|-------------|
| A, B, x, y, C, z | Raw solution |
| C_canonical, z_canonical | Normalized form (e.g. 21952^3 → 28^9) |
| gcd, a, b, c | After gcd extraction |
| m, alpha, beta, gamma | Residual signature (α, β, γ) |
| macro_family, sub_family | Taxonomy |
| regime | prime / mixed / even |
| geometry | none / G0 / G0d / Gv / Gvd / G1 |
| proof_status | identity_proven / sage_integral_points_pending / sage_certified |
| certificate | Algebraic certificate string |

---

## Requirements

- CUDA-capable GPU (tested on RTX 5060, sm_120)
- CUDA toolkit ≥ 12.0
- Python ≥ 3.10
- SageMath (optional, for Mordell curve certification — https://sagecell.sagemath.org)

---

## Hardware

Developed and tested on RTX 5060 (Blackwell, sm_120) under WSL2 on Windows. Typical scan throughput: ~2400 MKeys/s for (3,5,3) over [1, 1M]².
