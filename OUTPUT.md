========================================================================
BIBLE v4 — SUMMARY BY MACRO-FAMILY
========================================================================

  Residual signature : a^x·g^α + b^y·g^β = c^z·g^γ
  Regimes : ● prime  ◑ mixed  ○ even  ◈ both

  ────────────────────────────────────────────────────────────────────
  M1 — SYMMETRIC : 1577 solutions

    M1.1 — sym-pure : 530 solutions  ● prime
           995,328^7 +      995,328^7 =  286,654,464^5
      → M1.1·sym-pure(A=995328,k=7) : 2·995328^7=286654464^5
           986,078^5 +      986,078^5 = 12,308,225,596^3
      → M1.1·sym-pure(A=986078,k=5) : 2·986078^5=12308225596^3
           953,312^4 +      953,312^4 =  118,210,688^3
      → M1.1·sym-pure(A=953312,k=4) : 2·953312^4=118210688^3

    M1.2 — next-power : 602 solutions  ● prime
      signature : (α=0, β=0, γ=1)
           953,316^3 +      882,700^3 =       35,308^4
      → M1.2·next-power(a=27,b=25,k=3,S=35308) : (27S)^3+(25S)^3=S^4
           937,678^4 +      870,701^4 =       66,977^5
      → M1.2·next-power(a=14,b=13,k=4,S=66977) : (14S)^4+(13S)^4=S^5
           955,332^3 +      784,737^3 =       34,119^4
      → M1.2·next-power(a=28,b=23,k=3,S=34119) : (28S)^3+(23S)^3=S^4

    M1.3 — scaled-sum : 445 solutions  ◑ mixed
           948,992^3 +      868,112^3 =       35,048^4
      → M1.3·scaled-sum(a=352,b=322,n=3,g=2696,s=77000456)
           830,609^3 +      789,308^3 =       32,123^4
      → M1.3·scaled-sum(a=181,b=172,n=3,g=4589,s=11018189)
           865,059^5 +      731,973^5 = 8,855,941,698^3
      → M1.3·scaled-sum(a=13,b=11,n=5,g=66543,s=532344)

  ────────────────────────────────────────────────────────────────────
  M2 — RESIDUAL : 524 solutions

    M2.1 — direct-residual : 491 solutions  ◑ mixed
           660,022^3 +      990,033^4 =   98,673,289^3
      → M2.1·direct-residual(a=2,b=3,c=299,α=0,β=1,γ=0,g=330011)
           864,045^3 +      729,638^4 =   65,686,621^3
      → M2.1·direct-residual(a=45,b=38,c=3421,α=0,β=1,γ=0,g=19201)
           941,535^3 +      564,921^4 =   46,700,136^3
      → M2.1·direct-residual(a=5,b=3,c=248,α=0,β=1,γ=0,g=188307)

    M2.2 — mixed-residual : 33 solutions  ◑ mixed
           979,985^6 +      587,991^7 = 28,964,777,302,786^3
      → M2.2·mixed-residual(a=5,b=3,c=147781738,α=3,β=4,γ=0,g=195997)
           756,268^6 +      567,201^7 = 26,631,016,214,305^3
      → M2.2·mixed-residual(a=4,b=3,c=140854915,α=3,β=4,γ=0,g=189067)
           866,110^6 +      346,444^7 = 8,431,647,020,804^3
      → M2.2·mixed-residual(a=5,b=2,c=48675382,α=3,β=4,γ=0,g=173222)

  ────────────────────────────────────────────────────────────────────
  M3 — BINOMIAL : 2324 solutions

    M3.G0 — binomial-G0 : 1201 solutions  ○ even
      signature : (α=0, β=1, γ=0)
           986,062^3 +      493,031^4 =   38,949,449^3  [G0]
      → M3.G0·binom(a=2,b=1,r=2,g=493031,x=3,y=4,z=3)
           949,088^3 +      474,544^4 =   37,014,432^3  [G0]
      → M3.G0·binom(a=2,b=1,r=2,g=474544,x=3,y=4,z=3)
           913,920^4 +      456,960^5 =   11,880,960^4  [G0]
      → M3.G0·binom(a=2,b=1,r=2,g=456960,x=4,y=5,z=4)

    M3.G1 — binomial-G1 : 57 solutions  ◑ mixed
           518,616^5 +       24,696^6 = 3,354,408,288^3  [G1]
      → M3.G1·mordell-sporadic(a=21,k=11) : 3·21²+8=11³ — isolated point on Y²
           428,064^3 +       15,288^5 =    9,417,408^3  [G1]
      → M3.G1·binom(a=28,b=1,r=28,g=15288,x=3,y=5,z=3)
           975,912^3 +        6,216^5 =    2,169,384^3  [G1]
      → M3.G1·binom(a=157,b=1,r=157,g=6216,x=3,y=5,z=3)

    M3.R — binomial-off-corridor : 743 solutions  ◑ mixed
           985,950^6 +      492,975^7 = 19,198,923,699,375^3
      → M3.R·binom(a=2,b=1,r=2,g=492975,x=6,y=7,z=3)
           492,280^4 +      984,560^5 =   31,013,640^4
      → M3.R·binom(a=2,b=1,r=2,g=492280,x=4,y=5,z=4)
           948,976^6 +      474,488^7 = 17,560,831,247,232^3
      → M3.R·binom(a=2,b=1,r=2,g=474488,x=6,y=7,z=3)

    M3.S — binomial-sym : 323 solutions  ○ even
           999,999^3 +      999,999^4 =   99,999,900^3  [G0d]
      → M3.S·binom-sym(A=999999,dx=1)
           999,999^6 +      999,999^7 = 99,999,800,000,100^3
      → M3.S·binom-sym(A=999999,dx=1)
           999,999^6 +      999,999^7 =    9,999,990^6  [G0d]
      → M3.S·binom-sym(A=999999,dx=1)

  Total: 4425 solutions


========================================================================
GEOMETRY OF SOLUTIONS WITH RATIO A=rB
========================================================================

  G0d — RATIONAL DIAGONAL CORRIDOR  r=1  →  B = X^x - 1
  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─
  Case A=B: B=X^x-1, C=XB (e.g. 2B^x=C^x)
    (3,4,3)  99 solutions — r=1(99)
    (4,5,4)  30 solutions — r=1(30)
    (5,6,5)  14 solutions — r=1(14)
    (6,7,6)  9 solutions — r=1(9)
    (7,8,7)  6 solutions — r=1(6)
    (8,9,8)  4 solutions — r=1(4)
    (9,10,9)  3 solutions — r=1(3)

  G0  — RATIONAL CORRIDOR  (x, x+1, x)  →  B = X^x - r^x
  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─
  Universal parametrization : B=X^x-r^x, A=rB, C=X·B
    (3,4,3)  956 solutions — r=2(77), r=3(66), r=4(59), r=5(53), r=6(49) +63 more
    (4,5,4)  165 solutions — r=2(24), r=3(21), r=4(18), r=5(16), r=6(14) +15 more
    (5,6,5)  48 solutions — r=2(11), r=3(9), r=4(8), r=5(6), r=6(5) +5 more
    (6,7,6)  17 solutions — r=2(6), r=3(5), r=4(3), r=5(2), r=6(1)
    (7,8,7)  8 solutions — r=2(4), r=3(3), r=4(1)
    (8,9,8)  4 solutions — r=2(3), r=3(1)
    (9,10,9)  3 solutions — r=2(2), r=3(1)

  Gvd — VALUATION DIAGONAL CORRIDOR  r=1  →  2·B^x = C^(x+1)
  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─
  Case A=B: 2B^x=C^(x+1)
    (3,3,4)  26 solutions — r=1(26)
    (4,4,5)  13 solutions — r=1(13)
    (5,5,6)  8 solutions — r=1(8)
    (6,6,7)  6 solutions — r=1(6)
    (7,7,8)  5 solutions — r=1(5)
    (8,8,9)  4 solutions — r=1(4)
    (9,9,10)  3 solutions — r=1(3)

  Gv  — VALUATION CORRIDOR  (x, x, x+1)  →  (r^x+1)·B^x = C^(x+1)
  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─
  Structure: D=r^x+1, B=B₀·t^(x+1), C=C₀·t^x
    (3,3,4)  91 solutions — r=2(15), r=3(10), r=4(7), r=5(6), r=6(5) +29 more
    (4,4,5)  30 solutions — r=2(7), r=3(5), r=4(3), r=5(3), r=6(2) +9 more
    (5,5,6)  14 solutions — r=2(4), r=3(3), r=4(2), r=5(1), r=6(1) +3 more
    (6,6,7)  9 solutions — r=2(3), r=3(2), r=4(1), r=5(1), r=6(1) +1 more
    (7,7,8)  5 solutions — r=2(2), r=3(1), r=4(1), r=5(1)
    (8,8,9)  4 solutions — r=2(2), r=3(1), r=4(1)
    (9,9,10)  2 solutions — r=2(1), r=3(1)

  G1 — ELLIPTIC CORRIDOR  (3m, 3m+2, 3)  →  Y² = X³ - r^x  (Mordell)
  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─  ─
  Lift : A=rY, B=Y, C=X·Y^(x/3)

    (3,5,3)  →  Y²=X³-r^3
      r=     6  n=-216  (2 sol. GPU)
        ex: 168^3+28^5=280^3
      r=     7  n=-343  (4 sol. GPU)
        ex: 91^3+13^5=104^3
      r=    11  n=-1331  (1 sol. GPU)
        ex: 102564^3+9324^5=4130532^3
      r=    23  n=-12167  (1 sol. GPU)
        ex: 13524^3+588^5=41748^3
      r=    24  n=-13824  (2 sol. GPU)
        ex: 5376^3+224^5=8960^3
      r=    26  n=-17576  (3 sol. GPU)
        ex: 13182^3+507^5=32955^3
      r=    28  n=-21952  (4 sol. GPU)
        ex: 2912^3+104^5=3328^3
      r=    31  n=-29791  (1 sol. GPU)
        ex: 59582^3+1922^5=297910^3
      r=    38  n=-54872  (3 sol. GPU)
        ex: 13718^3+361^5=20577^3
      r=    42  n=-74088  (1 sol. GPU)
        ex: 133770^3+3185^5=691145^3
      r=    47  n=-103823  (1 sol. GPU)
        ex: 25803^3+549^5=40626^3
      r=    54  n=-157464  (2 sol. GPU)
        ex: 40824^3+756^5=68040^3
      r=    55  n=-166375  (2 sol. GPU)
        ex: 78375^3+1425^5=185250^3
      r=    63  n=-250047  (3 sol. GPU)
        ex: 22113^3+351^5=25272^3
      r=    92  n=-778688  (1 sol. GPU)
        ex: 432768^3+4704^5=1335936^3
      r=    96  n=-884736  (1 sol. GPU)
        ex: 172032^3+1792^5=286720^3
      r=   104  n=-1124864  (3 sol. GPU)
        ex: 18824^3+181^5=19005^3
      r=   110  n=-1331000  (1 sol. GPU)
        ex: 42680^3+388^5=44232^3
      r=   111  n=-1367631  (2 sol. GPU)
        ex: 151959^3+1369^5=202612^3
      r=   112  n=-1404928  (2 sol. GPU)
        ex: 93184^3+832^5=106496^3
      r=   118  n=-1643032  (1 sol. GPU)
        ex: 277890^3+2355^5=454515^3
      r=   119  n=-1685159  (1 sol. GPU)
        ex: 122451^3+1029^5=144060^3
      r=   124  n=-1906624  (1 sol. GPU)
        ex: 357492^3+2883^5=625611^3
      r=   140  n=-2744000  (1 sol. GPU)
        ex: 802620^3+5733^5=1886157^3
      r=   143  n=-2924207  (1 sol. GPU)
        ex: 749177^3+5239^5=1634568^3
      r=   150  n=-3375000  (1 sol. GPU)
        ex: 525000^3+3500^5=875000^3
      r=   152  n=-3511808  (2 sol. GPU)
        ex: 438976^3+2888^5=658464^3
      r=   157  n=-3869893  (1 sol. GPU)
        ex: 975912^3+6216^5=2169384^3
      r=   175  n=-5359375  (1 sol. GPU)
        ex: 284375^3+1625^5=325000^3
      r=   182  n=-6028568  (1 sol. GPU)
        ex: 989898^3+5439^5=1789431^3
      r=   188  n=-6644672  (1 sol. GPU)
        ex: 825696^3+4392^5=1300032^3
      r=   189  n=-6751269  (1 sol. GPU)
        ex: 351918^3+1862^5=404054^3
      r=   244  n=-14526784  (1 sol. GPU)
        ex: 907924^3+3721^5=1134905^3
      r=   252  n=-16003008  (1 sol. GPU)
        ex: 707616^3+2808^5=808704^3
      r=   416  n=-71991296  (1 sol. GPU)
        ex: 602368^3+1448^5=608160^3
    (5,6,3)  →  Y²=X³-r^5
      r=    21  n=-4084101  (1 sol. GPU)
        ex: 518616^5+24696^6=3354408288^3

  ┌─ SAGEMATH (sagecell.sagemath.org) ──────────────────────────────────
  │  # Commands to certify the complete list of integer points
  │  # (sage_integral_points_pending → run to certify)
  │
  │  EllipticCurve([0, -216]).integral_points()  # (3,5,3) r=6
  │  EllipticCurve([0, -343]).integral_points()  # (3,5,3) r=7
  │  EllipticCurve([0, -1331]).integral_points()  # (3,5,3) r=11
  │  EllipticCurve([0, -12167]).integral_points()  # (3,5,3) r=23
  │  EllipticCurve([0, -13824]).integral_points()  # (3,5,3) r=24
  │  EllipticCurve([0, -17576]).integral_points()  # (3,5,3) r=26
  │  EllipticCurve([0, -21952]).integral_points()  # (3,5,3) r=28
  │  EllipticCurve([0, -29791]).integral_points()  # (3,5,3) r=31
  │  EllipticCurve([0, -54872]).integral_points()  # (3,5,3) r=38
  │  EllipticCurve([0, -74088]).integral_points()  # (3,5,3) r=42
  │  EllipticCurve([0, -103823]).integral_points()  # (3,5,3) r=47
  │  EllipticCurve([0, -157464]).integral_points()  # (3,5,3) r=54
  │  EllipticCurve([0, -166375]).integral_points()  # (3,5,3) r=55
  │  EllipticCurve([0, -250047]).integral_points()  # (3,5,3) r=63
  │  EllipticCurve([0, -778688]).integral_points()  # (3,5,3) r=92
  │  EllipticCurve([0, -884736]).integral_points()  # (3,5,3) r=96
  │  EllipticCurve([0, -1124864]).integral_points()  # (3,5,3) r=104
  │  EllipticCurve([0, -1331000]).integral_points()  # (3,5,3) r=110
  │  EllipticCurve([0, -1367631]).integral_points()  # (3,5,3) r=111
  │  EllipticCurve([0, -1404928]).integral_points()  # (3,5,3) r=112
  │  EllipticCurve([0, -1643032]).integral_points()  # (3,5,3) r=118
  │  EllipticCurve([0, -1685159]).integral_points()  # (3,5,3) r=119
  │  EllipticCurve([0, -1906624]).integral_points()  # (3,5,3) r=124
  │  EllipticCurve([0, -2744000]).integral_points()  # (3,5,3) r=140
  │  EllipticCurve([0, -2924207]).integral_points()  # (3,5,3) r=143
  │  EllipticCurve([0, -3375000]).integral_points()  # (3,5,3) r=150
  │  EllipticCurve([0, -3511808]).integral_points()  # (3,5,3) r=152
  │  EllipticCurve([0, -3869893]).integral_points()  # (3,5,3) r=157
  │  EllipticCurve([0, -5359375]).integral_points()  # (3,5,3) r=175
  │  EllipticCurve([0, -6028568]).integral_points()  # (3,5,3) r=182
  │  EllipticCurve([0, -6644672]).integral_points()  # (3,5,3) r=188
  │  EllipticCurve([0, -6751269]).integral_points()  # (3,5,3) r=189
  │  EllipticCurve([0, -14526784]).integral_points()  # (3,5,3) r=244
  │  EllipticCurve([0, -16003008]).integral_points()  # (3,5,3) r=252
  │  EllipticCurve([0, -71991296]).integral_points()  # (3,5,3) r=416
  │  EllipticCurve([0, -4084101]).integral_points()  # (5,6,3) r=21
  └──────────────────────────────────────────────────────────────────────
julien@DESKTOP-96478O6:~/Cribble$