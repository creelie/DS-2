"""Exact checks for "The Aperiodic Jordan Algebras of Quasicrystals Are Not Jordan".

Points of Z[tau] are pairs (p, q) = p + q*tau, tau^2 = tau + 1; arithmetic is exact.
Run:  python3 verify_jordan.py      (standard library only)
"""
from fractions import Fraction as Fr
from collections import defaultdict
import random

def tmul(x):                      # tau * (p + q tau) = q + (p + q) tau
    p, q = x
    return (q, p + q)

def qa(x, y):                     # x |- y = tau^2 x - tau y
    a, b = tmul(tmul(x)), tmul(y)
    return (a[0] - b[0], a[1] - b[1])

def ge(x, y):                     # exact order on Z[tau] inside R: x >= y ?
    p, q = x[0] - y[0], x[1] - y[1]     # p + q tau >= 0 ?  tau = (1+sqrt5)/2
    s, t = 2 * p + q, q                  # 2(p+q tau) = s + t sqrt5
    if s >= 0 and t >= 0: return True
    if s <= 0 and t <= 0: return s == 0 and t == 0
    return s * s >= 5 * t * t if s > 0 else 5 * t * t >= s * s

def star_in(x, lo, hi, lo_open=True):   # is x* = p + q(1 - tau) in (lo, hi] ?  lo, hi integers
    p, q = x
    xs = (p + q, -q)                    # p + q - q tau as element of Z[tau]
    a = ge(xs, (lo, 0)) and not (lo_open and xs == (lo, 0))
    return a and ge((hi, 0), xs)

def prod(u, v, S=None):
    r = defaultdict(Fr)
    for a, ca in u.items():
        for b, cb in v.items():
            for k in (qa(a, b), qa(b, a)):
                if S is None or k in S:
                    r[k] += ca * cb / 2
    return {k: c for k, c in r.items() if c}

def sub(u, v):
    r = defaultdict(Fr, u)
    for k, c in v.items(): r[k] -= c
    return {k: c for k, c in r.items() if c}

checks = []
def check(name, cond):
    checks.append((name, bool(cond)))
    print(("PASS " if cond else "FAIL ") + name)

T = lambda n: {0: (1, 0), 1: (0, 1), 2: (1, 1), 3: (1, 2), 4: (2, 3), 5: (3, 5), 6: (5, 8)}[n]
# tau^n = F_{n-1} + F_n tau
check("tau powers: tau^2..tau^6 = 1+tau, 1+2tau, 2+3tau, 3+5tau, 5+8tau",
      all(tmul(T(n - 1)) == T(n) for n in range(1, 7)))
check("tau^6 - (tau^4 + tau^2) = 4tau+2 > 0  and  tau^4 - 2tau^2 = tau > 0  and  tau^6-(tau^4+2tau^2)=3tau+1",
      (T(6)[0] - T(4)[0] - T(2)[0], T(6)[1] - T(4)[1] - T(2)[1]) == (2, 4)
      and (T(4)[0] - 2 * T(2)[0], T(4)[1] - 2 * T(2)[1]) == (0, 1)
      and (T(6)[0] - T(4)[0] - 2 * T(2)[0], T(6)[1] - T(4)[1] - 2 * T(2)[1]) == (1, 3))

# Quasiaddition identities quoted from the papers
pts = [(p, q) for p in range(-4, 5) for q in range(-4, 5)]
check("x|-x = x, x|-(x|-y) = y|-x, flexible, translation invariant",
      all(qa(x, x) == x for x in pts) and
      all(qa(x, qa(x, y)) == qa(y, x) and qa(x, qa(y, x)) == qa(qa(x, y), x) for x in pts[:30] for y in pts[:30]))
check("x|-y = x  iff  y = x  (and y|-x = x iff y = x)",
      all((qa(x, y) == x) == (x == y) and (qa(y, x) == x) == (x == y) for x in pts for y in pts))

# Theorem A, Fibonacci example x = L_1 + L_{1+tau}
x = {(1, 0): Fr(1), (1, 1): Fr(1)}
x2 = prod(x, x); x3 = prod(x2, x)
D = sub(prod(x2, x2), prod(x3, x))
check("x^2 = L_{-tau} + L_1 + L_{1+tau} + L_{2+2tau}",
      x2 == {(0, -1): 1, (1, 0): 1, (1, 1): 1, (2, 2): 1})
check("coefficient of L_{9+13tau} in x^2x^2 - x^3x is -1/4", D.get((9, 13)) == Fr(-1, 4))
x3x = prod(x3, x)
check("9+13tau = 1 + tau^6 * tau is the largest point of x^3 x, with coefficient 1/4",
      (9, 13) in x3x and x3x[(9, 13)] == Fr(1, 4) and all(ge((9, 13), z) for z in x3x)
      and (1 + T(6)[1], T(6)[0] + T(6)[1]) == (9, 13))
check("every point of x^2x^2 is <= 1 + (tau^4+tau^2) tau = 5+7tau < 9+13tau",
      all(ge((5, 7), z) for z in prod(x2, x2)) and not ge((5, 7), (9, 13)))
check("Jordan identity fails for (x, y) = (L_1+L_{1+tau}, L_1+L_{1+tau})",
      sub(prod(prod(x, x), x2), prod(x, prod(x, x2))) != {})

# Theorem A characterisation, random tests on Z[tau]
random.seed(1)
ok_fail = True; ok_hold = True
for _ in range(150):
    k = random.randint(2, 4)
    S = random.sample(pts, k)
    u = {s: Fr(random.choice([-3, -2, -1, 1, 2, 3])) for s in S}
    u2 = prod(u, u)
    if not sub(prod(u2, u2), prod(prod(u2, u), u)): ok_fail = False
for _ in range(60):
    p = random.choice(pts); c = Fr(random.choice([-2, -1, 1, 3]))
    u = {p: c}; u2 = prod(u, u)
    y = {s: Fr(random.randint(-3, 3)) for s in random.sample(pts, 3)}
    y = {a: b for a, b in y.items() if b}
    if sub(prod(prod(u, y), u2), prod(u, prod(y, u2))): ok_hold = False
check("random x with |supp| >= 2: x^2x^2 != x^3x (150 cases)", ok_fail)
check("x = cL_p: Jordan identity holds for random y (60 cases)", ok_hold)

# Points of the Fibonacci chain F_{1,0}: window (0,1]
fib = [(1, 0), (1, 1), (2, 2), (2, 3), (3, 4), (0, -1), (-1, -2), (-1, -3), (-2, -4)]
check("1, 1+tau, 2+2tau, 2+3tau, 3+4tau, -tau, -1-2tau, -1-3tau, -2-4tau lie in F_{1,0}",
      all(star_in(z, 0, 1) for z in fib))
check("0 and tau do not lie in F_{1,0}", not star_in((0, 0), 0, 1) and not star_in((0, 1), 0, 1))
check("1, 1+tau, 2+2tau are consecutive in F_{1,0}",
      all(not star_in((p, q), 0, 1) for p in range(-3, 6) for q in range(-3, 6)
          if ge((p, q), (1, 0)) and ge((2, 2), (p, q)) and (p, q) not in [(1, 0), (1, 1), (2, 2)]))

# Theorem B: three-point truncation
S3 = [(1, 0), (1, 1), (2, 2)]
e1, e2, e3 = ({z: Fr(1)} for z in S3)
tab = (prod(e1, e1, S3) == e1 and prod(e2, e2, S3) == e2 and prod(e3, e3, S3) == e3 and
       prod(e1, e2, S3) == {(2, 2): Fr(1, 2)} and prod(e1, e3, S3) == {} and prod(e2, e3, S3) == {})
check("truncation to {1, 1+tau, 2+2tau}: e_i^2 = e_i, e1e2 = e3/2, e1e3 = e2e3 = 0", tab)
xt = {(1, 0): Fr(1), (1, 1): Fr(1)}
xt2 = prod(xt, xt, S3)
check("truncation: x^2x^2 - x^3x = e3 for x = e1+e2",
      sub(prod(xt2, xt2, S3), prod(prod(xt2, xt, S3), xt, S3)) == {(2, 2): Fr(1)})

# Proposition C: the subspace L_1 o J is not an ideal of J(F_{1,0})
w = prod(prod({(1, 0): Fr(1)}, {(1, 1): Fr(1)}), {(2, 3): Fr(1)})
check("(L_1 o L_{1+tau}) o L_{2+3tau} = 1/4(L_{-4-7tau}+L_{6+9tau}+L_{1+tau}+L_{3+4tau})",
      w == {(-4, -7): Fr(1, 4), (6, 9): Fr(1, 4), (1, 1): Fr(1, 4), (3, 4): Fr(1, 4)})
check("1+tau = 1|-y only for y = 0, and = y|-1 only for y = tau (searched |p|,|q| <= 30)",
      [y for y in [(p, q) for p in range(-30, 31) for q in range(-30, 31)] if qa((1, 0), y) == (1, 1)] == [(0, 0)]
      and [y for y in [(p, q) for p in range(-30, 31) for q in range(-30, 31)] if qa(y, (1, 0)) == (1, 1)] == [(0, 1)])

# augmentation is multiplicative
u = {(0, 0): Fr(2), (1, 3): Fr(-1)}; v = {(2, 1): Fr(3), (1, 1): Fr(5)}
check("augmentation is multiplicative", sum(prod(u, v).values()) == sum(u.values()) * sum(v.values()))

n = sum(ok for _, ok in checks)
print(f"{n} checks passed, {len(checks) - n} failed")
