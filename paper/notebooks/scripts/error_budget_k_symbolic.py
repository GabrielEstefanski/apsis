"""
error_budget_k_symbolic.py
==========================

Symbolic (sympy) Lindstedt-Poincare derivation, to second order in 1/c^2, of
the apsidal advance of a test particle under the EXACT 1PN force implemented
in crates/apsis-1pn/src/lib.rs::accumulate_force (G = 1, planar polar):

    radial:     r'' - r*phi'^2 = -mu/r^2 + (mu/(c^2 r^2)) * (4*mu/r - v^2 + 4*r'^2)
    azimuthal:  (1/r) * d/dt(r^2 phi') = (mu/(c^2 r^2)) * 4 * r' * (r*phi')
    v^2 = r'^2 + r^2 phi'^2

Target:

    Delta_omega(eps) = 6*pi*eps * (1 + k(e)*eps + O(eps^2)),   eps = mu/(c^2*p)

Conventions (they matter at order eps^2 -- k(e) is convention-dependent):

  * h0 is the exact integral of motion found in gate G0, referenced at u = 0:
        h0 = h * exp(q*mu*u/c^2),   h = r^2*phi',  u = 1/r,  q determined by G0.
    p = h0^2/mu is the (Newtonian-limit) semi-latus rectum, eps = mu/(c^2*p).
  * e is the amplitude of the zeroth-order Lindstedt solution
        u0(s) = (mu/h0^2)*(1 + e*cos(s)),  s = Omega*phi,
    with the standard normalization that u1, u2 carry NO homogeneous cos(s)
    component (particular + bounded solutions only).
  * Delta_omega here is the advance per radial period of the s-variable.
    It coincides with the geometric periapsis-to-periapsis advance because
    periapses are minima of r = maxima of u, and u is exactly 2*pi-periodic
    in s, so successive periapses are separated by Delta(phi) = 2*pi/Omega.

Method note (deliberate): the derivation works DIRECTLY on the force
equations above. No Hamiltonian/Lagrangian is introduced -- a Lagrangian
matching this force at O(1/c^2) differs from it at O(1/c^4), which is
exactly the order computed here.

Hard gates:
  G0: the exact integral h*exp(q*mu*u/c^2) = const is verified symbolically
      under the azimuthal EOM, with q determined (not assumed) by sympy.
  G1: the first-order advance must equal 6*pi*eps EXACTLY (symbolic).
  G2: u1, u2 strictly periodic -- order-by-order residuals vanish
      identically with pure trig-polynomial solutions (no secular terms).
  G3: k(e) finite at e = 0.

No commits. Run: .venv/Scripts/python.exe paper/notebooks/scripts/error_budget_k_symbolic.py
"""

import sys

import sympy as sp


def fail(msg):
    print(f"BLOCKED: {msg}")
    sys.exit(1)


# ════════════════════════════════════════════════════════════════════════════
# G0 — exact integral of the azimuthal equation (time domain)
# ════════════════════════════════════════════════════════════════════════════
t = sp.Symbol("t", positive=True)
mu, c = sp.symbols("mu c", positive=True)
q = sp.Symbol("q")

r = sp.Function("r", positive=True)(t)
f = sp.Function("phi")(t)
rd, phd = r.diff(t), f.diff(t)

# Azimuthal EOM exactly as specified: (1/r) d/dt(r^2 phi') = (mu/(c^2 r^2))*4*r'*(r*phi')
azim_eom = sp.Eq(sp.diff(r**2 * phd, t) / r, (mu / (c**2 * r**2)) * 4 * rd * (r * phd))
phidd_from_eom = sp.solve(azim_eom, f.diff(t, 2))
if len(phidd_from_eom) != 1:
    fail(f"G0: azimuthal EOM not uniquely solvable for phi'': {phidd_from_eom}")
phidd_from_eom = phidd_from_eom[0]

h_t = r**2 * phd
u_t = 1 / r
candidate = h_t * sp.exp(q * mu * u_t / c**2)
dIdt = candidate.diff(t).subs(f.diff(t, 2), phidd_from_eom)
dIdt = sp.simplify(dIdt)

q_solutions = sp.solve(sp.Eq(dIdt, 0), q)
if not q_solutions:
    fail(f"G0: no constant q makes d/dt[h*exp(q*mu*u/c^2)] vanish. Residual: {dIdt}")
q_val = q_solutions[0]
if not q_val.is_constant():
    fail(f"G0: solved q is not a constant: {q_val}")
residual_at_q = sp.simplify(dIdt.subs(q, q_val))
if residual_at_q != 0:
    fail(f"G0: residual nonzero at q={q_val}: {residual_at_q}")
print(f"[G0] PASS  d/dt[ h*exp(q*mu*u/c^2) ] == 0 under the azimuthal EOM, q = {q_val}")
print(f"      => h0 = h*exp({q_val}*mu*u/c^2) exact invariant (reference u=0)")

# ════════════════════════════════════════════════════════════════════════════
# Step 1 — reduce to the orbit ODE u'' + u = F(u, u'; c) in phi, exact h(u)
# ════════════════════════════════════════════════════════════════════════════
ph = sp.Symbol("varphi")
d = sp.Symbol("delta", positive=True)  # delta = 1/c^2
h0, e = sp.symbols("h0 e", positive=True)

U = sp.Function("u", positive=True)(ph)
Up, Upp = U.diff(ph), U.diff(ph, 2)
hU = h0 * sp.exp(-q_val * mu * d * U)  # exact h(u) from G0 (delta = 1/c^2)

# Chain rule d/dt = phi' * d/dphi = h*u^2 * d/dphi, with r = 1/u:
r_x = 1 / U
phd_x = hU * U**2                       # phi'
rd_x = -hU * Up                         # r'  = -(1/u^2) du/dt = -h u'
rdd_x = hU * U**2 * sp.diff(rd_x, ph)   # r'' = (h u^2 d/dphi) r'
v2_x = rd_x**2 + r_x**2 * phd_x**2

# Consistency: substitutions must satisfy the azimuthal EOM identically.
azim_res = U * (hU * U**2) * sp.diff(hU, ph) - (mu * U**2 * d) * 4 * rd_x * (r_x * phd_x)
if sp.simplify(azim_res) != 0:
    fail(f"reduction violates azimuthal EOM: residual {sp.simplify(azim_res)}")
print("[chk] azimuthal EOM satisfied identically by the (h(u), u(phi)) substitutions")

radial_res = rdd_x - r_x * phd_x**2 - (
    -mu * U**2 + mu * U**2 * d * (4 * mu * U - v2_x + 4 * rd_x**2)
)
upp_solutions = sp.solve(radial_res, Upp)
if len(upp_solutions) != 1:
    fail(f"radial EOM not uniquely solvable for u'': {upp_solutions}")
F_exact = sp.simplify(upp_solutions[0] + U)  # u'' + u = F_exact(u, u'; delta), EXACT
print("[eq ] u'' + u =", F_exact)

# Series in delta = 1/c^2, KEEPING terms through delta^2 (i.e. 1/c^4):
F_ser = sp.expand(F_exact.series(d, 0, 3).removeO())
F2_coeff = sp.simplify(F_ser.coeff(d, 2))
print("[eq ] series to 1/c^4:  u'' + u =", sp.collect(F_ser, d))
print(f"      (note: the explicit delta^2 source is {F2_coeff} — the exp(8*mu*u*delta)"
      f" and (1-4*mu*u*delta) factors cancel at that order)")

# Move to plain symbols for substitution into the Lindstedt ansatz.
Uu, Uv = sp.symbols("Uu Uv")
F_sym = F_ser.subs({Up: Uv, U: Uu})

# ════════════════════════════════════════════════════════════════════════════
# Step 2 — Lindstedt-Poincare: s = Omega*phi, two orders
# ════════════════════════════════════════════════════════════════════════════
s = sp.Symbol("s")
w1, w2 = sp.symbols("w1 w2")
a0, a2, a3 = sp.symbols("a0 a2 a3")
b0, b2, b3, b4 = sp.symbols("b0 b2 b3 b4")

B = mu / h0**2  # = mu/p = 1/p in mu=1 units
Om = 1 + w1 * d + w2 * d**2

u0e = B * (1 + e * sp.cos(s))
# bounded ansatz: no homogeneous cos(s) component (e-normalization convention);
# RHS is even in s at every order (u0 even), so no sin terms — asserted below.
u1e = a0 + a2 * sp.cos(2 * s) + a3 * sp.cos(3 * s)
u2e = b0 + b2 * sp.cos(2 * s) + b3 * sp.cos(3 * s) + b4 * sp.cos(4 * s)
ue = u0e + d * u1e + d**2 * u2e

# d/dphi = Omega d/ds. F_sym is already the delta-series (polynomial in
# delta); the full residual is therefore polynomial in delta and plain
# expand + coeff extraction is exact — no series() call (its internal
# routing was observed to introduce spurious rational-in-cos terms).
res = Om**2 * ue.diff(s, 2) + ue - F_sym.subs({Uu: ue, Uv: Om * ue.diff(s)})
res = sp.expand(res)

E0 = sp.simplify(res.coeff(d, 0))
if E0 != 0:
    fail(f"order delta^0 not satisfied by u0 = (mu/h0^2)(1+e*cos s): {E0}")
print("[O:0] PASS  u0 = (mu/h0^2)*(1 + e*cos s) solves the zeroth order exactly")


def _laurent_coeffs(expr, n_guard):
    """Exact Fourier coefficients via the exponential basis: rewrite
    cos/sin as exp, substitute z = exp(I*s), and read Laurent
    coefficients with sp.Poly. No simplify in the pipeline."""
    z = sp.Symbol("z_exp")
    ex = sp.expand(sp.expand(expr.rewrite(sp.exp), power_exp=True))
    ex = ex.subs(sp.exp(sp.I * s), z)
    ex = sp.expand(ex.subs(sp.exp(-sp.I * s), 1 / z))
    for a in sp.preorder_traversal(ex):
        if a.has(s):
            fail(f"exp-basis reduction left s-dependence: {a}")
    poly = sp.Poly(sp.expand(sp.together(ex) * z**n_guard), z)
    coeffs = {}
    for (p,), coef in poly.terms():
        n = p - n_guard
        if abs(n) > n_guard:
            fail(f"harmonic |n|={abs(n)} exceeds guard {n_guard}")
        coeffs[n] = sp.expand(coef)
    return coeffs


def solve_order(E, unknowns, n_max, label):
    """Project order-residual E on harmonics in the exp basis, kill the
    e^{±is} secular pair via the frequency unknown, solve the bounded
    coefficients, and assert the residual dies identically."""
    n_guard = n_max + 3
    cf = _laurent_coeffs(E, n_guard)
    for n in range(0, n_guard + 1):
        plus, minus = cf.get(n, sp.Integer(0)), cf.get(-n, sp.Integer(0))
        if sp.expand(plus - minus) != 0:
            fail(f"{label}: n={n} exp-coefficients not even in s ({plus} vs {minus})")
    eqs = [cf.get(n, sp.Integer(0)) for n in range(0, n_guard + 1)]
    secular = eqs[1]
    if not secular.has(unknowns[0]):
        fail(f"{label}: e^(is) secular coefficient does not contain {unknowns[0]}: {secular}")
    sol = sp.solve(eqs, unknowns, dict=True)
    if len(sol) != 1:
        fail(f"{label}: harmonic system not uniquely solvable: {sol}")
    sol = sol[0]
    leftover_cf = _laurent_coeffs(sp.expand(E.subs(sol)), n_guard)
    leftover = {n: v for n, v in leftover_cf.items() if sp.expand(v) != 0}
    if leftover:
        fail(f"{label}: residual after solving all harmonics: {leftover}")
    return sol


# ── order delta (1/c^2) ─────────────────────────────────────────────────────
E1 = sp.expand(res.coeff(d, 1))
sol1 = solve_order(E1, [w1, a0, a2, a3], n_max=3, label="O(1/c^2)")
w1_v = sp.simplify(sol1[w1])
u1_v = sp.simplify(u1e.subs(sol1))
print(f"[O:1] PASS  w1 = {w1_v},  u1(s) = {u1_v}")

# ── order delta^2 (1/c^4) ───────────────────────────────────────────────────
E2 = sp.expand(res.coeff(d, 2).subs(sol1))
sol2 = solve_order(E2, [w2, b0, b2, b3, b4], n_max=4, label="O(1/c^4)")
w2_v = sp.simplify(sol2[w2])
u2_v = sp.simplify(u2e.subs(sol2))
print(f"[O:2] PASS  w2 = {w2_v},  u2(s) = {u2_v}")

# G2: solutions are finite trigonometric polynomials with constant
# coefficients — strictly periodic iff the exp-basis reduction succeeds
# (any secular s-polynomial term would survive as bare s and fail it).
for _name, val in [("w1", w1_v), ("u1", u1_v), ("w2", w2_v), ("u2", u2_v)]:
    _laurent_coeffs(sp.expand(val), 8)
print("[G2 ] PASS  u1, u2 strictly periodic (finite harmonic polynomials)")

# ════════════════════════════════════════════════════════════════════════════
# Step 3/4 — apsidal advance per radial period, eps-expansion, k(e)
# ════════════════════════════════════════════════════════════════════════════
inv_Om = (1 / (1 + w1_v * d + w2_v * d**2)).series(d, 0, 3).removeO()
d_omega = sp.expand(2 * sp.pi * (inv_Om - 1))  # = 2*pi*(-w1*d + (w1^2 - w2)*d^2 + ...)

eps_over_d = mu**2 / h0**2  # eps = mu/(c^2*p) = mu^2*delta/h0^2

c1 = sp.simplify(d_omega.coeff(d, 1))
g1_residual = sp.simplify(c1 - 6 * sp.pi * eps_over_d)
if g1_residual != 0:
    fail(f"G1: first-order advance is {c1}, not 6*pi*mu^2/h0^2; residual {g1_residual}")
print("[G1 ] PASS  Delta_omega_1 = 6*pi*eps exactly (eps = mu/(c^2*p), p = h0^2/mu)")

c2 = sp.simplify(d_omega.coeff(d, 2))
k = sp.simplify(sp.cancel(c2 / (6 * sp.pi * eps_over_d**2)))
if k.has(mu) or k.has(h0) or k.has(d):
    fail(f"k(e) failed to close in e alone: {k}")

k0 = sp.simplify(k.subs(e, 0))
if not k0.is_finite:
    fail(f"G3: k(0) not finite: {k0}")
print(f"[G3 ] PASS  k(0) = {k0} finite")

print()
print("RESULT")
print("  Delta_omega = 6*pi*eps*(1 + k(e)*eps + O(eps^2)),  eps = mu/(c^2*p),  p = h0^2/mu")
print(f"  k(e) = {k}")
print(f"  k(e) factored = {sp.factor(k)}")
print(f"  k(0.0) = {sp.N(k.subs(e, 0), 12)}")
print(f"  k(0.2) = {sp.N(k.subs(e, sp.Rational(1, 5)), 12)}")
print(f"  k(0.4) = {sp.N(k.subs(e, sp.Rational(2, 5)), 12)}")
