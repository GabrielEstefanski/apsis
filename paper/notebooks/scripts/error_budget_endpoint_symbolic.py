"""
error_budget_endpoint_symbolic.py
=================================

The O(eps) structure of the *osculating* argument of periapsis omega
along the truncated-1PN flow, in the gate's measurement
convention (Newtonian e-vector, mu = 1).

Motivation. The Phase-B ensembles show a constant-in-N angle offset of
~ -1.5e-8 rad between the integrated advance and the first-order secular
formula. `System::integrate_until` exits at the first accepted step with
t >= t_end (overshoot up to one adaptive IAS15 step), so the endpoint
state is sampled at a small true anomaly nu past the N-th periapsis.
This script derives, blind (from the osculating definitions and the
equation of motion only — no perturbation-theory formulas recalled),
the closed form of d(omega_osc)/d(phi) at O(delta), its secular part,
its periodic part, and the endpoint-offset function

    Q(nu) = omega_osc(nu past periapsis) - omega_osc(periapsis)

which converts a measured endpoint true anomaly into a predicted angle
residual. Every claim is asserted symbolically (gates GB0..GB5).

Derivation route (definitional, no Gauss equations recalled)
------------------------------------------------------------
Planar motion, u = 1/r, ' = d/dphi, h = r^2 phidot. The osculating
(Newtonian, mu) elements of the instantaneous state satisfy

    e_osc cos(nu) = h^2 u / mu - 1  =: Aexpr
    e_osc sin(nu) = -h^2 u' / mu    =: Bexpr        (rdot = -h u')
    omega_osc     = phi - nu .

Differentiating along the flow and eliminating with the perturbed Binet
equation u'' + u = mu/h^2 + W and the angular-momentum equation
h' = S/(h u^3) gives an expression for omega' that is *identically zero*
for any Kepler arc (gate GB0) and is therefore purely O(delta): the
perturbing force enters explicitly, the O(delta) deformation of the
trajectory does not.

Force (crates/apsis-1pn, test-particle limit, delta = 1/c^2):

    a_pert = (mu delta u^2) [ (4 mu u - v^2) n_hat + 4 (n_hat . v) v ]

with components R (radial, outward) and S (transverse, along-track):

    R = mu delta u^2 (4 mu u - v^2 + 4 rdot^2)
    S = 4 mu delta u^2 rdot (h u)

Run:  python paper/notebooks/scripts/error_budget_endpoint_symbolic.py
"""

import sympy as sp

# ── Symbols ───────────────────────────────────────────────────────────────────

mu, delta, hp, e, nu, phi = sp.symbols("mu delta h_p e nu phi", positive=True)
# hp: angular momentum at the periapsis IC (the gate's IC fixes it:
# h_p^2 = mu * p with p = a (1 - e^2)).

# ── Zeroth-order Kepler arc through the periapsis IC ─────────────────────────
# omega0 = 0 (IC on +x axis), so nu = phi at leading order.

w = 1 + e * sp.cos(phi)                  # 1 + e cos nu
u0 = (mu / hp**2) * w
u0p = sp.diff(u0, phi)
v2_0 = hp**2 * (u0p**2 + u0**2)          # v^2 = h^2 (u'^2 + u^2)
rdot0 = -hp * u0p

# ── GB-h: angular momentum along the flow ─────────────────────────────────────
# dh/dphi = hdot / phidot = (r S) / (h u^2) = S / (h u^3).
# With S = 4 mu delta u^2 rdot (h u) and rdot = -h u':
#   dh/dphi = -4 mu delta h u'  ==>  h = h_p exp(-4 mu delta (u - u_p)).
# (Exact integral; A1's G0 gate found the same exponent independently.)

S0 = 4 * mu * delta * u0**2 * rdot0 * (hp * u0)
dh_dphi_flow = S0 / (hp * u0**3)
gb_h = sp.simplify(dh_dphi_flow - (-4 * mu * delta * hp * u0p))
assert gb_h == 0, f"GB-h FAIL: dh/dphi identity residual {gb_h}"
print("[GB-h] PASS  dh/dphi = -4 mu delta h u'  (h = h_p exp(-4 mu delta (u - u_p)))")

# ── omega' from the osculating definitions, exact in the perturbation ────────
# With Aexpr = h^2 u/mu - 1, Bexpr = -h^2 u'/mu, nu = atan2(B, A):
#   nu' = (B' A - A' B) / (A^2 + B^2)
# Substituting u'' = mu/h^2 - u + W and h' (flow):
#   B' = A - (h^2/mu) W - 2 h h' u'/mu ,   A' = -B + 2 h h' u/mu
#   omega' = 1 - nu'
#          = [ (h^2 W/mu + 2 h h' u'/mu) A + (2 h h' u/mu) B ] / e_osc^2
# Kepler limit: W = 0, h' = 0  ==>  omega' = 0 identically (GB0).
# We now *re-derive* that algebra with sympy from scratch (no shortcut),
# using generic functions u(phi), h(phi) and an unforced Binet residual.

uf = sp.Function("u")(phi)
hf = sp.Function("h")(phi)
Wf = sp.Function("W")(phi)      # Binet residual: u'' + u - mu/h^2
A_expr = hf**2 * uf / mu - 1
B_expr = -(hf**2) * sp.diff(uf, phi) / mu

nup = (sp.diff(B_expr, phi) * A_expr - sp.diff(A_expr, phi) * B_expr) / (
    A_expr**2 + B_expr**2
)
omega_prime_exact = sp.simplify(
    (1 - nup)
    .subs(sp.diff(uf, phi, 2), mu / hf**2 - uf + Wf)
    .doit()
)

# GB0: Kepler arc — W = 0, h' = const — omega' vanishes identically.
omega_prime_kepler = omega_prime_exact.subs(Wf, 0).subs(sp.Derivative(hf, phi), 0)
gb0 = sp.simplify(omega_prime_kepler)
assert gb0 == 0, f"GB0 FAIL: Kepler omega' = {gb0}"
print("[GB0] PASS  omega' == 0 on any Kepler arc (independent of the arc)")

# ── O(delta): substitute force terms and the zeroth-order arc ────────────────
# W = -R/(h^2 u^2) - S u'/(h^2 u^3); both evaluated on the Kepler arc
# (omega' is already O(delta), so zeroth-order substitution is exact at O(delta)).

R0 = mu * delta * u0**2 * (4 * mu * u0 - v2_0 + 4 * rdot0**2)
W0 = -R0 / (hp**2 * u0**2) - S0 * u0p / (hp**2 * u0**3)
hp_flow = -4 * mu * delta * hp * u0p     # h' on the arc (GB-h)

omega_prime = omega_prime_exact.subs(
    [(Wf, W0), (sp.Derivative(hf, phi), hp_flow), (sp.Derivative(uf, phi), u0p), (uf, u0), (hf, hp)]
)
omega_prime = sp.simplify(sp.expand_trig(sp.expand(omega_prime)))

# Strip O(delta^2): omega' is a polynomial in delta with zero constant term.
omega_prime_1 = sp.expand(omega_prime).coeff(delta, 1) * delta
eps_t = mu**2 * delta / hp**2            # eps-tilde = mu delta / p,  p = hp^2/mu

# GB1: closed form of omega'/eps at O(delta).
closed = -(3 / e) * sp.cos(phi) + 3 - 5 * sp.cos(2 * phi) + e * sp.cos(phi)
gb1 = sp.simplify(omega_prime_1 / eps_t - closed)
assert gb1 == 0, f"GB1 FAIL: omega'/eps - closed = {gb1}"
print("[GB1] PASS  omega'/eps = -(3/e) cos nu + 3 - 5 cos 2nu + e cos nu")

# GB2: secular part = 3 eps per radian  ==>  6 pi eps per orbit.
secular = sp.integrate(omega_prime_1, (phi, 0, 2 * sp.pi)) / (2 * sp.pi)
gb2 = sp.simplify(secular - 3 * eps_t)
assert gb2 == 0, f"GB2 FAIL: secular - 3 eps = {gb2}"
print("[GB2] PASS  <omega'> = 3 eps  ==>  Delta_omega = 6 pi eps / orbit "
      "(cross-checks A1 gate G1)")

# GB3: periodic part, odd in nu (antisymmetric about periapsis).
P_nu = sp.integrate(omega_prime_1 - secular, (phi, 0, nu))
P_closed = eps_t * (-(3 / e - e) * sp.sin(nu) - sp.Rational(5, 2) * sp.sin(2 * nu))
gb3 = sp.simplify(P_nu - P_closed)
assert gb3 == 0, f"GB3 FAIL: P(nu) - closed = {gb3}"
gb3_odd = sp.simplify(P_closed + P_closed.subs(nu, -nu))
assert gb3_odd == 0, f"GB3 FAIL: P not odd: {gb3_odd}"
print("[GB3] PASS  P(nu) = eps [ -(3/e - e) sin nu - (5/2) sin 2nu ]   (odd)")

# ── Endpoint-offset function Q(nu) ────────────────────────────────────────────
# omega_osc at successive periapsis passages increments by exactly the
# secular advance (P(0) = 0), so for an endpoint sampled at small true
# anomaly nu past the N-th periapsis the residual vs 6 pi eps N is
#   Q(nu) = 3 eps nu + P(nu).

Q_nu = 3 * eps_t * nu + P_closed

# GB4: slope at periapsis, with closed-form factorisation.
slope = sp.diff(Q_nu, nu).subs(nu, 0)
slope_closed = -eps_t * (3 - e) * (1 + e) / e
gb4 = sp.simplify(slope - slope_closed)
assert gb4 == 0, f"GB4 FAIL: Q'(0) - closed = {gb4}"
print("[GB4] PASS  Q'(0) = -eps (3 - e)(1 + e) / e")

# GB5: time chain — instantaneous rate at periapsis (for the overshoot
# story; the per-run H2 correction uses Q(nu_end) directly and does not
# need this mapping).  phidot(peri) = h_p u_p^2, u_p = (mu/h_p^2)(1+e).
phidot_peri = hp * ((mu / hp**2) * (1 + e)) ** 2
domega_dt_peri = sp.simplify(slope_closed * phidot_peri)
dodt_closed = -(mu**4 * delta / hp**5) * (3 - e) * (1 + e) ** 3 / e
gb5 = sp.simplify(domega_dt_peri - dodt_closed)
assert gb5 == 0, f"GB5 FAIL: domega/dt(peri) - closed = {gb5}"
print("[GB5] PASS  domega/dt|peri = -(mu^4 delta / h_p^5) (3 - e)(1 + e)^3 / e")

# ── Numeric evaluation for the gate scenario (Mercury, both conventions) ─────

print()
print("=" * 72)
print("MACHINE-READABLE SUMMARY")
print("=" * 72)
print("omega'/eps      = -(3/e) cos nu + 3 - 5 cos 2nu + e cos nu")
print("P(nu)/eps       = -(3/e - e) sin nu - (5/2) sin 2nu")
print("Q(nu)/eps       = 3 nu - (3/e - e) sin nu - (5/2) sin 2nu")
print("Q'(0)/eps       = -(3 - e)(1 + e)/e")
print()

A_lit = sp.Float("0.387098", 15)
E_lit = sp.Float("0.20563", 15)
for label, c_val in [("raw_c   c=10065.3201686", sp.Float("10065.3201686", 15))]:
    # The numeric c only scales eps; the harness reads eps per-row from the
    # CSV (eps = predicted/(6 pi N)), so only the e-dependent factors below
    # are consumed downstream. Printed for the notebook table.
    eps_val = 1 / (c_val**2 * A_lit * (1 - E_lit**2))
    slope_over_eps = -((3 - E_lit) * (1 + E_lit)) / E_lit
    print(f"{label}:")
    print(f"  eps(Mercury)        = {float(eps_val):.6e}")
    print(f"  Q'(0)/eps           = {float(slope_over_eps):.6f}")
    print(f"  Q'(0)               = {float(slope_over_eps * eps_val):.6e} rad/rad")
print()
print("All gates GB-h, GB0..GB5: PASS")
