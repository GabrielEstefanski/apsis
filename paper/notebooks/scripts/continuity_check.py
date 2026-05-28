"""Independent truncated-Plummer counter-test reproduction."""
import numpy as np
from scipy.integrate import solve_ivp

G_eff = 1.0
a = 1.0
e = 0.5
R_c = 1.0
alpha = 0.8
T_end = 60.0
T_orbit = 2.0 * np.pi * np.sqrt(a**3 / G_eff)
dt_for_yoshida = 1e-3  # canonical time units; matches apsis test

def rhs(t, y):
    x, yy, vx, vy = y[0], y[1], y[2], y[3]
    r2 = float(x*x + yy*yy)
    if r2 < R_c * R_c:
        f = r2**(-1.5)
    else:
        f = alpha * r2**(-1.5)
    return [float(vx), float(vy), float(-G_eff*f*x), float(-G_eff*f*yy)]

r_peri = a * (1 - e)
v_peri = float(np.sqrt(G_eff * (1 + e) / (a * (1 - e))))
y0 = [float(r_peri), 0.0, 0.0, v_peri]

t_eval = np.linspace(0, T_end, 30000)
sol = solve_ivp(rhs, (0, T_end), y0, t_eval=t_eval,
                rtol=1e-12, atol=1e-13, method='DOP853')

x = sol.y[0]; yy = sol.y[1]; vx = sol.y[2]; vy = sol.y[3]
r = np.sqrt(x*x + yy*yy)
v = np.sqrt(vx*vx + vy*vy)

signs = r - R_c
crossings = []
for i in range(len(signs) - 1):
    if signs[i] * signs[i+1] < 0:
        frac = -signs[i] / (signs[i+1] - signs[i])
        tc = sol.t[i] + frac * (sol.t[i+1] - sol.t[i])
        vc = v[i] + frac * (v[i+1] - v[i])
        direction = "out" if signs[i] < 0 else "in"
        crossings.append((tc, vc, direction))

delta_F = 1.0 - alpha
E0 = 0.5 * v_peri**2 - G_eff / r_peri

print(f"=== Truncated-Plummer counter-test (scipy reproduction) ===")
print(f"R_c={R_c}, alpha={alpha}, dF={delta_F:.3f}")
print(f"T_orbit={T_orbit:.4f}, dt(Yoshida4)={dt_for_yoshida:.4e} canonical")
print(f"E0 specific = {E0:.6f}  (|E0| = {abs(E0):.6f})")
print(f"Crossings found over T=60: {len(crossings)}")
print()
print(f"{'#':>3} {'t_cross':>9} {'v_cross':>9} {'dir':>4} {'|dE|_bound':>13} {'/|E0|':>12}")
for i, (tc, vc, d) in enumerate(crossings):
    b = delta_F * vc * dt_for_yoshida
    print(f"  {i+1:>1} {tc:>9.4f} {vc:>9.6f} {d:>4} {b:>13.4e} {b/abs(E0):>12.4e}")

if crossings:
    vlist = [c[1] for c in crossings]
    print()
    print(f"=== Summary ===")
    print(f"v_cross range [{min(vlist):.4f}, {max(vlist):.4f}], mean {np.mean(vlist):.4f}")
    print(f"|dE|/|E0| bound: [{delta_F*min(vlist)*dt_for_yoshida/abs(E0):.4e}, "
          f"{delta_F*max(vlist)*dt_for_yoshida/abs(E0):.4e}]")
    print(f"Paper-reported measured: [4.7e-6, 2.0e-4]")
