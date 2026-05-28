"""
Independent Plummer-softened Kepler integration to verify the
apsidal-precession closed form.

Integrates Sun + Mercury under U(r) = -GM/sqrt(r^2 + eps^2) with
RK45 (scipy), measures the cumulative drift of the periapsis angle
without modulo-2pi aliasing, and compares against the closed-form
prediction.
"""
import numpy as np
from scipy.integrate import solve_ivp

# Canonical units: G = M_sun = 1, AU.
GM = 1.0
a = 0.387098
e = 0.20563
eps = 0.02  # softening, AU
N_ORBITS = 50  # short run to keep wall-clock low

# Orbital period in canonical units.
T = 2.0 * np.pi * np.sqrt(a**3 / GM)

# Initial conditions: at periapsis along +x.
r_peri = a * (1 - e)
v_peri = np.sqrt(GM * (2.0 / r_peri - 1.0 / a))

# State: [x, y, vx, vy] of Mercury; Sun fixed at origin.
y0 = [r_peri, 0.0, 0.0, v_peri]

def rhs(t, y):
    x, yy, vx, vy = y
    r2 = x*x + yy*yy
    soft = (r2 + eps*eps)**1.5
    ax = -GM * x / soft
    ay = -GM * yy / soft
    return [vx, vy, ax, ay]

# Integrate, sampling at every orbit period.
t_eval = np.arange(0, N_ORBITS * T + T/2, T)
sol = solve_ivp(rhs, (0, t_eval[-1]), y0, t_eval=t_eval,
                rtol=1e-13, atol=1e-13, method='DOP853')

# Compute periapsis longitude at each sample.
# omega from eccentricity vector e_vec = v x L / GM - r_hat
omegas = []
for k in range(sol.y.shape[1]):
    x, yy, vx, vy = sol.y[:, k]
    r = np.sqrt(x*x + yy*yy)
    Lz = x*vy - yy*vx
    # Eccentricity vector components (planar):
    # e_vec = (v x L)/GM - r/|r|
    e_x = (vy * Lz) / GM - x / r
    e_y = (-vx * Lz) / GM - yy / r
    omega = np.arctan2(e_y, e_x)
    omegas.append(omega)

omegas = np.array(omegas)

# Unwrap to get the true cumulative drift.
omegas_unwrap = np.unwrap(omegas)
drift_per_orbit_unwrap = (omegas_unwrap[-1] - omegas_unwrap[0]) / N_ORBITS
drift_naive = ((omegas[-1] - omegas[0] + np.pi) % (2*np.pi)) - np.pi  # apsis-style wrap

# Closed-form prediction.
predicted = -3.0 * np.pi * eps**2 / (a**2 * (1 - e**2)**2)

print(f"=== Plummer Mercury, N={N_ORBITS} orbits, eps={eps} AU ===")
print(f"Per orbit:")
print(f"  Measured (true unwrap):   {drift_per_orbit_unwrap:.6e} rad/orbit")
print(f"  Predicted (closed-form):  {predicted:.6e} rad/orbit")
print(f"  Ratio meas/pred:          {drift_per_orbit_unwrap / predicted:.6f}")
print(f"")
print(f"Total over {N_ORBITS} orbits:")
print(f"  Measured cumulative:      {omegas_unwrap[-1] - omegas_unwrap[0]:.4f} rad")
print(f"  Predicted cumulative:     {predicted * N_ORBITS:.4f} rad")
print(f"")
print(f"Naive (apsis-style):")
print(f"  Naive drift end-init:     {drift_naive:.4f} rad (wraps to (-pi, pi])")
