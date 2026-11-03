import numpy as np
from scipy.integrate import solve_ivp


def raser_rhs(t, y, deltas, T1, T2, tau_rd, Mz_eq):
    n = len(deltas)
    Mx, My, Mz = y[:n], y[n:2 * n], y[2 * n:]
    k = 1.0 / tau_rd
    Sx, Sy = Mx.sum(), My.sum()                 # what the coil sees
    dMx = -deltas * My - Mx / T2 - k * Mz * Sx
    dMy = deltas * Mx - My / T2 - k * Mz * Sy
    dMz = (Mz_eq - Mz) / T1 + k * (Mx * Sx + My * Sy)
    return np.concatenate([dMx, dMy, dMz])


def simulate_raser(deltas, t_end, T1, T2, tau_rd, Mz_eq, n_points=20000, tilt=1e-3, seed=0):
    deltas = np.asarray(deltas, dtype=float)
    n = len(deltas)
    rng = np.random.default_rng(seed)
    y0 = np.concatenate([tilt * rng.standard_normal(n), tilt * rng.standard_normal(n), -np.ones(n)])
    t = np.linspace(0, t_end, n_points)
    sol = solve_ivp(raser_rhs, (0, t_end), y0, t_eval=t, method="DOP853",
                    args=(deltas, T1, T2, tau_rd, np.asarray(Mz_eq, dtype=float)),
                    rtol=1e-9, atol=1e-11)
    Mx, My, Mz = sol.y[:n], sol.y[n:2 * n], sol.y[2 * n:]
    return sol.t, Mx, My, Mz