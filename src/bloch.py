import numpy as np
from scipy.integrate import solve_ivp

def bloch_rhs(t, M, delta, T1, T2, M0=1.0):
    Mx, My, Mz = M
    dMx = -delta * My - Mx / T2
    dMy = delta * Mx - My / T2
    dMz = (M0 - Mz) / T1
    return [dMx, dMy, dMz]


def simulate(M_init, t_end, delta, T1, T2, n_points=4000):
    t = np.linspace(0, t_end, n_points)
    sol = solve_ivp(bloch_rhs, (0, t_end), M_init, t_eval=t, args=(delta, T1, T2), rtol=1e-08, atol=1e-10)
    return sol.t, sol.y
