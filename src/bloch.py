import numpy as np
from scipy.integrate import solve_ivp

def bloch_rhs(t, M, delta, T1, T2, M0=1.0, tau_rd=np.inf, Mz_eq=None):
    if Mz_eq is None:
        Mz_eq = M0
    Mx, My, Mz = M
    k = 1.0 / (M0 * tau_rd)
    dMx = -delta * My - Mx / T2 - k * Mx * Mz
    dMy = delta * Mx - My / T2 - k * My * Mz
    dMz = (Mz_eq - Mz) / T1 + k * (Mx**2 + My**2)
    return [dMx, dMy, dMz]


def simulate(M_init, t_end, delta, T1, T2, n_points=4000, tau_rd=np.inf, Mz_eq=None):
    t = np.linspace(0, t_end, n_points)
    sol = solve_ivp(bloch_rhs, (0, t_end), M_init, t_eval=t, args=(delta, T1, T2, 1.0, tau_rd, Mz_eq), rtol=1e-08, atol=1e-10)
    return sol.t, sol.y
