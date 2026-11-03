import numpy as np
from scipy.integrate import solve_ivp

def bloch_rhs(t, M, delta, T1, T2, M0=1.0, tau_rd=np.inf):
    """Right-hand side of the Bloch equations (no RF pulse).

    M = [Mx, My, Mz]. delta = offset from the rotating frame (rad/s).
    tau_rd = radiation-damping time (s); np.inf switches it off.
    """
    Mx, My, Mz = M
    k = 1.0 / (M0 * tau_rd)
    dMx = -delta * My - Mx / T2 - k * Mx * Mz
    dMy = delta * Mx - My / T2 - k * My * Mz
    dMz = (M0 - Mz) / T1 + k * (Mx**2 + My**2)
    return [dMx, dMy, dMz]


def simulate(M_init, t_end, delta, T1, T2, n_points=4000, tau_rd=np.inf):
    t = np.linspace(0, t_end, n_points)
    sol = solve_ivp(bloch_rhs, (0, t_end), M_init, t_eval=t, args=(delta, T1, T2, 1.0, tau_rd), rtol=1e-08, atol=1e-10)
    return sol.t, sol.y
