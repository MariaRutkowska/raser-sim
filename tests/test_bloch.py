import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.bloch import simulate


def test_matches_analytic_solution():
    delta, T1, T2 = 2 * np.pi * 5.0, 3.0, 1.0
    t, (Mx, My, Mz) = simulate([1.0, 0.0, 0.0], 5.0, delta, T1, T2)
    assert np.allclose(Mx, np.cos(delta * t) * np.exp(-t / T2), atol=1e-6)
    assert np.allclose(Mz, 1 - np.exp(-t / T1), atol=1e-6)




def test_radiation_damping_analytic():
    tau = 0.02
    t, (Mx, My, Mz) = simulate([1.0, 0.0, 0.0], 0.2, 0.0, 1e9, 1e9, tau_rd=tau)
    assert np.allclose(Mx, 1 / np.cosh(t / tau), atol=1e-6)
    assert np.allclose(Mz, np.tanh(t / tau), atol=1e-6)



def test_raser_steady_state():
    T1, T2, tau, Mz_eq = 5.0, 1.0, 0.2, -1.0
    t, (Mx, My, Mz) = simulate([1e-3, 0.0, -1.0], 200.0, 0.0, T1, T2,
                               tau_rd=tau, Mz_eq=Mz_eq)
    assert abs(Mz[-1] - (-tau / T2)) < 1e-4
    assert abs(np.hypot(Mx[-1], My[-1]) - np.sqrt(tau * (-tau / T2 - Mz_eq) / T1)) < 1e-4





def test_one_line_raser_matches_bloch():
    from src.raser2 import raser_rhs
    from src.bloch import bloch_rhs
    y = np.array([0.1, -0.05, -0.7])
    a = raser_rhs(0, y, np.array([3.0]), 5.0, 1.0, 0.2, np.array([-1.0]))
    b = bloch_rhs(0, y, 3.0, 5.0, 1.0, tau_rd=0.2, Mz_eq=-1.0)
    assert np.allclose(a, b)