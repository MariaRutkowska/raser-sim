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