"""Level-0: Kalman/RTS numerics on synthetic data with known truth (TDS § 30.1: z-STD 0.9-1.1)."""
import copy

import numpy as np

from ppp import estimator as est
from ppp import timesys as ts
from tests import synth


def _z(sol, truth):
    z, s, _ = est.state_series(sol, "ZWD")
    ok = np.isfinite(z)
    return (z[ok] - truth["zwd"][ok]) / s[ok], np.sqrt(np.mean((z[ok] - truth["zwd"][ok]) ** 2))


def test_fixed_mode_synthetic_truth():
    obs, arcs, truth = synth.make_obs()
    sol = est.solve(obs, copy.deepcopy(arcs), est.EstConfig(mode="fixed"), max_passes=0)
    zs, rmse = _z(sol, truth)
    assert 0.9 <= np.std(zs) <= 1.1, np.std(zs)
    assert rmse < 0.006
    assert 0.8 < sol.nis < 1.2                               # NIS/dof criterion [TDS § 6.7]
    assert sol.chol_failures == 0


def test_static_mode_recovers_coordinates():
    dX = np.array([0.012, -0.020, 0.031])
    obs, arcs, truth = synth.make_obs(dX=dX, seed=5)
    sol = est.solve(obs, copy.deepcopy(arcs), est.EstConfig(mode="static"), max_passes=0)
    X, P = est.coordinate_estimate(sol, obs)
    err = X - obs.X0 - dX
    assert np.all(np.abs(err) < 4 * np.sqrt(np.diag(P)) + 1e-3)
    zs, rmse = _z(sol, truth)
    assert 0.85 <= np.std(zs) <= 1.15


def test_smoother_better_than_forward_and_outputs():
    obs, arcs, truth = synth.make_obs(seed=7)
    sol = est.solve(obs, copy.deepcopy(arcs), est.EstConfig(mode="fixed"), max_passes=1)
    z, s, zf = est.state_series(sol, "ZWD")
    ok = np.isfinite(z)
    assert np.sqrt(np.mean((z[ok] - truth["zwd"][ok]) ** 2)) < np.sqrt(np.mean((zf[ok] - truth["zwd"][ok]) ** 2))
    fm = est.extract_5min(obs, sol, 0, ts.DAY_NS, (), None, None, (False, False))
    assert len(fm.t_utc) == 288
    good = np.isfinite(fm.zwd)
    assert good.sum() > 250
    assert "EDGE" in fm.conv[:3] and fm.conv[100] == "CONVERGED"


def test_outlier_screening_and_gap_segments():
    obs, arcs, truth = synth.make_obs(seed=11)
    jp = int(np.nonzero(obs.usable[1000])[0][0])
    jc = int(np.nonzero(obs.usable[1500])[0][0])
    obs.Lif[1000, jp] += 0.5                                 # gross phase outlier
    obs.Pif[1500, jc] += 60.0                                # gross code outlier
    obs.usable[1200:1400] = False                            # 100-min station-wide gap -> single window
    sol = est.solve(obs, copy.deepcopy(arcs), est.EstConfig(mode="fixed"), max_passes=2)
    assert sol.rejected[1000, jp] or sol.edit_counts["phase_outliers_split"] >= 1
    assert sol.rejected[1500, jc] or sol.edit_counts["code_outliers_removed"] >= 1
    zs, rmse = _z(sol, truth)
    assert rmse < 0.008
    obs2, arcs2, _ = synth.make_obs(seed=11)
    obs2.usable[1200:1400] = False                           # > window_break_s (3600 s) -> 2 segments
    sol2 = est.solve(obs2, copy.deepcopy(arcs2), est.EstConfig(mode="fixed"), max_passes=0)
    assert len(sol2.segments) == 2
