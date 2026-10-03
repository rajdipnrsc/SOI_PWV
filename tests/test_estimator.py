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


def test_forward_screening_rejects_phase_outlier():
    obs, arcs, truth = synth.make_obs(seed=17)
    k = 1500
    # a satellite whose ambiguity is mature (arc started > 1 h earlier); for a newborn ambiguity an outlier is
    # indistinguishable from the ambiguity itself
    j = next(jj for jj in np.nonzero(obs.usable[k])[0] if arcs.arc_meta[int(arcs.arc[k, jj])]["first"] < k - 120
             and obs.el[k, jj] > np.radians(30))
    obs.Lif[k, j] += 0.30
    sol = est.solve(obs, copy.deepcopy(arcs), est.EstConfig(mode="fixed"), max_passes=0)
    assert sol.rejected[k, j]                                # caught by IGG-III in the forward pass
    assert sol.rejected.sum() < 0.002 * obs.usable.sum()      # and almost nothing else


def test_snr_variance_factor():
    """SNR model (TDS § 5.4): 10^((ref - SNR)/10) per frequency, IF-combined; elevation model where SNR is missing."""
    a1, a2 = 2.546, -1.546
    el = np.radians(np.array([[90.0, 30.0, 30.0]]))
    s1 = np.array([[45.0, 35.0, np.nan]])
    s2 = np.array([[45.0, 35.0, 40.0]])
    g = est.snr_var_factor(el, s1, s2, a1, a2, ref=45.0)
    assert g[0, 0] == np.float64(1.0)                       # reference SNR = zenith elevation model
    assert g[0, 1] == np.float64(10.0) or abs(g[0, 1] - 10.0) < 1e-12
    assert abs(g[0, 2] - 4.0) < 1e-12                       # 1/sin^2(30 deg)
    assert est._sigma2(el[0, 1], (1.0, 2.0), g[0, 1]) == 1.0 + 4.0 * g[0, 1]


def test_snr_weighting_runs_on_synthetic():
    """With an SNR model equivalent to the elevation model, the solution is unchanged (same weights)."""
    obs, arcs, truth = synth.make_obs(ne=720)
    sol0 = est.solve(copy.deepcopy(obs), copy.deepcopy(arcs), est.EstConfig(mode="fixed"), max_passes=0)
    # SNR chosen so that 10^((45 - S)/10) = 1/sin^2(e) on both frequencies
    s = 45.0 + 20.0 * np.log10(np.sin(np.maximum(obs.el, np.radians(1.0))))
    obs.var_b = est.snr_var_factor(obs.el, s, s, 2.546, -1.546)
    sol1 = est.solve(obs, copy.deepcopy(arcs), est.EstConfig(mode="fixed"), max_passes=0)
    z0, _, _ = est.state_series(sol0, "ZWD")
    z1, _, _ = est.state_series(sol1, "ZWD")
    assert np.nanmax(np.abs(z0 - z1)) < 1e-9
