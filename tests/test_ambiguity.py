"""Stage-6 unit tests: LAMBDA, success rate, WL/NL resolution on synthetic data with known integers."""
import copy

import numpy as np

from ppp import ambiguity as amb
from ppp import estimator as est
from ppp import timesys as ts
from tests import synth


def test_lambda_finds_integer_ils():
    rng = np.random.default_rng(1)
    n = 6
    A = rng.normal(size=(n, n))
    Q = A @ A.T * 0.01 + np.eye(n) * 0.001
    true = rng.integers(-20, 20, n).astype(float)
    a = true + rng.multivariate_normal(np.zeros(n), Q)
    F, s, D = amb.lambda_ils(a, Q, 2)
    # brute force check of the best candidate in a small neighbourhood
    Qi = np.linalg.inv(Q)
    best = F[:, 0]
    fbest = (a - best) @ Qi @ (a - best)
    assert abs(fbest - s[0]) < 1e-6 * max(1, fbest)
    for _ in range(300):
        cand = np.round(a) + rng.integers(-1, 2, n)
        assert (a - cand) @ Qi @ (a - cand) >= fbest - 1e-9
    assert s[1] >= s[0]
    assert 0 < amb.bootstrap_success_rate(D) <= 1


def test_wl_nl_resolution_synthetic():
    obs, arcs, sel, truth = synth.make_ar()
    sol = est.solve(obs, copy.deepcopy(arcs), est.EstConfig(mode="fixed"), max_passes=0)
    extract = lambda s_: est.extract_5min(obs, s_, 0, ts.DAY_NS, (), None, None, (False, False))  # noqa: E731
    fm = extract(sol)
    res = amb.run_ar(obs, arcs, sol, sel, "fixed", fm, extract)
    m = res.metrics
    assert res.status in ("FIXED", "PARTIAL"), (res.status, m)
    assert m["n_WL_fixed"] >= 0.75 * m["n_WL_candidates"]
    assert m["fixing_ratio"] >= 0.8
    # every fixed SD must equal the true SD integer combination
    f1, f2 = synth.F1, synth.F2
    lam_nl = synth.C / (f1 + f2)
    k_wl = synth.C / (f1 - f2) * f2 / (f1 + f2)
    for a, b, val, sig in res.constraints:
        n1 = truth["ints"][a][0] - truth["ints"][b][0]
        nwl = (truth["ints"][a][0] - truth["ints"][a][1]) - (truth["ints"][b][0] - truth["ints"][b][1])
        assert abs(val - (lam_nl * n1 + k_wl * nwl)) < 1e-6
    # fixed ZWD not worse than float against truth
    zf, _, _ = est.state_series(sol, "ZWD")
    zx, _, _ = est.state_series(res.solution, "ZWD")
    ok = np.isfinite(zf) & np.isfinite(zx)
    rf = np.sqrt(np.mean((zf[ok] - truth["zwd"][ok]) ** 2))
    rx = np.sqrt(np.mean((zx[ok] - truth["zwd"][ok]) ** 2))
    assert rx <= rf * 1.1
