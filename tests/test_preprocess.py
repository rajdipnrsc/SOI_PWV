"""Level-0/3: slip detection (injected slips), false alarms, clock jumps, observable selection (TDS § 13, § 30.4)."""
import numpy as np

from ppp import preprocess as pp
from ppp import settings
from tests import synth

SLIPS = [(1, 0), (0, 1), (1, -1), (2, 1), (2, 2), (5, 5), (10, 9), (3, 0), (1, 1)]


def _detected(arcs, k, j):
    """Slip at epoch k detected if the arc id changes between k-1 and k (+/- 1 epoch tolerance)."""
    col = arcs.arc[:, j]
    for kk in (k - 1, k, k + 1):
        prev = col[:kk][col[:kk] >= 0]
        nxt = col[kk:][col[kk:] >= 0]
        if len(prev) and len(nxt) and prev[-1] != nxt[0]:
            return True
    return False


def run_injection(seed=3, spacing=80):
    """Inject every slip type repeatedly, separated by >= `spacing` epochs on each satellite, in one run."""
    rng = np.random.default_rng(seed)
    _, t, el, vis = synth.make_raw(seed=seed)
    ne, ns = el.shape
    slips = []
    kinds = SLIPS * 6
    rng.shuffle(kinds)
    for j in range(ns):
        k = 60
        while k < ne - 60 and kinds:
            if vis[k - 30:k + 30, j].all() and el[k, j] > np.radians(12):
                d1, d2 = kinds.pop()
                slips.append((k, j, d1, d2))
                k += spacing
            else:
                k += 7
    sel, t, el, vis = synth.make_raw(seed=seed, slips=slips)
    arcs = pp.detect_slips_and_arcs(sel, t, el, vis, np.radians(7.0))
    res = {s_: [0, 0] for s_ in SLIPS}
    for (k, j, d1, d2) in slips:
        res[(d1, d2)][1] += 1
        res[(d1, d2)][0] += _detected(arcs, k, j)
    return res


def test_injected_slips_detected():
    """All slips that change MW by >= 2 cycles or GF by >= 0.3 m are detected; single-cycle L1- or L2-only
    slips >= 80 %. (1,1)-type slips change neither MW nor GF materially and are not required (TDS § 30.4
    gate is defined for MW-detectable slips)."""
    tot = {"clear": [0, 0], "single": [0, 0]}
    for seed in (3, 4):
        res = run_injection(seed=seed)
        for (d1, d2), (hit, n) in res.items():
            dgf = abs(settings.C_LIGHT / synth.F1 * d1 - settings.C_LIGHT / synth.F2 * d2)
            if abs(d1 - d2) >= 2 or dgf >= 0.3:
                tot["clear"][0] += hit
                tot["clear"][1] += n
            elif (d1, d2) in ((1, 0), (0, 1)):
                tot["single"][0] += hit
                tot["single"][1] += n
    assert tot["clear"][0] == tot["clear"][1], tot
    assert tot["single"][0] >= 0.8 * tot["single"][1], tot


def test_no_false_alarms_quiet_ionosphere():
    sel, t, el, vis = synth.make_raw(seed=5)
    arcs = pp.detect_slips_and_arcs(sel, t, el, vis, np.radians(7.0))
    n_slip = arcs.reasons["MW"] + arcs.reasons["GF"]
    sat_days = vis.sum() * 30.0 / 86400.0
    assert n_slip / sat_days <= 1.0                     # <= 1 false alarm per satellite-day [TDS § 30.4]


def test_clock_jump_repair():
    sel, t, el, vis = synth.make_raw(seed=7)
    k0 = 1000
    sel.P1[k0:] += 1e-3 * settings.C_LIGHT               # 1 ms receiver clock jump on code only
    sel.P2[k0:] += 1e-3 * settings.C_LIGHT
    jumps = pp.detect_clock_jumps(sel, vis)
    assert len(jumps) == 1 and jumps[0][0] == k0 and jumps[0][2]
    Pif, Lif, _, _ = pp.if_combination(sel)
    d = Pif - Lif
    j = int(np.nonzero(vis[k0 - 1] & vis[k0])[0][0])
    assert abs(d[k0, j] - d[k0 - 1, j]) < 5.0            # consistency restored after repair


def test_observable_selection_stable():
    class O:
        pass
    o = O()
    o.epochs = np.arange(10, dtype=np.int64)
    o.sats = ["G01"]
    nan = np.nan
    o.obs = {"C1C": np.full((10, 1), 2e7), "L1C": np.full((10, 1), 1e8), "C2W": np.full((10, 1), 2e7),
             "L2W": np.array([[1e8]] * 9 + [[nan]]), "L2X": np.full((10, 1), 1e8), "C2X": np.full((10, 1), 2e7)}
    o.lli = {c: np.zeros((10, 1), dtype=np.int8) for c in o.obs}
    sel = pp.select_observables(o)
    assert set(sel.pair_names.values()) == {("C1C", "L1C", "C2W", "L2W")}   # no L2W<->L2X flipping
    assert np.isnan(sel.L2[9, 0])                                           # missing epoch becomes a gap
