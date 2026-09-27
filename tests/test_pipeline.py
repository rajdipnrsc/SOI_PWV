"""End-to-end run on the HYDE sample with TEST-ONLY broadcast-derived products (slow; --run-slow).

This exercises every stage of pwv_ppp.run() on real RINEX data. Broadcast orbits/clocks are ~100x worse than
CODE products, so only plausibility is checked, never accuracy.
"""
import os

import numpy as np
import pytest

from tests.conftest import SAMPLE

REF = os.environ.get("PWV_PPP_TEST_REF")        # folder with igs20*.atx and gpt3_1.grd (optional)


@pytest.mark.slow
@pytest.mark.skipif(not os.path.exists(os.path.join(SAMPLE, "HYDE015.24o")), reason="sample data absent")
@pytest.mark.skipif(not REF, reason="set PWV_PPP_TEST_REF to a folder with an igs20 ANTEX and gpt3_1.grd")
def test_end_to_end_broadcast(tmp_path):
    import glob

    import pwv_ppp
    from ppp import antenna as an
    from ppp import rinex
    from ppp import troposphere as tr
    from tests.tools import broadcast_products as bp
    obs = os.path.join(SAMPLE, "HYDE015.24o")
    nav = os.path.join(SAMPLE, "HYDE015.24n")
    atx = an.read_antex(sorted(glob.glob(os.path.join(REF, "igs20*.atx")))[-1])
    o = rinex.read_obs(obs, decimate_s=30)
    ps = bp.build(rinex.read_nav(nav), int(o.epochs[0]), int(o.epochs[-1]), atx)
    g = tr.read_gpt3(os.path.join(REF, "gpt3_1.grd"))
    args = pwv_ppp.parse_args([obs, nav, "--offline", "--out", str(tmp_path), "--quiet"])
    res = pwv_ppp.run(args, products_override=ps, tropo_override=([], None, g), atx_override=atx)
    fm = res["fm"]
    ztd = (fm.zhd0 + fm.zwd) * 1e3
    assert np.isfinite(ztd).sum() > 250
    assert 2200 < np.nanmedian(ztd) < 2700                       # physically plausible at 420 m
    for f in ("_ZTD.csv", "_PWV.csv", "_ZTD.nc", "_manifest.json", "_quicklook.png"):
        assert os.path.exists(os.path.join(str(tmp_path), "HYDE_2024015" + f))
