"""TEST-ONLY: build SP3/CLK-like product objects from a RINEX navigation file.

Used to exercise the full processing chain in environments where CODE products cannot be downloaded.
Broadcast orbits (~1 m) and clocks (~1-2 ns) are far worse than CODE products, so ZTD from such a run is
NOT a product and is never written as one (the run is labelled PRODUCTS_OVERRIDE / TEST_BROADCAST).
This module is not imported by the program.
"""
import numpy as np

from ppp import antenna as an
from ppp import orbclk
from ppp import products as prd
from ppp import timesys as ts


def build(nav, t0, t1, atx=None, sp3_step_s=300, clk_step_s=30):
    """nav: dict prn -> [GpsEph]; returns a products.ProductSet with sp3 + clk (no ERP/OSB/ATT)."""
    ep = np.arange(t0 - (t0 % (sp3_step_s * ts.NS)) - 6 * sp3_step_s * ts.NS, t1 + 7 * sp3_step_s * ts.NS,
                   sp3_step_s * ts.NS, dtype=np.int64)
    sats = sorted(nav)
    pos = np.full((len(ep), len(sats), 3), np.nan)
    clk = np.full((len(ep), len(sats)), np.nan)
    for j, prn in enumerate(sats):
        for i, t in enumerate(ep):
            e = orbclk.select_eph(nav[prn], int(t), max_age_h=4.0)
            if e is None:
                continue
            r, dt = orbclk.broadcast_pos_clk(e, int(t))
            if atx is not None:
                ent = an.satellite_entry(atx, prn, int(t))
                p = an.pco_if(ent) if ent is not None else None
                if p is not None:           # broadcast refers to the APC: move to CoM along the nadir axis
                    r = r + p[2] * r / np.linalg.norm(r)
            pos[i, j] = r
            clk[i, j] = dt
    sp3 = orbclk.SP3(ep, sats, pos, clk, np.zeros((len(ep), len(sats)), dtype=bool), "IGS20", "BRD", "GPS",
                     float(sp3_step_s), "d", ["broadcast"])
    sp3.frames = {"IGS20"}
    cep = np.arange(t0 - (t0 % (clk_step_s * ts.NS)) - 10 * clk_step_s * ts.NS, t1 + 11 * clk_step_s * ts.NS,
                    clk_step_s * ts.NS, dtype=np.int64)
    csats = {}
    for prn in sats:
        vals = []
        for t in cep:
            e = orbclk.select_eph(nav[prn], int(t), max_age_h=4.0)
            vals.append(np.nan if e is None else orbclk.broadcast_pos_clk(e, int(t), include_rel=False)[1])
        v = np.array(vals)
        ok = np.isfinite(v)
        csats[prn] = (cep[ok], v[ok])
    clk_p = orbclk.ClockProduct(csats, 3.04, "BRD", "", "IGS20", "GPS", float(clk_step_s), ["broadcast"])
    ps = prd.ProductSet(status="OK_FLOAT_ONLY", reason="TEST_BROADCAST", family="TEST_BROADCAST", tier="TEST")
    ps.sp3, ps.clk = sp3, clk_p
    ps.frame_label = "IGS20"
    ps.files = [{"type": "BRDC", "filename": "broadcast (test only)", "sha256": None}]
    return ps
