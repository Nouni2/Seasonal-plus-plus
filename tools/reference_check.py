"""Correctness of Spp against independent references (CONTRIBUTING.md, Change Verification).

Every check compares Spp, through spp_dump, with a reference that shares no code
with it, and passes when the difference is within the documented model
difference (wiki: Accuracy, Error budget). The program prints one line per check and PASS
or FAIL, writes every measurement to a JSON file, and exits with 0 or 1.

References
  IERS EOP 20 C04 and Leap_Second.dat       UT1 - UTC, TAI - UTC     (astropy-iers-data)
  ERFA, the IAU SOFA routines                time, sidereal time, nutation, precession (pyerfa)
  JPL DE421                                  Sun, 1900-2050           (skyfield, skyfield-data)
  JPL DE431, Swiss Ephemeris compression     Sun, 1000 BC-3000 AD     (pyswisseph + files)
  Swiss Ephemeris calendar routines          Julian/Gregorian calendar

Setup (once):
  python3 -m venv .venv
  .venv/bin/pip install skyfield skyfield-data pyerfa pyswisseph numpy astropy-iers-data
Run:
  .venv/bin/python -I tools/reference_check.py build/spp_dump [--ephe DIR] [--json FILE] [--only a,b]

The Swiss Ephemeris files (sepl*/semo*, 1200 BC-3600 AD) are downloaded into
--ephe (default: reference/ephe) on first use.
"""

import argparse
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

import erfa
import numpy as np
import swisseph as swe
from skyfield.api import load, load_file, wgs84
from skyfield.framelib import ecliptic_frame

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import hour_angle_cio as hc  # noqa: E402

import astropy_iers_data  # noqa: E402
import skyfield_data  # noqa: E402

IERS_DIR = Path(astropy_iers_data.__file__).parent / "data"
DE421 = Path(skyfield_data.__file__).parent / "data" / "de421.bsp"
EPHE_FILES = ["seplm12", "seplm06", "sepl_00", "sepl_06", "sepl_12", "sepl_18", "sepl_24", "sepl_30",
              "semom12", "semom06", "semo_00", "semo_06", "semo_12", "semo_18", "semo_24", "semo_30"]
EPHE_URL = "https://raw.githubusercontent.com/aloistr/swisseph/master/ephe/{}.se1"

ARCSEC = 3600.0
DAY_S = 86400.0
TT_TAI = 32.184

OBSERVERS = {
    "Paris": (48.8566, 2.3522, 35.0),
    "Quito": (-0.1807, -78.4678, 2850.0),
    "Singapore": (1.3521, 103.8198, 15.0),
    "Sydney": (-33.8688, 151.2093, 58.0),
    "Reykjavik": (64.1466, -21.9426, 20.0),
    "Tromso": (69.6492, 18.9553, 10.0),
    "Longyearbyen": (78.2232, 15.6267, 10.0),
    "McMurdo": (-77.8460, 166.6760, 10.0),
}
# Sites where the Sun can graze the condition, so that a small altitude difference
# becomes a large time difference: the horizon above about 66°, -6° above about 60°.
POLAR = {"Tromso", "Longyearbyen", "McMurdo"}
GRAZING_TWILIGHT = POLAR | {"Reykjavik"}

EXE = None
RESULTS = {}
CHECKS = []


# ----------------------------------------------------------------------------
# Plumbing

def check(name, value, limit, unit, exact=False):
    ok = (value == limit) if exact else (value <= limit)
    CHECKS.append({"check": name, "value": value, "limit": limit, "unit": unit, "pass": bool(ok)})
    rel = "==" if exact else "<="
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: {value:.6g} {unit} ({rel} {limit:g})", flush=True)


def run_spp(lines):
    out = subprocess.run([str(EXE)], input="\n".join(lines) + "\n", capture_output=True, text=True, check=True)
    rows = [r.split() for r in out.stdout.splitlines()]
    assert len(rows) == len(lines), (len(rows), len(lines))
    return rows


def fnum(x):
    return "nan" if x is None or not np.isfinite(x) else repr(float(x))


POS_COLS = ["tai_utc", "tt_utc", "dut1", "delta_t", "tt_day", "tt_frac", "ut1_day", "ut1_frac", "gmst", "gast",
            "dpsi", "deps", "eps0", "eps", "geo_lon", "geo_lat", "geo_r", "app_lon", "ra", "dec",
            "alt_geom", "az", "alt_refr", "sd"]


def spp_pos(jd_day, frac, dut1, obs):
    rows = run_spp([f"pos {int(d)} {fnum(f)} {fnum(u)} {obs[0]} {obs[1]} {obs[2]}" for d, f, u in zip(jd_day, frac, dut1)])
    data = np.array([[float(x) for x in r[1:]] for r in rows])
    return {c: data[:, i] for i, c in enumerate(POS_COLS)}


def spp_tt(jd_day, frac):
    rows = run_spp([f"tt {int(d)} {fnum(f)}" for d, f in zip(jd_day, frac)])
    return np.array([[float(x) for x in r[1:]] for r in rows])


def angdiff_deg(a, b):
    return (np.asarray(a) - np.asarray(b) + 180.0) % 360.0 - 180.0


def stats(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {"n": 0, "max_abs": 0.0}
    return {"n": int(x.size), "max_abs": float(np.max(np.abs(x))), "rms": float(np.sqrt(np.mean(x * x))),
            "mean": float(np.mean(x)), "min": float(np.min(x)), "max": float(np.max(x))}


def load_c04():
    rows = []
    for line in (IERS_DIR / "eopc04.1962-now").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        rows.append((int(p[0]), int(p[1]), int(p[2]), float(p[4]), float(p[7])))
    a = np.array(rows)
    return {"y": a[:, 0].astype(int), "m": a[:, 1].astype(int), "d": a[:, 2].astype(int), "mjd": a[:, 3],
            "dut1": a[:, 4]}


def load_leap_seconds():
    steps = []
    for line in (IERS_DIR / "Leap_Second.dat").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        steps.append((int(float(p[0])), int(p[4])))
    return np.array([m for m, _ in steps]), np.array([v for _, v in steps])


def ensure_ephe(path):
    path.mkdir(parents=True, exist_ok=True)
    for name in EPHE_FILES:
        f = path / f"{name}.se1"
        if not f.exists():
            print(f"  downloading {f.name}", flush=True)
            urllib.request.urlretrieve(EPHE_URL.format(name), f)
    swe.set_ephe_path(str(path))


C04 = None
LEAP_MJD = LEAP_VAL = None


def iers_dut1_at(mjd):
    return np.interp(mjd, C04["mjd"], C04["dut1"])


# ----------------------------------------------------------------------------
# Time scales

def section_time():
    mjd = C04["mjd"]
    s = spp_pos(mjd.astype(np.int64) + 2400000, np.full(mjd.shape, 0.5), [None] * mjd.size, OBSERVERS["Paris"])
    ref_tai = np.array([erfa.dat(*erfa.jd2cal(2400000.5, m)) for m in mjd])
    post72 = mjd >= 41317
    mism = int(np.sum(s["tai_utc"][post72] != ref_tai[post72]))
    future = np.arange(int(mjd[-1]) + 1, 61600)
    sf = spp_pos(future + 2400000, np.full(future.shape, 0.5), [None] * future.size, OBSERVERS["Paris"])
    ref_f = LEAP_VAL[np.searchsorted(LEAP_MJD, future, side="right") - 1]
    known = future < 61584
    mism += int(np.sum(sf["tai_utc"][known] != ref_f[known]))
    check("TAI - UTC vs IERS, every day 1972-01-01..2027-06-27 (mismatched days)", mism, 0, "days", exact=True)
    check("TAI - UTC unset after the IERS file expiry (days set)", int(np.sum(np.isfinite(sf["tai_utc"][~known]))),
          0, "days", exact=True)

    err = s["dut1"] - C04["dut1"]
    sel = post72
    RESULTS["dut1_vs_iers_s"] = stats(err[sel])
    check("estimated DUT1 vs IERS C04, every day from 1972", float(np.max(np.abs(err[sel]))), 0.010, "s")
    first = (C04["d"] == 1) & post72
    check("estimated DUT1 vs IERS C04, 1st of each month from 1972", float(np.max(np.abs(err[first]))), 1e-9, "s")

    # Across each leap second (estimated DUT1): ΔT = TT - UT1 continuous, and UT1
    # advancing by the inserted second plus the UTC interval (1 ms here).
    lines = []
    for m in LEAP_MJD[1:]:
        lines += [f"pos {int(m) + 2400000} {fnum(0.5 - 0.001 / DAY_S)} nan 0 0 0", f"pos {int(m) + 2400000} 0.5 nan 0 0 0"]
    rows = run_spp(lines)
    dt_jump, ut1_err = 0.0, 0.0
    for i in range(0, len(rows), 2):
        b0, b1 = rows[i], rows[i + 1]
        dt_jump = max(dt_jump, abs(float(b1[4]) - float(b0[4])))
        ut1 = ((float(b1[7]) - float(b0[7])) + (float(b1[8]) - float(b0[8]))) * DAY_S
        ut1_err = max(ut1_err, abs(ut1 - 1.001))
    check(f"ΔT jump across the {len(LEAP_MJD) - 1} leap seconds (estimated DUT1)", dt_jump, 1e-6, "s")
    check("UT1 advance across a leap second minus 1.001 s", ut1_err, 1e-6, "s")
    check("largest |estimated DUT1| from 1972 (the clamp is at 0.9 s)",
          float(np.max(np.abs(np.concatenate([s["dut1"][post72], sf["dut1"]])))), 0.9 - 1e-9, "s")


def section_delta_t_pieces():
    # Jumps of the ΔT model where one piece hands over to the next, and TT - UTC
    # across the end of the leap-second table.
    table_first = 2000.0 + (41317 - 51544.5) / 365.25          # 1972-01-01
    table_last = 2000.0 + (int(C04["mjd"][C04["d"] == 1][-1]) - 51544.5) / 365.25
    bounds = [-500.0, 500.0, 1600.0, 1700.0, 1800.0, 1860.0, 1900.0, 1920.0, 1941.0, 1961.0, table_first,
              table_last, 2050.0, 2150.0]
    lines = []
    for b in bounds:
        lines += [f"dt {fnum(np.nextafter(b, -np.inf))}", f"dt {fnum(b)}"]
    v = [float(r[1]) for r in run_spp(lines)]
    jumps = {f"{b:.4f}": round(v[2 * i + 1] - v[2 * i], 4) for i, b in enumerate(bounds)}
    RESULTS["delta_t_jumps_s"] = jumps
    for b, j in jumps.items():
        print(f"  info  ΔT jump at {b}: {j:+.4f} s", flush=True)
    rows = run_spp(["pos 2461584 0.49999999999 nan 0 0 0", "pos 2461584 0.5 nan 0 0 0"])
    step = abs(float(rows[1][2]) - float(rows[0][2]))
    check("TT - UTC step at the end of the leap-second table (2027-06-28 0h UTC)", step, 1e-9, "s")


def section_calendar():
    rng = np.random.default_rng(5)
    lines, ref = [], []
    for jd in rng.uniform(1355818.5, 2816787.5, 20000):
        cal = swe.GREG_CAL if jd >= 2299160.5 else swe.JUL_CAL
        y, m, d, _ = swe.revjul(float(np.floor(jd - 0.5) + 0.5), cal)
        h, mi = int(rng.integers(0, 24)), int(rng.integers(0, 60))
        sec = min(float(np.round(rng.uniform(0, 60), 3)), 59.999)
        lines.append(f"cal {y} {m} {d} {h} {mi} {sec!r}")
        ref.append((y, m, d, h, mi, sec, swe.julday(y, m, d, h + mi / 60.0 + sec / 3600.0, cal)))
    rows = run_spp(lines)
    jd_err, bad = [], 0
    for r, (y, m, d, h, mi, sec, jd_ref) in zip(rows, ref):
        jd_err.append((float(r[1]) + float(r[2]) - jd_ref) * DAY_S)
        if tuple(int(x) for x in r[3:8]) != (y, m, d, h, mi) or abs(float(r[8]) - sec) > 1e-6:
            bad += 1
    check("civil -> JD vs Swiss Ephemeris julday, 1000 BC-3000 AD", float(np.max(np.abs(jd_err))), 1e-4, "s")
    check("JD -> civil round trip mismatches", bad, 0, "dates", exact=True)


# ----------------------------------------------------------------------------
# Earth orientation

def section_earth():
    rng = np.random.default_rng(1)
    jd = np.concatenate([2451545.0 + rng.uniform(-36525 * 30, 36525 * 10, 4000),
                         2451545.0 + rng.uniform(-36525, 36525, 4000)])
    d = np.floor(jd).astype(np.int64)
    s = spp_pos(d, jd - d, [0.0] * jd.size, OBSERVERS["Paris"])
    tt1, tt2, ut1, ut2 = s["tt_day"], s["tt_frac"], s["ut1_day"], s["ut1_frac"]
    modern = np.abs((jd - 2451545.0) / 365.25) <= 100.0
    mas = ARCSEC * 1e3
    dpsi80, deps80 = erfa.nut80(tt1, tt2)
    dpsi06, deps06 = erfa.nut06a(tt1, tt2)
    g82 = angdiff_deg(s["gmst"], np.degrees(erfa.gmst82(ut1, ut2))) * mas
    g94 = angdiff_deg(s["gast"], np.degrees(erfa.gst94(ut1, ut2))) * mas
    n80 = np.maximum(np.abs(s["dpsi"] - np.degrees(dpsi80)), np.abs(s["deps"] - np.degrees(deps80))) * mas
    n06 = np.maximum(np.abs(s["dpsi"] - np.degrees(dpsi06)), np.abs(s["deps"] - np.degrees(deps06))) * mas
    check("IAU 1980 nutation vs ERFA nut80, 1000 BC-3000 AD", float(np.max(n80)), 0.01, "mas")
    check("IAU 1980 nutation vs ERFA nut06a (model difference), 1900-2100", float(np.max(n06[modern])), 25.0, "mas")
    check("GMST vs ERFA gmst82, 1900-2100", float(np.max(np.abs(g82[modern]))), 1.0, "mas")
    check("GMST vs ERFA gmst82, 1000 BC-3000 AD", float(np.max(np.abs(g82))), 20.0, "mas")
    check("GAST vs ERFA gst94 (1994 terms omitted), 1900-2100", float(np.max(np.abs(g94[modern]))), 4.0, "mas")


def section_gmst_exact():
    # IAU 1982 GMST (Meeus 12.4) in exact rational arithmetic, with the decimal
    # coefficients exact, at the exact UT1 instants Spp used; and the same
    # formula evaluated as written in double ("direct form").
    from fractions import Fraction as Fr
    rng = np.random.default_rng(17)
    jd = rng.uniform(1355818.5, 2816787.5, 3000)
    d = np.floor(jd).astype(np.int64)
    s = spp_pos(d, jd - d, [0.0] * jd.size, (0.0, 0.0, 0.0))
    a, r, c, dd = Fr("280.46061837"), Fr("360.98564736629"), Fr("0.000387933"), Fr(1, 38710000)
    split_err, direct_err = [], []
    for day, frac, gmst in zip(s["ut1_day"], s["ut1_frac"], s["gmst"]):
        D = Fr(int(day) - 2451545) + Fr(float(frac))
        T = D / 36525
        exact = (a + r * D + c * T * T - dd * T * T * T) % 360
        split_err.append(float((Fr(float(gmst)) - exact + 180) % 360 - 180) * ARCSEC)
        Dd = float(int(day) - 2451545) + float(frac)
        Td = Dd / 36525.0
        direct = (280.46061837 + 360.98564736629 * Dd + 0.000387933 * Td * Td - Td * Td * Td / 38710000.0) % 360.0
        direct_err.append(float((Fr(direct) - exact + 180) % 360 - 180) * ARCSEC)
    RESULTS["gmst_exact_arcsec"] = {"split": stats(split_err), "direct": stats(direct_err)}
    print(f"  info  GMST direct form vs exact: max {np.max(np.abs(direct_err)):.3g} arcsec", flush=True)
    check("GMST (split form) vs exact rational evaluation, 1000 BC-3000 AD", float(np.max(np.abs(split_err))),
          5e-7, "arcsec")


def section_hour_angle():
    rng = np.random.default_rng(11)
    years = np.arange(-1000.0, 3001.0, 50.0)
    jd = 2451545.0 + (years - 2000.0) * 365.25 + rng.uniform(0.0, 1.0, years.size)
    d = np.floor(jd).astype(np.int64)
    s = spp_pos(d, jd - d, [0.0] * jd.size, (0.0, 0.0, 0.0))

    def s_curve(lo, hi):
        e = np.linspace(lo, hi, int(abs(hi - lo) / 0.01) + 1)
        X, Y, Z = hc.mean_pole(e)
        g = (X * np.gradient(Y, e) - Y * np.gradient(X, e)) / (1.0 + Z)
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(e))])
        return e, hc.S0 - cum
    ef, sf = s_curve(2000.0, 3002.0)
    eb, sb = s_curve(2000.0, -1002.0)
    err = np.empty(jd.size)
    for i in range(jd.size):
        tt1, tt2 = s["tt_day"][i], s["tt_frac"][i]
        epj = erfa.epj(tt1, tt2)
        sl = np.interp(epj, ef, sf) if epj >= 2000.0 else np.interp(epj, eb[::-1], sb[::-1])
        (x, y, z, *_), _ = swe.calc(tt1 + tt2, swe.SUN, swe.FLG_SWIEPH | swe.FLG_EQUATORIAL | swe.FLG_XYZ
                                    | swe.FLG_J2000 | swe.FLG_ICRS | swe.FLG_NONUT)
        v = np.array([x, y, z]) / np.linalg.norm([x, y, z])
        X, Y = hc.true_cip(epj, tt1, tt2)
        vc = erfa.c2ixys(X, Y, sl) @ v
        h = np.degrees(erfa.era00(jd[i], 0.0)) - np.degrees(np.arctan2(vc[1], vc[0]))
        err[i] = angdiff_deg(s["gast"][i] - s["ra"][i], h) * ARCSEC
    RESULTS["hour_angle_vs_era_arcsec"] = {int(y): round(float(e), 3) for y, e in zip(years, err)}
    for lo, hi, lim in [(1600, 2100, 0.2), (1000, 3001, 1.5), (-1000, 1000, 30.0)]:
        sel = (years >= lo) & (years < hi)
        check(f"hour angle vs ERA + long-term precession (CIO), {lo}..{min(hi, 3000)}",
              float(np.max(np.abs(err[sel]))), lim, "arcsec")


# ----------------------------------------------------------------------------
# The Sun

def section_sun_de421():
    ts = load.timescale(builtin=True)
    eph = load_file(str(DE421))
    earth, sun = eph["earth"], eph["sun"]
    jd = np.arange(2415030.5, 2469800.5, 3.7) + 0.123456
    d = np.floor(jd).astype(np.int64)
    sp = spp_tt(d, jd - d)
    t = ts.tt_jd(d.astype(float), jd - d)
    app = earth.at(t).observe(sun).apparent(deflectors=())
    ra, dec, _ = app.radec(epoch="date")
    lat, lon, _ = app.frame_latlon(ecliptic_frame)
    dist = earth.at(t).observe(sun).distance().km
    d_ra = angdiff_deg(sp[:, 0], ra._degrees) * ARCSEC * np.cos(np.radians(dec.degrees))
    d_dec = (sp[:, 1] - dec.degrees) * ARCSEC
    d_lon = angdiff_deg(sp[:, 4], lon.degrees) * ARCSEC
    T = (jd - 2451545.0) / 36525.0
    slope, icept = np.polyfit(T, d_lon, 1)
    resid = d_lon - (icept + slope * T)
    RESULTS["sun_de421"] = {"separation_arcsec": stats(np.hypot(d_ra, d_dec)), "lon_drift_arcsec_per_cy": slope,
                            "lon_offset_j2000_arcsec": icept, "lon_residual_arcsec": stats(resid)}
    check("apparent Sun vs DE421 (IAU 2006/2000A), separation, 1900-2050",
          float(np.max(np.hypot(d_ra, d_dec))), 0.35, "arcsec")
    check("ecliptic longitude vs DE421 after the FK5 equinox drift, RMS", float(np.std(resid)), 0.02, "arcsec")
    check("FK5 equinox drift vs IAU 2006, |rate - 0.274|", abs(float(slope) - 0.274), 0.02, "arcsec/cy")
    check("ecliptic latitude vs DE421", float(np.max(np.abs(sp[:, 5] - lat.degrees)) * ARCSEC), 0.1, "arcsec")
    check("distance vs DE421", float(np.max(np.abs(sp[:, 3] * 1.495978707e8 - dist))), 15.0, "km")

    # Light time and aberration together, from the geometric to the apparent place,
    # against the constant-of-aberration term -20.4898"/R that Spp applies.
    glat, glon, _ = (sun - earth).at(t).frame_latlon(ecliptic_frame)
    d_lon_ref = angdiff_deg(lon.degrees, glon.degrees) * ARCSEC
    d_lat_ref = (lat.degrees - glat.degrees) * ARCSEC
    ab = -20.4898 / sp[:, 3] - d_lon_ref
    RESULTS["aberration_arcsec"] = {"longitude": stats(ab), "latitude_ref": stats(d_lat_ref)}
    print(f"  info  aberration term vs DE421: longitude max {np.max(np.abs(ab)):.4f}, "
          f"latitude max {np.max(np.abs(d_lat_ref)):.5f} arcsec", flush=True)
    check("light time + aberration vs -20.4898\"/R, longitude", float(np.max(np.abs(ab))), 0.012, "arcsec")


def section_sun_de431():
    jd = np.arange(1355818.5, 2816787.5, 36.7) + 0.271828
    d = np.floor(jd).astype(np.int64)
    sp = spp_tt(d, jd - d)
    ref = np.array([swe.calc(float(x), swe.SUN, swe.FLG_SWIEPH | swe.FLG_EQUATORIAL)[0][:2] for x in jd])
    years = 2000.0 + (jd - 2451545.0) / 365.25
    sep = np.hypot(angdiff_deg(sp[:, 0], ref[:, 0]) * np.cos(np.radians(ref[:, 1])), sp[:, 1] - ref[:, 1]) * ARCSEC
    RESULTS["sun_de431_separation_arcsec"] = stats(sep)
    modern = (years >= 1900) & (years < 2100)
    check("apparent Sun vs DE431, 1900-2100", float(np.max(sep[modern])), 0.5, "arcsec")
    check("apparent Sun vs DE431, 1000 BC-3000 AD", float(np.max(sep)), 4.0, "arcsec")


def section_topocentric():
    ts = load.timescale(builtin=True)
    eph = load_file(str(DE421))
    earth, sun = eph["earth"], eph["sun"]
    rng = np.random.default_rng(7)
    worst_alt = worst_az = 0.0
    for name, obs in OBSERVERS.items():
        mjd = rng.uniform(51544.0, float(C04["mjd"][-1]) - 1.0, 1000)
        step = np.abs(iers_dut1_at(np.floor(mjd) + 1) - iers_dut1_at(np.floor(mjd))) < 0.5   # no leap second
        mjd = mjd[step]
        dut1 = iers_dut1_at(mjd)
        jd = mjd + 2400000.5
        d = np.floor(jd).astype(np.int64)
        s = spp_pos(d, jd - d, dut1, obs)
        t = ts.tt_jd(s["tt_day"], s["tt_frac"])
        t.delta_t = s["tt_utc"] - dut1
        alt, az, _ = (earth + wgs84.latlon(obs[0], obs[1], elevation_m=obs[2])).at(t).observe(sun).apparent(
            deflectors=()).altaz()
        up = alt.degrees > -5.0
        worst_alt = max(worst_alt, float(np.max(np.abs(s["alt_geom"] - alt.degrees)[up])) * ARCSEC)
        worst_az = max(worst_az, float(np.max(np.abs(angdiff_deg(s["az"], az.degrees) * np.cos(alt.radians))[up]))
                       * ARCSEC)
    check("topocentric altitude vs DE421, 8 sites, 2000-2026", worst_alt, 0.1, "arcsec")
    check("topocentric azimuth (x cos alt) vs DE421, 8 sites, 2000-2026", worst_az, 0.1, "arcsec")


# ----------------------------------------------------------------------------
# Events

def saemundsson_apparent(h, pressure_mbar, temp_c=10.0):
    r = 1.02 / np.tan(np.radians(h + 10.3 / (h + 5.11))) + 0.0019279
    return h + r * (pressure_mbar / 1010.0) * (283.0 / (273.0 + temp_c)) / 60.0


def isa_pressure(elev_m):
    return 1010.0 * (1.0 - 2.25577e-5 * elev_m) ** 5.25588


def section_events(step_days=5):
    ts = load.timescale(builtin=True)
    eph = load_file(str(DE421))
    earth, sun = eph["earth"], eph["sun"]
    out = {}
    lim = {"transit": 0, "limb": 0, "limb_polar": 0, "twilight": 0, "twilight_polar": 0,
           "types": 0, "own_limb": 0, "own_twilight": 0, "own_transit": 0}

    for name, (lat, lon, elev) in OBSERVERS.items():
        topo = earth + wgs84.latlon(lat, lon, elevation_m=elev)
        pressure = isa_pressure(elev)
        dates = [tuple(int(v) for v in erfa.jd2cal(2400000.5, m)[:3]) + (m,)
                 for m in range(51544, int(C04["mjd"][-1]) - 2, step_days)]
        dut1 = iers_dut1_at(np.array([q[3] for q in dates]) + 0.5 - lon / 360.0)
        rows = run_spp([f"day {y} {m} {d} {lat} {lon} {elev} {fnum(u)}" for (y, m, d, _), u in zip(dates, dut1)])

        def tt_of(day, frac, u):
            mjd = np.floor(day + frac - 2400000.5)
            tai = LEAP_VAL[np.searchsorted(LEAP_MJD, mjd, side="right") - 1]
            t = ts.tt_jd(day.astype(float), frac + (TT_TAI + tai) / DAY_S)
            t.delta_t = TT_TAI + tai - u
            return t

        def apparent(day, frac, u):
            return topo.at(tt_of(day, frac, u)).observe(sun).apparent(deflectors=())

        def upper_limb(day, frac, u):
            alt, _, dist = apparent(day, frac, u).altaz()
            return saemundsson_apparent(alt.degrees + 959.63 / dist.au / ARCSEC, pressure)

        def almanac(day, frac, u):
            return apparent(day, frac, u).altaz()[0].degrees + 50.0 / 60.0

        def civil(day, frac, u):
            return apparent(day, frac, u).altaz()[0].degrees + 6.0

        def geo_hour_angle(day, frac, u):
            t = tt_of(day, frac, u)
            ra, _, _ = earth.at(t).observe(sun).apparent(deflectors=()).radec(epoch="date")
            return (t.gast + lon / 15.0 - ra.hours + 12.0) % 24.0 - 12.0

        def solve(day, frac, u, f):
            # Secant from Spp's instant to 1e-8 s; roots more than 2 h away are dropped.
            lost = np.zeros(frac.shape, dtype=bool)
            x0, x1 = frac.copy(), frac + 1.0 / DAY_S
            f0, f1 = f(day, x0, u), f(day, x1, u)
            for _ in range(40):
                x2 = x1 - f1 * (x1 - x0) / np.where(f1 != f0, f1 - f0, 1.0)
                lost |= ~np.isfinite(x2) | (np.abs(x2 - frac) > 2.0 / 24.0)
                x0, f0, x1 = x1, f1, np.where(lost, frac, x2)
                f1 = f(day, x1, u)
                if np.max(np.where(lost, 0.0, np.abs(x1 - x0))) * DAY_S < 1e-8:
                    break
            lost |= np.abs(x1 - x0) * DAY_S >= 1e-6
            return np.where(lost, np.nan, x1)

        typ = np.array([int(r[1]) for r in rows])
        tr_d = np.array([float(r[2]) for r in rows])
        tr_f = np.array([float(r[3]) for r in rows])
        res = {"transit_s": stats((tr_f - solve(tr_d, tr_f, dut1, geo_hour_angle)) * DAY_S)}
        lim["transit"] = max(lim["transit"], res["transit_s"]["max_abs"])

        # Transit in Spp's own pipeline: geocentric hour angle at the instant.
        p = spp_pos(tr_d, tr_f, dut1, (lat, lon, elev))
        lim["own_transit"] = max(lim["own_transit"],
                                 float(np.max(np.abs(angdiff_deg(p["gast"] + lon, p["ra"])))))

        for col, label in [(4, "rise"), (6, "set")]:
            have = np.array([r[col] != "nan" for r in rows])
            d = np.array([float(r[col]) for r in rows])[have]
            f = np.array([float(r[col + 1]) for r in rows])[have]
            u = dut1[have]
            if d.size == 0:
                continue
            ul = solve(d, f, u, upper_limb)
            al = solve(d, f, u, almanac)
            e_ul = (f - ul) * DAY_S
            res[f"{label}_vs_upper_limb_s"] = stats(e_ul)
            # Offset from the almanac convention (-0.8333°, geometric), and the offset the two
            # definitions predict with the reference geometry. Their difference is e_ul.
            res[f"{label}_vs_almanac_s"] = stats((f - al) * DAY_S)
            res[f"{label}_almanac_offset_predicted_s"] = stats((ul - al) * DAY_S)
            res[f"{label}_unresolved"] = int(np.sum(np.isnan(ul)))
            key = "limb_polar" if name in POLAR else "limb"
            lim[key] = max(lim[key], res[f"{label}_vs_upper_limb_s"]["max_abs"])
            a = run_spp([f"alt {int(x)} {fnum(y)} {fnum(z)} {lat} {lon} {elev}" for x, y, z in zip(d, f, u)])
            lim["own_limb"] = max(lim["own_limb"], max(abs(float(r[4])) for r in a) * ARCSEC)

            # Civil twilight: the center 6° below the horizon, geometric there.
            dirn = -1 if label == "rise" else 1
            sel = [i for i, h in enumerate(have) if h]
            cr = run_spp([f"cross {dates[i][0]} {dates[i][1]} {dates[i][2]} {lat} {lon} {elev} {fnum(dut1[i])} -6 {dirn}"
                          for i in sel])
            okc = np.array([r[1] != "nan" for r in cr])
            if np.any(okc):
                cd = np.array([float(r[1]) for r in cr])[okc]
                cf = np.array([float(r[2]) for r in cr])[okc]
                cu = dut1[np.array(sel)[okc]]
                tw = (cf - solve(cd, cf, cu, civil)) * DAY_S
                res[f"civil_{label}_s"] = stats(tw)
                key = "twilight_polar" if name in GRAZING_TWILIGHT else "twilight"
                lim[key] = max(lim[key], res[f"civil_{label}_s"]["max_abs"])
                a = run_spp([f"alt {int(x)} {fnum(y)} {fnum(z)} {lat} {lon} {elev}" for x, y, z in zip(cd, cf, cu)])
                lim["own_twilight"] = max(lim["own_twilight"], max(abs(float(r[3]) + 6.0) for r in a) * ARCSEC)

        # Day type: does the refracted upper limb stay above or below 0 within 12 h of transit?
        grid = np.linspace(-0.5, 0.5, 145)
        ulg = np.empty((len(rows), grid.size))
        for a0 in range(0, len(rows), 100):
            b0 = min(len(rows), a0 + 100)
            ulg[a0:b0] = upper_limb(np.repeat(tr_d[a0:b0], grid.size), (tr_f[a0:b0, None] + grid).ravel(),
                                    np.repeat(dut1[a0:b0], grid.size)).reshape(b0 - a0, grid.size)
        ref_type = np.where(ulg[:, 72] < 0.0, 2, np.where(np.all(ulg > 0.0, axis=1), 1, 0))
        res["day_type_disagreements"] = int(np.sum(ref_type != typ))
        lim["types"] += res["day_type_disagreements"]
        out[name] = res
        print(f"    events: {name} done", flush=True)

    RESULTS["events"] = out
    check("transit vs geocentric apparent hour angle = 0 (DE421)", lim["transit"], 0.005, "s")
    check("transit in Spp's own pipeline: |hour angle|", lim["own_transit"], 1e-9, "deg")
    check("sunrise/sunset vs upper limb definition (DE421 geometry), 0°-64°", lim["limb"], 0.05, "s")
    check("sunrise/sunset vs upper limb definition (DE421 geometry), 69°-78°", lim["limb_polar"], 1.0, "s")
    check("sunrise/sunset in Spp's own pipeline: |refracted upper limb|", lim["own_limb"], 1e-5, "arcsec")
    check("civil twilight vs DE421 (center at -6°), 0°-49°", lim["twilight"], 0.05, "s")
    check("civil twilight vs DE421 (center at -6°), 64°-78°", lim["twilight_polar"], 1.0, "s")
    check("civil twilight in Spp's own pipeline: |altitude + 6°|", lim["own_twilight"], 1e-5, "arcsec")
    check("day type disagreements with the upper limb definition", lim["types"], 0, "days", exact=True)


def section_instant_split():
    # Each split {day, frac} with frac in [-100, 100] against its normalized split,
    # {day + floor(frac), frac - floor(frac)}, which is the same instant exactly.
    rng = np.random.default_rng(13)
    raw, norm = [], []
    for _ in range(2000):
        day = int(rng.integers(1355818, 2816787))
        frac = float(rng.uniform(-100.0, 100.0))
        whole = float(np.floor(frac))
        raw.append(f"tt {day} {fnum(frac)}")
        norm.append(f"tt {day + int(whole)} {fnum(frac - whole)}")
    a, b = run_spp(raw), run_spp(norm)
    check("results that differ from the normalized split (2000 splits, frac in [-100, 100])",
          sum(1 for x, y in zip(a, b) if x != y), 0, "splits", exact=True)


def section_libm_dispatch():
    # The same requests with glibc's FMA/AVX2 math variants on (default) and off.
    import os
    rng = np.random.default_rng(3)
    lines = []
    for _ in range(20000):
        jd = rng.uniform(1355818.5, 2816787.5)
        d = int(np.floor(jd))
        lines.append(f"pos {d} {float(jd - d)!r} {float(rng.uniform(-0.9, 0.9))!r} "
                     f"{float(rng.uniform(-89, 89))!r} {float(rng.uniform(-180, 180))!r} 100")
    inp = "\n".join(lines) + "\n"
    env = dict(os.environ, GLIBC_TUNABLES="glibc.cpu.hwcaps=-AVX2,-FMA,-AVX2_Usable,-FMA_Usable")
    a = subprocess.run([str(EXE)], input=inp, capture_output=True, text=True, check=True).stdout.splitlines()
    b = subprocess.run([str(EXE)], input=inp, capture_output=True, text=True, check=True, env=env).stdout.splitlines()
    differ, worst = 0, 0.0
    angles = [POS_COLS.index(c) for c in ("gmst", "gast", "ra", "dec", "alt_geom", "az", "alt_refr")]
    for x, y in zip(a, b):
        if x != y:
            differ += 1
            p, q = x.split()[1:], y.split()[1:]
            worst = max(worst, max(abs(float(angdiff_deg(float(p[i]), float(q[i])))) for i in angles))
    RESULTS["libm_dispatch"] = {"positions_differing": differ, "of": len(lines), "max_angle_deg": worst}
    print(f"  info  glibc FMA/AVX2 variants off: {differ} of {len(lines)} positions differ, "
          f"max {worst:.3g} deg", flush=True)
    check("glibc FMA/AVX2 variants on vs off, largest angle difference", worst * ARCSEC, 1e-8, "arcsec")


SECTIONS = [("time", section_time), ("delta_t_pieces", section_delta_t_pieces), ("calendar", section_calendar), ("earth", section_earth),
            ("gmst_exact", section_gmst_exact), ("hour_angle", section_hour_angle), ("sun_de421", section_sun_de421), ("sun_de431", section_sun_de431),
            ("topocentric", section_topocentric), ("events", section_events), ("instant_split", section_instant_split),
            ("libm_dispatch", section_libm_dispatch)]


def main():
    global EXE, C04, LEAP_MJD, LEAP_VAL
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("spp_dump")
    ap.add_argument("--ephe", default=str(HERE.parent / "reference" / "ephe"))
    ap.add_argument("--json", default=None)
    ap.add_argument("--only", default=None, help="comma-separated sections: " + ",".join(n for n, _ in SECTIONS))
    args = ap.parse_args()
    EXE = Path(args.spp_dump).resolve()
    ensure_ephe(Path(args.ephe))
    C04 = load_c04()
    LEAP_MJD, LEAP_VAL = load_leap_seconds()
    only = set(args.only.split(",")) if args.only else None
    for name, fn in SECTIONS:
        if only and name not in only:
            continue
        print(f"[{name}]", flush=True)
        fn()
    if args.json:
        Path(args.json).write_text(json.dumps({"checks": CHECKS, "measurements": RESULTS}, indent=1))
    failed = [c for c in CHECKS if not c["pass"]]
    print(f"{len(CHECKS) - len(failed)}/{len(CHECKS)} checks passed")
    print("FAIL" if failed else "PASS")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
