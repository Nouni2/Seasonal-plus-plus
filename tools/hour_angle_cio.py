"""Hour angle of the Sun from the Earth Rotation Angle, without sidereal-time polynomials.

H = ERA(UT1) + longitude - alpha_CIO, where alpha_CIO is the Sun's right ascension
measured from the Celestial Intermediate Origin. The CIP comes from the long-term
precession of Vondrak et al. (2011) (ERFA ltpb) with IAU 2000A nutation; the CIO
locator s is integrated numerically from its kinematic definition
    s(t) = s0 - integral_{J2000}^{t} (X dY/dt - Y dX/dt) / (1 + Z) dt.
The apparent geocentric direction of the Sun is taken from the Swiss Ephemeris
files (compressed JPL DE431) in the GCRS (ICRS axes, no precession or nutation).

Used by reference_check.py. Standalone:
    python -I hour_angle_cio.py <spp_dump> <swiss_ephemeris_dir>
prints, per epoch, the hour-angle differences Spp - CIO and Swiss Ephemeris - CIO.
"""

import subprocess
import sys

import erfa
import numpy as np
import swisseph as swe
from pathlib import Path

S0 = 94e-6 / 3600.0 * np.pi / 180.0   # s at J2000.0, rad


def mean_pole(epj):
    rpb = erfa.ltpb(epj)
    return rpb[..., 2, 0], rpb[..., 2, 1], rpb[..., 2, 2]


def cio_locator(epj_target, step_years=0.01):
    n = max(2, int(abs(epj_target - 2000.0) / step_years) + 1)
    e = np.linspace(2000.0, epj_target, n)
    X, Y, Z = mean_pole(e)
    dX = np.gradient(X, e)
    dY = np.gradient(Y, e)
    integrand = (X * dY - Y * dX) / (1.0 + Z)
    return S0 - np.trapezoid(integrand, e)


def true_cip(epj, tt1, tt2):
    rpb = erfa.ltpb(epj)
    eps = np.arccos(np.clip(np.dot(erfa.ltpequ(epj), erfa.ltpecl(epj)), -1.0, 1.0))
    dpsi, deps = erfa.nut06a(tt1, tt2)
    rnpb = erfa.numat(eps, dpsi, deps) @ rpb
    return rnpb[2, 0], rnpb[2, 1]


def main():
    exe = Path(sys.argv[1]).resolve()
    swe.set_ephe_path(sys.argv[2])
    jds = [1355900.3, 1538400.3, 1721500.3, 1904100.3, 2086300.3, 2269000.3, 2342000.3, 2415200.3,
           2451545.3, 2488000.3, 2524600.3, 2633000.3, 2816000.3]
    print(f"{'year':>8} {'Spp-CIO':>10} {'Swiss-CIO':>10}   (arcsec of hour angle; UT1 given, longitude 0)")
    for jd in jds:
        d = int(np.floor(jd))
        f = jd - d
        out = subprocess.run([str(exe)], input=f"pos {d} {f!r} 0 0 0 0\n", capture_output=True, text=True,
                             check=True).stdout.split()
        tt1, tt2 = float(out[5]), float(out[6])
        gast, ra_spp = float(out[10]), float(out[19])
        epj = erfa.epj(tt1, tt2)

        (x, y, z, *_), flag = swe.calc(tt1 + tt2, swe.SUN,
                                       swe.FLG_SWIEPH | swe.FLG_EQUATORIAL | swe.FLG_XYZ | swe.FLG_J2000
                                       | swe.FLG_ICRS | swe.FLG_NONUT)
        v = np.array([x, y, z]) / np.linalg.norm([x, y, z])
        X, Y = true_cip(epj, tt1, tt2)
        s = cio_locator(epj)
        v_cirs = erfa.c2ixys(X, Y, s) @ v
        alpha_cio = np.degrees(np.arctan2(v_cirs[1], v_cirs[0]))
        era = np.degrees(erfa.era00(jd, 0.0))
        h_cio = era - alpha_cio

        (ra_swe, *_), _ = swe.calc(tt1 + tt2, swe.SUN, swe.FLG_SWIEPH | swe.FLG_EQUATORIAL)
        h_swe = swe.sidtime(jd) * 15.0 - ra_swe
        h_spp = gast - ra_spp
        w = lambda a: ((a + 180.0) % 360.0 - 180.0) * 3600.0
        year = 2000.0 + (jd - 2451545.0) / 365.25
        print(f"{year:8.1f} {w(h_spp - h_cio):10.2f} {w(h_swe - h_cio):10.2f}")


if __name__ == "__main__":
    main()
