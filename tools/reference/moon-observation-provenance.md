# Lunar observation references

Generated independently of Spp on 2026-11-03 by
`tools/lunar_sky_reference.py`. The generator consumes only public reference
services; it never runs or imports the Spp implementation.

`moon_horizons.csv` contains JPL Horizons API v1.2 observer tables for Moon
301 relative to geocenter 500@399 or geodetic sites on Earth. Horizons reported
DE441 and Earth-orientation file **eop.261102.p270129**. Site positions:

| Fixture | East longitude | Geodetic latitude | Elevation, metres |
| --- | ---: | ---: | ---: |
| paris | 2.3522 | 48.8566 | 35 |
| south | 151.2093 | -33.8688 | 58 |
| polar | -156.7886 | 71.2906 | 0 |

API: https://ssd.jpl.nasa.gov/api/horizons.api

The request uses observer ephemerides, UTC dates, airless apparent positions,
ICRF input reference, decimal-degree angles, extra precision and quantities
2,4,7,10,13,16,20,23,24,31,43,49. Apparent RA/Dec are true-of-date, not J2000 ICRF
angles. DUT1 is copied from quantity 49 so comparisons use matched Earth
rotation times. Lunar size is quantity 13, range quantity 20, illumination
quantity 10, elongation quantity 23 and **physical phase angle quantity 43
(`phi`)**. Quantity 24's approximate S-T-O angle is retained in the raw response
for inspection and is not substituted for phi. Quantity 16 supplies the
sub-solar sky bearing at 0.01-degree precision. These additional references
were fetched into the `*-disk43.json` cached responses; the CSV now has 16
columns, including the bearing. Both the 15 topocentric disk illuminations and
five geocentric illuminations are checked.

Five retained epochs are 2000-01-01 12:00, 2024-04-08 18:00,
2026-11-03 12:00, 2026-11-04 00:00 and 2026-11-24 12:00 UTC. The generator
also requests historical/future points for inspection but excludes them from
the fixtures where time/orientation policies differ. It retains unmodified
API responses in the explicitly supplied cache directory. Later regenerated
Earth-orientation predictions can change; review differences rather than
blindly replacing the regression data.

`moon_usno_phases.csv` contains all 198 principal phases from USNO API v4.0.1
for 2000, 2024, 2026 and 2050. The endpoint is
`https://aa.usno.navy.mil/api/moon/phases/year?year=YEAR`. It reports Universal
Time at minute precision. The fixture preserves that precision; it does not
add fictitious seconds. API version can be rechecked in the cached responses.

The Horizons comparison allows differing analytical/measured pole models and
relativistic phase conventions. USNO comparisons allow reference quantization,
time-scale and model differences. Neither dataset establishes a general
astronomical accuracy bound. See `docs/lunar-ephemeris.md` for measured maxima
and implementation limits.

## Matched-TT phase reference

`moon_horizons_phase_tt.csv` was independently generated on 2026-11-03 by
`tools/lunar_phase_reference.py`, without evaluating Spp. Its SHA-256 is
`5313d4b67c825005107faa321ac38b892013392da097be1980bf7a6cc986f5ff`.

The source is JPL Horizons API observer tables for targets **301 (Moon)** and
**10 (Sun)**, geocenter `500@399`, **TT input and output**, quantity **31**,
airless, ICRF reference, extra precision and decimal degrees. Both bodies
reported **DE441**. Quantity 31 is the apparent observer-centered IAU76/80
ecliptic-of-date longitude, including light time, aberration and gravitational
deflection. The raw response labels TT explicitly and prints longitudes at
1e-7-degree resolution. Requests are batched sequentially in groups of 32
epochs and unmodified responses are cached under their request-URL hashes.

The 62 neighborhoods cover all 50 principal phases in 2026 and the first four
each in 2000, 2024 and 2050. USNO's rounded times only select a neighborhood;
they are not reused as reference event times. Each neighborhood has four
samples at -180, 0, +90 and +180 seconds relative to the USNO neighborhood
plus 90 seconds, explicitly interpreted as TT. The benchmark uses the first,
second and fourth samples for a quadratic longitude-difference root and
checks it against the unused +90-second point. The maximum unused-point
residual was 0.000476 arcsecond, comparable to table quantization.

The Spp event is converted back to TT for comparison. This avoids treating
uncertain future Delta-T or UTC/UT1 differences as lunar-model error. With a
0.0001-second Spp numerical root tolerance, sampled DE405-vs-DE441 residuals
were at most 0.211037 seconds, RMS 0.081283 seconds. Cached/direct phase times
differed by at most 0.00001782 seconds in that run. This is sampled model
comparison, not a global error bound or a guarantee of future UTC predictions.

Primary contracts: [Horizons quantities](https://ssd.jpl.nasa.gov/horizons/manual.html#observer-table),
[API time scales and TLIST](https://ssd-api.jpl.nasa.gov/doc/horizons.html),
and [USNO phase definitions](https://aa.usno.navy.mil/faq/moon_phases).

## Extended numerical-force fixtures (2026-11-05)

`gravity_de440.csv` and `earth_tides.csv` are generated by
`tools/gravity_reference.py` from the pinned DE440s comments and published
Moon GM. The script uses 70-digit Decimal arithmetic: ordinary rational
Legendre polynomials define each potential, and five-point numerical
derivatives supply acceleration. It loads no Spp implementation. The 28
gravity cases cover degree 2 and full degree/order, equatorial/oblique/polar
directions, and near-body distances. Eighteen separate tidal-mode gradients
cover all three response orders with lunar- and solar-distance sources.

The production harmonic evaluator uses Cartesian solid-harmonic recurrences
with analytic gradients. Separate zonal and quadrupole paths are checked
against the same independent potential fixtures. Tidal modes are additionally
checked against their isotropic zero-lag limit about different poles.

The numerical benchmark now also evaluates the same 62 matched-TT Horizons
phase arcs described above. It constructs its own Moon trajectories for four
±190-day windows; source/reference correlation through the initial seed is
explicit. Sample maximum/RMS phase residuals are 0.268178/0.128456 seconds.
This differs from the ELP result; orbital and phase accuracy are separate
reported quantities.

Force ablation and short-lag prediction controls are described in
[lunar integration](../../docs/lunar-integration.md). The delay-source predictor
is local and bounded to ±0.25 day, rather than using an external Moon history.
The reference kernel is queried for truth only by the benchmark. Physical
forces consume the propagated state, EMB/planet trajectories and response
constants. Primary model sources are
[Park et al. (2021)](https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf) and
[Folkner et al. (2014)](https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf).
