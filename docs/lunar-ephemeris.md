# Lunar ephemeris, sky coordinates and phases

Spp evaluates the complete analytical ELP/MPP02 Moon model. It
provides geometric position and analytical velocity, apparent observer coordinates,
illumination, principal phase times and lunar horizon events, accessible through
the CLI. Surface/libration coordinates use an independent
[rotation provider](lunar-rotation.md). A separately selectable
[numerical provider](lunar-integration.md) propagates a finite window from a
kernel seed. Existing solar calculations keep their existing pipeline.

Hermite interpolation avoids repeated full-series evaluation in observations
and phase searches. The [Moon overview](moon.md) covers the complete API,
CLI, OpenGL application, model choices and measured limits.

## Try it yourself

Build the CLI with Ctrl+Shift+B in VS Code, then run **Tasks: Run Task → CLI: Run
Debug** and enter:

```text
moon 2026-11-03 12:00
```

Or use the project terminal:

```bash
cmake --build --preset debug --target spp_cli
./build/debug/spp moon 2026-11-03 12:00
./build/debug/spp moon 2026-11-04 12:00
./build/debug/spp moon 2026-11-03 12:00 --fit=llr
```

Change the date/time to explore the Moon's changing distance, speed and
geometric direction. The output also displays the full position and velocity
vectors and identifies the fitted set and coordinate frame.

The Moon command's calendar input is **TDB (ephemeris time)**, explicitly labeled
in the output and usage. It is not UTC or local clock time. Longitude/latitude
are geometric ecliptic coordinates in the J2000 frame, not altitude/azimuth
for an observer. The usual Spp historical-calendar convention applies. TDB
clock fields are validated; leap-second notation is not accepted.

### Altitude and azimuth for your location

```bash
./build/debug/spp moon sky 2026-11-03 12:00 48.8566 2.3522 35
./build/debug/spp moon sky 2026-11-03 18:00 48.8566 2.3522 35
./build/debug/spp moon sky 2026-11-03 12:00 48.8566 2.3522 35 --no-refraction
```

These examples use Paris. Replace latitude, longitude and optional elevation
(metres above the WGS84 ellipsoid) with your location. Longitude is positive
east. **Sky input is UTC**, not your local clock time. The date and clock are
strictly validated. The output identifies the fitted set and time scale.
ELP orbital commands accept `--fit=llr` or `--fit=de405`. Raw `moon rotation`
is independent of the orbital fit. Numerical selection has its own explicit
kernel, seed and propagation options.

Altitude is angular height above the horizon: +90 degrees is overhead, zero
is the horizon, and a negative value is below it. It refers to the Moon's
**center**. Azimuth is the compass direction: 0 north, 90 east, 180 south,
270 west. Refraction is enabled by default using the existing atmospheric
model, with standard pressure scaled to elevation. `--pressure=`, `--temp=`
and `--no-refraction` work as for the Sun. The unrefracted center altitude is
also printed; it still includes light time, parallax and aberration.

The sky output also includes apparent topocentric right ascension/declination
of date, observer distance, angular diameter, and **geocentric** phase and
illumination. A negative center altitude alone is not a rigorous limb-visibility
or moonrise test. Terrain, horizon dip and extinction are not modeled.
`--dut1=S` supplies UT1 minus UTC for Earth rotation. Without it Spp estimates
DUT1, explicitly labeled in the output.

### Illumination and upcoming phases

```bash
./build/debug/spp moon phase 2026-11-03 12:00
./build/debug/spp moon phase 2026-11-24 12:00
```

Phase input is also **UTC**. The command prints the lit fraction of the disk,
waxing/waning state, elongation from the Sun, phase angle, and principal phase
times during the following 35 days. Illumination uses spherical reflection
geometry at the retarded lunar position, including the finite solar distance
and incoming-light aberration at the Moon. The outgoing observer ray is
astrometric; apparent sky elongation and the longitude cycle are separate
quantities. Eclipses and lunar surface relief are excluded. It is not an
empirical estimate of brightness. The physical phase angle is approximately
zero for full Moon and 180 degrees for new Moon.

The phase cycle is apparent geocentric lunar longitude minus solar longitude,
wrapped to [0,360). Principal phases are roots at 0 (new), 90 (first quarter),
180 (full) and 270 (last quarter) degrees. They are not roots at exactly 0%,
50% or 100% illumination, because lunar latitude and finite solar distance
matter. Crescent/gibbous names describe the intervals; within 0.5 degrees of
a principal longitude the label says **Near** that phase. The event list gives
the actual modeled phase instants. No fixed 29.53-day cycle is substituted.

### Moonrise, moonset and transits

```bash
./build/debug/spp --help
./build/debug/spp moon day 2026-11-03 48.8566 2.3522 35
./build/debug/spp moon day 2026-11-03 48.8566 2.3522 35 --no-refraction
```

`moon day` searches the entire UTC day and prints all rise/set, upper/lower
meridian crossings and grazing contacts, plus hours above the horizon. A
missing rise or set is reported explicitly. Rise/set use the spherical upper
limb and the configured atmosphere at a level horizon. Elevation changes the
observer and default pressure; terrain and horizon dip are not inferred.
Upper transit means crossing the upper meridian; changing lunar declination
means it need not be the exact maximum of altitude.

The OpenGL app exposes the same day through **Today's sky** / **H**, including
an altitude curve, event tooltips and clickable event times. See
[moon-app.md](moon-app.md).

The separate **TDB-input** `moon rotation DATE TIME` command exposes the
IAU2009 pole/prime-meridian approximation. See [lunar-rotation.md](lunar-rotation.md)
for its independent provider, frame/rate contract and measured limits.

In VS Code, use **CLI: Run Debug** with any of the argument strings above,
omitting `./build/debug/spp`.

## State and time contract

`moon::geometric()` returns the Moon relative to Earth's center in the mean
dynamical ecliptic and equinox of J2000.0. It is not ICRF equatorial coordinates,
ecliptic-of-date coordinates, a topocentric state or an apparent direction.

- `position_km`: Cartesian x/y/z in kilometres.
- `velocity_km_per_day`: derivative in kilometres per TDB day, including the
  derivative of the model's moving-ecliptic-to-J2000 rotation.
- Input: `Instant<TDB>`, a split TDB Julian date. Equivalent splits normalize
  before evaluation. `Instant<TT>` is intentionally not accepted.
- Computational interval: -31 through +11 Julian centuries from J2000, inclusive
  (JD 1319270.0 through 2853320.0), encompassing Spp's intended chronology.
- Non-finite fractions/invalid fitted-set values throw `std::invalid_argument`;
  epochs outside the interval or unnormalizable extreme splits throw
  `std::out_of_range`.

The interval bounds prevent meaningless extrapolation and numerical overflow;
they are not an accuracy guarantee. The geometric command still consumes
actual TDB. Observer/phase functions use `Time::tdb()`, which derives TT using
Spp's existing leap-second/Delta-T policy and adds the geocentric seven-term
USNO Circular 179 equation 2.6 approximation. The stated TDB-TT approximation
error is about 10 microseconds over 1600-2200; outside that interval it is
extrapolated. It excludes the observer-dependent diurnal term. The conversion
preserves split-date precision and is available as `tdb_from_tt()` for typed TT
inputs. Future/ancient UTC conversions still inherit uncertainty in Delta-T
and future leap seconds. Do not relabel TT as TDB.

```cpp
#include <spp/moon.hpp>
#include <spp/model/elp_mpp02.hpp>

// J2000.0 expressed as a TDB Julian date.
const spp::Instant<spp::TDB> tdb{2451545, 0.0};
const auto state = spp::moon::geometric(tdb);
// Same untruncated position series, without evaluating unused velocity terms.
const auto position = spp::moon::position(tdb);

// Explicit parameter selection for reproducible comparisons.
const spp::model::elp_mpp02::Ephemeris llr{
    spp::model::elp_mpp02::Fit::LLR
};
const auto llr_state = spp::moon::geometric(tdb, llr);

const auto utc = spp::Time::from_utc(2026, 11, 3, 12, 0);
const spp::Observer observer{48.8566, 2.3522, 35.0};
const auto sky = spp::moon::observe(utc, observer);
const auto horizontal = spp::moon::horizontal(utc, observer);
const auto phase = spp::moon::phase(utc);
const auto next_phases = spp::moon::phases(utc, utc.shifted(35.0));

// All observation/phase operations can use a different ephemeris provider.
const auto llr_sky = spp::moon::observe(utc, observer, llr);
const auto llr_phase = spp::moon::phase(utc, llr);
const auto llr_events = spp::moon::phases(utc, utc.shifted(35.0), llr);

// Explicit cached selection: retain this instance across queries, per worker.
spp::model::elp_mpp02::CachedEphemeris fast_llr{spp::model::elp_mpp02::Fit::LLR};
const auto fast_sky = spp::moon::horizontal(utc, observer, fast_llr);
```

The default parameter set is DE405/DE406. Geometric `moon::geometric()` and
`moon::position()` still evaluate the full series directly. Default observation
and phase APIs share a bounded, per-thread interpolating provider, described
below; the `moon sky` and `moon phase` CLI commands use it too. Explicit
`Ephemeris` providers continue to request direct evaluation. Both parameter
sets retain the published adjustments. Provenance, redistribution
license and coefficient hashes are in [data/elp_mpp02](../data/elp_mpp02/README.md).

## Multiple bodies and implementations

Body API and model are separate. `moon::geometric(t, model)` uses a callable
provider returning the exact typed geometric state. Dispatch is resolved at
compile time; the model can be immutable or carry its own cache or propagation
state. The frame and origin are type parameters of `CartesianState`, so a
provider for a different origin/frame needs an explicit adapter.

This is the initial extension boundary, not an implementation of automatic
model selection. Future numerical dynamics, alternative lunar theories and
other moons can share state conventions and downstream mathematical operations
without requiring ELP-specific initialization. Selection policies still need
measured accuracy, date-range, workload and performance constraints.

`observe`, `coordinates`, `apparent`, `horizontal`, `phase`, `phases` and `horizon` accept the same
provider. They also accept a provider exposing
`GeocentricEclipticJ2000Position position(Instant<TDB>)` without requiring a
velocity. If both entry points exist, observations prefer `position()`; a
full-state-only callable still works by extracting its position. Frame/origin
compatibility is checked for both contracts. `moon::geometric()` keeps its
full-state contract and analytic velocity.

The compiled observation algorithms use a synchronous non-owning
reference to it, without allocating, copying its cache, or storing its lifetime.
The provider is evaluated at emission times as well as reception, so its valid
range must include the preceding light-travel interval. `apparent(Time)` returns
the geocentric apparent position; `observe().equatorial` is topocentric.

`phases(start,end,model,tolerance_s)` returns events in **[start,end)**, sorted
in time. Empty intervals return no events; reversed intervals and non-finite
or unsupported tolerances are errors. The default refinement tolerance is
0.01 seconds, permitted range 0.00001-60 seconds. This is numerical precision,
not an astronomical error bound. A one-day sampling bracket relies on the
Earth Moon's monotonically advancing synodic longitude; providers that violate
that assumption produce an error. It is not a general event finder for other
satellites. Illinois root refinement is shared with the existing solar solver.
Refinement uses offsets within the current one-day bracket, avoiding the loss
of subsecond precision from large floating-point offsets from the search start.

### Horizon search contract

```cpp
const auto start = spp::Time::from_utc(2026, 11, 3);
const spp::Observer paris{48.8566, 2.3522, 35};
const auto day = spp::moon::horizon(start, start.shifted(1), paris);
const auto center = spp::moon::horizon(start, start.shifted(1), paris,
    {.refraction=false}, {.altitude_deg=5, .upper_limb=false});
```

`horizon(start,end,observer[,model][,atmosphere][,options])` returns every
event in **[start,end)**, in time order, the above-horizon state immediately
after the start, and the duration above the target over the requested interval.
Empty intervals have no events and zero duration. No single-rise/single-set
assumption is made. `Coordinates` provides apparent equatorial/horizontal
coordinates, geometric center altitude, range, semidiameter and hour angle
without the additional geocentric phase/disk work of `observe`.

`HorizonOptions` selects the apparent horizon altitude (-80 to +80 degrees),
upper limb versus center, a 0.01-second numerical time tolerance (range
0.00001-60 seconds), and a 0.001-arcsecond angular contact tolerance (range
0.00001-1 arcsecond). These are solver controls, not physical prediction bounds.
Atmospheric pressure/temperature must give a monotonic terrestrial refraction
scale no greater than two; refraction is inverted once before the search.

The adaptive search operates on a smooth sine-altitude condition, using a
conservative terrestrial-lunar second-derivative bound to discard intervals
and certify monotonic crossings. It can subdivide same-sign endpoints to find
short visibility windows and stationary contacts. Contacts require an
interior stationary extremum or a stationary boundary contact, rather than
merely being close to the horizon. A sampled bound violation, non-finite
provider or unresolved continuous contact fails explicitly. This is a
terrestrial Moon algorithm: providers must obey its smooth-motion assumptions
and observer ranges of 250,000-600,000 km. It does not promise completeness
for arbitrary artificial satellites or pathological functions. Near-flat
contacts have angular uncertainty as well as numerical time tolerance.

`spp_moon_event_bench` couples timings with analytic synthetic cases,
one-minute dense lunar comparisons, direct/cache agreement and grazing,
short-window and interval-boundary gates. These numerical comparisons do not
independently establish physical moonrise timing accuracy.

### Work performed by each query

`moon::horizontal()` calculates only the observer's ray and horizontal
coordinates. It does not calculate a second geocentric ray, phase, illumination
or angular diameter. `moon::observe()` retains the complete result and shares
one reception-time position between the observer and phase calculations.

ELP's position-only evaluation retains all 35,901 terms and parameter
adjustments. Compile-time branches omit argument derivatives, differentiated
series, complementary trigonometric evaluations used for rates, and the
velocity/frame-rotation derivative. No coefficient truncation is used. A light-time iteration is
skipped only if both emission-time splits repeat exactly.

Horizontal-only queries require at most four provider position calls; complete
observations at most seven. With direct ELP each is a full-series evaluation;
with the default cache they usually share one interpolated interval. The previous
horizontal wrapper performed the complete observation, up to eight full-state
evaluations including velocities. These are workload counts, not measured
speedups. Alternative models that expose only a full state still calculate
whatever work their own callables require.

Principal-phase searches evaluate the longitude cycle without calculating
illumination at every sample/refinement. The `spp_bench` target, when built
with both `SPP_BUILD_BENCH=ON` and `SPP_BUILD_MOON=ON`, now includes matched
direct-state/position and complete-observation/horizontal cases, plus cached
repeat workloads and explicitly forced lunar-cache misses. Numerical and
performance results are recorded in [lunar-testing.md](lunar-testing.md).

### Interpolation policy

`model::elp_mpp02::CachedEphemeris` owns an immutable fitted-set selection and a
bounded cache of twenty-minute intervals. Cubic Hermite interpolation uses the
full analytical position and velocity at each endpoint. A fresh isolated
interval costs two full state evaluations; adjacent intervals can reuse a
shared endpoint. Hits use cubic arithmetic, including a differentiated cubic
for callers requesting velocity. This replaces repeated trigonometric series
evaluations at reception and nearby emission instants. The direct evaluator
and its analytical velocity remain available without interpolation.

The reusable `model::HermiteEphemeris<Model,Scale,IntervalsPerDay>` adapter
preserves the Cartesian origin/frame and input time-scale types. It owns its
provider and domain; no model-fit identifiers or body-specific physics are
embedded in its cache. Other bodies/models require their own interval/error
validation before adopting it. The generic adapter defaults to 288 intervals/day;
the ELP adapter selects 72/day (twenty minutes). It holds 256 intervals with
compact hash lookup and LRU replacement, allowing every slot to serve any date.
Storage is bounded, with no
runtime allocation or external trajectory data. Coefficients are generated
on demand from Spp's own full analytical model.

An instance has mutable cache state and must be confined to a worker or
externally synchronized. Default APIs use one instance per thread; explicit
instances isolate fitted sets and model state. `clear()` drops the segments
and shared endpoint memo. `source_evaluations()` reports actual source calls
since construction/clear. Equivalent split dates normalize before lookup;
partial edge intervals are clipped to the source domain and the exact final
endpoint is returned directly. Failed fills do not publish partial segments.

Interpolation deliberately adds a separate numerical approximation. New
checks enforce sampled position budgets of 1 mm during 2000-2100 and 0.1 m
across the full computational domain, alongside velocity, observation and
phase-time comparisons. The twenty-minute checks passed with sampled maximum
position differences of 0.807 mm in 2000-2100 and 0.028 m over the computational
domain, and phase-time differences below 0.000006 seconds. These are sample
results, not global error bounds. JPL/USNO checks exercise the cached
default path; the direct upstream comparison remains unchanged.

For isolated geometric queries, direct position evaluation avoids paying for
two endpoint states. For repeated observations or phase refinements, cached
states can be reused. The benchmark exposes both repeat and forced-miss costs;
no workload-independent speedup or automatic model-selection policy is claimed.

## Observation model

The observation pipeline performs these operations:

1. Convert UTC-bearing `Time` to TT, TDB and UT1, using each scale for its own purpose.
2. Rotate the geometric J2000 dynamical ecliptic state to date using IAU 2006
   P03 precession, with the J2000 frame-bias rotation removed. Apply Spp's
   existing IAU 1980 nutation and Laskar obliquity for true-of-date coordinates.
3. Compute the WGS84 observer vector and subtract it exactly, retaining lunar
   parallax and the changed topocentric distance. No solar parallax constant
   is reused for the Moon.
4. Iterate Newtonian light travel three times using the Moon's geocentric
   state plus Earth's heliocentric displacement. The Earth trajectory comes
   from the existing VSOP87D/FK5 geometric Sun, rotated back to fixed axes.
   Apply special-relativistic aberration using Earth's orbital velocity plus
   the observer's rotational velocity. Retarding only the geocentric Moon
   and adding ordinary annual aberration would be inconsistent.
5. Convert to local horizontal coordinates using UT1/GAST, then apply the
   shared atmospheric refraction model if requested.

This is a self-contained analytical observation model. It approximates an
inertial solar frame with heliocentric Earth states; it does not include solar
barycentric motion, gravitational deflection, measured polar motion or measured
celestial-pole offsets. P03 polynomials and the TDB approximation have more
limited precision domains than the geometric ephemeris's computational range.
There is no uniform sky-coordinate accuracy claim across the entire interval.
The Moon is modeled as a sphere of mean radius 1737.4 km for angular size.

`Observation.disk` supplies topocentric spherical illumination and orientation
for the phase viewer. Its bright-limb angle is measured eastward from celestial
north; subtract the parallactic angle to measure from the local zenith. The
geocentric almanac phase remains separately available in `Observation.phase`.
Undefined bright-limb/zenith directions have explicit flags. Stable `atan2`
separation and half-angle illumination avoid cancellation near alignment.

Illumination follows the physical `phi` convention, rather than Horizons'
approximate S-T-O angle. The incoming Sun direction uses lunar heliocentric
velocity, estimated from the already available reception/emission positions
and the centered Earth-velocity estimate. It needs no additional ephemeris
calls or velocity-only provider API. The short lunar secant and heliocentric
solar approximation are model limitations, tested through independent physical
phase/illumination fixtures. The bright-limb sky bearing retains the apparent
Sun/Moon formula (Meeus 48.5), checked against Horizons' sub-solar bearing.

Terrestrial observer inputs require finite coordinates, latitude within
[-90,90] and elevation between -12 km and 10,000 km. Atmosphere inputs require
finite nonnegative pressure and temperature above -273 degrees Celsius (the
shared refraction formula's singular point). These bounds
do not establish atmospheric-model accuracy at extreme elevations.

## Build

`SPP_BUILD_MOON=ON` is the default. ELP coefficients are compiled into the static
library; its geometric direct evaluation and cached state/position queries use no
runtime files, network or allocation. Event searches allocate their result and
search storage. Cached instances have per-instance mutable state;
the convenience observation APIs confine it to each thread.

```bash
cmake -S . -B build/moon -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/moon --target spp_lunar_check spp_moon_sky_check
ctest --test-dir build/moon --output-on-failure
```

`SPP_BUILD_VALIDATION` must be enabled for `spp_lunar_check` (it defaults to ON
when Spp is the top-level project). Set `-DSPP_BUILD_MOON=OFF` for a solar-only
build; the CLI then omits `moon`, and calls to lunar functions require a
different build of the library.
For cross-compilation without an emulator, `elp_gen` is built using
`SPP_HOST_CXX`, or supplied with `SPP_ELP_GEN_EXECUTABLE`, matching the existing
VSOP generator's host-tool scheme.

## Verification and limits

The [Moon speed and accuracy report](lunar-testing.md) provides one Release
configure/build/report workflow, a VS Code task, detailed numerical comparisons
and a Moon-only repeated benchmark. Performance verification includes direct
upstream, interpolation, JPL/USNO and renderer measurements. Consult that page
and the app documentation for the workload and limits of each measurement.

The committed CSV covers 72 TDB epochs for each fitted set, including J2000,
upstream example dates, both computational endpoints and seeded random dates
across the interval. It was produced by compiling the upstream evaluator
separately, not by Spp's evaluator. This verifies compatibility with that
implementation; it is not independent validation of the underlying theory
against JPL observations or modern DE440/DE441 ephemerides.

The checker also verifies analytical velocity by five-point differences of
positions, equivalent split dates, input errors, the default fitted set,
explicit stateful-provider dispatch and compile-time rejection of TT inputs
and incompatible state frames. Assertions use explicit failures in Release
builds as well as Debug builds.

The user-supplied preceding report passed position-only evaluation against the same
144 upstream states, full-state fallback, position-only providers, exact model
selection, observer-path consistency, reception-state reuse and call-count
limits. New interpolation checks additionally verify rates, disk orientation,
finite-distance topocentric illumination, bounds, and cache replacement.

To regenerate the fixtures from the pinned, unmodified upstream checkout:

```bash
python3 tools/lunar_reference.py /path/to/ElpMpp02 tools/reference/elp_mpp02_states.csv
```

The upstream directory must contain `ElpMpp02.cpp`, `ElpMpp02.h` and its 14
coefficient files. Regeneration needs Python and a C++20 compiler; ordinary Spp
builds need neither Python nor an upstream checkout. Inspect changes to the
fixtures before accepting them.

Before this follow-up, the observation checker used independent, committed Horizons DE441 observer
results at five epochs (2000, the 2024 eclipse and three November 2026 instants)
for Paris, Sydney and a high Arctic location, including positions below the
horizon. DUT1 is taken from each Horizons response. It separately checks five
geocentric apparent/phase cases. In the tested Release build the largest sky
direction difference was 1.016 arcseconds, equatorial difference 0.099 arcseconds,
range difference 0.022 km and illumination difference 0.0023 percentage points.
Horizons includes measured/predicted Earth-orientation corrections and uses
slightly different phase/illumination conventions. These sample maxima are
regression evidence, not bounds for untested dates, sites or atmospheres.

Before this follow-up, principal times were compared with USNO's minute-precision Universal Time tables
for 2000, 2024, 2026 and 2050. The checker validates all event kinds/order/counts
and longitude residuals. The largest difference was 41.25 seconds across 198
phases. Differences below a minute are required; reference quantization,
model differences and UT-versus-UTC prevent inferring subsecond accuracy from
that check.
The separate contracts exercise model selection, stateful providers, refraction,
observer elevation/poles, DUT1, invalid inputs, frame inverses and narrow/empty
phase intervals. Existing solar output is checked independently for regressions.

Regenerate observational fixtures without running Spp:

```bash
python3 tools/lunar_sky_reference.py --cache /tmp/lunar-reference-downloads
```

The cache retains unmodified API responses; inspect the coefficient/model,
Earth-orientation data and fixture changes before accepting regenerated values.
Future/historical Horizons sky rows are intentionally excluded where UTC/TT
or Earth-orientation policies differ. Those independent sky fixtures establish
neither physical moonrise timing nor millennia-wide apparent-position accuracy.
Horizon events have the separate numerical/performance coverage described above.

The additional `moon_horizons_phase_tt.csv` fixture uses 62 geocentric Moon/Sun
phase neighborhoods in **TT**. Its quadratic reference roots are checked
against a fourth independent sample per arc. The matched-time benchmark
removes UT/UTC/Delta-T ambiguity and measured about 0.21 s worst difference
against DE441, at the sampled 2000/2024/2026/2050 events. Its regeneration tool
does not run Spp:

```bash
python3 tools/lunar_phase_reference.py --cache /tmp/spp-moon-observation-reference
```

## Observation references

- [USNO Circular 179, equation 2.6](https://arxiv.org/pdf/astro-ph/0602086): TDB-TT approximation and its precision domain.
- [ERFA pfw06](https://github.com/liberfa/erfa/blob/master/src/pfw06.c): P03 Fukushima-Williams precession coefficients, attributed in `data/licenses/ERFA.txt`.
- [JPL aberration corrections](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/abcorr.html): light travel and observer motion.
- [USNO phase definitions](https://aa.usno.navy.mil/faq/moon_phases): apparent-longitude definitions of principal phases.
- [JPL Horizons manual](https://ssd.jpl.nasa.gov/horizons/manual.html): observer quantities, frames, phase conventions and Earth orientation.
- [USNO API](https://aa.usno.navy.mil/data/api): independent principal-phase reference tables.
