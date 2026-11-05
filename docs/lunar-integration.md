# Numerical lunar provider

Spp has two independently selectable Moon approaches: the self-contained
ELP/MPP02 analytical series, with its interpolation adapter, and a finite-window
numerical trajectory driven by prescribed planetary ephemerides. ELP remains
the default. Selection is explicit; there is no automatic fallback or reseeding.

The numerical model includes Earth/Moon attraction, the Sun and Mercury,
Venus, Mars, Jupiter, Saturn, Uranus, Neptune and Pluto system barycenters,
major-body first post-Newtonian corrections, Earth zonal gravity through degree
5, and lunar gravity through degree/order 6. Delayed Earth and lunar tides and
lunar spin deformation modify the degree-2 fields. Figure/point-mass reactions
with the Sun through Saturn, solar J2, photon pressure and solar Lense–Thirring
frame dragging are included. It uses DE440 mass/figure/response
constants. Planetary and Earth–Moon barycenter
trajectories come from the supplied SPK; Earth and Moon are derived from that
barycenter and the propagated separation at each force evaluation. The kernel's
Earth/Moon trajectories are not queried after seeding. This is a reduced dynamical model:
it omits asteroids, integrated physical libration and a fully coupled
solar-system integration.
It does not reproduce DE440's complete dynamical model.

## Local data and commands

The existing local kernel is available at
`context/Spp Moon/reference-data/de440s.bsp`, also accessible from Windows
through the project on drive A:. It is runtime data for the numerical option;
ordinary builds and default ELP queries do not need or download it.

From the project root on the NAS:

```bash
cmake --preset moon
cmake --build --preset moon --target spp_cli spp_moon
./build/moon-Linux/spp moon sky 2026-11-04 12:00 48.8566 2.3522 35 --model=integrated "--kernel=context/Spp Moon/reference-data/de440s.bsp" --seed=2026-11-03T12:00
./build/moon-Linux/spp moon phase 2026-11-04 12:00 --model=integrated "--kernel=context/Spp Moon/reference-data/de440s.bsp" --seed=2026-11-03T12:00
./build/moon-Linux/spp_moon 2026-11-04 12:00 48.8566 2.3522 35 --model=integrated "--kernel=context/Spp Moon/reference-data/de440s.bsp" --seed=2026-11-03T12:00 --events --grid
```

From a **local Windows MSYS2 UCRT64 terminal**, after installing the graphics
dependencies in [moon-app.md](moon-app.md):

```bash
cd /a/Documents/Code/C-C++/Seasonal-plus-plus
cmake --preset moon
cmake --build --preset moon --target spp_cli spp_moon
./build/moon-Windows/spp.exe moon sky 2026-11-04 12:00 48.8566 2.3522 35 --model=integrated "--kernel=context/Spp Moon/reference-data/de440s.bsp" --seed=2026-11-03T12:00
./build/moon-Windows/spp_moon.exe 2026-11-04 12:00 48.8566 2.3522 35 --model=integrated "--kernel=context/Spp Moon/reference-data/de440s.bsp" --seed=2026-11-03T12:00 --events --grid
```

These observer coordinates are Paris. Observation timestamps are UTC;
**`--seed` is TDB**. The seed names an actual initial state read once from the
kernel. In these examples the observation is a day after that state, so it
requires genuine propagation. `--span-days=45` is the CLI default and
`--span-days=90` the GUI default, in each direction from the seed. The GUI's
phase calendar searches beyond the displayed instant; it needs sufficient
window on both sides. Requests outside the propagated window fail explicitly.
The CLI reports step/force counts, and the GUI labels the stored trajectory.
`--fit` and `--direct` belong to ELP. In integrated GUI mode F/D do not change
the provider. Windows execution still requires verification on Windows.

## C++ API

```cpp
#include <spp/spp.hpp>
#include <spp/model/lunar_numerical.hpp>

const spp::model::spk::Kernel kernel("context/Spp Moon/reference-data/de440s.bsp");
const spp::Instant<spp::TDB> seed{2461348, 0.0}; // 2026-11-03 12:00 TDB
const spp::model::lunar_numerical::Ephemeris orbit(
    kernel, seed, {seed.day-45, seed.frac}, {seed.day+45, seed.frac});
const auto t = spp::Time::from_utc(2026, 11, 4, 12);
const spp::Observer paris{48.8566, 2.3522, 35};
const auto state = orbit(t.tdb());
const auto sky = spp::moon::observe(t, paris, orbit);
const auto events = spp::moon::horizon(t, t.shifted(1), paris, orbit);
const auto phases = spp::moon::phases(t, t.shifted(35), orbit);
```

The returned orbital frame/origin and units match ELP: geocentric mean
ecliptic/equinox J2000, kilometres and kilometres per TDB day. `icrf(t)` exposes
the typed geocentric ICRF state. A fixed frame-bias adapter performs the
conversion. `stats()` reports accepted/rejected steps and force evaluations;
`storage_bytes()` reports allocated dense-trajectory storage, excluding the
shared kernel. `environment_bytes()` reports the separate propagation-input
cache, and `environment_samples()` counts the exact time-only samples used to
construct it. Configuration and results belong to each instance.

`IntegrationOptions` controls relative, position and velocity error scales,
initial/minimum/maximum step and a hard step budget. The defaults use 1e-12
relative tolerance, 1e-7 km position and 1e-7 km/day velocity scales, with a
30-minute maximum step. These are local numerical controls, not physical
accuracy promises. `ForceOptions` exposes relativity, body figures, maximum
degrees, Earth/lunar tidal responses and external figure interactions.
Setting `earth_tides=false`, `lunar_tides=false`,
`external_figure_interactions=false`, `solar_oblateness=false`,
`solar_radiation_pressure=false`, `solar_frame_dragging=false`,
`earth_degree=2`, `moon_degree=2` selects
the earlier quadrupole family for comparisons. A second constructor accepts an explicit ICRF/geocentric
initial state and `Forces`, avoiding dependence on a kernel Moon seed.

### Propagation performance

The default `ForceOptions::environment_cache=true` prepares immutable,
degree-12 Chebyshev segments of at most four days for the prescribed planetary
states, pole/orientation matrices, delayed solar inputs, and the planet-only
Newtonian acceleration/potential needed by the 1PN equations. It fits residuals
about central samples to preserve small variations in large coordinates.
Every stage still evaluates the Moon-dependent forces from its own trial
position and velocity. All gravity degrees, tides, relativity and solar terms
remain enabled; the RK step limits and tolerances are unchanged.

This cache contains neither a kernel lunar trajectory nor a fitted propagated
Moon trajectory. It is shared read-only and does not depend on query order.
SPK segment-priority boundaries use direct evaluation, as do force requests
outside the prepared window. Exact endpoint samples preserve coverage errors.
Auxiliary storage is bounded at 4096 segments; longer construction windows
use direct inputs. Accuracy measurements below apply to the sampled DE440s
workloads, rather than guaranteeing interpolation accuracy for arbitrary SPKs.

To compare with direct inputs through the C++ API:

```cpp
spp::model::lunar_numerical::ForceOptions direct;
direct.environment_cache = false;
const spp::model::lunar_numerical::Ephemeris reference(
    kernel, seed, {seed.day-45, seed.frac}, {seed.day+45, seed.frac}, {}, direct);
```

A standalone `Forces` object evaluates inputs directly until
`forces.prepared(first,last)` returns a prepared copy. `Ephemeris` prepares
its force automatically. The benchmark compares cached/direct forces and
trajectories separately from comparisons against DE440s.

## Gravity and tidal response

The fixed lunar PA coefficients come from the pinned DE440s comments. A shared
Cartesian harmonic evaluator computes analytic gradients without a pole
singularity. Earth uses a faster axisymmetric Legendre evaluation. The Sun and
other inner/giant planets also exert forces on the extended Earth and Moon;
the resulting reactions affect their relative acceleration. Earth uses its
selected zonal degree for those interactions; distant lunar interactions use
degree 2. Higher lunar degrees at planetary distances are neglected.

Earth tides have separate long-period, diurnal and semidiurnal Love numbers and
orbital/rotational lags in `earth_tide`. The force derives from the three
trace-free quadrupole-potential modes. Equal Love numbers with zero lag reduce
to the standard isotropic equilibrium tide. The passive rotation in DE430's
equation 31 corresponds to a positive vector rotation by the rotational lag.

Lunar deformation uses the delayed Earth direction and prescribed PA spin rate
to update the degree-2 inertia/field coefficients. `moon_love`, `moon_lag_days`
and `moon_reference_spin_rad_day` control this response. The default reference
spin is IAU2009's mean prime-meridian rate. Rotation remains a separate
prescribed approximation; mantle/core torques are not integrated.

Delayed/advanced lunar positions are **local cubic predictions from the
propagated position, velocity, acceleration and Newtonian jerk**, with lags
limited to ±0.25 day. No historical/future kernel Moon state enters the force.
This approximates the lagged source, rather than solving a delay differential
system with exact state history. Sun-raised Earth tides query the Sun/EMB at
the lagged epoch, so kernel coverage must also include those lags.

The sampled lag-source error reached 301 m at the extreme ±0.25-day controls.
That is a small error in a tide-inducing direction about 384,000 km away.
The benchmark separately measures its force sensitivity; after the real source
mass scaling, the 90-day displacement scale is about 2.9 micrometres, below the
roughly 50-micrometre integration-refinement difference. This sampled estimate
is not a global delay-model error bound.

Solar terms are separately selectable. `solar_oblateness` adds the fixed
DE440 J2 field; `solar_frame_dragging` adds the spin term when `relativity`
is also enabled. `solar_radiation_pressure` uses an inverse-square spherical
approximation, with `earth_srp_fraction=2e-14` and `moon_srp_fraction=1.44e-13`.
Its Sun-to-body vector convention makes the acceleration repulsive. The
benchmark checks this direction, the polar J2/Lense–Thirring magnitudes and
frame-dragging velocity symmetry/work before the ephemeris comparisons.

The reusable `model::IntegratedOrbit<Force,Origin,Frame>` accepts a force
provider returning acceleration in km/TDB day² from `(Instant<TDB>, typedState)`.
It propagates forward/backward from the same initial state using Dormand–Prince
5(4), reuses the last stage between accepted steps, rejects excessive local
error and stops on non-finite states/forces, step underflow or budget exhaustion.
It stores quintic Hermite position segments from endpoint positions, velocities
and accelerations. Query velocity is the analytic derivative of that same
polynomial. Const queries use binary search and make no force evaluations.
The immutable trajectory can be shared by render and calendar workers.
Other bodies and different force or integration algorithms can supply the same
typed provider contract; numerical lunar assumptions are not built into the
generic integrator or the analytical evaluator.

## Kernel reader

`model::spk::Kernel` reads IEEE little/big-endian DAF/SPK type 2 and 3 segments
in frame 1 (ICRF, labeled J2000 in SPICE). Files are loaded once into immutable
shared memory. `state(target,center,t)` uses NAIF IDs, returns geometric relative
ICRF states in km and km/day, and never applies light time or aberration.
Its origin type is `origin::Relative`; runtime target/center IDs define the
origin. `segments()` provides names, frames, types and TDB coverage.
`states(targets,center,t,output)` evaluates a batch at one epoch, reusing shared
center segments with bounded local scratch storage. It preserves individual
call arithmetic and supports concurrent readers.

Center paths meet at their common ancestor before summation, avoiding a large
barycentric subtraction for Moon/Earth states. The last applicable segment
within the file has priority. Unsupported selected segment types/frames, bad
addresses/directories, cycles, malformed FTP signatures and missing coverage
raise errors. This reader does not implement the entire SPICE toolkit, multiple
kernel priority pools, UTC conversion or non-inertial frame kernels.

## Reproduce accuracy and performance measurements

```bash
cmake --preset moon "-DSPP_MOON_REFERENCE_KERNEL=$PWD/context/Spp Moon/reference-data/de440s.bsp"
cmake --build --preset moon --target spp_moon_report
```

The optional kernel enables independent CSPICE reader checks and comparisons
of propagated Moon states against DE440s. Without it the report still runs
the generic integrator's analytic controls; no download occurs. Standalone:

```bash
cmake --build --preset moon --target spp_spk_bench spp_integrator_bench
./build/moon-Linux/spp_spk_bench "context/Spp Moon/reference-data/de440s.bsp" tools/reference/spk_de440s.csv
./build/moon-Linux/spp_integrator_bench "context/Spp Moon/reference-data/de440s.bsp"
```

On the NAS (GCC 14.2 Release), seven seeds from 1900–2100, each sampled at
257 points over ±16 days, differed from DE440s by at most **0.041 m** and
0.000008 km/day. Each 32-day trajectory took about **23–34 ms** to construct
(1538 accepted steps, 9230 force evaluations). Half maximum step plus tighter
tolerances changed sampled positions by at most **0.000042 m**. These are
short-arc samples; longer propagation and dates outside those samples do not
inherit that bound. Seed/reference correlation also means this comparison
measures short-arc propagation from an accurate seed, not an independent
solution for the initial lunar orbit.

The GUI's larger **±90-day** 2026 window was sampled at 1441 instants. Maximum
position residual was **0.020164 m**. Three alternating before/after runs of
the same physical model measured median construction **467.028 → 138.703 ms**
(3.37× faster, 70.3% less time), including preparing the environment cache.
The trajectory allocates 2,752,512 bytes, plus **681,896 bytes** for the
environment and the unchanged shared 32,726,016-byte kernel. Resident observer
coordinates took about **3 µs**, complete observations **6–7 µs**, with no
force evaluations during queries. This is slower than a resident ELP lookup
but covers the entire propagated window without lunar-series cache misses.
The compared observer directions differed from direct kernel states by at most
0.0000021 arcsecond in the separate ±16-day 2026 arc.

The environment cache differed from direct inputs by at most **0.00004173 m**
in 2881 trajectory comparisons over that 180-day window, and by at most
0.00000230 m in seven 32-day windows. Separate force comparisons use 3619
dates/trial states at seven epochs, including interval boundaries and partial
segments: maximum 2.40e-11 km/day². The cache/direct trajectory gate is 0.1 mm;
these sampled errors do not establish a global interpolation bound.

Ablation uses the same seed, 180-day window, integration settings and 721
reference instants. On 2026-11-05, three sequential NAS runs measured these
median construction times with the optimized environment:

| Force configuration | Maximum position residual | Trajectory construction |
| --- | --- | --- |
| Earlier quadrupole family | 0.707542 m | 82.370 ms |
| Higher-order gravity | 0.599136 m | 94.583 ms |
| Add Earth tides | 0.292405 m | 108.346 ms |
| Add lunar tides/spin deformation | 0.186463 m | 110.251 ms |
| Add external figure reactions | 0.023888 m | 126.234 ms |
| Add solar J2 | 0.022226 m | 127.295 ms |
| Add photon pressure | 0.020497 m | 127.685 ms |
| Add solar frame dragging | 0.020164 m | 129.243 ms |

The tool prints the full sequence on each run. The more complete forces cost
time while constructing the trajectory; queries still use the same dense
output and make no force calls. Both paths make 51,854 force calls, with the
same tolerances and no rejections; the improvement reduces cost per call.
The final ablation row is warmer than the preceding standalone GUI-window
construction. These are workload-specific timings, rather than per-call
latency percentiles. Before/after logs are in
`context/Spp Moon/performance/`: the files
`2026-11-05-performance-before-1.txt` and
`2026-11-05-performance-after-1.txt`, followed by pairs 2 and 3. The machine-readable
summary is `2026-11-05-performance-summary.json` in that directory.

Independent matched-TT phase checks also cover all 62 pinned DE441 short arcs.
The integrated model's sampled maximum was **0.268178 s**, RMS **0.128456 s**.
Four ±190-day trajectories took a median **1.178 s**, down from **3.918 s**
before this optimization; finding the 62 events then took about 1.2 ms.
The earlier direct ELP comparison was 0.211037 s maximum. Phase
predictions include the shared solar/apparent-place pipeline, so improved lunar
orbital agreement does not imply an equally improved phase-time comparison.

Its optimized grid/daily OpenGL benchmark crossed midnight with 120 one-minute frames:
median 2.197 ms, mean 2.340 ms, p95 4.162 ms, maximum 5.632 ms, including one
daily schedule refresh. This was Mesa llvmpipe software rendering, not a Windows
GPU measurement. The same run reported **149.882 ms** of initial propagation,
including kernel loading, compared with the preceding run's 473.357 ms.
Rendering itself is unchanged; these separate runs do not establish a frame-time
improvement. The [application guide](moon-app.md) explains benchmark timing
boundaries and the separate background scheduling used interactively.

An analytic circular-orbit control over ±32 days had maximum position error
0.000046 m. Random stored-state queries measured about **0.08–0.10 µs** without
force calls. Independent CSPICE controls cover 180 states at endpoints,
record boundaries and fractional epochs. Maximum Moon/Earth reader difference
was 6.51e-8 m; all-body maximum 0.000955 m, dominated by very large outer-planet
vectors. Kernel loading took about 8 ms for 32,726,016 bytes; state queries
about 0.45 µs. Synthetic kernels cover both endians, both polynomial layouts,
derivatives, center chains and corrupt inputs.

## Provenance

- [DE440s kernel](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp), SHA256 `c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2`, coverage 1849-12-26 through 2150-01-22 TDB.
- [DE440 gravitational parameters](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc), SHA256 `924ddf4fb9ead9fe8a1aa55780bcabde40b09d00065d58226e24b68d8092f140`. Figure constants are in the DE440s DAF comments.
- [NAIF lunar frame definition](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_250416.tf): fixed ME/PA offset; IAU2009 rotation supplies its approximate time-dependent orientation.
- [Park et al. (2021), DE440/DE441](https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf), especially equations 27–29 for major-body relativity and figure/reaction forces. Spp uses a subset of that paper's model.
- [Folkner et al. (2014), DE430/DE431](https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf), equations 13 and 31–32 for tidal lag/rotation conventions; Park2021 equations 30–34 and 54 for the lunar inertia-to-gravity response.
- [LLNL SSAPy solar radiation model](https://software.llnl.gov/SSAPy/api/ssapy.accel.AccelSolRad.html): repulsive solar-line direction convention. The lunar/terrestrial fractions and solar spin/J2 parameters come from Park2021 and the pinned DE440s comments.
- [NAIF DAF](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html) and [SPK](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html) format specifications.

`tools/spk_reference.py` generates the checked-in CSV using independent CSPICE
via SpiceyPy; its header records tool versions and kernel hash. That Python
dependency is needed only for regeneration, not Spp builds or runtime.

`tools/gravity_reference.py` independently evaluates Legendre potentials with
70-digit Decimal arithmetic and finite differences. Its checked-in fixtures
cover 28 harmonic gradients (including poles) and 18 separate tidal-mode
gradients. The benchmark also checks the isotropic tidal limit, quadrupole and
zonal fast paths, numerical refinement, short-lag sensitivity and the ablation
above. Regeneration needs only Python's standard library:

```bash
python3 tools/gravity_reference.py "context/Spp Moon/reference-data/de440s.bsp"
```
