# Moon speed and accuracy report

The performance work is built and checked against direct ELP and the independent
reference data. The same commands let the user reproduce the accuracy gates and
timings. The [OpenGL app](moon-app.md) has a separate rendering benchmark.

## Run the complete report

For a NAS/Windows CPU comparison, run this identical command from the project
root in either Bash or PowerShell; no shell variables or previous session are
needed:

```text
cmake -P cmake/BenchmarkMoon.cmake
```

It configures portable Release settings (`SPP_NATIVE=OFF`, CAPI disabled),
builds the report executables and runs all ten stages with the existing DE440s
kernel. Windows builds under `%LOCALAPPDATA%/SPP/moon-build`; Linux uses
`build/moon-Linux`. Timestamped UTF-8 logs, including host/processor information,
are saved under `context/Spp Moon/performance/comparison/`. Compare matching
stage-8 median query times and stage-10 construction times; NAS time divided
by Windows time gives the Windows speedup. Kernel loading includes storage and
network costs. The command does not run the GUI or download reference data.

From the project root:

```bash
cmake --preset moon
cmake --build --preset moon --target spp_moon_report
```

In VS Code, **Tasks: Run Task → Moon: Accuracy and Speed (Moon Preset)** does the
same configure/build/report sequence. The existing Ctrl+Shift+B task is retained.

The report target builds current sources before executing the checks. It runs
nine stages (ten with an optional numerical reference kernel) and returns failure if any stage fails. It still attempts the
other stages, including the speed report, after an accuracy-stage failure.
The report runs only when explicitly requested. Existing CTest tests remain
separately available.
Reference data is already in the repository; no downloads are performed.
The numerical orbital comparisons additionally require an existing DE440s
kernel. Enable them with the configuration command in
[lunar-integration.md](lunar-integration.md). Without that option, the generic
integrator's analytic numerical/performance controls still run.

## Accuracy output

The report prints maximum absolute error, RMS error, mean residual, units,
sample count, the worst-case date and the acceptance limit for each metric.

| Reference | Coverage | What is checked |
| --- | --- | --- |
| Pinned upstream ELP/MPP02 evaluator | 144 states; both fits; computational endpoints and seeded interior dates | Cartesian position-only and full-state positions, analytic velocity, and agreement between evaluator paths |
| Direct full ELP/MPP02 | Both fits; seeded intervals throughout the computational domain, denser 2000-2100 coverage, segment edges and lunar observation cases | Error introduced by interpolation, cached velocity, observer/phase agreement, and cache contracts |
| JPL Horizons DE441 | 15 topocentric and 5 geocentric cases; 2000, 2024, 2026 | Sky/equatorial directions, range, angular size, elongation, physical illumination/phase angle and apparent bright-limb bearing |
| USNO | 198 principal phases in 2000, 2024, 2026 and 2050 | Event count/order, date/time residuals, and separate longitude-root residuals |
| JPL Horizons DE441, matched TT | 62 principal phases: all 50 in 2026 and four each in 2000/2024/2050; four independent short-arc samples per event | Physical phase-time residual without UT/UTC/Delta-T mismatch, cached/direct timing agreement, and search cost |
| Adaptive horizon-search benchmark | Both fits; Paris, Sydney, equator, Arctic, Antarctic and both poles; synthetic close crossings, stationary contacts and interval boundaries | Counts/directions, one-minute dense comparison, direct/cache agreement, visibility duration and timings |
| Earth-cache benchmark | 256 resident segments, old modulo collisions and repeated eviction/index rebuilds | No unnecessary resident refill, bitwise state stability and memory/lookup cost |
| Lunar rotation | Eight independent evaluations of NAIF's IAU/PCK series; seven independent DE441 orientation cases, 1900-2100 | Axis/pole agreement, analytic rotation rates, state-velocity transforms, inverse/composition and speed; see [lunar-rotation.md](lunar-rotation.md) |
| Surface/libration | Thirteen independent geocentric/topocentric DE441 cases; fixed reference rotations isolate geometry | Sub-observer/subsolar coordinates, pole bearing, view matrix and replaceable rotation provider |
| Numerical integration | Analytic circular orbit over ±32 days; adaptive rejection, backwards propagation, dense derivatives and step budget | Numerical errors, force counts, trajectory construction and query cost |
| Optional DE440s kernel | 180 independent CSPICE states, synthetic type-2/3 kernels in both endians; seven ±16-day lunar arcs, 1900–2100 | Kernel interpolation, malformed inputs, short-arc propagation residuals and numerical refinement |
| Extended numerical forces | 28 independent Decimal harmonic gradients, 18 tidal gradients, 2026 ±90-day ablation, short-lag controls and 62 matched-TT phase arcs | Higher gravity, response modes, external reactions, orbital/phase residuals and construction/query cost |
| Numerical environment cache | 3619 cached/direct force comparisons at seven epochs; seven ±16-day trajectories and a 2026 ±90-day trajectory; endpoints, partial segments, out-of-window and missing-coverage cases | Force residual below 1e-10 km/day²; cached/direct trajectory difference below 0.1 mm; construction time, sample counts and auxiliary memory |
| GUI workers (separate CTest/benchmark) | 2048 idle shutdowns, 128 active cancellations, 128 coalescing/error cycles and a superseded uninterruptible result | Actual UI-worker lifecycle, cancellation and latest-result publication; `spp_moon_worker_bench` needs no graphics context |

The JPL/USNO end-to-end comparisons use the default DE405/DE406 fit. Upstream
state compatibility and the speed measurements cover both fitted sets.
Default observation/phase calls now use the interpolated provider, so the
independent JPL/USNO checks exercise the production cached path as well.

The separate `spp_lunar_cache_check` stage requires sampled interpolation
position errors below **1 mm during 2000-2100** and **0.1 m across the full
computational domain**, per Cartesian component. Velocity limits are 0.0001
and 0.02 km/day respectively. Distant-date rounding in direct ELP's time and
argument polynomials affects the comparison too; it is reported separately.
Cached-vs-direct observations must agree within 0.00001 arcsecond in sky
direction/cycle, 1 mm in range, and 0.02 seconds for the tested phase events.
These are regression budgets, not global interpolation-error guarantees.
The direct upstream limits are unchanged. Measured sampled errors are recorded
below; interpolation/disk checks are kept separate from physical model residuals.

Synthetic cubic/derivative checks also cover normalized splits, partial first
and last intervals, cache eviction, shared endpoint reuse, failed fills,
fit isolation, invalid dates, and light-time queries without new source calls.

ELP agreement verifies implementation of the same theory; it is not an
independent measurement of lunar orbit accuracy. JPL sky cases use airless
coordinates and each reference's DUT1. Earth-orientation and phase conventions
differ slightly from Spp. Illumination error is printed in **percentage points**.
USNO event references have minute precision and report UT. Small numerical
root residuals must not be interpreted as equally small physical timing errors.
Maximum sampled errors are not guaranteed limits for untested dates/locations.

The illumination checker now uses Horizons **quantity 43 (`phi`)**, not the
approximate quantity 24 (S-T-O). It compares both the geocentric phase and
the topocentric disk: 20 illumination cases in total. Correcting incoming-light
aberration and using the astrometric outgoing ray reduced the sampled worst
illumination error to **0.0000195 percentage points**. The old apparent-triangle
implementation had 0.002205 percentage points on its five geocentric cases.
Physical angle residuals are at most about 0.204 arcsecond against `phi`, whose
reference resolution is 0.0001 degree. Bright-limb reference bearings have
only 0.01-degree resolution; the sampled worst difference is about 54.4
arcseconds and does not establish subarcsecond orientation accuracy.

The matched-TT phase benchmark uses separate Moon and Sun quantity-31 tables
from Horizons, geocentric and airless. It finds each reference phase with a
quadratic short-arc fit and checks it against a fourth unused sample. USNO
only locates the neighborhood; its rounded times are not the accuracy reference.
Spp's result is compared in TT, including for 2050, so different predictions
of future leap seconds or Delta-T do not appear as phase-model errors. The
largest sampled direct-vs-DE441 difference was **0.212 seconds**, RMS about
**0.082 seconds**. Horizons includes relativistic/frame conventions absent
from Spp. This is independent sample evidence, not a global subsecond promise.

```bash
cmake --build --preset moon --target spp_moon_phase_bench spp_moon_event_bench spp_earth_cache_bench
./build/moon-Linux/spp_moon_phase_bench tools/reference/moon_horizons_phase_tt.csv
./build/moon-Linux/spp_moon_event_bench
./build/moon-Linux/spp_earth_cache_bench
```

To see calculated and reference values for each sky case and each phase event:

```bash
./build/release/spp_moon_sky_check tools/reference/moon_horizons.csv tools/reference/moon_usno_phases.csv --details
```

NaN/infinite residuals fail explicitly. Accuracy results are never inferred
from elapsed benchmark time, and a speed improvement cannot hide an accuracy
failure. The model/provider/refraction/interval checks also remain enabled.

## Speed output

The Moon-only benchmark measures both DE405/DE406 and LLR parameter fits, using
nearby one-minute instants and random dates during 2000-2026 with seed 7.
It retains the preceding direct full-series timings and adds interpolated
timings for the same inputs and corrections:

- Full geometric position and velocity.
- Position-only evaluation.
- Observer altitude/azimuth.
- Complete observation, including geocentric phase.
- Phase/illumination.
- A 35-day principal-phase search, with returned event count.

The cached provider owns 256 twenty-minute cubic segments, with hash lookup
and LRU replacement across all slots. A new isolated interval evaluates two full ELP states; an adjacent
interval reuses its common endpoint. Subsequent queries within a resident
interval evaluate only a polynomial. Both parameter fits have separate instances.

**Repeat** rows warm up the workload before timing; the cache has bounded
capacity, so a large or unlucky random workload can still miss. The printed
**ELP calls/query** value is from the last measured batch and exposes these
misses. **Forced miss** rows clear the lunar segments and endpoint memo before
every query, including that clearing cost. The position-only miss can be slower
than direct position evaluation because it constructs two full endpoint states;
coordinate/observation queries can reuse them across their light-time calls.
Near a segment boundary those calls can require a second segment.

Each point-query case has one untimed warmup, followed by five timed batches of
128 queries. The interval case performs one 35-day search per batch.
Input generation and input-TDB conversion are outside the timers; observation
functions' own required time conversions remain part of their cost. The
benchmark includes every returned numeric field in a checksum consumed after
timing, so optimization cannot silently turn the full-state/observation cases
into position/altitude-only calculations. File reads and printing are outside
the measured loops.

Best/median/worst values are **batch-average microseconds per query**, not
individual-query latency percentiles. Throughput is derived from the median.
Search timing is milliseconds per 35-day request, reported for direct ELP,
a repeated cached workload, and an empty lunar cache at the start of each
search. The last case clears outside the timer and includes all subsequent
segment construction. Search source-evaluation counts are printed too.
Earth's existing cache and the CPU/OS caches are not flushed in these cases;
an empty lunar cache is not a claim of a cold machine. Large random workloads
can also refill Earth's bounded cache, independently of the lunar model.
Compiler and CMake build configuration are printed. Use Release for performance
comparisons and retain the hardware/workload details when comparing machines.
There is no arbitrary speed pass/fail threshold.

The previous direct report measured roughly 3.9-4.8 ms per coordinate query
and 6.8-7.7 ms per full observation on the user's NAS. A twenty-minute cache
run measured approximately 2.6 microseconds for resident coordinate queries
and 2.5 ms for forced lunar misses. A 35-day principal-phase search measured
about 0.17 ms repeated or 102 ms starting with an empty lunar cache, versus
285 ms direct. These results preceded the cache-index improvement described
below; they are not universal latency guarantees.

The first five-minute implementation failed its 2000-2100 velocity budget
(0.000351 km/day versus 0.0001). Twenty-minute intervals reduce the sensitivity
of polynomial derivatives to rounding in endpoint positions while retaining
the positional budget. The recorded difficult dates remain explicit regression
cases. Sample maxima after this adjustment were 0.807 mm position and
0.0000746 km/day velocity in 2000-2100; full-domain position was 0.028 m.

The benchmark also exposed unnecessary misses in the original four-way
index: the 128-date repeat workload performed 0.312 direct evaluations/query
despite fitting within total capacity. The index now uses compact open-addressed
lookup and LRU replacement across all 256 slots. A synthetic 256-interval
regression checks that hash collisions cannot trigger these avoidable refills.

After the index change, an isolated GCC 14.2 Release DE405 run (128 queries,
three measured batches) produced these median batch means on the NAS:

| Coordinate workload | Direct | Cached repeat | Forced lunar miss |
| --- | --- | --- | --- |
| Nearby one-minute dates | 3788 us | 2.386 us | 2334 us |
| Random 2000-2026 dates | 4609 us | 798 us | 3126 us |

Both repeat cases performed **zero** direct lunar source evaluations. The
remaining random-date cost includes refilling Earth's existing 64-slot cache;
the lunar cache alone cannot eliminate that work. The repeated 35-day search
was 0.165 ms, empty-lunar-cache search 97.179 ms, and direct search 271.855 ms.
The collision/reindex/eviction regression and numerical interpolation/disk
checks passed after the index change. These workloads are not machine-wide
cold-cache measurements and do not establish global accuracy bounds.

The subsequent Earth-cache change keeps the same fit, interval length,
coefficients and Clenshaw arithmetic, and replaces the old direct-mapped
64 slots with 256 slots and global LRU/hash lookup. On this build it uses
**273,496 bytes per thread**, versus about 66 KiB for the old storage. Resident
Earth queries, including bitwise verification in the benchmark, took about
0.115 microseconds with zero refills across all 256 slots, including dates
deliberately colliding under the old modulo index.

A fresh DE405 comparison (128 dates, three batches) measured repeated random
coordinates at **819.412 us before / 2.606 us after**, with zero lunar source
evaluations in both cases. Nearby repeated coordinates measured 2.766 / 2.518
us. These are separate runs and include normal machine variation; the random
speedup is explained by eliminating avoidable Earth rebuilds. Working sets
larger than either cache can still miss; genuinely cold full-series fills
remain millisecond operations.

The event benchmark's 14 ordinary daily workloads measured about 0.08-0.16 ms
for a warm search and 44-69 ms with an empty lunar cache. Adaptive and dense
times differed by at most 0.004600 s. Grazing, same-sign-endpoint visibility
windows, stationary boundary visibility and half-open intervals have separate
gates. Neither those timings nor self-comparisons establish physical
moonrise accuracy under real atmospheric conditions.

For speed alone, or a larger measurement:

```bash
./build/release/spp_bench --moon-only
./build/release/spp_bench --moon-only --samples=512 --rounds=9 --fit=de405
```

`--samples` accepts 16-8192; `--rounds` accepts 3-31; `--fit` accepts `both`,
`de405` or `llr`. Without `--moon-only`, the existing solar benchmarks also run.
