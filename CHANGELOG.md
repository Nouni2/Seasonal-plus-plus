# Spp 1.0.1 — changelog and implementation notes

| | |
|---|---|
| Base | `159c970` (Spp 1.0.0, GitHub `Nouni2/Seasonal-plus-plus`) |
| Date | 2026-10-27 |
| Version in `CMakeLists.txt` | still `1.0.0`; this changelog describes the unreleased 1.0.1 |
| Built and checked with | GCC 14.2 (Release, Debug, UBSan), Clang 18 + libstdc++ 14, GCC 14.2 cross to aarch64 Linux; Linux x86-64 |
| Audit | [audit.md](audit.md) |

The README's own changelog (section *Changelog → 1.0.1*) is the short, user-facing version. This file
gives, for every change, the problem as measured, what was changed and where, how it was verified,
and its effect on results, classified as the Constitution requires (Article II.4):
**[B]** bit-identical, **[R]** results change.

---

## Summary

| Issue | Change | Class | Verified by |
|---|---|---|---|
| #1 BUG-09 | Sunrise/sunset solve the refracted upper limb = 0 | [R] events | `reference_check.py` (DE421), `spp_pipeline_check`, `spp_day_scan` |
| #2 BUG-02 | Monthly observed ΔT from IERS EOP 20 C04, 1972-01 to 2026-09 | [R] estimated DUT1, TT after 2027-06-28 | `reference_check.py` (IERS) |
| #3 BUG-01 | Every `Instant` normalized at the API boundary | [B] for normalized input | `reference_check.py` (2,000 splits), 66,000-request bit comparison |
| #4 TST-01b | Correctness checks against independent references; event specification | tools, docs | — |
| #5 TST-01a | Extended-precision reference of the whole pipeline | tool | — |
| #7 BUG-03 | Day type at the poles near the equinoxes | [R] ±0.1° of the poles only | `spp_day_scan` (88,330 days) |
| #8 BUG-05 | Invalid dates and non-finite input rejected | behaviour on invalid input | UBSan run, CLI |
| #9 BUG-04 | Angle wraps stay in [0, 360) | [B] except at the wrap point | boundary probe |
| #10 BUG-06 | ISO 8601 expanded years | output format | CLI, probe |
| #11 BUG-07 | Held DUT1 taken at the table boundary; TT continuous | [R] after 2027-06-28 (0.3 ms) | `reference_check.py` |
| #12 BUG-08 | ΔT joins documented and measured | docs | `reference_check.py` |
| #26 API-03 | `crossing` precondition documented and asserted | debug builds only | debug run |
| #27 API-04 | One position pipeline (`src/observe.cpp`) | [B] | 66,000-request bit comparison |
| #28 API-05 | DUT1 clamp shown inactive; leap-second continuity checked | docs, checks | `reference_check.py` |
| #29 BLD-02 | Cross-compiling builds or emulates `vsop_gen` | build | aarch64 cross build |
| #30 BLD-03 | `-fexcess-precision` role stated (x87 only) | docs | object-code comparison |
| #31 BLD-04 | `libstdc++exp` chosen by a link test | build | GCC and Clang builds |
| #32 CLI-01 | `spp` rejects `inf`/`nan` | CLI | CLI |
| #33 CLM-01 | Warm and cold performance published with context | docs, `spp_bench` | `spp_bench` |
| #34 CLM-02 | Floating-point reproducibility claim corrected and measured | docs, check | `reference_check.py`, aarch64 run |
| #35–#41 CLM-03…10 | Every README figure has a reproducing tool; wrong figures corrected | docs, tools | `spp_readme_numbers`, `spp_cache_check`, `reference_check.py` |

Not started: #13–#25 (milestones 1.1 Free speed, 1.2 Sub-0.1 µs, 1.3 Events and batch). See
[handoff.md](handoff.md).

---

## Correctness fixes

### #1 BUG-09 — sunrise and sunset solve the documented condition [R]

**Problem.** The README defines sunrise and sunset as the upper limb touching the horizon. 1.0.0
refracted the Sun's *center* and solved for "refracted center altitude = −semidiameter". Refraction
grows by 2.7′ over the 16′ between the upper limb and the center near the horizon, so at 1.0.0's
instants the refracted upper limb was 161″–167″ below the horizon (sea level). Measured against the
documented definition solved with JPL DE421 geometry, 1.0.0's sunrise was early (and sunset late) by
8.5–9.5 s at Quito (2,850 m), 13–15 s at Singapore, 20–25 s in Paris, 30–90 s at Reykjavík and up to
18 minutes on days next to polar day. The day type disagreed with the definition on 12 days at
Tromsø and Longyearbyen over 2000–2026.

**Change** (`src/events.cpp`, `include/spp/events.hpp`):
- `upper_limb_altitude_from(g, obs, atm)`: geometric topocentric altitude of the center (the same
  pipeline with refraction off) plus the semidiameter 959.63″/R *at that instant*, then refracted
  with the atmosphere at the observer's elevation.
- `crossing_impl` takes a `Point` (`Center` or `UpperLimb`). For the center the objective is
  computed exactly as before, so `crossing`, twilights, transit and `altitude_deg` are bit-identical.
  For the upper limb the analytic first guess subtracts the semidiameter at transit.
- `day()` classifies and solves on the upper-limb condition with target 0.
- New public functions: `events::upper_limb_altitude_deg`, `events::upper_limb_crossing`. The
  constant `upper_limb_deg` is kept; its comment now says it is the geometric center altitude for use
  with refraction off.
- `tools/figures.cpp` draws the daylight band from `Day::rise`/`Day::set` instead of a center crossing
  at −16′. All figures were regenerated.

**Verified.**
- Spp's instants against the definition solved with DE421 geometry and the same refraction model
  (`reference_check.py`): within 11 ms at 0°–64°, 0.44 s at 69°–78° (grazing days).
- The condition in Spp's own pipeline: |refracted upper limb| ≤ 6.4 × 10⁻⁶″ (acceptance: 10⁻⁵″).
- Day type: 0 disagreements at 8 sites, 2000–2026.
- Extended precision (`spp_pipeline_check`): sunrise/sunset within 0.74 µs (|φ| ≤ 60°).
- Center crossings and transits: bit-identical to 1.0.0 on 3,000 random days with DUT1 given.

**Effect.** Sunrise later and sunset earlier by the figures above. Offsets from the −0.8333° almanac
convention shrink from 15–48 s to 1.3–5.2 s. README examples, tables, gallery numbers and figures were
re-measured (`spp_readme_numbers`).

### #2 BUG-02 — observed ΔT from the primary source [R]

**Problem.** The annual 2005–2026 table was wrong by up to 0.075 s in 2023–2025. Measured against
IERS EOP 20 C04 on every day since 1972, the estimated DUT1 was off by up to 0.12 s (1.8″ of hour
angle): annual sampling misses the ±0.03 s seasonal variation, and the 1961–2005 polynomials are off
by up to 0.12 s themselves. 111 of 261 months since 2005 failed the issue's 20 ms criterion.

**Change.**
- `tools/delta_t_table.py` generates `src/model/delta_t_observed.hpp` from `eopc04.1962-now` and
  `Leap_Second.dat`: ΔT = 32.184 s + (TAI − UTC) − (UT1 − UTC) on the 1st of every month at 0h UTC,
  1972-01 to 2026-09, 657 values with the source's 7 decimals, keyed by MJD.
- `src/model/espenak_meeus.cpp` interpolates the table linearly in MJD. The year→MJD conversion can
  round just outside the table at its ends; the lookup clamps to the end intervals. The 1986–2005
  polynomial is no longer used; the 1961 piece ends at 1972-01-01. The bridge to 2050 starts from the
  last entry (69.1815 s at 2026-09-01).
- README §2.5/§2.6 document the source, the derivation and the refresh command.

**Verified.** Estimated DUT1 equals IERS on the 1st of every month (1.5 × 10⁻¹³ s) and is within
6.7 ms on every day from 1972. 1 January values match USNO's figures quoted in #2 (2024: 69.1752 s).

**Effect.** TT is unchanged from 1972 to 2027-06-28 (exact from leap seconds). Hour angles, altitudes,
azimuths and events change by up to 1.8″ or 0.12 s where DUT1 is estimated. TT after 2027-06-28 is
earlier than in 1.0.0 by 0.002–0.29 s up to 2050 and 0.01 s after (0.012″ of solar longitude at most).

### #3 BUG-01 — `Instant` splits [B for normalized input]

**Problem.** `Instant` is an aggregate; a hand-built `{day, frac}` with frac outside [0, 1) picked the
wrong 16-day Chebyshev segment and extrapolated its fit: up to 179° in right ascension.

**Change.** `src/core/split.hpp` (new) holds `detail::normalize` (moved from `time.cpp`) and
`detail::normalized(Instant)`. Every public function taking an `Instant` normalizes first:
`earth::nutation`, `mean_obliquity_deg`, `gmst_deg` (and through it `gast_deg`, `lst_deg`),
`sun::geometric` (and through it `apparent_ecliptic`, `apparent`), `sun::semidiameter_deg`,
`tai_minus_utc`, and `detail::earth_state`. A fast path returns at once when frac is already in
[0, 1), so the cost is one comparison; without it `std::floor` (a library call on baseline x86-64)
made `sun::apparent` 15 % slower. Non-finite or astronomically large fractions are left as they are
(results NaN, no undefined behaviour, see #8).

**Verified.** 2,000 random splits with frac in [−100, 100] give results identical to their
normalized split (1.0.0: 1,823 differ). 66,000 requests (positions, Sun states, altitudes, days)
are bit-identical to the build before the change.

### #7 BUG-03 — day type at the poles [R, within about 0.1° of a pole]

**Problem.** `PolarNight` was declared when the limb was below the horizon at transit. At ±90° near
an equinox the change of declination over 12 h exceeds the daily swing of altitude, so the maximum can
fall at an edge of the window: 2026-09-25 at 90° N and 2026-09-20 at 90° S were wrongly polar night.

**Change** (`src/events.cpp`). The edge altitudes are now computed for every day. If the limb is below
the horizon at transit but above it at an edge, the day is `Normal`, with a sunset between the edge and
transit (Sun up at the start) or a sunrise between transit and the edge (Sun up at the end), refined
by the same Illinois method (extracted into `refine()`, bit-identical for existing paths). Duration is
the time above the horizon within the window.

**Verified.** `spp_day_scan` (new): every day of 2026 at 60°–90° in 0.25° steps, both hemispheres,
88,330 days, against the condition sampled every 2 minutes: 0 inconsistencies. The two dates now give
a sunset at 03:43:08 UTC (90° N) and a sunrise at 20:25:07 UTC (90° S).

### #8 BUG-05 — validation [behaviour on invalid input]

**Problem.** Invalid dates rolled over silently (2026-02-30 → 2026-03-02, month 13, the 1582 gap) and
NaN reached an integer conversion (undefined behaviour).

**Change.**
- `is_valid(CivilDate)` (new, `calendar.hpp/.cpp`): month, day within the month (Julian leap rule
  before 1582-10-15, Gregorian after), not 1582-10-05…14.
- `Time::from_civil/from_utc` throw `std::invalid_argument` for invalid dates and non-finite seconds;
  hours, minutes and seconds still carry over (the library relies on 24:00 for longitude −180°).
- The `Time` constructor throws for a non-finite DUT1, a fraction with |frac| ≥ 2⁵², or a day beyond
  ±2⁶²; `shifted` throws for a non-finite or huge shift.
- `events::day`, `transit`, `crossing`, `upper_limb_crossing` throw for non-finite observer coordinates
  or targets (the approximate-noon conversion cast the longitude to `int`).
- `detail::normalize` leaves non-finite or huge fractions untouched (found by UBSan: a hand-built
  `Instant` with an infinite fraction overflowed the day).
- `spp` reports these errors and exits with status 2; `spp_bench` draws only valid dates.

**Verified.** A probe with 26 invalid inputs built with `-fsanitize=undefined
-fno-sanitize-recover=all` (library and probe): each throws or returns NaN, no sanitizer report.

### #9 BUG-04 — angle wraps [B except exactly at the wrap]

`norm360` returns 0 when the wrap rounds to 360, +0 for −0, and keeps NaN. `src/sky.cpp` uses one
`wrap_two_pi` for right ascension, azimuth and the parallax wrap. A probe of −1e-20, −1e-14, −0.0, 360,
720 − 1e-13 and others returns values in [0, 360) (1.0.0 returned 360 three times).

### #10 BUG-06 — ISO 8601 years [output format]

`iso_date(CivilDate)` (new) writes years outside 0…9999 with a sign and at least four digits
(`-0005-01-01`, `-0100-03-01`); used by `Time::iso_utc`, `spp year` and `spp_spiral`.

### #11 BUG-07 — continuous TT at the end of the leap-second table [R after 2027-06-28]

The held DUT1 was evaluated 1.5 days before the boundary, so TT − UTC stepped by 0.506 ms. It is now
the last exact TT − UTC minus ΔT at the boundary instant (2027-06-28 0h UTC): TT − UTC is equal to all
17 digits on both sides. Held value −0.0271 s.

### #12 BUG-08 — ΔT joins [docs]

The comment claiming continuity at 1600, 500 and −500 is corrected. The published pieces are kept as
published (Article V forbids undocumented corrections); the jumps are measured
(`reference_check.py`, section `delta_t_pieces`) and listed in README §2.6: up to −0.251 s at 1600,
−0.021 s at the new 1972 table start, 0 at 2026-09, 2050 and 2150. A jump changes TT only: 0.01″ of
solar longitude, about 0.7 ms on an event time at the equator.

---

## API, build and tools

### #26 API-03 — `crossing` precondition

Documented in `events.hpp` and the README; `crossing` and `upper_limb_crossing` assert
|hour angle| < 0.01° (2.4 s) in debug builds (checked: aborts for transit + 0.1 day, silent for real
transits; the debug-build `spp_day_scan` passes).

### #27 API-04 — one position pipeline [B]

`src/observe.hpp/.cpp` (new): `detail::Geocentric`, `detail::geocentric(Time)`, `detail::lst_deg`,
`detail::observed`. `sun::horizontal` and the event solvers both use them. 66,000 requests are
bit-identical to the build before; `spp_bench` before/after on the same machine shows no change
(0.18–0.24 µs per `sun::apparent` both ways, run-to-run noise).

### #28 API-05 — DUT1 clamp

With the monthly table, the largest estimated |DUT1| from 1972 is 0.81 s, so the ±0.9 s clamp never
engages. New checks across all 27 leap seconds: ΔT continuous to 4 × 10⁻¹¹ s, and UT1 advances by the
inserted second plus the UTC interval to 4 × 10⁻¹¹ s. README §2.5 states where the clamp would act.

### #29 BLD-02 — cross-compiling

`CMakeLists.txt`: with `CMAKE_CROSSCOMPILING_EMULATOR`, `vsop_gen` is built for the target and CMake
runs it through the emulator; without one, it is built for the host with `SPP_HOST_CXX` (default
`c++`) or taken from `SPP_VSOP_GEN_EXECUTABLE`. Checked with `aarch64-linux-gnu-g++-14` both ways: the
generated table is byte-identical to the native one, and the aarch64 `spp` (under QEMU) prints the
same sunrise and sunset as x86-64.

### #30 BLD-03 — `-fexcess-precision=standard`

All 13 library sources compile to identical x86-64 machine code with and without the flag (GCC 14.2,
`objdump` comparison). The flag is kept for x87 targets; `SppFlags.cmake` and the README say so.

### #31 BLD-04 — `libstdc++exp`

The configure step test-links a `std::print` program, then with `stdc++exp`, and links it only if
needed, for any compiler (the checks compile as C++23; the library stays C++20). GCC 14.2 and Clang 18
build on Linux without it; the fallback path was exercised by forcing the first check to fail.

### #32 CLI-01

`num()` rejects non-finite values: `inf`, `-inf`, `infinity`, `nan`, `-nan`, `1e999` → exit 2.

---

## Published figures and claims

### #33 CLM-01, #35 CLM-03 — performance

`spp_bench` measures a cold path (random instants and dates over 1000 BC–3000 AD, timed on first use)
and prints compiler and build type. README Performance states machine, compiler and build, warm and
cold: `sun::apparent` 0.18 µs, `altitude_deg` 0.47 µs warm / 0.91 ms cold, `day` 7.9 µs / 0.97 ms,
`spp year` 29 ms. The segment fill is about 0.9 ms here, not 0.3 ms. Banner and Features say "warm
cache".

### #34 CLM-02 — floating-point reproducibility

The policy now says what the flags guarantee (no compiler contraction or excess precision) and what
they do not (libm differences). Measured: on x86-64, turning off glibc's FMA/AVX2 variants changes 320
of 20,000 positions by at most 5.7 × 10⁻¹⁴° (`reference_check.py`, `libm_dispatch`); x86-64 with the
default variants and aarch64 (QEMU) are bit-identical on the same 20,000 positions.

### #36 CLM-04 — "reentrant" → "thread-safe", with the reason.

### #37 CLM-05 — aberration

Light time plus aberration against DE421 differ from −20.4898″/R by up to 0.010″ in longitude, almost
all a monthly term of 0.0093″ amplitude (Earth about the Earth–Moon barycentre); 0.0006″ left after
fitting it; latitude below 0.001″. README error budget and §7 corrected.

### #38 CLM-06 — 1972 step: 42.250 s → 42.184 s, 0.066 s (`spp_readme_numbers`).

### #39 CLM-07 — GMST

`reference_check.py`, section `gmst_exact`, evaluates the formula in exact rational arithmetic
(Python `Fraction`, decimal coefficients exact) at 3,000 random instants: split form 1.9 × 10⁻⁷″,
direct double form 3.6 × 10⁻⁴″ on that sample.

### #40 CLM-08 — the cache offset

`spp_cache_check` rebuilds the fit without the offset: worst L error 1.9 × 10⁻⁵″ instead of
7.5 × 10⁻⁶″, ×2.6 or 1.4 bits, not "about 10 bits". Comment and README corrected.

### #41 CLM-10 — figures without a tool

- `spp_readme_numbers` (new) prints the event, sidereal-time and time-scale figures of the README.
- `spp_cache_check` measures 32-day segments: B 4.1 × 10⁻⁹″ instead of 2.9 × 10⁻¹¹″.
- Altitude residuals: `reference_check.py` and `spp_pipeline_check`.
- GMST rational check: `reference_check.py`.

---

## Validation tools (#4, #5)

| Tool | Checks | Result on the final build |
|---|---|---|
| `tools/reference_check.py` + `spp_dump` | 40 checks against IERS, ERFA, JPL DE421/DE431, Swiss Ephemeris calendar, CIO/ERA hour angle, event definitions, splits, libm dispatch | 40/40 PASS |
| `spp_pipeline_check` | Whole pipeline vs 80-bit re-implementation, 1000 BC–3000 AD | PASS: positions ≤ 1.03 × 10⁻⁵″, events ≤ 0.74 µs (≤ 3.8 µs above 60°) |
| `spp_cache_check` | Cache vs 80-bit series; fit variants | PASS |
| `spp_day_scan` | Day type and events vs 2-minute sampling, 60°–90° | PASS, 88,330 days |
| `spp_readme_numbers` | README figures | matches README |

The event specification (condition, code, verification) is the table in README §13.1.

---

## Files

New: `CHANGELOG.md`, `audit.md`, `handoff.md`, `src/core/split.hpp`, `src/observe.hpp`,
`src/observe.cpp`, `src/model/delta_t_observed.hpp`, `tools/delta_t_table.py`,
`tools/reference_check.py`, `tools/hour_angle_cio.py`, `tools/spp_dump.cpp`,
`tools/pipeline_check.cpp`, `tools/day_scan.cpp`, `tools/readme_numbers.cpp`.

Changed: `CMakeLists.txt`, `cmake/SppFlags.cmake`, `README.md`, `include/spp/core/{angles,calendar,time}.hpp`,
`include/spp/events.hpp`, `include/spp/model/espenak_meeus.hpp`, `src/core/{calendar,time}.cpp`,
`src/{earth,earth_state,events,sky,sun}.cpp`, `src/model/espenak_meeus.cpp`, `apps/spp/main.cpp`,
`apps/spiral/main.cpp`, `tools/{bench,cache_check,figures}.cpp`, `docs/images/*.svg`.

Public API additions (source-compatible): `is_valid`, `iso_date`, `events::upper_limb_altitude_deg`,
`events::upper_limb_crossing`. Behaviour changes on existing API: the `Time` factories and event
solvers throw `std::invalid_argument` on invalid input.
