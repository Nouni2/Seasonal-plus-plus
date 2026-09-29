<p align="center">
  <img src="docs/images/banner.svg" alt="Spp (Seasonal++): high-precision solar astronomy in C++20" width="100%">
</p>

<p align="center">
  <a href="#building"><img src="https://img.shields.io/badge/C%2B%2B-20-00599C?style=flat-square&logo=cplusplus&logoColor=white" alt="C++20"></a>
  <a href="#building"><img src="https://img.shields.io/badge/CMake-3.28%2B-064F8C?style=flat-square&logo=cmake&logoColor=white" alt="CMake 3.28+"></a>
  <a href="#3-heliocentric-earth-vsop87d"><img src="https://img.shields.io/badge/ephemeris-VSOP87D-E07B39?style=flat-square" alt="VSOP87D ephemeris"></a>
  <a href="#floating-point-policy"><img src="https://img.shields.io/badge/fp64-strict%20IEEE%20754-6E40C9?style=flat-square" alt="Strict IEEE 754 fp64"></a>
  <a href="#features"><img src="https://img.shields.io/badge/dependencies-none-2EA043?style=flat-square" alt="No dependencies"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-GPL--3.0-blue?style=flat-square" alt="License: GPL-3.0-or-later"></a>
</p>

<p align="center">
  <b>Where the Sun is, and when it rises, transits and sets: anywhere, from 1000 BC to 3000 AD.</b><br>
  Full VSOP87D ephemeris · strict IEEE 754 fp64 · ~0.4 µs per sky position · events solved to 1 µs
</p>

Spp is a C++20 library for high-precision solar positional astronomy. It
computes where the Sun is and when it rises, transits and sets: apparent
geocentric coordinates, topocentric refracted altitude and azimuth, sidereal
time, nutation, ΔT, and event instants for any place and date from 1000 BC to
3000 AD.

**The library is the engine.** Everything Spp can compute is available through
its C++ API (`include/spp/`). The `spp` executable is a small demonstration
front end that prints a few fixed results. It is not the engine, and what it
does not print says nothing about what the engine can compute. See
[What the `spp` executable is and is not](#what-the-spp-executable-is-and-is-not).

It is built for large batch workloads, such as every day of a year across a
grid of latitudes. It never trades accuracy for speed: every shortcut is
measured against an extended-precision evaluation of the same theory and must
match it to within double-precision rounding.

```
$ spp day 2026-06-21 48.8566 2.3522 35
type      NORMAL
transit   2026-06-21T11:52:24Z
sunrise   2026-06-21T03:46:28Z
sunset    2026-06-21T19:58:21Z
duration  16.1981 h
```

---

## Contents

- [Gallery](#gallery)
- [Features](#features)
- [Building](#building)
- [Command-line tool](#command-line-tool) (a demonstration front end, not the engine)
- [Azimuth spirals](#azimuth-spirals)
- [Library](#library)
- [Conventions](#conventions)
- [Theory](#theory)
  1. [Pipeline](#1-pipeline)
  2. [Time](#2-time)
  3. [Heliocentric Earth: VSOP87D](#3-heliocentric-earth-vsop87d)
  4. [Geometric geocentric Sun](#4-geometric-geocentric-sun)
  5. [Nutation](#5-nutation)
  6. [Obliquity of the ecliptic](#6-obliquity-of-the-ecliptic)
  7. [Aberration and the apparent longitude](#7-aberration-and-the-apparent-longitude)
  8. [Ecliptic to equatorial](#8-ecliptic-to-equatorial)
  9. [Sidereal time and hour angle](#9-sidereal-time-and-hour-angle)
  10. [Topocentric parallax](#10-topocentric-parallax)
  11. [Horizontal coordinates](#11-horizontal-coordinates)
  12. [Atmospheric refraction](#12-atmospheric-refraction)
  13. [Events](#13-events)
- [Numerical methods](#numerical-methods)
- [Accuracy](#accuracy)
- [Performance](#performance)
- [Accuracy tools](#accuracy-tools)
- [Project layout](#project-layout)
- [Limitations](#limitations)
- [References](#references)
- [License](#license)

---

## Gallery

Every figure in this README is computed by the library and drawn by
`spp_figures` (see [Accuracy tools](#accuracy-tools)).

**Daylight and twilight through the year in Paris.** Every boundary is an event
solved by `spp::events`: sunrise and sunset, then civil (−6°), nautical (−12°)
and astronomical (−18°) twilight. Around the June solstice astronomical
twilight lasts all night, because at 48.86° N the Sun never gets 18° below the
horizon.

![Daylight and twilight through the year in Paris, computed with Spp](docs/images/daylight.svg)

**Day length through the year, from the equator to 80° N.** The curves meet
just above 12 h, a few days before the March equinox and after the September
one: refraction and the Sun's radius make the day at the equinox slightly
longer than 12 h (12.12 h at the equator, 12.17 h in Paris on 20 March).
Above the Arctic Circle the curves reach the polar day (24 h) and polar night
(0 h) plateaus.

![Day length through the year by latitude](docs/images/day_length.svg)

**The Sun's path across the sky over Paris** at the June solstice, the March
equinox and the December solstice. On these dates solar noon falls a few
minutes before 12:00 UTC (11:52, 11:58 and 11:49): Paris lies 2.35° east of
Greenwich, 9 minutes of time, and the equation of time adds or removes a few
more.

![The Sun's path across the sky over Paris](docs/images/sun_paths.svg)

**The analemma.** The Sun photographed at the same clock time every day traces
a figure eight. The vertical extent is the change in declination over the year;
the horizontal swing is the equation of time, from the eccentricity of the
Earth's orbit and the tilt of its axis.

![The analemma seen from Greenwich](docs/images/analemma.svg)

**Where the Sun rises and sets.** The azimuth of sunrise and sunset over the
year as spirals, the same view as the interactive `spp_spiral` viewer. In
Paris, sunrise swings between 51.5° (northeast, June) and 125.9° (southeast,
December), and sunset between 234.1° and 308.5°.

![Sunrise and sunset azimuth spirals](docs/images/spirals.svg)

---

## Features

- **Solar ephemeris.** The full VSOP87D Earth series (2,425 terms) is compiled
  into the binary as exact hex-float literals, so there are no data files at
  runtime.
- **Apparent place.** Nutation, annual aberration and true obliquity, giving
  right ascension and declination referred to the true equator and equinox of
  date.
- **Earth rotation.** GMST, GAST and local apparent sidereal time, with optional
  DUT1.
- **Time scales.** UTC, UT1 and TT are distinct C++ types, so mixing them is a
  compile error. ΔT comes from observed values and long-term fits. Dates are
  stored as a split Julian Date for full precision over millennia.
- **Observer corrections.** Topocentric parallax on the WGS84 ellipsoid with
  elevation, and atmospheric refraction scaled for pressure and temperature.
- **Events.** Meridian transit, sunrise and sunset (upper limb), civil, nautical
  and astronomical twilight, or any target altitude. Days are classified as
  normal, polar day or polar night.
- **Speed.** About 0.4 µs for a refracted topocentric altitude and about 7 µs for
  a complete sunrise/transit/sunset day, on one core.
- **Portable standard C++.** The library has no dependencies beyond the
  standard library. Only the optional `spp_spiral` viewer uses OpenGL.

---

## Building

### Requirements

- A C++20 compiler; the executables use C++23 for `std::print`. Developed with
  GCC 15.2 (MSYS2 UCRT64 on Windows).
- CMake ≥ 3.28 and Ninja.
- Optional, for `spp_spiral`: GLFW 3, GLEW, GLM and OpenGL 3.3. The target is
  built automatically when they are found, and skipped otherwise.

### Build

```sh
cmake --preset release
cmake --build --preset release
./build/release/spp day 2026-06-21 48.8566 2.3522
```

| Preset | Output | Purpose |
|---|---|---|
| `release` | `build/release/` | Optimized |
| `debug` | `build/debug/` | Debug information, no optimization |

CMake options:

| Option | Default | Effect |
|---|---|---|
| `SPP_NATIVE` | `OFF` | Adds `-march=native` |
| `SPP_BUILD_APPS` | `ON` | Builds `spp` and, when GLFW/GLEW/GLM are found, `spp_spiral` |
| `SPP_BUILD_VALIDATION` | `ON` | Builds the `spp_cache_check` accuracy tool |
| `SPP_BUILD_BENCH` | `OFF` | Builds the `spp_bench` benchmark |
| `SPP_BUILD_FIGURES` | `OFF` | Builds `spp_figures`, which renders the README figures |

Other CMake projects link the library as `spp::spp`.

The build first compiles `tools/vsop_gen`, which turns `data/VSOP87D.ear` into
`build/<preset>/generated/spp/vsop87d_earth_terms.inc` (see
[VSOP87D tables](#vsop87d-tables)).

On Windows with MSYS2, the executables link against the UCRT64 runtime DLLs,
so `C:\msys64\ucrt64\bin` must be on `PATH` to run them.

### Floating-point policy

All targets are compiled with:

- `-fno-fast-math`: strict IEEE 754 semantics.
- `-ffp-contract=off`: no implicit fused multiply-add, so results do not depend
  on the target CPU.
- `-fexcess-precision=standard`: no hidden extended-precision intermediates.

All computation is in `double`.

---

## Command-line tool

### What the `spp` executable is and is not

**It is** a thin demonstration front end (`apps/spp/main.cpp`, about 180
lines). It parses arguments, calls a handful of library functions and prints
the results. It contains no astronomy of its own. It exists for quick one-off
queries and for checking results by eye.

**It is not** the engine, and not a complete interface to it. It offers three
fixed outputs (`sun`, `day`, `year`) and exposes only a small part of the API.
For example, it prints sunrise and sunset *times* but not the Sun's azimuth at
those times. It prints no twilight instants, no altitude curves and no events
for custom target altitudes. The library computes all of these.

**Rule:** if the executable doesn't print something, that is a limitation of
the executable, never of Spp. To get anything else, write a C++ program that
links the `spp` library and calls the API directly (see [Library](#library)).
Any quantity the API provides can be computed at any instant, combined with
others, and fed to whatever consumes it (files, plots, visualizations).
Extending the executable is optional and never required to use the engine.

### Commands

```
spp sun  YYYY-MM-DD HH:MM[:SS] [LAT LON [ELEV]]   Sun at a UTC instant
spp day  YYYY-MM-DD LAT LON [ELEV]                transit, sunrise, sunset, day length
spp year YYYY LAT LON [ELEV]                      one CSV row per day

options:
  --dut1=S          UT1 - UTC in seconds (default 0)
  --pressure=MBAR   surface pressure for refraction (default 1010)
  --temp=C          surface temperature for refraction (default 10)
  --no-refraction   geometric altitudes
```

Latitude is in degrees, positive North. Longitude is in degrees, positive East.
Elevation is in metres. All times are UTC. Negative (astronomical) years are
accepted, e.g. `-0100-03-01`.

### `spp sun`

At the March equinox of 2026, λ passes through 0°:

```
$ spp sun 2026-03-20 14:46 48.8566 2.3522 35
time        2026-03-20T14:46:00Z  (ΔT 69.27 s)
ecliptic    λ 0.00007180°  β 0.00009035°  R 0.995917594 AU
equatorial  α 0.00002994° (0.00000200 h)  δ 0.00011145°
nutation    Δψ 6.2175"  Δε 9.1309"  ε 23.43841881°
sidereal    GMST 39.64780049°  GAST 39.64938507°
horizontal  alt 29.298672° (geometric 29.268708°)  az 230.093725°
```

Without a location, only the geocentric lines are printed.

### `spp day`

Tromsø at the June solstice, above the Arctic Circle:

```
$ spp day 2026-06-21 69.6492 18.9553
type      POLAR_DAY
transit   2026-06-21T10:45:59Z
sunrise   -
sunset    -
duration  24.0000 h
```

Events are reported in UTC. For Sydney, sunrise falls on the previous UTC date:

```
$ spp day 2026-12-21 -33.8688 151.2093 58
type      NORMAL
transit   2026-12-21T01:53:01Z
sunrise   2026-12-20T18:40:18Z
sunset    2026-12-21T09:05:45Z
duration  14.4240 h
```

### `spp year`

```
$ spp year 2026 48.8566 2.3522
date,type,transit_utc,sunrise_utc,sunset_utc,duration_h
2026-01-01,NORMAL,11:54:09,07:43:36,16:04:54,8.354994
2026-01-02,NORMAL,11:54:37,07:43:32,16:05:54,8.372716
...
```

For 1582, the ten days removed by the Gregorian reform (October 5–14) are
skipped.

---

## Azimuth spirals

`spp_spiral` plots where the Sun rises and sets over a whole year, as two 3D
spirals in cylindrical coordinates. It is a separate program built directly on
the library. It computes the Sun's azimuth at each sunrise and sunset with the
public API, which the `spp` executable does not print.

```
spp_spiral YEAR LAT LON [ELEV]
```

![Sunrise and sunset azimuth spirals for Paris](docs/images/spirals.svg)

It opens two windows, **sunrise** and **sunset**. In each one, every day of
the year is a point:

| Coordinate | Meaning |
|---|---|
| Angle | Azimuth of the Sun at the event: 0° = North (red), 90° = East |
| Height | Day of the year: 1 January at the bottom, 31 December at the top |
| Radius | 1 (no meaning) |

- **Points** are colored by season: blue around January, orange around July.
- **Guides:** faint rings mark the 1st of each month, and N/E/S/W letters mark
  both ends of the cylinder.
- **Gaps:** days without the event, such as polar day or polar night, leave a
  gap in the spiral.
- **Console summary:** before the windows open, the program prints the azimuth
  range of each event and the dates where the extremes occur.

| Control | Action |
|---|---|
| Left drag | Orbit |
| Mouse wheel | Zoom |
| `T` | View from above (North up, East right) |
| `R` | Reset the view |
| `Esc` | Close the window |

How each azimuth is computed:

1. `events::day` gives the sunrise and sunset instants.
2. At each instant, `sun::horizontal` gives the Sun's topocentric azimuth: the
   apparent place, sidereal time, parallax and horizontal coordinates, with the
   same nutation throughout.

Refraction changes only the altitude, not the azimuth.

---

## Library

Link against the `spp` target and include `spp/spp.hpp`, or individual headers.
The library is C++20; this example prints with `std::println`, which needs
C++23.

```cpp
#include <print>
#include "spp/spp.hpp"

int main() {
    // An instant in UTC, and an observer (longitude positive East).
    const auto t = spp::Time::from_utc(2026, 3, 20, 14, 46, 0.0);
    const spp::Observer paris{48.8566, 2.3522, 35.0};

    // Apparent geocentric position of the Sun (uses Terrestrial Time).
    const spp::Equatorial eq = spp::sun::apparent(t.tt());
    std::println("RA {:.6f} h  Dec {:+.6f} deg", eq.ra_hours(), eq.dec_deg);

    // The Sun in the observer's sky: refracted altitude and azimuth.
    const spp::Horizontal sky = spp::sun::horizontal(t, paris);
    std::println("altitude {:+.6f} deg  azimuth {:.6f} deg", sky.alt_deg, sky.az_deg);

    // Sunrise, sunset and day length. Next to polar day or night either can be
    // missing, so check before use.
    const auto day = spp::events::day({2026, 6, 21}, paris);
    if (day.rise && day.set)
        std::println("sunrise {}  sunset {}  {:.4f} h", day.rise->iso_utc(), day.set->iso_utc(),
                     day.duration_hours);

    // Civil twilight: the Sun's center 6 deg below the horizon, either side of transit.
    const auto dusk = spp::events::crossing(day.transit, paris, {}, spp::events::civil_twilight_deg,
                                            spp::events::Direction::After);
    if (dusk) std::println("civil dusk {}", dusk->iso_utc());
}
```

```
RA 0.000002 h  Dec +0.000111 deg
altitude +29.298672 deg  azimuth 230.093725 deg
sunrise 2026-06-21T03:46:28Z  sunset 2026-06-21T19:58:21Z  16.1981 h
civil dusk 2026-06-21T20:40:32Z
```

### API overview

The public API has two layers. The object namespaces (`sun`, `earth`, `sky`,
`events`) are what most code uses. `spp::model::<theory>` exposes the
underlying theories directly, each named after its source, so which model
produced a number is always explicit.

#### `spp::` — core types (`spp/core/`)

| Symbol | Description |
|---|---|
| `Time` | An instant held in UTC, with its DUT1 |
| `Time::from_utc(y, mo, d, h, mi, s, dut1)` | From a civil UTC date and time |
| `Time::from_civil(CivilDateTime, dut1)` | Same, from a struct |
| `Time::from_jd(day, frac, dut1)` | From a split Julian Date (normalized) |
| `t.utc()`, `t.ut1()`, `t.tt()` | The instant as `Instant<UTC/UT1/TT>` |
| `t.delta_t()`, `t.dut1()` | ΔT = TT − UT1 and DUT1 = UT1 − UTC, seconds |
| `t.decimal_year()` | The year value used to evaluate ΔT |
| `t.shifted(days)` | Another instant, same DUT1 |
| `t.to_civil()`, `t.iso_utc()` | Civil UTC date and time; ISO 8601 rounded to the second |
| `Instant<Scale>` | `{day, frac}` split Julian Date; `days_since_j2000()`, `centuries_since_j2000()` |
| `CivilDate`, `CivilDateTime` | Calendar dates (Julian before 1582-10-15, Gregorian after) |
| `is_gregorian(CivilDate)` | Whether a date is on or after 1582-10-15 |
| `Ecliptic` | `lon_deg`, `lat_deg`, `dist_au` |
| `Equatorial` | `ra_deg`, `dec_deg`, `dist_au`, `ra_hours()` |
| `Horizontal` | `alt_deg`, `az_deg` |
| `Observer` | `lat_deg`, `lon_deg` (East +), `elev_m` |
| `Atmosphere` | `pressure_mbar` (1010), `temp_c` (10), `refraction` (true) |
| `norm360`, `wrap180` | Angle normalization |
| Constants (`spp/core/constants.hpp`) | `deg2rad`, `rad2deg`, `jd_j2000_day`, `days_per_century`, WGS84, solar constants |

#### `spp::sun` — the Sun from the geocenter (`spp/sun.hpp`)

| Function | Returns |
|---|---|
| `geometric(Instant<TT>)` | Geometric ecliptic position ([§4](#4-geometric-geocentric-sun)) |
| `apparent_ecliptic(Instant<TT>, Nutation)` | Plus nutation in longitude and aberration ([§7](#7-aberration-and-the-apparent-longitude)) |
| `apparent(Instant<TT> [, Nutation])` | Apparent RA and Dec, true equator and equinox of date ([§8](#8-ecliptic-to-equatorial)) |
| `horizontal(Time, Observer [, Atmosphere])` | The Sun in the observer's sky: topocentric azimuth and refracted altitude ([§10](#10-topocentric-parallax)–[§12](#12-atmospheric-refraction)) |

#### `spp::earth` — Earth orientation (`spp/earth.hpp`)

| Function | Returns |
|---|---|
| `nutation(Instant<TT>)` | `Nutation{dpsi_deg, deps_deg, mean_obliquity_deg, true_obliquity_deg}` ([§5](#5-nutation), [§6](#6-obliquity-of-the-ecliptic)) |
| `mean_obliquity_deg(Instant<TT>)` | ε₀ |
| `delta_t(Time)` | ΔT, seconds ([§2.6](#26-δt)) |
| `gmst_deg(Instant<UT1>)` | Greenwich mean sidereal time ([§9](#9-sidereal-time-and-hour-angle)) |
| `equation_of_equinoxes_deg(Nutation)` | Δψ cos ε |
| `gast_deg(Instant<UT1>, Nutation)`, `gast_deg(Time)` | Greenwich apparent sidereal time |
| `lst_deg(Instant<UT1>, Nutation, lon)`, `lst_deg(Time, lon)` | Local apparent sidereal time |

#### `spp::sky` — the observer's sky (`spp/sky.hpp`)

| Function | Returns |
|---|---|
| `to_equatorial(Ecliptic, true_obliquity_deg)` | Ecliptic → equatorial ([§8](#8-ecliptic-to-equatorial)) |
| `to_horizontal(Equatorial, lat_deg, lst_deg)` | Equatorial → horizontal, geometric ([§11](#11-horizontal-coordinates)) |
| `parallax(Equatorial, Observer, lst_deg)` | Geocentric → topocentric ([§10](#10-topocentric-parallax)) |
| `refraction(geometric_alt_deg, Atmosphere)` | Apparent (refracted) altitude ([§12](#12-atmospheric-refraction)) |
| `observe(Equatorial, Observer, lst_deg, Atmosphere)` | Parallax, horizontal coordinates and refraction in one call: what the observer sees |

#### `spp::events` — events in time (`spp/events.hpp`)

| Symbol | Description |
|---|---|
| `altitude_deg(Time, Observer, Atmosphere)` | Full pipeline: apparent, topocentric, refracted altitude of the Sun's center |
| `transit(approx, lon, cfg)` | Meridian transit near an approximate instant ([§13.3](#133-transit)) |
| `crossing(transit, obs, atm, target_alt_deg, dir, cfg)` | Instant the Sun's center reaches a target altitude before or after transit; `nullopt` if it never does within 12 h ([§13.4](#134-crossings)) |
| `day(date, obs, atm, dut1, cfg)` | `Day{date, type, transit, rise, set, duration_hours}` ([§13.2](#132-one-day)) |
| `upper_limb_deg` | −16′: sunrise/sunset |
| `civil_twilight_deg`, `nautical_twilight_deg`, `astronomical_twilight_deg` | −6°, −12°, −18° |
| `DayType` | `Normal`, `PolarDay`, `PolarNight` |
| `Direction` | `Before` (rising side), `After` (setting side) |
| `SolverConfig` | `time_tolerance_s` (1e-6), `max_iter` (60) |

#### `spp::model` — the theories (`spp/model/`)

| Namespace | Contents |
|---|---|
| `vsop87d` | `earth(t)`: heliocentric L, B, R of the Earth, t in Julian millennia of TT ([§3](#3-heliocentric-earth-vsop87d)) |
| `iau1980` | `nutation(T)`: Δψ, Δε in arcseconds ([§5](#5-nutation)) |
| `laskar` | `mean_obliquity_deg(T)` ([§6](#6-obliquity-of-the-ecliptic)) |
| `meeus` | Civil date ↔ Julian Date ([§2.2](#22-civil-date-to-julian-date), [§2.3](#23-julian-date-to-civil-date)); GMST ([§9](#9-sidereal-time-and-hour-angle)) |
| `espenak_meeus` | `delta_t_seconds(decimal_year)` ([§2.6](#26-δt)) |
| `saemundsson` | `apparent_altitude_deg(h, pressure, temperature)` ([§12](#12-atmospheric-refraction)) |

### Thread safety

All functions are reentrant. The ephemeris cache (see
[Chebyshev earth state](#chebyshev-earth-state)) is `thread_local`, so threads
never share or lock it. Each thread that computes positions holds 67 KB of
cache.

---

## Conventions

| Quantity | Convention |
|---|---|
| Angles | Degrees unless the name says otherwise (`_rad`, `_arcsec`, `ra_hours()`) |
| Longitude | Positive East, −180 … +180 |
| Latitude | Geodetic, positive North |
| Azimuth | From North, increasing East (N = 0°, E = 90°, S = 180°, W = 270°) |
| Hour angle | Positive West of the meridian |
| Right ascension | Degrees in `ra_deg`; hours via `ra_hours()` |
| Elevation | Metres, applied on the WGS84 ellipsoid |
| Distance | Astronomical units |

### Time

- **Julian Dates** are stored split into an integer day and a fraction in
  [0, 1), where fraction 0 is noon ([§2.1](#21-julian-dates-and-split-storage)).
- **UTC** is the input and output scale. **UT1** = UTC + DUT1 drives Earth
  rotation (sidereal time). **TT** = UT1 + ΔT drives the ephemeris.
- **DUT1** defaults to 0, which makes UT1 equal to UTC. |UT1 − UTC| < 0.9 s,
  which bounds the resulting rotation error to 13.5″. Pass the IERS value for
  sub-arcsecond hour angles.
- **Calendar.** The Julian calendar is used before 1582-10-15 and the Gregorian
  from that date on. Years use astronomical numbering: 1 BC is year 0, 2 BC is
  year −1.

### Events

- **Sunrise and sunset** are the instants the apparent (refracted) altitude of
  the Sun's center equals −16′, i.e. the upper limb touches the horizon for a
  fixed semidiameter of 16′. See [§13.1](#131-definitions) for how this
  relates to the usual −0.8333° almanac convention.
- **Transit** is the instant the hour angle of the apparent Sun is zero.
- **`events::day(date, …)`** searches around local solar noon of `date`, taken
  as 12:00 − longitude/15 h UTC. Rise and set are searched within 12 h either
  side of that transit, so they can fall on the neighbouring UTC date.
- **Day type** is decided from the altitudes at transit and 12 h before and
  after it. Next to polar day, a normal day can have a sunset but no sunrise
  (or the reverse) ([§13.2](#132-one-day)).

---

## Theory

This section documents every physical model and algorithm in the engine, with
the exact formulas, constants and operation choices used by the code.
Throughout:

- φ is the observer's geodetic latitude, λ<sub>obs</sub> the longitude (East
  positive), h<sub>e</sub> the height above the ellipsoid.
- λ, β are ecliptic longitude and latitude; α, δ right ascension and
  declination; h, A altitude and azimuth; H the hour angle.
- T is Julian centuries of TT since J2000.0, and t = T/10 Julian millennia.

### 1. Pipeline

```
UTC ──(+DUT1)──► UT1 ──(+ΔT)──► TT
                  │               ├─► VSOP87D Earth L, B, R ─► geometric Sun λ, β, R      §3, §4
                  │               ├─► nutation Δψ, Δε ; obliquity ε₀, ε = ε₀ + Δε        §5, §6
                  │               └─► apparent λ = λ + Δψ − aberration ─► α, δ            §7, §8
                  └─► GMST ─► + Δψ cos ε ─► GAST ─► + λobs ─► LST ─► H = LST − α        §9
                                                                        │
            α, δ, H, observer (WGS84) ─► parallax ─► α′, δ′ ─► h, A ─► refraction ─► h_app   §10–§12
                                                                                        │
                                                        root finding in time ─► events   §13
```

Every stage exists as a public function (see [API overview](#api-overview)),
so any intermediate quantity can be obtained on its own.

### 2. Time

#### 2.1 Julian Dates and split storage

The Julian Date (JD) counts days continuously from 4713 BC January 1, 12:00
(Julian proleptic calendar). A JD near the present is about 2.46 × 10⁶; stored
in one `double`, which lies in [2²¹, 2²²), the spacing between adjacent values
would be 2⁻³¹ day ≈ 40 µs (a rounding error of up to 2⁻³² day ≈ 20 µs).

Spp stores every instant as an integer day plus a fraction:

$$
\mathrm{JD} = d + f,\qquad d \in \mathbb{Z},\quad f \in [0, 1),
$$

where f = 0 is 12:00. The fraction's resolution is at most 2⁻⁵³ day ≈ 10 ps,
whatever the date.

Construction normalizes any (d, f):

1. k = trunc(f); d ← d + k; f ← f − k.
2. If f < 0: d ← d − 1; f ← f + 1.

Adding an offset (DUT1, ΔT) renormalizes the same way: a fraction ≥ 1 moves
⌊f⌋ days into d, and a negative fraction borrows ⌈|f|⌉ days.

Differences are formed from the split parts, e.g. days since J2000.0:

$$
D = (d - 2451545) + f .
$$

The integer subtraction is exact, so D carries the full precision of f.

#### 2.2 Civil date to Julian Date

Meeus, Algorithm 7.1, returning the integer JD of the date's noon.

Given year Y, month M (1–12), day D, and time of day h:m:s:

1. If M ≤ 2: Y ← Y − 1, M ← M + 12. January and February are counted as months
   13 and 14 of the previous year, so the leap day falls at the end of the
   counting year.
2. If the date is on or after 1582-10-15 (Gregorian calendar):
   A = ⌊Y/100⌋, B = 2 − A + ⌊A/4⌋. Before that date (Julian calendar), B = 0.
   B is the accumulated difference between the two calendars' leap-year rules.
3. Integer JD at noon:

$$
d = \lfloor 365.25\,(Y + 4716) \rfloor + \lfloor 30.6001\,(M + 1) \rfloor + D + B - 1524 .
$$

4. Fraction: f = h/24 + m/1440 + s/86400 − 0.5, then normalized
   ([§2.1](#21-julian-dates-and-split-storage)).

The term ⌊30.6001 (M + 1)⌋ reproduces the cumulative month lengths of the
March-based year. The factor 30.6001 instead of 30.6 guards against 30.6 × k
landing just below an integer in floating point. The Gregorian test uses the
date as given, before the month shift.

#### 2.3 Julian Date to civil date

Meeus, Chapter 7, inverse. From (d, f):

1. g = f + 0.5; Z = d + ⌊g⌋; F = g − ⌊g⌋. This moves the day boundary to
   civil midnight.
2. If Z < 2299161 (1582-10-15): A = Z. Otherwise
   α = ⌊(Z − 1867216.25)/36524.25⌋ and A = Z + 1 + α − ⌊α/4⌋.
3. B = A + 1524, C = ⌊(B − 122.1)/365.25⌋, D = ⌊365.25 C⌋,
   E = ⌊(B − D)/30.6001⌋.
4. Day with fraction: B − D − ⌊30.6001 E⌋ + F. The integer part is the day of
   the month, and the remainder the fraction of the civil day.
5. Month = E − 1 if E < 14, else E − 13. Year = C − 4716 if month > 2, else
   C − 4715.

The time of day is split from s = fraction × 86400 as hour = ⌊s/3600⌋,
minute = ⌊(s mod 3600)/60⌋ and second = s − 3600·hour − 60·minute. Floor
division and modulo are derived from `fmod`, so they stay consistent with each
other at exact multiples.

`iso_utc()` rounds to the nearest second by adding 0.5 s to the JD fraction
before converting, so a time such as 23:59:59.7 correctly rolls over to the
next date.

#### 2.4 Epochs

J2000.0 is JD 2451545.0 TT (2000-01-01 12:00 TT). The time arguments of the
theories are

$$
T = \frac{\mathrm{JD} - 2451545}{36525}\ \text{(Julian centuries)},\qquad
t = \frac{T}{10}\ \text{(Julian millennia)},
$$

computed from the split form as T = D / 36525, with D from
[§2.1](#21-julian-dates-and-split-storage). Unless noted, T is on the TT scale.

#### 2.5 Time scales

| Scale | Definition in Spp | Drives |
|---|---|---|
| UTC | Input and output | Civil time |
| UT1 | UTC + DUT1 | Earth rotation (§9) |
| TT | UT1 + ΔT = UTC + DUT1 + ΔT | Ephemeris, nutation, obliquity (§3–§8) |

In code the conversions add (DUT1)/86400 or (ΔT + DUT1)/86400 to the UTC
fraction and renormalize. `Instant<UTC>`, `Instant<UT1>` and `Instant<TT>` are
distinct types, so passing a UTC instant where TT is required does not compile.

Sensitivity to time errors:

- The Sun moves 360°/365.2422 d ≈ 0.0411″ per second of time in longitude, so
  a 1 s error in TT moves the computed Sun by only 0.04″.
- The Earth rotates 15.04″ per second of UT1, so a 1 s error in UT1 shifts
  hour angles, altitudes and azimuths by up to 15″. With DUT1 = 0 the error is
  |UT1 − UTC| < 0.9 s, i.e. under 13.5″.

**Leap seconds.** UTC is treated as a uniform scale. A leap second
(23:59:60) has no distinct representation and maps onto 00:00:00 of the next
day. TT is derived as UT1 + ΔT rather than through TAI, so leap seconds do not
enter TT except through the caller's DUT1.

#### 2.6 ΔT

ΔT = TT − UT1 is evaluated at the decimal year

$$
y = 2000 + 100\,T_{\mathrm{UTC}},
$$

counted in Julian years from 2000-01-01 12:00 UTC. It differs from the
calendar-year fraction by under a day, which is negligible for ΔT. The pieces
are evaluated in this order; the first one whose range contains y is used.

**Observed values, 2005 ≤ y < 2026.** Annual values at the start of each year
(IERS/USNO), interpolated linearly between consecutive years:

| Year | ΔT (s) | Year | ΔT (s) | Year | ΔT (s) |
|---|---|---|---|---|---|
| 2005 | 64.69 | 2013 | 66.91 | 2021 | 69.36 |
| 2006 | 64.85 | 2014 | 67.28 | 2022 | 69.29 |
| 2007 | 65.15 | 2015 | 67.64 | 2023 | 69.18 |
| 2008 | 65.46 | 2016 | 68.10 | 2024 | 69.10 |
| 2009 | 65.78 | 2017 | 68.59 | 2025 | 69.07 |
| 2010 | 66.07 | 2018 | 68.97 | 2026 | 69.06 (provisional) |
| 2011 | 66.32 | 2019 | 69.22 | | |
| 2012 | 66.60 | 2020 | 69.36 | | |

**Near future, 2026 ≤ y < 2050.** A linear bridge from the last value to the
Espenak–Meeus estimate of 93 s at 2050:

$$
\Delta T = 69.06 + (93 - 69.06)\,\frac{y - 2026}{24}.
$$

This bridge rises by about 1 s per year, whereas observed ΔT has been flat or
slightly falling since 2021 (see the table): the Earth's rotation has recently
sped up, against the long-term trend the fits assume. Since
ΔT = 32.184 s + (TAI − UTC) − (UT1 − UTC), with TAI − UTC = 37 s since 2017 and
|UT1 − UTC| < 0.9 s, the current value stays near 69.1 s until a leap second
is introduced, while the bridge already gives about 69.8 s by late 2026. The
resulting error in the Sun's position is under 0.03″ (0.0411″ per second of
ΔT, [§2.5](#25-time-scales)); hour angles are unaffected, since they depend on
UT1, not ΔT.

**Prediction, 2050 ≤ y < 2150.** With u = (y − 1820)/100:

$$
\Delta T = -20 + 32u^2 - 0.5628\,(2150 - y).
$$

**Historical fits (Espenak & Meeus).** With the variable shown:

| Range | Variable | ΔT (s) |
|---|---|---|
| 1986 ≤ y < 2005 | t = y − 2000 | 63.86 + 0.3345t − 0.060374t² + 0.0017275t³ + 0.000651814t⁴ + 0.00002373599t⁵ |
| 1961 ≤ y < 1986 | t = y − 1975 | 45.45 + 1.067t − t²/260 − t³/718 |
| 1941 ≤ y < 1961 | t = y − 1950 | 29.07 + 0.407t − t²/233 + t³/2547 |
| 1920 ≤ y < 1941 | t = y − 1920 | 21.20 + 0.84493t − 0.076100t² + 0.0020936t³ |
| 1900 ≤ y < 1920 | t = y − 1900 | −2.79 + 1.494119t − 0.0598939t² + 0.0061966t³ − 0.000197t⁴ |
| 1860 ≤ y < 1900 | t = y − 1860 | 7.62 + 0.5737t − 0.251754t² + 0.01680668t³ − 0.0004473624t⁴ + t⁵/233174 |
| 1800 ≤ y < 1860 | t = y − 1800 | 13.72 − 0.332447t + 0.0068612t² + 0.0041116t³ − 0.00037436t⁴ + 0.0000121272t⁵ − 0.0000001699t⁶ + 0.000000000875t⁷ |
| 1700 ≤ y < 1800 | t = y − 1700 | 8.83 + 0.1603t − 0.0059285t² + 0.00013336t³ − t⁴/1174000 |
| 1600 ≤ y < 1700 | t = y − 1600 | 120 − 0.9808t − 0.01532t² + t³/7129 |
| 500 ≤ y < 1600 | u = (y − 1000)/100 | 1574.2 − 556.01u + 71.23472u² + 0.319781u³ − 0.8503463u⁴ − 0.005050998u⁵ + 0.0083572073u⁶ |
| −500 ≤ y < 500 | u = y/100 | 10583.6 − 1014.41u + 33.78311u² − 5.952053u³ − 0.1798452u⁴ + 0.022174192u⁵ + 0.0090316521u⁶ |
| y < −500 or y ≥ 2150 | u = (y − 1820)/100 | −20 + 32u² (long-term tidal parabola) |

In the 1800–1860 fit the large terms nearly cancel, and the t⁷ term is
essential: without it the fit diverges by thousands of seconds by 1860.

ΔT cannot be overridden through the API; see [Limitations](#limitations).

### 3. Heliocentric Earth: VSOP87D

VSOP87 (Bretagnon & Francou 1988) is a semi-analytical theory of the planets.
Version D gives heliocentric spherical coordinates referred to the dynamical
ecliptic and equinox **of date**, so general precession is built into the
series. The body used is the Earth itself (not the Earth–Moon barycenter), so
the Moon's perturbation of the Earth is included.

Each coordinate X ∈ {L, B, R} is a Poisson series in t (Julian millennia of
TT):

$$
X(t) = \sum_{n=0}^{5} t^{\,n} \sum_{i} A_{n,i} \cos\!\left(B_{n,i} + C_{n,i}\,t\right),
$$

with L and B in radians and R in astronomical units. Spp evaluates the full
series, with no truncation:

| Block | t⁰ | t¹ | t² | t³ | t⁴ | t⁵ | Total |
|---|---|---|---|---|---|---|---|
| L (longitude) | 559 | 341 | 142 | 22 | 11 | 5 | 1,080 |
| B (latitude) | 184 | 99 | 49 | 11 | 5 | — | 348 |
| R (radius vector) | 526 | 292 | 139 | 27 | 10 | 3 | 997 |
| | | | | | | | **2,425** |

- L is returned unwrapped: its t¹ block contains the mean motion
  (≈ 6283.08 rad per millennium), so L grows to about 2 × 10⁴ rad at 3000 AD.
  It is reduced to [0°, 360°) only when converted to the Sun's longitude.
- The fastest term in the Earth series has C = 327,574 rad per millennium, a
  period of 7.0 days, and an amplitude between 10⁻¹⁰ and 10⁻⁸ rad. Every term
  with amplitude above 10⁻⁸ rad has a period of at least 9.4 days. These
  figures set the Chebyshev segment design ([Numerical methods](#chebyshev-earth-state)).

#### VSOP87D tables

The data file `data/VSOP87D.ear` holds, per term, index columns and the
coefficients S, K, A, B, C. Only A, B and C are needed for the reduced form
above. At build time `tools/vsop_gen`:

1. Parses each block header (variable 1/2/3 = L/B/R, power of t, term count)
   and checks the counts, duplicate blocks and malformed lines.
2. Parses A, B, C with `std::from_chars`, which rounds correctly: each value is
   the double nearest to the decimal in the file.
3. Writes them as hex-float literals, which the compiler reads back
   bit-for-bit. No decimal round trip can alter a coefficient.

### 4. Geometric geocentric Sun

In the heliocentric ecliptic frame, the Earth is at

$$
\mathbf{r}_\oplus = R\,(\cos B\cos L,\ \cos B\sin L,\ \sin B).
$$

The Sun seen from the Earth is at −**r**<sub>⊕</sub>. Negating a unit vector
maps (L, B) to (L + 180°, −B), so the geometric geocentric position is

$$
\lambda = L + 180^\circ \pmod{360^\circ},\qquad \beta = -B,\qquad \Delta = R .
$$

Spp uses this closed form directly. Building the rectangular vector, negating
it and converting back with atan2 gives the same angles, with more rounding.

β is tiny (the Earth's heliocentric latitude on the ecliptic of date stays
below 1.3″ over 1000 BC–3000 AD), but it is kept throughout.

### 5. Nutation

Nutation is the short-period motion of the Earth's rotation axis, caused
mainly by the Moon's torque on the equatorial bulge. It shifts the true
equinox along the ecliptic by Δψ (nutation in longitude) and tilts the true
equator by Δε (nutation in obliquity).

Spp uses the IAU 1980 theory, truncated to its 20 largest terms (Meeus,
Table 22.A). The five fundamental luni-solar arguments, in degrees, with T in
Julian centuries of TT, each reduced modulo 360°:

| Argument | Meaning | Value (°) |
|---|---|---|
| D | Mean elongation of the Moon from the Sun | 297.85036 + 445267.111480T − 0.0019142T² + T³/189474 |
| M | Mean anomaly of the Sun | 357.52772 + 35999.050340T − 0.0001603T² − T³/300000 |
| M′ | Mean anomaly of the Moon | 134.96298 + 477198.867398T + 0.0086972T² + T³/56250 |
| F | Moon's argument of latitude | 93.27191 + 483202.017538T − 0.0036825T² + T³/327270 |
| Ω | Longitude of the Moon's ascending node | 125.04452 − 1934.136261T + 0.0020708T² + T³/450000 |

Each term i has an argument θ<sub>i</sub> = d·D + m·M + m′·M′ + f·F + ω·Ω, and

$$
\Delta\psi = \sum_i (S_i + S'_i\,T)\,\sin\theta_i,\qquad
\Delta\varepsilon = \sum_i (C_i + C'_i\,T)\,\cos\theta_i ,
$$

with coefficients in units of 0.0001″:

| d | m | m′ | f | ω | S | S′ | C | C′ |
|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 1 | −171996 | −174.2 | 92025 | 8.9 |
| −2 | 0 | 0 | 2 | 2 | −13187 | −1.6 | 5736 | −3.1 |
| 0 | 0 | 0 | 2 | 2 | −2274 | −0.2 | 977 | −0.5 |
| 0 | 0 | 0 | 0 | 2 | 2062 | 0.2 | −895 | 0.5 |
| 0 | 1 | 0 | 0 | 0 | 1426 | −3.4 | 54 | −0.1 |
| 0 | 0 | 1 | 0 | 0 | 712 | 0.1 | −7 | 0 |
| −2 | 1 | 0 | 2 | 2 | −517 | 1.2 | 224 | −0.6 |
| 0 | 0 | 0 | 2 | 1 | −386 | −0.4 | 200 | 0 |
| 0 | 0 | 1 | 2 | 2 | −301 | 0 | 129 | −0.1 |
| −2 | −1 | 0 | 2 | 2 | 217 | −0.5 | −95 | 0.3 |
| −2 | 0 | 1 | 0 | 0 | −158 | 0 | 0 | 0 |
| −2 | 0 | 0 | 2 | 1 | 129 | 0.1 | −70 | 0 |
| 0 | 0 | −1 | 2 | 2 | 123 | 0 | −53 | 0 |
| 2 | 0 | 0 | 0 | 0 | 63 | 0 | 0 | 0 |
| 0 | 0 | 1 | 0 | 1 | 63 | 0.1 | −33 | 0 |
| 2 | 0 | −1 | 2 | 2 | −59 | 0 | 26 | 0 |
| 0 | 0 | −1 | 0 | 1 | −58 | −0.1 | 32 | 0 |
| 0 | 0 | 1 | 2 | 1 | −51 | 0 | 27 | 0 |
| −2 | 0 | 2 | 0 | 0 | 48 | 0 | 0 | 0 |
| 0 | 0 | −2 | 2 | 1 | 46 | 0 | −24 | 0 |

- The dominant term, with argument Ω, has the 18.6-year period of the lunar
  node and an amplitude of 17.2″ in Δψ and 9.2″ in Δε.
- The next largest, 2F − 2D + 2Ω, has a half-year period (1.3″ in Δψ).
- The fastest term in the table, M′ + 2F + Ω, has a period of 9.1 days.
- The terms omitted from the full IAU 1980 series are each below about
  0.004″; their combined effect is a few hundredths of an arcsecond.

### 6. Obliquity of the ecliptic

The mean obliquity ε₀, the angle between the mean equator and the ecliptic of
date, follows Laskar (1986), valid for |U| < 1 (±10,000 years), with U = T/100:

$$
\begin{aligned}
\varepsilon_0 = {} & 84381.448'' - 4680.93''U - 1.55''U^2 + 1999.25''U^3 - 51.38''U^4 - 249.67''U^5 \\
& - 39.05''U^6 + 7.12''U^7 + 27.87''U^8 + 5.79''U^9 + 2.45''U^{10} .
\end{aligned}
$$

Spp evaluates it in Horner form. The true obliquity, the angle between the
true equator and the ecliptic, is

$$
\varepsilon = \varepsilon_0 + \Delta\varepsilon .
$$

At J2000.0, ε₀ = 23° 26′ 21.448″.

### 7. Aberration and the apparent longitude

Light from the Sun takes about 8.3 minutes to reach the Earth, and the Earth
moves at about 30 km/s across the line of sight. The resulting displacement of
the Sun's apparent direction is annual aberration.

In the heliocentric frame the Sun is at rest, so light-time and aberration are
the same effect seen in two frames. It is applied once.

The displacement in longitude is the Earth's velocity perpendicular to the
Sun's direction, divided by c. For a Keplerian orbit that component is h/R,
with h the constant specific angular momentum, so the 1/R dependence is exact.
The classical constant of aberration, κ = 20.49552″, is defined with the
factor 1/√(1 − e²); writing h in terms of a, e and the period then gives,
with a = 1 AU and no approximation in e,

$$
\Delta\lambda_{\text{ab}} = -\,\kappa\,\frac{1 - e^2}{R} = -\,\frac{20.4898''}{R}
\qquad(\text{Meeus 25.10}),
$$

since 20.49552″ × (1 − 0.0167086²) = 20.4898″. What is neglected is the effect
in latitude, which is proportional to β, and the Sun's own motion about the
barycenter of the solar system.

The apparent ecliptic position is then

$$
\lambda_{\text{app}} = \lambda + \Delta\psi - \frac{20.4898''}{R} \pmod{360^\circ},
\qquad \beta_{\text{app}} = \beta .
$$

Adding Δψ moves the longitude from the mean to the true equinox of date.

### 8. Ecliptic to equatorial

The equatorial frame of date is the ecliptic frame rotated about the equinox
direction (the x axis) by the true obliquity ε. A unit vector
(cos β cos λ, cos β sin λ, sin β) becomes

$$
x' = x,\qquad y' = y\cos\varepsilon - z\sin\varepsilon,\qquad z' = y\sin\varepsilon + z\cos\varepsilon ,
$$

and reading off the angles gives

$$
\sin\delta = \sin\beta\cos\varepsilon + \cos\beta\sin\varepsilon\sin\lambda,
\qquad
\alpha = \operatorname{atan2}\!\left(\sin\lambda\cos\varepsilon - \tan\beta\sin\varepsilon,\ \cos\lambda\right).
$$

The second follows from y′/x′ after dividing both by cos β. sin δ is clamped
to [−1, 1] before the arcsine, and α is reduced to [0°, 360°). The distance is
unchanged.

### 9. Sidereal time and hour angle

**GMST.** The Greenwich mean sidereal time is the IAU 1982 expression in the
form of Meeus (12.4), with D days and T = D/36525 centuries of **UT1** since
J2000.0:

$$
\theta_0 = 280.46061837 + 360.98564736629\,D + 0.000387933\,T^2 - \frac{T^3}{38710000}\qquad(\text{degrees}).
$$

**Exact reduction.** Evaluated as written, 360.98564736629 × D reaches about
4 × 10⁸° at 1000 BC. A double resolves only about 6 × 10⁻⁸° there, and the
decimal coefficient itself is stored with a relative error of about 10⁻¹⁶.
Spp removes the large part exactly:

1. Split D = D<sub>int</sub> + f, with D<sub>int</sub> = d − 2451545 an
   integer and f the JD fraction.
2. Write 360.98564736629 = 360 + r, r = 0.98564736629. Then
   360·D<sub>int</sub> is a whole number of turns and drops out modulo 360°.
3. Split r = r<sub>head</sub> + r<sub>tail</sub>, where r<sub>head</sub> is r
   rounded to 24 significant bits, so r<sub>head</sub>·D<sub>int</sub> is exact
   for |D<sub>int</sub>| < 2²⁹. Its reduction modulo 360° with `fmod` is also
   exact.
4. θ₀ = 280.46061837 + 360f + fmod(r<sub>head</sub>·D<sub>int</sub>, 360) +
   r<sub>tail</sub>·D<sub>int</sub> + r·f + cubic terms, then reduced to
   [0°, 360°).

The largest intermediate stays below about 10³°. Compared with the formula
evaluated in exact rational arithmetic over 1000 BC–3000 AD, the direct form
errs by up to 3.8 × 10⁻⁴″ and the split form by 1.9 × 10⁻⁷″.

**Apparent and local sidereal time.** The equation of the equinoxes
converts mean to apparent sidereal time:

$$
\mathrm{GAST} = \mathrm{GMST} + \Delta\psi\cos\varepsilon ,\qquad
\mathrm{LST} = \mathrm{GMST} + \lambda_{\text{obs}} + \Delta\psi\cos\varepsilon \pmod{360^\circ}.
$$

The Δψ and ε used here are the same values used for the apparent place
(§5–§7), computed once per instant. The complementary terms added to the
equation of the equinoxes in 1994 (0.00264″ sin Ω + 0.000063″ sin 2Ω) are
not included.

**Hour angle.** H = LST − α, wrapped to [−180°, 180°). It is positive west of
the meridian.

### 10. Topocentric parallax

The observer is not at the Earth's center. Seen from the surface, the Sun is
displaced by up to its horizontal parallax, about 8.8″.

**Observer position.** On the WGS84 ellipsoid (a = 6,378,137 m,
f = 1/298.257223563, b = a(1 − f)), the point at geodetic latitude φ lies on
the meridian ellipse at (a cos u, b sin u), where u is the reduced latitude:

$$
\tan u = (1 - f)\tan\varphi .
$$

Moving a height h<sub>e</sub> along the ellipsoid normal adds
h<sub>e</sub>(cos φ, sin φ). In units of a, the observer's distance from the
axis and from the equatorial plane are

$$
\rho\cos\varphi' = \cos u + \frac{h_e}{a}\cos\varphi,\qquad
\rho\sin\varphi' = (1 - f)\sin u + \frac{h_e}{a}\sin\varphi .
$$

**Horizontal parallax.** The Sun's equatorial horizontal parallax at distance
R (AU) is π = 8.794″/R, and sin π is the Earth's equatorial radius over the
Sun's distance.

**Shift.** Subtracting the observer's vector from the Sun's geocentric vector
(of length 1/sin π Earth radii) and converting back to angles gives
(Meeus 40.2, 40.3)

$$
\Delta\alpha = \operatorname{atan2}\!\left(-\rho\cos\varphi'\sin\pi\sin H,\ \cos\delta - \rho\cos\varphi'\sin\pi\cos H\right),
$$

$$
\delta' = \operatorname{atan2}\!\left((\sin\delta - \rho\sin\varphi'\sin\pi)\cos\Delta\alpha,\ \cos\delta - \rho\cos\varphi'\sin\pi\cos H\right),
$$

and α′ = α + Δα, reduced to [0°, 360°). The distance is kept.

Because Δα is proportional to sin H, it vanishes at H = 0: the topocentric and
geocentric transits are the same instant. The parallax shift acts almost
entirely on altitude, lowering the Sun by about π cos h.

### 11. Horizontal coordinates

The local frame is obtained by rotating the hour-angle frame about the
east–west axis by 90° − φ. With H′ = LST − α′ (wrapped to [−180°, 180°)), the
Sun's direction has components

$$
\begin{aligned}
\text{north} &= \sin\delta'\cos\varphi - \cos\delta'\cos H'\sin\varphi,\\
\text{east} &= -\cos\delta'\sin H',\\
\text{up} &= \sin\delta'\sin\varphi + \cos\delta'\cos H'\cos\varphi,
\end{aligned}
$$

so

$$
\sin h = \sin\varphi\sin\delta' + \cos\varphi\cos\delta'\cos H',\qquad
A = \operatorname{atan2}(\text{east},\ \text{north}) \pmod{360^\circ}.
$$

A is measured from North through East. sin h is clamped to [−1, 1]. The
latitude is geodetic, so "up" is the ellipsoid normal, which is the direction a
level instrument defines. This is the geometric (airless) topocentric position.

### 12. Atmospheric refraction

The atmosphere bends light downward, so a body appears higher than its
geometric altitude. The effect is about half a degree at the horizon and falls
to about 1′ at 45°.

**As implemented.** Given a geometric altitude h<sub>g</sub> (degrees),
pressure P (mbar) and temperature T<sub>C</sub> (°C):

1. If h<sub>g</sub> < −2°, no refraction is applied.
2. The standard-atmosphere refraction is scaled by

$$
k = \frac{P}{1010}\cdot\frac{283}{273 + T_C}.
$$

3. The Saemundsson (1986) expression (Meeus 16.4) gives the refraction as a
   function of the true (geometric) altitude:

$$
R(h_g) = k \cdot \frac{1.02'}{\tan\!\left(h_g + \dfrac{10.3^\circ}{h_g + 5.11^\circ}\right)},
\qquad h_a = h_g + R(h_g).
$$

   It is evaluated once, directly; no iteration is involved.

`Atmosphere::refraction = false` returns the geometric altitude.

With the default atmosphere (1010 mbar, 10 °C):

| Geometric altitude | Refraction |
|---|---|
| 0° | 28.98′ |
| 5° | 9.67′ |
| 10° | 5.41′ |
| 20° | 2.74′ |
| 45° | 1.01′ |

An apparent altitude of 0° corresponds to 34.43′ of refraction, in agreement
with the standard 34′ horizon refraction used by almanacs.

Refraction does not change the azimuth.

### 13. Events

#### 13.1 Definitions

| Event | Condition on the Sun's center |
|---|---|
| Sunrise, sunset | Apparent (refracted, topocentric) altitude = −16′ |
| Civil twilight | Apparent altitude = −6° |
| Nautical twilight | Apparent altitude = −12° |
| Astronomical twilight | Apparent altitude = −18° |
| Any event | Apparent altitude = any target, via `events::crossing` |
| Transit | Hour angle of the apparent Sun = 0 |

- **−16′** is the Sun's semidiameter, so the event is the upper limb touching
  the horizon. The semidiameter is fixed; the real value varies between about
  15.7′ (July) and 16.3′ (January) with the Earth–Sun distance.
- **Twilights** lie below −2°, where no refraction is applied
  ([§12](#12-atmospheric-refraction)), so they are geometric altitudes.
- **Transit** uses the geocentric hour angle. It is identical to the
  topocentric one, because the parallax shift in right ascension vanishes at
  H = 0 ([§10](#10-topocentric-parallax)).

**Relation to the almanac convention.** Almanacs define sunrise and sunset
as the geometric altitude of the center at −0.8333°: the 34′ of refraction at
the horizon plus the 16′ semidiameter, both taken as fixed. Spp instead
refracts the center itself, which sits 16′ below the horizon, where the
refraction is larger. With the default atmosphere its condition corresponds
to a geometric altitude of −0.8956° (37.73′ of refraction). So Spp's sunrise
comes earlier and its sunset later than the almanac convention, symmetrically
(measured, sea level, 2026):

| Place and date | Sunrise | Sunset |
|---|---|---|
| Equator, March equinox | −14.9 s | +14.9 s |
| 48.86° N (Paris), June solstice | −29.0 s | +29.0 s |
| 48.86° N (Paris), December solstice | −28.1 s | +28.1 s |
| 60° N, June solstice | −52.4 s | +52.4 s |

The difference grows with latitude because the Sun crosses the horizon at a
shallower angle. To reproduce the almanac convention exactly, solve for the
geometric altitude −0.8333° with refraction disabled:
`events::crossing(transit, obs, {.refraction = false}, -50.0 / 60.0, dir)`.

#### 13.2 One day

`events::day(date, obs, atm, dut1)`:

1. **First guess.** Local solar noon is taken as 12:00 − λ<sub>obs</sub>/15 h
   UTC on the given date, split into whole hours, minutes and seconds.
2. **Transit** is refined from it ([§13.3](#133-transit)).
3. **Classification.** h<sub>noon</sub> is the apparent altitude at transit,
   and h<sub>before</sub>, h<sub>after</sub> the apparent altitudes 12 h
   before and after it, approximating the two lower culminations.
   - h<sub>noon</sub> < −16′: `PolarNight`, duration 0 h.
   - h<sub>before</sub> > −16′ and h<sub>after</sub> > −16′: `PolarDay`,
     duration 24 h.
   - Otherwise `Normal`.
4. **Rise and set** are the crossings of −16′ before and after transit
   ([§13.4](#134-crossings)).
5. **Duration** is the time above the horizon within the 24 h centered on
   transit: set − rise, in hours, computed from the split JDs. Rise and set
   may fall on the UTC dates on either side of `date`.

**Days next to polar day.** On the last day before polar day, the Sun rises
but does not set again within 12 h of transit; on the first day after it, the
Sun sets but was already up 12 h before transit. Such a day is `Normal` with
only one of `rise` and `set`, and the missing end of the window counts as
daylight. For example, at 70° N on 2026-07-27 the Sun sets at 23:45 UTC with
no sunrise, and the duration is 23.64 h.

#### 13.3 Transit

The hour angle of the Sun advances at

$$
\frac{dH}{dt} = 360.9856^\circ/\text{day} - \frac{d\alpha}{dt} \approx 360^\circ/\text{day},
$$

since α increases by about 1° per day (0.9°–1.1° over the year); over 2026
the rate stays between 359.876° and 360.089° per day. Newton's method with
the fixed slope 360°/day,

$$
t \leftarrow t - \frac{H(t)}{360}\ \text{days},
$$

shrinks the error by a factor of about 3,000 or more per step, because the
slope is off by at most 0.035 %. The iteration stops when |H| < 10⁻⁹°
(≈ 0.2 µs of time), or when a step no longer reduces |H|, which is the
rounding noise floor, or after `max_iter` steps. From the first guess (off by
up to about 16 minutes, the equation of time) it typically converges in 3 or 4
steps.

#### 13.4 Crossings

`events::crossing` finds when the apparent altitude reaches a target
h<sub>t</sub> on one side of transit. With f(t) = h(t) − h<sub>t</sub>,
f ≥ 0 at transit is required; otherwise there is no crossing and the result
is empty.

**1. Analytic first guess.** From the altitude equation of
[§11](#11-horizontal-coordinates), the hour angle at which a body of
declination δ reaches altitude h₀ satisfies

$$
\cos H_0 = \frac{\sin h_0 - \sin\varphi\sin\delta}{\cos\varphi\cos\delta}.
$$

Here δ is the Sun's declination at transit, and the target is converted to a
geometric altitude by h₀ = h<sub>t</sub> − r, where r = refraction(h<sub>t</sub>)
− h<sub>t</sub> estimates the refraction at the event. If |cos H₀| < 1, the
guess is t<sub>transit</sub> ± H₀/360 days, with H₀ in degrees: minus before
transit, plus after.

**2. Bracket.** f is evaluated at the guess.

- If f ≥ 0 (not yet past the crossing), the search steps away from transit by
  2 minutes, doubling the step each time, until f < 0. It stops if it would go
  more than half a day from transit.
- If f < 0 (already past), it steps back toward transit the same way until
  f ≥ 0, and uses transit itself if it gets there, since f ≥ 0 is known there.

The guess is usually within a few minutes, so this takes one or two
evaluations.

**3. Fallback.** When there is no usable guess (|cos H₀| ≥ 1, near the polar
limits) or the walk fails, the bracket is found by scanning outward from
transit in 1-hour steps, up to 12 hours.

**4. Refinement.** The Illinois method (a modified regula falsi) keeps a
bracket [a, b] with f(a) ≥ 0 > f(b):

$$
x = \frac{a\,f(b) - b\,f(a)}{f(b) - f(a)} .
$$

x replaces the endpoint with the same sign of f. When the same endpoint is
replaced twice in a row, the function value stored at the other endpoint is
halved. This prevents the one-sided stagnation of plain regula falsi and gives
superlinear convergence (order about 1.44) while always keeping the root
bracketed. If x falls outside the bracket, the midpoint is used instead. The
iteration stops when the bracket is narrower than
`SolverConfig::time_tolerance_s` (default 1 µs) or after `max_iter`
iterations, and returns the bracket's midpoint.

Measured on random dates and observers, the refined instants satisfy the event
condition to better than 10⁻⁵″ in altitude.

---

## Numerical methods

### Chebyshev earth state

Evaluating VSOP87D (2,425 cosines) and the nutation series (20 sines and
cosines, 5 arguments) dominates the cost of a solar position. All five
quantities, L, B, R, Δψ and Δε, are smooth functions of time, so Spp replaces
them with Chebyshev interpolants over short segments.

**Segments.** Time is divided into 16-day segments aligned on J2000.0: segment
k covers days [16k, 16k + 16) from JD 2451545. For an instant (d, f), the
segment is k = ⌊(d − 2451545)/16⌋, taken from the integer day alone since
f ∈ [0, 1). The local coordinate is formed from the split time:

$$
x = \frac{2}{16}\big((d - 2451545 - 16k) + f\big) - 1 \in [-1, 1).
$$

**Fit.** Each segment is built on first use from N + 1 = 25 direct
evaluations at the Chebyshev nodes

$$
x_j = \cos\frac{\pi (j + \tfrac12)}{N + 1},\qquad j = 0,\dots,N,
$$

and each quantity's degree-N (N = 24) interpolant has coefficients

$$
c_m = \frac{2}{N + 1}\sum_{j=0}^{N} g(x_j)\cos\frac{\pi m (j + \tfrac12)}{N + 1},\qquad c_0 \leftarrow \tfrac12 c_0 ,
$$

so that $g(x) \approx \sum_{m=0}^{N} c_m T_m(x)$.

**Offset.** Each quantity is fitted as the residual g = X − X<sub>mid</sub>
from its value at the middle node. VSOP87D longitude reaches about
2 × 10⁴ rad, and fitting it directly would cost about 10 bits of precision
to rounding in the coefficients. L advances about 0.28 rad over a segment, so
the residual from the middle node is at most about 0.14 rad. The offset is
added back as the last operation, with a single rounding.

**Evaluation.** Clenshaw's recurrence, run for the five quantities in the
same loop so that their independent dependency chains overlap:

$$
b_{N+1} = b_{N+2} = 0,\quad b_m = 2x\,b_{m+1} - b_{m+2} + c_m\ (m = N,\dots,1),\quad
g(x) = x\,b_1 - b_2 + c_0 .
$$

**Why degree 24 over 16 days is enough.** For a single term cos(ω t) over a
segment of half-length h, the Chebyshev coefficients decay like the Bessel
functions J<sub>m</sub>(ωh), which fall off faster than exponentially once
m > ωh. Taking the worst cases for h = 8 days:

- The fastest VSOP87D Earth term (period 7.0 days, amplitude ≤ 10⁻⁸ rad) has
  ωh ≈ 7.2. Its coefficient beyond degree 24 is of order 10⁻¹⁹ rad.
- Terms above 10⁻⁸ rad have periods of at least 9.4 days (ωh ≤ 5.4) and
  decay even faster.
- The fastest nutation term (9.1 days, 0.005″) is negligible after degree 24.

The measured error confirms it (see [Accuracy](#accuracy)). The fit is exactly
as close to the true series as direct double-precision evaluation, and the
closer of the two varies at random. Longer segments (32 days) were measured to
start adding error in B.

**Cache.** Segments are kept in a 64-slot direct-mapped table, slot
k mod 64, of 1,048 bytes each (67 KB), declared `thread_local`. That covers
1,024 consecutive days (about 2.8 years) without eviction. The pipeline
requests the state of the same instant twice (for nutation, then for the
Sun), so the most recent instant is also memoized.

### Other optimizations

- **Split GMST** removes the 10⁸° intermediate exactly
  ([§9](#9-sidereal-time-and-hour-angle)).
- **Horner form** for the obliquity polynomial, instead of ten `pow` calls
  ([§6](#6-obliquity-of-the-ecliptic)).
- **Direct inversion** λ = L + 180°, β = −B instead of a round trip through
  rectangular coordinates ([§4](#4-geometric-geocentric-sun)).
- **Shared nutation.** Δψ and ε are computed once per instant and used for
  both the apparent place and the equation of the equinoxes.
- **Solver.** The analytic first guess brackets each crossing within minutes,
  and the superlinear refinement converges to 1 µs ([§13.4](#134-crossings)).

None of these changes a result by more than double-precision rounding; the
split GMST is more accurate than the direct form.

---

## Accuracy

### Error budget

| Source | Size | Affects |
|---|---|---|
| VSOP87D theory (complete series) | About 1″ over several thousand years around J2000 (authors' estimate) | Solar longitude, latitude, distance |
| Dynamical vs FK5 equinox | 0.09″, not applied ([Limitations](#limitations)) | Solar longitude |
| Nutation: IAU 1980 truncated to 20 terms | A few hundredths of an arcsecond | Apparent place, sidereal time |
| Aberration: none in latitude, no barycentric motion of the Sun | Below 0.01″ | Apparent place |
| Equation of the equinoxes without the 1994 terms | Below 0.003″ | Sidereal time |
| ΔT, observed era (2005–2025) | About 0.2 s → 0.01″ in position | Solar position |
| ΔT, linear bridge after 2026 ([§2.6](#26-δt)) | About 0.7 s high by late 2026 → 0.03″ | Solar position |
| ΔT, future and pre-telescopic | Seconds to hours → 0.041″ per second in position | Solar position |
| DUT1 = 0 (default) | Up to 0.9 s → up to 13.5″ | Hour angle, altitude, azimuth, event times |
| Sunrise/sunset definition vs the almanac convention ([§13.1](#131-definitions)) | Sunrise 15–52 s earlier (sunset later) at 0°–60° latitude | Rise/set |
| Fixed 16′ semidiameter | ±0.3′ → a few seconds | Rise/set |
| Real atmosphere vs model | Often minutes near the horizon | Rise/set, low altitudes |
| Numerical error | See below | Everything |

### Numerical error

Numerical error is far below the model error:

- **Cache.** Against an extended-precision (80-bit) evaluation of the same
  series over 1000 BC–3000 AD, both the cache and direct double evaluation are
  within:

  | Quantity | Direct (double) | Cached |
  |---|---|---|
  | L | 1.2e-5″ | 1.2e-5″ |
  | B | 3.3e-11″ | 3.3e-11″ |
  | R | 1.1e-13 AU | 9.0e-14 AU |
  | Δψ | 1.8e-10″ | 2.1e-10″ |
  | Δε | 8.4e-11″ | 8.2e-11″ |

  The remaining error of both is the rounding of phases B + C·t, which reach
  about 10⁶ rad at the ends of the range. The larger figure for L reflects
  that L itself reaches 2 × 10⁴ rad.
- **GMST** is within 1.9e-7″ of the formula evaluated in exact rational
  arithmetic ([§9](#9-sidereal-time-and-hour-angle)).
- **Events.** Transits converge to below 10⁻⁹° of hour angle (~0.2 µs). Rise,
  set and twilight instants converge to 1 µs, and their altitude residuals are
  below 10⁻⁵″.

---

## Performance

Measured with `spp_bench` on one core (release build). Each measurement runs
once before it is timed, so all of them see the ephemeris cache already
filled for the instants they use:

| Operation | Time |
|---|---|
| `sun::apparent` | 0.15 µs |
| `events::altitude_deg` (includes the apparent Sun) | 0.35 µs |
| `events::day` (transit, rise, set) | 6.6 µs |
| `spp year` (365 days, including process start-up) | ~60 ms |

Filling one 16-day cache segment costs 25 direct series evaluations, about
0.3 ms. How the rest is achieved is described in
[Numerical methods](#numerical-methods).

---

## Accuracy tools

| Target | What it checks |
|---|---|
| `spp_cache_check` | The Chebyshev cache and direct double evaluation, both against the same series evaluated in 80-bit extended precision, at random instants over 1000 BC–3000 AD, densely inside segments, and on both sides of segment boundaries. The cache passes when it is no further from the truth than direct evaluation |
| `spp_bench` | Throughput of `sun::apparent`, `events::altitude_deg` and `events::day` over a year × latitude grid |
| `spp_figures` | Renders the figures in `docs/images/` as SVG from the library's results: daylight and twilight, day length by latitude, analemma, sun paths, azimuth spirals |

```sh
./build/release/spp_cache_check
cmake --preset release -DSPP_BUILD_BENCH=ON && cmake --build --preset release && ./build/release/spp_bench
cmake --preset release -DSPP_BUILD_FIGURES=ON && cmake --build --preset release && ./build/release/spp_figures docs/images
```

---

## Project layout

```
Seasonal++/
├── CMakeLists.txt, CMakePresets.json
├── README.md, LICENSE
├── .editorconfig, .gitattributes, .gitignore
├── cmake/SppFlags.cmake          compiler flags and floating-point policy
├── data/VSOP87D.ear              VSOP87D Earth series (source data)
├── tools/
│   ├── vsop_gen.cpp              .ear → exact C++ tables, run at build time
│   ├── cache_check.cpp           cache accuracy against extended precision
│   ├── bench.cpp                 benchmark
│   └── figures.cpp               renders the README figures (SVG)
├── docs/images/                  README figures
├── include/spp/
│   ├── spp.hpp                   umbrella header
│   ├── core/                     Time, Instant, calendar, angles, constants, types
│   ├── sun.hpp  earth.hpp  sky.hpp  events.hpp
│   └── model/                    vsop87d, iau1980, laskar, meeus, espenak_meeus, saemundsson
├── src/
│   ├── core/  model/             implementations
│   ├── sun.cpp  earth.cpp  sky.cpp  events.cpp
│   └── earth_state.{hpp,cpp}     Chebyshev earth-state cache
└── apps/
    ├── spp/main.cpp              command-line tool (demonstration front end)
    ├── spiral/main.cpp           sunrise/sunset azimuth spirals (OpenGL)
    └── common/                   helpers shared by the viewer and the figures
```

---

## Limitations

- **The Sun only.** Moon and planets are not implemented. The `spp::model`
  layout is designed to take them, e.g. further VSOP87 bodies.
- **Nutation.** The truncated IAU 1980 series limits the apparent place to a
  few hundredths of an arcsecond. IAU 2000A/B would be the upgrade path.
- **Frame.** Positions are referred to the dynamical ecliptic and equinox of
  date (VSOP87D). The small rotation to the FK5 system (Meeus 32.3:
  ΔL = −0.09033″, ΔB = 0.03916″ (cos L′ − sin L′) with
  L′ = L − 1.397°T − 0.00031°T²) is not applied. There is no ICRS/GCRS
  pipeline, precession-nutation matrix model or CIO-based sidereal time.
- **Sidereal time** omits the 1994 complementary terms of the equation of the
  equinoxes (below 0.003″).
- **ΔT** always comes from the built-in model; the API has no way to supply
  another value. The observed table ends with a provisional value for 2026,
  and later years are predictions; the 2026–2050 bridge is already about
  0.7 s above the observed value ([§2.6](#26-δt)). **DUT1** can be supplied
  per call.
- **Refraction** is a standard-atmosphere model. It applies no refraction
  below −2° geometric altitude, and real horizon refraction varies from night
  to night.
- **Semidiameter** is fixed at 16′ instead of following the Earth–Sun
  distance.
- **Day classification** samples the altitude 12 h before and after transit
  instead of finding the exact lower culminations ([§13.2](#132-one-day)).
- **Rise/set with a raised horizon** (mountains, observer height above
  terrain) is not modelled; elevation only affects parallax.

---

## References

- P. Bretagnon, G. Francou, *Planetary theories in rectangular and spherical
  variables: VSOP87 solutions*, Astronomy & Astrophysics 202, 309 (1988).
- J. Meeus, *Astronomical Algorithms*, 2nd ed., Willmann-Bell (1998):
  chapters 7 (Julian Day), 12 (sidereal time), 13 (transformations),
  16 (refraction), 22 (nutation and obliquity), 25 (solar coordinates),
  32 (VSOP87 and FK5), 40 (parallax).
- J. Laskar, *Secular terms of classical planetary theories using the results
  of general theory*, Astronomy & Astrophysics 157, 59 (1986).
- P. K. Seidelmann (ed.), *Explanatory Supplement to the Astronomical Almanac*
  (1992): IAU 1980 nutation, IAU 1982 GMST.
- F. Espenak, J. Meeus, *Five Millennium Canon of Solar Eclipses*, NASA
  TP-2006-214141 (2006): polynomial expressions for ΔT.
- T. Saemundsson, *Astronomical refraction*, Sky and Telescope 72, 70 (1986).
- IERS / USNO: observed values of ΔT and DUT1.
- NIMA TR8350.2, *Department of Defense World Geodetic System 1984*: WGS84
  ellipsoid.

---

## License

GNU General Public License v3.0 or later. See [LICENSE](LICENSE).
