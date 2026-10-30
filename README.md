<p align="center">
  <img src="docs/images/banner.svg" alt="Spp (Seasonal++): high-precision solar astronomy in C++20" width="50%">
</p>

<p align="center">
  <a href="https://github.com/Nouni2/Seasonal-plus-plus/releases/tag/v1.0.1"><img src="https://img.shields.io/badge/version-1.0.1-F2A93B?style=flat-square" alt="Version 1.0.1"></a>
  <a href="#building"><img src="https://img.shields.io/badge/C%2B%2B-20-00599C?style=flat-square&logo=cplusplus&logoColor=white" alt="C++20"></a>
  <a href="#building"><img src="https://img.shields.io/badge/CMake-3.28%2B-064F8C?style=flat-square&logo=cmake&logoColor=white" alt="CMake 3.28+"></a>
  <a href="#accuracy--precision"><img src="https://img.shields.io/badge/ephemeris-VSOP87D-E07B39?style=flat-square" alt="VSOP87D ephemeris"></a>
  <a href="#floating-point-policy"><img src="https://img.shields.io/badge/fp64-strict%20IEEE%20754-6E40C9?style=flat-square" alt="Strict IEEE 754 fp64"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-GPL--3.0-blue?style=flat-square" alt="License: GPL-3.0-or-later"></a>
  <a href="http://github.com/Nouni2/Seasonal-plus-plus/wiki"><img src="https://img.shields.io/badge/docs-wiki-0E8A16?style=flat-square" alt="Wiki"></a>
</p>

<p align="center">
  <b>Where the Sun is, and when it rises, transits and sets: anywhere, from 1000 BC to 3000 AD.</b><br>
  Full VSOP87D ephemeris · strict IEEE 754 fp64 · ~0.5 µs per sky position (warm cache) · events solved to 1 µs
</p>

Spp is a C++20 library for solar positional astronomy: apparent geocentric coordinates, topocentric refracted altitude and azimuth, sidereal time, nutation, ΔT, and solar event instants from 1000 BC to 3000 AD.

The library is the engine. All computation is exposed through `include/spp/`. The `spp` CLI front-end prints a subset of these results.

Built for batch evaluation across dates and observer grids. Calculations are verified against an 80-bit extended-precision reference to within double-precision rounding.

```powershell
$ spp day 2026-06-21 48.8566 2.3522 35
type      NORMAL
transit   2026-06-21T11:52:24Z
sunrise   2026-06-21T03:46:57Z
sunset    2026-06-21T19:57:52Z
duration  16.1820 h
```

---

## Features

- **Coordinates**: Apparent geocentric equatorial (RA, Dec, distance) and topocentric horizontal (true and refracted altitude, azimuth).
- **Events**: Root-finders for transit, twilight boundaries (civil, nautical, astronomical), and sunrise/sunset (refracted upper limb touching horizon). Polar day and polar night handling.
- **Astronomical Models**:
  - Earth heliocentric orbit: **VSOP87D** (compiled in, 2,425 terms).
  - Nutation: **IAU 1980** (106 periodic terms).
  - Sidereal time: **IAU 1982** GMST and GAST in split precision.
  - Earth orientation & ΔT: **IERS EOP 20 C04** (1972–2026 monthly table) and Espenak–Meeus polynomials.
  - Atmospheric refraction: Sæmundsson model with local elevation, temperature, and pressure scaling.
- **Throughput**: ~0.5 µs per position on warm Chebyshev cache; event solutions converge to 1 µs.

---

## Building

Requires a C++20 conforming compiler (GCC 14+, Clang 18+, or MSVC 2022) and CMake 3.28+.

```bash
# Configure and build the library and command-line tools
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Using presets (Ninja + GCC/UCRT64):

```bash
cmake --preset release
cmake --build --preset release
```

### CMake Integration

Spp installs standard CMake package files:

```cmake
find_package(Spp 1.0 CONFIG REQUIRED)
target_link_libraries(your_app PRIVATE spp::spp)
```

---

## C++ Example

```cpp
#include <iostream>
#include <spp/spp.hpp>

int main() {
    using namespace spp;

    // Define observer: Paris (48.8566° N, 2.3522° E, 35 m elevation)
    const Observer obs{
        .lat_deg = 48.8566,
        .lon_deg = 2.3522,
        .elev_m  = 35.0
    };

    // Construct a UTC timestamp
    const Time t = Time::from_utc(2026, 6, 21, 12, 0, 0.0);

    // 1. Topocentric solar coordinates (apparent position with atmospheric refraction)
    const Horizontal sky = sun::horizontal(t, obs);
    std::cout << "Altitude: " << sky.alt_deg << "°\n";
    std::cout << "Azimuth:  " << sky.az_deg  << "°\n";

    // 2. Solar events for the date
    const events::Day day = events::day({2026, 6, 21}, obs);
    if (day.type == events::DayType::Normal && day.rise && day.set) {
        std::cout << "Transit: " << day.transit.iso_utc() << "\n";
        std::cout << "Sunrise: " << day.rise->iso_utc()    << "\n";
        std::cout << "Sunset:  " << day.set->iso_utc()     << "\n";
    }

    return 0;
}
```

---

## Accuracy & Precision

### Floating-Point Policy
Computations run in IEEE 754 binary64 (`double`). With GCC and Clang, Spp builds with `-fno-fast-math` and `-ffp-contract=off`, so the compiler never fuses `a*b+c` into an FMA. MSVC builds use `/fp:precise` (untested).

### Verification
- **Solar Position**: Topocentric altitude and azimuth within 0.064″ of JPL DE421 at 8 sites, 78° S to 78° N (2000–2026).
- **Nutation & Sidereal Time**: IAU 1980 nutation within 0.0008 mas of ERFA `nut80` and IAU 1982 GMST within 0.5 mas of ERFA `gmst82` (1900–2100).
- **Solar Events**: Sunrise/sunset instants agree with DE421 geometry to ≤ 11 ms (latitudes 0°–64°) and ≤ 6.4 × 10⁻⁶″ of refracted upper limb altitude.
- **Full Pipeline**: Verified against an 80-bit extended-precision evaluation (`tools/pipeline_check.cpp`).

---

## Performance

| Operation | Warm cache | Cold cache |
|---|---|---|
| `sun::apparent` | 0.18 µs | |
| `events::altitude_deg` (includes the apparent Sun) | 0.48 µs | 0.48 ms |
| `events::day` (transit, rise, set) | 7.76 µs | 0.52 ms |

Warm: each query's 16-day Chebyshev segment is already filled (consecutive instants, or a year of days at 19 latitudes for `events::day`). Cold: random instants or dates over 1000 BC–3000 AD, timed on first use, so nearly every query fills a new segment.

*Measured with `spp_bench` (`tools/bench.cpp`) on one thread: AMD Ryzen 7 7840HS, Windows 11 x86-64, GCC 15.2.0, `release` CMake preset.*

---

## Limitations

- **Solar body only**: Computes solar coordinates and events. Lunar and planetary ephemerides are out of scope.
- **Platform `libm`**: While compiler code generation is deterministic, transcendental functions (`sin`, `cos`, `atan2`) depend on the host C runtime and may differ by the last ULP.
- **Pre-1600 ΔT**: Historical ΔT uncertainty dominates timing precision prior to 1600 AD.

---

## Extended Documentation (Wiki)

In-depth astronomical derivations, tutorials, and full API references are hosted on the [Project Wiki](http://github.com/Nouni2/Seasonal-plus-plus/wiki):

- [Getting Started & Installation Guide](http://github.com/Nouni2/Seasonal-plus-plus/wiki/Getting-Started)
- [C++ API Reference](http://github.com/Nouni2/Seasonal-plus-plus/wiki/Library-API)
- [Solar Position Pipeline](http://github.com/Nouni2/Seasonal-plus-plus/wiki/Solar-Position)
- [Event Solvers & Day Types](http://github.com/Nouni2/Seasonal-plus-plus/wiki/Events)
- [Time Scales & ΔT Calibration](http://github.com/Nouni2/Seasonal-plus-plus/wiki/Time-and-Delta-T)
- [Validation Tools & Audit Reproduction](http://github.com/Nouni2/Seasonal-plus-plus/wiki/Validation-Tools)
- [Known Issues & Workarounds](http://github.com/Nouni2/Seasonal-plus-plus/wiki/Known-Issues)
- [Roadmap](http://github.com/Nouni2/Seasonal-plus-plus/wiki/Roadmap)

---

## Contributing

Please see [CONTRIBUTING.md](CONTRIBUTING.md) for precision policies, development standards, and test verification requirements.

---

## License

Spp is licensed under the [GNU General Public License v3.0 or later](LICENSE).

