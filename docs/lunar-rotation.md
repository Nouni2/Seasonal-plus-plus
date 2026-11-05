# Lunar orientation

The rotation provider is separate from the lunar orbit provider. It evaluates
the IAU2009 trigonometric approximation of the Moon's Mean Earth/Polar Axis
frame, named `IAU_MOON` in NAIF's text PCK. Orbital ELP fits do not change it.
It supplies the ICRF pole, prime meridian, their analytic angular rates, and a
typed rotating-frame transform with its time derivative.

This is an approximate orientation model. It is not the high-precision DE441
physical-libration solution. Those data can be supplied by a future rotation
provider returning the same frame contract. The GUI's **Libration grid** / **L**
option shows this body-fixed grid and sub-observer coordinates. It has no terrain map.

## Explore it

```bash
cmake --preset moon
cmake --build --preset moon --target spp_cli
./build/moon-Linux/spp moon rotation 2026-11-03 12:00
```

For the Windows build, from the same project root in MSYS2 UCRT64:

```bash
./build/moon-Windows/spp.exe moon rotation 2026-11-03 12:00
./build/moon-Windows/spp.exe moon surface 2026-11-03 12:00
./build/moon-Windows/spp.exe moon surface 2026-11-03 12:00 48.8566 2.3522 35
```

This command takes **TDB** calendar fields. It prints the pole's ICRF right
ascension and declination, prime-meridian angle W, angular rates in degrees
per TDB day, and the body axes' ICRF components. Its output identifies the
approximation. `--fit` belongs to orbital commands and is rejected here.

The raw rotation is evaluated at the supplied physical instant. An observer
application must evaluate it at lunar reflection time, accounting for light
travel. It must also transform mean-ecliptic ELP vectors explicitly into ICRF.
The `surface` commands take **UTC** and perform those operations, geocentrically
or for the supplied observer (the example is Paris).

## C++ API and frames

```cpp
#include <spp/spp.hpp>
#include <spp/model/iau_moon.hpp>

const auto time = spp::Time::from_utc(2026, 11, 3, 12);
const spp::model::iau_moon::Rotation model;
const auto orientation = model(time.tdb());
const auto& rotation = orientation.rotation;
```

`orientation` contains `pole_ra_deg`, `pole_dec_deg`, `prime_meridian_deg`
and their `*_rate_deg_day` values. `rotation` has type
`FrameRotation<frame::Icrf, frame::MoonMeanEarth>`. TT inputs and incompatible
state frames are rejected at compile time.

`FrameRotation` preserves the origin tag, rotates typed Cartesian positions,
and transforms a state as **p' = R p**, **v' = R v + (dR/dt) p**. Rates and
state velocities use TDB days. Inversion and composition retain frame tags
and rotation derivatives. All operands must describe the same epoch; origin
translation is a separate operation. The IAU model evaluates without runtime
files, allocation or shared state.

The explicit fixed adapter `frames::mean_ecliptic_j2000_to_icrf()` includes the
J2000 frame bias and obliquity. Its inverse converts back; neither changes origin.

## Libration and surface geometry

```cpp
spp::model::elp_mpp02::CachedEphemeris orbit;
spp::model::iau_moon::Rotation spin;
const spp::Observer paris{48.8566, 2.3522, 35};
const auto surface = spp::moon::surface(time, paris, orbit, spin);
const auto observation = spp::moon::observe(time, paris, orbit, spin);
```

`Surface` provides `subobserver_lon_deg`, `subobserver_lat_deg`,
`subsolar_lon_deg`, `subsolar_lat_deg`, `pole_position_angle_deg` and
`pole_defined`. Longitudes are east-positive in [-180,180); latitudes are
planetocentric. The sub-observer coordinates describe libration of the point
at the apparent disk center. The subsolar coordinates use the incoming
aberrated Sun at reflection. The pole angle is measured east from celestial north.

`body_to_view` projects body-fixed unit directions into screen right,
celestial north and towards the viewer. It is an orthographic sphere orientation,
including the conformal orientation change from observer aberration. It omits
terrain, differential surface light travel and atmospheric image distortion.

The rotation provider is evaluated once at reflection TDB and may return a
typed `BodyRotation` (by value or reference) or an object containing `.rotation`.
It is borrowed synchronously and never copied. Invalid/non-finite or reflected
frames are rejected. An independent provider can replace the IAU model without
changing orbital evaluation. Default overloads use IAU rotation and the cached
ELP orbit. `Observation.surface` shares the observation's ray calculations;
it adds no orbital provider calls. The lean coordinate API omits surface work.

Its numerical evaluation bounds are -31 through +11 Julian centuries from
J2000, matching the current lunar engine. These bounds are not an accuracy
claim for ancient or future orientation. Non-finite/unusable split dates and
out-of-domain inputs fail explicitly.

## References and measured limits

Coefficients come from [NAIF pck00011.tpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc),
the unchanged IAU2009 lunar pole/prime-meridian series. The PCK's SHA-256 is
`3dff7b1dbeceaa01f25467767d3fa25816051c85d162d1edf04acb310ee28bb1`.
The reference generator independently parses that source, evaluates its
phases with decimal arithmetic, and constructs the axes from a Cartesian pole
basis. Eight epochs include the numerical endpoints. Separate finite
differences check matrix rates and moving-state velocities; inverse and
composition checks include the velocity terms.

Seven independent [Horizons](https://ssd.jpl.nasa.gov/horizons/manual.html)
cases from 1900 through 2100 use DE441's **MOON_ME** rotation. Its ICRF pole,
astrometric sight line and sub-observer longitude recover an independent
rotation matrix; sub-observer latitude supplies a consistency control.
Reception TT is shifted by the reported one-way light time and converted to
TDB before evaluating the model.

The sampled IAU/DE441 difference was at most **0.002807 degree** (about
10.1 arcseconds) for the complete axes, and 7.70 arcseconds for the pole.
This is not a uniform physical bound. Implementation differences against the
separately evaluated PCK were below 5e-10 degree; matrix-rate differences were
below 1.8e-12 per day. Model disagreement and numerical agreement are separate.

The NAS GCC 14.2 Release benchmark measured about **0.78 microsecond** per
orientation plus derivative, for 4096 dates and nine batches. No speed
threshold is used. Reproduce the benchmark or regenerate only the references:

```bash
cmake --build --preset moon --target spp_moon_rotation_bench
./build/moon-Linux/spp_moon_rotation_bench tools/reference/moon_iau_rotation.csv tools/reference/moon_rotation_horizons.csv
python3 tools/lunar_rotation_reference.py --cache /tmp/spp-moon-observation-reference
```

Normal builds and checks use the committed CSVs without Python or network.
Reference regeneration does not run Spp and preserves unmodified downloads
under their request hashes. It replaces a fixture only after successful
generation. The complete Moon report includes these rotation gates.

Thirteen independent Horizons surface cases cover geocentric observations
from 1900–2100 and topocentric observations from Paris and Sydney. Sampled
IAU/DE441 maximum differences are 0.002769 degree for the sub-observer point,
0.002556 degree for the subsolar point and 0.001114 degree for pole bearing.
Injecting independently reconstructed DE441 rotations isolates ray/frame
geometry: its maximum geocentric observer/Sun residual is **0.283 arcsecond**.
The default resident surface operation measured about **3.2 microseconds/query**
without lunar cache refills. These sampled results are not global error bounds.
