#pragma once

namespace spp {

// Ecliptic of date: longitude λ, latitude β, radius vector.
struct Ecliptic {
    double lon_deg;
    double lat_deg;
    double dist_au = 1.0;
};

// Equatorial of date: right ascension α, declination δ.
struct Equatorial {
    double ra_deg;
    double dec_deg;
    double dist_au = 1.0;

    double ra_hours() const { return ra_deg / 15.0; }
};

// Local horizon. Azimuth is measured from North, increasing East.
struct Horizontal {
    double alt_deg;
    double az_deg;
};

struct Observer {
    double lat_deg;        // geodetic, positive North
    double lon_deg;        // positive East
    double elev_m = 0.0;   // above the WGS84 ellipsoid
};

struct Atmosphere {
    double pressure_mbar = 1010.0;
    double temp_c = 10.0;
    bool refraction = true;
};

}  // namespace spp
