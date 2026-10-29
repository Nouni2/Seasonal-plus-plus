/*
 * spp_c.h — the C interface of Spp (Seasonal++), for use from other languages
 * through a shared library (spp.dll on Windows, spp_capi.so elsewhere).
 *
 * The interface carries no C++ types and lets no exception cross it. Every
 * function returns an int32_t status (SPP_OK or an SPP_E_* code) and writes
 * its results through pointers.
 *
 * Instants are UTC, as integer milliseconds since 1970-01-01T00:00:00Z (Unix
 * time: every day has 86,400 s, as in Spp's UTC). Results are rounded to the
 * nearest millisecond. The scale is part of every name (*_utc_ms). Supported
 * years: -1000 to 3000 (astronomical numbering; Julian calendar before
 * 1582-10-15, Gregorian from then on).
 *
 * Physical assumptions, fixed in this version of the interface:
 *   - observer on the WGS84 ellipsoid at the given elevation;
 *   - refraction: Spp's standard atmosphere (1010 mbar scaled to the elevation,
 *     10 degrees C);
 *   - DUT1 = UT1 - UTC estimated by Spp from its observed Delta T table.
 *
 * Definitions are those of the Spp README (section 13.1): sunrise and sunset
 * are the instants the refracted upper limb of the Sun is at 0 degrees; civil
 * dawn and dusk the instants its center is at -6 degrees; transit the instant
 * its geocentric apparent hour angle is 0.
 *
 * All functions are thread-safe.
 */

#ifndef SPP_C_H
#define SPP_C_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32)
#  if defined(SPP_C_BUILD)
#    define SPP_C_API __declspec(dllexport)
#  else
#    define SPP_C_API __declspec(dllimport)
#  endif
#  define SPP_C_CALL __cdecl
#elif defined(__GNUC__)
#  define SPP_C_API __attribute__((visibility("default")))
#  define SPP_C_CALL
#else
#  define SPP_C_API
#  define SPP_C_CALL
#endif

/* Version of this interface: the layout of the structures and the meaning of
 * the functions. It changes whenever either changes. */
#define SPP_C_ABI_VERSION 1

/* Status codes. */
#define SPP_OK 0
#define SPP_E_ARGUMENT 1 /* null pointer, wrong struct_size, invalid date, non-finite or out-of-range value */
#define SPP_E_RANGE 2    /* an instant or date outside the years -1000 to 3000 */
#define SPP_E_INTERNAL 3 /* any other failure inside the library */

/* Day types (spp::events::DayType). */
#define SPP_DAY_NORMAL 0      /* not polar day or night; sunrise or sunset may still be missing (see below) */
#define SPP_DAY_POLAR_DAY 1   /* the upper limb stays above the horizon over the 24 h window */
#define SPP_DAY_POLAR_NIGHT 2 /* the upper limb stays below the horizon over the 24 h window */

/* An observer. Latitude geodetic, positive North, in [-90, 90]; longitude
 * positive East, in [-180, 180]; elevation in metres above the WGS84 ellipsoid,
 * in [-1000, 11000]. */
typedef struct spp_observer {
    double lat_deg;
    double lon_deg;
    double elev_m;
} spp_observer;

/* One solar day: the 24 hours centered on the transit that Spp finds from
 * 12:00 - longitude/15 h UTC on the given UTC date (README section 13.2).
 * Events can fall on the UTC date before or after it.
 *
 * Each event has a presence flag (0 or 1); its time field is 0 when the flag
 * is 0. A Normal day can lack sunrise or sunset: next to polar day the Sun is
 * already up at the start (or still up at the end) of the window, and within
 * about 0.1 degree of a pole it can set before transit or rise after it.
 * Civil dawn and dusk are absent when the Sun's center does not reach -6
 * degrees within 12 h of transit (it stays above, or stays below).
 *
 * duration_hours is the time the upper limb is above the horizon within the
 * window: 24 on polar day, 0 on polar night. */
typedef struct spp_solar_day {
    int32_t struct_size; /* set by the caller to sizeof(spp_solar_day) */
    int32_t day_type;    /* SPP_DAY_* */
    int32_t has_rise;
    int32_t has_set;
    int32_t has_civil_dawn;
    int32_t has_civil_dusk;
    int64_t transit_utc_ms;
    int64_t rise_utc_ms;
    int64_t set_utc_ms;
    int64_t civil_dawn_utc_ms;
    int64_t civil_dusk_utc_ms;
    double duration_hours;
} spp_solar_day;

/* SPP_C_ABI_VERSION of the library actually loaded. */
SPP_C_API int32_t SPP_C_CALL spp_c_abi_version(void);

/* Library version and source revision, e.g. "1.0.0+f609065". Static storage. */
SPP_C_API const char* SPP_C_CALL spp_c_library_version(void);

/* sizeof(spp_solar_day) as compiled into the library, for layout checks. */
SPP_C_API int32_t SPP_C_CALL spp_c_solar_day_size(void);

/* The solar day for a UTC civil date (year, month 1-12, day 1-31). */
SPP_C_API int32_t SPP_C_CALL spp_solar_day_utc(const spp_observer* observer, int32_t year, int32_t month,
                                                int32_t day, spp_solar_day* out);

/* Apparent (refracted) altitude of the Sun's upper limb, in degrees, at a UTC
 * instant. Positive: the Sun is up, by the same definition as sunrise and
 * sunset. */
SPP_C_API int32_t SPP_C_CALL spp_upper_limb_altitude(const spp_observer* observer, int64_t utc_ms,
                                                     double* out_alt_deg);

/* The same at count instants start_utc_ms + i * step_ms, i = 0 .. count - 1,
 * written to out_alt_deg[0 .. count - 1]. count in [0, 1000000]. On an error
 * nothing is guaranteed about out_alt_deg. */
SPP_C_API int32_t SPP_C_CALL spp_upper_limb_altitudes(const spp_observer* observer, int64_t start_utc_ms,
                                                      int64_t step_ms, int32_t count, double* out_alt_deg);

#ifdef __cplusplus
}
#endif

#endif /* SPP_C_H */
