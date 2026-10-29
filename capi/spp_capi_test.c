/*
 * spp_capi_test: checks the C interface (spp/spp_c.h) through the shared
 * library, in C.
 *
 *   - Layout and version: SPP_C_ABI_VERSION and sizeof(spp_solar_day) agree
 *     with the library.
 *   - Error codes for invalid input: null pointers, wrong struct_size, invalid
 *     dates, non-finite and out-of-range values, out-of-range instants.
 *   - Reference results: solar days for fixed places and dates, against values
 *     computed with Spp's C++ API on Linux x86-64 (GCC 14.2) at the revision
 *     the library was built from. Events must agree within 1 ms (the results
 *     are rounded to the millisecond, and another platform's math library can
 *     move one across a rounding boundary), durations within 1e-6 h.
 *   - Definitions: at every reported sunrise and sunset the upper limb is at
 *     the horizon (|altitude| <= 1e-5 degrees after rounding to the
 *     millisecond), rising at sunrise and setting at sunset; events are in
 *     order; the batch and single altitude functions agree exactly.
 *
 *   spp_capi_test           runs the checks, prints PASS or FAIL, exits with 0 or 1
 *   spp_capi_test --print   prints the reference table as computed by this library
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spp/spp_c.h"

typedef struct {
    const char* name;
    spp_observer obs;
    int32_t year, month, day;
    int32_t day_type, has_rise, has_set, has_dawn, has_dusk;
    int64_t transit, rise, set, dawn, dusk;
    double duration_hours;
} reference_case;

/* Generated with `spp_capi_test --print`; see the header comment. */
static const reference_case kCases[] = {
#include "spp_capi_test_cases.inc"
};

static const int kCaseCount = (int)(sizeof(kCases) / sizeof(kCases[0]));

static int failures = 0;

static void fail(const char* what, const char* detail) {
    ++failures;
    printf("FAIL  %s: %s\n", what, detail);
}

static void expect_status(const char* what, int32_t got, int32_t want) {
    char buf[96];
    if (got != want) {
        snprintf(buf, sizeof buf, "status %d, expected %d", (int)got, (int)want);
        fail(what, buf);
    }
}

static spp_solar_day solar_day(const spp_observer* obs, int32_t y, int32_t m, int32_t d, int32_t* status) {
    spp_solar_day r;
    memset(&r, 0, sizeof r);
    r.struct_size = (int32_t)sizeof r;
    *status = spp_solar_day_utc(obs, y, m, d, &r);
    return r;
}

static void print_table(void) {
    int i;
    for (i = 0; i < kCaseCount; ++i) {
        const reference_case* c = &kCases[i];
        int32_t status;
        const spp_solar_day r = solar_day(&c->obs, c->year, c->month, c->day, &status);
        printf("    {\"%s\", {%.4f, %.4f, %.1f}, %d, %d, %d, %d, %d, %d, %d, %d,\n", c->name, c->obs.lat_deg,
               c->obs.lon_deg, c->obs.elev_m, (int)c->year, (int)c->month, (int)c->day, (int)r.day_type,
               (int)r.has_rise, (int)r.has_set, (int)r.has_civil_dawn, (int)r.has_civil_dusk);
        printf("     %lldLL, %lldLL, %lldLL, %lldLL, %lldLL, %.9f},\n", (long long)r.transit_utc_ms,
               (long long)r.rise_utc_ms, (long long)r.set_utc_ms, (long long)r.civil_dawn_utc_ms,
               (long long)r.civil_dusk_utc_ms, r.duration_hours);
        if (status != SPP_OK) printf("    /* status %d */\n", (int)status);
    }
}

static void check_event(const char* name, const char* label, int32_t has, int64_t got, int32_t has_ref, int64_t ref) {
    char buf[160];
    if (has != has_ref) {
        snprintf(buf, sizeof buf, "%s present %d, expected %d", label, (int)has, (int)has_ref);
        fail(name, buf);
    } else if (has && llabs((long long)(got - ref)) > 1) {
        snprintf(buf, sizeof buf, "%s %lld ms, expected %lld ms", label, (long long)got, (long long)ref);
        fail(name, buf);
    }
}

static double altitude(const spp_observer* obs, int64_t ms) {
    double a = NAN;
    if (spp_upper_limb_altitude(obs, ms, &a) != SPP_OK) fail("upper limb altitude", "status not SPP_OK");
    return a;
}

static void check_definitions(const reference_case* c, const spp_solar_day* r) {
    char buf[160];
    if (r->has_rise) {
        const double a = altitude(&c->obs, r->rise_utc_ms);
        if (fabs(a) > 1e-5) {
            snprintf(buf, sizeof buf, "upper limb at sunrise %.3g deg", a);
            fail(c->name, buf);
        }
        if (!(altitude(&c->obs, r->rise_utc_ms - 60000) < altitude(&c->obs, r->rise_utc_ms + 60000)))
            fail(c->name, "upper limb not rising at sunrise");
    }
    if (r->has_set) {
        const double a = altitude(&c->obs, r->set_utc_ms);
        if (fabs(a) > 1e-5) {
            snprintf(buf, sizeof buf, "upper limb at sunset %.3g deg", a);
            fail(c->name, buf);
        }
        if (!(altitude(&c->obs, r->set_utc_ms - 60000) > altitude(&c->obs, r->set_utc_ms + 60000)))
            fail(c->name, "upper limb not setting at sunset");
    }
    if (r->has_civil_dawn && r->has_rise && !(r->civil_dawn_utc_ms < r->rise_utc_ms)) fail(c->name, "dawn after sunrise");
    if (r->has_rise && r->rise_utc_ms < r->transit_utc_ms && r->has_set && r->set_utc_ms > r->transit_utc_ms &&
        !(r->rise_utc_ms < r->set_utc_ms))
        fail(c->name, "sunrise after sunset");
    if (r->has_civil_dusk && r->has_set && !(r->set_utc_ms < r->civil_dusk_utc_ms)) fail(c->name, "dusk before sunset");
    if (r->day_type == SPP_DAY_POLAR_DAY && (r->has_rise || r->has_set || r->duration_hours != 24.0))
        fail(c->name, "polar day with an event or not 24 h");
    if (r->day_type == SPP_DAY_POLAR_NIGHT && (r->has_rise || r->has_set || r->duration_hours != 0.0))
        fail(c->name, "polar night with an event or not 0 h");
    if (r->day_type == SPP_DAY_NORMAL && !r->has_rise && !r->has_set) fail(c->name, "normal day without any event");
}

static void check_reference(void) {
    int i;
    char buf[160];
    for (i = 0; i < kCaseCount; ++i) {
        const reference_case* c = &kCases[i];
        int32_t status;
        const spp_solar_day r = solar_day(&c->obs, c->year, c->month, c->day, &status);
        if (status != SPP_OK) {
            expect_status(c->name, status, SPP_OK);
            continue;
        }
        if (r.day_type != c->day_type) {
            snprintf(buf, sizeof buf, "day type %d, expected %d", (int)r.day_type, (int)c->day_type);
            fail(c->name, buf);
        }
        check_event(c->name, "transit", 1, r.transit_utc_ms, 1, c->transit);
        check_event(c->name, "sunrise", r.has_rise, r.rise_utc_ms, c->has_rise, c->rise);
        check_event(c->name, "sunset", r.has_set, r.set_utc_ms, c->has_set, c->set);
        check_event(c->name, "civil dawn", r.has_civil_dawn, r.civil_dawn_utc_ms, c->has_dawn, c->dawn);
        check_event(c->name, "civil dusk", r.has_civil_dusk, r.civil_dusk_utc_ms, c->has_dusk, c->dusk);
        if (fabs(r.duration_hours - c->duration_hours) > 1e-6) {
            snprintf(buf, sizeof buf, "duration %.9f h, expected %.9f h", r.duration_hours, c->duration_hours);
            fail(c->name, buf);
        }
        check_definitions(c, &r);
    }
}

static void check_errors(void) {
    const spp_observer ok = {48.8566, 2.3522, 35.0};
    spp_observer bad;
    spp_solar_day r;
    double a, buf[4];
    const int64_t jun2026 = 1782043200000LL; /* 2026-06-21T12:00:00Z */

    memset(&r, 0, sizeof r);
    r.struct_size = (int32_t)sizeof r;
    expect_status("null observer", spp_solar_day_utc(NULL, 2026, 6, 21, &r), SPP_E_ARGUMENT);
    expect_status("null result", spp_solar_day_utc(&ok, 2026, 6, 21, NULL), SPP_E_ARGUMENT);
    r.struct_size = 64;
    expect_status("wrong struct_size", spp_solar_day_utc(&ok, 2026, 6, 21, &r), SPP_E_ARGUMENT);
    r.struct_size = (int32_t)sizeof r;
    expect_status("2026-02-30", spp_solar_day_utc(&ok, 2026, 2, 30, &r), SPP_E_ARGUMENT);
    expect_status("month 13", spp_solar_day_utc(&ok, 2026, 13, 1, &r), SPP_E_ARGUMENT);
    expect_status("1582-10-10 (reform gap)", spp_solar_day_utc(&ok, 1582, 10, 10, &r), SPP_E_ARGUMENT);
    expect_status("year 3001", spp_solar_day_utc(&ok, 3001, 1, 1, &r), SPP_E_RANGE);
    expect_status("year -1001", spp_solar_day_utc(&ok, -1001, 1, 1, &r), SPP_E_RANGE);

    bad = ok;
    bad.lat_deg = 91.0;
    expect_status("latitude 91", spp_solar_day_utc(&bad, 2026, 6, 21, &r), SPP_E_ARGUMENT);
    bad = ok;
    bad.lon_deg = NAN;
    expect_status("longitude NaN", spp_solar_day_utc(&bad, 2026, 6, 21, &r), SPP_E_ARGUMENT);
    bad = ok;
    bad.elev_m = INFINITY;
    expect_status("elevation inf", spp_upper_limb_altitude(&bad, jun2026, &a), SPP_E_ARGUMENT);
    bad = ok;
    bad.elev_m = 20000.0;
    expect_status("elevation 20 km", spp_upper_limb_altitude(&bad, jun2026, &a), SPP_E_ARGUMENT);

    expect_status("null altitude", spp_upper_limb_altitude(&ok, jun2026, NULL), SPP_E_ARGUMENT);
    expect_status("instant after 3000", spp_upper_limb_altitude(&ok, 40000000000000LL, &a), SPP_E_RANGE);
    expect_status("instant before -1000", spp_upper_limb_altitude(&ok, -100000000000000LL, &a), SPP_E_RANGE);
    expect_status("count -1", spp_upper_limb_altitudes(&ok, jun2026, 3600000, -1, buf), SPP_E_ARGUMENT);
    expect_status("count too large", spp_upper_limb_altitudes(&ok, jun2026, 3600000, 1000001, buf), SPP_E_ARGUMENT);
    expect_status("count 0", spp_upper_limb_altitudes(&ok, jun2026, 3600000, 0, NULL), SPP_OK);
    expect_status("null buffer", spp_upper_limb_altitudes(&ok, jun2026, 3600000, 4, NULL), SPP_E_ARGUMENT);
    expect_status("step overflow", spp_upper_limb_altitudes(&ok, jun2026, 9000000000000000000LL, 4, buf), SPP_E_RANGE);
    expect_status("last instant after 3000",
                  spp_upper_limb_altitudes(&ok, jun2026, 400LL * 365 * 86400000, 4, buf), SPP_E_RANGE);
}

static void check_batch(void) {
    const spp_observer obs = {69.6492, 18.9553, 10.0};
    const int64_t start = 1789862400000LL; /* 2026-09-20T00:00:00Z */
    double batch[48];
    int i;
    expect_status("batch", spp_upper_limb_altitudes(&obs, start, 3600000, 48, batch), SPP_OK);
    for (i = 0; i < 48; ++i) {
        const double one = altitude(&obs, start + (int64_t)i * 3600000);
        if (memcmp(&one, &batch[i], sizeof one) != 0) {
            fail("batch", "differs from single calls");
            break;
        }
    }
}

int main(int argc, char** argv) {
    if (argc > 1 && strcmp(argv[1], "--print") == 0) {
        print_table();
        return 0;
    }
    printf("spp_capi %s, interface %d\n", spp_c_library_version(), (int)spp_c_abi_version());
    if (spp_c_abi_version() != SPP_C_ABI_VERSION) fail("interface version", "header and library differ");
    if (spp_c_solar_day_size() != (int32_t)sizeof(spp_solar_day)) fail("spp_solar_day", "size differs from library");
    check_errors();
    check_reference();
    check_batch();
    printf("%d reference cases, %d failures\n", kCaseCount, failures);
    printf("%s\n", failures ? "FAIL" : "PASS");
    return failures ? 1 : 0;
}
