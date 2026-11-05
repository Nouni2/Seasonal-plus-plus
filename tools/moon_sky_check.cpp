// Independent lunar observation/phase references and public API contracts.
#include "spp/spp.hpp"
#include "spp/model/elp_mpp02.hpp"
#include "core/frames.hpp"
#include "error_stats.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <format>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

void require(bool condition, const char* why) {
    if (!condition) throw std::runtime_error(why);
}
template<class Exception, class F> void throws(F f) {
    try { f(); }
    catch (const Exception&) { return; }
    throw std::runtime_error("expected exception not raised");
}
double elapsed(const spp::Time& a, const spp::Time& b) {
    const auto u = a.utc(), v = b.utc();
    return (static_cast<double>(u.day-v.day) + u.frac-v.frac)*spp::seconds_per_day;
}
std::vector<std::string> columns(std::string line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    for (std::string s; std::getline(stream, s, ',');) result.push_back(s);
    return result;
}

void horizons(const char* path, bool details) {
    std::ifstream file(path);
    require(file.good(), "cannot open Horizons fixtures");
    double sky_error = 0.0, eq_error = 0.0, range_error = 0.0, diameter_error = 0.0;
    double illumination_error = 0.0, elongation_error = 0.0, phase_error = 0.0;
    int sky_count = 0, phase_count = 0;
    spp::check::ErrorStats sky, eqs, ranges, diameters, illumination, elongations, phase_angles, bright_limb;
    for (std::string line; std::getline(file, line);) {
        if (line.empty() || line.front() == '#') continue;
        const auto c = columns(line);
        require(c.size() == 16, "invalid Horizons fixture row");
        std::array<double,15> v{};
        for (std::size_t i=0; i<v.size(); ++i) v[i] = std::stod(c[i+1]);
        const auto day = static_cast<std::int64_t>(v[0]);
        const auto t = spp::Time::from_jd(day, v[0]-static_cast<double>(day), v[4]);
        const auto context = std::format("{} {} DUT1={:+.5f}s", c[0], t.iso_utc(), v[4]);
        auto light=[&](double fraction,double angle) {
            illumination_error=std::max(illumination_error,std::abs(fraction-v[11]));
            phase_error=std::max(phase_error,std::abs(angle-v[13])*3600);
            illumination.add((fraction-v[11])*100,context);
            phase_angles.add((angle-v[13])*3600,context);
        };
        spp::Equatorial eq{};
        if (c[0] == "geocenter") {
            eq = spp::moon::apparent(t);
            const auto p = spp::moon::phase(t);
            light(p.illuminated_fraction,p.phase_angle_deg);
            elongation_error = std::max(elongation_error, std::abs(p.elongation_deg-v[12])*3600);
            elongations.add((p.elongation_deg-v[12])*3600, context);
            if (details)
                std::printf("%s\n  illumination: Spp %.7f%%  JPL %.7f%%  delta %+.7f percentage points\n",
                    context.c_str(), p.illuminated_fraction*100, v[11]*100,
                    (p.illuminated_fraction-v[11])*100);
            ++phase_count;
        } else {
            const auto o = spp::moon::observe(t, {v[1],v[2],v[3]}, {.refraction=false});
            const auto h = spp::moon::horizontal(t, {v[1],v[2],v[3]}, {.refraction=false});
            light(o.disk.illuminated_fraction,o.disk.phase_angle_deg);
            require(o.disk.bright_limb_defined,"fixture limb direction unexpectedly undefined");
            bright_limb.add(spp::wrap180(o.disk.bright_limb_deg-v[14])*3600,context);
            require(h.alt_deg == o.horizontal.alt_deg && h.az_deg == o.horizontal.az_deg,
                    "horizontal-only path differs at a Horizons reference epoch");
            eq = o.equatorial;
            const double daz = spp::wrap180(o.horizontal.az_deg-v[7])*std::cos(v[8]*spp::deg2rad);
            sky_error = std::max(sky_error, std::hypot(daz, o.horizontal.alt_deg-v[8])*3600);
            range_error = std::max(range_error, std::abs(o.distance_km-v[9]));
            diameter_error = std::max(diameter_error, std::abs(o.semidiameter_deg*7200-v[10]));
            sky.add(std::hypot(daz,o.horizontal.alt_deg-v[8])*3600, context);
            ranges.add((o.distance_km-v[9])*1000, context);
            diameters.add(o.semidiameter_deg*7200-v[10], context);
            if (details)
                std::printf("%s\n  altitude: Spp %+.9f  JPL %+.9f deg  delta %+.6f arcsec\n"
                    "  azimuth:  Spp %.9f  JPL %.9f deg  delta %+.6f arcsec\n"
                    "  distance: Spp %.6f  JPL %.6f km  delta %+.6f m\n",
                    context.c_str(), h.alt_deg,v[8],(h.alt_deg-v[8])*3600,
                    h.az_deg,v[7],spp::wrap180(h.az_deg-v[7])*3600,
                    o.distance_km,v[9],(o.distance_km-v[9])*1000);
            ++sky_count;
        }
        const double dra = spp::wrap180(eq.ra_deg-v[5])*std::cos(v[6]*spp::deg2rad);
        eq_error = std::max(eq_error, std::hypot(dra,eq.dec_deg-v[6])*3600);
        eqs.add(std::hypot(dra,eq.dec_deg-v[6])*3600, context);
    }
    std::printf("JPL DE441: %d sky / %d geocentric phase epochs\n"
        "  max sky %.6f arcsec, equatorial %.6f arcsec, range %.6f km, diameter %.6f arcsec\n"
        "  max illumination %.9f fraction, elongation %.6f arcsec, phase angle %.6f arcsec\n",
        sky_count, phase_count, sky_error, eq_error, range_error, diameter_error,
        illumination_error, elongation_error, phase_error);
    std::puts("Independent JPL DE441 residuals (airless, matched DUT1):");
    bool pass = sky.report("sky direction", "arcsec", 1.5);
    pass = eqs.report("equatorial direction", "arcsec", 0.5) && pass;
    pass = ranges.report("observer range", "m", 100.0) && pass;
    pass = diameters.report("angular diameter", "arcsec", 0.1) && pass;
    pass = illumination.report("illuminated fraction", "percentage points", 0.00002) && pass;
    pass = elongations.report("solar elongation", "arcsec", 0.5) && pass;
    pass = phase_angles.report("physical phase (phi)", "arcsec", 0.5) && pass;
    pass = bright_limb.report("bright limb (0.01 deg ref)","arcsec",60.0) && pass;
    std::puts("Sky direction uses delta-azimuth*cos(altitude) and delta-altitude.\n"
              "JPL includes EOP corrections absent from Spp; these are sampled residuals, not global bounds.");
    require(sky_count == 15 && phase_count == 5, "Horizons fixture cases missing");
    // Horizons includes measured/predicted pole motion and celestial-pole
    // corrections; Spp uses its analytical Earth orientation, without EOP.
    require(sky_error < 1.5 && eq_error < 0.5, "lunar apparent coordinates disagree with JPL");
    require(range_error < 0.1 && diameter_error < 0.1, "lunar distance/angular size disagree with JPL");
    // Illu% and phi use the physical illumination convention; S-T-O is an
    // approximate angle with an additional outgoing aberration correction.
    require(illumination_error < 2.0e-7 && elongation_error < 0.5 && phase_error < 0.5,
            "lunar phase geometry disagrees with JPL");
    require(pass, "JPL residual limits exceeded");
}

void usno(const char* path, bool details) {
    std::ifstream file(path);
    require(file.good(), "cannot open USNO fixtures");
    std::vector<spp::moon::PhaseEvent> expected;
    std::vector<int> years;
    for (std::string line; std::getline(file,line);) {
        if (line.empty() || line.front() == '#') continue;
        const auto c=columns(line);
        require(c.size()==6, "invalid USNO fixture row");
        const int y=std::stoi(c[0]);
        if (years.empty() || years.back()!=y) years.push_back(y);
        expected.push_back({static_cast<spp::moon::PrincipalPhase>(std::stoi(c[5])),
            spp::Time::from_utc(y,std::stoi(c[1]),std::stoi(c[2]),std::stoi(c[3]),std::stoi(c[4]))});
    }
    double error=0.0, residual=0.0;
    std::size_t index=0;
    spp::check::ErrorStats times, roots;
    for (int y:years) {
        const auto start=spp::Time::from_utc(y,1,1), end=spp::Time::from_utc(y+1,1,1);
        const auto found=spp::moon::phases(start,end);
        double year_error=0.0;
        for (const auto& e:found) {
            require(index<expected.size(), "extra principal phase found");
            require(e.phase==expected[index].phase, "principal phase ordering wrong");
            const double dt=elapsed(e.time,expected[index].time);
            const double angle=spp::wrap180(spp::moon::phase(e.time).cycle_deg
                -90.0*static_cast<int>(e.phase));
            year_error=std::max(year_error,std::abs(dt));
            residual=std::max(residual,std::abs(angle));
            const auto reference=expected[index].time.to_civil();
            const auto context=std::format("{} Spp {} / USNO {} {:02}:{:02} UT",
                spp::moon::name(e.phase),e.time.iso_utc(),spp::iso_date(reference.date()),reference.hour,reference.minute);
            times.add(dt,context);
            roots.add(angle*3600,context);
            if (details) std::printf("%s  delta %+.6f s  longitude residual %+.9f arcsec\n",context.c_str(),dt,angle*3600);
            ++index;
        }
        std::printf("USNO %d: %zu phases, max difference %.3f s (reference at minute precision)\n",y,found.size(),year_error);
        error=std::max(error,year_error);
    }
    require(years.size()==4 && index==expected.size() && index>190, "USNO phase cases missing");
    std::puts("Independent USNO event times and separate numerical root residuals:");
    const bool times_ok=times.report("principal phase time","s",60.0);
    const bool roots_ok=roots.report("longitude root residual","arcsec",3.0e-6*3600);
    std::puts("USNO times have minute precision and use UT; root precision is not physical timing accuracy.");
    require(error<60.0, "principal phase times disagree with USNO");
    require(residual<3.0e-6, "principal phase longitude root not converged");
    require(times_ok && roots_ok, "USNO time/root residual limits exceeded");
}

struct Alternative {
    int calls=0;
    spp::moon::State operator()(spp::Instant<spp::TDB>) {
        ++calls;
        return {{400000,0,0},{0,0,0}};
    }
};

struct PositionOnly {
    const spp::model::elp_mpp02::Ephemeris model;
    int calls = 0;
    spp::moon::Position position(spp::Instant<spp::TDB> t) {
        ++calls;
        return model.position(t);
    }
};
static_assert(spp::moon::PositionSource<PositionOnly>);
static_assert(!spp::moon::Ephemeris<PositionOnly>);

struct CountedModel {
    const spp::model::elp_mpp02::Ephemeris model;
    int position_calls = 0, state_calls = 0, reception_calls = 0;
    spp::Instant<spp::TDB> reception;

    explicit CountedModel(spp::Instant<spp::TDB> t) : reception(t) {}
    spp::moon::Position position(spp::Instant<spp::TDB> t) {
        ++position_calls;
        if (t.day == reception.day && t.frac == reception.frac) ++reception_calls;
        return model.position(t);
    }
    spp::moon::State operator()(spp::Instant<spp::TDB> t) {
        ++state_calls;
        return model(t);
    }
};

struct StateOnly {
    const spp::model::elp_mpp02::Ephemeris model;
    spp::moon::State operator()(spp::Instant<spp::TDB> t) const { return model(t); }
};

struct InvalidPosition {
    spp::moon::Position position(spp::Instant<spp::TDB>) const {
        return {{std::numeric_limits<double>::quiet_NaN(),0,0}};
    }
};

void observation_paths() {
    const auto t = spp::Time::from_utc(2026,11,3,12,0,0.0,0.0);
    const spp::Observer obs{48.8566,2.3522,35.0};
    CountedModel coordinates(t.tdb()), complete(t.tdb());
    const auto horizontal = spp::moon::horizontal(t,obs,coordinates,{.refraction=false});
    const auto observation = spp::moon::observe(t,obs,complete,{.refraction=false});
    require(coordinates.state_calls == 0 && complete.state_calls == 0,
            "observer calculations evaluated unused model velocities");
    require(coordinates.position_calls <= 4 && coordinates.position_calls < complete.position_calls,
            "horizontal-only calculation did extra geocentric phase work");
    require(complete.position_calls <= 7 && complete.reception_calls == 1,
            "full observation did not share its reception-time position");
    require(horizontal.alt_deg == observation.horizontal.alt_deg
            && horizontal.az_deg == observation.horizontal.az_deg,
            "horizontal and full observation use different observer calculations");

    PositionOnly positions;
    const auto only = spp::moon::observe(t,obs,positions,{.refraction=false});
    require(positions.calls > 0 && only.distance_km == observation.distance_km
            && only.phase.illuminated_fraction == observation.phase.illuminated_fraction,
            "position-only provider not supported by the full observation");
    const StateOnly states;
    const auto fallback = spp::moon::observe(t,obs,states,{.refraction=false});
    require(std::abs(fallback.horizontal.alt_deg-only.horizontal.alt_deg) < 1e-9
            && std::abs(spp::wrap180(fallback.horizontal.az_deg-only.horizontal.az_deg)) < 1e-9
            && std::abs(fallback.distance_km-only.distance_km) < 1e-6,
            "position-only and full-state observation paths disagree");
    require(spp::moon::apparent(t,positions).ra_deg == spp::moon::apparent(t,states).ra_deg
            && spp::moon::phase(t,positions).cycle_deg == spp::moon::phase(t,states).cycle_deg,
            "position-only provider not supported by apparent/phase APIs");
    const InvalidPosition invalid;
    throws<std::invalid_argument>([&]{spp::moon::horizontal(t,obs,invalid);});
    throws<std::invalid_argument>([&]{spp::moon::phase(t,invalid);});
    const auto begin = spp::Time::from_utc(2026,11,1);
    const auto expected = spp::moon::phases(begin,begin.shifted(10.0),states);
    const auto found = spp::moon::phases(begin,begin.shifted(10.0),positions);
    require(!found.empty() && found.size() == expected.size()
            && found.front().phase == expected.front().phase
            && std::abs(elapsed(found.front().time,expected.front().time)) < 1e-6,
            "position-only provider not supported by phase search");
}

void contracts() {
    const auto t=spp::Time::from_utc(2026,11,3,12,0,0.0,0.0);
    const spp::Observer obs{48.8566,2.3522,35.0};
    const auto a=spp::moon::observe(t,obs,{.refraction=false});
    const auto b=spp::moon::observe(t,obs);
    require(a.horizontal.alt_deg==a.geometric_alt_deg && b.horizontal.alt_deg>a.horizontal.alt_deg,
            "lunar refraction not applied separately");
    require(a.distance_km==b.distance_km && a.horizontal.az_deg==b.horizontal.az_deg,
            "refraction changed azimuth/range");
    const auto h=spp::moon::horizontal(t,obs,{.refraction=false});
    require(h.alt_deg==a.horizontal.alt_deg && h.az_deg==a.horizontal.az_deg, "horizontal wrapper differs");
    require(a.phase.illuminated_fraction==b.phase.illuminated_fraction && !a.phase.waxing
            && std::string(a.phase.name())=="Waning Crescent", "phase/refraction contract wrong");
    const spp::model::elp_mpp02::Ephemeris llr{spp::model::elp_mpp02::Fit::LLR};
    require(spp::moon::observe(t,obs,llr).distance_km!=a.distance_km, "observation ignores selected fit");
    Alternative model;
    const auto alternate=spp::moon::observe(t,obs,model);
    require(model.calls>=2 && std::abs(alternate.distance_km-a.distance_km)>10000,
            "observation ignores stateful provider");
    require(spp::moon::phase(t,model).cycle_deg!=a.phase.cycle_deg, "phase ignores provider");
    require(spp::moon::apparent(t,model).ra_deg!=a.equatorial.ra_deg, "apparent ignores provider");
    throws<std::runtime_error>([&]{spp::moon::phases(t,t.shifted(35),model);});
    for (double lat:{-90.0,90.0}) {
        const auto pole=spp::moon::observe(t,{lat,0,0});
        require(std::isfinite(pole.horizontal.alt_deg) && std::isfinite(pole.horizontal.az_deg),
                "observer pole singularity");
    }
    const auto equator=spp::sky::observer_position_km({0,0,0},0);
    require(equator[0]==spp::wgs84_a_m/1000 && equator[1]==0 && equator[2]==0, "observer units wrong");
    const auto high=spp::moon::observe(t,{obs.lat_deg,obs.lon_deg,5000},{.refraction=false});
    require(high.distance_km!=a.distance_km, "observer elevation ignored");
    const auto rotated=spp::moon::observe(spp::Time::from_utc(2026,11,3,12,0,0.0,0.5),obs,{.refraction=false});
    require(std::abs(rotated.horizontal.az_deg-a.horizontal.az_deg)>1e-5, "DUT1 ignored");
    throws<std::invalid_argument>([&]{spp::moon::observe(t,{91,0,0});});
    throws<std::invalid_argument>([&]{spp::moon::observe(t,{0,std::numeric_limits<double>::quiet_NaN(),0});});
    throws<std::invalid_argument>([&]{spp::moon::observe(t,obs,{.pressure_mbar=-1.0});});
    throws<std::invalid_argument>([&]{spp::moon::observe(t,obs,{.temp_c=-273.15});});
    throws<std::invalid_argument>([&]{spp::moon::observe(t,obs,{.temp_c=-273.0});});
    throws<std::invalid_argument>([&]{spp::moon::phases(t,t.shifted(-1));});
    throws<std::invalid_argument>([&]{spp::moon::phases(t,t.shifted(1),0);});
    require(spp::moon::phases(t,t).empty(), "empty phase interval not empty");
    const auto events=spp::moon::phases(t,t.shifted(35));
    const auto first=events.front().time;
    require(spp::moon::phases(first.shifted(-1e-4),first.shifted(1e-4)).size()==1,
            "narrow interval missed phase");
    require(spp::moon::phases(first.shifted(1e-4),first.shifted(1e-3)).empty(), "duplicate phase past root");
    const auto tt=t.tt();
    const auto td=t.tdb();
    require(std::abs((static_cast<double>(td.day-tt.day)+td.frac-tt.frac)*spp::seconds_per_day)<0.002,
            "TDB correction outside expected magnitude");
    const auto split=spp::tdb_from_tt({tt.day-1,tt.frac+1});
    const auto canonical=spp::tdb_from_tt({tt.day,std::fmod(tt.frac+1,1.0)});
    require(split.day==canonical.day && split.frac==canonical.frac, "TDB split normalization wrong");
    throws<std::invalid_argument>([]{spp::tdb_from_tt({2451545,std::numeric_limits<double>::infinity()});});
    throws<std::out_of_range>([&]{spp::moon::phases(t,spp::Time::from_utc(5000,1,1));});
    for (double century:{-30.0,-1.0,0.0,0.26,10.0}) {
        const auto q=spp::Time::from_jd(2451545,century*36525,0);
        const auto n=spp::earth::nutation(q.tt());
        const spp::detail::Vector v{300000,-200000,10000};
        const auto w=spp::detail::ecliptic_j2000_to_equatorial(v,q.tt(),n);
        const auto r=spp::detail::equatorial_to_ecliptic_j2000(w,q.tt(),n);
        for (int i=0;i<3;++i) require(std::abs(v[i]-r[i])<1e-8, "frame inverse inconsistent");
    }
    std::puts("Lunar observer/model/time/interval contracts: PASS");
}

}  // namespace

int main(int argc,char** argv) {
    try {
        const bool details=argc==4 && std::string_view(argv[3])=="--details";
        if (argc!=3 && !details)
            throw std::runtime_error("usage: spp_moon_sky_check horizons.csv usno_phases.csv [--details]");
        bool failed=false;
        auto check=[&](const char* name,auto body) {
            try { body(); }
            catch(const std::exception& e) {
                std::fprintf(stderr,"FAIL [%s]: %s\n",name,e.what());
                failed=true;
            }
        };
        check("JPL",[&]{horizons(argv[1],details);});
        check("USNO",[&]{usno(argv[2],details);});
        check("contracts",[]{contracts();});
        check("observation paths",[]{observation_paths();});
        if (failed) return 1;
        std::puts("PASS");
        return 0;
    } catch(const std::exception& e) {
        std::fprintf(stderr,"FAIL: %s\n",e.what());
        return 1;
    }
}
