// spp: command-line front end for the Spp library.
//
//   spp sun  YYYY-MM-DD HH:MM[:SS] [LAT LON [ELEV]]   apparent Sun (UTC input)
//   spp day  YYYY-MM-DD LAT LON [ELEV]                transit, sunrise, sunset
//   spp year YYYY LAT LON [ELEV]                      one CSV row per day
//   spp moon YYYY-MM-DD HH:MM[:SS] [--fit=de405|llr]  geometric Moon (TDB input)
//   spp moon sky YYYY-MM-DD HH:MM[:SS] LAT LON [ELEV]  lunar sky (UTC input)
//   spp moon phase YYYY-MM-DD HH:MM[:SS]              phases (UTC input)
//
// Options: --dut1=S (default: estimated)  --pressure=MBAR (default: standard at
// the elevation)  --temp=C  --no-refraction
// Longitude is positive East; elevation in metres.

#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <print>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>
#include <variant>

#include "spp/spp.hpp"
#if SPP_CLI_MOON
#include "spp/model/elp_mpp02.hpp"
#include "spp/model/iau_moon.hpp"
#include "spp/model/lunar_numerical.hpp"
#include "spp/model/meeus.hpp"
#endif

namespace {

struct Options {
    std::optional<double> dut1;   // estimated when not given
    spp::Atmosphere atm;
};

[[noreturn]] void usage(int status=2) {
    std::println(status?stderr:stdout,
                 "usage:\n"
                 "  spp sun  YYYY-MM-DD HH:MM[:SS] [LAT LON [ELEV]]\n"
                 "  spp day  YYYY-MM-DD LAT LON [ELEV]\n"
                 "  spp year YYYY LAT LON [ELEV]\n"
#if SPP_CLI_MOON
                 "  spp moon YYYY-MM-DD HH:MM[:SS] [--fit=de405|llr]  (TDB input)\n"
                 "  spp moon sky YYYY-MM-DD HH:MM[:SS] LAT LON [ELEV] (UTC input)\n"
                 "  spp moon phase YYYY-MM-DD HH:MM[:SS]              (UTC input)\n"
                 "  spp moon surface YYYY-MM-DD HH:MM[:SS] [LAT LON [ELEV]] (UTC input)\n"
                 "  spp moon day YYYY-MM-DD LAT LON [ELEV]             (UTC day)\n"
                 "  spp moon rotation YYYY-MM-DD HH:MM[:SS]           (TDB input, IAU approximation)\n"
                 "  orbital moon commands accept --fit=de405|llr\n"
                 "  --model=integrated --kernel=PATH --seed=YYYY-MM-DDTHH:MM[:SS] (TDB)\n"
                 "    optional --span-days=45: propagate this far each side of the seed\n"
#endif
                 "options: --dut1=S --pressure=MBAR --temp=C --no-refraction\n"
                 "  DUT1 is estimated and pressure follows the elevation unless given");
    std::exit(status);
}

double num(std::string_view s) {
    double v{};
    auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc{} || p != s.data() + s.size() || !std::isfinite(v)) {
        std::println(stderr, "spp: not a finite number: '{}'", s);
        std::exit(2);
    }
    return v;
}

int integer(std::string_view s) {
    int v{};
    auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc{} || p != s.data() + s.size()) {
        std::println(stderr, "spp: not an integer: '{}'", s);
        std::exit(2);
    }
    return v;
}

// Splits "a<sep>b<sep>c" into integers/numbers; leading '-' allowed on the first field.
std::vector<std::string_view> fields(std::string_view s, char sep) {
    std::vector<std::string_view> out;
    std::size_t start = s.starts_with('-') ? 1 : 0;
    std::size_t from = 0;
    for (std::size_t i = start; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == sep) {
            out.push_back(s.substr(from, i - from));
            from = i + 1;
        }
    }
    return out;
}

spp::CivilDate parse_date(std::string_view s) {
    auto f = fields(s, '-');
    if (f.size() != 3) usage();
    return {integer(f[0]), integer(f[1]), integer(f[2])};
}

spp::Observer parse_observer(const std::vector<std::string_view>& a, std::size_t i) {
    return {num(a[i]), num(a[i + 1]), a.size() > i + 2 ? num(a[i + 2]) : 0.0};
}

// Full ISO timestamp: events can fall on the neighbouring UTC date.
std::string when(const std::optional<spp::Time>& t) {
    if (!t) return "-";
    return t->iso_utc();
}

const char* type_name(spp::events::DayType t) {
    switch (t) {
        case spp::events::DayType::Normal: return "NORMAL";
        case spp::events::DayType::PolarDay: return "POLAR_DAY";
        case spp::events::DayType::PolarNight: return "POLAR_NIGHT";
    }
    return "?";
}

int cmd_sun(const std::vector<std::string_view>& a, const Options& o) {
    if (a.size() != 2 && a.size() != 4 && a.size() != 5) usage();
    const auto d = parse_date(a[0]);
    const auto hms = fields(a[1], ':');
    if (hms.size() < 2 || hms.size() > 3) usage();
    const auto t = spp::Time::from_utc(d.year, d.month, d.day, integer(hms[0]), integer(hms[1]),
                                       hms.size() == 3 ? num(hms[2]) : 0.0, o.dut1);

    const auto tt = t.tt();
    const auto n = spp::earth::nutation(tt);
    const auto ecl = spp::sun::apparent_ecliptic(tt, n);
    const auto eq = spp::sun::apparent(tt, n);

    std::println("time        {}  (ΔT {:.2f} s, DUT1 {:+.3f} s{})", t.iso_utc(), t.delta_t(), t.dut1(),
                 t.dut1_estimated() ? " estimated" : "");
    std::println("ecliptic    λ {:.8f}°  β {:.8f}°  R {:.9f} AU", ecl.lon_deg, ecl.lat_deg, ecl.dist_au);
    std::println("equatorial  α {:.8f}° ({:.8f} h)  δ {:.8f}°", eq.ra_deg, eq.ra_hours(), eq.dec_deg);
    std::println("nutation    Δψ {:.4f}\"  Δε {:.4f}\"  ε {:.8f}°", n.dpsi_deg * 3600.0,
                 n.deps_deg * 3600.0, n.true_obliquity_deg);
    std::println("sidereal    GMST {:.8f}°  GAST {:.8f}°", spp::earth::gmst_deg(t.ut1()),
                 spp::earth::gast_deg(t.ut1(), n));

    if (a.size() >= 4) {
        const auto obs = parse_observer(a, 2);
        const auto hor = spp::sun::horizontal(t, obs, o.atm);
        const auto geo = spp::sun::horizontal(t, obs, {.refraction = false});
        std::println("horizontal  alt {:.6f}° (geometric {:.6f}°)  az {:.6f}°", hor.alt_deg, geo.alt_deg, hor.az_deg);
    }
    return 0;
}

int cmd_day(const std::vector<std::string_view>& a, const Options& o) {
    if (a.size() != 3 && a.size() != 4) usage();
    const auto d = spp::events::day(parse_date(a[0]), parse_observer(a, 1), o.atm, o.dut1);
    std::println("type      {}", type_name(d.type));
    std::println("transit   {}", d.transit.iso_utc());
    std::println("sunrise   {}", d.rise ? d.rise->iso_utc() : "-");
    std::println("sunset    {}", d.set ? d.set->iso_utc() : "-");
    std::println("duration  {:.4f} h", d.duration_hours);
    return 0;
}

#if SPP_CLI_MOON
void print_phase(const spp::moon::Phase& p) {
    std::println("phase       {} ({})", p.name(), p.waxing ? "waxing" : "waning");
    std::println("illuminated {:.4f}%", p.illuminated_fraction * 100.0);
    std::println("cycle       {:.6f} deg (0 new, 90 first quarter, 180 full, 270 last quarter)", p.cycle_deg);
    std::println("elongation  {:.6f} deg from the Sun", p.elongation_deg);
    std::println("phase angle {:.6f} deg (0 fully illuminated)", p.phase_angle_deg);
}
void print_surface(const spp::moon::Surface& s) {
    std::println("orientation IAU 2009 approximation, at reflection (Mean Earth/Polar Axis)");
    std::println("subobserver lon {:+.6f} deg  lat {:+.6f} deg (east-positive)",s.subobserver_lon_deg,s.subobserver_lat_deg);
    std::println("subsolar    lon {:+.6f} deg  lat {:+.6f} deg (east-positive)",s.subsolar_lon_deg,s.subsolar_lat_deg);
    if(s.pole_defined)std::println("north pole  {:.6f} deg (east from celestial north)",s.pole_position_angle_deg);
    else std::println("north pole  sky bearing undefined");
}

spp::CivilDate moon_date(std::string_view value) {
    const auto date=parse_date(value);
    if(date.year< -2000 || date.year>4000)
        throw std::out_of_range("moon: date outside the ephemeris interval");
    if(!spp::is_valid(date))throw std::invalid_argument("moon: invalid calendar date");
    return date;
}

// Runtime selection belongs at the CLI boundary. The library's observation
// algorithms remain independently parameterized by an ephemeris provider.
class MoonProvider {
    using Direct=spp::model::elp_mpp02::Ephemeris;
    using Cached=spp::model::elp_mpp02::CachedEphemeris;
    using Integrated=spp::model::lunar_numerical::Ephemeris;
    std::variant<Direct,Cached,Integrated> model_;
    spp::model::elp_mpp02::Fit fit_;
    std::optional<spp::Instant<spp::TDB>> seed_;double span_;
public:
    MoonProvider(spp::model::elp_mpp02::Fit fit,bool direct,bool integrated,
        std::string_view kernel,std::optional<spp::Instant<spp::TDB>> seed,double span):model_(Direct(fit)),fit_(fit),seed_(seed),span_(span) {
        if(integrated) {
            if(kernel.empty() || !seed)throw std::invalid_argument("integrated Moon requires --kernel and explicit --seed (TDB)");
            model_.emplace<Integrated>(spp::model::spk::Kernel(std::string(kernel)),*seed,
                spp::Instant<spp::TDB>{seed->day,seed->frac-span},spp::Instant<spp::TDB>{seed->day,seed->frac+span});
        } else if(!direct)model_.emplace<Cached>(fit);
    }
    spp::moon::State operator()(spp::Instant<spp::TDB> t) {
        return std::visit([&](auto& model){return model(t);},model_);
    }
    spp::moon::Position position(spp::Instant<spp::TDB> t) {
        return std::visit([&](auto& model){return model.position(t);},model_);
    }
    void describe()const {
        if(const auto* integrated=std::get_if<Integrated>(&model_)) {
            const auto& s=integrated->stats();
            std::println("model       numerical Moon, DE440 constants; kernel seed and prescribed planetary trajectories");
            std::println("forces      major bodies + 1PN + Earth degree 5 + Moon degree 6 + tidal deformation and figure reactions");
            std::println("solar       J2, radiation pressure and spin frame dragging");
            std::println("seed        JD {} + {:.12f} TDB; propagated +/-{:g} days",seed_->day,seed_->frac,span_);
            std::println("propagation {} accepted / {} rejected steps, {} force evaluations",s.accepted_steps,s.rejected_steps,s.force_evaluations);
        } else std::println("model       ELP/MPP02, {} fit",fit_==spp::model::elp_mpp02::Fit::LLR?"LLR":"DE405/DE406");
    }
};
spp::Instant<spp::TDB> moon_seed(std::string_view text) {
    const auto fields_t=fields(text,'T');
    if(fields_t.size()!=2)throw std::invalid_argument("--seed requires YYYY-MM-DDTHH:MM[:SS] in TDB");
    const auto date=moon_date(fields_t[0]);const auto hms=fields(fields_t[1],':');
    if(hms.size()!=2 && hms.size()!=3)throw std::invalid_argument("invalid --seed clock");
    const int hour=integer(hms[0]),minute=integer(hms[1]);const double second=hms.size()==3?num(hms[2]):0;
    if(hour<0 || hour>23 || minute<0 || minute>59 || second<0 || second>=60)throw std::invalid_argument("invalid --seed clock");
    const auto jd=spp::model::meeus::jd_from_civil({date.year,date.month,date.day,hour,minute,second});
    return {jd.day,jd.frac};
}
int cmd_moon(std::vector<std::string_view> a, const Options& o) {
    using spp::model::elp_mpp02::Fit;
    Fit fit = Fit::DE405;
    bool fit_given=false;
    bool integrated=false,orbit_options=false;std::string_view kernel;
    std::optional<spp::Instant<spp::TDB>> seed;double span=45;
    for (auto i = a.begin(); i != a.end();) {
        const auto option=*i;
        if(option.starts_with("--model=")) {
            if(option!="--model=elp" && option!="--model=integrated")throw std::invalid_argument("choose --model=elp or --model=integrated");
            integrated=option=="--model=integrated";orbit_options=true;i=a.erase(i);continue;
        }
        if(option.starts_with("--kernel=")){kernel=option.substr(9);orbit_options=true;i=a.erase(i);continue;}
        if(option.starts_with("--seed=")){seed=moon_seed(option.substr(7));orbit_options=true;i=a.erase(i);continue;}
        if(option.starts_with("--span-days=")){span=num(option.substr(12));orbit_options=true;i=a.erase(i);continue;}
        if (!i->starts_with("--fit=")) { ++i; continue; }
        fit_given=true;
        if (*i == "--fit=llr") fit = Fit::LLR;
        else if (*i == "--fit=de405") fit = Fit::DE405;
        else
            throw std::invalid_argument("moon: choose --fit=de405 or --fit=llr");
        i = a.erase(i);
    }
    if(!std::isfinite(span) || span<=0 || span>365)throw std::invalid_argument("--span-days must be positive and at most 365");
    if(integrated && fit_given)throw std::invalid_argument("--fit selects ELP parameters and cannot be combined with --model=integrated");
    if(!integrated && (!kernel.empty() || seed || span!=45))throw std::invalid_argument("kernel/seed/span options require --model=integrated");
    if (a.empty()) usage();
    if(a[0]=="day") {
        a.erase(a.begin());
        if(a.size()!=3 && a.size()!=4)usage();
        const auto date=moon_date(a[0]);const auto obs=parse_observer(a,1);
        const auto start=spp::Time::from_utc(date.year,date.month,date.day,0,0,0,o.dut1);
        MoonProvider model(fit,false,integrated,kernel,seed,span);
        const auto found=spp::moon::horizon(start,start.shifted(1),obs,model,o.atm);
        std::println("Moon - horizon events for {} (UTC day)",spp::iso_date(date));
        model.describe();
        std::println("at midnight {} the horizon",found.above_at_start?"above":"below");
        std::println("above       {:.6f} hours (upper limb, selected atmosphere)",found.duration_hours);
        bool rise=false,set=false;
        for(const auto& e:found.events) {
            rise=rise || e.kind==spp::moon::HorizonKind::Rise;
            set=set || e.kind==spp::moon::HorizonKind::Set;
            std::println("{:15} {}  center alt {:+.6f} deg  az {:.6f} deg",spp::moon::name(e.kind),
                e.time.iso_utc(),e.center_altitude_deg,e.azimuth_deg);
        }
        if(!rise)std::println("Moonrise        none in this UTC day");
        if(!set)std::println("Moonset         none in this UTC day");
        return 0;
    }
    const bool sky = a[0] == "sky", phase = a[0] == "phase", rotation = a[0] == "rotation", surface = a[0]=="surface";
    if(rotation && fit_given)throw std::invalid_argument("moon rotation: the IAU orientation model is independent of the orbital --fit option");
    if(rotation && orbit_options)throw std::invalid_argument("moon rotation is independent of orbital model/kernel options");
    if (sky || phase || rotation || surface) a.erase(a.begin());
    if(sky ? (a.size()!=4 && a.size()!=5) : surface ? (a.size()!=2 && a.size()!=4 && a.size()!=5) : a.size()!=2)usage();
    const auto date = moon_date(a[0]);
    const auto hms = fields(a[1], ':');
    if (hms.size() != 2 && hms.size() != 3) usage();
    const int hour = integer(hms[0]), minute = integer(hms[1]);
    const double second = hms.size() == 3 ? num(hms[2]) : 0.0;
    if (hour < 0 || hour > 23 || minute < 0 || minute > 59 || second < 0.0 || second >= 60.0)
        throw std::invalid_argument("moon: time must be 00:00:00 through 23:59:59.999...");

    if (sky || phase || surface) {
        MoonProvider model(fit,false,integrated,kernel,seed,span);
        const auto t = spp::Time::from_utc(date.year, date.month, date.day, hour, minute, second, o.dut1);
        std::println("Moon — {}", sky ? "observer's sky" : surface ? "surface orientation" : "geocentric phases");
        std::println("time        {} (UTC)", t.iso_utc());
        model.describe();
        if(surface) {
            if(a.size()>=4)print_surface(spp::moon::surface(t,parse_observer(a,2),model));
            else print_surface(spp::moon::surface(t,model));
        } else if (sky) {
            const auto obs = parse_observer(a, 2);
            const auto result = spp::moon::observe(t, obs, model, o.atm);
            std::println("observer    lat {:+.6f} deg  lon {:+.6f} deg  elevation {:.1f} m",
                         obs.lat_deg, obs.lon_deg, obs.elev_m);
            std::println("DUT1        {:+.3f} s{}", t.dut1(), t.dut1_estimated() ? " (estimated)" : "");
            std::println("altitude    {:+.6f} deg (geometric {:+.6f} deg, center{})",
                         result.horizontal.alt_deg, result.geometric_alt_deg,
                         o.atm.refraction ? ", with refraction" : ", no refraction");
            std::println("azimuth     {:.6f} deg (0 north, 90 east, 180 south, 270 west)", result.horizontal.az_deg);
            std::println("equatorial  RA {:.8f} deg  Dec {:+.8f} deg (topocentric apparent, of date)",
                         result.equatorial.ra_deg, result.equatorial.dec_deg);
            std::println("distance    {:.3f} km from observer", result.distance_km);
            std::println("diameter    {:.4f} arcmin", result.semidiameter_deg * 120.0);
            std::println("disk lit    {:.6f}% (topocentric)", result.disk.illuminated_fraction*100);
            if(result.disk.bright_limb_defined)
                std::println("bright limb {:.6f} deg (east from celestial north)",result.disk.bright_limb_deg);
            if(result.disk.zenith_defined)
                std::println("parallactic {:+.6f} deg (celestial north to local zenith)",result.disk.parallactic_angle_deg);
            print_surface(result.surface);
            print_phase(result.phase);
        } else {
            print_phase(spp::moon::phase(t, model));
            std::println("\nPrincipal phases during the next 35 days (UTC):");
            for (const auto& event : spp::moon::phases(t, t.shifted(35.0), model))
                std::println("  {:14} {}", spp::moon::name(event.phase), event.time.iso_utc());
        }
        return 0;
    }

    // This converts calendar fields to a Julian date only. The input fields
    // explicitly represent TDB; no UTC-to-TDB approximation or relabeling.
    const auto jd = spp::model::meeus::jd_from_civil(
        {date.year, date.month, date.day, hour, minute, second});
    if(rotation) {
        const auto result=spp::model::iau_moon::Rotation{}({jd.day,jd.frac});
        std::println("Moon - IAU 2009 rotation (approximate Mean Earth/Polar Axis model)");
        std::println("time        {} {:02}:{:02}:{:09.6f} TDB",spp::iso_date(date),hour,minute,second);
        std::println("frame       ICRF to lunar Mean Earth/Polar Axis (IAU_MOON)");
        std::println("pole        RA {:.8f} deg  Dec {:+.8f} deg",result.pole_ra_deg,result.pole_dec_deg);
        std::println("meridian    W {:.8f} deg",result.prime_meridian_deg);
        std::println("pole rate   RA {:+.9f}  Dec {:+.9f} deg/TDB day",result.pole_ra_rate_deg_day,result.pole_dec_rate_deg_day);
        std::println("spin rate   {:+.9f} deg/TDB day",result.prime_meridian_rate_deg_day);
        for(std::size_t i=0;i<3;++i)
            std::println("body {}      {:+.9f} {:+.9f} {:+.9f} (ICRF components)",
                static_cast<char>('X'+i),result.rotation.matrix[i][0],result.rotation.matrix[i][1],result.rotation.matrix[i][2]);
        return 0;
    }
    MoonProvider model(fit,true,integrated,kernel,seed,span);
    const auto state = spp::moon::geometric(spp::Instant<spp::TDB>{jd.day, jd.frac}, model);
    const auto& p = state.position_km;
    const auto& v = state.velocity_km_per_day;
    const double distance = std::hypot(p[0], p[1], p[2]);
    const double speed = std::hypot(v[0], v[1], v[2]) / spp::seconds_per_day;
    const double longitude = spp::norm360(std::atan2(p[1], p[0]) * spp::rad2deg);
    const double latitude = std::atan2(p[2], std::hypot(p[0], p[1])) * spp::rad2deg;

    std::println("Moon — geometric position from Earth's center");
    std::println("time        {} {:02}:{:02}:{:09.6f} TDB", spp::iso_date(date), hour, minute, second);
    model.describe();
    std::println("frame       mean ecliptic and equinox of J2000");
    std::println("distance    {:.3f} km", distance);
    std::println("speed       {:.6f} km/s", speed);
    std::println("longitude   {:.8f} deg", longitude);
    std::println("latitude    {:+.8f} deg", latitude);
    std::println("position    x {:+.6f}  y {:+.6f}  z {:+.6f} km", p[0], p[1], p[2]);
    std::println("velocity    x {:+.9f}  y {:+.9f}  z {:+.9f} km/s",
                 v[0] / spp::seconds_per_day, v[1] / spp::seconds_per_day, v[2] / spp::seconds_per_day);
    return 0;
}
#endif

bool leap(int y, bool gregorian) {
    return gregorian ? (y % 4 == 0 && (y % 100 != 0 || y % 400 == 0)) : (y % 4 == 0);
}

int cmd_year(const std::vector<std::string_view>& a, const Options& o) {
    if (a.size() != 3 && a.size() != 4) usage();
    const int year = integer(a[0]);
    const auto obs = parse_observer(a, 1);
    static constexpr int kDays[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    std::println("date,type,transit_utc,sunrise_utc,sunset_utc,duration_h");
    for (int m = 1; m <= 12; ++m) {
        const int n = kDays[m - 1] + (m == 2 && leap(year, year > 1582) ? 1 : 0);
        for (int dd = 1; dd <= n; ++dd) {
            const spp::CivilDate date{year, m, dd};
            if (year == 1582 && m == 10 && dd > 4 && dd < 15) continue;   // dropped by the reform
            const auto d = spp::events::day(date, obs, o.atm, o.dut1);
            std::println("{},{},{},{},{},{:.6f}", spp::iso_date(date), type_name(d.type), when(d.transit),
                         when(d.rise), when(d.set), d.duration_hours);
        }
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    Options o;
    std::vector<std::string_view> args;
    for (int i = 1; i < argc; ++i) {
        std::string_view s = argv[i];
        if(s=="--help" || s=="-h")usage(0);
        else if (s.starts_with("--dut1=")) o.dut1 = num(s.substr(7));
        else if (s.starts_with("--pressure=")) o.atm.pressure_mbar = num(s.substr(11));
        else if (s.starts_with("--temp=")) o.atm.temp_c = num(s.substr(7));
        else if (s == "--no-refraction") o.atm.refraction = false;
        else args.push_back(s);
    }
    if (args.empty()) usage();

    const std::string_view cmd = args.front();
    args.erase(args.begin());
    try {
        if (cmd == "sun") return cmd_sun(args, o);
        if (cmd == "day") return cmd_day(args, o);
        if (cmd == "year") return cmd_year(args, o);
#if SPP_CLI_MOON
        if (cmd == "moon") return cmd_moon(args, o);
#endif
    } catch (const std::system_error&) {
        return 1;   // output closed early, e.g. piped into `head`
    } catch (const std::invalid_argument& e) {
        std::println(stderr, "spp: {}", e.what());
        return 2;
    } catch (const std::out_of_range& e) {
        std::println(stderr, "spp: {}", e.what());
        return 2;
    } catch (const std::exception& e) {
        std::println(stderr,"spp: {}",e.what());
        return 1;
    }
    usage();
}
