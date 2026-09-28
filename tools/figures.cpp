// spp_figures: renders the README figures as SVG, computed with the Spp library.
//
//   spp_figures [OUTPUT_DIR]      (default: docs/images)
//
// Figures:
//   daylight.svg     night, twilights and daylight across a year (Paris)
//   day_length.svg   day length through the year at several latitudes
//   analemma.svg     the Sun at 12:00 UTC every day of the year (Greenwich)
//   sun_paths.svg    the Sun's path across the sky at solstices and equinox (Paris)
//   spirals.svg      sunrise and sunset azimuths as 3D spirals (Paris)

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <optional>
#include <print>
#include <string>
#include <utility>
#include <vector>

#include "spp/spp.hpp"

namespace {

constexpr int kYear = 2026;
const spp::Observer kParis{48.8566, 2.3522, 35.0};
const spp::Observer kGreenwich{51.4769, -0.0005, 46.0};

// ---------------------------------------------------------------------------
// SVG
// ---------------------------------------------------------------------------

constexpr const char* kFont = "Segoe UI, Helvetica, Arial, sans-serif";
constexpr const char* kInk = "#1f2328";
constexpr const char* kMuted = "#656d76";
constexpr const char* kGrid = "#e6e9ec";
constexpr const char* kAxis = "#8c959f";

class Svg {
public:
    Svg(int width, int height) : width_(width), height_(height) {}

    template <class... Args>
    void add(std::format_string<Args...> fmt, Args&&... args) {
        body_ += std::format(fmt, std::forward<Args>(args)...);
        body_ += '\n';
    }

    void save(const std::filesystem::path& path) const {
        std::ofstream out(path, std::ios::binary);
        out << std::format(R"(<svg xmlns="http://www.w3.org/2000/svg" width="{0}" height="{1}" viewBox="0 0 {0} {1}" font-family="{2}">)",
                           width_, height_, kFont)
            << '\n'
            << body_ << "</svg>\n";
        std::println("wrote {}", path.string());
    }

    int width() const { return width_; }
    int height() const { return height_; }

private:
    int width_;
    int height_;
    std::string body_;
};

// Rounded card behind a figure, with a title and a subtitle.
void card(Svg& svg, const char* fill, const char* stroke, const char* title_color, const std::string& title,
          const std::string& subtitle) {
    svg.add(R"(<rect x="0.5" y="0.5" width="{}" height="{}" rx="12" fill="{}" stroke="{}"/>)", svg.width() - 1,
            svg.height() - 1, fill, stroke);
    svg.add(R"(<text x="24" y="34" font-size="17" font-weight="600" fill="{}">{}</text>)", title_color, title);
    svg.add(R"(<text x="24" y="54" font-size="12.5" fill="{}">{}</text>)", kMuted, subtitle);
}

// Maps data coordinates into a pixel rectangle (y grows upward in data).
struct Frame {
    double left, top, width, height;
    double x_min, x_max, y_min, y_max;

    double x(double v) const { return left + (v - x_min) / (x_max - x_min) * width; }
    double y(double v) const { return top + (y_max - v) / (y_max - y_min) * height; }
    double right() const { return left + width; }
    double bottom() const { return top + height; }
};

using Points = std::vector<std::pair<double, double>>;

std::string polyline(const Points& pts, const Frame& f) {
    std::string s;
    for (const auto& [x, y] : pts) s += std::format("{:.1f},{:.1f} ", f.x(x), f.y(y));
    return s;
}

void clip_to(Svg& svg, const Frame& f, const char* id) {
    svg.add(R"(<clipPath id="{}"><rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}"/></clipPath>)", id, f.left,
            f.top, f.width, f.height);
}

void y_grid(Svg& svg, const Frame& f, double from, double to, double step, const char* unit) {
    for (double v = from; v <= to + 1e-9; v += step) {
        svg.add(R"(<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="{}"/>)", f.left, f.y(v), f.right(),
                f.y(v), kGrid);
        svg.add(R"(<text x="{:.1f}" y="{:.1f}" font-size="11" fill="{}" text-anchor="end">{:g}{}</text>)", f.left - 8,
                f.y(v) + 4, kMuted, v, unit);
    }
}

void x_grid(Svg& svg, const Frame& f, double from, double to, double step, const char* unit) {
    for (double v = from; v <= to + 1e-9; v += step) {
        svg.add(R"(<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="{}"/>)", f.x(v), f.top, f.x(v),
                f.bottom(), kGrid);
        svg.add(R"(<text x="{:.1f}" y="{:.1f}" font-size="11" fill="{}" text-anchor="middle">{:g}{}</text>)", f.x(v),
                f.bottom() + 18, kMuted, v, unit);
    }
}

void frame_border(Svg& svg, const Frame& f) {
    svg.add(R"(<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}" fill="none" stroke="{}"/>)", f.left, f.top,
            f.width, f.height, kAxis);
}

void axis_label_y(Svg& svg, const Frame& f, const std::string& text) {
    svg.add(R"svg(<text transform="translate({:.1f},{:.1f}) rotate(-90)" font-size="12" fill="{}" text-anchor="middle">{}</text>)svg",
            f.left - 44, f.top + f.height / 2, kInk, text);
}

void axis_label_x(Svg& svg, const Frame& f, const std::string& text) {
    svg.add(R"(<text x="{:.1f}" y="{:.1f}" font-size="12" fill="{}" text-anchor="middle">{}</text>)",
            f.left + f.width / 2, f.bottom() + 40, kInk, text);
}

// ---------------------------------------------------------------------------
// Calendar and astronomy helpers
// ---------------------------------------------------------------------------

struct Year {
    std::vector<spp::CivilDate> dates;   // every civil date, index = day of year - 1
    std::array<int, 12> month_start{};   // day of year of the 1st of each month
};

Year year_dates(int year) {
    Year y;
    const auto jan1 = spp::Time::from_utc(year, 1, 1);
    for (auto t = jan1;; t = t.shifted(1.0)) {
        const auto c = t.to_civil();
        if (c.year != year) break;
        y.dates.push_back(c.date());
        if (c.day == 1) y.month_start[static_cast<std::size_t>(c.month - 1)] = static_cast<int>(y.dates.size());
    }
    return y;
}

// Month grid lines and centered month names along a day-of-year axis.
void month_axis(Svg& svg, const Frame& f, const Year& y) {
    static constexpr std::array<const char*, 12> names{"Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
    const int days = static_cast<int>(y.dates.size());
    for (std::size_t m = 0; m < 12; ++m) {
        const double start = y.month_start[m];
        const double end = m + 1 < 12 ? y.month_start[m + 1] : days + 1;
        svg.add(R"(<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="{}"/>)", f.x(start), f.top, f.x(start),
                f.bottom(), kGrid);
        svg.add(R"(<text x="{:.1f}" y="{:.1f}" font-size="11" fill="{}" text-anchor="middle">{}</text>)",
                f.x((start + end) / 2), f.bottom() + 18, kMuted, names[m]);
    }
}

double hours_since(const spp::Time& t, const spp::Time& ref) {
    return (static_cast<double>(t.utc().day - ref.utc().day) + (t.utc().frac - ref.utc().frac)) * 24.0;
}

// Topocentric position of the Sun, with the altitude refracted.
spp::Horizontal sun_sky(const spp::Time& t, const spp::Observer& obs) {
    const auto tt = t.tt();
    const auto n = spp::earth::nutation(tt);
    const auto eq = spp::sun::apparent(tt, n);
    const double lst = spp::earth::lst_deg(t.ut1(), n, obs.lon_deg);
    const auto hor = spp::sky::to_horizontal(spp::sky::parallax(eq, obs, lst), obs.lat_deg, lst);
    return {spp::sky::refraction(hor.alt_deg, {}), hor.az_deg};
}

// Blue around January, orange around July (as in the spp_spiral viewer).
std::string season_color(int doy, int days) {
    const double s = 0.5 - 0.5 * std::cos(2.0 * spp::pi * (doy - 1) / days);
    const auto mix = [s](double a, double b) { return static_cast<int>(std::lround(255.0 * (a + (b - a) * s))); };
    return std::format("#{:02x}{:02x}{:02x}", mix(0.30, 1.00), mix(0.60, 0.62), mix(1.00, 0.20));
}

// ---------------------------------------------------------------------------
// 1. Daylight and twilight across the year
// ---------------------------------------------------------------------------

// Hours (UTC, from the date's midnight) between which the Sun's center is above
// `target` on the day of `d`. Spans the whole day when it never sets below it.
std::pair<double, double> above(const spp::events::Day& d, const spp::Time& midnight, double target) {
    using spp::events::Direction;
    const double tr = hours_since(d.transit, midnight);
    if (spp::events::altitude_deg(d.transit, kParis, {}) < target) return {tr, tr};
    const auto dawn = spp::events::crossing(d.transit, kParis, {}, target, Direction::Before);
    const auto dusk = spp::events::crossing(d.transit, kParis, {}, target, Direction::After);
    return {dawn ? hours_since(*dawn, midnight) : tr - 12.0, dusk ? hours_since(*dusk, midnight) : tr + 12.0};
}

void render_daylight(const std::filesystem::path& dir) {
    const Year y = year_dates(kYear);
    const int days = static_cast<int>(y.dates.size());

    struct Level {
        double target;
        const char* color;
        const char* name;
    };
    const std::array<Level, 4> levels{{
        {spp::events::astronomical_twilight_deg, "#23406e", "Astronomical twilight"},
        {spp::events::nautical_twilight_deg, "#3a6aa6", "Nautical twilight"},
        {spp::events::civil_twilight_deg, "#79a7d9", "Civil twilight"},
        {spp::events::upper_limb_deg, "#fbd97a", "Daylight"},
    }};

    std::array<Points, 4> dawn, dusk;
    Points noon;
    for (int i = 0; i < days; ++i) {
        const auto& date = y.dates[static_cast<std::size_t>(i)];
        const auto midnight = spp::Time::from_utc(date.year, date.month, date.day);
        const auto d = spp::events::day(date, kParis);
        noon.push_back({i + 1.0, hours_since(d.transit, midnight)});
        for (std::size_t l = 0; l < levels.size(); ++l) {
            const auto [a, b] = above(d, midnight, levels[l].target);
            dawn[l].push_back({i + 1.0, a});
            dusk[l].push_back({i + 1.0, b});
        }
    }

    Svg svg(920, 540);
    card(svg, "#ffffff", "#d0d7de", kInk, "Daylight and twilight through the year",
         std::format("Paris (48.86° N, 2.35° E), {}. Hours in UTC. Every boundary is an event solved by spp::events.",
                     kYear));
    const Frame f{70, 76, 820, 380, 1, static_cast<double>(days + 1), 0, 24};
    clip_to(svg, f, "daylight-clip");
    svg.add(R"(<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}" fill="#0f1f3d"/>)", f.left, f.top, f.width,
            f.height);
    svg.add(R"svg(<g clip-path="url(#daylight-clip)">)svg");
    for (std::size_t l = 0; l < levels.size(); ++l) {
        Points band = dusk[l];
        band.insert(band.end(), dawn[l].rbegin(), dawn[l].rend());
        svg.add(R"(<polygon points="{}" fill="{}"/>)", polyline(band, f), levels[l].color);
    }
    svg.add(R"(<polyline points="{}" fill="none" stroke="#c26a00" stroke-width="1.6" stroke-dasharray="6 4"/>)",
            polyline(noon, f));
    svg.add("</g>");

    for (int h = 0; h <= 24; h += 3)
        svg.add(R"(<text x="{:.1f}" y="{:.1f}" font-size="11" fill="{}" text-anchor="end">{:02d}:00</text>)",
                f.left - 8, f.y(h) + 4, kMuted, h);
    month_axis(svg, f, y);
    frame_border(svg, f);
    axis_label_y(svg, f, "Time of day (UTC)");

    // Legend
    double x = f.left;
    const double ly = f.bottom() + 44;
    auto swatch = [&](const char* color, const char* name, bool dashed) {
        if (dashed)
            svg.add(R"(<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="#c26a00" stroke-width="1.6" stroke-dasharray="6 4"/>)",
                    x, ly - 4, x + 22, ly - 4);
        else
            svg.add(R"(<rect x="{:.1f}" y="{:.1f}" width="22" height="12" rx="2" fill="{}" stroke="#d0d7de"/>)", x,
                    ly - 10, color);
        svg.add(R"(<text x="{:.1f}" y="{:.1f}" font-size="12" fill="{}">{}</text>)", x + 30, ly, kInk, name);
        x += 44 + 7.0 * static_cast<double>(std::char_traits<char>::length(name));
    };
    swatch("#0f1f3d", "Night", false);
    for (const auto& l : levels) swatch(l.color, l.name, false);
    swatch("", "Solar noon", true);

    svg.save(dir / "daylight.svg");
}

// ---------------------------------------------------------------------------
// 2. Day length by latitude
// ---------------------------------------------------------------------------

void render_day_length(const std::filesystem::path& dir) {
    const Year y = year_dates(kYear);
    const int days = static_cast<int>(y.dates.size());
    const std::array<double, 9> lats{0, 15, 30, 45, 55, 60, 65, 70, 80};
    const std::array<const char*, 9> colors{"#d4572a", "#e0822f", "#d9a62e", "#8aa83a", "#3f9b6e",
                                            "#2b8a9e", "#3070b8", "#4f55b8", "#6b3fa0"};

    Svg svg(920, 520);
    card(svg, "#ffffff", "#d0d7de", kInk, "Day length through the year, by latitude",
         std::format("Sunrise to sunset (upper limb, refracted), {}. Polar day counts 24 h, polar night 0 h.", kYear));
    const Frame f{70, 76, 700, 370, 1, static_cast<double>(days + 1), 0, 24};
    y_grid(svg, f, 0, 24, 3, " h");
    month_axis(svg, f, y);

    for (std::size_t k = 0; k < lats.size(); ++k) {
        const spp::Observer obs{lats[k], 0.0, 0.0};
        // A normal day with a missing rise or set has no defined length
        // (see Limitations in the README): leave a gap instead of plotting 0 h.
        std::vector<Points> runs(1);
        for (int i = 0; i < days; ++i) {
            const auto d = spp::events::day(y.dates[static_cast<std::size_t>(i)], obs);
            if (d.type == spp::events::DayType::Normal && (!d.rise || !d.set)) {
                runs.emplace_back();
                continue;
            }
            runs.back().push_back({i + 1.0, d.duration_hours});
        }
        for (const auto& run : runs)
            if (run.size() > 1)
                svg.add(R"(<polyline points="{}" fill="none" stroke="{}" stroke-width="2" stroke-linejoin="round"/>)",
                        polyline(run, f), colors[k]);
        const double ly = f.top + 14 + 26.0 * static_cast<double>(k);
        svg.add(R"(<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="{}" stroke-width="3"/>)",
                f.right() + 28, ly - 4, f.right() + 52, ly - 4, colors[k]);
        svg.add(R"(<text x="{:.1f}" y="{:.1f}" font-size="12" fill="{}">{:g}° N</text>)", f.right() + 60, ly, kInk,
                lats[k]);
    }
    frame_border(svg, f);
    axis_label_y(svg, f, "Day length");
    svg.save(dir / "day_length.svg");
}

// ---------------------------------------------------------------------------
// 3. Analemma
// ---------------------------------------------------------------------------

void render_analemma(const std::filesystem::path& dir) {
    const Year y = year_dates(kYear);
    const int days = static_cast<int>(y.dates.size());

    Svg svg(560, 660);
    card(svg, "#ffffff", "#d0d7de", kInk, "Analemma",
         std::format("The Sun at 12:00 UTC every day of {}, seen from Greenwich.", kYear));
    const Frame f{78, 76, 440, 510, 170, 190, 10, 65};
    x_grid(svg, f, 170, 190, 5, "°");
    y_grid(svg, f, 10, 65, 5, "°");

    std::vector<spp::Horizontal> sky;
    for (const auto& d : y.dates) sky.push_back(sun_sky(spp::Time::from_utc(d.year, d.month, d.day, 12), kGreenwich));
    for (int i = 0; i < days; ++i) {
        const auto& h = sky[static_cast<std::size_t>(i)];
        svg.add(R"(<circle cx="{:.1f}" cy="{:.1f}" r="3.2" fill="{}"/>)", f.x(h.az_deg), f.y(h.alt_deg),
                season_color(i + 1, days));
    }

    // Solstices and equinoxes of 2026 (dates of the events in UTC).
    struct Mark {
        int month, day;
        const char* label;
        double dx;
    };
    for (const auto& m : std::array<Mark, 4>{{{3, 20, "Mar 20", 12}, {6, 21, "Jun 21", 12}, {9, 23, "Sep 23", -12},
                                              {12, 21, "Dec 21", 12}}}) {
        const auto h = sun_sky(spp::Time::from_utc(kYear, m.month, m.day, 12), kGreenwich);
        svg.add(R"(<circle cx="{:.1f}" cy="{:.1f}" r="5.5" fill="none" stroke="{}" stroke-width="1.5"/>)",
                f.x(h.az_deg), f.y(h.alt_deg), kInk);
        svg.add(R"(<text x="{:.1f}" y="{:.1f}" font-size="12" fill="{}" text-anchor="{}">{}</text>)",
                f.x(h.az_deg) + m.dx, f.y(h.alt_deg) + 4, kInk, m.dx > 0 ? "start" : "end", m.label);
    }
    frame_border(svg, f);
    axis_label_y(svg, f, "Altitude");
    axis_label_x(svg, f, "Azimuth (180° = due South)");
    svg.save(dir / "analemma.svg");
}

// ---------------------------------------------------------------------------
// 4. Sun paths
// ---------------------------------------------------------------------------

void render_sun_paths(const std::filesystem::path& dir) {
    Svg svg(920, 470);
    card(svg, "#ffffff", "#d0d7de", kInk, "The Sun's path across the sky",
         std::format("Paris, {}. Refracted altitude against azimuth; dots every hour (UTC).", kYear));
    const Frame f{70, 76, 820, 300, 30, 330, -2, 70};
    x_grid(svg, f, 30, 330, 30, "°");
    y_grid(svg, f, 0, 70, 10, "°");
    for (auto [az, name] : std::array<std::pair<double, const char*>, 3>{{{90, "East"}, {180, "South"}, {270, "West"}}})
        svg.add(R"(<text x="{:.1f}" y="{:.1f}" font-size="11" font-weight="600" fill="{}" text-anchor="middle">{}</text>)",
                f.x(az), f.top - 6, kMuted, name);
    svg.add(R"(<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}" fill="#eef1f4"/>)", f.left, f.y(0),
            f.width, f.y(-2) - f.y(0));

    struct Path {
        int month, day;
        const char* color;
        const char* name;
    };
    const std::array<Path, 3> paths{{{6, 21, "#e0822f", "June solstice"},
                                     {3, 20, "#3f9b6e", "March equinox"},
                                     {12, 21, "#3070b8", "December solstice"}}};
    for (std::size_t k = 0; k < paths.size(); ++k) {
        const auto& p = paths[k];
        const auto midnight = spp::Time::from_utc(kYear, p.month, p.day);
        Points pts;
        for (int m = 0; m <= 24 * 60; m += 4) {
            const auto h = sun_sky(midnight.shifted(m / 1440.0), kParis);
            if (h.alt_deg > -2.0) pts.push_back({h.az_deg, h.alt_deg});
        }
        svg.add(R"(<polyline points="{}" fill="none" stroke="{}" stroke-width="2.4" stroke-linejoin="round"/>)",
                polyline(pts, f), p.color);
        for (int hour = 0; hour < 24; ++hour) {
            const auto h = sun_sky(midnight.shifted(hour / 24.0), kParis);
            if (h.alt_deg < 0.0) continue;
            svg.add(R"(<circle cx="{:.1f}" cy="{:.1f}" r="3.4" fill="#ffffff" stroke="{}" stroke-width="1.8"/>)",
                    f.x(h.az_deg), f.y(h.alt_deg), p.color);
            if (k == 0)
                svg.add(R"(<text x="{:.1f}" y="{:.1f}" font-size="10" fill="{}" text-anchor="middle">{}h</text>)",
                        f.x(h.az_deg), f.y(h.alt_deg) - 8, kMuted, hour);
        }
        const double lx = f.left + 230.0 * static_cast<double>(k);
        svg.add(R"(<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="{}" stroke-width="3"/>)", lx,
                f.bottom() + 46, lx + 24, f.bottom() + 46, p.color);
        svg.add(R"(<text x="{:.1f}" y="{:.1f}" font-size="12" fill="{}">{}</text>)", lx + 32, f.bottom() + 50, kInk,
                p.name);
    }
    frame_border(svg, f);
    axis_label_y(svg, f, "Altitude");
    svg.save(dir / "sun_paths.svg");
}

// ---------------------------------------------------------------------------
// 5. Sunrise and sunset azimuth spirals
// ---------------------------------------------------------------------------

struct Vec3 {
    double x, y, z;
};

Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3 cross(Vec3 a, Vec3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
Vec3 normalize(Vec3 a) {
    const double n = std::sqrt(dot(a, a));
    return {a.x / n, a.y / n, a.z / n};
}

// Same geometry as spp_spiral: x = East, -z = North, y = day of year.
constexpr double kHeight = 3.0;

Vec3 cyl(double az_deg, double h, double r = 1.0) {
    const double a = az_deg * spp::deg2rad;
    return {r * std::sin(a), h, -r * std::cos(a)};
}

// Orthographic view from the viewer's default camera (yaw 0.7, pitch 0.35).
struct Camera {
    Vec3 target{0.0, kHeight / 2.0, 0.0};
    Vec3 forward, right, up;
    double cx, cy, scale;

    Camera(double center_x, double center_y, double s) : cx(center_x), cy(center_y), scale(s) {
        const double yaw = 0.7, pitch = 0.35;
        forward = normalize({-std::cos(pitch) * std::sin(yaw), -std::sin(pitch), -std::cos(pitch) * std::cos(yaw)});
        right = normalize(cross(forward, {0.0, 1.0, 0.0}));
        up = cross(right, forward);
    }
    double sx(Vec3 p) const { return cx + scale * dot(p - target, right); }
    double sy(Vec3 p) const { return cy - scale * dot(p - target, up); }
    double depth(Vec3 p) const { return dot(p - target, forward); }
};

struct Segment {
    Vec3 a, b;
    std::string color;
    double width;
    double base_opacity;
};

void spiral_panel(Svg& svg, const Camera& cam, const std::vector<std::pair<int, double>>& samples, int days,
                  const Year& y, const char* title) {
    const double top = kHeight * (days - 1) / std::max(days - 1, 1);
    const auto height = [&](int doy) { return kHeight * (doy - 1) / std::max(days - 1, 1); };
    std::vector<Segment> segs;

    constexpr int n = 96;
    auto ring = [&](double h, const char* color, double opacity) {
        for (int i = 0; i < n; ++i)
            segs.push_back({cyl(360.0 * i / n, h), cyl(360.0 * (i + 1) / n, h), color, 1.0, opacity});
    };
    ring(0.0, "#6b7280", 0.9);
    ring(top, "#6b7280", 0.9);
    for (std::size_t m = 1; m < 12; ++m) ring(height(y.month_start[m]), "#3a4150", 0.8);
    for (int az = 0; az < 360; az += 90) {
        const char* color = az == 0 ? "#d94c4c" : "#4b5563";
        segs.push_back({cyl(az, 0.0), cyl(az, top), color, 1.0, 0.9});
        segs.push_back({{0.0, 0.0, 0.0}, cyl(az, 0.0, 1.15), color, 1.2, 0.9});
    }
    segs.push_back({{0.0, 0.0, 0.0}, {0.0, top, 0.0}, "#6b7280", 1.0, 0.9});
    for (std::size_t i = 1; i < samples.size(); ++i) {
        if (samples[i].first != samples[i - 1].first + 1) continue;
        segs.push_back({cyl(samples[i - 1].second, height(samples[i - 1].first)),
                        cyl(samples[i].second, height(samples[i].first)), season_color(samples[i].first, days), 2.6,
                        1.0});
    }

    // Painter's order, far to near; farther segments fade.
    std::ranges::sort(segs, [&](const Segment& s, const Segment& t) {
        return cam.depth(s.a) + cam.depth(s.b) > cam.depth(t.a) + cam.depth(t.b);
    });
    for (const auto& s : segs) {
        const double d = (cam.depth(s.a) + cam.depth(s.b)) / 2.0;
        const double fade = std::clamp(0.62 - 0.28 * d, 0.3, 1.0);
        svg.add(R"(<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="{}" stroke-width="{:g}" stroke-linecap="round" opacity="{:.2f}"/>)",
                cam.sx(s.a), cam.sy(s.a), cam.sx(s.b), cam.sy(s.b), s.color, s.width, s.base_opacity * fade);
    }

    for (auto [c, az] : std::array<std::pair<const char*, int>, 4>{{{"N", 0}, {"E", 90}, {"S", 180}, {"W", 270}}}) {
        const Vec3 p = cyl(az, 0.0, 1.36);
        svg.add(R"(<text x="{:.1f}" y="{:.1f}" font-size="14" font-weight="600" fill="{}" text-anchor="middle">{}</text>)",
                cam.sx(p), cam.sy(p) + 5, az == 0 ? "#f07171" : "#c9d1d9", c);
    }
    svg.add(R"(<text x="{:.1f}" y="92" font-size="15" font-weight="600" fill="#e6edf3" text-anchor="middle">{}</text>)",
            cam.cx, title);
}

void render_spirals(const std::filesystem::path& dir) {
    const Year y = year_dates(kYear);
    const int days = static_cast<int>(y.dates.size());

    std::vector<std::pair<int, double>> rise, set;
    for (int i = 0; i < days; ++i) {
        const auto d = spp::events::day(y.dates[static_cast<std::size_t>(i)], kParis);
        if (d.rise) rise.push_back({i + 1, sun_sky(*d.rise, kParis).az_deg});
        if (d.set) set.push_back({i + 1, sun_sky(*d.set, kParis).az_deg});
    }

    Svg svg(920, 640);
    card(svg, "#0d1117", "#30363d", "#e6edf3", "Sunrise and sunset azimuth through the year",
         std::format("Paris, {}. Angle = azimuth of the Sun at the event (N in red), height = day of the year "
                     "(1 January at the bottom), as in spp_spiral.",
                     kYear));
    spiral_panel(svg, Camera(240, 360, 132), rise, days, y, "Sunrise");
    spiral_panel(svg, Camera(680, 360, 132), set, days, y, "Sunset");
    svg.save(dir / "spirals.svg");
}

}  // namespace

int main(int argc, char** argv) {
    const std::filesystem::path dir = argc > 1 ? argv[1] : "docs/images";
    std::filesystem::create_directories(dir);
    render_daylight(dir);
    render_day_length(dir);
    render_analemma(dir);
    render_sun_paths(dir);
    render_spirals(dir);
    return 0;
}
