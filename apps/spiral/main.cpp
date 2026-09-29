// spp_spiral: sunrise and sunset azimuths over a year, as two 3D spirals.
//
//   spp_spiral YEAR LAT LON [ELEV]
//
// Two windows, one per event. Each day of the year is a point in cylindrical
// coordinates: angle = azimuth of the Sun at the event (0 = North, 90 = East),
// height = day of the year (1 January at the bottom), radius = 1.
//
// Mouse: left-drag to orbit, wheel to zoom. Keys: T top view, R reset view,
// Esc close the window.

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <print>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "season_color.hpp"
#include "spp/spp.hpp"

namespace {

constexpr float kPi = 3.14159265358979f;
constexpr float kHeight = 3.0f;   // length of the year axis, in units of the radius

// ---------------------------------------------------------------------------
// Astronomy: everything below comes from the Spp library
// ---------------------------------------------------------------------------

struct Sample {
    int doy;               // day of the year, 1 = 1 January
    spp::CivilDate date;
    double az_deg;         // azimuth of the Sun at the event, 0 = North, 90 = East
};

struct YearData {
    int days = 0;
    std::vector<int> month_start;   // day of the year of the 1st of each month
    std::vector<Sample> rise, set;
};

YearData compute(int year, const spp::Observer& obs) {
    YearData data;
    const auto jan1 = spp::Time::from_utc(year, 1, 1);
    for (int m = 1; m <= 12; ++m)
        data.month_start.push_back(static_cast<int>(spp::Time::from_utc(year, m, 1).utc().day - jan1.utc().day) + 1);

    // Walk the civil calendar day by day (handles leap years and the 1582 reform).
    for (auto t = jan1;; t = t.shifted(1.0)) {
        const auto c = t.to_civil();
        if (c.year != year) break;
        const int doy = ++data.days;
        const auto day = spp::events::day(c.date(), obs);
        if (day.rise) data.rise.push_back({doy, c.date(), spp::sun::horizontal(*day.rise, obs).az_deg});
        if (day.set) data.set.push_back({doy, c.date(), spp::sun::horizontal(*day.set, obs).az_deg});
    }
    return data;
}

void print_summary(const char* name, const std::vector<Sample>& s, const YearData& data, int year) {
    if (s.empty()) {
        std::println("{}: none in {}", name, year);
        return;
    }
    const auto [lo, hi] = std::minmax_element(s.begin(), s.end(),
                                              [](const Sample& a, const Sample& b) { return a.az_deg < b.az_deg; });
    auto iso = [](const spp::CivilDate& d) { return std::format("{:04d}-{:02d}-{:02d}", d.year, d.month, d.day); };
    std::println("{}: {} days, azimuth {:.2f} deg ({}) to {:.2f} deg ({})", name, s.size(), lo->az_deg,
                 iso(lo->date), hi->az_deg, iso(hi->date));
    if (static_cast<int>(s.size()) < data.days)
        std::println("  {} days without {}", data.days - static_cast<int>(s.size()), name);
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

struct Vertex {
    glm::vec3 pos;
    glm::vec3 col;
};

// World axes: x = East, -z = North, y = up (the day-of-year axis).
glm::vec3 cyl(double az_deg, float h, float r = 1.0f) {
    const float a = static_cast<float>(az_deg) * (kPi / 180.0f);
    return {r * std::sin(a), h, -r * std::cos(a)};
}

float height(int doy, int days) {
    return kHeight * static_cast<float>(doy - 1) / static_cast<float>(std::max(days - 1, 1));
}

glm::vec3 season_color(int doy, int days) {
    const auto [r, g, b] = common::season_rgb(doy, days);
    return {static_cast<float>(r), static_cast<float>(g), static_cast<float>(b)};
}

const glm::vec3 kNorth(0.85f, 0.30f, 0.30f);
const glm::vec3 kGuide(0.42f, 0.44f, 0.48f);
const glm::vec3 kFaint(0.20f, 0.21f, 0.24f);

void segment(std::vector<Vertex>& v, glm::vec3 a, glm::vec3 b, glm::vec3 col) {
    v.push_back({a, col});
    v.push_back({b, col});
}

void ring(std::vector<Vertex>& v, float h, float r, glm::vec3 col) {
    constexpr int n = 128;
    for (int i = 0; i < n; ++i)
        segment(v, cyl(360.0 * i / n, h, r), cyl(360.0 * (i + 1) / n, h, r), col);
}

// Cylinder guides: rings at the ends and on the 1st of each month, the four
// cardinal lines along the year axis, the central axis and compass spokes.
std::vector<Vertex> guides(const YearData& data) {
    std::vector<Vertex> v;
    const float top = height(data.days, data.days);
    ring(v, 0.0f, 1.0f, kGuide);
    ring(v, top, 1.0f, kGuide);
    for (std::size_t m = 1; m < data.month_start.size(); ++m) ring(v, height(data.month_start[m], data.days), 1.0f, kFaint);

    for (int az = 0; az < 360; az += 90) {
        const glm::vec3 col = az == 0 ? kNorth : kGuide;
        segment(v, cyl(az, 0.0f), cyl(az, top), az == 0 ? kNorth : kFaint);
        segment(v, glm::vec3(0.0f), cyl(az, 0.0f, 1.15f), col);
        segment(v, glm::vec3(0.0f, top, 0.0f), cyl(az, top, 1.15f), col);
    }
    segment(v, glm::vec3(0.0f), glm::vec3(0.0f, top, 0.0f), kGuide);
    return v;
}

// Stroke glyphs for the compass letters, in a unit box.
using Stroke = std::array<float, 4>;

std::vector<Stroke> glyph(char c) {
    switch (c) {
        case 'N': return {{0, 0, 0, 1}, {0, 1, 1, 0}, {1, 0, 1, 1}};
        case 'E': return {{1, 1, 0, 1}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0.5f, 0.75f, 0.5f}};
        case 'S': return {{1, 1, 0, 1}, {0, 1, 0, 0.5f}, {0, 0.5f, 1, 0.5f}, {1, 0.5f, 1, 0}, {1, 0, 0, 0}};
        case 'W': return {{0, 1, 0.25f, 0}, {0.25f, 0, 0.5f, 0.6f}, {0.5f, 0.6f, 0.75f, 0}, {0.75f, 0, 1, 1}};
        default: return {};
    }
}

// Compass letters facing the camera, around both ends of the cylinder.
std::vector<Vertex> labels(float top, glm::vec3 right, glm::vec3 up) {
    std::vector<Vertex> v;
    constexpr float size = 0.14f;
    constexpr std::array<std::pair<char, int>, 4> letters{{{'N', 0}, {'E', 90}, {'S', 180}, {'W', 270}}};
    for (float h : {0.0f, top}) {
        for (auto [c, az] : letters) {
            const glm::vec3 at = cyl(az, h, 1.32f);
            const glm::vec3 col = c == 'N' ? kNorth : glm::vec3(0.85f);
            for (auto [x0, y0, x1, y1] : glyph(c)) {
                auto p = [&](float x, float y) { return at + size * ((x - 0.5f) * 0.7f * right + (y - 0.5f) * up); };
                segment(v, p(x0, y0), p(x1, y1), col);
            }
        }
    }
    return v;
}

// ---------------------------------------------------------------------------
// OpenGL
// ---------------------------------------------------------------------------

constexpr const char* kVertexShader = R"(#version 330 core
layout(location = 0) in vec3 pos;
layout(location = 1) in vec3 col;
uniform mat4 mvp;
uniform float point_size;
out vec3 v_col;
void main() {
    gl_Position = mvp * vec4(pos, 1.0);
    gl_PointSize = point_size;
    v_col = col;
}
)";

constexpr const char* kFragmentShader = R"(#version 330 core
in vec3 v_col;
uniform bool round_points;
out vec4 frag;
void main() {
    if (round_points) {
        vec2 c = gl_PointCoord * 2.0 - 1.0;
        if (dot(c, c) > 1.0) discard;
    }
    frag = vec4(v_col, 1.0);
}
)";

GLuint compile(GLenum type, const char* src) {
    const GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok = 0;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        std::array<char, 2048> log{};
        glGetShaderInfoLog(s, static_cast<GLsizei>(log.size()), nullptr, log.data());
        std::println(stderr, "spp_spiral: shader error: {}", log.data());
        std::exit(1);
    }
    return s;
}

GLuint program() {
    const GLuint vs = compile(GL_VERTEX_SHADER, kVertexShader);
    const GLuint fs = compile(GL_FRAGMENT_SHADER, kFragmentShader);
    const GLuint p = glCreateProgram();
    glAttachShader(p, vs);
    glAttachShader(p, fs);
    glLinkProgram(p);
    GLint ok = 0;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) {
        std::array<char, 2048> log{};
        glGetProgramInfoLog(p, static_cast<GLsizei>(log.size()), nullptr, log.data());
        std::println(stderr, "spp_spiral: link error: {}", log.data());
        std::exit(1);
    }
    glDeleteShader(vs);
    glDeleteShader(fs);
    return p;
}

struct Mesh {
    GLuint vao = 0;
    GLuint vbo = 0;
    GLsizei count = 0;

    void create() {
        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        const auto offset = [](std::size_t o) { return reinterpret_cast<const void*>(static_cast<std::uintptr_t>(o)); };
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), offset(offsetof(Vertex, pos)));
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), offset(offsetof(Vertex, col)));
        glEnableVertexAttribArray(1);
    }
    void upload(const std::vector<Vertex>& v) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(v.size() * sizeof(Vertex)), v.data(), GL_DYNAMIC_DRAW);
        count = static_cast<GLsizei>(v.size());
    }
    void destroy() {
        glDeleteBuffers(1, &vbo);
        glDeleteVertexArrays(1, &vao);
    }
};

// One window: one event (sunrise or sunset) and its own camera.
struct Plot {
    GLFWwindow* win = nullptr;
    GLuint prog = 0;
    Mesh guides, curve, labels;
    std::vector<GLint> firsts;     // one line strip per run of consecutive days
    std::vector<GLsizei> counts;
    float top = kHeight;

    float yaw = 0.0f, pitch = 0.0f, dist = 0.0f;
    bool dragging = false;
    double last_x = 0.0, last_y = 0.0;

    void reset_view() {
        yaw = 0.7f;
        pitch = 0.35f;
        dist = 6.5f;
    }
};

Plot& plot_of(GLFWwindow* w) { return *static_cast<Plot*>(glfwGetWindowUserPointer(w)); }

void install_callbacks(Plot& p) {
    glfwSetWindowUserPointer(p.win, &p);
    glfwSetMouseButtonCallback(p.win, [](GLFWwindow* w, int button, int action, int) {
        if (button != GLFW_MOUSE_BUTTON_LEFT) return;
        Plot& q = plot_of(w);
        q.dragging = action == GLFW_PRESS;
        glfwGetCursorPos(w, &q.last_x, &q.last_y);
    });
    glfwSetCursorPosCallback(p.win, [](GLFWwindow* w, double x, double y) {
        Plot& q = plot_of(w);
        if (q.dragging) {
            q.yaw -= static_cast<float>(x - q.last_x) * 0.008f;
            q.pitch = std::clamp(q.pitch + static_cast<float>(y - q.last_y) * 0.008f, -1.55f, 1.55f);
        }
        q.last_x = x;
        q.last_y = y;
    });
    glfwSetScrollCallback(p.win, [](GLFWwindow* w, double, double dy) {
        Plot& q = plot_of(w);
        q.dist = std::clamp(q.dist * std::pow(0.9f, static_cast<float>(dy)), 1.5f, 30.0f);
    });
    glfwSetKeyCallback(p.win, [](GLFWwindow* w, int key, int, int action, int) {
        if (action != GLFW_PRESS) return;
        Plot& q = plot_of(w);
        if (key == GLFW_KEY_ESCAPE) glfwSetWindowShouldClose(w, GLFW_TRUE);
        if (key == GLFW_KEY_R) q.reset_view();
        if (key == GLFW_KEY_T) {   // from above: North up, East right
            q.yaw = 0.0f;
            q.pitch = 1.55f;
        }
    });
}

void setup(Plot& p, const YearData& data, const std::vector<Sample>& samples) {
    glfwMakeContextCurrent(p.win);
    glfwSwapInterval(0);
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::println(stderr, "spp_spiral: cannot load OpenGL functions");
        std::exit(1);
    }
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_PROGRAM_POINT_SIZE);

    p.prog = program();
    p.guides.create();
    p.curve.create();
    p.labels.create();
    p.top = height(data.days, data.days);
    p.guides.upload(guides(data));

    // Consecutive days form one strip; days without the event split the curve.
    std::vector<Vertex> v;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (i == 0 || samples[i].doy != samples[i - 1].doy + 1) {
            p.firsts.push_back(static_cast<GLint>(v.size()));
            p.counts.push_back(0);
        }
        const Sample& s = samples[i];
        v.push_back({cyl(s.az_deg, height(s.doy, data.days)), season_color(s.doy, data.days)});
        ++p.counts.back();
    }
    p.curve.upload(v);

    p.reset_view();
    install_callbacks(p);
}

void render(Plot& p) {
    glfwMakeContextCurrent(p.win);
    int fw = 0, fh = 0;
    glfwGetFramebufferSize(p.win, &fw, &fh);
    if (fw == 0 || fh == 0) return;   // minimized

    glViewport(0, 0, fw, fh);
    glClearColor(0.07f, 0.08f, 0.10f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    const glm::vec3 target(0.0f, p.top / 2.0f, 0.0f);
    const glm::vec3 eye = target + p.dist * glm::vec3(std::cos(p.pitch) * std::sin(p.yaw), std::sin(p.pitch),
                                                      std::cos(p.pitch) * std::cos(p.yaw));
    const glm::mat4 view = glm::lookAt(eye, target, glm::vec3(0.0f, 1.0f, 0.0f));
    const glm::mat4 proj = glm::perspective(glm::radians(40.0f), static_cast<float>(fw) / static_cast<float>(fh),
                                            0.05f, 100.0f);
    const glm::mat4 mvp = proj * view;

    glUseProgram(p.prog);
    glUniformMatrix4fv(glGetUniformLocation(p.prog, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
    const GLint point_size = glGetUniformLocation(p.prog, "point_size");
    const GLint round_points = glGetUniformLocation(p.prog, "round_points");

    glUniform1f(point_size, 1.0f);
    glUniform1i(round_points, 0);
    glBindVertexArray(p.guides.vao);
    glDrawArrays(GL_LINES, 0, p.guides.count);

    const glm::vec3 forward = glm::normalize(target - eye);
    const glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0.0f, 1.0f, 0.0f)));
    const glm::vec3 up = glm::cross(right, forward);
    p.labels.upload(labels(p.top, right, up));
    glBindVertexArray(p.labels.vao);
    glDrawArrays(GL_LINES, 0, p.labels.count);

    glBindVertexArray(p.curve.vao);
    if (!p.firsts.empty())
        glMultiDrawArrays(GL_LINE_STRIP, p.firsts.data(), p.counts.data(), static_cast<GLsizei>(p.firsts.size()));
    glUniform1f(point_size, 5.0f);
    glUniform1i(round_points, 1);
    glDrawArrays(GL_POINTS, 0, p.curve.count);

    glfwSwapBuffers(p.win);
}

void destroy(Plot& p) {
    glfwMakeContextCurrent(p.win);
    p.guides.destroy();
    p.curve.destroy();
    p.labels.destroy();
    glDeleteProgram(p.prog);
    glfwDestroyWindow(p.win);
    p.win = nullptr;
}

// ---------------------------------------------------------------------------
// Command line
// ---------------------------------------------------------------------------

[[noreturn]] void usage() {
    std::println(stderr, "usage: spp_spiral YEAR LAT LON [ELEV]");
    std::exit(2);
}

template <class T>
T parse(std::string_view s) {
    T v{};
    auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc{} || p != s.data() + s.size()) usage();
    return v;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) usage();
    const int year = parse<int>(argv[1]);
    const spp::Observer obs{parse<double>(argv[2]), parse<double>(argv[3]), argc == 5 ? parse<double>(argv[4]) : 0.0};

    const YearData data = compute(year, obs);
    print_summary("sunrise", data.rise, data, year);
    print_summary("sunset", data.set, data, year);

    glfwSetErrorCallback([](int, const char* msg) { std::println(stderr, "spp_spiral: GLFW: {}", msg); });
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 8);

    const std::string place = std::format("{} at {:.4f}, {:.4f}, {:.0f} m", year, obs.lat_deg, obs.lon_deg, obs.elev_m);
    const std::array<std::pair<std::string, const std::vector<Sample>*>, 2> events{{
        {"Spp - sunrise azimuth - " + place, &data.rise},
        {"Spp - sunset azimuth - " + place, &data.set},
    }};

    constexpr int size = 760;
    std::array<Plot, 2> plots;
    for (std::size_t i = 0; i < plots.size(); ++i) {
        plots[i].win = glfwCreateWindow(size, size, events[i].first.c_str(), nullptr, nullptr);
        if (!plots[i].win) {
            glfwTerminate();
            return 1;
        }
        glfwSetWindowPos(plots[i].win, 60 + static_cast<int>(i) * (size + 20), 80);
        setup(plots[i], data, *events[i].second);
    }

    // Redraw on input only: the scene is static.
    while (true) {
        bool open = false;
        for (Plot& p : plots) {
            if (!p.win) continue;
            if (glfwWindowShouldClose(p.win)) {
                destroy(p);
                continue;
            }
            open = true;
            render(p);
        }
        if (!open) break;
        glfwWaitEvents();
    }

    glfwTerminate();
    return 0;
}
