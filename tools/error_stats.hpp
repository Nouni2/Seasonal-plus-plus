#pragma once

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace spp::check {

// Reporting only: signed residuals, with explicit rejection of NaNs. These
// statistics never participate in the astronomy calculation being checked.
struct ErrorStats {
    std::size_t count = 0;
    double sum = 0.0;
    double squares = 0.0;
    double max_abs = 0.0;
    std::string worst_case;

    void add(double error, const std::string& context) {
        if (!std::isfinite(error) || !std::isfinite(error*error))
            throw std::runtime_error("non-finite residual: " + context);
        if (count == 0 || std::abs(error) > max_abs) {
            max_abs = std::abs(error);
            worst_case = context;
        }
        ++count;
        sum += error;
        squares += error*error;
    }
    double rms() const { return count ? std::sqrt(squares/static_cast<double>(count)) : 0.0; }
    double mean() const { return count ? sum/static_cast<double>(count) : 0.0; }
    bool passes(double limit) const { return count > 0 && max_abs < limit; }

    bool report(const char* label, const char* unit, double limit) const {
        const bool pass = passes(limit);
        std::printf("  %-24s max %.9g  RMS %.9g  mean %+.9g %s\n"
                    "    limit < %.9g %s; n=%zu; %s; worst: %s\n",
                    label, max_abs, rms(), mean(), unit, limit, unit, count,
                    pass ? "PASS" : "FAIL", worst_case.c_str());
        return pass;
    }
};

}  // namespace spp::check
