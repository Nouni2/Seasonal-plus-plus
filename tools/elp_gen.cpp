// ELP/MPP02 data -> exact C++ coefficient literals. No network or Python is
// required to build Spp. Source coefficients are vendored under data/elp_mpp02.
#include <array>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

namespace {

struct Table { const char* file; const char* name; int count; bool main; };
constexpr Table tables[] = {
    {"elp_main.long", "main_long", 1023, true},
    {"elp_main.lat", "main_lat", 918, true},
    {"elp_main.dist", "main_dist", 704, true},
    {"elp_pert.longT0", "pert_long0", 11314, false},
    {"elp_pert.longT1", "pert_long1", 1199, false},
    {"elp_pert.longT2", "pert_long2", 219, false},
    {"elp_pert.longT3", "pert_long3", 2, false},
    {"elp_pert.latT0", "pert_lat0", 6462, false},
    {"elp_pert.latT1", "pert_lat1", 516, false},
    {"elp_pert.latT2", "pert_lat2", 52, false},
    {"elp_pert.distT0", "pert_dist0", 12115, false},
    {"elp_pert.distT1", "pert_dist1", 1165, false},
    {"elp_pert.distT2", "pert_dist2", 210, false},
    {"elp_pert.distT3", "pert_dist3", 2, false}
};

template<class T>
T number(std::istringstream& row) {
    std::string token;
    if (!(row >> token)) throw std::runtime_error("missing coefficient");
    T value{};
    const auto [end, error] = std::from_chars(token.data(), token.data() + token.size(), value);
    if (error != std::errc{} || end != token.data() + token.size())
        throw std::runtime_error("invalid coefficient: " + token);
    if constexpr (std::is_floating_point_v<T>)
        if (!std::isfinite(value)) throw std::runtime_error("non-finite coefficient");
    return value;
}

void end_of_row(std::istringstream& row) {
    std::string extra;
    if (row >> extra) throw std::runtime_error("unexpected extra coefficient: " + extra);
}

std::string literal(double value) {
    std::array<char, 64> buffer{};
    const bool negative = std::signbit(value);
    const auto [end, error] = std::to_chars(buffer.data(), buffer.data() + buffer.size(),
        negative ? -value : value, std::chars_format::hex);
    if (error != std::errc{}) throw std::runtime_error("cannot format coefficient");
    return std::string(negative ? "-0x" : "0x") + std::string(buffer.data(), end);
}

void emit(const std::filesystem::path& directory, const Table& table, std::ostream& output) {
    std::ifstream input(directory / table.file);
    if (!input) throw std::runtime_error("cannot open input");
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("missing term count");
    std::istringstream header(line);
    if (number<int>(header) != table.count) throw std::runtime_error("wrong term count");
    end_of_row(header);
    output << "inline constexpr " << (table.main ? "MainTerm " : "PerturbationTerm ")
           << table.name << "[] = {\n";
    for (int i = 0; i < table.count; ++i) {
        if (!std::getline(input, line)) throw std::runtime_error("truncated table");
        std::istringstream row(line);
        output << "    {{";
        for (int j = 0; j < (table.main ? 4 : 13); ++j) {
            if (j) output << ", ";
            output << number<int>(row);
        }
        output << (table.main ? "}, {" : "}, ");
        for (int j = 0; j < (table.main ? 7 : 2); ++j) {
            if (j) output << ", ";
            output << literal(number<double>(row));
        }
        output << (table.main ? "}},\n" : "},\n");
        end_of_row(row);
    }
    while (std::getline(input, line)) {
        std::istringstream row(line);
        end_of_row(row);
    }
    if (input.bad()) throw std::runtime_error("input read failure");
    output << "};\n\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::runtime_error("usage: elp_gen <data-directory> <output.inc>");
        // Validate and generate in memory before opening the output file.
        std::ostringstream output;
        output << "// Generated from data/elp_mpp02 by elp_gen. Do not edit.\n"
                  "namespace spp::model::elp_mpp02::detail::data {\n\n";
        for (const Table& table : tables) {
            try { emit(argv[1], table, output); }
            catch (const std::exception& error) {
                throw std::runtime_error(std::string(table.file) + ": " + error.what());
            }
        }
        output << "}  // namespace spp::model::elp_mpp02::detail::data\n";
        std::ofstream file(argv[2], std::ios::binary);
        file << output.str();
        file.close();
        if (!file) throw std::runtime_error("cannot write output");
        std::cout << "elp_gen: 35901 terms written to " << argv[2] << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "elp_gen: " << error.what() << '\n';
        return 1;
    }
}
