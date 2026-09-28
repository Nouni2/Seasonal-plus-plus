cmake_minimum_required(VERSION 3.28)
project(Spp VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

option(SPP_NATIVE "Tune for the build machine (-march=native)" OFF)
option(SPP_BUILD_BENCH "Build the benchmark executable" OFF)
option(SPP_EXACT "Reference mode: direct series evaluation, textbook formulas (30-40x slower)" OFF)

include(cmake/SppFlags.cmake)

# ------------------------------------------------------------------------------
# VSOP87D tables: data/VSOP87D.ear is compiled into the library at build time.
# ------------------------------------------------------------------------------
add_executable(vsop_gen tools/vsop_gen.cpp)
spp_apply_flags(vsop_gen)

set(SPP_GENERATED_DIR ${CMAKE_BINARY_DIR}/generated)
set(SPP_VSOP_INC ${SPP_GENERATED_DIR}/spp/vsop87d_earth_terms.inc)

add_custom_command(
    OUTPUT ${SPP_VSOP_INC}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${SPP_GENERATED_DIR}/spp
    COMMAND vsop_gen ${CMAKE_SOURCE_DIR}/data/VSOP87D.ear ${SPP_VSOP_INC}
    DEPENDS vsop_gen ${CMAKE_SOURCE_DIR}/data/VSOP87D.ear
    COMMENT "Generating VSOP87D Earth tables"
    VERBATIM
)

# ------------------------------------------------------------------------------
# Library
# ------------------------------------------------------------------------------
add_library(spp STATIC
    src/core/time.cpp
    src/core/calendar.cpp
    src/sun.cpp
    src/earth.cpp
    src/sky.cpp
    src/events.cpp
    src/earth_state.cpp
    src/model/vsop87d.cpp
    src/model/iau1980.cpp
    src/model/laskar.cpp
    src/model/meeus.cpp
    src/model/espenak_meeus.cpp
    src/model/saemundsson.cpp
    ${SPP_VSOP_INC}
)
target_include_directories(spp
    PUBLIC  ${CMAKE_SOURCE_DIR}/include
    PRIVATE ${SPP_GENERATED_DIR}
)
spp_apply_flags(spp)
if(SPP_EXACT)
    target_compile_definitions(spp PUBLIC SPP_EXACT)
endif()

# ------------------------------------------------------------------------------
# Apps
# ------------------------------------------------------------------------------
add_executable(spp_cli apps/spp/main.cpp)
# The library is C++20; executables use C++23 for std::print.
set_target_properties(spp_cli PROPERTIES OUTPUT_NAME spp CXX_STANDARD 23)
target_link_libraries(spp_cli PRIVATE spp stdc++exp)   # stdc++exp: std::print on GCC 15
spp_apply_flags(spp_cli)

# Sunrise/sunset azimuth spirals (OpenGL); built when GLFW, GLEW and GLM are available.
find_package(glfw3 CONFIG QUIET)
find_package(glew CONFIG QUIET)
find_package(glm CONFIG QUIET)
find_package(OpenGL QUIET)
if(TARGET glfw AND TARGET GLEW::glew AND TARGET glm::glm AND TARGET OpenGL::GL)
    add_executable(spp_spiral apps/spiral/main.cpp)
    set_target_properties(spp_spiral PROPERTIES CXX_STANDARD 23)
    target_link_libraries(spp_spiral PRIVATE spp glfw GLEW::glew glm::glm OpenGL::GL stdc++exp)
    spp_apply_flags(spp_spiral)
else()
    message(STATUS "spp_spiral not built: GLFW, GLEW, GLM or OpenGL not found")
endif()

if(SPP_BUILD_BENCH)
    add_executable(spp_bench bench/bench_pipeline.cpp)
    set_target_properties(spp_bench PROPERTIES CXX_STANDARD 23)
    target_link_libraries(spp_bench PRIVATE spp stdc++exp)
    spp_apply_flags(spp_bench)
endif()

# Accuracy of the Chebyshev cache against direct series evaluation
add_executable(spp_cache_check validation/cache_check.cpp)
set_target_properties(spp_cache_check PROPERTIES CXX_STANDARD 23)
target_include_directories(spp_cache_check PRIVATE ${CMAKE_SOURCE_DIR}/src ${SPP_GENERATED_DIR})
target_link_libraries(spp_cache_check PRIVATE spp stdc++exp)
spp_apply_flags(spp_cache_check)
