# Contributing to Seasonal++

Thank you for your interest in contributing to Seasonal++.
Spp is a **precision-first** solar positional astronomy engine so all contributions must respect these rules:

### 1. Numerical Policy
* **IEEE 754 Binary64**: All computations use `double`. Fast-math and implicit FMA contraction are forbidden: GCC and Clang builds use `-fno-fast-math` and `-ffp-contract=off`, MSVC builds use `/fp:precise` (see `cmake/SppFlags.cmake`). Code generation must remain deterministic across target compilers.
* **Model Implementation**: Algorithms live under `spp::model`, cite primary literature, and contain no undocumented fudge factors.
* **Typed Time**: Time scales (UTC, UT1, TT) are distinct types (`Instant<Scale>`). Instants preserve split Julian Date precision. The C interface cannot carry C++ types: there, every instant is UTC, an `int64_t` of milliseconds since 1970-01-01T00:00Z whose name ends in `utc_ms`, converted to and from `Time` at the boundary only.

### 2. Architecture & Scope
* **Engine First**: The library holds all astronomical calculations: public headers in `include/spp/`, implementation in `src/`. Front-end apps (`apps/`, `tools/`) are thin I/O wrappers. The C interface (`include/spp/spp_c.h`, `capi/`) is one too: it converts types, validates input and calls the C++ API, and exposes no result the C++ API lacks.
* **Standalone**: Zero external runtime dependencies. Astronomical data tables are embedded as compile-time constants.

### 3. Change Verification
* **`[B]` Bit-Identical**: Refactorings and optimizations must produce bit-identical results against golden reference tests.
* **`[R]` Result Changes**: Algorithmic or model fixes must be verified against independent references (SOFA/ERFA, JPL ephemerides) and extended precision (`tools/pipeline_check.cpp`), accompanied by an entry in the wiki [Changelog](https://github.com/Nouni2/Seasonal-plus-plus/wiki/Changelog).

### 4. Integrity
* **Reproducible Claims**: Any figure in documentation or benchmarks must be directly reproducible via repository tools (`tools/`).
* **Clean Compilation**: Code must compile without warnings under `-Wall -Wextra -Wpedantic -Wshadow -Wdouble-promotion`.
