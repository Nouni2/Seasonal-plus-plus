# Contributing to Seasonal++

Thank you for your interest in contributing to Seasonal++.
Spp is a **precision-first** solar positional astronomy engine so all contributions must respect these rules:

### 1. Numerical Policy
* **IEEE 754 Binary64**: All computations use `double`. Fast-math (`-ffast-math`) and implicit contraction (`-ffp-contract=off`) are forbidden. Code generation must remain deterministic across target compilers.
* **Model Implementation**: Algorithms live under `spp::model`, cite primary literature, and contain no undocumented fudge factors.
* **Typed Time**: Time scales (UTC, UT1, TT) are distinct types (`Instant<Scale>`). Instants preserve split Julian Date precision.

### 2. Architecture & Scope
* **Engine First**: All astronomical calculations reside in `include/spp/`. Front-end apps (`apps/`, `tools/`) are thin I/O wrappers.
* **Standalone**: Zero external runtime dependencies. Astronomical data tables are embedded as compile-time constants.

### 3. Change Verification
* **`[B]` Bit-Identical**: Refactorings and optimizations must produce bit-identical results against golden reference tests.
* **`[R]` Result Changes**: Algorithmic or model fixes must be verified against independent references (SOFA/ERFA, JPL ephemerides) and extended precision (`tools/pipeline_check.cpp`), accompanied by a release changelog entry.

### 4. Integrity
* **Reproducible Claims**: Any figure in documentation or benchmarks must be directly reproducible via repository tools (`tools/`).
* **Clean Compilation**: Code must compile without warnings under `-Wall -Wextra -Wpedantic -Wshadow -Wdouble-promotion`.
