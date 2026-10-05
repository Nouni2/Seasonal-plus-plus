# The Spp Constitution

| | |
|---|---|
| **Version** | 1.0 |
| **Ratified** | 2026-09-30, with Spp 1.0.0 |
| **Applies to** | Everything in this repository: library, tools, executables, build, documentation and wiki |

## Preamble

Spp (Seasonal++) exists to compute where the Sun is, and when it rises, transits and sets, as
precisely as its theories allow and as fast as that precision permits. This document states the
principles that make it that, and that no change may break. Code, documentation and process all
answer to it.

Each article gives the **rule**, the **reason**, and how compliance is **checked**.

---

## Article I — Purpose and scope

1. Spp is a **calculation engine** for solar positional astronomy: the Sun's position, Earth
   orientation, observer corrections and solar events.
2. Spp computes **the Sun only**. Other bodies are outside its scope for the time being.
3. Spp has **no web layer**: no server, no network interface, no HTML or JSON front end. It is a library
   plus thin command-line and visual tools.
4. Time zones, calendars other than Julian/Gregorian, and presentation of results belong to the programs
   that use Spp, not to Spp.

**Reason.** A narrow scope is what allows the depth of precision the project promises.
**Check.** A change that adds a body, a service or a presentation layer requires an amendment first
(Article XIV).

## Article II — Precision first

1. **Speed never costs accuracy.** No optimization may make any result less accurate.
2. All computation is in IEEE 754 `double`. **Fast-math, and every one of its sub-flags, is forbidden.
   Floating-point contraction (implicit FMA) is forbidden.** The flags live in `cmake/SppFlags.cmake`
   and are private to the library.
3. Every shortcut (caches, fits, reformulations) is **measured against an extended-precision evaluation**
   of the same theory, and must match it to within double-precision rounding.
4. Every change is classified before it lands:
   - **[B] bit-identical**: results must not change at all;
   - **[R] rounding order changes**: the error against the extended-precision reference must be no worse
     than before (max, RMS and bias), on the full supported date range.
5. The limits of the platform are stated, not hidden: the math library (`sin`, `cos`, …) is not
   correctly rounded and differs between platforms.

**Reason.** Precision is the project's reason to exist; once lost it cannot be recovered silently.
**Check.** `spp_cache_check` must pass; [R] changes need the end-to-end reference comparison.

## Article III — Correct definitions

1. Every physical condition the engine solves for (sunrise, sunset, twilight, transit, day type) has a
   **written specification**, and the code is reviewed against it.
2. Solving an equation precisely does not make it the right equation. **Correctness is established
   against independent references** (published almanacs, SOFA/ERFA, observatory data), not only against
   the engine's own higher-precision copy.
3. A difference from an independent reference must be explained by a documented, **correct** definition;
   it is never waved through as "a different convention".

**Reason.** The 1.0.0 audit found a sunrise condition computed to a microsecond and wrong by tens of
seconds. Self-consistency checks could not see it.
**Check.** Each event definition exists in the documentation and in a test against an independent
source.

## Article IV — The library is the engine

1. Everything Spp can compute is available through its C++ API (`include/spp/`).
2. The `spp` executable and the other tools are **thin front ends**: they parse input, call the API and
   print. They contain no astronomy of their own.
3. What an executable does not print says nothing about what the engine can compute. A missing output
   is solved with a program that uses the library, never by moving computation into a tool.

**Reason.** One implementation, reachable from everywhere, cannot drift from itself.
**Check.** No file under `apps/` or `tools/` implements a model, a transform or a solver. The one
exception is validation: a checking tool (such as `spp_cache_check`) may re-evaluate a model
independently, in higher precision, to verify the engine.

## Article V — Explicit models

1. Every number comes from a **named theory**, in its own namespace under `spp::model`, named after its
   source (`vsop87d`, `iau1980`, `laskar`, `meeus`, `espenak_meeus`, `saemundsson`).
2. Theories are used **in full** unless a truncation is documented with its error.
3. No undocumented empirical corrections, tuning constants or "fudge factors".

**Reason.** Which model produced a number must always be known, so that its error can be known.
**Check.** Every constant in `src/` is traceable to a cited source or derived in the documentation.

## Article VI — Time is typed

1. **UTC, UT1 and TT are distinct types.** Mixing time scales must not compile.
2. Instants are stored as a **split Julian Date** (integer day + fraction), so precision does not depend
   on the date.
3. UTC is the scale of input and output; UT1 drives the Earth's rotation; TT drives the ephemeris.

**Reason.** A 69-second scale confusion is the easiest large error to make and the hardest to see.
**Check.** New APIs that take instants take `Instant<Scale>` or `Time`, never a bare `double`.

## Article VII — Standalone

1. The library depends on **nothing but the C++ standard library**. Optional visual tools may use
   graphics libraries; the library never does.
2. The library needs **no data files at run time**: source data is compiled in, exactly (hex-float
   literals).
3. Spp is described on its own terms. Its documentation, comments and file names do not refer to other
   implementations it may resemble or have been derived from.

**Reason.** A self-contained engine is reproducible, portable and easy to embed.
**Check.** `target_link_libraries(spp …)` lists no external library; no file I/O in `src/`.

## Article VIII — Honest claims

1. **Every number published** in the README, the wiki or a release note is **reproducible** with a tool
   or an example in the repository, and is re-checked when the code it describes changes.
2. **Known defects are published** with their impact and a workaround, until they are fixed.
3. **Performance figures state their context**: machine, compiler, build type, and warm or cold cache.
4. Measurements are reported as measured; estimates are labelled as estimates.

**Reason.** A precision library is only as trustworthy as its statements about itself.
**Check.** A release is not published while a published figure is known to be wrong.

## Article IX — Performance, within Article II

1. Performance is pursued **only within the bounds of Article II**.
2. Optimizations are **measured, not assumed**: every performance change is benchmarked before and
   after, on the same machine.
3. The published figures cover both the warm path (consecutive instants) and the cold path (a new
   16-day segment).

**Reason.** Speed matters for batch work, but a fast wrong answer is worthless.
**Check.** `spp_bench` results accompany every [R] or performance change.

## Article X — Data provenance

1. Every data table (the VSOP87D series, the nutation terms, the leap-second table, the observed ΔT
   values) **cites its exact source**.
2. Tables that go stale (leap seconds, observed ΔT) state the date up to which they are valid, and the
   project documents how to refresh them.
3. Values are copied from the **primary source**, never from a summary or a secondary copy.

**Reason.** The 1.0.0 audit found an observed ΔT table wrong by up to 0.075 s because it had never been
checked against its source.
**Check.** Each table's comment names its source and date; the refresh procedure is in the Development
documentation.

## Article XI — Stable interface

1. Spp follows **Semantic Versioning**. The public API is the headers under `include/spp/`.
2. A release that breaks the public API raises the major version.
3. A fix that **changes results** without changing the API is released as a patch, with a **prominent
   changelog entry** stating what changes and by how much.

**Reason.** Programs that depend on Spp must be able to upgrade knowingly.
**Check.** Every release has a changelog entry; result-changing fixes are marked as such.

## Article XII — Documentation

1. The **README** is the complete reference for a release. The **wiki** is the extended documentation:
   tutorials, cookbook, API reference, theory, known issues and roadmap.
2. Every code example in the documentation **compiles and runs** against the release it documents, and
   its printed output is the real output.
3. Command examples are written for **PowerShell**: `.\spp day 2026-06-21 48.8566 2.3522`.

**Reason.** Documentation that does not run misleads exactly the readers who trust it.
**Check.** Examples are rebuilt and re-run for every release.

## Article XIII — Process

1. Work is tracked as **issues**, one per finding or task, grouped into **milestones** that follow the
   roadmap, and labelled by severity and kind.
2. Releases are **tagged**, built from the tag, packaged, and published with release notes that include
   known issues.
3. A release is not published while an issue labelled `S1-critical` is open against it, unless the release
   notes state it prominently with its workaround.

**Reason.** A visible plan and a repeatable release are what make the other articles enforceable.
**Check.** Every release has a tag, a package, notes and a milestone.

---

## Article XIV — Amendments

1. This constitution changes only through a commit that **edits this file**, states the **reason** in its
   message, and **raises the version** (minor for a clarification, major for a changed principle).
2. An amendment takes effect for the next change after it is committed; it does not retroactively
   excuse earlier violations.

## Article XV — Supremacy

1. Where a guideline, comment, habit or document conflicts with this constitution, **the constitution
   prevails**.
2. Existing code that violates it is a defect to be tracked as an issue and fixed, not a precedent.

---

*Known deviations at ratification (Spp 1.0.0), tracked as issues:* the sunrise/sunset definition
(Article III, issue BUG-09), the unchecked observed ΔT table (Article X, BUG-02), published claims
without a reproducing tool (Article VIII, CLM-10), and the absence of an independent end-to-end
reference (Article III, TST-01).
