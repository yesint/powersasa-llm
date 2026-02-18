# Legacy C++ Code

This folder keeps the C++ implementation history separate from the Rust crate at repository root.

## Layout

- `initial/`
  - Snapshot of the original two headers from commit `912340cd059dd000862d85834d8cef198cb351cd`:
    - `power_diagram.h`
    - `power_sasa.h`
- `current/`
  - The latest C++ codebase, tests, and test data after LLM-assisted refactoring and parity work.

## What changed from `initial` to `current`

The `current` code is not a new algorithm. It is a structural refactor of the same Power-SASA method documented in `POWER_SASA_ALGORITHM.md`.

Main change themes (from commit history and `SESSION_HANDOFF.md`):

- Modernization and cleanup:
  - macro-to-typed constants and helpers
  - stronger type usage (`std::array`, `size_t`, aliases)
  - removal of dead code and obsolete preprocessor branches
  - safer ownership with `std::unique_ptr`
- ID-first data model migration:
  - introduced mirrored IDs for cells/vertices/neighbours/generator refs
  - moved hot-path reads from pointer-first to validated ID access helpers
  - centralized mutation helpers to keep pointer/ID mirrors synchronized
- Numerical and robustness hardening:
  - centralized tolerance/guard logic
  - stricter range checks and invariant checks in insertion/rebuild paths
  - preserved known-safe pointer fallback paths where full ID-only conversion caused regressions
- Test and validation improvements:
  - modernized regression harness (`sasatest`)
  - maintained golden-data checks and benchmark/testing utilities

In short, `current` aims to preserve behavior while making the code easier to reason about and safer to evolve.
