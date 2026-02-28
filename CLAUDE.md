# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```sh
cargo build
cargo test
cargo test sasatest_matches_golden_912340c   # run single test
cargo clippy
cargo check
cargo fmt
```

## Architecture

`powersasa` is a Rust library crate (no binary) that computes solvent-accessible surface area (SASA) and molecular volume using a 3D weighted power (Laguerre/Voronoi) diagram. The code is an LLM-assisted translation of the C++ reference in `legacy/`.

### Source layout (`src/`)

- **`lib.rs`** — re-exports the three public items: `PowerDiagram`, `PowerSasa`, `SasaError`.
- **`power_diagram.rs`** — low-level 3D weighted power diagram (incremental insertion). Key types:
  - `PowerDiagramParams<Scalar>` / `PowerDiagramRuntimeParams` — build configuration.
  - `PowerDiagram<Scalar>` — owns cells, vertices, zero-crossing caches, bounding cube.
  - `GeneratorRef` / `GeneratorKind` — internal handles distinguishing atom generators from boundary-side generators.
- **`power_sasa.rs`** — higher-level orchestration. `PowerSasa<Scalar>` wraps a `PowerDiagram` and adds contour-integral machinery for per-atom SASA and volume (plus optional gradients).
  - Constructor: `PowerSasa::new(coords_iter, radii_iter, with_sasa, with_dsasa, with_vol, with_dvol)`
  - Computation: `ps.calc_all() -> Result<(), SasaError>`
  - Results: `ps.per_atom_sasa()`, `ps.per_atom_vol()`, `ps.dsasa()`, `ps.dvol()`
- **`error.rs`** — `SasaError` enum (via `thiserror`), covering all geometry edge cases.

All structs are generic over `Scalar: nalgebra::RealField + Copy`, supporting both `f32` and `f64`. Tests use `f32`.

### Per-atom algorithm (`calc_sasa_single` / `calc_all`)

1. **Build power diagram** — atoms modelled as spheres (input radii already include probe radius, conventionally `vdW + 0.14`).
2. **Derive topology caches** — `FillAllNeighbours`, `FillAllMyVertices`, `FillAllZeroPoints`.
3. **For each atom**, using local cell topology:
   - compute intersection-circle geometry for each neighbor (`costheta`, `sintheta`, direction `e`).
   - collect contour vertices from zero-power crossings (`myZeroPoints`) and fully-covered edges.
   - order vertices angularly on each neighbor circle (`Get_Ang` / `Get_Next`).
   - walk contour cycles and accumulate spherical integrals for SASA and volume.
   - add closed-form cap contributions for isolated single circles.
4. **Finalize** — `Sasa[i] = RAD² * sasa_ia`; `Vol[i]` from three additive terms.

Scratch arrays (`np`, `nt`, `vx`, `br_c`, `br_p`, `ang`, `next`, `volnb`, …) are pre-allocated per-atom at construction time.

Numeric guards: `DRAD2()` / `DANG()` tolerances, `acos` clamping, `SasaError` variants for all degenerate states, contour-walk overflow at `K_MAX_COUNT = 100`.

### Tests (`tests/`)

`tests/sasatest.rs` runs golden-data integration tests. The manifest `tests/testdata/sasa_cases/golden_912340c.txt` lists input files with expected total SASA and total volume; tolerance is `1e-3`. Input files contain whitespace-separated `x y z radius` lines (comment lines start with `#`); the test adds probe radius `0.14` to each atom.

### Legacy C++ reference (`legacy/`)

Not part of the Rust build. Use as algorithmic ground truth when debugging numerical discrepancies:
- `legacy/POWER_SASA_ALGORITHM.md` — detailed step-by-step algorithm documentation.
- `legacy/initial/` — original C++ headers from commit `912340c`.
- `legacy/current/` — refactored C++ with tests and benchmark.

## Key dependencies

| Crate | Role |
|-------|------|
| `nalgebra` | `Vector3<Scalar>`, `RealField` trait |
| `thiserror` | Derive macros for `SasaError` |
