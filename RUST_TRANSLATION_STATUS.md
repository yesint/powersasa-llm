# Rust PowerSASA Translation Status

This document is a session handoff for the ongoing literal C++ -> Rust translation of:
- `power_diagram.h`
- `power_sasa.h`

Target Rust files:
- `rust/src/power_diagram.rs`
- `rust/src/power_sasa.rs`
- regression runner: `rust/src/bin/sasatest.rs`

Date context: current working session up to commit `a60b140`.

## 1. Goal and Constraints

Primary goal is a literal translation with near 1:1 logic/dataflow correspondence to C++.

Key constraints followed:
- `nalgebra::Vector3<Scalar>` is used.
- iterator-style APIs use `impl Iterator<Item = ...>`.
- stage-by-stage commits were made continuously.
- avoid long-lived temporary stubs; progressively replace with literal logic.

## 2. Current Regression State

Command used:
- `cargo run --bin sasatest` from `rust/`

Current totals (latest):
- Protein (`../protein_coords.txt`)
  - SASA: `800.452942` expected `144.811996`
  - Volume: `147.065170` expected `56.747002`
- Subset cases 01..05 remain relatively close but still failing tolerance.

Interpretation:
- Global topology/contour behavior is still non-equivalent, despite large portions of insertion logic now ported.
- Not a simple “no neighbors” failure mode (neighbor density was checked during debugging and is not the dominant issue).

## 3. What Has Been Ported (High-Level)

### 3.1 `PowerDiagram` major blocks now present in Rust

Implemented in `rust/src/power_diagram.rs`:
- Core data types and helper accessors:
  - `Vertex::is_corner`, `is_connected`, `resolved_*`, `powerdiff3d`, endpoint slot helper.
- Cube build path:
  - corner creation, side generators, endpoint wiring.
  - slot permutation (`swap_vertex_link_slots`) matching C++ `buildCube` ordering step.
- `insert_first` behavior:
  - corner power init from first site.
  - point-0 corner ownership seeding into `my_vertices_ids`.
- Insertion preparation and replacement traversal:
  - `prepare_insertion` with stabilization loop.
  - `find_replaced_vertex` (1-hop/2-hop/3-hop/full scan style).
  - `finite_replaced`, `replace_check`, corner/finite replaced propagation.
- Insertion materialization pieces:
  - `try_to_build_vertex_on_edge`.
  - `init_new_vertex_from_replaced`.
  - `create_finite_vertices_from_replaced`.
  - connection via planes (`connect_new_finites_among_themselves_3d`, `register_vertex_for_connection_3d`).
  - `update_unused`, representative reassignment.
- `build_vertices` retry/rollback loop:
  - includes rollback/dampening path and identical-point-style handling logic scaffold.
  - now takes `(n_points, from)` semantics and callers updated.
- End-of-build compaction:
  - unused-slot compaction with endpoint rewiring (`compact_unused_vertices`, `move_vertex_network_update_only`).
- Cache fills:
  - `fill_all_my_vertices*`, `fill_all_neighbours*`, `fill_all_zero_points*` are present.

### 3.2 `PowerSasa` alignment items

Implemented/fixed in `rust/src/power_sasa.rs`:
- `calc_sasa_single` main contour traversal pipeline is present.
- `nnb == 0` owner gate aligned with C++ (full sphere only if first vertex owner matches atom).

## 4. What Was Explicitly Removed as Non-Literal

Removed during translation to reduce heuristic drift:
- post-build brute-force corner owner recomputation (not in C++ insertion path semantics).
- overlap-based neighbor fallback in `fill_all_neighbours` (non-literal emergency path).
- `prepare_insertion` mutations of `bond_to_id` and non-C++ fallback return logic.

## 5. Recent Commit Trail (Most Relevant)

Recent sequence around current state:
- `a60b140` step53: strict, failure-propagating finite connection registration.
- `4c1cd09` step52: corrected first-neighbor control flow in `find_replaced_vertex`.
- `cea47cb` step51: fixed signed rollback arithmetic for vertex count delta.
- `38939e2` step50: `build_vertices(n_points, from)` semantics.
- `f392911` step49: C++ ownership gate for `nnb == 0` in `PowerSasa`.
- `740ff93` step48: `prepare_insertion` returns stabilized hint id directly.
- `eefe5b5` step47: removed non-literal neighbour overlap fallback.
- `7a5656f` step46: removed non-literal corner-owner recomputation.
- `bfd529d` step45: removed non-literal bond-to mutations.
- `95369ee` step44: insert-first corner ownership seeding.
- `15c9f84` step43: unused compaction + endpoint rewiring.
- `6aab2e8` step42: retry/rollback insertion loop + failure semantics.

## 6. Remaining Gaps (Actionable)

The translation is substantial but still behavior-divergent. Highest-priority remaining equivalence gaps:

1. `build_vertices` rollback state transitions are still approximate in corner cases.
- File: `rust/src/power_diagram.rs`
- Function: `build_vertices`
- Risk: small deviation in rollback bookkeeping yields large contour errors.

2. `do_insertion` lifecycle is not yet fully equivalent to C++ invariants.
- File: `rust/src/power_diagram.rs`
- Functions: `do_insertion`, `create_finite_vertices_from_replaced`, `connect_new_finites_among_themselves_3d`, `update_unused`
- Risk: malformed finite graph produces inflated exposed surface/volume.

3. `add_more` / `revert` remain simplified compared to full C++ branch logic.
- File: `rust/src/power_diagram.rs`
- Functions: `add_more`, `revert`
- Missing detail: full reallocation/rebuild branch parity and mirror restoration behavior.

4. Full invariant/mirror maintenance layer from C++ is not fully ported as active checks.
- C++ equivalents include `sync_*` / `validate_*` helpers.
- Rust has partial helpers but no complete invariant framework.

5. Dead fallback/helper remnants indicate non-final control flow.
- `find_replaced_vertex_corner`, `n_virtual_generators`, `power_err_scaled_epsilon` currently unused.
- These should either be integrated exactly where C++ uses equivalent paths or removed after parity is proven.

## 7. Recommended Next Session Start Plan

1. Lock down `build_vertices` + `do_insertion` parity first.
- Work from C++ blocks in this order:
  - `buildVertices` retry/rollback internals
  - `CreateFiniteVerticesFromReplaced`
  - `ConnectNewFinitesAmongThemselves3D`
  - `UpdateUnused`
  - `AssignRepresentativeVerticesToCells`

2. Add temporary invariant assertions (Rust-only debug checks) after each insertion.
- Validate endpoint symmetry and generator slot consistency.
- Validate each `my_vertices_ids` entry points to a connected vertex or expected corner case.

3. Re-run `sasatest` after each large parity block, keep committing in stages.

4. Only after insertion topology converges, revisit any remaining `PowerSasa` divergences.

## 8. Fast Resume Commands

From repo root:

```bash
cd rust
cargo check
cargo run --bin sasatest
```

Commit history focus:

```bash
git log --oneline -30
```

## 9. Notes for Next Agent

- Do not reintroduce heuristic fallbacks that were already removed.
- Treat C++ insertion topology as the ground truth and prioritize exact state transitions over local numeric tuning.
- Large protein error magnitude indicates topology-level mismatch, not a tiny numeric epsilon issue.
