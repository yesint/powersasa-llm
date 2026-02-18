# PowerSASA Modernization Handoff

## Session Snapshot
- Branch: `main`
- HEAD: `dc39e02`
- Regression gate: pass (`make -C build` then `./build/sasatest`)
- Working tree: clean

## Required Test Gate (after every change)
1. `make -C build`
2. `./build/sasatest`

If `./build/sasatest` exits non-zero, treat as regression.

## What Was Completed In This Session
- Continued ID-first migration automatically across insertion/recovery/revert hot paths.
- Converted remaining insertion-front reads from pointer helper lookups to ID/data-driven access:
  - replaced `involved_front_ptr()` call sites with `involved_id_at(0)` / direct insertion cell usage.
- Removed dead helper:
  - `involved_front_ptr()`.
- Removed additional pointer-fallback read branches where IDs are already validated:
  - switched involved loops to `cell_at(id)` in:
    - revert zero cleanup
    - involved-neighbour refresh/reset loops
    - insertion recovery representative reassignment
    - connection/reserve/representative helper paths
    - nested `vertex` insertion/replace-check paths
- Kept known-safe exception path unchanged:
  - `SetInvolvedPersistingVisitedToZero()` still uses `involved_ptr_at(...)` for index `>=1` to include side generators.
- Mutation synchronization remains helper-only for cell mirrors:
  - `bondToId`, `neighboursIds`, `myVerticesIds` writes go through helper APIs.

## New Commits (latest first)
- `dc39e02` Drop redundant null check in involved visited reset
- `ea9c2e0` Use validated involved IDs directly in read paths
- `5c20e82` Remove unused involved front pointer helper
- `0fa5948` Use insertion cell directly in buildVertices recovery loop
- `1ab9a29` Use involved front IDs in insertion and vertex paths
- `66a70b6` Use insertion cell directly in numerical fallback loop
- `b59a563` Use insertion cell directly in replaced/involved seed check
- `f1ebd25` Use helper to erase cell vertices by ID in UpdateUnused
- `6a5ec09` Use involved point IDs in revert and neighbour resets
- `4485202` Use involved point IDs in insertion recovery reads

## Current Technical State
- Hot-path migration plan items are functionally complete for this tier:
  - temp-work containers are ID/ref based (`unused`/`Invalids` IDs, `ReplacedIds`, `InvolvedRefs`).
  - conversion helpers are centralized and heavily used (`cell_at`, `vertex_at`, `*_id_at`, setter helpers).
  - targeted read-path fallback cleanup was completed where IDs are guaranteed.
  - mutation writes for mirrored cell links are helper-centric.
- Remaining pointer model is still structural in algorithms (`cellPtr`/`vertexPtr` throughout geometry/topology ops), but no regression from migration work.

### Confirmed Migration Blockers (keep as-is unless re-validated)
- Immediate endpoint-ID writes in live rewiring can regress geometry:
  - signature: `PowerSasa: odd number of crossing...`
  - safe pattern remains: update pointer now, defer ID write (`kInvalidId`) and sync later.
- `SetInvolvedPersistingVisitedToZero()` must still clear via full involved refs (`involved_ptr_at(...)` for index `>=1`):
  - prior point-ID-only conversion caused segfault (`EXIT:139`), indicating non-point involved refs participate.

## Next Session Options
1. Start deeper pointer-retirement work (step 6): progressively replace pointer-heavy internal APIs with ID-based wrappers where practical.
2. Expand regression coverage beyond `sasatest` before additional structural pointer removals.
3. Keep current state as stable migration checkpoint and move to adjacent files (`power_sasa.h`) with the same method.
