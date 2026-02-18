# Power SASA Algorithm and Code Structure

## Scope
This document explains how this codebase computes solvent-accessible surface area (SASA) and molecular volume using a 3D weighted power diagram.

Primary implementation files:
- `power_sasa.h`
- `power_diagram.h`
- `sasatest.cpp` (example and regression check)

## High-Level Idea
The implementation models each atom as a sphere (radius already includes probe radius in `sasatest.cpp`: `radius + 0.14`).

1. Build a weighted power diagram of all atoms.
2. For each atom, use local topology from the diagram (neighbor list, cell vertices, zero-crossing points).
3. Reconstruct exposed spherical contours on that atom.
4. Integrate contour geometry to obtain:
- per-atom SASA
- per-atom volume contribution
- optional gradients (`DSasa`, `DVol`)

## Core Components

### `POWER_DIAGRAM::PowerDiagram<PDFloat, PDCoord, 3>` (`power_diagram.h`)
Provides geometric/topological primitives used by SASA.

Important nested structures:
- `cell`:
  - `position`, `r`, `r2`
  - `neighboursIds`
  - `myVerticesIds`
  - `myZeroPoints`
- `vertex`:
  - `position`, `powerValue`
  - `generatorRefs` (point/side generators)
  - `endPointIds` (graph connectivity)
- `zeroPoint`:
  - edge parameter `pos` where power crosses zero
  - `fromId`, `branch`
  - generating refs (`generatorRefs[3]`)

Important public workflows:
- `create(...)`: build initial clipped diagram.
- `recalculate(...)`: rebuild diagram for updated full coordinates/radii.
- `addMore(...)`: incremental insertion of atoms.
- `revert()`: rollback to pre-insertion snapshot.
- `FillAllMyVertices()`, `FillAllNeighbours()`, `FillAllZeroPoints()`: derive cached topology used by SASA.

### `POWERSASA::PowerSasa<PDFloat, PDCoord>` (`power_sasa.h`)
High-level SASA/volume evaluator on top of `PowerDiagram`.

Main API:
- `calc_sasa_single(iatom)`
- `calc_sasa_all()`
- `update_coords(...)`
- `add_more(...)` + `revert()`
- getters: `getSasa()`, `getVol()`, `getDSasa()`, `getDVol()`

Feature flags in constructor:
- `withSasa`
- `withDSasa`
- `withVol`
- `withDVol`

Scratch arrays (`np`, `nt`, `vx`, `br_c`, `br_p`, `ang`, `next`, `volnb`, etc.) are resized dynamically by `Resize_*` helpers.

## Per-Atom Algorithm (`calc_sasa_single`)
For atom `i` with center `pos` and radius `RAD` (`RAD2 = RAD*RAD`):

### 1. Initialize and fetch local neighborhood
- Load `atom.neighboursIds`.
- Reset outputs for this atom.
- Early exits:
  - `nnb == 0`: full isolated sphere (SASA `= 4*pi*RAD2`, volume `= 4/3*pi*RAD^3` if owned).

### 2. Validate local power state
- For all vertices in `atom.myVerticesIds`:
  - if `abs(powerValue) < tol_pow`, throw exception (unstable near-zero classification).
  - track whether atom is fully covered (`powerValue <= 0` everywhere).
- If fully covered and volume is not requested, return.

### 3. Compute neighbor-circle geometry
For each neighbor `j`:
- `rel_pos = neighbor.position - pos`
- `dist = |rel_pos|`
- classify interactions:
  - total cover: `dist <= nb_RAD - RAD` => return (atom inaccessible)
  - non-contributing: `dist >= RAD + nb_RAD` or `dist <= RAD - nb_RAD`
- otherwise, compute intersection-circle parameters:
  - `costheta[j] = 0.5 * (dist + (RAD2 - nb_RAD2)/dist) / RAD`
  - `sintheta[j] = sqrt(1 - costheta[j]^2)`
  - `e[j] = rel_pos / dist`

If no contributing neighbors remain, return full sphere values.

### 4. Register contour vertices
Two sources are collected into `vx[]` (normalized on atom sphere):

1. `myZeroPoints`: explicit zero-power crossings along diagram edges.
2. fully covered negative-power edges between local power vertices.

Each contour vertex is attached to two neighbor circles via:
- `br_c[v][0/1]` (neighbor ids)
- `br_p[v][0/1]` (position index within that neighbor circle list)

For volume integration, additional per-neighbor triangle-fan accumulators are tracked via `knot[]`, `fknot[]`, `volnb[]`.

### 5. Order vertices on each neighbor circle
For each neighbor circle with `np[j] > 0`:
- `Get_Ang(...)`: compute angular parameter `phi` around the circle.
- `Get_Next(...)`: rank/order angles with tolerance handling to build successor links.

This produces cyclic traversal links `next[j][k]`.

### 6. Traverse contours and integrate
For each unvisited contour vertex:
- determine traversal direction via sign of `(e[ic2] x e[ic1]) . vx[start]`
- walk contour cycle using `next` links and bridge tables.

During traversal:
- angular increment: `phi = ang[next] - ang[current]` (wrapped to `[0, 2*pi)`)
- spherical term:
  - `co = (e1.e2 - costheta1*costheta2)/(sintheta1*sintheta2)` (clamped to `[-1,1]`)
  - accumulate `sasa_ia += phi*costheta2 - acos(co)` plus per-contour `2*pi` offset
- derivatives (optional): update `DSasa_parts` and `DVol` using local arc/cross-product terms
- volume auxiliaries:
  - `scone = sin(theta)^2 * (phi - sin(phi))`
  - accumulate `vol2` and per-neighbor `volnb`

### 7. Add isolated cap contributions (“single circles”)
For neighbors with `np[j]==0 && nt[j]==0`:
- test cap-center visibility against all other neighbors by power comparison.
- if visible, add closed-form cap contribution:
  - area: `2*pi*(1 + costheta[j])`
  - volume/derivative terms using same `scone` form.

### 8. Finalize atom outputs
- `Sasa[i] = RAD2 * sasa_ia`
- `DSasa[i] = -sum_j DSasa_parts[i][j]`
- `Vol[i] = RAD*RAD2*sasa_ia/3 + RAD*RAD2*vol2/6 + vol3/6`
- reset temporary `visitedAs` flags on neighbors

## Numeric Stability and Error Guards
`power_sasa.h` includes multiple protections:
- tolerance functions:
  - `DRAD2()` for power/radius comparisons
  - `DANG()` for angular sorting
- `clamp_unit_interval()` before `acos`
- near-collinear handling in `Ang_About`
- explicit exceptions for ambiguous/invalid states:
  - axis too short
  - unresolved angle ties
  - odd circle crossing counts
  - invalid neighbor/vertex links
  - contour walk overflow (`kMaxCount`)

## Data and Coordinate Conventions
- `PowerDiagram` stores point positions shifted by `center` (bounding-box center).
- `calc_sasa_single` uses those shifted coordinates consistently.
- In `sasatest.cpp`, radii are already solvent-expanded (`vdW + probe`).
- Units are inherited from input (example prints `nm^2`, `nm^3`).

## Execution Path in `sasatest.cpp`
1. Read atom coordinates/radii from `protein_coords.txt`.
2. Add probe radius `0.14`.
3. Construct `PowerSasa<float, Vec3<float>>(coords, weights, withSasa=1, withDSasa=0, withVol=1, withDVol=0)`.
4. Call `calc_sasa_all()`.
5. Sum per-atom arrays from `getSasa()` and `getVol()`.
6. Compare against golden values:
- total SASA: `144.812`
- total volume: `56.747`

## Practical Notes for Future Sessions
- If coordinates/radii change for all atoms, prefer `update_coords(...)` + `calc_sasa_all()`.
- If exploring local perturbations/insertion sensitivity, use `add_more(...)`, compute, then `revert()`.
- If enabling gradients, construct with `withDSasa` and/or `withDVol`.
- For debugging topology issues, inspect per-cell:
  - `neighboursIds`
  - `myVerticesIds`
  - `myZeroPoints`

## `bond_to` Parameter (Constructor Guidance)
`PowerSasa` has a constructor that accepts `bond_to` and forwards it to `PowerDiagram::create(...)`.

Purpose:
- `bond_to` is an insertion/topology hint for power-diagram construction.
- It does not directly change SASA/volume equations; it helps the diagram code pick a good representative start for local insertion/update steps.
- Without explicit `bond_to`, the fallback constructor builds a chain (`bond_to[i] = i-1`).

How to choose `bond_to` for a real molecule:
1. Use a consistent atom order (typically input/PDB or residue order).
2. Build a spanning-tree parent map from the molecular bond graph.
3. For each atom `i > 0`, choose `bond_to[i]` as one bonded parent atom with index `< i`.
4. For multiple disconnected fragments, root each fragment and assign parents within that fragment.
5. If no bond graph is available, use the nearest already-inserted atom in space as a fallback heuristic.

Implementation constraint:
- `bond_to[i]` should usually point to an already added atom (`< i`). Out-of-range or not-yet-valid links are sanitized internally.

## Known Constraints
- Code is hardwired to 3D usage for SASA (`PowerDiagram<...,3>`).
- Several operations assume valid topology from `FillAll*` caches.
- Errors near degeneracies are handled by throwing exceptions rather than silently continuing.

## Licensing Note
Both headers include SASA license text from Karlsruhe Institute of Technology. Verify your use is compliant before redistribution or commercial usage.
