# powersasa

`powersasa` is a Rust crate for solvent-accessible surface area (SASA) and molecular volume computation based on a **3D weighted power diagram** (Laguerre/Voronoi) formulation.

## Implemented algorithm

This crate implements the Power-SASA method as described in the [original paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21844):

- Atoms are modeled as spheres (typically van der Waals radius + probe radius).
- A weighted power diagram is built for the atom set
- per-atom exposed spherical contours are reconstructed from local diagram topology.
- Contour integrals are used to compute:
  - per-atom and total SASA
  - per-atom and total volume

## Provenance

The Rust code is an **LLM-assisted translation** (GPT-5.2-Codex) of the legacy C++ implementation. Originally C++ code was taken from [SIMONA](https://www.int.kit.edu/1636.php) packege. Now this code is open sourced and available on [GitHub](https://github.com/multiscale-modelling/pdsasa), so original restrictive licence **do not apply**. Original C++ code (stored in `legacy/initial`) was first heavily refactored with Codex to get rid of the pointer-based logic. The final C++ code is stored in `legacy/current` Then the code was converted to Rust.

LLM translation preserved the original logic as much as possible, so all credits to the implementation of the algorithm should still be attibuted to original authors.

## Project layout

- `src/`: Rust implementation
- `tests/`: integration tests and golden test data
- `legacy/`: C++ history and reference snapshots

## License

LGPLv3 (see `LICENSE`).
