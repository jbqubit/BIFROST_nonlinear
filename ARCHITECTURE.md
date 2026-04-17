# BIFROST Architecture

This document describes the high-level architecture of the repository, with emphasis on the Julia port used for modernized polarization-transfer and DGD simulation. 

# structure
This is the project file structure. This is a high level schematic. Do not update it to reflect the location of all files. 

```text
.
├── .cursor                          [1]
├── AGENT.md                         [2]
├── ARCHITECTURE.md
├── README.md                        [3]
├── TODO.md                          [17]
├── fibers.py                        [4L]
├── test_fibers.py                   [5L]
├── Example-Fibers.ipynb             [6L]
├── brillouin.py, demo_wdm_with_raman.py, raman*.py [16L]
├── Papers                           [7]
├── julia-port                       [8]
│   ├── README.md                    [16]
│   ├── Project.toml
│   ├── Manifest.toml
│   ├── material-properties.jl       [9]
│   ├── fiber-cross-section.jl       [10]
│   ├── fiber-path.jl                [11]
│   ├── path-integral.jl             [12]
│   ├── fiber-path-plot.jl           [13]
│   ├── demo.jl                      [14]
│   └── test                         [15]
│       ├── runtests.jl
│       ├── test_fiber_path.jl
│       ├── test_material_properties.jl
│       ├── test_fiber_cross_section.jl
│       ├── test_paddle_transfer.jl
│       ├── test_dgd.jl
│       └── test_path_integral.jl
└── *.py / *.jl supporting scripts
```

- [2] Notes for tooling workflows (for agents).
- [3] Primary project overview and scientific context (for humans).
- [4L] Legacy Python implementation for birefringence simulation.
- [5L] Python regression and behavior checks.
- [6L] Interactive notebook example for exploration and demos.
- [7L] Research references and source material.
- [8] Active Julia refactor and solver architecture.
- [9] Material models and refractive-index behavior.
- [10] Cross-sectional physics and birefringence calculations.
- [11] Fiber authoring/specification layer and source assembly.
- [12] Adaptive propagation engine and DGD sensitivity solver.
- [13] 3D geometry and visualization output pipeline.
- [14] End-to-end runnable examples over composed fibers.
- [15] Julia test harness and physics-oriented unit coverage.
- [16] User guide for the julia port (for humans and agents).
- [17] 

## Getting Started
If you are a new agent working with this repository read the following.
- AGENT.md
- README.md is to be treated as user onboarding and system-level design.

Do not read the following unless motivated by a workflow.
- All the files marked [xL] as these are legacy files for the python. 

Comment
- The files marked [xL] are authoritative as regards physics. 

## Architectural Intent

- Separate **physics modeling**, **material & artifact specifications** and **numerical
solving** into distinct modules.
- Keep propagation methods stable under non-commuting generator terms through Lie-group style exponential stepping.
- Support continuous/function-valued fiber definitions instead of only fixed pre-sliced segment grids.
- Make extensibility explicit: new birefringence mechanisms plug in as typed sources with shared contracts.

## Layered Design

1. **Material layer** (`material-properties.jl`)
   - Encodes intrinsic optical material properties and optional spectral derivatives.
   - No dependency on fiber-path assembly.

2. **Cross-section layer** (`fiber-cross-section.jl`)
   - Encodes transverse step-index fiber geometry and local birefringence response laws.
   - Bridges material properties into physically meaningful local response coefficients.

3. **Fiber specification layer** (`fiber-path.jl`)
   - This layer considers 3D extrusions of the fiber cross-section specified in the 
     previous layer.
   - Defines `FiberSpec`, `BendSegment`, `TwistSegment`, and compiled source types (`BendSource`, `TwistSource`).
   - Enforces coverage and breakpoint validity over fiber domain.
   - Assembles fiber-level `K(s)` and `Kω(s)` through source contribution composition.

4. **Propagation layer** (`path-integral.jl`)
   - Solves `dJ/ds = K(s)J` with adaptive step-doubling exponential midpoint integration.
   - Solves coupled sensitivity system for `G = ∂ωJ` and derives DGD from `J` and `G`.
   - Uses breakpoint-aware interval decomposition to avoid integrating across discontinuities.

5. **Presentation layer** (`fiber-path-plot.jl`, `demo.jl`)
   - Generates visual diagnostics and runnable examples.
   - Demonstrates typical composition + propagation workflow.

## Runtime Flow
0. See demo.jl as an example of how to setup a simluation. 
1. Build a `FiberSpec` with domain, cross-section, and wavelength.
2. Author bend/twist segments (`bend!`, `twist!`) with scalar or function-valued profiles.
3. Compile to a `Fiber` (source objects + coverage + breakpoints).
4. Propagate with `propagate_fiber` for Jones output or `propagate_fiber_sensitivity` for DGD.
5. Post-process outputs for diagnostics, plots, and regression checks.

## Contracts and Invariants

- Each source must cover the full fiber domain; uncovered intervals are validation failures.
- Breakpoints are normalized and globally merged before piecewise propagation.
- Numerical tolerances (`rtol`, `atol`, step controls) are explicit API inputs, not hidden globals.
- Global phase-insensitive error metrics are used in adaptive acceptance checks for physically meaningful convergence behavior.

## Testing Strategy

- Test entrypoint: `julia-port/test/runtests.jl`.
- Emphasis areas:
  - source composition and domain validation,
  - material/cross-section calculations,
  - propagation behavior and DGD computation,
  - paddle-like path construction.
- Current status note: the committed test harness references `test_path_integral_sources.jl`, while ongoing work is transitioning toward `test_path_integral.jl`. Keep test wiring aligned as files evolve.

## Extension Guidance

- Add new birefringence mechanisms by introducing a new `AbstractBirefringenceSource` subtype plus:
  - `generator_K_contribution`,
  - `generator_Kω_contribution`,
  - coverage and breakpoint declarations.
- Keep source logic pure and local; avoid coupling source behavior to solver internals.
- Prefer adding validation tests for new source invariants before tuning solver parameters.
- Preserve separation between lossless Jones propagation assumptions and any future loss/attenuation model.

## Operational Best Practices

- Keep architecture docs schematic; avoid drifting into exhaustive file inventory.
- When implementing a new feature hold off on updating docs, tests and interfaces 
  until things have settled. 
- When a feature looks to be settled ask the user if it is time to update docs, tests, and interfaces when changing source contracts or solver semantics.
- Favor deterministic demos and tests so numerical regressions are quickly detectable.
