# Project Structure

This is a high-level schematic. Do not update it to reflect every file.

```text
.
├── .cursor                          [1]
├── AGENTS.md                        [2]
├── ARCHITECTURE.md
├── README.md                        [3]
├── TODO.md                          [17]
├── fibers.py                        [4L]
├── test_fibers.py                   [5L]
├── Example-Fibers.ipynb             [6L]
├── brillouin.py, demo_wdm_with_raman.py, raman*.py [16L]
├── Papers                           [7L]
├── julia-port                       [8]
│   ├── README.md                    [18]
│   ├── Project.toml
│   ├── Manifest.toml
│   ├── material-properties.jl       [9]
│   ├── path-geometry.jl             [10]
│   ├── path-integral.jl             [11]
│   ├── fiber-cross-section.jl       [12]
│   ├── fiber-path.jl                [13]
│   ├── fiber-path-plot.jl           [14]
│   ├── demo*.jl                     [15]
│   └── test                         [19]
├── output                           [20]
└── *.py / *.jl supporting scripts
```

- [2] Notes for tooling workflows.
- [3] Primary project overview and scientific context.
- [4L] Legacy Python implementation for birefringence simulation.
- [5L] Python regression and behavior checks.
- [6L] Interactive notebook example for exploration and demos.
- [7L] Research references and source material.
- [8] Active Julia refactor and solver architecture.
- [9] Standalone material models and refractive-index behavior.
- [10] Standalone path construction and differential geometry.
- [11] Generic adaptive propagation for callable Jones generators.
- [12] Cross-sectional fiber optics and local birefringence responses.
- [13] Path-backed fiber assembly and generator construction.
- [14] 3D geometry and propagation visualization pipeline.
- [15] Runnable Julia demos; visual demos write HTML files to `output/`.
- [16L] Legacy Raman and Brillouin Python scripts.
- [17] TODO list for humans. Starting TODO items requires user authorization.
- [18] User guide for the Julia port.
- [19] Julia tests.
- [20] Output of demo methods and generated visual artifacts.

Files marked `[L]` are legacy files for the old Python implementation. Do not
read them unless a specific workflow requires it. They are authoritative for
legacy behavior and must not be modified without explicit user authorization.

## Architectural Intent

- Separate material physics, path geometry, fiber assembly, and numerical
  propagation.
- Keep the core propagation API usable with any callable `K(s)` and `Kω(s)`.
- Support continuous/function-valued geometry and twist rather than only fixed
  pre-sliced segment grids.
- Keep lossless Jones propagation isolated from any future gain/loss model.
- Preserve MCM compatibility on uncertainty-carrying code paths.

## Standalone Building Blocks

These files are intentionally useful on their own:

| File | Standalone role |
| --- | --- |
| `material-properties.jl` | Material constants and spectra; no path or fiber geometry. |
| `path-geometry.jl` | Three-dimensional path construction and geometric queries; no optics. |
| `path-integral.jl` | Adaptive propagation for callable `K(s)` and `Kω(s)` generators. |

The fiber-specific layers combine those pieces:

| File | How it extends the standalone pieces |
| --- | --- |
| `fiber-cross-section.jl` | Adds step-index fiber optics and birefringence responses. |
| `fiber-path.jl` | Binds path geometry to a cross section and assembles bend/twist `K` and `Kω`. |

## Layered Design

0. **Geometry layer** (`path-geometry.jl`, `path-geometry-connector.jl`,
   `path-geometry-plot.jl`)

   - Builds and queries three-dimensional paths.
   - Provides straight, bend, catenary, helix, `JumpBy`, and `JumpTo` authoring.
   - Resolves material twist metadata into path-coordinate twist runs.
   - Resolves `JumpBy` and `JumpTo` into G2 quintic connectors.

1. **Material layer** (`material-properties.jl`)

   - Encodes intrinsic optical material properties.
   - Provides spectral responses and derivatives needed by DGD calculations.

2. **Cross-section layer** (`fiber-cross-section.jl`)

   - Encodes transverse step-index fiber geometry.
   - Converts material properties into guided-index, dispersion, nonlinearity,
     and local birefringence response coefficients.

3. **Fiber assembly layer** (`fiber-path.jl`, `fiber-path-meta.jl`,
   `fiber-path-modify.jl`)

   - Binds a built `PathSpecCached` to a `FiberCrossSection` and `T_ref_K`.
   - Keeps operating wavelength as a per-query argument rather than `Fiber`
     state.
   - Assembles fiber-level bend and twist generators `K(s)` and `Kω(s)`.
   - Interprets per-segment metadata such as `Nickname`, `MCMadd`, and
     `MCMmul`.
   - Applies meta-driven path perturbations and thermal length scaling through
     `modify(fiber)`.

4. **Propagation layer** (`path-integral.jl`)

   - Solves `dJ/ds = K(s)J` with adaptive step-doubling exponential midpoint
     integration.
   - Solves the coupled sensitivity system for `G = ∂ωJ`.
   - Derives DGD from `J` and `G`.
   - Uses breakpoint-aware interval decomposition to avoid integrating across
     discontinuities.
   - Uses phase-insensitive error metrics for Jones propagation.

5. **Presentation layer** (`fiber-path-plot.jl`, `demo*.jl`)

   - Generates visual diagnostics and runnable examples.
   - Keeps visual demos as human-inspected outputs rather than reusable library
     code.

## Runtime Flow

0. See `julia-port/demo-smallest.jl` for the smallest runnable example.
1. Build a `PathSpecBuilder` with path primitives and optional metadata.
2. Compile it with `build(...)`, producing a `PathSpecCached`.
3. Bind it into `Fiber(path; cross_section, T_ref_K)`.
4. Propagate with `propagate_fiber(fiber; λ_m=...)` for Jones output.
5. Use `propagate_fiber_sensitivity(fiber; λ_m=...)` when DGD is needed.
6. Post-process outputs for diagnostics, plots, demos, and regression checks.

## Contracts and Invariants

- Path breakpoints are normalized and globally merged before piecewise
  propagation.
- The propagator must not step across path segment or twist-run boundaries.
- Numerical tolerances (`rtol`, `atol`, step controls) are explicit API inputs,
  not hidden globals.
- Global phase-insensitive error metrics are used in adaptive acceptance checks.
- `path-integral.jl` assumes lossless Jones propagation.
- MCM-compatible paths must avoid scalar coercions and particle-dependent
  conditionals on uncertainty-carrying values.

## Testing Strategy

- Test entrypoint: `julia --project=. julia-port/test/runtests.jl`.
- Emphasis areas:
  - path geometry and connector invariants,
  - material and cross-section calculations,
  - fiber assembly and breakpoint behavior,
  - propagation behavior and DGD computation,
  - MCM compatibility,
  - paddle-like path construction.
- Tests should follow the taxonomy in `AGENTS.md`: `T-PHYSICS`,
  `T-VALIDATION`, `T-SIM-REGRESSION`, `T-GUARDRAIL`, and visual demos.
- Validation against published fiber data and direct legacy `fiber.py`
  comparisons is still limited.

## Extension Guidance

- Add new path shapes by implementing the `AbstractPathSegment` local geometry
  interface in `path-geometry.jl`.
- Add new per-segment annotations by extending the `AbstractMeta` vocabulary and
  keeping interpretation in the consuming layer.
- Add new fiber-level birefringence mechanisms by extending generator assembly
  in `fiber-path.jl` and adding guardrail tests first.
- Keep solver changes in `path-integral.jl` deliberate; step controller, error
  metric, and exponential formulas are core numerical contracts.
- Preserve separation between lossless Jones propagation and any future
  attenuation, gain, or polarization-dependent loss model.

## Operational Best Practices

- Keep architecture docs schematic; avoid drifting into exhaustive file
  inventory.
- When a feature is still moving, delay broad docs/interface updates until the
  design has settled.
- When a feature settles, update docs and tests together.
- Favor deterministic demos and tests so numerical regressions are quickly
  detectable.
