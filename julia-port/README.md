# Julia port: path-integral birefringence propagation

J. Britton, updated 2026-05-04

This directory contains the current Julia prototype for calculating the
end-to-end polarization transfer function (PTF) and differential group delay
(DGD) of a single-mode optical fiber link. It is a refactoring of the
`fiber.py` polarization model from the BIFROST paper into a path-backed Julia
implementation.

The major architectural change is that the optical fiber is not represented as
a pre-sliced list of Jones matrices. Instead, the code builds a continuous
centerline path, binds it to a transverse fiber cross section, and integrates a
local Jones generator:

```math
\frac{dJ}{ds}=K(s,\omega)J,\qquad J(s_0)=I.
```

Propagation uses an adaptive exponential-midpoint, Lie-group style integrator.
The adaptive controller never steps across path breakpoints, and its error
metric is insensitive to physically irrelevant global Jones phase.

## Quick Start

From the repository root:

```bash
cd julia-port
julia --project=. test/runtests.jl
```

The demo files write standalone HTML artifacts under `../output/`:

```bash
julia --project=. demo-smallest.jl
julia --project=. demo1.jl
julia --project=. demo2.jl
julia --project=. demo3mcm.jl
julia --project=. demo3benchmark.jl
```

The code is currently organized as include-based Julia scripts rather than a
packaged module. Most tests and demos include only the files they need.

## Building Blocks

These files are intentionally useful on their own:

| File | Standalone role |
| --- | --- |
| `material-properties.jl` | Material constants and spectra; no path or fiber geometry. |
| `path-geometry.jl` | Three-dimensional path construction and geometric queries, with no optics. |
| `path-integral.jl` | Generic adaptive propagation for callable `K(s)` and `Kω(s)` generators. |

The fiber layers combine and specialize those pieces:

| File | How it extends the standalone pieces |
| --- | --- |
| `fiber-cross-section.jl` | Adds step-index fiber optics and birefringence responses. |
| `fiber-path.jl` | Binds path geometry to a cross section and assembles bend/twist `K` and `Kω`. |

## Current Model

High-level authoring is path based:

1. Build geometry with `PathSpecBuilder`.
2. Freeze and place it with `build(...)`, producing a `PathSpecCached`.
3. Bind that path to a `FiberCrossSection` with `Fiber(path; cross_section,
   T_ref_K)`.
4. Propagate at a requested wavelength with `propagate_fiber(fiber; λ_m=...)`.

The current path primitives include:

- `StraightSegment`
- `BendSegment`
- `CatenarySegment`
- `HelixSegment`
- `JumpBy`
- `JumpTo`

`JumpBy` and `JumpTo` are authoring conveniences. At build time they are
resolved into a G2 quintic connector implemented in
`path-geometry-connector.jl`.

Material twist is attached as per-segment metadata using `Twist`. A twist run
starts on the segment carrying the annotation and continues until the next
twist annotation or the end of the path. Twist rates may be constant or
callable functions of run-local arc length.

The optical `Fiber` stores:

- the built `PathSpecCached`,
- the `FiberCrossSection`,
- a reference temperature `T_ref_K`,
- the fiber domain `[s_start, s_end]`.

The operating wavelength is not stored on `Fiber`. It is supplied per query to
`generator_K(fiber, λ_m)`, `generator_Kω(fiber, λ_m)`, `propagate_fiber`, and
`propagate_fiber_sensitivity`.

## Example

The smallest runnable example is kept in `demo-smallest.jl`:

```bash
julia --project=. demo-smallest.jl
```

```julia
include("material-properties.jl")
include("fiber-cross-section.jl")
include("path-geometry.jl")
include("path-integral.jl")

xs = FiberCrossSection(
    GermaniaSilicaGlass(0.036),
    GermaniaSilicaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like",
)

spec = PathSpecBuilder()
straight!(spec; length = 0.5, meta = [Nickname("lead-in")])
bend!(spec; radius = 0.05, angle = pi / 2, meta = [Nickname("90 deg bend")])
straight!(spec; length = 0.5, meta = [Nickname("lead-out")])

path = build(spec)
fiber = Fiber(path; cross_section = xs, T_ref_K = 297.15)

J, stats = propagate_fiber(fiber; λ_m = 1550e-9, rtol = 1e-9, verbose = false)
```

For DGD:

```julia
J, G, stats = propagate_fiber_sensitivity(
    fiber;
    λ_m = 1550e-9,
    rtol = 1e-9,
    verbose = false,
)

dgd = output_dgd(J, G)
```

For Monte Carlo Measurements (MCM) paths, prefer `output_dgd_2x2(J, G)` over
`output_dgd(J, G)` because it avoids `eigvals`.

## Generator Formulation

The original BIFROST-style sliced approach calculates

```math
J_{\mathrm{total}}=\prod_i J_i,
```

with matrix order matching the order light encounters along the fiber. That
approach is simple, but it is difficult to attach a meaningful error bound when
linear birefringence, twist, and other non-commuting terms vary along the path.

The Julia implementation instead assembles a local generator:

```math
K(s,\omega)=K_{\mathrm{bend}}(s,\omega)+K_{\mathrm{twist}}(s,\omega).
```

The bending contribution comes from path curvature. For a local bend radius
`R(s)`, the implemented perturbation uses the bending birefringence response
from `fiber-cross-section.jl`; in the simplest stress model the magnitude
scales like `1/R(s)^2`.

The twist contribution uses the total frame twist rate:

```math
\tau_{\mathrm{path}}(s)=\tau_{\mathrm{geom}}(s)+\tau_{\mathrm{material}}(s).
```

Here `geometric_torsion(path, s)` comes from the centerline, while
`material_twist(path, s)` comes from resolved `Twist` metadata.

The same decomposition exists for the frequency derivative:

```math
K_\omega(s,\omega)
=K_{\mathrm{bend},\omega}(s,\omega)+K_{\mathrm{twist},\omega}(s,\omega).
```

That keeps ordinary Jones propagation and DGD sensitivity propagation aligned:
both use the same `Fiber`, wavelength, breakpoint partition, and adaptive
integration strategy.

## Propagation

The exponential midpoint step is

```math
J_{n+1}=\exp\!\left(hK(s_n+h/2)\right)J_n.
```

It is useful here because the solution of a constant-coefficient matrix ODE is
exactly an exponential. Step by step, the method preserves the multiplicative
structure of Jones propagation. Under the lossless assumption, after removing
common phase, the Jones matrices live in `SU(2)`.

The adaptive controller uses step doubling:

- take one full step of size `h`,
- take two half steps of size `h/2`,
- compare the two results using `phase_insensitive_error`,
- accept or reject the step,
- update `h` with a cubic-root controller because the estimate scales as
  `O(h^3)`.

Path breakpoints come from the built path:

```julia
fiber_breakpoints(fiber) = breakpoints(fiber.path)
```

Those breakpoints include path segment boundaries and resolved twist-run
boundaries. `propagate_fiber` calls `propagate_piecewise`, which integrates
independently over each smooth interval.

The 2x2 Jones exponential uses a closed form based on Cayley-Hamilton. The
implementation also factors out small numerical trace drift before applying the
traceless formula:

```math
\exp(A)=\exp(\operatorname{tr}(A)/2)
\left[\cosh(\mu)I+\operatorname{sinhc}(\mu)\tilde A\right],
\qquad \mu^2=-\det(\tilde A).
```

Here `sinhc(mu) = sinh(mu) / mu`, with a Taylor branch near zero.

## DGD Sensitivity Propagation

The finite-difference DGD estimate used by the legacy implementation has the
form

```math
\partial_\omega J
\approx \frac{J(\omega+\Delta\omega)-J(\omega)}{\Delta\omega}.
```

The Julia propagator instead integrates the sensitivity matrix
`G = partial_omega J` directly:

```math
\frac{dJ}{ds}=KJ,\qquad J(s_0)=I,
```

```math
\frac{dG}{ds}=K_\omega J+KG,\qquad G(s_0)=0.
```

At the output, the PMD generator is

```math
H_{\mathrm{PMD}}=-iJ^{-1}G.
```

`output_dgd(J, G)` returns the eigenvalue spread of that generator. For 2x2
MCM-valued matrices, `output_dgd_2x2(J, G)` computes the same spread using a
closed-form Hermitian 2x2 formula, avoiding `LinearAlgebra.eigvals`.

The coupled sensitivity step is implemented using a closed-form Frechet
derivative of the 2x2 exponential:

```math
\exp\!\left(h\begin{bmatrix}K & K_\omega\\0 & K\end{bmatrix}\right)
=
\begin{bmatrix}E & F\\0 & E\end{bmatrix}.
```

This is implemented by `exp_block_upper_triangular_2x2`. Avoiding generic 4x4
`LinearAlgebra.exp` is important for MCM compatibility and is also faster for
ordinary `Float64` cases.

## Monte Carlo Measurements Compatibility

Several files are written to lift through
`MonteCarloMeasurements.Particles`:

- `material-properties.jl`
- `fiber-cross-section.jl`
- `path-geometry.jl`
- `fiber-path.jl`
- `fiber-path-modify.jl`
- `path-integral.jl`

Important conventions:

- avoid `::Real` annotations on uncertain-input slots,
- avoid `Float64(...)` coercions on paths that may carry `Particles`,
- avoid conditionals that would need to branch independently per particle,
- use `MonteCarloMeasurements.unsafe_comparisons(true)` in MCM tests,
- reduce ensemble-wide adaptive decisions through `scalar_reduce`.

`scalar_reduce` currently uses the maximum particle value when reducing an
MCM-valued scalar error metric. That makes the adaptive controller conservative:
the whole ensemble takes one step size selected by the worst particle. A
`pmean`-style reduction is the documented performance compromise if worst-case
step counts become too high.

Per-segment uncertainty and annotations live in the `meta` vector:

- `Nickname(label)` labels a segment for visual diagnostics,
- `MCMadd(symbol, distribution)` applies additive perturbations,
- `MCMmul(symbol, distribution)` applies multiplicative perturbations.

`fiber-path-modify.jl` interprets those annotations. For example, `:T_K`
metadata is converted into thermal length scaling using the cladding material
CTE at the fiber reference temperature.

## File Overview

- `material-properties.jl`: material refractive index, thermo-optic behavior,
  CTE, and nonlinear index helpers.
- `fiber-cross-section.jl`: step-index cross-section quantities, guided index,
  dispersion, nonlinear coefficient, and perturbative birefringence responses.
- `path-geometry.jl`: path authoring, placement, differential geometry,
  material twist resolution, sampling, and global path diagnostics.
- `path-geometry-connector.jl`: quintic G2 connector used by `JumpBy` and
  `JumpTo`.
- `path-geometry-plot.jl`: path plotting and HTML helpers.
- `fiber-path-meta.jl`: concrete per-segment metadata vocabulary.
- `fiber-path.jl`: `Fiber`, bend/twist generator assembly, and fiber-level
  diagnostics.
- `fiber-path-modify.jl`: meta-driven path perturbation and thermal scaling.
- `fiber-path-plot.jl`: fiber and propagation visualization helpers.
- `path-integral.jl`: Jones propagation, sensitivity propagation, DGD, and
  MCM-aware exponential formulas.
- `demo1.jl`: path geometry, segment labels, helix, modification, and
  adaptive-step visual demos.
- `demo2.jl`: `JumpBy` and `JumpTo` connector demos.
- `demo3mcm.jl`: MCM temperature/PTF demos.
- `demo3benchmark.jl`: MCM propagation benchmark demos.
- `test/`: unit and regression tests.

## Tests

The test orchestrator is:

```bash
julia --project=. test/runtests.jl
```

It currently includes:

- `test_path_geometry.jl`
- `test_fiber_path.jl`
- `test_fiber_path_modify.jl`
- `test_material_properties.jl`
- `test_mcm_compatability.jl`
- `test_paddle_transfer.jl`
- `test_dgd.jl`
- `test_fiber_cross_section.jl`
- `test_path_integral.jl`

Tests are intended to fall into the taxonomy described in the repository
`AGENTS.md`: `T-PHYSICS`, `T-VALIDATION`, `T-SIM-REGRESSION`,
`T-GUARDRAIL`, and visual demos. In practice, the strongest current coverage
is physics-motivated and guardrail testing. Validation against published fiber
data and direct legacy `fiber.py` comparisons is still limited.

MCM test blocks must use `MonteCarloMeasurements.unsafe_comparisons(true)`.
Under unsafe comparisons, structural checks such as sorted breakpoints are
reduced through particle means rather than failing on particle-valued booleans.

## Known Limitations

- `path-integral.jl` assumes lossless Jones propagation. Do not introduce
  gain, loss, or polarization-dependent loss there without a separate design.
- The model is intended for single-mode, weakly guiding, nearly circular fibers
  in the wavelength and temperature ranges described in the repository
  guidance.
- External stress, cladding noncircularity, non-concentric cores, nonlinear
  scattering, electric/magnetic effects, and polarization-dependent loss are
  not modeled.
- `modify(fiber)` handles geometry and metadata perturbations, including
  thermal scaling, but twist remapping through modification remains a known
  caveat area.
- The Julia port has substantial internal tests, but it is not yet a validated
  replacement for the legacy Python model or for published fiber data.

## Appendix: Cayley-Hamilton 2x2 Exponential

For any 2x2 matrix `A`, Cayley-Hamilton gives

```math
A^2-\operatorname{tr}(A)A+\det(A)I=0.
```

If `A` is traceless, then

```math
A^2=-\det(A)I.
```

Define `mu^2 = -det(A)`. Then

```math
A^2=\mu^2 I,\qquad A^3=\mu^2 A,\qquad A^4=\mu^4 I,\ldots
```

The exponential series splits into even and odd powers:

```math
e^A
=\sum_{k=0}^{\infty}\frac{A^{2k}}{(2k)!}
 +\sum_{k=0}^{\infty}\frac{A^{2k+1}}{(2k+1)!}.
```

Using `A^(2k) = mu^(2k) I` and `A^(2k+1) = mu^(2k) A`:

```math
e^A
=
\left(\sum_{k=0}^{\infty}\frac{\mu^{2k}}{(2k)!}\right)I
+
\left(\sum_{k=0}^{\infty}\frac{\mu^{2k}}{(2k+1)!}\right)A.
```

Therefore

```math
e^A=\cosh(\mu)I+\frac{\sinh(\mu)}{\mu}A,
\qquad \mu^2=-\det(A).
```

At `mu = 0`, interpret `sinh(mu) / mu -> 1`, so `e^A = I + A`.
