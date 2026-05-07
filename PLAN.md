# Plan: Separate shrinkage from static path geometry; put reference temperature on Fiber

## Context

`path-geometry.jl` is meant to hold a **static, immutable representation of geometry** — a space curve with segments, spinning overlays, and Frenet-frame queries. Today it also carries `shrinkage`: a per-segment scalar that rescales arc length and curvature. Shrinkage is a *fiber* phenomenon (thermal / drawing contraction), not a property of an abstract space curve, so it violates the file's stated intent.

Two related concerns follow from the same observation:

1. A fiber's physical length and geometry are **temperature-dependent**. The shrinkage layer needs a reference temperature so it has something to shrink relative to.
2. A `FiberCrossSection`'s optical properties are temperature-dependent. Today `T_K` is passed per call into `bending_birefringence` / `twisting_birefringence`, with no notion of a reference T at which core/cladding diameters are specified. That is a latent correctness gap when thermal expansion of the glass is eventually modeled.

The outcome we want:

- `path-geometry.jl` contains only pure geometry. No temperature, no `shrinkage` field, no `shrinkage` kwarg on builders, no `_apply_shrinkage_override`, no `nominal_arc_length` vs `arc_length` split.
- `FiberCrossSection` stays purely optical/geometric with no temperature field. `T_K` continues to be passed per call into `bending_birefringence` / `twisting_birefringence`.
- **Temperature lives on `Fiber`**, as a single scalar reference temperature captured at Fiber construction. This `T_ref_K` is the shared reference for both the geometric shrinkage layer and (eventually) thermal expansion of the cross section. A path has no temperature — it's just geometry. A cross section has no temperature — it's just optics.
- A new `fiber-path-shrinkage.jl` holds all shrinkage logic: a pure geometric transform `shrink(path::Path, α) → Path`. The caller supplies `α` directly; the file is purely geometric and does not know about temperature.

## Approach

### 1. Strip shrinkage out of path-geometry.jl

Remove from [path-geometry.jl](julia-port/path-geometry.jl):

- `shrinkage` field from every segment struct: `StraightSegment`, `BendSegment`, `CatenarySegment`, `HelixSegment`, `JumpBy`, `JumpTo`, `HermiteConnector`.
- `shrinkage` kwarg from every outer constructor and every `PathSpec` builder (`straight!`, `bend!`, `helix!`, `catenary!`, `jumpby!`, `jumpto!`).
- The `nominal_arc_length` interface — with no shrinkage, `arc_length` and `nominal_arc_length` coincide. Delete `nominal_arc_length` methods and all call sites inside `path-geometry.jl` (notably `_resolve_overlay`, which uses `ps.s_offset_nom`).
- `_apply_shrinkage_override` (L901–910) and the `shrinkage` kwarg on `build(spec; shrinkage=...)` (L948, L961–967).
- `segment_shrinkage` accessor (L92).
- `PlacedSegment.s_offset_nom` field — no longer meaningful when nominal == effective.
- SpinningOverlay handling in `_resolve_overlay` simplifies: `α = 1` everywhere, so `s_eff == s_nom` and the two offsets collapse.

All segment arithmetic that currently multiplies by `seg.shrinkage` (e.g. BendSegment `R = shrinkage * radius`) drops the factor.

### 2. Add `T_ref_K` to Fiber (not Path, not FiberCrossSection)

`Path` stays purely geometric — no temperature field. `FiberCrossSection` stays purely optical — no temperature field, and `bending_birefringence` / `twisting_birefringence` signatures are unchanged (still take an operating `T_K` per call).

Temperature enters the model **only at the `Fiber` level**, where geometry and cross section are bound together. Add a `T_ref_K::Float64` field to [Fiber](julia-port/fiber-path.jl) (L101–109):

```julia
struct Fiber{P,FTK}
    path::P
    cross_section::FiberCrossSection
    λ_m::Float64
    T_K::FTK               # operating temperature profile T(s)
    T_ref_K::Float64       # reference T for path geometry + cross-section dimensions
    s_start::Float64
    s_end::Float64
end
```

Add a `T_ref_K` kwarg to the `Fiber` constructor, defaulting to 297.15 K. Document in the docstring that this single value is the reference temperature for:

- **the path geometry** — the temperature at which `path`'s segment lengths/radii are valid. Future callers that want to model thermal contraction will build an effective path via `shrink(path, α)` with `α` derived from `(T_operating - T_ref_K)`.
- **the cross-section dimensions** — the temperature at which `core_diameter_m` and `cladding_diameter_m` are valid. A future thermal-expansion layer would derive effective diameters from `T_K(s) - T_ref_K`.

Neither derivation is implemented in this refactor; the field just captures the reference point so future work has somewhere to anchor. The operating `T_K(s)` profile continues to drive the per-call `bending_birefringence(..., T_K, ...)` arguments exactly as today.

### 3. New file: fiber-path-shrinkage.jl

Create [julia-port/fiber-path-shrinkage.jl](julia-port/fiber-path-shrinkage.jl). Responsibility: given a `Path` and a geometric shrinkage factor supplied by the caller, produce a new *effective* `Path` whose segment dimensions are scaled.

**Public API (decided):**

- `shrink(path::Path, α::Real) → Path` — uniform scaling.
- `shrink(path::Path, α_by_index::Dict{Int,<:Real}) → Path` — per-segment (matches today's `build(...; shrinkage=Dict)` override).

The caller supplies `α` directly. This file is **purely geometric** — it does not know about glass thermal expansion or `T_ref_K`. A future layer can convert `(T_operating, T_ref, material)` → `α`, but that lives elsewhere.

Under the hood, `shrink` rebuilds each segment with scaled dimensional fields (e.g. `BendSegment(radius * α, angle, axis_angle)` — angle is preserved, radius scales), re-places them via the same frame-advancing logic as `build`, and reconciles spinning overlays. The overlay reconciliation is the logic currently in `_resolve_overlay` (L912–936): overlays authored in reference coordinates are mapped into effective coordinates after scaling. That mapping moves from `build` into `shrink`.

Key semantic note: after this refactor, `build(spec)` always produces a reference-T `Path` where effective == nominal. `shrink` is the *only* way to introduce a divergence.

### 4. Update builders and callers

**Drop the `shrinkage` kwarg from every `PathSpec` builder entirely** (`straight!`, `bend!`, `helix!`, `catenary!`, `jumpby!`, `jumpto!`). The only way to introduce shrinkage is `shrink(build(spec), α)`. This forces call sites to restructure but keeps the authoring layer unambiguous.

- [fiber-path.jl](julia-port/fiber-path.jl): adds `T_ref_K` to `Fiber` (§2). No shrinkage changes — `Fiber` continues to operate on whatever `Path` it's handed (a reference-T path by convention).
- [test_path_geometry.jl](julia-port/test/test_path_geometry.jl): tests at L24–25, 28, 54–57, 67, 72, 135–141, 168, 170, 195, 198, 323–329, 334, 339, 525–528 move to a new `test_fiber_path_shrinkage.jl` and are rewritten against `shrink(build(spec), α)` instead of constructor/`build` kwargs.
- [test_mcm_compatability.jl](julia-port/test/test_mcm_compatability.jl): MCM tests at L267–268, 282–283, 324, 327, 333 that exercise uncertain (`Particles`) shrinkage move to the shrinkage test file. `Particles` must flow through `shrink` — the existing `promote`-based constructors already support it, so `shrink` rebuilding a segment with `α::Particles` should Just Work; verify in tests.

### 5. Out of scope

- Applying shrinkage inside `Fiber` (eager or lazy). `Fiber` continues to see a reference-T path for now.
- Thermal expansion of glass (changing effective `core_diameter_m` as a function of `T_K(s) - T_ref_K`). We only add the `T_ref_K` field to `Fiber` and document intent.
- A `T_operating → α` material model for path shrinkage.
- Changing `bending_birefringence` / `twisting_birefringence` signatures.
- Adding a temperature field to `Path` or `FiberCrossSection` — both stay pure.
- Anything in `path-integral.jl`.

## Critical files

- [julia-port/path-geometry.jl](julia-port/path-geometry.jl) — strip shrinkage; no temperature field added.
- [julia-port/fiber-cross-section.jl](julia-port/fiber-cross-section.jl) — **no changes**. Cross section stays purely optical.
- [julia-port/fiber-path.jl](julia-port/fiber-path.jl) — add `T_ref_K::Float64` field and constructor kwarg to `Fiber`.
- [julia-port/fiber-path-shrinkage.jl](julia-port/fiber-path-shrinkage.jl) — new. Pure geometric transform `shrink(path, α)`.
- [julia-port/test/test_path_geometry.jl](julia-port/test/test_path_geometry.jl) — remove shrinkage tests.
- [julia-port/test/test_mcm_compatability.jl](julia-port/test/test_mcm_compatability.jl) — move MCM shrinkage tests out.
- [julia-port/test/test_fiber_path.jl](julia-port/test/test_fiber_path.jl) — add a small test that `Fiber` round-trips `T_ref_K` and defaults it to 297.15.
- New: `julia-port/test/test_fiber_path_shrinkage.jl`.

## Verification

1. `julia --project=julia-port -e 'include("julia-port/test/test_path_geometry.jl")'` — passes with the shrinkage-free geometry.
2. `julia --project=julia-port -e 'include("julia-port/test/test_fiber_path_shrinkage.jl")'` — new tests cover uniform and per-index shrinkage, arc-length scaling, curvature invariance, spinning-overlay remapping, and MCM `Particles` compatibility.
3. `julia --project=julia-port -e 'include("julia-port/test/test_fiber_path.jl")'` and `test_dgd.jl`, `test_paddle_transfer.jl`, `test_mcm_compatability.jl` — all still pass. Most don't touch shrinkage; the MCM shrinkage cases have moved.
4. Spot-check: construct a `Fiber` without passing `T_ref_K`, confirm it defaults to 297.15; construct with `T_ref_K = 300.0`, confirm the field round-trips. `Path` and `FiberCrossSection` remain temperature-free.
5. `demo.jl` still runs end-to-end.

## Decisions locked in

- **Shrink API**: `shrink(path::Path, α) → Path` only (no spec-level variant).
- **α source**: Caller supplies `α` directly. `fiber-path-shrinkage.jl` is purely geometric.
- **Fiber shrinkage**: Out of scope. `Fiber` keeps operating on reference-T paths.
- **Builders**: Drop the `shrinkage` kwarg from every `PathSpec` builder. Shrinkage is introduced only via `shrink(...)`.
- **Temperature location**: `T_ref_K` lives on `Fiber` only — not on `Path`, not on `FiberCrossSection`. A path is pure geometry; a cross section is pure optics. Temperature is a property of the bound `Fiber` that combines them.
