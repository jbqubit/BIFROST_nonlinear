# Plan: Subpath refactor, Pass 4 — consumer migration to green runtests

## Context

Passes 1–3 landed the new `SubpathBuilder` / `Subpath` / `SubpathBuilt` /
`PathBuilt` architecture and the `Fiber` API surface that consumes them.
Pass 3 also implemented `conserve_path_length=true` and added a standalone
`test_fiber_path_pass3.jl` that's green.

What's still red on this branch:

- `julia-port/fiber-path-plot.jl` — works on `Fiber` and on sampled outputs;
  has no direct `PathSpec`/`JumpTo` references but currently fails to load
  in context because of dependency chains.
- `julia-port/demo-smallest.jl`, `demo.jl`, `demo1.jl`, `demo2.jl`,
  `demo3mcm.jl`, `demo3benchmark.jl` — all still on the old API.
- `julia-port/test/test_fiber_path.jl`, `test_fiber_path_modify.jl`,
  `test_mcm_compatability.jl`, `test_paddle_transfer.jl`, `test_dgd.jl`,
  `test_path_integral.jl` — still on the old API or otherwise blocked.
- `runtests.jl` — red.

Pass 4 migrates these consumers in **one pass** and brings `runtests.jl`
back to green. Per Pass 3's planning Q&A:

- **Delete `demo.jl`** (near-duplicate of `demo1.jl`). Port its unique
  `demo_helix_mcm_spinning` into `demo1.jl` first.
- **Strip the 6 geometry-only demos** out of `demo1.jl` (they already live
  in `demo-path-geometry.jl`). `demo1.jl` becomes the modify + adaptive
  step demo file.
- **Drop the 6 old `:T_K`-on-JumpTo testsets** in `test_fiber_path_modify.jl`
  (lines 327–415). The same physics is covered by Pass 3's
  `conserve_path_length` tests; do not migrate the obsolete authoring.

---

## Goals

1. `runtests.jl` runs green end-to-end.
2. All `demo*.jl` files run end-to-end and produce the HTML/SVG artifacts
   under `julia-port/output/`.
3. `demo1.jl` is the canonical "modify" demo file (no geometry-only demos).
   `demo.jl` is removed.

---

## Critical files

| File | Change |
|------|--------|
| `julia-port/fiber-path-plot.jl` | Verify it loads under the migrated stack; fix any incidental old-API references (esp. inside the `PlotRuntime` module's sampling path). |
| `julia-port/demo-smallest.jl` | Migrate to new API. |
| `julia-port/demo.jl` | **Delete** after porting `demo_helix_mcm_spinning` into `demo1.jl`. |
| `julia-port/demo1.jl` | Strip geometry demos; migrate modify + adaptive-step demos; absorb `demo_helix_mcm_spinning`. |
| `julia-port/demo2.jl` | Restructure each multi-jumpto path as a multi-Subpath `PathBuilt`. Migrate kwargs (`destination=` → `point=`, `tangent=` → `incoming_tangent=`, `curvature_out=` → `incoming_curvature=`). Translate the 3 thermal-anchor demos to `conserve_path_length=true`. |
| `julia-port/demo3mcm.jl` | Migrate path build + `Fiber()` + `modify()` to new API. |
| `julia-port/demo3benchmark.jl` | Migrate `Fiber()` + `modify()` calls. |
| `julia-port/test/test_fiber_path.jl` | Migrate 15 testsets to new builder API. |
| `julia-port/test/test_fiber_path_modify.jl` | Migrate ~21 testsets. **Delete** the 6 old JumpTo-`:T_K` testsets (lines 327–415). |
| `julia-port/test/test_mcm_compatability.jl` | Migrate any path-build references; mostly compatible. |
| `julia-port/test/test_paddle_transfer.jl` | Verify it still loads (no path-build code; depends on `propagate_fiber`). |
| `julia-port/test/test_dgd.jl`, `test_path_integral.jl` | Pure math tests; verify they still load. |
| `julia-port/test/runtests.jl` | Add `include("test_fiber_path_pass3.jl")`. Reorder if needed. |

Reused (no changes):

- `path-geometry.jl`, `path-geometry-meta.jl`, `path-geometry-connector.jl`
  (Pass 1).
- `path-geometry-plot.jl`, `demo-path-geometry.jl` (Pass 2).
- `fiber-path.jl`, `fiber-path-modify.jl`, `test_fiber_path_pass3.jl`
  (Pass 3).

API mapping reference (used throughout the migration):

| Old | New |
|---|---|
| `PathSpecBuilder()` | `SubpathBuilder()` + `start!(sb)` |
| `straight!(spec; …)` etc. | unchanged on the new builder |
| terminal `jumpto!(spec; destination=p, tangent=t, curvature_out=k, min_bend_radius=r)` | `jumpto!(sb; point=p, incoming_tangent=t, incoming_curvature=k, min_bend_radius=r, conserve_path_length=…)` |
| interior `jumpby!(spec; …)` | unchanged |
| `build(spec)` → `PathSpecCached` | `build(sb)` → `SubpathBuilt` (single subpath) or `build([sub1, sub2, …])` → `PathBuilt` |
| `path.spec.s_start` | `0.0` |
| `path.s_end` | `arc_length(path)` or `s_end(path)` |
| `path.placed_segments` | `b.placed_segments` (interior only); for "all placed including connector" use `_all_placed_segs(b)` from path-geometry.jl, or iterate `b.placed_segments` then `b.jumpto_placed` |
| `modify(fiber)` returning a path | `modify(fiber).path` |
| `:T_K` on JumpTo (thermal anchor) | `conserve_path_length=true` on the Subpath that ends at that destination |

---

## Implementation

### Phase 1 — `fiber-path-plot.jl`

Smoke-test load order: `julia --project=. -e 'include("julia-port/fiber-path-plot.jl"); println("OK")'`. The file has no direct `PathSpec`/`JumpTo` references in its top-level surface, but its `PlotRuntime` submodule samples from a `Fiber`. Walk through:

1. Confirm it loads. Fix any include chains pulling in deleted types.
2. `sample_fiber_centerline` and `sample_fiber_input` — verify they accept
   a `Fiber{<:SubpathBuilt}` and a `Fiber{<:PathBuilt}` and return the
   expected named-tuple structure. Adjust if they reach into
   `path.placed_segments`/`path.spec.s_start`/`path.s_end` directly.
3. `write_fiber_input_plot3d`, `render_poincare_sphere`, `render_pol_circle`,
   `collect_adaptive_steps`, `write_adaptive_steps_plot` — verify they
   still produce HTML output given a `Fiber` built from a `SubpathBuilt`.
4. If any function referenced `path.spec.s_start`/`path.s_end` directly,
   replace with `0.0` / `arc_length(f.path)` (or `f.s_start`/`f.s_end`).

### Phase 2 — `demo-smallest.jl`

Smallest demo. Pattern:

```julia
include("material-properties.jl")
include("fiber-cross-section.jl")
include("path-geometry.jl")
include("path-integral.jl")

xs = FiberCrossSection(...)

sb = SubpathBuilder(); start!(sb)
straight!(sb; length = 0.5, meta = [Nickname("lead-in")])
bend!(sb;     radius = 0.05, angle = π/2, meta = [Nickname("90 deg bend")])
straight!(sb; length = 0.5, meta = [Nickname("lead-out")])
# Seal at the natural exit. Compute the natural exit point/tangent
# analytically: after a quarter bend at R=0.05 starting on +z and a
# trailing straight of 0.5, end is (0.05 + 0.5, 0, 0.5 + 0.05) with
# tangent (1, 0, 0).
jumpto!(sb; point = (0.55, 0.0, 0.55), incoming_tangent = (1.0, 0.0, 0.0))

fiber = Fiber(build(sb); cross_section = xs, T_ref_K = 297.15)

J, stats = propagate_fiber(fiber; λ_m = 1550e-9, rtol = 1e-9, verbose = false)

println("J ="); display(J)
println("intervals = ", length(stats))
```

### Phase 3 — `demo.jl` deletion + `demo_helix_mcm_spinning` port

1. Read `demo.jl`'s `demo_helix_mcm_spinning` function body.
2. Append a migrated equivalent to `demo1.jl`, using the new builder API.
3. `git rm julia-port/demo.jl`.

### Phase 4 — `demo1.jl` migration

#### 4a. Strip geometry-only demos

Delete these 6 functions from `demo1.jl` (already in `demo-path-geometry.jl`):

- `demo_path_geometry`
- `demo_path_geometry_segment_labels`
- `demo_path_geometry_helix_0`
- `demo_path_geometry_helix_pi_3`
- `demo_path_geometry_helix_2pi_3`
- `demo_path_geometry_jumps_min_radius`

#### 4b. Migrate modify demos

12 modify demos + `demo_adaptive_step_doubling` + `demo_helix_mcm_spinning`
(ported from `demo.jl`). For each:

1. Convert `_build_modify_variant` (and its helix variant) helpers from
   `PathSpecBuilder` to `SubpathBuilder` + `start!` + `jumpto!`. The
   inverted-U baseline becomes a single Subpath sealed at its natural
   exit; the modified version is `modify(fiber).path` (one extra
   `.path` access).
2. Update `_sample_segment_xyz(path, seg_index; n)` to accept
   `Union{SubpathBuilt, PathBuilt}` and use `_all_placed_segs(path)` (or
   iterate `path.placed_segments` then `path.jumpto_placed`).
3. Update `path.spec.s_start` → `0.0` and `path.s_end` → `arc_length(path)`.

#### 4c. `demo_all` index

Update the `demo1.html` index so links still resolve (the geometry-only
entries are gone; entries for the 12 modify demos + adaptive-step +
helix-mcm-spinning remain).

### Phase 5 — `demo2.jl` migration

The structurally biggest demo migration. demo2 has ~20 demo functions;
each builds a paddle-style path that interleaves `straight!`/`jumpto!`/
`jumpby!` calls. In the new architecture each `jumpto!` seals a Subpath,
so the typical 4-jumpto paddle path becomes a 5-Subpath `PathBuilt`.

Reuse the pattern from `demo_path_geometry_jumps_min_radius` in
`demo-path-geometry.jl` as the migration template:

```julia
sb1 = SubpathBuilder(); start!(sb1)
straight!(sb1; …)
jumpto!(sb1; point = pin1, incoming_tangent = …, min_bend_radius = …)

sb2 = SubpathBuilder()
start!(sb2; point = pin1, outgoing_tangent = …)
straight!(sb2; …)
jumpto!(sb2; point = pin2, incoming_tangent = …)
# … etc.

p = build([Subpath(sb1), Subpath(sb2), …])
```

#### 5a. Geometry-only jump demos (15 functions)

Migrate kwargs: `destination=` → `point=`, `tangent=` → `incoming_tangent=`,
`curvature_out=` → `incoming_curvature=`. Restructure into multi-Subpath
PathBuilts.

#### 5b. Modify jump demos (3 functions)

- `demo_modify_jumpby_drift_2d`: pure JumpBy interior; stays a single
  Subpath. Migrate `Fiber(modify(f).path; …)` access patterns.
- `demo_modify_jumpto_anchor_2d`: the old version exercised
  destination-pinning under upstream `:T_K`. New: same physics via the
  Subpath's terminal connector being lab-frame-pinned.
  `conserve_path_length=false` (default) — the destination stays put but
  the connector arc length grows naturally.
- `demo_modify_jumpto_anchor_thermal_2d`: the old version exercised
  thermal expansion of the JumpTo connector itself. New: set
  `conserve_path_length=true` on the Subpath's terminal jumpto so the
  connector absorbs the length change instead of growing.

#### 5c. SVG/HTML helpers

`_jump_row_svg`, `_jump_row_html`, `_modify_overlay_svg` consume
`path.placed_segments` and similar fields. Update to handle the new
shape (interior segments + terminal connector, possibly across multiple
Subpaths).

#### 5d. `demo_all` index

Update `demo2.html` index entries.

### Phase 6 — `demo3mcm.jl` and `demo3benchmark.jl`

#### 6a. `demo3mcm.jl`

`_mcm_demo_fiber()` and `demo_mcm_temperature_ptf*()`:

- Replace `PathSpecBuilder()` + `build()` with the new builder pattern.
- Replace `modify(fiber)` (returning a path) with
  `modify(fiber).path` where the path is reused.
- Update plot code if it touches `placed_segments` directly.

#### 6b. `demo3benchmark.jl`

Same pattern, smaller surface (~10 lines).

### Phase 7 — Test migration

#### 7a. `test/test_fiber_path.jl` (15 testsets)

Mechanical migration:

1. `if !isdefined(Main, :PathSpecCached)` → `if !isdefined(Main, :SubpathBuilt)`.
2. `PathSpecBuilder()` → `SubpathBuilder()` + `start!(sb)`.
3. Append `jumpto!(sb; point = …, incoming_tangent = …)` before
   `build(sb)`.
4. `path.spec.s_start` → `0.0`; `path.s_end` → `arc_length(path)`.
5. `path.placed_segments` semantics: interior only. Tests asserting on
   "the last segment" of a path that previously included a JumpTo segment
   should now read `path.jumpto_placed.segment` for the connector.

#### 7b. `test/test_fiber_path_modify.jl` (21 testsets after deletion)

1. Delete the 6 testsets at lines 327–415 (old `:T_K`-on-JumpTo).
2. Migrate the remaining 21 testsets (radius/length/angle MCM
   perturbations on Straight/Bend/Catenary/Helix segments + JumpBy):
   - Same builder migration as 7a.
   - Replace `modify(fiber)` (path) with `modify(fiber).path` where the
     test then iterates segments.
   - `path.placed_segments[i].segment` semantics still apply for
     interior segments.

#### 7c. `test/test_mcm_compatability.jl` (14 testsets)

Most testsets exercise material-properties / fiber-cross-section / path
integral. The one testset that builds a path and runs `propagate_fiber`
needs the standard builder migration. Other testsets pass through.

#### 7d. `test/test_paddle_transfer.jl`, `test/test_dgd.jl`, `test/test_path_integral.jl`

Verify each loads cleanly under the new stack. test_paddle_transfer
depends on `propagate_fiber` only and should be untouched. test_dgd and
test_path_integral are pure math; should be untouched.

### Phase 8 — `runtests.jl`

Add `test_fiber_path_pass3.jl` to the include list. Final order:

```julia
include("test_path_geometry_connector.jl")
include("test_path_geometry.jl")
include("test_fiber_path.jl")
include("test_fiber_path_modify.jl")
include("test_fiber_path_pass3.jl")        # NEW
include("test_material_properties.jl")
include("test_fiber_cross_section.jl")
include("test_mcm_compatability.jl")
include("test_paddle_transfer.jl")
include("test_dgd.jl")
include("test_path_integral.jl")
```

---

## Verification

```bash
# Pass 4 acceptance — every command must succeed end to end:
julia --project=. julia-port/test/runtests.jl                        # GREEN

julia --project=. julia-port/demo-smallest.jl                        # prints J
julia --project=. julia-port/demo-path-geometry.jl                   # 7 HTML
julia --project=. -e 'include("julia-port/demo1.jl"); demo_all()'    # demo1.html + N HTML
julia --project=. -e 'include("julia-port/demo2.jl"); demo_all()'    # demo2.html + N HTML
julia --project=. -e 'include("julia-port/demo3mcm.jl"); demo3mcm_all()'        # demo3mcm.html
julia --project=. -e 'include("julia-port/demo3benchmark.jl"); demo3benchmark_all()'  # demo3benchmark.html

# Sanity: demo.jl is gone.
test ! -f julia-port/demo.jl && echo "demo.jl removed"
```

Expected outputs in `julia-port/../output/`:

- `demo-path-geometry-index.html` and the seven path-geometry HTMLs (Pass 2).
- `demo1.html` index plus 12 modify HTMLs + adaptive-step + helix-mcm-spinning.
- `demo2.html` index plus 18 jump HTMLs.
- `demo3mcm.html` plus the MCM temperature PTF artifacts.
- `demo3benchmark.html` plus benchmark artifacts.

After Pass 4 the only deferred items remaining are:

- Future: parallelize `propagate_fiber` across Subpaths.
- Future: `AbstractMeta.resolved` flag for end-to-end meta processing
  certification (see TODO.md).
- Optional: deeper `path-geometry-plot.jl` rework (the Pass 1 minimal
  load-fix is sufficient for this pass).
