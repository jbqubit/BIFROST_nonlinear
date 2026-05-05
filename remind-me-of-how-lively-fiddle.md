# Plan: Subpath refactor, Pass 1 — path-geometry*.jl, hard cutover

## Context

This is the first of several passes that replace
`PathSpecBuilder`/`PathSpec`/`PathSpecCached` with a Subpath-centric
architecture. This pass touches **only** `path-geometry*.jl` (and renames
`fiber-path-meta.jl` → `path-geometry-meta.jl`).

**Hard cutover.** `PathSpec`, `PathSpecBuilder`, `PathSpecCached`, and
`JumpTo` are deleted in this pass. (`JumpBy` is **kept** as a regular interior
segment.) The rest of the codebase
(`fiber-path.jl`, `fiber-path-modify.jl`, `fiber-path-plot.jl`, demos, and
all tests outside `test_path_geometry*.jl`) will fail to load on this branch
until subsequent passes migrate them. This is intentional — we want
`path-geometry*.jl` clean and self-consistent in one shot rather than
carrying half-implemented duplicate types.

## Architectural Reframing

These principles drive every design decision in this pass:

1. **`SubpathBuilder` is the authoring entry point.** A Subpath specification
   begins with `start!(...)` and ends with `jumpto!(...)`. Between those two
   sealing calls the user appends interior segments (`straight!`, `bend!`,
   `helix!`, `catenary!`, `jumpby!`).

2. **A `Subpath` contains only user-supplied information.** Nothing is computed
   inside a `Subpath`. The terminal connector geometry is solved at
   `build()` time.

3. **Each `Subpath` is fully independent of all others.**

4. **`SubpathBuilt` is the result of `build(sp::Subpath)`.** It holds the
   `Subpath` plus its computed placement: the placed interior segments, the
   resolved terminal `QuinticConnector`, and the resolved twist runs.

5. **Each `SubpathBuilt` is fully independent of all others.** Its queries
   resolve from its own contents alone. The Subpath's `start_point` and
   `jumpto_point` are the only globally-anchored values it references.

6. **`PathBuilt` is a distinct ordered container of `SubpathBuilt`s.** It is
   the first layer in the stack that knows about Subpath ordering, and is
   therefore the natural site for any cross-Subpath reconciliation.

7. **There is no `Path` / `PathSpec` / `PathSpecCached` type.** The old
   distinctions are replaced by `Subpath`/`SubpathBuilt` (single) plus
   `PathBuilt` (ordered bag).

8. **Independence enables future parallelism in `propagate_fiber`.** Building
   geometry isn't the cost; ODE integration is. Subpath independence lets a
   later pass parallelize propagation across Subpaths.

---

## Phase 1 — Mechanical renames

### 1a. `target_arc_length` → `target_path_length`

- `path-geometry-connector.jl` — function signature, docstring, internal
  variable names, error messages.
- `path-geometry.jl` — internal call sites (most go away when JumpTo dispatch
  is removed; the kwarg lives only on `_build_quintic_connector`).

### 1b. `fiber-path-meta.jl` → `path-geometry-meta.jl`

- Rename the file.
- `path-geometry.jl` includes `path-geometry-meta.jl` directly (legitimate
  now that the meta types are part of the geometry layer).
- Update all `include("fiber-path-meta.jl")` references — those files won't
  build until subsequent passes anyway, but the include path needs to be
  correct.

---

## Phase 2 — New types replace old types

### 2a. `SubpathBuilder`

Two sealing calls delineate the builder lifecycle:

1. `start!(builder; ...)` — seals the start state. Must be called before any
   interior segments. Calling after a segment, or twice, throws.
2. `jumpto!(builder; ...)` — seals the end. Must be called after `start!`;
   further `straight!`/`bend!`/etc. or a second `jumpto!` after this throws.

```julia
mutable struct SubpathBuilder
    # Top-level annotation (e.g. Nickname; reserved for future plotting use)
    meta::Vector{AbstractMeta}
    # Start state (sealed by start!)
    start_point::Union{Nothing, NTuple{3, Float64}}              # default (0,0,0)
    start_outgoing_tangent::Union{Nothing, NTuple{3, Float64}}   # default (0,0,1)
    start_outgoing_curvature::Union{Nothing, NTuple{3, Float64}} # default (0,0,0)
    # Interior segments (straight, bend, helix, catenary, jumpby)
    segments::Vector{AbstractPathSegment}
    # Terminal connector spec (sealed by jumpto!)
    jumpto_point::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_tangent::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_curvature::Union{Nothing, NTuple{3, Float64}}
    jumpto_min_bend_radius::Union{Nothing, Float64}
    jumpto_conserve_path_length::Bool
end
```

No `jumpto_meta` field — the terminal connector is structural, not
meta-bearing.

**`start!(builder; point=(0,0,0), outgoing_tangent=(0,0,1),
        outgoing_curvature=(0,0,0))`**:
- Throws if `builder.start_point` is not Nothing (already sealed).
- Throws if `builder.segments` is non-empty.

**`jumpto!(builder; point, incoming_tangent=nothing,
         incoming_curvature=nothing, min_bend_radius=nothing,
         conserve_path_length=false)`**:
- Throws if `builder.start_point` is Nothing ("call start! first").
- Throws if `builder.jumpto_point` is not Nothing (already sealed).
- After this call, `start!`/`straight!`/`bend!`/etc. throw.

Existing segment functions (`straight!`, `bend!`, `helix!`, `catenary!`,
`jumpby!`) gain checks: throw if `start!` hasn't sealed the start, or if
`jumpto!` already sealed the end.

### 2b. `Subpath`

Validation lives in the constructor:

```julia
struct Subpath
    meta::Vector{AbstractMeta}
    start_point::NTuple{3, Float64}
    start_outgoing_tangent::NTuple{3, Float64}
    start_outgoing_curvature::NTuple{3, Float64}
    segments::Vector{AbstractPathSegment}   # interior, authored only
    jumpto_point::NTuple{3, Float64}
    jumpto_incoming_tangent::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_curvature::Union{Nothing, NTuple{3, Float64}}
    jumpto_min_bend_radius::Union{Nothing, Float64}
    jumpto_conserve_path_length::Bool
end

Subpath(b::SubpathBuilder)   # validates: start_point and jumpto_point not Nothing
```

A `Subpath` carries authored inputs only. Geometry queries on a `Subpath`
(arc_length, curvature, position, etc.) throw with "call build(subpath)
first" — analogous to today's throws on unresolved `JumpBy` segments.

### 2c. `SubpathBuilt`

```julia
struct SubpathBuilt
    subpath::Subpath
    placed_segments::Vector{PlacedSegment}      # interior segments only
    jumpto_quintic_connector::QuinticConnector  # the resolved terminal connector
    resolved_twists::Vector{ResolvedTwistRate}
    pending_continuous_first_twist::Bool        # true if first twist anchor has
                                                # is_continuous=true and was left
                                                # unresolved for PathBuilt to handle
end
```

- The terminal connector is **not** in `placed_segments`; it lives in its
  own field `jumpto_quintic_connector`. This matches the conceptual structure:
  `placed_segments` contains the interior; the jumpto is a distinct
  built-in feature of the Subpath.
- `s_end` is **not** stored — it's derived on demand as the sum of placed
  segment arc lengths plus the terminal connector arc length.
- `pending_continuous_first_twist` lets the PathBuilt build know it must
  resolve this Subpath's first twist anchor against the prior Subpath's
  terminal twist phase. See §2i.

### 2d. `PathBuilt`

A distinct struct, ordered container only. Stores no derived global
quantities — those are computed on demand to avoid consistency hazards:

```julia
struct PathBuilt
    subpaths::Vector{SubpathBuilt}
end

# Computed on demand. Lightweight.
s_offsets(p::PathBuilt) = ...   # cumulative starts; s_offsets(p)[1] == 0
s_end(p::PathBuilt)     = ...   # total global arc length

arc_length(p::PathBuilt) = s_end(p)
```

`s_offsets` and `s_end` are functions, not stored fields. They never go
stale because they're recomputed from the `subpaths` vector each time.
Cost is O(N) per call where N is the number of Subpaths — typically tiny.

### 2e. `build(sub::Subpath) → SubpathBuilt`

**Coordinate convention.** v1 stores `PlacedSegment.origin` and
`PlacedSegment.frame` in **global coordinates**, matching today's `build()`
exactly. The only difference is that `pos` is initialized at
`collect(sub.start_point)` rather than `zeros(3)`. Each `PlacedSegment.frame`
is the global frame at that segment's start, evolving forward as segments
are added. Subpath independence is preserved logically (each Subpath's
geometry is a function of its own inputs alone — `start_point`,
`start_outgoing_*`, segments, `jumpto_*`); independence does **not** require
storage to be Subpath-local.

The build inlines the old `_resolve_at_placement(::JumpTo)` logic for the
terminal connector. Since `JumpTo` is no longer a segment, no dispatch is
needed.

```julia
function build(sub::Subpath)
    # Initial state. The Subpath's start_point and start_outgoing_tangent are
    # global anchors held on the Subpath struct; SubpathBuilt's PlacedSegments
    # work in the same coordinate convention as today's build().
    pos         = collect(sub.start_point)
    T_frame     = _safe_normalize(collect(sub.start_outgoing_tangent))
    # Derive N_frame, B_frame from T_frame via Gram-Schmidt (as today)
    K_in_global = collect(sub.start_outgoing_curvature)
    placed      = PlacedSegment[]

    s_eff = 0.0  # local arc-length accumulator for placement records

    for seg_orig in sub.segments
        frame      = hcat(N_frame, B_frame, T_frame)
        seg_placed = _resolve_at_placement(seg_orig, pos, frame, K_in_global)
        push!(placed, PlacedSegment(seg_placed, s_eff, copy(pos), copy(frame)))
        # ... advance pos / frame / K_in_global / s_eff ...
    end

    # Resolve the terminal connector inline. Stored on its own field.
    frame    = hcat(N_frame, B_frame, T_frame)
    p1_local = frame' * (collect(sub.jumpto_point) .- pos)
    chord    = norm(p1_local)
    t_hat_out = isnothing(sub.jumpto_incoming_tangent) ?
        (chord > 1e-15 ? p1_local ./ chord : [0., 0., 1.]) :
        normalize(frame' * normalize(collect(sub.jumpto_incoming_tangent)))
    K0_local = frame' * K_in_global
    K1_local = isnothing(sub.jumpto_incoming_curvature) ? zeros(3) :
        frame' * collect(sub.jumpto_incoming_curvature)
    jumpto_quintic_connector = _build_quintic_connector(
        p1_local, t_hat_out, K0_local, K1_local;
        min_bend_radius = sub.jumpto_min_bend_radius)

    # Resolve twists local to this Subpath. Defer the first anchor if it is
    # is_continuous=true; the PathBuilt build handles cross-Subpath
    # continuity.
    resolved, pending_first = _resolve_twists_subpath_local(placed)

    return SubpathBuilt(sub, placed, jumpto_quintic_connector,
                       resolved, pending_first)
end

build(b::SubpathBuilder) = build(Subpath(b))
```

Note (`collect`): `start_point`, `start_outgoing_tangent`, and
`start_outgoing_curvature` are stored as `NTuple{3, Float64}` to match
JumpTo's existing field convention. The `collect(...)` calls convert tuples
to `Vector{Float64}` once, so subsequent `frame * v` arithmetic works with
plain matrices. (StaticArrays would let us drop the conversion entirely;
that's a separate refactor.)

`conserve_path_length=true` is a no-op in `build()` itself — there is no
baseline to compare against. It is consumed by the modify pipeline in the
later `_modified_rebuild` pass.

### 2f. Stitching builds into a `PathBuilt`

`build(::Vector{SubpathBuilt})` is the **first layer in the stack that knows
Subpath ordering**. v1 does three things:

1. **Conformity check between adjacent Subpaths.** For each consecutive
   pair `(N-1, N)` of Subpaths, validate:
   - `subpaths[N-1].subpath.jumpto_point == subpaths[N].subpath.start_point`
   - `subpaths[N-1].subpath.jumpto_incoming_tangent ==
      subpaths[N].subpath.start_outgoing_tangent`
   - `subpaths[N-1].subpath.jumpto_incoming_curvature ==
      subpaths[N].subpath.start_outgoing_curvature`

   Each check uses approximate equality (e.g. `≈` with a small tolerance).
   The check enforces the policy that Subpaths are **independent but
   stackable in a particular order**: independence does not relax C2
   continuity — it just says continuity is enforced *at stitch time*, not
   built into the Subpath specs.

   Mismatches throw `ArgumentError` with the offending pair index and
   field name.

2. **Cross-Subpath twist continuity** (when applicable). For each Subpath_N
   with `pending_continuous_first_twist == true`, look up Subpath_{N-1}'s
   terminal twist phase and resolve Subpath_N's first twist anchor with
   that phi_0. If `N == 1`, throw — there is no prior Subpath to inherit
   from.

3. Construct the `PathBuilt`.

```julia
function build(builts::Vector{SubpathBuilt})
    # 1. Conformity check
    for i in 2:length(builts)
        _check_subpath_conformity(builts[i-1].subpath, builts[i].subpath, i)
    end
    # 2. Cross-Subpath twist continuity
    builts = _resolve_pending_continuous_twists(builts)
    # 3. Construct
    return PathBuilt(builts)
end

build(subpaths::Vector{Subpath}) = build([build(sp) for sp in subpaths])
build(spb::SubpathBuilt)         = build(SubpathBuilt[spb])
```

The list comprehension in `build(::Vector{Subpath})` is the parallelization
seam (Threads.@threads or similar) when `propagate_fiber` arrives.

### 2g. JumpTo is deleted; JumpBy stays as a regular interior segment

`JumpTo` (struct, constructor, `_resolve_at_placement` method, the `@eval`
unresolved-error throws specific to it) is removed. The terminal connector
solve lives only in `build(::Subpath)` — there is no longer a JumpTo
segment dispatch.

`JumpBy` is unchanged from today's behavior. It remains an authored interior
segment (via `jumpby!`) that resolves to a `QuinticConnector` at placement
time inside `build()`'s segment loop. It does not seal the builder, does not
anchor a Subpath boundary, and is not treated as terminal.

The `@eval` block that defines unresolved-error stubs for `JumpBy` is
preserved (it stays as a guardrail for any code path that tries to query
arc_length/curvature on an authored-but-unresolved JumpBy).

### 2h. Twist scoping per Subpath

`_resolve_twists_subpath_local(placed)` (port of the existing
`_resolve_twists`) walks the placed segments of one Subpath. It differs in
one place: when the first twist anchor of a Subpath has
`is_continuous=true`, it does **not** throw. It leaves that anchor
unresolved (returning `pending_first=true`) and the PathBuilt build
resolves it from the prior Subpath's terminal phase.

If no such pending anchor exists, `pending_first=false` and the resolved
twists are complete after `build(sp::Subpath)`.

The "Subpath_1 with `is_continuous=true` on its first anchor" case throws
inside the PathBuilt build (not inside `build(sp::Subpath)`), with a clear
message that the first Subpath has no prior phase to inherit.

A standalone `build(sp::Subpath)` that produces a `SubpathBuilt` with
`pending_continuous_first_twist=true` is **valid** but not yet complete:
querying material twist on it raises until the PathBuilt build resolves it.
(See the TODO entry on `AbstractMeta.resolved` for a future
generalization.)

---

## Phase 3 — Query interface

### 3a. `SubpathBuilt` queries

Port of the existing `PathSpecCached` query surface, adapted for
`SubpathBuilt`'s structure (placed_segments + separate
jumpto_quintic_connector):

```julia
arc_length(b::SubpathBuilt) = sum(arc_length(ps.segment) for ps in b.placed_segments) +
                               arc_length(b.jumpto_quintic_connector)

arc_length(::SubpathBuilt, s1, s2) = s2 - s1

function _find_placed_segment(b::SubpathBuilt, s)
    cum = 0.0
    for ps in b.placed_segments
        seg_len = arc_length(ps.segment)
        if s <= cum + seg_len + 1e-12
            return ps, s - cum
        end
        cum += seg_len
    end
    # Past all interior segments → terminal connector
    s_local = s - cum
    return _terminal_placed_segment(b), clamp(s_local, 0.0,
        arc_length(b.jumpto_quintic_connector))
end
```

`_terminal_placed_segment(b)` synthesizes a `PlacedSegment` for the terminal
connector with origin/frame computed from the end-state of the last interior
placed segment (or from the Subpath start state if there are no interior
segments).

The point/vector queries (`curvature`, `geometric_torsion`,
`material_twist`, `position`, `tangent`, `normal`, `binormal`, `frame`,
`sample`, `sample_uniform`, `breakpoints`) follow the existing patterns,
adapted for the new structure.

### 3b. `PathBuilt` queries — single glue + generated forwards

A single `_find_subpath` glue method does the routing; an `@eval` loop
generates forwarding overloads for every "query at s" function:

```julia
function _find_subpath(p::PathBuilt, s)
    n      = length(p.subpaths)
    offs   = s_offsets(p)
    for i in 1:n
        local_end = offs[i] + arc_length(p.subpaths[i])
        if s <= local_end + 1e-12 || i == n
            return p.subpaths[i], s - offs[i]
        end
    end
    error("s out of path bounds")
end

# Forward all "point query at s" methods through _find_subpath in one shot.
for f in (:curvature, :geometric_torsion, :material_twist,
          :position, :tangent, :normal, :binormal, :frame)
    @eval function $f(p::PathBuilt, s::Real)
        sb, s_local = _find_subpath(p, s)
        return $f(sb, s_local)
    end
end
```

The few queries that don't fit the routing pattern get one-liners:

```julia
arc_length(p::PathBuilt)        = s_end(p)
arc_length(::PathBuilt, s1, s2) = s2 - s1

breakpoints(p::PathBuilt) = normalize_breakpoints(reduce(vcat,
    [breakpoints(p.subpaths[i]) .+ s_offsets(p)[i] for i in 1:length(p.subpaths)]))

sample(p::PathBuilt, s_values) = [frame(p, s) for s in s_values]

function sample_uniform(p::PathBuilt; n::Int = 256)
    return sample(p, range(0.0, s_end(p); length = n))
end
```

This pattern means:
- Adding a new "point query at s" function on `SubpathBuilt` requires only
  appending its name to the `@eval` tuple to gain `PathBuilt` support.
- The "doesn't fit" set is small and explicit.
- `_find_subpath` is the single global → local s translation.

### 3c. `path-geometry-plot.jl`

This file currently consumes `PathSpecCached`. Update it minimally so it
loads (its public API takes `SubpathBuilt` or `PathBuilt`); deeper plot
rework is deferred to a later pass.

---

## Phase 4 — Test migration

Hard-cutover tests migrate from `PathSpecBuilder`/`PathSpec`/`PathSpecCached`
to `SubpathBuilder`/`Subpath`/`SubpathBuilt`/`PathBuilt`.

### 4a. New guardrail tests in `test_path_geometry.jl`

- `start!` lifecycle: throws on second call, throws after a segment has
  been added.
- `jumpto!` lifecycle: throws before `start!`, throws on second call.
- Segment funcs throw if called before `start!` or after `jumpto!`.
- `Subpath(b)` throws if `start!` not called.
- `Subpath(b)` throws if `jumpto!` not called.
- Geometry queries on a `Subpath` (`arc_length`, `curvature`, etc.) throw
  with a clear message.
- `build(::Vector{Subpath})` rejects mismatched `jumpto_point[N-1]` vs
  `start_point[N]` with `ArgumentError`. Same for tangent and curvature.
- Stitch a pre-built `Vector{SubpathBuilt}` into a `PathBuilt` and verify
  it is identical (via deep field comparison) to `build(::Vector{Subpath})`
  applied to the same Subpaths.
- Twist continuity:
  - `is_continuous=true` first anchor in Subpath_1 throws at PathBuilt
    build time.
  - `is_continuous=true` first anchor in Subpath_{N≥2} resolves to the
    correct phi_0 inherited from Subpath_{N-1}'s terminal phase.
  - A standalone `build(sp::Subpath)` produces a `SubpathBuilt` with
    `pending_continuous_first_twist=true` when applicable; querying
    material twist on it raises.
- `start_outgoing_curvature` non-default: a Subpath with non-zero
  `start_outgoing_curvature` followed by an interior segment correctly
  threads the incoming curvature into the segment placement.

### 4b. Migration of existing tests in `test_path_geometry.jl`

Migrate every test using `PathSpecBuilder()`/`build(spec)` to the new API.
**Constraint: underlying logical assertions must remain non-trivial.** Do
not relax tolerances or weaken invariants to make migrated tests pass.

Mechanics:
1. `PathSpecBuilder()` → `SubpathBuilder()` + `start!()` (with the implicit
   defaults).
2. Interior segment calls (`straight!`, `bend!`, `helix!`, `catenary!`)
   are unchanged.
3. Append a terminal `jumpto!`. For tests that already had a terminal
   `JumpTo`, migrate kwargs: `destination=` → `point=`, `tangent=` →
   `incoming_tangent=`, `curvature_out=` → `incoming_curvature=`. For
   tests without a prior `JumpTo`, choose `point` along the natural exit
   ray and `incoming_tangent` equal to the natural exit tangent so the
   connector adds a small, predictable terminus.
4. `build(spec)` works directly on `SubpathBuilder`.
5. Update query assertions: where a test asserted a property of the whole
   path's arc length, restructure so the assertion targets a specific
   `placed_segment` length (or the `jumpto_quintic_connector` arc length)
   rather than the total. Where a test queried at a specific `s` inside the
   original geometry, the assertion stays unchanged (the connector lies
   beyond that s).

Tests that explicitly exercised JumpTo behavior (endpoint matches
destination, min_bend_radius validation, infeasibility throws, etc.)
migrate straightforwardly — they were already terminal-jumpto tests.

Tests using `JumpBy` (interior relative jump) keep their `jumpby!` calls
intact. JumpBy remains an interior segment in the new architecture. Only
the surrounding authoring boilerplate (`PathSpecBuilder` →
`SubpathBuilder` + `start!` + ... + `jumpto!`) changes.

`test_path_geometry_connector.jl` does **not** migrate to the new builder
API; it tests `_build_quintic_connector` directly. The
`target_arc_length` → `target_path_length` rename is its only change.

### 4c. Other tests

`test_fiber_path*.jl` and `test_mcm_compatability.jl` will fail to load on
this branch because they include `fiber-path.jl`, which still expects
`PathSpec*`. They're migrated in the next pass when `fiber-path.jl` and
related are updated.

This is the accepted cost of hard cutover. The branch's `runtests.jl`
will be red until the downstream pass lands.

```bash
julia --project=. julia-port/test/runtests.jl  # red until downstream pass lands
julia --project=. julia-port/test/test_path_geometry.jl
julia --project=. julia-port/test/test_path_geometry_connector.jl
```

The path-geometry and connector tests pass on their own at the end of this
pass.

---

## Out of scope (deferred to subsequent passes)

- **`Fiber` architecture rework.** A future `Fiber` is built from `Subpath`
  parts and holds an internal `PathBuilt` containing placed `SubpathBuilt`s.
  `SubpathBuilt` and `PathBuilt` are intended to be implementation details
  of `Fiber` from the end-user's perspective.
- Migration of `fiber-path.jl` to consume `SubpathBuilt`/`PathBuilt`.
- Migration of `_modified_rebuild` in `fiber-path-modify.jl`.
- Implementation of `conserve_path_length` semantics (lives in the modify
  pipeline).
- Migration of all `demo*.jl` files (each demo gets ceremonial
  `SubpathBuilder` / `start!` / `jumpto!`; existing `jumpby!` calls survive
  unchanged).
- Migration of `test_fiber_path*.jl`, `test_mcm_compatability.jl`.
- Plot updates in `path-geometry-plot.jl` (beyond a minimal load-fix) and
  `fiber-path-plot.jl`.
- Parallelization of `propagate_fiber` (the actual end goal of the
  independence invariant).
- `AbstractMeta.resolved` flag — added as a TODO entry; revisited when
  the downstream meta-processing audit happens.
