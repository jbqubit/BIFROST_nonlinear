"""
fiber-path-modify.jl

Meta-driven segment modifier for a `Fiber`. Walks every placed segment,
applies any `MCMadd` / `MCMmul` annotations (via `MCMcombine`) to the
matching scalar fields of that segment type, and rebuilds the path with
the perturbed segments re-placed.

# jumpto! terminal connector: lab-frame anchoring

The terminal connector (authored via `jumpto!` on a `SubpathBuilder`) **conserves
its lab-frame destination** under `modify(fiber)`: the endpoint stays fixed
regardless of what `:T_K` or field-level perturbations any upstream segment carries.
On long fibers with many meta-annotated segments this anchoring is essential. Without
it, each segment's perturbation propagates into every downstream position; the
accumulated drift compounds along the path and the geometry can swing well outside the
author's intent.

Note: in Phase 2, `:T_K` meta on the terminal connector's `jumpto_meta` is not
applied. Length-constrained connector behavior under `:T_K` is pending Phase 3.

# Field-level symbols (direct)

| Segment           | Scalar symbols                               |
|-------------------|----------------------------------------------|
| `StraightSegment` | `:length`                                    |
| `BendSegment`     | `:radius`, `:angle`, `:axis_angle`           |
| `CatenarySegment` | `:a`, `:length`, `:axis_angle`               |
| `HelixSegment`    | `:radius`, `:pitch`, `:turns`, `:axis_angle` |

All segment types also accept `:T_K`; see below.

# `:T_K` sugar (indirect)

`:T_K` is not a field; it expands to a thermal length-scaling of each
segment's T-sensitive fields, with
`α_lin = cte(fiber.cross_section.cladding_material, fiber.T_ref_K)` and
`ΔT = MCMcombine(T_ref, seg, :T_K) - T_ref`:

- `StraightSegment` → `:length`
- `BendSegment`     → `:radius`
- `CatenarySegment` → `:a`, `:length`
- `HelixSegment`    → `:radius`, `:pitch`
- `JumpBy`          → resolved connector geometry scales by `τ` (chord
  and arc both expand; `delta` is fiber-relative)
- `QuinticConnector` (already-resolved, e.g. inside a `JumpBy`)
  → polynomial coeffs and `s_table` scale by `τ`

Angles (`:angle`, `:axis_angle`) and counts (`:turns`) do not respond to
`:T_K`.

# Composition order (per scalar field)

1. `base_T = baseline * (1 + α_lin · ΔT)` if the field is T-sensitive and
   any `:T_K` MCM is present; otherwise `base_T = baseline`.
2. `final = MCMcombine(base_T, seg, :<field>)` (multiplicative first, then
   additive).

With no matching MCM entries `MCMcombine` returns its input unchanged, so
segments without annotations incur no arithmetic.

# Twist overlays

STUB: twist overlay remapping was removed pending the per-segment-meta twist
refactor. `modify(fiber)` no longer carries twist data through.
"""

if !@isdefined(Fiber)
    include("fiber-path.jl")
end
if !@isdefined(MCMadd)
    include("path-geometry-meta.jl")
end

# ----------------------------
# Public API
# ----------------------------

"""
    modify(fiber::Fiber) → SubpathCached

Return a new `SubpathCached` whose segments have been perturbed according
to the `MCMadd` / `MCMmul` / `:T_K` annotations on each segment's `meta`
vector. See this file's module docstring for the per-field and T_K semantics.
"""
function modify(fiber::Fiber)
    T_ref = fiber.T_ref_K
    α_lin = cte(fiber.cross_section.cladding_material, T_ref)
    return _modified_rebuild(fiber.path, T_ref, α_lin)
end

# ----------------------------
# Per-segment modification
# ----------------------------

# Compute the T_K length-scaling factor for this segment: (1 + α_lin·ΔT),
# returning the baseline type's multiplicative identity when no :T_K MCM is
# present (so the calling code can skip arithmetic via identity).
function _T_K_factor(seg, T_ref, α_lin)
    T_op = MCMcombine(T_ref, seg, :T_K)
    T_op === T_ref && return nothing
    return 1 + α_lin * (T_op - T_ref)
end

# Apply T_K scaling (if any) then the field-level MCM for `sym`.
_apply(base, seg, sym::Symbol, factor) =
    MCMcombine(isnothing(factor) ? base : base * factor, seg, sym)

# StraightSegment: length is T-sensitive.
function _modify_segment(seg::StraightSegment, T_ref, α_lin)
    τ = _T_K_factor(seg, T_ref, α_lin)
    new_length = _apply(seg.length, seg, :length, τ)
    return StraightSegment(new_length; meta = seg.meta)
end

# BendSegment: radius is T-sensitive; angle and axis_angle are not.
function _modify_segment(seg::BendSegment, T_ref, α_lin)
    τ = _T_K_factor(seg, T_ref, α_lin)
    new_radius     = _apply(seg.radius, seg, :radius, τ)
    new_angle      = MCMcombine(seg.angle, seg, :angle)
    new_axis_angle = MCMcombine(seg.axis_angle, seg, :axis_angle)
    return BendSegment(new_radius, new_angle, new_axis_angle; meta = seg.meta)
end

# CatenarySegment: a and length are T-sensitive; axis_angle is not.
function _modify_segment(seg::CatenarySegment, T_ref, α_lin)
    τ = _T_K_factor(seg, T_ref, α_lin)
    new_a          = _apply(seg.a, seg, :a, τ)
    new_length     = _apply(seg.length, seg, :length, τ)
    new_axis_angle = MCMcombine(seg.axis_angle, seg, :axis_angle)
    return CatenarySegment(new_a, new_length, new_axis_angle; meta = seg.meta)
end

# HelixSegment: radius and pitch are T-sensitive; turns and axis_angle are not.
function _modify_segment(seg::HelixSegment, T_ref, α_lin)
    τ = _T_K_factor(seg, T_ref, α_lin)
    new_radius     = _apply(seg.radius, seg, :radius, τ)
    new_pitch      = _apply(seg.pitch, seg, :pitch, τ)
    new_turns      = MCMcombine(seg.turns, seg, :turns)
    new_axis_angle = MCMcombine(seg.axis_angle, seg, :axis_angle)
    return HelixSegment(new_radius, new_pitch, new_turns, new_axis_angle;
                        meta = seg.meta)
end

# QuinticConnector: only :T_K is supported (no field-level MCMs). Scale the
# quintic polynomial coefficients and arc-length table by (1 + α_lin·ΔT).
# Geometric scaling r(u) → τ·r(u) implies all six coefficient rows scale by τ
# and the parameter-speed handle λ scales by τ too.
function _modify_segment(seg::QuinticConnector, T_ref, α_lin)
    τ = _T_K_factor(seg, T_ref, α_lin)
    isnothing(τ) && return seg
    # `lambda` is a Float64 metadata handle that records the search outcome
    # at build() time. Field-by-field rescaling only needs to touch the
    # geometry-bearing arrays; τ may carry MCM Particles which would not fit
    # the Float64 field anyway.
    return QuinticConnector(seg.a .* τ, seg.lambda, seg.s_table .* τ;
                            meta = seg.meta)
end

# Raw JumpBy / JumpTo should not appear inside a built Path (build() resolves
# them into QuinticConnector). Preserve the historical fallback.
_modify_segment(seg::AbstractPathSegment, ::Any, ::Any) =
    error("modify: cannot modify segment of type $(typeof(seg)); " *
          "expected a segment produced by build()")

# ----------------------------
# Internal: reassemble a path with replaced segments
# ----------------------------

function _modify_authored_segment(seg::JumpBy, ::Any, ::Any)
    return seg
end

_modify_authored_segment(seg::AbstractPathSegment, T_ref, α_lin) =
    _modify_segment(seg, T_ref, α_lin)

function _modified_rebuild(path::SubpathCached, T_ref, α_lin)
    sub = path.subpath

    # Re-place the authored segments, matching build()'s frame-advancing loop.
    # JumpBy is resolved after upstream modifications so G2 connector endpoints
    # inherit the modified incoming curvature. Connector-local :T_K meta is then
    # applied to the resolved QuinticConnector.
    pos         = collect(sub.origin_point)
    T_frame, N_frame, B_frame = _initial_frame(collect(sub.origin_outgoing_tangent))
    s_eff       = 0.0
    placed      = PlacedSegment[]
    K_in_global = collect(sub.origin_outgoing_curvature)

    for seg_orig in sub.segments
        frame = hcat(N_frame, B_frame, T_frame)

        seg_authored = _modify_authored_segment(seg_orig, T_ref, α_lin)
        seg_new = if seg_authored isa JumpBy
            # JumpBy delta is fiber-relative; chord-and-arc both scale by τ.
            seg_resolved = _resolve_at_placement(seg_authored, pos, frame, K_in_global)
            _modify_segment(seg_resolved, T_ref, α_lin)
        else
            _resolve_at_placement(seg_authored, pos, frame, K_in_global)
        end

        push!(placed, PlacedSegment(seg_new, s_eff, copy(pos), copy(frame)))

        pos_end_local         = end_position_local(seg_new)
        (T_end_l, N_end_l, _) = end_frame_local(seg_new)

        pos     = pos + frame * pos_end_local
        T_frame = _safe_normalize(frame * T_end_l)
        N_end_g = frame * N_end_l
        N_frame = _safe_normalize(N_end_g - dot(N_end_g, T_frame) * T_frame)
        B_frame = cross(T_frame, N_frame)

        L_seg = arc_length(seg_new)
        κ_end = curvature(seg_new, L_seg)
        K_in_global = κ_end .* N_frame

        s_eff = s_eff + L_seg
    end

    # Terminal connector — re-resolve from the subpath's jumpto fields.
    if !isnothing(sub.jumpto_point)
        frame    = hcat(N_frame, B_frame, T_frame)
        p1_local = frame' * (collect(sub.jumpto_point) .- pos)
        chord    = norm(p1_local)
        t_hat_out = if isnothing(sub.jumpto_incoming_tangent)
            chord > 1e-15 ? p1_local ./ chord : [0.0, 0.0, 1.0]
        else
            normalize(frame' * normalize(collect(sub.jumpto_incoming_tangent)))
        end
        K0_local = frame' * K_in_global
        K1_local = isnothing(sub.jumpto_incoming_curvature) ? zeros(3) :
                   frame' * collect(sub.jumpto_incoming_curvature)
        connector_new = _build_quintic_connector(p1_local, t_hat_out, K0_local, K1_local;
                            min_bend_radius = sub.jumpto_min_bend_radius,
                            meta            = sub.jumpto_meta)
        push!(placed, PlacedSegment(connector_new, s_eff, copy(pos), copy(frame)))
        s_eff += arc_length(connector_new)
    end

    has_twist = any(ps -> any(m -> m isa Twist, segment_meta(ps.segment)), placed)
    resolved_twists = has_twist ?
        _resolve_twists(placed, Float64(_qc_nominalize(s_eff))) :
        ResolvedTwistRate[]
    return SubpathCached(sub, placed, s_eff, resolved_twists)
end

function _replace_segments(path::SubpathCached, new_segments::Vector{<:AbstractPathSegment})
    @assert length(new_segments) == length(path.placed_segments)
    # Rebuild from the original subpath's origin, with new segments replacing the authored ones.
    # Note: the terminal connector (if any) is re-resolved from sub.jumpto_point in build().
    sub = path.subpath
    # Build a new SubpathBuilder with the same origin and new segments, then build.
    b = SubpathBuilder()
    b.origin_point             = sub.origin_point
    b.origin_outgoing_tangent  = sub.origin_outgoing_tangent
    b.origin_outgoing_curvature = sub.origin_outgoing_curvature
    # Strip the terminal connector out of new_segments: last segment is the connector
    # if the original subpath had a jumpto_point.
    n_segs = isnothing(sub.jumpto_point) ? length(new_segments) : length(new_segments) - 1
    for seg in new_segments[1:n_segs]
        push!(b.segments, seg)
    end
    if !isnothing(sub.jumpto_point)
        b.jumpto_point             = sub.jumpto_point
        b.jumpto_incoming_tangent  = sub.jumpto_incoming_tangent
        b.jumpto_incoming_curvature = sub.jumpto_incoming_curvature
        b.jumpto_min_bend_radius   = sub.jumpto_min_bend_radius
        b.jumpto_conserve_path_length = sub.jumpto_conserve_path_length
        b.jumpto_meta              = sub.jumpto_meta
    end
    return build(b)
end
