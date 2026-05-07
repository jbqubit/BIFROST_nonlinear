"""
fiber-path-modify.jl

Meta-driven Subpath modifier for a `Fiber`. Walks every authored interior
segment, applies any `MCMadd` / `MCMmul` annotations (via `MCMcombine`) to
the matching scalar fields, and rebuilds the Subpath(s) with the perturbed
inputs.

# Subpath terminal connector: lab-frame anchoring

Each `Subpath` ends with a structural terminal connector pinned at its
authored `jumpto_point` (a lab-frame invariant). Under `modify(fiber)` the
terminal connector is re-solved to land at the same `jumpto_point`
regardless of upstream perturbations, so the connector absorbs whatever
positional drift the modified interior segments would have introduced.
This re-pinning is what keeps long fibers with many meta-annotated
segments from drifting outside the author's intent.

If the Subpath's `jumpto_conserve_path_length == true`, the terminal
connector's arc length is additionally constrained so that the Subpath's
**total arc length** matches the baseline (pre-modify) value. This is the
mechanism for thermal-anchor and similar slack-conservation use cases.
The connector solves for an arc length of
`L_baseline - sum(modified interior arc lengths)`. The chord (= chord to
`jumpto_point`) and tangent constraints are unchanged; only the arc-length
search engages, which can change the connector's local shape (curvature
spikes possible).

!!! warning "min_bend_radius and birefringence"
    The slack a `conserve_path_length` connector absorbs shows up as
    **curvature variation**, even when the total path length is correctly
    preserved. Set `jumpto_min_bend_radius` on each conserve-path-length
    Subpath large enough that the resulting peak curvature stays well
    below the regime where bend-induced birefringence is consequential —
    otherwise length-conserving shape variations alone can dominate the
    polarization budget.

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
- `QuinticConnector` (already-resolved interior connector inside a
  `JumpBy`) → polynomial coeffs and `s_table` scale by `τ`

The Subpath's terminal connector is **structural**, not a segment. It
does not carry `meta` and so does not respond to `:T_K` directly. To
absorb thermal expansion of the Subpath's interior at a fixed lab-frame
endpoint, use `jumpto_conserve_path_length=true` on the Subpath.

Angles (`:angle`, `:axis_angle`) and counts (`:turns`) do not respond to
`:T_K`.

# Composition order (per scalar field)

1. `base_T = baseline * (1 + α_lin · ΔT)` if the field is T-sensitive and
   any `:T_K` MCM is present; otherwise `base_T = baseline`.
2. `final = MCMcombine(base_T, seg, :<field>)` (multiplicative first, then
   additive).

With no matching MCM entries `MCMcombine` returns its input unchanged, so
segments without annotations incur no arithmetic.

# Spinning overlays

Per-segment `Spinning` meta is carried through `modify(fiber)`: the rebuilt
Subpath's spinning runs are re-resolved against the perturbed segment arc
lengths via `_resolve_spinning_subpath_local`, so spinning anchors track the
modified geometry. Spinning `rate` and `phi_0` are not themselves perturbed by
MCM annotations.
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
    modify(fiber::Fiber) → Fiber

Return a new `Fiber` whose path has been perturbed according to the
`MCMadd` / `MCMmul` / `:T_K` annotations on each interior segment's `meta`
vector, plus any `jumpto_conserve_path_length=true` constraint on each
Subpath. See this file's module docstring for the per-field and `:T_K`
semantics.
"""
function modify(fiber::Fiber)
    T_ref = fiber.T_ref_K
    α_lin = cte(fiber.cross_section.cladding_material, T_ref)
    new_path = _modified_rebuild(fiber.path, T_ref, α_lin)
    return Fiber(new_path; cross_section = fiber.cross_section, T_ref_K = T_ref)
end

# ----------------------------
# Per-segment modification
# ----------------------------

# Compute the T_K length-scaling factor for this segment: (1 + α_lin·ΔT),
# returning `nothing` when no `:T_K` MCM is present (so the calling code
# can skip arithmetic via identity).
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
    return QuinticConnector(seg.a .* τ, seg.lambda, seg.s_table .* τ;
                            meta = seg.meta)
end

# Authored interior segments other than JumpBy go through the field-level
# modification dispatch above.
_modify_authored_segment(seg::JumpBy, ::Any, ::Any) = seg
_modify_authored_segment(seg::AbstractPathSegment, T_ref, α_lin) =
    _modify_segment(seg, T_ref, α_lin)

# ----------------------------
# Subpath rebuild
# ----------------------------

# Rebuild a single Subpath under meta perturbations. Uses the baseline
# `SubpathBuilt` only to read the baseline total arc length when the
# Subpath has `jumpto_conserve_path_length=true`. Otherwise the result is
# a function of the new Subpath's authored inputs alone.
function _modified_rebuild_subpath(sp_built::SubpathBuilt, T_ref, α_lin)
    sub = sp_built.subpath

    # Apply meta perturbations to each authored interior segment.
    new_interior = AbstractPathSegment[
        _modify_authored_segment(seg, T_ref, α_lin) for seg in sub.segments
    ]

    new_sub = Subpath(
        deepcopy(sub.meta),
        sub.start_point,
        sub.start_outgoing_tangent,
        sub.start_outgoing_curvature,
        new_interior,
        sub.jumpto_point,
        sub.jumpto_incoming_tangent,
        sub.jumpto_incoming_curvature,
        sub.jumpto_min_bend_radius,
        sub.jumpto_conserve_path_length,
    )

    # Walk the new interior segments with build()'s frame-advancing loop.
    pos = collect(new_sub.start_point)
    T_frame, N_frame, B_frame = _initial_frame_from_tangent(new_sub.start_outgoing_tangent)
    K_in_global = collect(new_sub.start_outgoing_curvature)
    placed = PlacedSegment[]
    s_eff = 0.0

    for seg_orig in new_sub.segments
        frame = hcat(N_frame, B_frame, T_frame)
        # JumpBy resolves at placement time. Other segments are passed through
        # _resolve_at_placement (which is identity for non-Jump segments).
        seg_placed = if seg_orig isa JumpBy
            seg_resolved = _resolve_at_placement(seg_orig, pos, frame, K_in_global)
            # Apply :T_K (and any other MCM) to the resolved QuinticConnector.
            _modify_segment(seg_resolved, T_ref, α_lin)
        else
            _resolve_at_placement(seg_orig, pos, frame, K_in_global)
        end
        push!(placed, PlacedSegment(seg_placed, s_eff, copy(pos), copy(frame)))

        pos_end_local         = end_position_local(seg_placed)
        (T_end_l, N_end_l, _) = end_frame_local(seg_placed)

        pos     = pos + frame * pos_end_local
        T_frame = _safe_normalize(frame * T_end_l)
        N_end_g = frame * N_end_l
        N_frame = _safe_normalize(N_end_g - dot(N_end_g, T_frame) * T_frame)
        B_frame = cross(T_frame, N_frame)

        L_seg = arc_length(seg_placed)
        κ_end = curvature(seg_placed, L_seg)
        K_in_global = κ_end .* N_frame

        s_eff = s_eff + L_seg
    end

    # Resolve the terminal connector. Optionally constrain its arc length
    # so the Subpath's total arc length matches the baseline.
    frame = hcat(N_frame, B_frame, T_frame)
    p1_local  = frame' * (collect(new_sub.jumpto_point) .- pos)
    chord     = norm(p1_local)
    t_hat_out = isnothing(new_sub.jumpto_incoming_tangent) ?
        (chord > 1e-15 ? p1_local ./ chord : [0.0, 0.0, 1.0]) :
        _safe_normalize(frame' * _safe_normalize(collect(new_sub.jumpto_incoming_tangent)))
    K0_local = frame' * K_in_global
    K1_local = isnothing(new_sub.jumpto_incoming_curvature) ?
        zeros(eltype(K0_local), 3) :
        frame' * collect(new_sub.jumpto_incoming_curvature)

    target = nothing
    if new_sub.jumpto_conserve_path_length
        L_baseline = Float64(_qc_nominalize(arc_length(sp_built)))
        s_eff_nom  = Float64(_qc_nominalize(s_eff))
        target = max(L_baseline - s_eff_nom, 0.0)
    end

    connector = _build_quintic_connector(
        p1_local, t_hat_out, K0_local, K1_local;
        min_bend_radius    = new_sub.jumpto_min_bend_radius,
        target_path_length = target,
    )
    jumpto_placed = PlacedSegment(connector, s_eff, copy(pos), copy(frame))
    s_end_eff = s_eff + arc_length(connector)

    resolved, pending = _resolve_spinning_subpath_local(
        placed, connector, Float64(_qc_nominalize(s_eff)),
        Float64(_qc_nominalize(s_end_eff)),
    )

    return SubpathBuilt(new_sub, placed, connector, jumpto_placed, resolved, pending)
end

# Top-level dispatch: SubpathBuilt or PathBuilt.
_modified_rebuild(path::SubpathBuilt, T_ref, α_lin) =
    _modified_rebuild_subpath(path, T_ref, α_lin)

function _modified_rebuild(path::PathBuilt, T_ref, α_lin)
    new_subs = [_modified_rebuild_subpath(sp, T_ref, α_lin) for sp in path.subpaths]
    # Re-stitch via build(::Vector{SubpathBuilt}). This re-runs the
    # endpoint-conformity check and any cross-Subpath spinning resolution.
    return build(new_subs)
end
