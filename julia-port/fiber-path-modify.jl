"""
fiber-path-modify.jl

Meta-driven segment modifier for a `Fiber`. Walks every placed segment,
applies any `MCMadd` / `MCMmul` annotations (via `MCMcombine`) to the
matching scalar fields of that segment type, and rebuilds the path with
the perturbed segments re-placed.

# Field-level symbols (direct)

| Segment           | Scalar symbols                               |
|-------------------|----------------------------------------------|
| `StraightSegment` | `:length`                                    |
| `BendSegment`     | `:radius`, `:angle`, `:axis_angle`           |
| `CatenarySegment` | `:a`, `:length`, `:axis_angle`               |
| `HelixSegment`    | `:radius`, `:pitch`, `:turns`, `:axis_angle` |

# `:T_K` sugar (indirect)

`:T_K` is not a field; it expands to a thermal length-scaling of each
segment's T-sensitive fields, with
`α_lin = cte(fiber.cross_section.cladding_material, fiber.T_ref_K)` and
`ΔT = MCMcombine(T_ref, seg, :T_K) - T_ref`:

- `StraightSegment` → `:length`
- `BendSegment`     → `:radius`
- `CatenarySegment` → `:a`, `:length`
- `HelixSegment`    → `:radius`, `:pitch`
- `QuinticConnector`→ polynomial coeffs and `s_table`

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
    include("fiber-path-meta.jl")
end

# ----------------------------
# Public API
# ----------------------------

"""
    modify(fiber::Fiber) → PathSpecCached

Return a new `PathSpecCached` whose segments have been perturbed according
to the `MCMadd` / `MCMmul` / `:T_K` annotations on each segment's `meta`
vector. See this file's module docstring for the per-field and T_K semantics.
"""
function modify(fiber::Fiber)
    T_ref = fiber.T_ref_K
    α_lin = cte(fiber.cross_section.cladding_material, T_ref)
    new_segments = AbstractPathSegment[]
    for ps in fiber.path.placed_segments
        push!(new_segments, _modify_segment(ps.segment, T_ref, α_lin))
    end
    return _replace_segments(fiber.path, new_segments)
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

function _replace_segments(path::PathSpecCached, new_segments::Vector{<:AbstractPathSegment})
    n_seg = length(path.placed_segments)
    @assert length(new_segments) == n_seg

    # Re-place the segments, replicating build()'s frame-advancing loop.
    pos     = zeros(3)
    N_frame = [1.0, 0.0, 0.0]
    B_frame = [0.0, 1.0, 0.0]
    T_frame = [0.0, 0.0, 1.0]
    s_eff   = Float64(path.spec.s_start)
    placed  = PlacedSegment[]

    for seg_new in new_segments
        frame = hcat(N_frame, B_frame, T_frame)
        push!(placed, PlacedSegment(seg_new, s_eff, copy(pos), copy(frame)))

        pos_end_local         = end_position_local(seg_new)
        (T_end_l, N_end_l, _) = end_frame_local(seg_new)

        pos     = pos + frame * pos_end_local
        T_frame = _safe_normalize(frame * T_end_l)
        N_end_g = frame * N_end_l
        N_frame = _safe_normalize(N_end_g - dot(N_end_g, T_frame) * T_frame)
        B_frame = cross(T_frame, N_frame)

        s_eff = s_eff + arc_length(seg_new)
    end

    new_spec = PathSpec(collect(new_segments), path.spec.s_start)
    # Only resolve twists if any segment carries one. Skipping the resolve when
    # there are no anchors keeps MCM-valued (Particles) modify paths working;
    # _resolve_twists requires a Float64 s_end which is incompatible with
    # Particles-typed cumulative arc length.
    has_twist = any(ps -> any(m -> m isa Twist, segment_meta(ps.segment)), placed)
    resolved_twists = has_twist ? _resolve_twists(placed, Float64(s_eff)) : ResolvedTwistRate[]
    return PathSpecCached(new_spec, placed, s_eff, resolved_twists)
end
