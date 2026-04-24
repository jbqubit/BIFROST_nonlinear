"""
Shrinkage transform for `Path` objects.

`path-geometry.jl` holds static, geometries. Shrinkage
(due to eg thermal contraction or drawing contraction) rescales segment dimensions. This
file provides `shrink(path, α) → Path` as a pure geometric transform: given an
input `Path` and a per-segment (or uniform) scaling factor `α`, produce a new
`Path` whose segment dimensions have been multiplied by `α`, with the placement
(origin and frame of each segment) recomputed from the scaled geometry.

The caller supplies `α` directly. A higher layer — not this file — is
responsible for deriving `α` from `(T_operating, T_ref_K, material)`.

# Semantics

- `BendSegment`: `radius → α·radius`. The swept `angle` is preserved. Arc length
  scales by `α`.
- `StraightSegment`, `CatenarySegment`: `length → α·length`. For catenary the
  parameter `a` also scales (`a → α·a`), preserving the catenary shape under
  scaling.
- `HelixSegment`: `radius → α·radius`, `pitch → α·pitch`, `turns` preserved.
- `HermiteConnector` (resolved form of `JumpBy`/`JumpTo`): scaled by rebuilding
  from the scaled Hermite polynomial coefficients.
- Twist overlays: the *overlay arc-length coordinates* are rescaled so that the
  same physical segment portion is covered. A constant twist rate τ (rad/m)
  stays constant under shrinkage; a function rate `τ(s)` is composed with the
  inverse coordinate map, `τ_new(s') = τ_old(s'/α)`.
"""

if !isdefined(Main, :Path)
    include("path-geometry.jl")
end

# ----------------------------
# Per-segment scaling
# ----------------------------

_scale_segment(seg::StraightSegment, α::Real) =
    StraightSegment(seg.length * α; nickname = seg.nickname)

_scale_segment(seg::BendSegment, α::Real) =
    BendSegment(seg.radius * α, seg.angle, seg.axis_angle; nickname = seg.nickname)

_scale_segment(seg::HelixSegment, α::Real) =
    HelixSegment(seg.radius * α, seg.pitch * α, seg.turns, seg.axis_angle;
                 nickname = seg.nickname)

_scale_segment(seg::CatenarySegment, α::Real) =
    CatenarySegment(seg.a * α, seg.length * α, seg.axis_angle; nickname = seg.nickname)

function _scale_segment(seg::HermiteConnector, α::Real)
    a0 = seg.a0 .* α
    a1 = seg.a1 .* α
    a2 = seg.a2 .* α
    a3 = seg.a3 .* α
    return HermiteConnector(a0, a1, a2, a3, seg.s_table .* α)
end

# Unplaceable segments (raw JumpBy/JumpTo) should not appear inside a built Path,
# but define a fallback for completeness.
_scale_segment(seg::AbstractPathSegment, ::Real) =
    error("shrink: cannot scale segment of type $(typeof(seg)); " *
          "expected a segment produced by build()")

# ----------------------------
# Public API
# ----------------------------

"""
    shrink(path::Path, α::Real) → Path

Return a new `Path` with every segment's dimensions scaled by the uniform
factor `α`. Twist overlay coordinates are rescaled so that the same physical
portion of each segment is covered; constant rates are preserved, function
rates are composed with the inverse coordinate map.
"""
function shrink(path::Path, α::Real)
    return _shrink_with(path, _ -> α)
end

"""
    shrink(path::Path, α_by_index::Dict{Int,<:Real}) → Path

Return a new `Path` with per-segment scaling. Segments whose 1-based index is
not in the dict are passed through with α = 1.0.
"""
function shrink(path::Path, α_by_index::AbstractDict{Int,<:Real})
    return _shrink_with(path, i -> get(α_by_index, i, 1.0))
end

# ----------------------------
# Internal: rebuild a scaled Path
# ----------------------------

function _shrink_with(path::Path, α_of_index)
    # 1. Scale each segment. α can be Float64 or Particles; keep element type
    #    flexible rather than coercing to Float64.
    n_seg = length(path.placed_segments)
    α_per_segment = [α_of_index(i) for i in 1:n_seg]
    scaled_segments = [_scale_segment(path.placed_segments[i].segment, α_per_segment[i])
                       for i in 1:n_seg]

    # 2. Re-place the scaled segments, replicating the frame-advancing logic of
    #    `build`. s_eff is initialised as 0.0; if any α is non-Float64 (e.g.
    #    Particles), promotion occurs through arc_length(seg_new).
    pos     = zeros(3)
    N_frame = [1.0, 0.0, 0.0]
    B_frame = [0.0, 1.0, 0.0]
    T_frame = [0.0, 0.0, 1.0]
    s_eff   = 0.0
    placed  = PlacedSegment[]

    s_old_starts = [ps.s_offset_eff for ps in path.placed_segments]
    s_new_starts = Vector{Any}(undef, n_seg)

    for (i, seg_new) in enumerate(scaled_segments)
        s_new_starts[i] = s_eff

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

    # 3. Rescale overlays. For each ResolvedTwistRate we find which segment it
    #    lives in (by construction it lives in exactly one — see
    #    `_resolve_overlay` in path-geometry.jl) and apply that segment's α.
    function _find_segment_index(s_old)
        for i in 1:n_seg
            seg_old_end = s_old_starts[i] + arc_length(path.placed_segments[i].segment)
            if s_old <= seg_old_end + 1e-12
                return i
            end
        end
        return n_seg
    end

    new_overlays = ResolvedTwistOverlay[]
    for ov in path.resolved_overlays
        new_rates = ResolvedTwistRate[]
        for r in ov.rates
            i = _find_segment_index(r.s_eff_start)
            αi = α_per_segment[i]
            s_old_off = r.s_eff_start - s_old_starts[i]
            s_old_len = r.s_eff_end   - r.s_eff_start
            s_new_start = s_new_starts[i] + αi * s_old_off
            s_new_end   = s_new_start + αi * s_old_len
            rate_old    = r.rate
            s_new_seg_start = s_new_starts[i]
            s_old_seg_start = s_old_starts[i]
            rate_new = let rate_old = rate_old, αi = αi,
                           s_new_seg_start = s_new_seg_start,
                           s_old_seg_start = s_old_seg_start
                s_new -> rate_old(s_old_seg_start + (s_new - s_new_seg_start) / αi)
            end
            push!(new_rates, ResolvedTwistRate(s_new_start, s_new_end, rate_new))
        end
        push!(new_overlays, ResolvedTwistOverlay(new_rates))
    end

    return Path(0.0, s_eff, placed, new_overlays)
end
