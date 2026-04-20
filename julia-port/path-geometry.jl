"""
path-geometry.jl

Three-dimensional path geometry for smooth space curves.

This file  describes the shape of a 3D curve and its associated differential 
geometry.

# Assembly

Two approaches are supported and can be freely mixed:

(1) Sliding-frame approach: each segment is specified relative to the frame left by
    the previous segment. The tangent direction at the start of each new segment is
    exactly the tangent at the end of the previous one, so continuity is structural.

        spec = PathSpec()
        straight!(spec; length=1.0)
        bend!(spec; radius=0.05, angle=π/2) [ASK1]
        twist!(spec; s_start=0.3, length=0.8, rate=2π)          # constant
        twist!(spec; s_start=0.3, length=0.8, rate=s->sin(s))  # function

(2) Endpoint approach: specify the displacement (jumpby!) or absolute destination
    (jumpto!) and an optional outgoing tangent. The connecting segment is an Euler
    spiral (clothoid). The incoming tangent is always the current frame tangent.

        jumpby!(spec; delta=(0.0, 0.0, 0.5))
        jumpto!(spec; destination=(1.0, 0.0, 0.5), tangent=(0.0, 1.0, 0.0))

# Shrinkage

Each segment carries a `shrinkage` scalar (default 1.0) that uniformly rescales its
metric: arc lengths scale by α, curvature scales by 1/α, angles are preserved.
Shrinkage can be overridden at build time for parametric studies.

Twist overlays conserve turn count under shrinkage: if a nominal length L contains
N turns, an effective length α·L still contains N turns, so the effective twist rate
changes as τ(s) = 2π·N·α(s) / ∫α ds over the overlay interval.

# Interface

    arc_length(seg_or_path)
    arc_length(path, s1, s2)
    curvature(seg_or_path, s)
    geometric_torsion(seg_or_path, s)
    material_twist(path, s)
    position(path, s)
    tangent(path, s)
    normal(path, s)
    binormal(path, s)
    frame(path, s)
    start_point(path), end_point(path)
    start_tangent(path), end_tangent(path)
    path_length(path)
    cartesian_distance(path, s1, s2)
    bounding_box(path)
    total_turning_angle(path)
    total_torsion(path)
    total_material_twist(path)
    writhe(path)
    sample(path, s_values)
    sample_uniform(path; n)

# JSON

    path_from_json(filename) → PathSpec
    path_to_json(spec, filename)
"""

using LinearAlgebra

# -----------------------------------------------------------------------
# Abstract segment type
# -----------------------------------------------------------------------

abstract type AbstractPathSegment end

# Required interface for each concrete segment (local arc-length s ∈ [0, arc_length(seg)]):
#   arc_length(seg)                   → Float64
#   nominal_arc_length(seg)           → Float64  (arc length at shrinkage = 1)
#   curvature(seg, s)                 → Float64  (κ, 1/m)
#   geometric_torsion(seg, s)         → Float64  (τ_geom, rad/m)
#   position_local(seg, s)            → Vector{Float64} length 3
#   tangent_local(seg, s)             → Vector{Float64} unit, length 3
#   normal_local(seg, s)              → Vector{Float64} unit, length 3
#   binormal_local(seg, s)            → Vector{Float64} unit, length 3
#   end_position_local(seg)           → Vector{Float64} length 3
#   end_frame_local(seg)              → (T, N, B) each Vector{Float64} length 3

segment_shrinkage(seg::AbstractPathSegment) = seg.shrinkage

# -----------------------------------------------------------------------
# TwistOverlay
# -----------------------------------------------------------------------

"""
    TwistOverlay(s_start, length, rate)

Material twist applied over a nominal arc-length interval.
`s_start` and `length` are in nominal (shrinkage = 1) arc-length coordinates.
`rate` is a function `rate(s_eff) → Float64` giving the twist rate in rad/m
of effective arc-length.  Use `twist!` to construct overlays — it handles
conversion from scalar or function inputs.
"""
struct TwistOverlay
    s_start::Float64
    length::Float64
    rate::Function
end

# -----------------------------------------------------------------------
# StraightSegment
# -----------------------------------------------------------------------

struct StraightSegment <: AbstractPathSegment
    length::Float64
    shrinkage::Float64
end

StraightSegment(length::Real; shrinkage::Real = 1.0) =
    StraightSegment(Float64(length), Float64(shrinkage))

arc_length(seg::StraightSegment)         = seg.length * seg.shrinkage
nominal_arc_length(seg::StraightSegment) = seg.length
curvature(::StraightSegment, ::Real)     = 0.0
geometric_torsion(::StraightSegment, ::Real) = 0.0

position_local(::StraightSegment, s::Real)  = [0.0, 0.0, Float64(s)]
tangent_local(::StraightSegment, ::Real)        = [0.0, 0.0, 1.0]
normal_local(::StraightSegment, ::Real)         = [1.0, 0.0, 0.0]
binormal_local(::StraightSegment, ::Real)       = [0.0, 1.0, 0.0]
end_position_local(seg::StraightSegment)        = [0.0, 0.0, arc_length(seg)]
end_frame_local(::StraightSegment) = ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])

# -----------------------------------------------------------------------
# BendSegment  (circular arc)
# -----------------------------------------------------------------------

"""
    BendSegment(radius, angle, axis_angle; shrinkage)

Circular arc of radius `radius` (m) sweeping `angle` (rad) in the plane whose
inward normal is at `axis_angle` (rad) from the local N-axis.

In the local frame (local z = incoming tangent, local x = incoming normal,
local y = incoming binormal), the inward normal direction is:
    n̂ = cos(axis_angle)·x̂ + sin(axis_angle)·ŷ

Shrinkage scales the radius and arc length but preserves the swept angle.
Curvature κ = 1 / (shrinkage · radius).
"""
struct BendSegment <: AbstractPathSegment
    radius::Float64
    angle::Float64       # total angle swept (rad), preserved under shrinkage
    axis_angle::Float64  # orientation of inward normal in transverse plane (rad)
    shrinkage::Float64

    function BendSegment(radius::Real, angle::Real, axis_angle::Real = 0.0;
                         shrinkage::Real = 1.0)
        @assert radius > 0 "BendSegment: radius must be positive"
        new(Float64(radius), Float64(angle), Float64(axis_angle), Float64(shrinkage))
    end
end

arc_length(seg::BendSegment)         = seg.shrinkage * seg.radius * abs(seg.angle)
nominal_arc_length(seg::BendSegment) = seg.radius * abs(seg.angle)
curvature(seg::BendSegment, ::Real)  = 1.0 / (seg.shrinkage * seg.radius)
geometric_torsion(::BendSegment, ::Real) = 0.0

function position_local(seg::BendSegment, s::Real)
    R   = seg.shrinkage * seg.radius
    θ   = s / R
    φ   = seg.axis_angle
    n̂   = [cos(φ), sin(φ), 0.0]
    return R * (1 - cos(θ)) * n̂ + [0.0, 0.0, R * sin(θ)]
end

function tangent_local(seg::BendSegment, s::Real)
    R = seg.shrinkage * seg.radius
    θ = s / R
    φ = seg.axis_angle
    return [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)]
end

function normal_local(seg::BendSegment, s::Real)
    R = seg.shrinkage * seg.radius
    θ = s / R
    φ = seg.axis_angle
    return [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)]
end

function binormal_local(seg::BendSegment, ::Real)
    φ = seg.axis_angle
    return [-sin(φ), cos(φ), 0.0]   # constant for circular arc (zero torsion)
end

function end_position_local(seg::BendSegment)
    R = seg.shrinkage * seg.radius
    θ = seg.angle
    φ = seg.axis_angle
    n̂ = [cos(φ), sin(φ), 0.0]
    return R * (1 - cos(θ)) * n̂ + [0.0, 0.0, R * sin(θ)]
end

function end_frame_local(seg::BendSegment)
    θ = seg.angle
    φ = seg.axis_angle
    T = [sin(θ) * cos(φ), sin(θ) * sin(φ),  cos(θ)]
    N = [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)]
    B = [-sin(φ), cos(φ), 0.0]
    return (T, N, B)
end

# -----------------------------------------------------------------------
# CatenarySegment
# -----------------------------------------------------------------------

"""
    CatenarySegment(a, length, axis_angle; shrinkage)

A catenary curve in the plane whose horizontal direction is at `axis_angle`
from the local N-axis.  `a` (m) is the catenary parameter (a = T₀/(ρg),
the ratio of horizontal tension to weight per unit length).  The curve
starts with tangent along the local z-axis (vertical catenary vertex) and
curves horizontally.

Under shrinkage α: effective parameter a → α·a, arc length → α·length,
curvature κ(s) = a_eff / (a_eff² + s²).  Angles are preserved.
"""
struct CatenarySegment <: AbstractPathSegment
    a::Float64
    length::Float64
    axis_angle::Float64
    shrinkage::Float64

    function CatenarySegment(a::Real, length::Real, axis_angle::Real = 0.0;
                             shrinkage::Real = 1.0)
        @assert a > 0      "CatenarySegment: a must be positive"
        @assert length > 0 "CatenarySegment: length must be positive"
        new(Float64(a), Float64(length), Float64(axis_angle), Float64(shrinkage))
    end
end

arc_length(seg::CatenarySegment)         = seg.length * seg.shrinkage
nominal_arc_length(seg::CatenarySegment) = seg.length
geometric_torsion(::CatenarySegment, ::Real) = 0.0

function curvature(seg::CatenarySegment, s::Real)
    a = seg.a * seg.shrinkage
    return a / (a^2 + s^2)
end

function position_local(seg::CatenarySegment, s::Real)
    a = seg.a * seg.shrinkage
    φ = seg.axis_angle
    n̂ = [cos(φ), sin(φ), 0.0]
    horiz = a * (sqrt(1 + (s / a)^2) - 1)
    vert  = a * asinh(s / a)
    return horiz * n̂ + [0.0, 0.0, vert]
end

function tangent_local(seg::CatenarySegment, s::Real)
    a = seg.a * seg.shrinkage
    φ = seg.axis_angle
    q = sqrt(1 + (s / a)^2)
    return [(s / a) / q * cos(φ), (s / a) / q * sin(φ), 1.0 / q]
end

function normal_local(seg::CatenarySegment, s::Real)
    # N = dT/ds / |dT/ds|, derived analytically: N = [n̂_horiz/q, -s/a/q] normalised
    a = seg.a * seg.shrinkage
    φ = seg.axis_angle
    q = sqrt(1 + (s / a)^2)
    return [cos(φ) / q, sin(φ) / q, -(s / a) / q]
end

function binormal_local(seg::CatenarySegment, s::Real)
    return cross(tangent_local(seg, s), normal_local(seg, s))
end

function end_position_local(seg::CatenarySegment)
    return position_local(seg, arc_length(seg))
end

function end_frame_local(seg::CatenarySegment)
    s = arc_length(seg)
    return (tangent_local(seg, s), normal_local(seg, s), binormal_local(seg, s))
end

# -----------------------------------------------------------------------
# HelixSegment  (stub)
# -----------------------------------------------------------------------

"""
    HelixSegment(radius, pitch, turns, axis_angle; shrinkage)

A helix wound around a cylinder whose axis lies in the transverse plane at
`axis_angle` from the local N-axis.  `radius` is the helix radius (m),
`pitch` is the axial advance per turn (m), `turns` is the number of turns.

The fiber starts tangent along the local z-axis.  Because the helix tangent
is not generally along z in the natural parametrization, the segment applies
an internal coordinate rotation so the entry tangent is always aligned with
the incoming sliding frame.

    κ = R / (R² + (pitch/2π)²)
    τ_geom = (pitch/2π) / (R² + (pitch/2π)²)
    arc_length = turns · 2π · √(R² + (pitch/2π)²) · shrinkage

TODO: implement position_local, tangent_local, normal_local, binormal_local,
      end_frame_local.  The implementation requires choosing a consistent
      initial phase so that tangent(0) = ẑ in local coords.  See README.md
      for the relevant coordinate rotation derivation.
"""
struct HelixSegment <: AbstractPathSegment
    radius::Float64
    pitch::Float64
    turns::Float64
    axis_angle::Float64
    shrinkage::Float64

    function HelixSegment(radius::Real, pitch::Real, turns::Real,
                          axis_angle::Real = 0.0; shrinkage::Real = 1.0)
        @assert radius > 0 "HelixSegment: radius must be positive"
        @assert turns  > 0 "HelixSegment: turns must be positive"
        new(Float64(radius), Float64(pitch), Float64(turns),
            Float64(axis_angle), Float64(shrinkage))
    end
end

function _helix_h(seg::HelixSegment)
    seg.pitch / (2π)          # axial advance per radian
end

function arc_length(seg::HelixSegment)
    h = _helix_h(seg)
    return seg.turns * 2π * sqrt(seg.radius^2 + h^2) * seg.shrinkage
end

nominal_arc_length(seg::HelixSegment) = arc_length(seg) / seg.shrinkage

function curvature(seg::HelixSegment, ::Real)
    R = seg.radius * seg.shrinkage
    h = _helix_h(seg) * seg.shrinkage
    return R / (R^2 + h^2)
end

function geometric_torsion(seg::HelixSegment, ::Real)
    R = seg.radius * seg.shrinkage
    h = _helix_h(seg) * seg.shrinkage
    return h / (R^2 + h^2)
end

position_local(::HelixSegment, ::Real) =
    error("HelixSegment: position_local not yet implemented")
tangent_local(::HelixSegment, ::Real) =
    error("HelixSegment: tangent_local not yet implemented")
normal_local(::HelixSegment, ::Real) =
    error("HelixSegment: normal_local not yet implemented")
binormal_local(::HelixSegment, ::Real) =
    error("HelixSegment: binormal_local not yet implemented")
end_position_local(::HelixSegment) =
    error("HelixSegment: end_position_local not yet implemented")
end_frame_local(::HelixSegment) =
    error("HelixSegment: end_frame_local not yet implemented")

# -----------------------------------------------------------------------
# JumpBy and JumpTo  (stubs)
# -----------------------------------------------------------------------

"""
    JumpBy(delta, tangent_out, shrinkage)

Connects the current position to current_position + shrinkage·delta using an
Euler spiral (clothoid).  The incoming tangent is the current sliding frame
tangent.  `tangent_out` is the desired outgoing tangent direction; if nothing,
the connector computes an implicit tangent that minimizes curvature variation
(requires a nonlinear solve).

The Euler spiral respects the incoming curvature κ_in for C² continuity.

TODO: implement the Euler spiral solver.
"""
struct JumpBy <: AbstractPathSegment
    delta::NTuple{3, Float64}
    tangent_out::Union{Nothing, NTuple{3, Float64}}
    shrinkage::Float64
end

function JumpBy(delta; tangent_out = nothing, shrinkage::Real = 1.0)
    d = (Float64(delta[1]), Float64(delta[2]), Float64(delta[3]))
    t = isnothing(tangent_out) ? nothing :
        (Float64(tangent_out[1]), Float64(tangent_out[2]), Float64(tangent_out[3]))
    JumpBy(d, t, Float64(shrinkage))
end

"""
    JumpTo(destination, tangent_out, shrinkage)

Connects the current position to the fixed lab-frame `destination` using an
Euler spiral.  Shrinkage changes the arc length of the connector (the
connector absorbs the geometry change) but does not move the destination.

TODO: implement the Euler spiral solver.
"""
struct JumpTo <: AbstractPathSegment
    destination::NTuple{3, Float64}
    tangent_out::Union{Nothing, NTuple{3, Float64}}
    shrinkage::Float64
end

function JumpTo(destination; tangent_out = nothing, shrinkage::Real = 1.0)
    d = (Float64(destination[1]), Float64(destination[2]), Float64(destination[3]))
    t = isnothing(tangent_out) ? nothing :
        (Float64(tangent_out[1]), Float64(tangent_out[2]), Float64(tangent_out[3]))
    JumpTo(d, t, Float64(shrinkage))
end

for T in (JumpBy, JumpTo)
    @eval begin
        arc_length(::$T) = error($(string(T)) * ": not yet implemented (requires Euler spiral solver)")
        nominal_arc_length(::$T) = error($(string(T)) * ": not yet implemented")
        curvature(::$T, ::Real) = error($(string(T)) * ": not yet implemented")
        geometric_torsion(::$T, ::Real) = 0.0
        position_local(::$T, ::Real) = error($(string(T)) * ": not yet implemented")
        tangent_local(::$T, ::Real) = error($(string(T)) * ": not yet implemented")
        normal_local(::$T, ::Real) = error($(string(T)) * ": not yet implemented")
        binormal_local(::$T, ::Real) = error($(string(T)) * ": not yet implemented")
        end_position_local(::$T) = error($(string(T)) * ": not yet implemented")
        end_frame_local(::$T) = error($(string(T)) * ": not yet implemented")
    end
end

# -----------------------------------------------------------------------
# PathSpec  (mutable authoring struct)
# -----------------------------------------------------------------------

mutable struct PathSpec
    segments::Vector{AbstractPathSegment}
    twist_overlays::Vector{TwistOverlay}
    PathSpec() = new(AbstractPathSegment[], TwistOverlay[])
end

function straight!(spec::PathSpec; length::Real, shrinkage::Real = 1.0)
    push!(spec.segments, StraightSegment(length; shrinkage))
    return spec
end

function bend!(spec::PathSpec; radius::Real, angle::Real, axis_angle::Real = 0.0,
               shrinkage::Real = 1.0)
    push!(spec.segments, BendSegment(radius, angle, axis_angle; shrinkage))
    return spec
end

function helix!(spec::PathSpec; radius::Real, pitch::Real, turns::Real,
                axis_angle::Real = 0.0, shrinkage::Real = 1.0)
    push!(spec.segments, HelixSegment(radius, pitch, turns, axis_angle; shrinkage))
    return spec
end

function catenary!(spec::PathSpec; a::Real, length::Real, axis_angle::Real = 0.0,
                   shrinkage::Real = 1.0)
    push!(spec.segments, CatenarySegment(a, length, axis_angle; shrinkage))
    return spec
end

function jumpby!(spec::PathSpec; delta, tangent = nothing, shrinkage::Real = 1.0)
    push!(spec.segments, JumpBy(delta; tangent_out = tangent, shrinkage))
    return spec
end

function jumpto!(spec::PathSpec; destination, tangent = nothing, shrinkage::Real = 1.0)
    push!(spec.segments, JumpTo(destination; tangent_out = tangent, shrinkage))
    return spec
end

"""
    twist!(spec; s_start, length, rate)

Add a material twist overlay.  `s_start` and `length` are in nominal
(shrinkage = 1) arc-length coordinates.

`rate` is the material twist rate in rad/m of effective arc-length:
- `Real` — converted to a constant-rate function `_ -> rate`
- `Function` — `rate(s_eff)` evaluated at effective arc-length `s_eff`
"""
function twist!(spec::PathSpec; s_start::Real, length::Real, rate)
    r = rate isa Function ? rate : let τ = Float64(rate); _ -> τ; end
    push!(spec.twist_overlays, TwistOverlay(Float64(s_start), Float64(length), r))
    return spec
end

# -----------------------------------------------------------------------
# PlacedSegment and resolved twist data used by Path
# -----------------------------------------------------------------------

struct PlacedSegment
    segment::AbstractPathSegment
    s_offset_eff::Float64       # cumulative effective arc-length at segment start
    s_offset_nom::Float64       # cumulative nominal arc-length at segment start
    origin::Vector{Float64}     # global start position
    frame::Matrix{Float64}      # 3×3, columns [N_global | B_global | T_global]
                                # transforms local vectors → global: v_g = frame * v_l
end

struct ResolvedTwistRate
    s_eff_start::Float64
    s_eff_end::Float64
    rate::Function                     # τ_mat(s_eff) in rad/m
end

struct ResolvedTwistOverlay
    rates::Vector{ResolvedTwistRate}
end

# -----------------------------------------------------------------------
# Path  (immutable built struct)
# -----------------------------------------------------------------------

struct Path
    s_start::Float64
    s_end::Float64
    placed_segments::Vector{PlacedSegment}
    resolved_overlays::Vector{ResolvedTwistOverlay}
end

# -----------------------------------------------------------------------
# build()
# -----------------------------------------------------------------------

function _apply_shrinkage_override(seg::S, override) where {S <: AbstractPathSegment}
    isnothing(override) && return seg
    α = Float64(override)
    # Reconstruct with new shrinkage — uses positional fields common to all segments.
    # This works because every concrete segment stores shrinkage as the last field.
    flds = fieldnames(S)
    vals = [f === :shrinkage ? α : getfield(seg, f) for f in flds]
    return S(vals...)
end

function _resolve_overlay(overlay::TwistOverlay,
                          placed::Vector{PlacedSegment})
    nom_start = overlay.s_start
    nom_end   = overlay.s_start + overlay.length

    rates = ResolvedTwistRate[]
    for ps in placed
        seg = ps.segment
        seg_nom_start = ps.s_offset_nom
        seg_nom_end   = ps.s_offset_nom + nominal_arc_length(seg)
        overlap_nom_start = max(nom_start, seg_nom_start)
        overlap_nom_end   = min(nom_end,   seg_nom_end)
        overlap_nom_end <= overlap_nom_start && continue

        α = seg.shrinkage
        eff_in_seg_start = (overlap_nom_start - seg_nom_start) * α
        eff_in_seg_end   = (overlap_nom_end   - seg_nom_start) * α
        s_eff_start = ps.s_offset_eff + eff_in_seg_start
        s_eff_end   = ps.s_offset_eff + eff_in_seg_end

        push!(rates, ResolvedTwistRate(s_eff_start, s_eff_end, overlay.rate))
    end

    return ResolvedTwistOverlay(rates)
end

"""
    build(spec; shrinkage = nothing) → Path

Compile a `PathSpec` into an immutable `Path`.  Optional `shrinkage` argument:
- `nothing`  — use each segment's authored shrinkage value (default)
- `Float64`  — uniform override applied to all segments
- `Dict{Int,Float64}` — override by segment index (1-based)
"""
function build(spec::PathSpec; shrinkage = nothing)
    pos     = zeros(3)
    N_frame = [1.0, 0.0, 0.0]
    B_frame = [0.0, 1.0, 0.0]
    T_frame = [0.0, 0.0, 1.0]

    s_eff = 0.0
    s_nom = 0.0
    placed = PlacedSegment[]

    for (i, seg_orig) in enumerate(spec.segments)
        # Apply build-time shrinkage override
        seg = if isnothing(shrinkage)
            seg_orig
        elseif shrinkage isa Real
            _apply_shrinkage_override(seg_orig, shrinkage)
        elseif shrinkage isa Dict && haskey(shrinkage, i)
            _apply_shrinkage_override(seg_orig, shrinkage[i])
        else
            seg_orig
        end

        frame = hcat(N_frame, B_frame, T_frame)   # columns: [N | B | T]
        push!(placed, PlacedSegment(seg, s_eff, s_nom, copy(pos), copy(frame)))

        # Advance position and frame
        pos_end_local         = end_position_local(seg)
        (T_end_l, N_end_l, _) = end_frame_local(seg)

        pos     = pos + frame * pos_end_local
        T_frame = normalize(frame * T_end_l)
        N_end_g = frame * N_end_l
        # Re-orthogonalise to prevent floating-point drift
        N_frame = normalize(N_end_g - dot(N_end_g, T_frame) * T_frame)
        B_frame = cross(T_frame, N_frame)

        s_eff += arc_length(seg)
        s_nom += nominal_arc_length(seg)
    end

    resolved = [_resolve_overlay(ov, placed) for ov in spec.twist_overlays]

    return Path(0.0, s_eff, placed, resolved)
end

# -----------------------------------------------------------------------
# Segment lookup helpers
# -----------------------------------------------------------------------

function _find_placed_segment(path::Path, s::Real)
    sf = Float64(s)
    n = length(path.placed_segments)
    for i in 1:n
        ps = path.placed_segments[i]
        s_end = ps.s_offset_eff + arc_length(ps.segment)
        if sf <= s_end + 1e-12 || i == n
            s_local = clamp(sf - ps.s_offset_eff, 0.0, arc_length(ps.segment))
            return ps, s_local
        end
    end
    error("s = $s out of path bounds [$(path.s_start), $(path.s_end)]")
end

function _local_to_global(ps::PlacedSegment, v_local::Vector{Float64})
    return ps.frame * v_local
end

# -----------------------------------------------------------------------
# Differential geometry interface on Path
# -----------------------------------------------------------------------

arc_length(path::Path) = path.s_end - path.s_start

function arc_length(::Path, s1::Real, s2::Real)
    @assert s2 >= s1 "arc_length: require s2 >= s1"
    return Float64(s2) - Float64(s1)
end

function curvature(path::Path, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return curvature(ps.segment, s_local)
end

function geometric_torsion(path::Path, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return geometric_torsion(ps.segment, s_local)
end

function material_twist(path::Path, s::Real)
    sf = Float64(s)
    τ = 0.0
    for ov in path.resolved_overlays
        for r in ov.rates
            if r.s_eff_start <= sf <= r.s_eff_end
                τ += r.rate(sf)
            end
        end
    end
    return τ
end

function position(path::Path, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return ps.origin + _local_to_global(ps, position_local(ps.segment, s_local))
end

function tangent(path::Path, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return _local_to_global(ps, tangent_local(ps.segment, s_local))
end

function normal(path::Path, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return _local_to_global(ps, normal_local(ps.segment, s_local))
end

function binormal(path::Path, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return _local_to_global(ps, binormal_local(ps.segment, s_local))
end

function frame(path::Path, s::Real)
    T = tangent(path, s)
    N = normal(path, s)
    B = binormal(path, s)
    κ = curvature(path, s)
    τ = geometric_torsion(path, s)
    m = material_twist(path, s)
    return (; position = position(path, s), tangent = T, normal = N, binormal = B,
              curvature = κ, geometric_torsion = τ, material_twist = m)
end

# -----------------------------------------------------------------------
# Endpoint access
# -----------------------------------------------------------------------

start_point(path::Path)   = position(path, path.s_start)
end_point(path::Path)     = position(path, path.s_end)
start_tangent(path::Path) = tangent(path, path.s_start)
end_tangent(path::Path)   = tangent(path, path.s_end)

# -----------------------------------------------------------------------
# Path measures
# -----------------------------------------------------------------------

path_length(path::Path) = arc_length(path)

function cartesian_distance(path::Path, s1::Real, s2::Real)
    return norm(position(path, s2) - position(path, s1))
end

function bounding_box(path::Path; n::Int = 512)
    ss = range(path.s_start, path.s_end; length = n)
    pts = [position(path, s) for s in ss]
    lo = minimum(reduce(hcat, pts); dims = 2) |> vec
    hi = maximum(reduce(hcat, pts); dims = 2) |> vec
    return (; lo, hi)
end

function total_turning_angle(path::Path)
    # ∫κ ds  — for analytic segments computed exactly where possible
    total = 0.0
    for ps in path.placed_segments
        seg = ps.segment
        if seg isa StraightSegment
            # κ = 0
        elseif seg isa BendSegment
            total += abs(seg.angle)   # ∫κ ds = (1/R_eff) * R_eff*θ = θ (preserved)
        elseif seg isa HelixSegment
            total += curvature(seg, 0.0) * arc_length(seg)  # constant κ
        else
            # Numerical fallback for CatenarySegment and others
            n = 64
            ss = range(0.0, arc_length(seg); length = n + 1)
            h = ss[2] - ss[1]
            total += h * sum(curvature(seg, s) for s in ss)
        end
    end
    return total
end

function total_torsion(path::Path)
    total = 0.0
    for ps in path.placed_segments
        seg = ps.segment
        if seg isa HelixSegment
            total += geometric_torsion(seg, 0.0) * arc_length(seg)  # constant τ
        elseif seg isa (StraightSegment) || seg isa BendSegment || seg isa CatenarySegment
            # τ_geom = 0 for these
        else
            n = 64
            ss = range(0.0, arc_length(seg); length = n + 1)
            h = ss[2] - ss[1]
            total += h * sum(geometric_torsion(seg, s) for s in ss)
        end
    end
    return total
end

function total_material_twist(path::Path; n_quad::Int = 128)
    total = 0.0
    for ov in path.resolved_overlays
        for r in ov.rates
            ss = range(r.s_eff_start, r.s_eff_end; length = n_quad + 1)
            h  = (r.s_eff_end - r.s_eff_start) / n_quad
            total += h * sum(r.rate(s) for s in ss)
        end
    end
    return total
end

"""
    writhe(path; n) → Float64

Writhe of the path: the double integral
    Wr = (1/4π) ∫∫ (r'(s₁)×r'(s₂)) · (r(s₁)−r(s₂)) / |r(s₁)−r(s₂)|³ ds₁ ds₂

Computed numerically by sampling the path at `n` points.

# TODO
Implement Berry phase polarization contribution from writhe.  For a closed or
nearly-closed fiber loop the geometric (Berry) phase accumulated by the
polarization state equals π·Wr (in rad), where Wr is the writhe of the fiber
centerline.  See: Berry (1984) Proc. R. Soc. A 392, and Ross (1984) for the
fiber optics context.  This contribution should be added to the output of the
polarization propagator in path-integral.jl whenever the fiber forms a loop.
"""
function writhe(path::Path; n::Int = 256)
    ss = collect(range(path.s_start, path.s_end; length = n))
    rs = [position(path, s) for s in ss]
    ts = [tangent(path, s)  for s in ss]
    ds = (path.s_end - path.s_start) / (n - 1)

    Wr = 0.0
    for i in 1:n, j in 1:n
        i == j && continue
        r_ij = rs[i] - rs[j]
        d = norm(r_ij)
        d < 1e-14 * (path.s_end - path.s_start) && continue
        Wr += dot(cross(ts[i], ts[j]), r_ij) / d^3
    end
    return Wr * ds^2 / (4π)
end

# -----------------------------------------------------------------------
# Sampling
# -----------------------------------------------------------------------

function sample(path::Path, s_values)
    return [frame(path, s) for s in s_values]
end

function sample_uniform(path::Path; n::Int = 256)
    ss = range(path.s_start, path.s_end; length = n)
    return sample(path, ss)
end

# -----------------------------------------------------------------------
# JSON I/O
# -----------------------------------------------------------------------

# JSON segment format (one entry per segment in order):
#
#   {"type": "straight",  "length": 1.0, "shrinkage": 1.0}
#   {"type": "bend",      "radius": 0.05, "angle": 1.5708,
#                         "axis_angle": 0.0, "shrinkage": 1.0}
#   {"type": "catenary",  "a": 0.1, "length": 0.5,
#                         "axis_angle": 0.0, "shrinkage": 1.0}
#   {"type": "helix",     "radius": 0.03, "pitch": 0.01, "turns": 3.0,
#                         "axis_angle": 0.0, "shrinkage": 1.0}
#   {"type": "jumpby",    "delta": [0.5, 0.0, 0.3],
#                         "tangent": null, "shrinkage": 1.0}
#   {"type": "jumpto",    "destination": [1.0, 0.0, 0.5],
#                         "tangent": [0.0, 1.0, 0.0], "shrinkage": 1.0}
#
# Twist overlays follow the segments:
#
#   {"type": "twist",  "s_start": 0.5, "length": 1.0, "turns": 3}
#
# Full file structure:
#   { "path": [ ...segment dicts... ] }

"""
    path_from_json(filename) → PathSpec

Load a PathSpec from a JSON file.  Requires the JSON3 package.

TODO: add JSON3 to Project.toml and implement this function.
"""
function path_from_json(::AbstractString)
    error("path_from_json: JSON3 not yet wired; add JSON3 to Project.toml")
end

"""
    path_to_json(spec, filename)

Save a PathSpec to a JSON file.  Requires the JSON3 package.

TODO: add JSON3 to Project.toml and implement this function.
"""
function path_to_json(::PathSpec, ::AbstractString)
    error("path_to_json: JSON3 not yet wired; add JSON3 to Project.toml")
end
