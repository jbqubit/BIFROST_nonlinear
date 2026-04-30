"""
path-geometry.jl

Three-dimensional path geometry for smooth space curves.


# Three layers

- `SubpathBuilder` (mutable) — authoring target for the bang-DSL below.
- `Subpath` (immutable) — frozen snapshot: `segments` + origin state + optional terminal connector.
- `SubpathCached` (immutable) — derived layout: `subpath` + `placed_segments` + `s_end`.

`build(builder_or_subpath) → SubpathCached` runs the placement loop. Queries
(`position`, `tangent`, etc.) take a `SubpathCached`.

A `PathCached` holds multiple `SubpathCached` values concatenated end-to-end.

# Assembly

Build a subpath by calling functions on a `SubpathBuilder` in lifecycle order:

    spec = SubpathBuilder()
    origin!(spec; point=(0,0,0), outgoing_tangent=(0,0,1))   # optional; defaults shown
    straight!(spec; length=1.0)
    bend!(spec; radius=0.05, angle=π/2, axis_angle=0.0)
    helix!(spec; radius=0.03, pitch=0.01, turns=2.0)
    catenary!(spec; a=0.1, length=0.05)
    jumpby!(spec; delta=(0.0, 0.0, 0.5))                     # interior jump (not terminal)
    jumpto!(spec; destination=(1.0, 0.0, 0.5),               # seals terminal connector
            tangent=(0.0, 1.0, 0.0))
    path = build(spec)                                        # → SubpathCached

Call order rules:

- `origin!` must come before any segments; omitting it defaults the start to
  the origin pointing along +z.
- `straight!`, `bend!`, `helix!`, `catenary!`, and `jumpby!` add interior
  segments in sliding-frame order: each segment starts where the previous
  one ended.
- `jumpto!` seals the builder with a G2-continuous quintic connector to an
  absolute lab-frame destination. No further segments may be added afterward.
- `build(spec)` compiles the builder to an immutable `SubpathCached`.

To combine multiple subpaths into a single `PathCached` with a shared global
arc-length coordinate:

    sp1 = SubpathBuilder()
    straight!(sp1; length=1.0)
    jumpto!(sp1; destination=(0.0, 0.0, 2.0))

    sp2 = SubpathBuilder()
    origin!(sp2; point=(0.0, 0.0, 2.0), outgoing_tangent=(0.0, 0.0, 1.0))
    bend!(sp2; radius=0.05, angle=π/2)

    path = build([freeze(sp1), freeze(sp2)])   # → PathCached

`freeze(builder) → Subpath` snapshots the builder without compiling geometry;
`build(subpaths::Vector{Subpath}) → PathCached` compiles all subpaths and
assembles global arc-length offsets.

# Material twist

Material twist is attached as per-segment meta via `Twist <: AbstractMeta`.
A `Twist` placed in a segment's `meta` vector starts a twist run at that
segment's start and continues until the next `Twist`-bearing segment, or
until the path ends.

    spec = SubpathBuilder()
    straight!(spec; length = 1.0,
              meta = [Twist(; rate = 2π, phi_0 = 0.0)])           # constant rate
    straight!(spec; length = 1.0,
              meta = [Twist(; rate = s -> sin(s), is_continuous = true)])
    build(spec)

`rate` may be a `Real` (constant rad/m) or a `Function` `rate(s_local)` of
run-local arc length (`s_local = 0` at the start of the run). With
`is_continuous = true` the resolver computes `phi_0` from the prior run's
accumulated phase.


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
    total_material_twist(path; s_start, s_end, rtol, atol)
    total_frame_rotation(path; s_start, s_end, rtol, atol)
    writhe(path)
    sample(path, s_values)
    sample_uniform(path; n)
"""

using LinearAlgebra
using QuadGK

# -----------------------------------------------------------------------
# Abstract segment type
# -----------------------------------------------------------------------

abstract type AbstractPathSegment end

# Required interface for each concrete segment (local arc-length s ∈ [0, arc_length(seg)]).
# Return element type T follows the segment's own type parameter (Float64 for
# deterministic segments, Particles for MCM-valued fields):
#   arc_length(seg)                   → T
#   curvature(seg, s)                 → T          (κ, 1/m)
#   geometric_torsion(seg, s)         → T          (τ_geom, rad/m)
#   position_local(seg, s)            → Vector{T} length 3
#   tangent_local(seg, s)             → Vector{T} unit, length 3
#   normal_local(seg, s)              → Vector{T} unit, length 3
#   binormal_local(seg, s)            → Vector{T} unit, length 3
#   end_position_local(seg)           → Vector{T} length 3
#   end_frame_local(seg)              → (T, N, B) each Vector{T} length 3
# JumpBy resolves at build() time into a parametric QuinticConnector{T}
# that supports MCM Particles via nominalized-scalar branching.
# jumpto! on a SubpathBuilder creates a terminal connector resolved in build().

# -----------------------------------------------------------------------
# AbstractMeta — per-segment annotations
# -----------------------------------------------------------------------
# Every segment carries a `meta::Vector{AbstractMeta}` bag. path-geometry.jl
# is deliberately ignorant of what's inside; downstream layers define their
# own `AbstractMeta` subtypes (see path-geometry-meta.jl) and decide how to act
# on them.

abstract type AbstractMeta end

segment_meta(seg::AbstractPathSegment) =
    :meta ∈ fieldnames(typeof(seg)) ? seg.meta : AbstractMeta[]

# Hook for checking MCM perturbation meta.  Returns false by default; path-geometry-meta.jl
# adds true-returning methods for MCMadd and MCMmul so the check works without coupling
# path-geometry.jl to those concrete types.
_has_mcm_perturbation(::AbstractMeta) = false

# -----------------------------------------------------------------------
# Twist  (per-segment material-twist annotation)
# -----------------------------------------------------------------------

"""
    Twist(; rate, phi_0 = 0.0, is_continuous = false)

Material-twist annotation attached to a segment's `meta` vector. The twist run
begins at the host segment's start (effective arc length `s_offset_eff`) and
extends until the next segment carrying a `Twist`, or until the path ends.

- `rate` may be a `Real` (constant rad/m) or a `Function` `rate(s_local)` of
  run-local arc length where `s_local = 0` at the start of the twist run.
- `phi_0` is the absolute initial phase (rad) at the start of the run when
  `is_continuous = false`.
- `is_continuous = true` means the resolver computes `phi_0` from the prior
  twist run's accumulated phase. In that case `phi_0` must be left at its
  default `0.0`.

At most one `Twist` may appear in any single segment's `meta`. The first
`Twist` on a path must have `is_continuous = false`.
"""
struct Twist <: AbstractMeta
    rate::Union{Float64, Function}
    phi_0::Float64
    is_continuous::Bool

    function Twist(; rate, phi_0 = 0.0, is_continuous::Bool = false)
        if is_continuous && phi_0 != 0.0
            throw(ArgumentError(
                "Twist: do not specify phi_0 when is_continuous=true (phase is carried over from prior run)"))
        end
        r = rate isa Function ? rate : Float64(rate)
        new(r, Float64(phi_0), is_continuous)
    end
end

"""
    ResolvedTwistRate(s_eff_start, s_eff_end, rate, phi_0)

A single twist run resolved to absolute path coordinates. `rate` is called
(when a `Function`) with run-local arc length `s_local = s - s_eff_start`.
`phi_0` is the absolute phase (rad) at `s_eff_start`.
"""
struct ResolvedTwistRate
    s_eff_start::Float64
    s_eff_end::Float64
    rate::Union{Float64, Function}
    phi_0::Float64
end

# -----------------------------------------------------------------------
# Quadrature helper
# -----------------------------------------------------------------------
# Integrate a twist rate over a run-local interval. Constant rates take the
# analytic branch; function rates use QuadGK adaptive Gauss–Kronrod, which
# subdivides automatically for oscillatory integrands.

_integrate_rate(rate::Float64, a::Float64, b::Float64;
                rtol::Float64 = 1e-8, atol::Float64 = 0.0) = rate * (b - a)

function _integrate_rate(rate::Function, a::Float64, b::Float64;
                         rtol::Float64 = 1e-8, atol::Float64 = 0.0)
    val, _err = QuadGK.quadgk(rate, a, b; rtol = rtol, atol = atol)
    return val
end

# -----------------------------------------------------------------------
# StraightSegment
# -----------------------------------------------------------------------

struct StraightSegment{T} <: AbstractPathSegment
    length::T
    meta::Vector{AbstractMeta}
end

function StraightSegment(length; meta = AbstractMeta[])
    StraightSegment{typeof(length)}(length, Vector{AbstractMeta}(meta))
end

# A StraightSegment with negative `length` is treated as walking backward
# along the local tangent: `arc_length` is `|length|`, position advances as
# `sign(length)·s`, and the end-frame tangent/normal carry the same sign so
# downstream segments inherit a consistent, right-handed frame matching the
# direction of motion. Using `sign` (not a conditional) keeps this branch-
# free and MCM/Particles-compatible — `sign` broadcasts elementwise over a
# Particles length.

arc_length(seg::StraightSegment)         = abs(seg.length)
curvature(seg::StraightSegment, _)       = zero(seg.length)
geometric_torsion(seg::StraightSegment, _) = zero(seg.length)

position_local(seg::StraightSegment, s)   = [zero(s), zero(s), sign(seg.length) * s]
tangent_local(seg::StraightSegment, _)    = [zero(seg.length), zero(seg.length), sign(seg.length)]
normal_local(seg::StraightSegment, _)     = [sign(seg.length), zero(seg.length), zero(seg.length)]
binormal_local(seg::StraightSegment, _)   = [zero(seg.length), one(seg.length), zero(seg.length)]
end_position_local(seg::StraightSegment)  = [zero(seg.length), zero(seg.length), seg.length]
function end_frame_local(seg::StraightSegment)
    sgn = sign(seg.length)
    T = [zero(seg.length), zero(seg.length), sgn]
    N = [sgn, zero(seg.length), zero(seg.length)]
    B = [zero(seg.length), one(seg.length), zero(seg.length)]
    return (T, N, B)
end

# -----------------------------------------------------------------------
# BendSegment  (circular arc)
# -----------------------------------------------------------------------

"""
    BendSegment(radius, angle, axis_angle)

Circular arc of radius `radius` (m) sweeping `angle` (rad) in the plane whose
inward normal is at `axis_angle` (rad) from the local N-axis.

In the local frame (local z = incoming tangent, local x = incoming normal,
local y = incoming binormal), the inward normal direction is:
    n̂ = cos(axis_angle)·x̂ + sin(axis_angle)·ŷ

Curvature κ = 1 / radius.
"""
struct BendSegment{T} <: AbstractPathSegment
    radius::T
    angle::T       # total angle swept (rad)
    axis_angle::T  # orientation of inward normal in transverse plane (rad)
    meta::Vector{AbstractMeta}

    function BendSegment(radius, angle, axis_angle = 0.0;
                         meta = AbstractMeta[])
        @assert radius > 0 "BendSegment: radius must be positive"
        r, a, x = promote(radius, angle, axis_angle)
        new{typeof(r)}(r, a, x, Vector{AbstractMeta}(meta))
    end
end

arc_length(seg::BendSegment)         = seg.radius * abs(seg.angle)
curvature(seg::BendSegment, _)       = one(seg.radius) / seg.radius
geometric_torsion(seg::BendSegment, _) = zero(seg.radius)

function position_local(seg::BendSegment, s)
    R   = seg.radius
    θ   = s / R
    φ   = seg.axis_angle
    n̂   = [cos(φ), sin(φ), zero(φ)]
    return R * (1 - cos(θ)) * n̂ + [zero(R), zero(R), R * sin(θ)]
end

function tangent_local(seg::BendSegment, s)
    R = seg.radius
    θ = s / R
    φ = seg.axis_angle
    return [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)]
end

function normal_local(seg::BendSegment, s)
    R = seg.radius
    θ = s / R
    φ = seg.axis_angle
    return [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)]
end

function binormal_local(seg::BendSegment, _)
    φ = seg.axis_angle
    return [-sin(φ), cos(φ), zero(φ)]   # constant for circular arc (zero torsion)
end

function end_position_local(seg::BendSegment)
    R = seg.radius
    θ = seg.angle
    φ = seg.axis_angle
    n̂ = [cos(φ), sin(φ), zero(φ)]
    return R * (1 - cos(θ)) * n̂ + [zero(R), zero(R), R * sin(θ)]
end

function end_frame_local(seg::BendSegment)
    θ = seg.angle
    φ = seg.axis_angle
    T = [sin(θ) * cos(φ), sin(θ) * sin(φ),  cos(θ)]
    N = [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)]
    B = [-sin(φ), cos(φ), zero(φ)]
    return (T, N, B)
end

# -----------------------------------------------------------------------
# CatenarySegment
# -----------------------------------------------------------------------

"""
    CatenarySegment(a, length, axis_angle)

A catenary curve in the plane whose horizontal direction is at `axis_angle`
from the local N-axis.  `a` (m) is the catenary parameter (a = T₀/(ρg),
the ratio of horizontal tension to weight per unit length).  The curve
starts with tangent along the local z-axis (vertical catenary vertex) and
curves horizontally.  Curvature κ(s) = a / (a² + s²).
"""
struct CatenarySegment{T} <: AbstractPathSegment
    a::T
    length::T
    axis_angle::T
    meta::Vector{AbstractMeta}

    function CatenarySegment(a, length, axis_angle = 0.0;
                             meta = AbstractMeta[])
        @assert a > 0      "CatenarySegment: a must be positive"
        @assert length > 0 "CatenarySegment: length must be positive"
        av, L, x = promote(a, length, axis_angle)
        new{typeof(av)}(av, L, x, Vector{AbstractMeta}(meta))
    end
end

arc_length(seg::CatenarySegment)         = seg.length
geometric_torsion(seg::CatenarySegment, _) = zero(seg.a)

function curvature(seg::CatenarySegment, s)
    a = seg.a
    return a / (a^2 + s^2)
end

function position_local(seg::CatenarySegment, s)
    a = seg.a
    φ = seg.axis_angle
    n̂ = [cos(φ), sin(φ), zero(φ)]
    horiz = a * (sqrt(1 + (s / a)^2) - 1)
    vert  = a * asinh(s / a)
    return horiz * n̂ + [zero(a), zero(a), vert]
end

function tangent_local(seg::CatenarySegment, s)
    a = seg.a
    φ = seg.axis_angle
    q = sqrt(1 + (s / a)^2)
    return [(s / a) / q * cos(φ), (s / a) / q * sin(φ), one(q) / q]
end

function normal_local(seg::CatenarySegment, s)
    # N = dT/ds / |dT/ds|, derived analytically: N = [n̂_horiz/q, -s/a/q] normalised
    a = seg.a
    φ = seg.axis_angle
    q = sqrt(1 + (s / a)^2)
    return [cos(φ) / q, sin(φ) / q, -(s / a) / q]
end

function binormal_local(seg::CatenarySegment, s)
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
    HelixSegment(radius, pitch, turns, axis_angle)

A helix whose entry tangent is ẑ (the incoming sliding-frame tangent), ensuring
continuity with the prior segment.  `axis_angle` (rad) selects which transverse
direction n̂ = cos(axis_angle)·x̂ + sin(axis_angle)·ŷ the helix curves toward.

The helix axis â = (h·ẑ + R·n̂) / ℓ' is tilted from the transverse plane by
arctan(h/R) toward ẑ, where h = pitch/(2π) and ℓ' = √(R²+h²).  This tilt is
the geometric consequence of demanding tangent(0) = ẑ: a zero-pitch helix
reduces to a circular arc (BendSegment) in the n̂ direction.

    κ(s) = R / (R² + h²)          (constant)
    τ_geom(s) = h / (R² + h²)     (constant)
    arc_length = turns · 2π · √(R² + h²)

Local-frame basis vectors:
    n̂  = [cos(axis_angle), sin(axis_angle), 0]     (toward helix axis)
    r̂₀ = [-sin(axis_angle), cos(axis_angle), 0]    (outward radial at s=0)
    ê_φ = (R·ẑ - h·n̂) / ℓ'                        (tangential at s=0, ⊥ axis)
"""
struct HelixSegment{T} <: AbstractPathSegment
    radius::T
    pitch::T
    turns::T
    axis_angle::T
    meta::Vector{AbstractMeta}

    function HelixSegment(radius, pitch, turns,
                          axis_angle = 0.0; meta = AbstractMeta[])
        @assert radius > 0 "HelixSegment: radius must be positive"
        @assert turns  > 0 "HelixSegment: turns must be positive"
        r, p, n, x = promote(radius, pitch, turns, axis_angle)
        new{typeof(r)}(r, p, n, x, Vector{AbstractMeta}(meta))
    end
end

function _helix_h(seg::HelixSegment)
    seg.pitch / (2π)          # axial advance per radian
end

function arc_length(seg::HelixSegment)
    h = _helix_h(seg)
    return seg.turns * 2π * sqrt(seg.radius^2 + h^2)
end

function curvature(seg::HelixSegment, _)
    R = seg.radius
    h = _helix_h(seg)
    return R / (R^2 + h^2)
end

function geometric_torsion(seg::HelixSegment, _)
    R = seg.radius
    h = _helix_h(seg)
    return h / (R^2 + h^2)
end

function _helix_basis(seg::HelixSegment)
    R  = seg.radius
    h  = _helix_h(seg)
    φ_a = seg.axis_angle
    ℓ′  = sqrt(R^2 + h^2)
    z0 = zero(R); z1 = one(R)
    ẑ   = [z0, z0, z1]
    n̂   = [cos(φ_a), sin(φ_a), z0]
    â   = (h .* ẑ .+ R .* n̂) ./ ℓ′
    r̂₀  = [-sin(φ_a), cos(φ_a), z0]   # ẑ × n̂
    ê_φ = (R .* ẑ .- h .* n̂) ./ ℓ′  # â × r̂₀, tangential at s=0
    return R, h, ℓ′, â, r̂₀, ê_φ
end

function position_local(seg::HelixSegment, s)
    R, h, ℓ′, â, r̂₀, ê_φ = _helix_basis(seg)
    φ = s / ℓ′
    return h .* â .* φ .+ R .* (cos(φ) - 1) .* r̂₀ .+ R .* sin(φ) .* ê_φ
end

function tangent_local(seg::HelixSegment, s)
    R, h, ℓ′, â, r̂₀, ê_φ = _helix_basis(seg)
    φ = s / ℓ′
    return (h .* â .- R .* sin(φ) .* r̂₀ .+ R .* cos(φ) .* ê_φ) ./ ℓ′
end

function normal_local(seg::HelixSegment, s)
    _, _, ℓ′, _, r̂₀, ê_φ = _helix_basis(seg)
    φ = s / ℓ′
    return -(cos(φ) .* r̂₀ .+ sin(φ) .* ê_φ)
end

function binormal_local(seg::HelixSegment, s)
    R, h, ℓ′, â, r̂₀, ê_φ = _helix_basis(seg)
    φ = s / ℓ′
    # B = T × N; expanding in {â, r̂₀, ê_φ} orthonormal frame:
    # B = (R·â + h·sin φ·r̂₀ - h·cos φ·ê_φ) / ℓ′
    return (R .* â .+ h .* sin(φ) .* r̂₀ .- h .* cos(φ) .* ê_φ) ./ ℓ′
end

function end_position_local(seg::HelixSegment)
    return position_local(seg, arc_length(seg))
end

function end_frame_local(seg::HelixSegment)
    s = arc_length(seg)
    return (tangent_local(seg, s), normal_local(seg, s), binormal_local(seg, s))
end

"""
    JumpBy(delta, tangent_out, curvature_out, min_bend_radius)

Connects the current position to `current_position + delta` using a quintic
G2 Hermite connector. The incoming tangent and curvature vector are inherited
from the prior segment automatically. `tangent_out` is the desired outgoing
tangent direction (relative-frame); if nothing, falls back to the chord
direction. `curvature_out` is the desired outgoing curvature vector dT/ds in
the relative frame; defaults to zero (G1 outgoing match).

`min_bend_radius` (metres, nothing = unconstrained) sets a lower bound on the
radius of curvature of the connector. The handle scale λ is extended beyond
the chord default (with bisection refinement) to keep κ ≤ 1/R_min.
"""
struct JumpBy <: AbstractPathSegment
    delta::NTuple{3, Float64}
    tangent_out::Union{Nothing, NTuple{3, Float64}}
    curvature_out::Union{Nothing, NTuple{3, Float64}}
    min_bend_radius::Union{Nothing, Float64}
    meta::Vector{AbstractMeta}
end

function JumpBy(delta; tangent_out = nothing, curvature_out = nothing,
                min_bend_radius = nothing, meta = AbstractMeta[])
    d = (Float64(delta[1]), Float64(delta[2]), Float64(delta[3]))
    t = isnothing(tangent_out) ? nothing :
        (Float64(tangent_out[1]), Float64(tangent_out[2]), Float64(tangent_out[3]))
    k = isnothing(curvature_out) ? nothing :
        (Float64(curvature_out[1]), Float64(curvature_out[2]), Float64(curvature_out[3]))
    r = isnothing(min_bend_radius) ? nothing : Float64(min_bend_radius)
    JumpBy(d, t, k, r, Vector{AbstractMeta}(meta))
end

# JumpBy is context-dependent: geometry is resolved at build() time
# into a QuinticConnector.  Calling these methods on the raw struct is unsupported.
arc_length(::JumpBy)             = error("JumpBy: call build() to resolve jump geometry")
curvature(::JumpBy, ::Real)      = error("JumpBy: call build() to resolve jump geometry")
geometric_torsion(::JumpBy, ::Real) = 0.0
position_local(::JumpBy, ::Real) = error("JumpBy: call build() to resolve jump geometry")
tangent_local(::JumpBy, ::Real)  = error("JumpBy: call build() to resolve jump geometry")
normal_local(::JumpBy, ::Real)   = error("JumpBy: call build() to resolve jump geometry")
binormal_local(::JumpBy, ::Real) = error("JumpBy: call build() to resolve jump geometry")
end_position_local(::JumpBy)     = error("JumpBy: call build() to resolve jump geometry")
end_frame_local(::JumpBy)        = error("JumpBy: call build() to resolve jump geometry")

# -----------------------------------------------------------------------
# QuinticConnector  (resolved form of JumpBy / JumpTo)
# -----------------------------------------------------------------------

include(joinpath(@__DIR__, "path-geometry-connector.jl"))

# Resolve JumpBy / JumpTo to a QuinticConnector at build() time.
#
# `K_in_global` is the prior segment's terminal curvature vector in the
# global frame, computed by `build()`. It is rotated into the new segment's
# local frame here, mirroring the tangent-rotation pattern.

function _resolve_at_placement(seg::JumpBy, pos::AbstractVector, frame_mat::AbstractMatrix,
                                K_in_global::AbstractVector)
    p1_local  = collect(seg.delta)
    chord     = norm(p1_local)
    t_hat_out = isnothing(seg.tangent_out) ?
        (chord > 1e-15 ? p1_local ./ chord : [0.0, 0.0, 1.0]) :
        normalize(collect(seg.tangent_out))
    K0_local = frame_mat' * K_in_global
    K1_local = isnothing(seg.curvature_out) ? zeros(eltype(K0_local), 3) :
                                              collect(seg.curvature_out)
    return _build_quintic_connector(p1_local, t_hat_out, K0_local, K1_local;
                                    min_bend_radius = seg.min_bend_radius,
                                    meta = seg.meta)
end

# Default fallback for non-Jump segments (ignores K_in_global).
_resolve_at_placement(seg::AbstractPathSegment, ::AbstractVector, ::AbstractMatrix,
                      ::AbstractVector) = seg

# -----------------------------------------------------------------------
# PlacedSegment  (derived layout — lives inside PathSpecCached)
# -----------------------------------------------------------------------

struct PlacedSegment
    segment::AbstractPathSegment
    s_offset_eff::Real          # cumulative arc-length at segment start
    origin::AbstractVector      # global start position (length 3)
    frame::AbstractMatrix       # 3×3, columns [N_global | B_global | T_global]
                                # transforms local vectors → global: v_g = frame * v_l
end

# -----------------------------------------------------------------------
# _check_jumpto_meta
# -----------------------------------------------------------------------
# Validates that no MCMadd/MCMmul entries appear in the terminal connector
# Default no-op. Overridden in a later phase by path-geometry-meta.jl to
# prohibit MCMadd/MCMmul on the terminal connector meta.
_check_jumpto_meta(::AbstractVector{<:AbstractMeta}) = nothing

# -----------------------------------------------------------------------
# _initial_frame — Gram-Schmidt from an outgoing tangent
# -----------------------------------------------------------------------
function _initial_frame(T::AbstractVector)
    T_n = _safe_normalize(T)
    ref = abs(dot(T_n, [1.0, 0.0, 0.0])) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
    N_raw = ref .- dot(ref, T_n) .* T_n
    N_n = _safe_normalize(N_raw)
    B_n = cross(T_n, N_n)
    return T_n, N_n, B_n
end

# -----------------------------------------------------------------------
# SubpathBuilder / Subpath / SubpathCached / PathCached
# -----------------------------------------------------------------------

mutable struct SubpathBuilder
    meta::Vector{AbstractMeta}
    origin_point::Union{Nothing, NTuple{3, Float64}}
    origin_outgoing_tangent::Union{Nothing, NTuple{3, Float64}}
    origin_outgoing_curvature::Union{Nothing, NTuple{3, Float64}}
    segments::Vector{AbstractPathSegment}
    jumpto_point::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_tangent::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_curvature::Union{Nothing, NTuple{3, Float64}}
    jumpto_min_bend_radius::Union{Nothing, Float64}
    jumpto_conserve_path_length::Bool
    jumpto_meta::Vector{AbstractMeta}

    SubpathBuilder(; meta::AbstractVector{<:AbstractMeta} = AbstractMeta[]) =
        new(Vector{AbstractMeta}(meta), nothing, nothing, nothing,
            AbstractPathSegment[], nothing, nothing, nothing,
            nothing, false, AbstractMeta[])
end

struct Subpath
    meta::Vector{AbstractMeta}
    origin_point::NTuple{3, Float64}
    origin_outgoing_tangent::NTuple{3, Float64}
    origin_outgoing_curvature::NTuple{3, Float64}
    segments::Vector{AbstractPathSegment}
    jumpto_point::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_tangent::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_curvature::Union{Nothing, NTuple{3, Float64}}
    jumpto_min_bend_radius::Union{Nothing, Float64}
    jumpto_conserve_path_length::Bool
    jumpto_meta::Vector{AbstractMeta}

    function Subpath(b::SubpathBuilder)
        _check_jumpto_meta(b.jumpto_meta)
        op = isnothing(b.origin_point) ? (0.0, 0.0, 0.0) : b.origin_point
        ot = isnothing(b.origin_outgoing_tangent) ? (0.0, 0.0, 1.0) : b.origin_outgoing_tangent
        ok = isnothing(b.origin_outgoing_curvature) ? (0.0, 0.0, 0.0) : b.origin_outgoing_curvature
        new(b.meta, op, ot, ok,
            deepcopy(b.segments),
            b.jumpto_point, b.jumpto_incoming_tangent, b.jumpto_incoming_curvature,
            b.jumpto_min_bend_radius, b.jumpto_conserve_path_length, b.jumpto_meta)
    end
end

struct SubpathCached
    subpath::Subpath
    placed_segments::Vector{PlacedSegment}
    s_end::Real
    resolved_twists::Vector{ResolvedTwistRate}
end

struct PathCached
    meta::Vector{AbstractMeta}
    subpaths::Vector{SubpathCached}
    s_offsets::Vector{Float64}
    s_end::Float64
end

# Convenience: s_start is always 0 for a SubpathCached (local arc-length domain).
s_start(::SubpathCached) = 0.0
s_start(::PathCached) = 0.0

"""
    freeze(builder::SubpathBuilder) → Subpath

Snapshot a builder into an immutable `Subpath`.
"""
freeze(b::SubpathBuilder) = Subpath(b)

"""
    origin!(builder; point=(0,0,0), outgoing_tangent=(0,0,1), outgoing_curvature=(0,0,0))

Seal the start state of `builder`. Must be called before any segments are added.
Calling it after segments have been added, or calling it twice, throws `ArgumentError`.
"""
function origin!(builder::SubpathBuilder;
                 point = (0.0, 0.0, 0.0),
                 outgoing_tangent = (0.0, 0.0, 1.0),
                 outgoing_curvature = (0.0, 0.0, 0.0))
    !isnothing(builder.origin_point) && throw(ArgumentError(
        "origin!: already called on this builder"))
    !isempty(builder.segments) && throw(ArgumentError(
        "origin!: must be called before any segments are added"))
    builder.origin_point =
        (Float64(point[1]), Float64(point[2]), Float64(point[3]))
    builder.origin_outgoing_tangent =
        (Float64(outgoing_tangent[1]), Float64(outgoing_tangent[2]),
         Float64(outgoing_tangent[3]))
    builder.origin_outgoing_curvature =
        (Float64(outgoing_curvature[1]), Float64(outgoing_curvature[2]),
         Float64(outgoing_curvature[3]))
    return builder
end

function straight!(spec::SubpathBuilder; length,
                   meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    !isnothing(spec.jumpto_point) && throw(ArgumentError(
        "straight!: cannot add segments after jumpto! has been called"))
    push!(spec.segments, StraightSegment(length; meta))
    return spec
end

function bend!(spec::SubpathBuilder; radius::Real, angle::Real, axis_angle::Real = 0.0,
               meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    !isnothing(spec.jumpto_point) && throw(ArgumentError(
        "bend!: cannot add segments after jumpto! has been called"))
    push!(spec.segments, BendSegment(radius, angle, axis_angle; meta))
    return spec
end

function helix!(spec::SubpathBuilder; radius::Real, pitch::Real, turns::Real,
                axis_angle::Real = 0.0, meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    !isnothing(spec.jumpto_point) && throw(ArgumentError(
        "helix!: cannot add segments after jumpto! has been called"))
    push!(spec.segments, HelixSegment(radius, pitch, turns, axis_angle; meta))
    return spec
end

function catenary!(spec::SubpathBuilder; a::Real, length::Real, axis_angle::Real = 0.0,
                   meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    !isnothing(spec.jumpto_point) && throw(ArgumentError(
        "catenary!: cannot add segments after jumpto! has been called"))
    push!(spec.segments, CatenarySegment(a, length, axis_angle; meta))
    return spec
end

function jumpby!(spec::SubpathBuilder; delta, tangent = nothing,
                 curvature_out = nothing, min_bend_radius = nothing,
                 meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    !isnothing(spec.jumpto_point) && throw(ArgumentError(
        "jumpby!: cannot add segments after jumpto! has been called"))
    push!(spec.segments, JumpBy(delta; tangent_out = tangent,
                                curvature_out = curvature_out,
                                min_bend_radius, meta))
    return spec
end

"""
    jumpto!(builder; jumpto_point, ...)

Seal the terminal connector of `builder`. Calling it twice throws.
Adding segments after this call throws.

Accepts `destination` as a backward-compatible alias for `jumpto_point`,
and `tangent`/`curvature_out` as aliases for `jumpto_incoming_tangent`/
`jumpto_incoming_curvature`.
"""
function jumpto!(builder::SubpathBuilder;
                 jumpto_point = nothing, destination = nothing,
                 jumpto_incoming_tangent = nothing, tangent = nothing,
                 jumpto_incoming_curvature = nothing, curvature_out = nothing,
                 jumpto_min_bend_radius = nothing, min_bend_radius = nothing,
                 jumpto_conserve_path_length::Bool = false,
                 meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    !isnothing(builder.jumpto_point) && throw(ArgumentError(
        "jumpto!: terminal connector already set (called twice)"))
    _check_jumpto_meta(meta)
    point = !isnothing(jumpto_point) ? jumpto_point : destination
    isnothing(point) && throw(ArgumentError(
        "jumpto!: jumpto_point (or destination) is required"))
    incoming_tangent   = !isnothing(jumpto_incoming_tangent) ? jumpto_incoming_tangent : tangent
    incoming_curvature = !isnothing(jumpto_incoming_curvature) ? jumpto_incoming_curvature :
                         curvature_out
    R_min = !isnothing(jumpto_min_bend_radius) ? jumpto_min_bend_radius : min_bend_radius
    builder.jumpto_point = (Float64(point[1]), Float64(point[2]), Float64(point[3]))
    if !isnothing(incoming_tangent)
        t = incoming_tangent
        builder.jumpto_incoming_tangent = (Float64(t[1]), Float64(t[2]), Float64(t[3]))
    end
    if !isnothing(incoming_curvature)
        k = incoming_curvature
        builder.jumpto_incoming_curvature = (Float64(k[1]), Float64(k[2]), Float64(k[3]))
    end
    builder.jumpto_min_bend_radius = isnothing(R_min) ? nothing : Float64(R_min)
    builder.jumpto_conserve_path_length = jumpto_conserve_path_length
    builder.jumpto_meta = Vector{AbstractMeta}(meta)
    return builder
end

# -----------------------------------------------------------------------
# build()
# -----------------------------------------------------------------------

"""
    build(builder_or_subpath) → SubpathCached

Compile a `SubpathBuilder` or `Subpath` into an immutable `SubpathCached`.
"""
_safe_normalize(v::AbstractVector) = v ./ sqrt(sum(abs2, v))

build(b::SubpathBuilder) = build(Subpath(b))

function build(sub::Subpath)
    pos         = collect(sub.origin_point)
    T_frame, N_frame, B_frame = _initial_frame(collect(sub.origin_outgoing_tangent))
    K_in_global = collect(sub.origin_outgoing_curvature)
    s_eff       = 0.0
    placed      = PlacedSegment[]

    isempty(sub.segments) && isnothing(sub.jumpto_point) &&
        error("build: Subpath contains no segments and no jumpto! connector")

    for seg_orig in sub.segments
        frame      = hcat(N_frame, B_frame, T_frame)   # columns: [N | B | T]
        seg_placed = _resolve_at_placement(seg_orig, pos, frame, K_in_global)
        push!(placed, PlacedSegment(seg_placed, s_eff, copy(pos), copy(frame)))

        # Advance position and frame
        pos_end_local         = end_position_local(seg_placed)
        (T_end_l, N_end_l, _) = end_frame_local(seg_placed)

        pos     = pos + frame * pos_end_local
        T_frame = _safe_normalize(frame * T_end_l)
        N_end_g = frame * N_end_l
        # Re-orthogonalise to prevent floating-point drift
        N_frame = _safe_normalize(N_end_g - dot(N_end_g, T_frame) * T_frame)
        B_frame = cross(T_frame, N_frame)

        # Update K_in_global for the next segment: κ at end times the unit
        # normal at the end, both expressed globally.
        L_seg = arc_length(seg_placed)
        κ_end = curvature(seg_placed, L_seg)
        K_in_global = κ_end .* N_frame

        s_eff += L_seg
    end

    # Terminal connector (from jumpto!)
    if !isnothing(sub.jumpto_point)
        frame     = hcat(N_frame, B_frame, T_frame)
        p1_local  = frame' * (collect(sub.jumpto_point) .- pos)
        chord     = norm(p1_local)
        t_hat_out = if isnothing(sub.jumpto_incoming_tangent)
            chord > 1e-15 ? p1_local ./ chord : [0.0, 0.0, 1.0]
        else
            normalize(frame' * normalize(collect(sub.jumpto_incoming_tangent)))
        end
        K0_local = frame' * K_in_global
        K1_local = isnothing(sub.jumpto_incoming_curvature) ? zeros(3) :
                   frame' * collect(sub.jumpto_incoming_curvature)
        connector = _build_quintic_connector(p1_local, t_hat_out, K0_local, K1_local;
                        min_bend_radius = sub.jumpto_min_bend_radius,
                        meta            = sub.jumpto_meta)
        push!(placed, PlacedSegment(connector, s_eff, copy(pos), copy(frame)))
        s_eff += arc_length(connector)
    end

    # Twist runs are anchored at the *nominal* arc-length boundaries. Twist
    # data (rate, phi_0) is deterministic by design; nominalizing the anchors
    # via _qc_nominalize lets MCM-typed segments (Particles s_eff) coexist
    # with Twist meta on the same path. _qc_nominalize is identity on Real
    # and uses MonteCarloMeasurements.pmean when Particles are present.
    has_twist = any(ps -> any(m -> m isa Twist, segment_meta(ps.segment)), placed)
    resolved_twists = has_twist ?
        _resolve_twists(placed, Float64(_qc_nominalize(s_eff))) :
        ResolvedTwistRate[]
    return SubpathCached(sub, placed, s_eff, resolved_twists)
end

function build(subpaths::Vector{Subpath};
               meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    isempty(subpaths) && error("build: empty subpath vector")
    cached = [build(sp) for sp in subpaths]
    s_offsets = Vector{Float64}(undef, length(cached))
    s_acc = 0.0
    for (i, sc) in enumerate(cached)
        s_offsets[i] = s_acc
        s_acc += Float64(_qc_nominalize(sc.s_end))
    end
    return PathCached(Vector{AbstractMeta}(meta), cached, s_offsets, s_acc)
end

# -----------------------------------------------------------------------
# _resolve_twists
# -----------------------------------------------------------------------
# Walk placed segments in order, collect Twist anchors from each segment's
# meta, and emit one ResolvedTwistRate per anchor. Each run extends from the
# anchor's segment start to the next anchor's segment start (or to s_end).
# is_continuous=true anchors take their phi_0 from the prior run's accumulated
# phase.

function _resolve_twists(placed::Vector{PlacedSegment}, s_end::Float64)
    # Collect (s_eff_start, twist) anchors with at-most-one-Twist-per-segment validation.
    anchors = Tuple{Float64, Twist}[]
    for ps in placed
        twists_here = Twist[]
        for m in segment_meta(ps.segment)
            m isa Twist && push!(twists_here, m)
        end
        if length(twists_here) > 1
            throw(ArgumentError(
                "Path build: segment at s_offset_eff = $(ps.s_offset_eff) carries " *
                "$(length(twists_here)) Twist meta entries; at most one is permitted"))
        end
        if !isempty(twists_here)
            push!(anchors, (Float64(_qc_nominalize(ps.s_offset_eff)), twists_here[1]))
        end
    end

    isempty(anchors) && return ResolvedTwistRate[]

    n = length(anchors)
    out = Vector{ResolvedTwistRate}(undef, n)
    prev_phi_0 = 0.0
    prev_run_length = 0.0
    prev_rate::Union{Float64, Function} = 0.0

    for i in 1:n
        s_start, tw = anchors[i]
        s_run_end = (i < n) ? anchors[i + 1][1] : s_end

        if tw.is_continuous
            if i == 1
                throw(ArgumentError(
                    "Path build: first Twist on the path has is_continuous=true; " *
                    "the first run has no prior phase to continue from"))
            end
            phi_0 = prev_phi_0 + _integrate_rate(prev_rate, 0.0, prev_run_length)
        else
            phi_0 = tw.phi_0
        end

        out[i] = ResolvedTwistRate(s_start, s_run_end, tw.rate, phi_0)

        prev_phi_0 = phi_0
        prev_run_length = s_run_end - s_start
        prev_rate = tw.rate
    end

    return out
end

# -----------------------------------------------------------------------
# Segment lookup helpers
# -----------------------------------------------------------------------

function _find_placed_segment(path::SubpathCached, s)
    n = length(path.placed_segments)
    for i in 1:n
        ps = path.placed_segments[i]
        seg_len = arc_length(ps.segment)
        s_end = ps.s_offset_eff + seg_len
        if s <= s_end + 1e-12 || i == n
            s_local = clamp(s - ps.s_offset_eff, zero(seg_len), seg_len)
            return ps, s_local
        end
    end
    error("s = $s out of path bounds [0.0, $(path.s_end)]")
end

function _local_to_global(ps::PlacedSegment, v_local::AbstractVector)
    return ps.frame * v_local
end

# -----------------------------------------------------------------------
# Differential geometry interface on SubpathCached
# -----------------------------------------------------------------------

arc_length(path::SubpathCached) = path.s_end

function arc_length(::SubpathCached, s1, s2)
    @assert s2 >= s1 "arc_length: require s2 >= s1"
    return s2 - s1
end

function curvature(path::SubpathCached, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return curvature(ps.segment, s_local)
end

function geometric_torsion(path::SubpathCached, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return geometric_torsion(ps.segment, s_local)
end

"""
    material_twist(path, s)

Material twist rate (rad/m) at effective arc length `s`, summed over all
resolved twist runs that contain `s`. Runs are disjoint by construction; the
sum exists only as a robustness guard.
"""
function material_twist(path::SubpathCached, s)
    τ = zero(s isa AbstractFloat ? s : Float64(s))
    for r in path.resolved_twists
        if r.s_eff_start <= s <= r.s_eff_end
            τ += r.rate isa Function ? r.rate(s - r.s_eff_start) : r.rate
        end
    end
    return τ
end

function position(path::SubpathCached, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return ps.origin + _local_to_global(ps, position_local(ps.segment, s_local))
end

function tangent(path::SubpathCached, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return _local_to_global(ps, tangent_local(ps.segment, s_local))
end

function normal(path::SubpathCached, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return _local_to_global(ps, normal_local(ps.segment, s_local))
end

function binormal(path::SubpathCached, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return _local_to_global(ps, binormal_local(ps.segment, s_local))
end

function frame(path::SubpathCached, s::Real)
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

start_point(path::SubpathCached)   = position(path, 0.0)
end_point(path::SubpathCached)     = position(path, path.s_end)
start_tangent(path::SubpathCached) = tangent(path, 0.0)
end_tangent(path::SubpathCached)   = tangent(path, path.s_end)

# -----------------------------------------------------------------------
# Path measures
# -----------------------------------------------------------------------

path_length(path::SubpathCached) = arc_length(path)

function cartesian_distance(path::SubpathCached, s1::Real, s2::Real)
    return norm(position(path, s2) - position(path, s1))
end

function bounding_box(path::SubpathCached; n::Int = 512)
    s0 = 0.0
    s1 = Float64(_qc_nominalize(path.s_end))
    ss = range(s0, s1; length = n)
    pts = [position(path, s) for s in ss]
    lo = minimum(reduce(hcat, pts); dims = 2) |> vec
    hi = maximum(reduce(hcat, pts); dims = 2) |> vec
    return (; lo, hi)
end

function total_turning_angle(path::SubpathCached)
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

function total_torsion(path::SubpathCached)
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

"""
    total_material_twist(path; s_start, s_end, rtol = 1e-8, atol = 0.0) → Float64

Integrated material twist ``∫ τ_{\\mathrm{mat}}(s) \\, ds`` over effective arc length from
`s_start` to `s_end` (defaults: full path). Both endpoints must lie in `[0, path.s_end]`.

Constant twist rates integrate analytically; `Function` rates use adaptive
Gauss–Kronrod (QuadGK) at the supplied `rtol`/`atol`.
"""
function total_material_twist(
    path::SubpathCached;
    s_start::Real = 0.0,
    s_end::Real   = path.s_end,
    rtol::Real    = 1e-8,
    atol::Real    = 0.0,
)
    s_lo = Float64(_qc_nominalize(s_start))
    s_hi = Float64(_qc_nominalize(s_end))
    if s_lo > s_hi
        throw(ArgumentError(
            "total_material_twist: require s_start ≤ s_end; got s_start=$(s_lo), s_end=$(s_hi)"))
    end
    ps0 = 0.0
    ps1 = Float64(_qc_nominalize(path.s_end))
    if !(ps0 - 1e-12 <= s_lo <= ps1 + 1e-12) || !(ps0 - 1e-12 <= s_hi <= ps1 + 1e-12)
        throw(ArgumentError(
            "total_material_twist: require path.s_start ≤ s ≤ path.s_end for both endpoints; " *
            "got [$(s_lo), $(s_hi)] m vs path domain [$(ps0), $(ps1)] m"))
    end
    s_lo == s_hi && return 0.0

    total = 0.0
    rtolf = Float64(rtol)
    atolf = Float64(atol)
    for r in path.resolved_twists
        a = max(r.s_eff_start, s_lo)
        b = min(r.s_eff_end, s_hi)
        b <= a && continue
        total += _integrate_rate(r.rate, a - r.s_eff_start, b - r.s_eff_start;
                                 rtol = rtolf, atol = atolf)
    end
    return total
end

"""
    total_frame_rotation(path; s_start, s_end, rtol = 1e-8, atol = 0.0) → Float64

Total rotation of the polarization reference frame over effective arc length from `s_start`
to `s_end`, integrating both contributions:

    dψ/ds = τ_geom(s) + Ω_material(s)

where `τ_geom` is the geometric torsion of the centerline (nonzero for helices; zero for
straight segments and circular bends) and `Ω_material` is the applied material twist rate
from `Twist` meta annotations. Returns the integral in radians.

See `path-geometry.md` for a discussion of how this differs from `total_torsion` and
`total_material_twist` individually.
"""
function total_frame_rotation(
    path::SubpathCached;
    s_start::Real = 0.0,
    s_end::Real   = path.s_end,
    rtol::Real    = 1e-8,
    atol::Real    = 0.0,
)
    s_lo = Float64(_qc_nominalize(s_start))
    s_hi = Float64(_qc_nominalize(s_end))
    if s_lo > s_hi
        throw(ArgumentError(
            "total_frame_rotation: require s_start ≤ s_end; got s_start=$(s_lo), s_end=$(s_hi)"))
    end
    ps0 = 0.0
    ps1 = Float64(_qc_nominalize(path.s_end))
    if !(ps0 - 1e-12 <= s_lo <= ps1 + 1e-12) || !(ps0 - 1e-12 <= s_hi <= ps1 + 1e-12)
        throw(ArgumentError(
            "total_frame_rotation: require path.s_start ≤ s ≤ path.s_end for both endpoints; " *
            "got [$(s_lo), $(s_hi)] m vs path domain [$(ps0), $(ps1)] m"))
    end
    s_lo == s_hi && return 0.0

    # Geometric torsion: integrate segment-by-segment, with closed-form
    # shortcuts where available.
    #
    # The endpoints/clipping bounds (s_lo, s_hi, seg_s_*_nom) use nominalized
    # arc lengths — they're deterministic interval delimiters. The
    # *integration weight* (b_loc - a_loc) is computed from the segment's
    # full Particles arc_length so length-uncertainty propagates through the
    # τ_geom contribution. For full-path queries frac_a=0, frac_b=1 so
    # b_loc - a_loc == L_seg exactly.
    #
    # The QuadGK fallback branch passes Particles bounds when the segment is
    # MCM-typed; this requires `MonteCarloMeasurements.unsafe_comparisons(true)`
    # in the calling scope. Currently only `QuinticConnector` falls through
    # here (Helix/Straight/Bend/Catenary all hit closed-form branches).
    τ_total = 0.0
    rtolf = Float64(rtol)
    atolf = Float64(atol)
    for ps in path.placed_segments
        seg = ps.segment
        seg_s_start_nom = Float64(_qc_nominalize(ps.s_offset_eff))
        seg_s_end_nom   = seg_s_start_nom +
                          Float64(_qc_nominalize(arc_length(seg)))
        a_nom = max(seg_s_start_nom, s_lo)
        b_nom = min(seg_s_end_nom,   s_hi)
        b_nom <= a_nom && continue

        L_seg  = arc_length(seg)
        seg_span_nom = seg_s_end_nom - seg_s_start_nom
        # Avoid 0/0 for zero-length segments (skipped above by b <= a).
        frac_a = (a_nom - seg_s_start_nom) / seg_span_nom
        frac_b = (b_nom - seg_s_start_nom) / seg_span_nom
        a_loc = frac_a * L_seg
        b_loc = frac_b * L_seg

        if seg isa StraightSegment || seg isa BendSegment || seg isa CatenarySegment
            # τ_geom = 0
        elseif seg isa HelixSegment
            τ_total += geometric_torsion(seg, zero(L_seg)) * (b_loc - a_loc)
        else
            val, _ = QuadGK.quadgk(s -> geometric_torsion(seg, s), a_loc, b_loc;
                                   rtol = rtolf, atol = atolf)
            τ_total += val
        end
    end

    Ω_total = total_material_twist(path; s_start = s_lo, s_end = s_hi,
                                   rtol = rtolf, atol = atolf)
    return τ_total + Ω_total
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
function writhe(path::SubpathCached; n::Int = 256)
    s0 = 0.0
    s1 = Float64(_qc_nominalize(path.s_end))
    ss = collect(range(s0, s1; length = n))
    rs = [position(path, s) for s in ss]
    ts = [tangent(path, s)  for s in ss]
    ds = (s1 - s0) / (n - 1)

    Wr = 0.0
    for i in 1:n, j in 1:n
        i == j && continue
        r_ij = rs[i] - rs[j]
        d = norm(r_ij)
        d < 1e-14 * (s1 - s0) && continue
        Wr += dot(cross(ts[i], ts[j]), r_ij) / d^3
    end
    return Wr * ds^2 / (4π)
end

# -----------------------------------------------------------------------
# Sampling
# -----------------------------------------------------------------------

"""
    Sample

One evaluated point on a `Path`: arc-length coordinate `s` plus all Frenet frame quantities.

Fields: `s`, `position`, `tangent`, `normal`, `binormal`, `curvature`,
`geometric_torsion`, `material_twist`.
"""
struct Sample
    s                 :: Real
    position          :: AbstractVector
    tangent           :: AbstractVector
    normal            :: AbstractVector
    binormal          :: AbstractVector
    curvature         :: Real
    geometric_torsion :: Real
    material_twist    :: Real
end

"""
    PathSample

Dense samples of a `Path` over `[s_start, s_end]`.

Fields:
- `samples :: Vector{Sample}` — one entry per arc-length station.
- `s_start`, `s_end` — the arc-length interval that was sampled.
- `n` — number of sample points.
"""
struct PathSample
    samples :: Vector{Sample}
    s_start :: Float64
    s_end   :: Float64
    n       :: Int
end

"""
    _segment_total_angle(seg) → Float64

Total angle (rad) swept by the tangent vector along `seg`. This is the geometric
quantity that governs how many sample points a segment needs.

- `BendSegment`:     the literal swept arc angle.
- `CatenarySegment`: `arctan(L_eff / a_eff)`, the angle from the vertical entry
  tangent to the exit tangent.
- `HelixSegment`:    `turns · 2π` — each full turn is 2π of tangent rotation.
- `StraightSegment` and connectors: 0 (no tangent rotation).
"""
function _segment_total_angle(seg::BendSegment)
    return abs(seg.angle)
end

function _segment_total_angle(seg::CatenarySegment)
    return atan(arc_length(seg) / seg.a)
end

function _segment_total_angle(seg::HelixSegment)
    return seg.turns * 2π
end

function _segment_total_angle(seg::QuinticConnector)
    L = arc_length(seg)
    L < 1e-15 && return 0.0
    n = 16
    ss = range(0.0, L; length = n + 1)
    h  = L / n
    return h * (sum(curvature(seg, s) for s in ss) - curvature(seg, 0.0)/2 - curvature(seg, L)/2)
end

_segment_total_angle(::AbstractPathSegment) = 0.0

_budget_scalar(x::AbstractFloat) = Float64(x)
_budget_scalar(x::Integer) = Float64(x)
function _budget_scalar(x)
    if hasfield(typeof(x), :particles)
        return Float64(maximum(getfield(x, :particles)))
    end
    return Float64(x)
end

"""
    _segment_point_budget(ps, path, s_lo, s_hi, fidelity) → Int

Adaptive point budget for the portion of `PlacedSegment` `ps` that falls within
`[s_lo, s_hi]`.

Geometric rule (applied to the clipped fraction of the segment):
  `ceil(fidelity · total_angle / (2π) · 32)`, minimum 2.

Twist rule: compute `|∫ τ_mat ds|` over the same interval using
`total_material_twist` at `n_quad = 32`, then apply the same formula with
that integrated angle.

The point budget is the maximum of the two rules.
"""
function _segment_point_budget(
    ps::PlacedSegment,
    path,
    s_lo::Float64,
    s_hi::Float64,
    fidelity::Float64,
)
    seg = ps.segment
    seg_s_start = ps.s_offset_eff
    seg_s_end   = ps.s_offset_eff + arc_length(seg)

    # Clipped interval for this segment
    a = max(s_lo, seg_s_start)
    b = min(s_hi, seg_s_end)
    b <= a && return 2

    # Fraction of the segment that is sampled
    seg_len = seg_s_end - seg_s_start
    frac = seg_len > 0.0 ? (b - a) / seg_len : 1.0

    # Geometric budget
    geom_angle  = _budget_scalar(_segment_total_angle(seg) * frac)
    geom_budget = max(2, ceil(Int, fidelity * geom_angle / (2π) * 32))

    # Twist budget: integrate |τ_mat| over [a, b] at loose tolerance — this is
    # only sizing a sampling budget, so 1e-3 relative is plenty.
    twist_total  = total_material_twist(path; s_start = a, s_end = b, rtol = 1e-3)
    twist_angle  = abs(_budget_scalar(twist_total))
    twist_budget = max(2, ceil(Int, fidelity * twist_angle / (2π) * 32))

    return max(geom_budget, twist_budget)
end

"""
    sample_path(path, s1, s2; fidelity = 1.0) → PathSample

Sample `path` over `[s1, s2]` with an adaptive point distribution driven by `fidelity`.

Points are allocated per segment using the heuristic:

    n_seg = max(geometric_budget, twist_budget)

where each budget is `ceil(fidelity · Δφ / (2π) · 32)` (minimum 2), with:

- **geometric**: `Δφ` = total tangent-rotation angle of the segment (bend angle,
  catenary turning angle, or `turns·2π` for a helix; 0 for straight segments).
- **twist**: `|∫ τ_mat ds|` over the segment's portion of `[s1, s2]`, computed
  via exact quadrature.

`fidelity = 1.0` gives ~32 points per full 2π of rotation or twist. Straight
segments with no twist get 2 points (endpoints only). Junction points between
segments are always included and shared (no duplicates).
"""
function sample_path(path::SubpathCached, s1::Real, s2::Real;
                     fidelity::Float64 = 1.0)
    @assert s2 > s1    "sample_path: require s2 > s1"
    @assert fidelity > 0.0 "sample_path: fidelity must be positive"

    s_lo = Float64(_qc_nominalize(s1))
    s_hi = Float64(_qc_nominalize(s2))

    # Collect arc-length stations segment by segment, sharing junction points.
    # Segment edges are deterministic — nominalize Particles cumulative offsets
    # so range() and clipping comparisons stay Float64.
    all_s = Float64[]
    for ps in path.placed_segments
        seg_s_start = Float64(_qc_nominalize(ps.s_offset_eff))
        seg_s_end   = seg_s_start + Float64(_qc_nominalize(arc_length(ps.segment)))

        a = max(s_lo, seg_s_start)
        b = min(s_hi, seg_s_end)
        b <= a && continue

        n_seg = _segment_point_budget(ps, path, s_lo, s_hi, fidelity)
        seg_ss = collect(range(a, b; length = n_seg))

        if isempty(all_s)
            append!(all_s, seg_ss)
        else
            # Drop the first point if it duplicates the last accumulated point.
            start_idx = (seg_ss[1] ≈ all_s[end]) ? 2 : 1
            append!(all_s, @view seg_ss[start_idx:end])
        end
    end

    # Guarantee at least the two endpoints even for degenerate paths.
    if isempty(all_s)
        all_s = [s_lo, s_hi]
    elseif length(all_s) == 1
        push!(all_s, s_hi)
    end

    n = length(all_s)
    samples = Vector{Sample}(undef, n)
    for i in eachindex(all_s)
        fr = frame(path, all_s[i])
        samples[i] = Sample(
            all_s[i],
            fr.position,
            fr.tangent,
            fr.normal,
            fr.binormal,
            fr.curvature,
            fr.geometric_torsion,
            fr.material_twist,
        )
    end
    return PathSample(samples, s_lo, s_hi, n)
end

function normalize_breakpoints(breakpoints::AbstractVector{<:Real})
    return sort(unique(copy(breakpoints)))
end

function path_segment_breakpoints(path::SubpathCached)
    points = Float64[0.0]
    for ps in path.placed_segments
        push!(points, Float64(_qc_nominalize(ps.s_offset_eff)))
        push!(points, Float64(_qc_nominalize(
            ps.s_offset_eff + arc_length(ps.segment))))
    end
    push!(points, Float64(_qc_nominalize(path.s_end)))
    return normalize_breakpoints(points)
end

function path_twist_breakpoints(path::SubpathCached)
    points = Float64[0.0, Float64(_qc_nominalize(path.s_end))]
    for r in path.resolved_twists
        push!(points, r.s_eff_start)
        push!(points, r.s_eff_end)
    end
    return normalize_breakpoints(points)
end

function breakpoints(path::SubpathCached)
    return normalize_breakpoints(
        vcat(path_segment_breakpoints(path), path_twist_breakpoints(path)))
end

function sample(path::SubpathCached, s_values)
    return [frame(path, s) for s in s_values]
end

function sample_uniform(path::SubpathCached; n::Int = 256)
    ss = range(0.0, path.s_end; length = n)
    return sample(path, ss)
end

# -----------------------------------------------------------------------
# PathCached — multi-subpath container queries
# -----------------------------------------------------------------------

function _find_subpath(path::PathCached, s::Real)
    n = length(path.subpaths)
    for i in 1:n
        s_hi = path.s_offsets[i] + Float64(_qc_nominalize(path.subpaths[i].s_end))
        if s <= s_hi + 1e-12 || i == n
            s_local = clamp(s - path.s_offsets[i],
                            zero(s - path.s_offsets[i]),
                            path.subpaths[i].s_end)
            return path.subpaths[i], s_local
        end
    end
    error("s = $s out of PathCached bounds [0, $(path.s_end)]")
end

arc_length(path::PathCached) = path.s_end
arc_length(::PathCached, s1, s2) = s2 - s1

function curvature(path::PathCached, s::Real)
    sp, sl = _find_subpath(path, s); return curvature(sp, sl)
end
function geometric_torsion(path::PathCached, s::Real)
    sp, sl = _find_subpath(path, s); return geometric_torsion(sp, sl)
end
function material_twist(path::PathCached, s)
    sp, sl = _find_subpath(path, s); return material_twist(sp, sl)
end
function position(path::PathCached, s::Real)
    sp, sl = _find_subpath(path, s); return position(sp, sl)
end
function tangent(path::PathCached, s::Real)
    sp, sl = _find_subpath(path, s); return tangent(sp, sl)
end
function normal(path::PathCached, s::Real)
    sp, sl = _find_subpath(path, s); return normal(sp, sl)
end
function binormal(path::PathCached, s::Real)
    sp, sl = _find_subpath(path, s); return binormal(sp, sl)
end
function frame(path::PathCached, s::Real)
    sp, sl = _find_subpath(path, s); return frame(sp, sl)
end

start_point(path::PathCached)   = position(path, 0.0)
end_point(path::PathCached)     = position(path, path.s_end)
start_tangent(path::PathCached) = tangent(path, 0.0)
end_tangent(path::PathCached)   = tangent(path, path.s_end)

path_length(path::PathCached) = path.s_end
cartesian_distance(path::PathCached, s1::Real, s2::Real) =
    norm(position(path, s2) - position(path, s1))

function bounding_box(path::PathCached; n::Int = 512)
    ss = range(0.0, path.s_end; length = n)
    pts = [position(path, s) for s in ss]
    lo = minimum(reduce(hcat, pts); dims = 2) |> vec
    hi = maximum(reduce(hcat, pts); dims = 2) |> vec
    return (; lo, hi)
end

function path_segment_breakpoints(path::PathCached)
    points = Float64[0.0]
    for (i, sc) in enumerate(path.subpaths)
        off = path.s_offsets[i]
        for bp in path_segment_breakpoints(sc)
            push!(points, off + bp)
        end
    end
    push!(points, path.s_end)
    return normalize_breakpoints(points)
end

function path_twist_breakpoints(path::PathCached)
    points = Float64[0.0, path.s_end]
    for (i, sc) in enumerate(path.subpaths)
        off = path.s_offsets[i]
        for bp in path_twist_breakpoints(sc)
            push!(points, off + bp)
        end
    end
    return normalize_breakpoints(points)
end

function breakpoints(path::PathCached)
    return normalize_breakpoints(
        vcat(path_segment_breakpoints(path), path_twist_breakpoints(path)))
end

function sample(path::PathCached, s_values)
    return [frame(path, s) for s in s_values]
end

function sample_uniform(path::PathCached; n::Int = 256)
    ss = range(0.0, path.s_end; length = n)
    return sample(path, ss)
end

# -----------------------------------------------------------------------
# Backward-compatibility aliases (Phase 2 → Phase 4 migration shim)
# -----------------------------------------------------------------------
# These aliases preserve the pre-Phase-2 public names so that existing
# call sites (tests, demos, fiber-path.jl) continue to work without
# change during the transition to SubpathBuilder / Subpath / SubpathCached.
const PathSpecBuilder = SubpathBuilder
const PathSpecCached  = SubpathCached
