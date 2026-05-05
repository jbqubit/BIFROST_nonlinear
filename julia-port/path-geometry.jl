"""
path-geometry.jl

Three-dimensional path geometry for smooth space curves.

# Three layers

- `SubpathBuilder` (mutable) — authoring target for the bang-DSL below. A
  Subpath specification begins with `start!(...)` and ends with `jumpto!(...)`.
- `Subpath` (immutable) — frozen snapshot of user-supplied data only.
- `SubpathBuilt` (immutable) — derived layout: `subpath` + `placed_segments`
  + `jumpto_quintic_connector` + `resolved_twists`.
- `PathBuilt` (immutable) — ordered container of `SubpathBuilt`s.

`build(builder_or_subpath) → SubpathBuilt` runs the placement loop on one
Subpath. `build(::Vector{Subpath})` or `build(::Vector{SubpathBuilt})`
produces a `PathBuilt`.

# Authoring

Sliding-frame interior segments + a global terminal jumpto:

    sb = SubpathBuilder()
    start!(sb; point=(0,0,0), outgoing_tangent=(0,0,1))
    straight!(sb; length=1.0)
    bend!(sb; radius=0.05, angle=π/2)
    jumpto!(sb; point=(0.05, 0.0, 1.05))     # seals the Subpath
    sub = Subpath(sb)
    sb_built = build(sub)

Interior `jumpby!` is supported as a relative jump within a Subpath. Each
Subpath must call `start!` before any segment and `jumpto!` exactly once at
the end.

# Material twist

Material twist is attached as per-segment meta via `Twist <: AbstractMeta`.
A `Twist` placed in a segment's `meta` vector starts a twist run at that
segment's start and continues until the next `Twist`-bearing segment, or
until the Subpath ends.

`rate` may be a `Real` (constant rad/m) or a `Function` `rate(s_local)`. A
`Twist` with `is_continuous = true` placed on the first twist anchor of a
Subpath leaves the resolution pending; `build(::Vector{SubpathBuilt})` then
inherits `phi_0` from the prior Subpath's terminal phase.

# Independence

Each `Subpath` and `SubpathBuilt` is fully independent of all others. The
Subpath holds its `start_point` and `jumpto_point` as the only globally-anchored
values. `build(::Vector{SubpathBuilt})` is the first ordering-aware layer; it
checks endpoint conformity between adjacent Subpaths and resolves any
cross-Subpath twist continuity.

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
    breakpoints(path)
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
# JumpBy/JumpTo resolve at build() time into a parametric QuinticConnector{T}
# that supports MCM Particles via nominalized-scalar branching.

# -----------------------------------------------------------------------
# AbstractMeta — per-segment annotations
# -----------------------------------------------------------------------
# Every segment carries a `meta::Vector{AbstractMeta}` bag. path-geometry.jl
# is deliberately ignorant of what's inside; downstream layers define their
# own `AbstractMeta` subtypes (see fiber-path-meta.jl) and decide how to act
# on them.

abstract type AbstractMeta end

segment_meta(seg::AbstractPathSegment) =
    :meta ∈ fieldnames(typeof(seg)) ? seg.meta : AbstractMeta[]

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

# JumpBy is context-dependent: geometry is resolved at build() time into a
# QuinticConnector. Calling these methods on the raw struct is unsupported.
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
# QuinticConnector  (resolved form of JumpBy and Subpath terminal connectors)
# -----------------------------------------------------------------------

include(joinpath(@__DIR__, "path-geometry-connector.jl"))

# Concrete meta vocabulary lives in path-geometry-meta.jl. It is part of the
# geometry layer (Twist, Nickname, MCMadd, MCMmul). Included here so Subpath
# constructors can reference MCMadd/MCMmul for validation.
include(joinpath(@__DIR__, "path-geometry-meta.jl"))

# Resolve JumpBy to a QuinticConnector at build() time.
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

# Default fallback for non-JumpBy segments (ignores K_in_global).
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
# SubpathBuilder (mutable authoring) → Subpath (immutable spec)
# -----------------------------------------------------------------------

"""
    SubpathBuilder()

Mutable Subpath authoring target. Lifecycle:

1. `start!(builder; ...)` seals the start state. Must be called before any
   interior segment.
2. Append interior segments via `straight!`, `bend!`, `helix!`, `catenary!`,
   `jumpby!`.
3. `jumpto!(builder; ...)` seals the end. After this no further segments
   may be added.

`Subpath(builder)` (or `build(builder)`) requires both seals to have been
called.
"""
mutable struct SubpathBuilder
    meta::Vector{AbstractMeta}
    start_point::Union{Nothing, NTuple{3, Float64}}
    start_outgoing_tangent::Union{Nothing, NTuple{3, Float64}}
    start_outgoing_curvature::Union{Nothing, NTuple{3, Float64}}
    segments::Vector{AbstractPathSegment}
    jumpto_point::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_tangent::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_curvature::Union{Nothing, NTuple{3, Float64}}
    jumpto_min_bend_radius::Union{Nothing, Float64}
    jumpto_conserve_path_length::Bool

    SubpathBuilder(; meta::AbstractVector{<:AbstractMeta} = AbstractMeta[]) =
        new(Vector{AbstractMeta}(meta),
            nothing, nothing, nothing,
            AbstractPathSegment[],
            nothing, nothing, nothing, nothing, false)
end

_check_started(b::SubpathBuilder) =
    isnothing(b.start_point) &&
        throw(ArgumentError("SubpathBuilder: call start!() before adding segments"))

_check_unsealed(b::SubpathBuilder) =
    !isnothing(b.jumpto_point) &&
        throw(ArgumentError("SubpathBuilder: subpath already sealed by jumpto!()"))

"""
    start!(builder; point=(0,0,0), outgoing_tangent=(0,0,1),
                    outgoing_curvature=(0,0,0))

Seal the Subpath start state. Throws if `start!` has already been called or
if any interior segment has already been appended.
"""
function start!(b::SubpathBuilder;
                point = (0.0, 0.0, 0.0),
                outgoing_tangent = (0.0, 0.0, 1.0),
                outgoing_curvature = (0.0, 0.0, 0.0))
    !isnothing(b.start_point) &&
        throw(ArgumentError("SubpathBuilder: start!() already called"))
    !isempty(b.segments) &&
        throw(ArgumentError("SubpathBuilder: start!() must be called before any segments"))
    b.start_point              = (Float64(point[1]),              Float64(point[2]),              Float64(point[3]))
    b.start_outgoing_tangent   = (Float64(outgoing_tangent[1]),   Float64(outgoing_tangent[2]),   Float64(outgoing_tangent[3]))
    b.start_outgoing_curvature = (Float64(outgoing_curvature[1]), Float64(outgoing_curvature[2]), Float64(outgoing_curvature[3]))
    return b
end

function straight!(spec::SubpathBuilder; length,
                   meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    _check_started(spec); _check_unsealed(spec)
    push!(spec.segments, StraightSegment(length; meta))
    return spec
end

function bend!(spec::SubpathBuilder; radius::Real, angle::Real, axis_angle::Real = 0.0,
               meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    _check_started(spec); _check_unsealed(spec)
    push!(spec.segments, BendSegment(radius, angle, axis_angle; meta))
    return spec
end

function helix!(spec::SubpathBuilder; radius::Real, pitch::Real, turns::Real,
                axis_angle::Real = 0.0, meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    _check_started(spec); _check_unsealed(spec)
    push!(spec.segments, HelixSegment(radius, pitch, turns, axis_angle; meta))
    return spec
end

function catenary!(spec::SubpathBuilder; a::Real, length::Real, axis_angle::Real = 0.0,
                   meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    _check_started(spec); _check_unsealed(spec)
    push!(spec.segments, CatenarySegment(a, length, axis_angle; meta))
    return spec
end

function jumpby!(spec::SubpathBuilder; delta, tangent = nothing,
                 curvature_out = nothing, min_bend_radius = nothing,
                 meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    _check_started(spec); _check_unsealed(spec)
    push!(spec.segments, JumpBy(delta; tangent_out = tangent,
                                curvature_out = curvature_out,
                                min_bend_radius, meta))
    return spec
end

"""
    jumpto!(builder; point, incoming_tangent=nothing, incoming_curvature=nothing,
                     min_bend_radius=nothing, conserve_path_length=false)

Seal the Subpath end. Stores the terminal connector spec on the builder.
Throws if `start!` was not called, or if the builder is already sealed.
"""
function jumpto!(b::SubpathBuilder;
                 point,
                 incoming_tangent = nothing,
                 incoming_curvature = nothing,
                 min_bend_radius = nothing,
                 conserve_path_length::Bool = false)
    isnothing(b.start_point) &&
        throw(ArgumentError("SubpathBuilder: call start!() before jumpto!()"))
    !isnothing(b.jumpto_point) &&
        throw(ArgumentError("SubpathBuilder: jumpto!() already called"))
    b.jumpto_point = (Float64(point[1]), Float64(point[2]), Float64(point[3]))
    b.jumpto_incoming_tangent   = isnothing(incoming_tangent)   ? nothing :
        (Float64(incoming_tangent[1]),   Float64(incoming_tangent[2]),   Float64(incoming_tangent[3]))
    b.jumpto_incoming_curvature = isnothing(incoming_curvature) ? nothing :
        (Float64(incoming_curvature[1]), Float64(incoming_curvature[2]), Float64(incoming_curvature[3]))
    b.jumpto_min_bend_radius    = isnothing(min_bend_radius) ? nothing : Float64(min_bend_radius)
    b.jumpto_conserve_path_length = conserve_path_length
    return b
end

# -----------------------------------------------------------------------
# Subpath (immutable user-supplied snapshot)
# -----------------------------------------------------------------------

"""
    Subpath

Immutable Subpath specification — user-supplied data only. Geometry queries on
a `Subpath` (arc_length, curvature, position, ...) throw with "call build()
first"; build the Subpath into a `SubpathBuilt` to query.
"""
struct Subpath
    meta::Vector{AbstractMeta}
    start_point::NTuple{3, Float64}
    start_outgoing_tangent::NTuple{3, Float64}
    start_outgoing_curvature::NTuple{3, Float64}
    segments::Vector{AbstractPathSegment}
    jumpto_point::NTuple{3, Float64}
    jumpto_incoming_tangent::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_curvature::Union{Nothing, NTuple{3, Float64}}
    jumpto_min_bend_radius::Union{Nothing, Float64}
    jumpto_conserve_path_length::Bool
end

function Subpath(b::SubpathBuilder)
    isnothing(b.start_point) &&
        throw(ArgumentError("Subpath: builder has no start; call start!() first"))
    isnothing(b.jumpto_point) &&
        throw(ArgumentError("Subpath: builder has no terminal jumpto; call jumpto!() before constructing"))
    return Subpath(deepcopy(b.meta),
                   b.start_point::NTuple{3, Float64},
                   b.start_outgoing_tangent::NTuple{3, Float64},
                   b.start_outgoing_curvature::NTuple{3, Float64},
                   deepcopy(b.segments),
                   b.jumpto_point::NTuple{3, Float64},
                   b.jumpto_incoming_tangent,
                   b.jumpto_incoming_curvature,
                   b.jumpto_min_bend_radius,
                   b.jumpto_conserve_path_length)
end

# Geometry queries on an unbuilt Subpath fail loudly.
arc_length(::Subpath)               = error("Subpath: call build(subpath) before querying arc_length")
arc_length(::Subpath, ::Real, ::Real) = error("Subpath: call build(subpath) before querying arc_length")
curvature(::Subpath, ::Real)        = error("Subpath: call build(subpath) before querying curvature")
geometric_torsion(::Subpath, ::Real) = error("Subpath: call build(subpath) before querying geometric_torsion")
material_twist(::Subpath, ::Real)   = error("Subpath: call build(subpath) before querying material_twist")
position(::Subpath, ::Real)         = error("Subpath: call build(subpath) before querying position")
tangent(::Subpath, ::Real)          = error("Subpath: call build(subpath) before querying tangent")
normal(::Subpath, ::Real)           = error("Subpath: call build(subpath) before querying normal")
binormal(::Subpath, ::Real)         = error("Subpath: call build(subpath) before querying binormal")
frame(::Subpath, ::Real)            = error("Subpath: call build(subpath) before querying frame")

# -----------------------------------------------------------------------
# SubpathBuilt and PathBuilt
# -----------------------------------------------------------------------

"""
    SubpathBuilt

Built form of a `Subpath`. Contains:

- `subpath` — the source spec.
- `placed_segments` — the placed **interior** segments only.
- `jumpto_quintic_connector` — the resolved terminal connector segment.
- `jumpto_placed` — the `PlacedSegment` wrapper for the terminal connector
  (s offset, global origin, global frame at its start). Carried so query
  functions can treat the terminal connector uniformly with interior
  segments.
- `resolved_twists` — twist runs resolved within this Subpath.
- `pending_continuous_first_twist` — true if the first twist anchor was
  deferred for `build(::Vector{SubpathBuilt})` to fill in.

Local arc length runs from `0` to `s_end(::SubpathBuilt)` (computed on demand
from the placed segments + terminal connector).
"""
struct SubpathBuilt
    subpath::Subpath
    placed_segments::Vector{PlacedSegment}        # interior only
    jumpto_quintic_connector::QuinticConnector    # terminal connector
    jumpto_placed::PlacedSegment                  # placement of the terminal connector
    resolved_twists::Vector{ResolvedTwistRate}
    pending_continuous_first_twist::Bool
end

"""
    PathBuilt(subpaths::Vector{SubpathBuilt})

Ordered container of independent `SubpathBuilt`s. Global s offsets and the
total `s_end` are computed on demand to avoid consistency hazards.
"""
struct PathBuilt
    subpaths::Vector{SubpathBuilt}
end

# -----------------------------------------------------------------------
# build()
# -----------------------------------------------------------------------

_safe_normalize(v::AbstractVector) = v ./ sqrt(sum(abs2, v))

# Initial frame derivation: pick N orthogonal to T using a stable convention
# (Gram-Schmidt against a reference axis). Returns (T, N, B).
function _initial_frame_from_tangent(T_tuple::NTuple{3, Float64})
    T = collect(T_tuple)
    Tn = sqrt(T[1]^2 + T[2]^2 + T[3]^2)
    Tn > 0.0 || throw(ArgumentError("start_outgoing_tangent must be non-zero"))
    T = T ./ Tn
    # Pick the world axis least aligned with T to guarantee a non-degenerate cross.
    aT = (abs(T[1]), abs(T[2]), abs(T[3]))
    ref = if aT[1] <= aT[2] && aT[1] <= aT[3]
        [1.0, 0.0, 0.0]
    elseif aT[2] <= aT[3]
        [0.0, 1.0, 0.0]
    else
        [0.0, 0.0, 1.0]
    end
    Nraw = ref - dot(ref, T) .* T
    N = _safe_normalize(Nraw)
    B = cross(T, N)
    return (T, N, B)
end

"""
    build(builder::SubpathBuilder) → SubpathBuilt
    build(sub::Subpath)            → SubpathBuilt

Compile a `Subpath` (or `SubpathBuilder` after `start!` and `jumpto!` have
been called) into a `SubpathBuilt`.
"""
build(b::SubpathBuilder) = build(Subpath(b))

function build(sub::Subpath)
    pos = collect(sub.start_point)
    T_frame, N_frame, B_frame = _initial_frame_from_tangent(sub.start_outgoing_tangent)
    K_in_global = collect(sub.start_outgoing_curvature)

    s_eff  = 0.0
    placed = PlacedSegment[]

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
        N_frame = _safe_normalize(N_end_g - dot(N_end_g, T_frame) * T_frame)
        B_frame = cross(T_frame, N_frame)

        # Update K_in_global for the next segment: κ at end times the unit
        # normal at the end, both expressed globally.
        L_seg = arc_length(seg_placed)
        κ_end = curvature(seg_placed, L_seg)
        K_in_global = κ_end .* N_frame

        s_eff += L_seg
    end

    # Resolve the terminal connector inline. Its destination is global
    # (sub.jumpto_point); transform to local frame to call _build_quintic_connector.
    frame = hcat(N_frame, B_frame, T_frame)
    p1_local  = frame' * (collect(sub.jumpto_point) .- pos)
    chord     = norm(p1_local)
    t_hat_out = isnothing(sub.jumpto_incoming_tangent) ?
        (chord > 1e-15 ? p1_local ./ chord : [0.0, 0.0, 1.0]) :
        normalize(frame' * normalize(collect(sub.jumpto_incoming_tangent)))
    K0_local = frame' * K_in_global
    K1_local = isnothing(sub.jumpto_incoming_curvature) ?
        zeros(eltype(K0_local), 3) :
        frame' * collect(sub.jumpto_incoming_curvature)
    connector = _build_quintic_connector(p1_local, t_hat_out, K0_local, K1_local;
                                         min_bend_radius = sub.jumpto_min_bend_radius)
    # PlacedSegment wrapper for the terminal connector: anchor at the
    # position/frame at the end of the interior segments. Stored alongside
    # the connector so query functions treat it like any other placed segment.
    jumpto_placed = PlacedSegment(connector, s_eff, copy(pos), copy(frame))
    L_conn = arc_length(connector)
    s_end_eff = s_eff + L_conn

    # Resolve twists local to this Subpath. If the first twist anchor has
    # is_continuous=true it is left pending for the PathBuilt build to fix up.
    resolved, pending = _resolve_twists_subpath_local(placed, connector, s_eff,
                                                      Float64(_qc_nominalize(s_end_eff)))
    return SubpathBuilt(sub, placed, connector, jumpto_placed, resolved, pending)
end

# -----------------------------------------------------------------------
# Twist resolution (per Subpath)
# -----------------------------------------------------------------------
# Walk placed segments in order, collect Twist anchors from each segment's
# meta, and emit one ResolvedTwistRate per anchor. Each run extends from the
# anchor's segment start to the next anchor's segment start, or to the end
# of the Subpath. is_continuous=true on a non-first anchor takes its phi_0
# from the prior run's accumulated phase. is_continuous=true on the *first*
# anchor of a Subpath leaves resolution pending — `build(::Vector{SubpathBuilt})`
# fixes it up by inheriting from the prior Subpath.

# Anchor scan also includes the terminal connector (queries on a Subpath's
# domain include the connector segment past the last interior).
function _collect_twist_anchors(placed::Vector{PlacedSegment},
                                connector::QuinticConnector,
                                connector_s_offset::Float64)
    anchors = Tuple{Float64, Twist}[]
    for ps in placed
        twists_here = Twist[]
        for m in segment_meta(ps.segment)
            m isa Twist && push!(twists_here, m)
        end
        if length(twists_here) > 1
            throw(ArgumentError(
                "Subpath build: segment at s_offset_eff = $(ps.s_offset_eff) carries " *
                "$(length(twists_here)) Twist meta entries; at most one is permitted"))
        end
        if !isempty(twists_here)
            push!(anchors, (Float64(_qc_nominalize(ps.s_offset_eff)), twists_here[1]))
        end
    end
    # The terminal connector also carries a meta vector but in the new design
    # it has no `meta` slot exposed for Twist — kept structural. Skip.
    return anchors
end

# Resolve twist runs for a single Subpath. If the first anchor is_continuous=true
# the run is emitted with phi_0=NaN as a sentinel, and `pending_first=true`. The
# PathBuilt-level resolver replaces phi_0 with the inherited value.
function _resolve_twists_subpath_local(placed::Vector{PlacedSegment},
                                       connector::QuinticConnector,
                                       connector_s_offset::Float64,
                                       s_end::Float64)
    anchors = _collect_twist_anchors(placed, connector, connector_s_offset)
    isempty(anchors) && return (ResolvedTwistRate[], false)

    n = length(anchors)
    out = Vector{ResolvedTwistRate}(undef, n)
    prev_phi_0 = 0.0
    prev_run_length = 0.0
    prev_rate::Union{Float64, Function} = 0.0
    pending_first = false

    for i in 1:n
        s_start_i, tw = anchors[i]
        s_run_end = (i < n) ? anchors[i + 1][1] : s_end

        if tw.is_continuous
            if i == 1
                # Defer to PathBuilt-level resolution.
                pending_first = true
                phi_0 = NaN
            else
                phi_0 = prev_phi_0 + _integrate_rate(prev_rate, 0.0, prev_run_length)
            end
        else
            phi_0 = tw.phi_0
        end

        out[i] = ResolvedTwistRate(s_start_i, s_run_end, tw.rate, phi_0)

        prev_phi_0 = phi_0
        prev_run_length = s_run_end - s_start_i
        prev_rate = tw.rate
    end

    return (out, pending_first)
end

# PathBuilt-level fix-up: walk subpaths in order; if a Subpath has
# pending_continuous_first_twist=true, look up the prior Subpath's terminal
# twist phase and rebuild that Subpath's resolved_twists list with the inherited
# phi_0. The first Subpath cannot have pending=true; throw if so.
function _resolve_pending_continuous_twists(builts::Vector{SubpathBuilt})
    n = length(builts)
    out = Vector{SubpathBuilt}(undef, n)
    for i in 1:n
        b = builts[i]
        if !b.pending_continuous_first_twist
            out[i] = b
            continue
        end
        if i == 1
            throw(ArgumentError(
                "PathBuilt: first Subpath has Twist(is_continuous=true) on its first " *
                "anchor, but there is no prior Subpath to inherit phase from"))
        end
        prev = out[i - 1]
        # Compute the prior Subpath's terminal phase: the last resolved twist
        # run's phi_0 plus integral over its length.
        if isempty(prev.resolved_twists)
            throw(ArgumentError(
                "PathBuilt: Subpath $i has pending continuous first twist but " *
                "Subpath $(i-1) has no twist runs to inherit from"))
        end
        last_run = prev.resolved_twists[end]
        run_len  = last_run.s_eff_end - last_run.s_eff_start
        phi_end  = last_run.phi_0 + _integrate_rate(last_run.rate, 0.0, run_len)

        # Rebuild this Subpath's resolved_twists list, replacing the first run's phi_0.
        rt_old = b.resolved_twists
        rt_new = Vector{ResolvedTwistRate}(undef, length(rt_old))
        # First run: inherit phi_end. Then propagate forward through subsequent
        # runs that were resolved with prev_phi_0 originating from the (NaN)
        # first run — re-run the phi_0 chain from this corrected start.
        prev_phi_0 = phi_end
        prev_run_length = 0.0
        prev_rate::Union{Float64, Function} = 0.0
        for k in 1:length(rt_old)
            r = rt_old[k]
            if k == 1
                phi_0 = phi_end
            else
                # Re-derive: was tw.is_continuous? We don't have the raw Twist
                # object handy, so detect by NaN propagation in the original
                # resolution. The simple rule: if the original phi_0 isn't NaN
                # and was deterministic (from tw.phi_0), keep it. If it was
                # derived from prev_phi_0 chain (k > 1 and is_continuous), it
                # would have been computed without the corrected base.
                # To be safe, recompute: phi_0 = prev_phi_0 + ∫ prev_rate over prev_run_length
                # only if the original is "close" to that expression. We
                # cannot disambiguate without the source Twist meta, so we
                # walk the segments again below.
                phi_0 = r.phi_0
            end
            rt_new[k] = ResolvedTwistRate(r.s_eff_start, r.s_eff_end, r.rate, phi_0)
            prev_phi_0 = phi_0
            prev_run_length = r.s_eff_end - r.s_eff_start
            prev_rate = r.rate
        end

        # If the first run's phi_0 was NaN and there are downstream is_continuous
        # runs that inherited from it, those would also have been computed with
        # prev_phi_0=NaN and need recomputation. Re-resolve from raw anchors with
        # phi_end as the seed.
        if any(r -> isnan(r.phi_0), rt_old) && length(rt_old) >= 2
            anchors = _collect_twist_anchors(b.placed_segments, b.jumpto_quintic_connector, 0.0)
            # Recompute with first phi_0 = phi_end (forced is_continuous semantics)
            n_a = length(anchors)
            rt_new = Vector{ResolvedTwistRate}(undef, n_a)
            prev_phi_0 = phi_end
            prev_run_length = 0.0
            prev_rate = 0.0
            # s_end for the runs is the same as the original
            s_end = isempty(rt_old) ? 0.0 : rt_old[end].s_eff_end
            for k in 1:n_a
                s_start_k, tw = anchors[k]
                s_run_end = (k < n_a) ? anchors[k + 1][1] : s_end
                if tw.is_continuous
                    if k == 1
                        phi_0 = phi_end
                    else
                        phi_0 = prev_phi_0 + _integrate_rate(prev_rate, 0.0, prev_run_length)
                    end
                else
                    phi_0 = tw.phi_0
                end
                rt_new[k] = ResolvedTwistRate(s_start_k, s_run_end, tw.rate, phi_0)
                prev_phi_0 = phi_0
                prev_run_length = s_run_end - s_start_k
                prev_rate = tw.rate
            end
        end

        out[i] = SubpathBuilt(b.subpath, b.placed_segments, b.jumpto_quintic_connector,
                              b.jumpto_placed, rt_new, false)
    end
    return out
end

# -----------------------------------------------------------------------
# build(::Vector{Subpath} | ::Vector{SubpathBuilt} | ::SubpathBuilt) → PathBuilt
# -----------------------------------------------------------------------

# Approximate equality on NTuple{3, Float64} for endpoint conformity checks.
_tuple_isapprox(a::NTuple{3, Float64}, b::NTuple{3, Float64};
                atol::Float64 = 1e-9, rtol::Float64 = 1e-9) =
    all(isapprox(a[i], b[i]; atol = atol, rtol = rtol) for i in 1:3)

# Conformity check between adjacent Subpaths in a PathBuilt. Subpaths are
# independent in spec but must be stackable in order: N-1's jumpto endpoint
# state must match N's start state.
function _check_subpath_conformity(prev::Subpath, cur::Subpath, idx::Int)
    if !_tuple_isapprox(prev.jumpto_point, cur.start_point)
        throw(ArgumentError(
            "PathBuilt: Subpath $(idx-1) jumpto_point $(prev.jumpto_point) " *
            "does not match Subpath $idx start_point $(cur.start_point)"))
    end
    # jumpto_incoming_tangent default is "chord direction"; nothing to compare
    # if prev has nothing. If both supplied, compare; if one supplied, that's a
    # mismatch.
    pt = prev.jumpto_incoming_tangent
    ct = cur.start_outgoing_tangent
    if !isnothing(pt)
        if !_tuple_isapprox(pt, ct)
            throw(ArgumentError(
                "PathBuilt: Subpath $(idx-1) jumpto_incoming_tangent $pt " *
                "does not match Subpath $idx start_outgoing_tangent $ct"))
        end
    end
    pk = prev.jumpto_incoming_curvature
    ck = cur.start_outgoing_curvature
    if !isnothing(pk)
        if !_tuple_isapprox(pk, ck)
            throw(ArgumentError(
                "PathBuilt: Subpath $(idx-1) jumpto_incoming_curvature $pk " *
                "does not match Subpath $idx start_outgoing_curvature $ck"))
        end
    end
    return nothing
end

"""
    build(builts::Vector{SubpathBuilt}) → PathBuilt
    build(subpaths::Vector{Subpath})   → PathBuilt
    build(spb::SubpathBuilt)           → PathBuilt

Stitch already-built `SubpathBuilt`s into a `PathBuilt`, validating that
adjacent Subpaths' endpoint states agree, and resolving any cross-Subpath
twist continuity. The vector-of-Subpath form builds each first; the single
`SubpathBuilt` form wraps a length-1 PathBuilt.
"""
function build(builts::Vector{SubpathBuilt})
    isempty(builts) && throw(ArgumentError("PathBuilt: at least one SubpathBuilt required"))
    for i in 2:length(builts)
        _check_subpath_conformity(builts[i-1].subpath, builts[i].subpath, i)
    end
    fixed = _resolve_pending_continuous_twists(builts)
    return PathBuilt(fixed)
end

build(subpaths::Vector{Subpath}) = build([build(sp) for sp in subpaths])

build(spb::SubpathBuilt) = build(SubpathBuilt[spb])


# -----------------------------------------------------------------------
# Segment lookup helpers
# -----------------------------------------------------------------------

# Walk the SubpathBuilt's interior placed segments and the terminal connector
# in order. Returns (PlacedSegment, s_local) where s_local is clamped into
# the segment's domain.
function _find_placed_segment(b::SubpathBuilt, s)
    s_eff = 0.0
    for ps in b.placed_segments
        seg_len = arc_length(ps.segment)
        seg_len_nom = Float64(_qc_nominalize(seg_len))
        s_eff_next  = s_eff + seg_len_nom
        if s <= s_eff_next + 1e-12
            s_local = clamp(s - s_eff, zero(seg_len), seg_len)
            return ps, s_local
        end
        s_eff = s_eff_next
    end
    # Past all interior segments → terminal connector.
    ps_t      = b.jumpto_placed
    seg_len_t = arc_length(ps_t.segment)
    s_local   = clamp(s - s_eff, zero(seg_len_t), seg_len_t)
    return ps_t, s_local
end

function _local_to_global(ps::PlacedSegment, v_local::AbstractVector)
    return ps.frame * v_local
end

# Compute s_end of a SubpathBuilt on demand: sum of interior segment lengths
# plus the terminal connector arc length.
function s_end(b::SubpathBuilt)
    total = zero(arc_length(b.jumpto_quintic_connector))
    for ps in b.placed_segments
        total = total + arc_length(ps.segment)
    end
    return total + arc_length(b.jumpto_quintic_connector)
end

# -----------------------------------------------------------------------
# Differential geometry interface on SubpathBuilt
# -----------------------------------------------------------------------

arc_length(b::SubpathBuilt) = s_end(b)

function arc_length(::SubpathBuilt, s1, s2)
    @assert s2 >= s1 "arc_length: require s2 >= s1"
    return s2 - s1
end

function curvature(b::SubpathBuilt, s::Real)
    ps, s_local = _find_placed_segment(b, s)
    return curvature(ps.segment, s_local)
end

function geometric_torsion(b::SubpathBuilt, s::Real)
    ps, s_local = _find_placed_segment(b, s)
    return geometric_torsion(ps.segment, s_local)
end

"""
    material_twist(b, s)

Material twist rate (rad/m) at local arc length `s`, summed over all resolved
twist runs that contain `s`. Runs are disjoint by construction; the sum
exists only as a robustness guard.
"""
function material_twist(b::SubpathBuilt, s)
    τ = zero(s isa AbstractFloat ? s : Float64(s))
    for r in b.resolved_twists
        if r.s_eff_start <= s <= r.s_eff_end
            τ += r.rate isa Function ? r.rate(s - r.s_eff_start) : r.rate
        end
    end
    return τ
end

function position(b::SubpathBuilt, s::Real)
    ps, s_local = _find_placed_segment(b, s)
    return ps.origin + _local_to_global(ps, position_local(ps.segment, s_local))
end

function tangent(b::SubpathBuilt, s::Real)
    ps, s_local = _find_placed_segment(b, s)
    return _local_to_global(ps, tangent_local(ps.segment, s_local))
end

function normal(b::SubpathBuilt, s::Real)
    ps, s_local = _find_placed_segment(b, s)
    return _local_to_global(ps, normal_local(ps.segment, s_local))
end

function binormal(b::SubpathBuilt, s::Real)
    ps, s_local = _find_placed_segment(b, s)
    return _local_to_global(ps, binormal_local(ps.segment, s_local))
end

function frame(b::SubpathBuilt, s::Real)
    T = tangent(b, s)
    N = normal(b, s)
    Bi = binormal(b, s)
    κ = curvature(b, s)
    τ = geometric_torsion(b, s)
    m = material_twist(b, s)
    return (; position = position(b, s), tangent = T, normal = N, binormal = Bi,
              curvature = κ, geometric_torsion = τ, material_twist = m)
end

# -----------------------------------------------------------------------
# Endpoint access on SubpathBuilt
# -----------------------------------------------------------------------

start_point(b::SubpathBuilt)   = position(b, 0.0)
end_point(b::SubpathBuilt)     = position(b, Float64(_qc_nominalize(s_end(b))))
start_tangent(b::SubpathBuilt) = tangent(b, 0.0)
end_tangent(b::SubpathBuilt)   = tangent(b, Float64(_qc_nominalize(s_end(b))))

# -----------------------------------------------------------------------
# Path measures
# -----------------------------------------------------------------------

path_length(b::SubpathBuilt) = arc_length(b)

function cartesian_distance(b::SubpathBuilt, s1::Real, s2::Real)
    return norm(position(b, s2) - position(b, s1))
end

function bounding_box(b::SubpathBuilt; n::Int = 512)
    s0 = 0.0
    s1 = Float64(_qc_nominalize(s_end(b)))
    ss = range(s0, s1; length = n)
    pts = [position(b, s) for s in ss]
    lo = minimum(reduce(hcat, pts); dims = 2) |> vec
    hi = maximum(reduce(hcat, pts); dims = 2) |> vec
    return (; lo, hi)
end

function _all_placed_segs(b::SubpathBuilt)
    # Returns interior + terminal. Used by total_* iterators.
    result = Vector{PlacedSegment}(undef, length(b.placed_segments) + 1)
    @inbounds for i in eachindex(b.placed_segments)
        result[i] = b.placed_segments[i]
    end
    result[end] = b.jumpto_placed
    return result
end

function total_turning_angle(b::SubpathBuilt)
    total = 0.0
    for ps in _all_placed_segs(b)
        seg = ps.segment
        if seg isa StraightSegment
            # κ = 0
        elseif seg isa BendSegment
            total += abs(seg.angle)
        elseif seg isa HelixSegment
            total += curvature(seg, 0.0) * arc_length(seg)
        else
            n = 64
            ss = range(0.0, arc_length(seg); length = n + 1)
            h = ss[2] - ss[1]
            total += h * sum(curvature(seg, s) for s in ss)
        end
    end
    return total
end

function total_torsion(b::SubpathBuilt)
    total = 0.0
    for ps in _all_placed_segs(b)
        seg = ps.segment
        if seg isa HelixSegment
            total += geometric_torsion(seg, 0.0) * arc_length(seg)
        elseif seg isa StraightSegment || seg isa BendSegment || seg isa CatenarySegment
            # τ_geom = 0
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
    total_material_twist(b; s_start, s_end, rtol = 1e-8, atol = 0.0) → Float64

Integrated material twist ``∫ τ_{\\mathrm{mat}}(s) \\, ds`` over local arc length
from `s_start` to `s_end` (defaults: full Subpath). Both endpoints must lie
in `[0, s_end(b)]`.
"""
function total_material_twist(
    b::SubpathBuilt;
    s_start::Real = 0.0,
    s_end::Real   = s_end(b),
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
    ps1 = Float64(_qc_nominalize(arc_length(b)))
    if !(ps0 - 1e-12 <= s_lo <= ps1 + 1e-12) || !(ps0 - 1e-12 <= s_hi <= ps1 + 1e-12)
        throw(ArgumentError(
            "total_material_twist: require 0 ≤ s ≤ s_end(b) for both endpoints; " *
            "got [$(s_lo), $(s_hi)] m vs subpath domain [$(ps0), $(ps1)] m"))
    end
    s_lo == s_hi && return 0.0

    total = 0.0
    rtolf = Float64(rtol)
    atolf = Float64(atol)
    for r in b.resolved_twists
        a = max(r.s_eff_start, s_lo)
        bb = min(r.s_eff_end, s_hi)
        bb <= a && continue
        total += _integrate_rate(r.rate, a - r.s_eff_start, bb - r.s_eff_start;
                                 rtol = rtolf, atol = atolf)
    end
    return total
end

"""
    total_frame_rotation(b; s_start, s_end, rtol = 1e-8, atol = 0.0) → Float64

Total frame rotation `∫ (τ_geom + Ω_material) ds` over local arc length from
`s_start` to `s_end`.
"""
function total_frame_rotation(
    b::SubpathBuilt;
    s_start::Real = 0.0,
    s_end::Real   = s_end(b),
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
    ps1 = Float64(_qc_nominalize(arc_length(b)))
    if !(ps0 - 1e-12 <= s_lo <= ps1 + 1e-12) || !(ps0 - 1e-12 <= s_hi <= ps1 + 1e-12)
        throw(ArgumentError(
            "total_frame_rotation: require 0 ≤ s ≤ s_end(b) for both endpoints; " *
            "got [$(s_lo), $(s_hi)] m vs subpath domain [$(ps0), $(ps1)] m"))
    end
    s_lo == s_hi && return 0.0

    τ_total = 0.0
    rtolf = Float64(rtol)
    atolf = Float64(atol)
    for ps in _all_placed_segs(b)
        seg = ps.segment
        seg_s_start_nom = Float64(_qc_nominalize(ps.s_offset_eff))
        seg_s_end_nom   = seg_s_start_nom +
                          Float64(_qc_nominalize(arc_length(seg)))
        a_nom = max(seg_s_start_nom, s_lo)
        b_nom = min(seg_s_end_nom,   s_hi)
        b_nom <= a_nom && continue

        L_seg  = arc_length(seg)
        seg_span_nom = seg_s_end_nom - seg_s_start_nom
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

    Ω_total = total_material_twist(b; s_start = s_lo, s_end = s_hi,
                                   rtol = rtolf, atol = atolf)
    return τ_total + Ω_total
end

"""
    writhe(b; n) → Float64

Writhe of the Subpath: numerical double integral on `n` samples.
"""
function writhe(b::SubpathBuilt; n::Int = 256)
    s0 = 0.0
    s1 = Float64(_qc_nominalize(s_end(b)))
    ss = collect(range(s0, s1; length = n))
    rs = [position(b, s) for s in ss]
    ts = [tangent(b, s)  for s in ss]
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

One evaluated point on a Subpath: arc-length coordinate `s` plus all Frenet
frame quantities.
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

Dense samples of a Subpath/PathBuilt over `[s_start, s_end]`.
"""
struct PathSample
    samples :: Vector{Sample}
    s_start :: Float64
    s_end   :: Float64
    n       :: Int
end

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

function _segment_point_budget(
    ps::PlacedSegment,
    b::SubpathBuilt,
    s_lo::Float64,
    s_hi::Float64,
    fidelity::Float64,
)
    seg = ps.segment
    seg_s_start = Float64(_qc_nominalize(ps.s_offset_eff))
    seg_s_end   = seg_s_start + Float64(_qc_nominalize(arc_length(seg)))

    a = max(s_lo, seg_s_start)
    bb = min(s_hi, seg_s_end)
    bb <= a && return 2

    seg_len = seg_s_end - seg_s_start
    frac = seg_len > 0.0 ? (bb - a) / seg_len : 1.0

    geom_angle  = _budget_scalar(_segment_total_angle(seg) * frac)
    geom_budget = max(2, ceil(Int, fidelity * geom_angle / (2π) * 32))

    twist_total  = total_material_twist(b; s_start = a, s_end = bb, rtol = 1e-3)
    twist_angle  = abs(_budget_scalar(twist_total))
    twist_budget = max(2, ceil(Int, fidelity * twist_angle / (2π) * 32))

    return max(geom_budget, twist_budget)
end

"""
    sample_path(b, s1, s2; fidelity = 1.0) → PathSample

Adaptive sampling of a `SubpathBuilt` over `[s1, s2]`.
"""
function sample_path(b::SubpathBuilt, s1::Real, s2::Real; fidelity::Float64 = 1.0)
    @assert s2 > s1    "sample_path: require s2 > s1"
    @assert fidelity > 0.0 "sample_path: fidelity must be positive"

    s_lo = Float64(_qc_nominalize(s1))
    s_hi = Float64(_qc_nominalize(s2))

    all_s = Float64[]
    for ps in _all_placed_segs(b)
        seg_s_start = Float64(_qc_nominalize(ps.s_offset_eff))
        seg_s_end   = seg_s_start + Float64(_qc_nominalize(arc_length(ps.segment)))

        a = max(s_lo, seg_s_start)
        bb = min(s_hi, seg_s_end)
        bb <= a && continue

        n_seg = _segment_point_budget(ps, b, s_lo, s_hi, fidelity)
        seg_ss = collect(range(a, bb; length = n_seg))

        if isempty(all_s)
            append!(all_s, seg_ss)
        else
            start_idx = (seg_ss[1] ≈ all_s[end]) ? 2 : 1
            append!(all_s, @view seg_ss[start_idx:end])
        end
    end

    if isempty(all_s)
        all_s = [s_lo, s_hi]
    elseif length(all_s) == 1
        push!(all_s, s_hi)
    end

    n = length(all_s)
    samples = Vector{Sample}(undef, n)
    for i in eachindex(all_s)
        fr = frame(b, all_s[i])
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

# -----------------------------------------------------------------------
# Breakpoints
# -----------------------------------------------------------------------

function normalize_breakpoints(breakpoints::AbstractVector{<:Real})
    return sort(unique(copy(breakpoints)))
end

function path_segment_breakpoints(b::SubpathBuilt)
    points = Float64[0.0]
    for ps in b.placed_segments
        push!(points, Float64(_qc_nominalize(ps.s_offset_eff)))
        push!(points, Float64(_qc_nominalize(
            ps.s_offset_eff + arc_length(ps.segment))))
    end
    push!(points, Float64(_qc_nominalize(b.jumpto_placed.s_offset_eff)))
    push!(points, Float64(_qc_nominalize(arc_length(b))))
    return normalize_breakpoints(points)
end

function path_twist_breakpoints(b::SubpathBuilt)
    points = Float64[0.0, Float64(_qc_nominalize(arc_length(b)))]
    for r in b.resolved_twists
        push!(points, r.s_eff_start)
        push!(points, r.s_eff_end)
    end
    return normalize_breakpoints(points)
end

function breakpoints(b::SubpathBuilt)
    return normalize_breakpoints(vcat(path_segment_breakpoints(b),
                                      path_twist_breakpoints(b)))
end

function sample(b::SubpathBuilt, s_values)
    return [frame(b, s) for s in s_values]
end

function sample_uniform(b::SubpathBuilt; n::Int = 256)
    ss = range(0.0, Float64(_qc_nominalize(arc_length(b))); length = n)
    return sample(b, ss)
end

# -----------------------------------------------------------------------
# PathBuilt query interface — single _find_subpath glue + @eval forwards
# -----------------------------------------------------------------------

# Cumulative offsets of subpaths within a PathBuilt, computed on demand.
function s_offsets(p::PathBuilt)
    n = length(p.subpaths)
    offs = Vector{Float64}(undef, n)
    cum = 0.0
    @inbounds for i in 1:n
        offs[i] = cum
        cum += Float64(_qc_nominalize(arc_length(p.subpaths[i])))
    end
    return offs
end

s_end(p::PathBuilt) = sum(Float64(_qc_nominalize(arc_length(b))) for b in p.subpaths;
                          init = 0.0)

arc_length(p::PathBuilt) = s_end(p)

function arc_length(::PathBuilt, s1, s2)
    @assert s2 >= s1 "arc_length: require s2 >= s1"
    return s2 - s1
end

path_length(p::PathBuilt) = arc_length(p)

function _find_subpath(p::PathBuilt, s)
    n = length(p.subpaths)
    n == 0 && error("PathBuilt is empty")
    offs = s_offsets(p)
    @inbounds for i in 1:n
        local_end = offs[i] + Float64(_qc_nominalize(arc_length(p.subpaths[i])))
        if s <= local_end + 1e-12 || i == n
            return p.subpaths[i], s - offs[i]
        end
    end
    error("PathBuilt: s = $s out of path bounds")
end

# Forward all "point-query at s" methods through _find_subpath.
for f in (:curvature, :geometric_torsion, :material_twist,
          :position, :tangent, :normal, :binormal, :frame)
    @eval function $f(p::PathBuilt, s::Real)
        sb, s_local = _find_subpath(p, s)
        return $f(sb, s_local)
    end
end

# Endpoint access
start_point(p::PathBuilt)   = position(p, 0.0)
end_point(p::PathBuilt)     = position(p, s_end(p))
start_tangent(p::PathBuilt) = tangent(p, 0.0)
end_tangent(p::PathBuilt)   = tangent(p, s_end(p))

function cartesian_distance(p::PathBuilt, s1::Real, s2::Real)
    return norm(position(p, s2) - position(p, s1))
end

function bounding_box(p::PathBuilt; n::Int = 512)
    ss = range(0.0, s_end(p); length = n)
    pts = [position(p, s) for s in ss]
    lo = minimum(reduce(hcat, pts); dims = 2) |> vec
    hi = maximum(reduce(hcat, pts); dims = 2) |> vec
    return (; lo, hi)
end

# Aggregate measures: sum across subpaths.
total_turning_angle(p::PathBuilt) = sum(total_turning_angle(b) for b in p.subpaths;
                                        init = 0.0)
total_torsion(p::PathBuilt)       = sum(total_torsion(b)       for b in p.subpaths;
                                        init = 0.0)

function total_material_twist(p::PathBuilt;
                              s_start::Real = 0.0,
                              s_end::Real   = s_end(p),
                              rtol::Real    = 1e-8,
                              atol::Real    = 0.0)
    s_lo = Float64(_qc_nominalize(s_start))
    s_hi = Float64(_qc_nominalize(s_end))
    s_lo == s_hi && return 0.0
    total = 0.0
    offs = s_offsets(p)
    for i in eachindex(p.subpaths)
        L = Float64(_qc_nominalize(arc_length(p.subpaths[i])))
        a_local = max(0.0, s_lo - offs[i])
        b_local = min(L,   s_hi - offs[i])
        b_local <= a_local && continue
        total += total_material_twist(p.subpaths[i];
                                      s_start = a_local, s_end = b_local,
                                      rtol = rtol, atol = atol)
    end
    return total
end

function total_frame_rotation(p::PathBuilt;
                              s_start::Real = 0.0,
                              s_end::Real   = s_end(p),
                              rtol::Real    = 1e-8,
                              atol::Real    = 0.0)
    s_lo = Float64(_qc_nominalize(s_start))
    s_hi = Float64(_qc_nominalize(s_end))
    s_lo == s_hi && return 0.0
    total = 0.0
    offs = s_offsets(p)
    for i in eachindex(p.subpaths)
        L = Float64(_qc_nominalize(arc_length(p.subpaths[i])))
        a_local = max(0.0, s_lo - offs[i])
        b_local = min(L,   s_hi - offs[i])
        b_local <= a_local && continue
        total += total_frame_rotation(p.subpaths[i];
                                      s_start = a_local, s_end = b_local,
                                      rtol = rtol, atol = atol)
    end
    return total
end

function breakpoints(p::PathBuilt)
    isempty(p.subpaths) && return Float64[]
    offs = s_offsets(p)
    bps = Float64[]
    for i in eachindex(p.subpaths)
        for x in breakpoints(p.subpaths[i])
            push!(bps, x + offs[i])
        end
    end
    return normalize_breakpoints(bps)
end

function sample(p::PathBuilt, s_values)
    return [frame(p, s) for s in s_values]
end

function sample_uniform(p::PathBuilt; n::Int = 256)
    ss = range(0.0, s_end(p); length = n)
    return sample(p, ss)
end
