"""
path-geometry.jl

Three-dimensional path geometry for smooth space curves.


# Three layers

- `PathSpecBuilder` (mutable) — authoring target for the bang-DSL below.
- `PathSpec` (immutable) — frozen snapshot: `segments` + `s_start`.
- `PathSpecCached` (immutable) — derived layout: `spec` + `placed_segments` + `s_end`.

`build(builder_or_spec) → PathSpecCached` runs the placement loop. Queries
(`position`, `tangent`, etc.) take a `PathSpecCached`.

# Assembly

Two approaches are supported and can be freely mixed on a `PathSpecBuilder`:

(1) Sliding-frame approach: each segment is specified relative to the frame left by
    the previous segment. The tangent direction at the start of each new segment is
    exactly the tangent at the end of the previous one, so continuity is structural.

        spec = PathSpecBuilder()
        straight!(spec; length=1.0)
        bend!(spec; radius=0.05, angle=π/2)

(2) Endpoint approach: specify the displacement (jumpby!) or absolute destination
    (jumpto!) and an optional outgoing tangent. The connecting segment is an Euler
    spiral (clothoid). The incoming tangent is always the current frame tangent.

        jumpby!(spec; delta=(0.0, 0.0, 0.5))
        jumpto!(spec; destination=(1.0, 0.0, 0.5), tangent=(0.0, 1.0, 0.0))

# Material twist

Material twist is attached as per-segment meta via `Twist <: AbstractMeta`.
A `Twist` placed in a segment's `meta` vector starts a twist run at that
segment's start and continues until the next `Twist`-bearing segment, or
until the path ends.

    spec = PathSpecBuilder()
    straight!(spec; length = 1.0,
              meta = [Twist(; rate = 2π, phi_0 = 0.0)])           # constant rate
    straight!(spec; length = 1.0,
              meta = [Twist(; rate = s -> sin(s), is_continuous = true)])
    build(spec)

`rate` may be a `Real` (constant rad/m) or a `Function` `rate(s_local)` of
run-local arc length (`s_local = 0` at the start of the run). With
`is_continuous = true` the resolver computes `phi_0` from the prior run's
accumulated phase.

# Shrinkage

`path-geometry.jl` is intentionally free of shrinkage. A `Path` represents pure,
temperature-independent geometry. Thermal contraction and related length-scaling
are applied by `fiber-path-shrinkage.jl` as a `shrink(path, α) → Path` transform.

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
# HermiteConnector (and its JumpBy/JumpTo precursors) is Float64-only; see its
# docstring.

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

# -----------------------------------------------------------------------
# JumpBy and JumpTo  (stubs)
# -----------------------------------------------------------------------

"""
    JumpBy(delta, tangent_out, min_bend_radius)

Connects the current position to current_position + delta using a smooth function.
The incoming tangent is the current sliding frame tangent.  `tangent_out` is the
desired outgoing tangent direction; if nothing, the connector computes an implicit
tangent that minimizes curvature variation (requires a nonlinear solve).

`min_bend_radius` (metres, nothing = unconstrained) sets a lower bound on the
radius of curvature of the Hermite connector.  The tangent handle length is
extended beyond the chord default when necessary to keep κ ≤ 1/R_min.
Currently a stub.

"""
struct JumpBy <: AbstractPathSegment
    delta::NTuple{3, Float64}
    tangent_out::Union{Nothing, NTuple{3, Float64}}
    min_bend_radius::Union{Nothing, Float64}
    meta::Vector{AbstractMeta}
end

function JumpBy(delta; tangent_out = nothing, min_bend_radius = nothing,
                meta = AbstractMeta[])
    d = (Float64(delta[1]), Float64(delta[2]), Float64(delta[3]))
    t = isnothing(tangent_out) ? nothing :
        (Float64(tangent_out[1]), Float64(tangent_out[2]), Float64(tangent_out[3]))
    r = isnothing(min_bend_radius) ? nothing : Float64(min_bend_radius)
    JumpBy(d, t, r, Vector{AbstractMeta}(meta))
end

"""
    JumpTo(destination, tangent_out, min_bend_radius)

Connects the current position to the fixed lab-frame `destination` using a smooth
function.

`min_bend_radius` (metres, nothing = unconstrained) sets a lower bound on the
radius of curvature of the Hermite connector.  The tangent handle length is
extended beyond the chord default when necessary to keep κ ≤ 1/R_min.
Currently a stub.
"""
struct JumpTo <: AbstractPathSegment
    destination::NTuple{3, Float64}
    tangent_out::Union{Nothing, NTuple{3, Float64}}
    min_bend_radius::Union{Nothing, Float64}
    meta::Vector{AbstractMeta}
end

function JumpTo(destination; tangent_out = nothing, min_bend_radius = nothing,
                meta = AbstractMeta[])
    d = (Float64(destination[1]), Float64(destination[2]), Float64(destination[3]))
    t = isnothing(tangent_out) ? nothing :
        (Float64(tangent_out[1]), Float64(tangent_out[2]), Float64(tangent_out[3]))
    r = isnothing(min_bend_radius) ? nothing : Float64(min_bend_radius)
    JumpTo(d, t, r, Vector{AbstractMeta}(meta))
end

# JumpBy and JumpTo are context-dependent: geometry is resolved at build() time
# into a HermiteConnector.  Calling these methods on the raw structs is unsupported.
for T in (JumpBy, JumpTo)
    @eval begin
        arc_length(::$T)             = error($(string(T)) * ": call build() to resolve jump geometry")
        curvature(::$T, ::Real)      = error($(string(T)) * ": call build() to resolve jump geometry")
        geometric_torsion(::$T, ::Real) = 0.0
        position_local(::$T, ::Real) = error($(string(T)) * ": call build() to resolve jump geometry")
        tangent_local(::$T, ::Real)  = error($(string(T)) * ": call build() to resolve jump geometry")
        normal_local(::$T, ::Real)   = error($(string(T)) * ": call build() to resolve jump geometry")
        binormal_local(::$T, ::Real) = error($(string(T)) * ": call build() to resolve jump geometry")
        end_position_local(::$T)     = error($(string(T)) * ": call build() to resolve jump geometry")
        end_frame_local(::$T)        = error($(string(T)) * ": call build() to resolve jump geometry")
    end
end

# -----------------------------------------------------------------------
# HermiteConnector  (resolved form of JumpBy / JumpTo)
# -----------------------------------------------------------------------

# Returns a unit vector perpendicular to t (arbitrary orientation).
function _perp_unit(t::Vector{Float64})
    ref = abs(t[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
    n   = ref .- dot(ref, t) .* t
    nn  = norm(n)
    return nn > 1e-12 ? n ./ nn : [1.0, 0.0, 0.0]
end

"""
    HermiteConnector

Internal segment produced by resolving a `JumpBy` or `JumpTo` at `build()` time.
Implements a cubic Hermite spline in local frame coordinates, connecting
(0,0,0) with incoming tangent ẑ to `p1_local` with outgoing tangent `t_hat_out`.
Tangent vectors are scaled by the chord length.

A 256-point Gauss-quadrature arc-length table enables arc-length reparameterisation
of all interface methods.  Curvature and geometric torsion are computed analytically
from the polynomial derivatives.

!!! note "MCM compatibility"
    `HermiteConnector`, `JumpBy`, and `JumpTo` are Float64-only. Their internals
    (Newton inversion of the arc-length table, `_perp_unit` fallback, tuple
    arithmetic) are not lifted to `Particles`. A `PathSpecBuilder` that mixes
    an uncertain `BendSegment`/`HelixSegment`/etc. with a `jumpby!`/`jumpto!`
    will fail at `build()`.
"""
struct HermiteConnector <: AbstractPathSegment
    a0        :: NTuple{3,Float64}   # P(t) = a0 + a1 t + a2 t² + a3 t³, t ∈ [0,1]
    a1        :: NTuple{3,Float64}
    a2        :: NTuple{3,Float64}
    a3        :: NTuple{3,Float64}
    s_table   :: Vector{Float64}     # s_table[i] = arc-length at t = (i-1)/(n-1)
    meta      :: Vector{AbstractMeta}
end

HermiteConnector(a0, a1, a2, a3, s_table; meta = AbstractMeta[]) =
    HermiteConnector(a0, a1, a2, a3, s_table, Vector{AbstractMeta}(meta))

const _HC_GAUSS4_NODES   = (-0.8611363115940526, -0.3399810435848563,
                              0.3399810435848563,  0.8611363115940526)
const _HC_GAUSS4_WEIGHTS = ( 0.3478548451374538,  0.6521451548625461,
                              0.6521451548625461,  0.3478548451374538)

function _hc_dp(a1::NTuple{3,Float64}, a2::NTuple{3,Float64}, a3::NTuple{3,Float64}, t::Float64)
    t2 = t * t
    return (a1[1] + 2*a2[1]*t + 3*a3[1]*t2,
            a1[2] + 2*a2[2]*t + 3*a3[2]*t2,
            a1[3] + 2*a2[3]*t + 3*a3[3]*t2)
end

function _hc_speed(a1, a2, a3, t::Float64)
    dp = _hc_dp(a1, a2, a3, t)
    return sqrt(dp[1]^2 + dp[2]^2 + dp[3]^2)
end

function _hc_ddp(a2::NTuple{3,Float64}, a3::NTuple{3,Float64}, t::Float64)
    return (2*a2[1] + 6*a3[1]*t, 2*a2[2] + 6*a3[2]*t, 2*a2[3] + 6*a3[3]*t)
end

function _hc_build_table(a1, a2, a3, n::Int)
    s = zeros(n)
    for i in 2:n
        t0 = (i - 2) / (n - 1)
        t1 = (i - 1) / (n - 1)
        tm, th = (t0 + t1) / 2, (t1 - t0) / 2
        s[i] = s[i-1] + th * (
            _HC_GAUSS4_WEIGHTS[1] * _hc_speed(a1, a2, a3, tm + th * _HC_GAUSS4_NODES[1]) +
            _HC_GAUSS4_WEIGHTS[2] * _hc_speed(a1, a2, a3, tm + th * _HC_GAUSS4_NODES[2]) +
            _HC_GAUSS4_WEIGHTS[3] * _hc_speed(a1, a2, a3, tm + th * _HC_GAUSS4_NODES[3]) +
            _HC_GAUSS4_WEIGHTS[4] * _hc_speed(a1, a2, a3, tm + th * _HC_GAUSS4_NODES[4]))
    end
    return s
end

function _hc_t_from_s(seg::HermiteConnector, s_target::Float64)
    s_table = seg.s_table
    n = length(s_table)
    L = s_table[end]
    L < 1e-15 && return 0.0
    sc = clamp(s_target, 0.0, L)

    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >> 1
        s_table[mid] <= sc ? (lo = mid) : (hi = mid)
    end
    t0 = (lo - 1) / (n - 1)
    t1 = (hi - 1) / (n - 1)
    s0 = s_table[lo]
    ds = s_table[hi] - s0

    t = ds > 1e-15 ? t0 + (sc - s0) / ds * (t1 - t0) : t0

    # Newton refinement: s(t) ≈ s0 + (t-t0)*speed_at_midpoint
    a1, a2, a3 = seg.a1, seg.a2, seg.a3
    for _ in 1:2
        spd_mid = _hc_speed(a1, a2, a3, (t0 + t) / 2)
        spd_mid < 1e-15 && break
        s_est = s0 + (t - t0) * spd_mid
        spd   = _hc_speed(a1, a2, a3, t)
        spd   < 1e-15 && break
        t = clamp(t - (s_est - sc) / spd, t0, t1)
    end
    return t
end

# Compute peak curvature of a Hermite connector (in t-space) by sampling n_check
# uniformly-spaced t values.  At least 3 interior points are always included so
# that the mid-curve bulge is never missed.
function _hc_peak_curvature(a1, a2, a3; n_check::Int = 128)
    κ_max = 0.0
    for i in 0:n_check
        t  = i / n_check
        dp    = _hc_dp(a1, a2, a3, t)
        dp_v  = [dp[1], dp[2], dp[3]]
        dp_n  = norm(dp_v)
        dp_n < 1e-15 && continue
        ddp   = _hc_ddp(a2, a3, t)
        ddp_v = [ddp[1], ddp[2], ddp[3]]
        κ = norm(cross(dp_v, ddp_v)) / dp_n^3
        κ > κ_max && (κ_max = κ)
    end
    return κ_max
end

# Assemble Hermite coefficients from handle scale h and unit tangents.
function _hc_coeffs(p1_local, t_hat_out, h)
    v0  = h .* [0.0, 0.0, 1.0]
    v1  = h .* t_hat_out
    a1v = v0
    a2v = 3 .* p1_local .- 2 .* v0 .- v1
    a3v = 2 .* (.-p1_local) .+ v0 .+ v1
    a1 = (a1v[1], a1v[2], a1v[3])
    a2 = (a2v[1], a2v[2], a2v[3])
    a3 = (a3v[1], a3v[2], a3v[3])
    return a1, a2, a3
end

function _build_hermite_connector(p1_local::Vector{Float64}, t_hat_out::Vector{Float64};
                                   n_table::Int = 256,
                                   min_bend_radius::Union{Nothing,Float64} = nothing,
                                   meta = AbstractMeta[])
    chord = norm(p1_local)
    h = chord   # default: chord-proportioned handles

    # if !isnothing(min_bend_radius)
    #     κ_limit = 1.0 / min_bend_radius

    #     _κ(hv) = _hc_peak_curvature(_hc_coeffs(p1_local, t_hat_out, hv)...)

    #     # Scan κ(h) over an exponential ladder to find any h where κ(h) ≤ κ_limit.
    #     # κ(h) can be non-monotone for inflecting geometries (e.g. anti-parallel
    #     # tangents with transverse chord), so a simple two-point probe is insufficient.
    #     h0 = max(chord, 1e-12)
    #     h_scan = h0
    #     κ_prev = _κ(h_scan)
    #     h_hi_bracket = nothing   # smallest h where κ(h) ≤ κ_limit

    #     for _ in 1:64
    #         h_next = 2.0 * h_scan
    #         κ_next = _κ(h_next)
    #         if κ_next <= κ_limit
    #             h_hi_bracket = h_next
    #             break
    #         end
    #         h_scan = h_next
    #         κ_prev = κ_next
    #     end

    #     if isnothing(h_hi_bracket)
    #         throw(ArgumentError(
    #             "min_bend_radius=$(min_bend_radius) m is infeasible for this jump: " *
    #             "peak curvature could not be brought below $(round(κ_limit;digits=2)) m⁻¹ " *
    #             "within practical handle lengths (h_max=$(round(h_scan;digits=4)) m). " *
    #             "Geometry: chord=$(round(chord*1e3;digits=2)) mm, " *
    #             "outgoing tangent $(round.(t_hat_out;digits=3))."))
    #     end

    #     if _κ(h0) > κ_limit
    #         # Bisect in [h0, h_hi_bracket]: h0 violates, h_hi_bracket satisfies.
    #         # κ(h) may be non-monotone so this finds *a* valid h, not necessarily
    #         # the minimum; that is acceptable for the connector length optimisation.
    #         h_lo = h0
    #         h_hi = h_hi_bracket
    #         for _ in 1:64
    #             h_mid = (h_lo + h_hi) / 2
    #             if _κ(h_mid) > κ_limit
    #                 h_lo = h_mid
    #             else
    #                 h_hi = h_mid
    #             end
    #             h_hi - h_lo < 1e-10 * h_hi && break
    #         end
    #         h = h_hi
    #     end
    # end

    a1, a2, a3 = _hc_coeffs(p1_local, t_hat_out, h)
    a0 = (0.0, 0.0, 0.0)
    return HermiteConnector(a0, a1, a2, a3, _hc_build_table(a1, a2, a3, n_table);
                            meta = meta)
end

arc_length(seg::HermiteConnector) = seg.s_table[end]

function curvature(seg::HermiteConnector, s::Real)
    t      = _hc_t_from_s(seg, Float64(s))
    dp     = _hc_dp(seg.a1, seg.a2, seg.a3, t)
    dp_v   = [dp[1], dp[2], dp[3]]
    ddp    = _hc_ddp(seg.a2, seg.a3, t)
    ddp_v  = [ddp[1], ddp[2], ddp[3]]
    dp_n   = norm(dp_v)
    dp_n < 1e-15 && return 0.0
    return norm(cross(dp_v, ddp_v)) / dp_n^3
end

function geometric_torsion(seg::HermiteConnector, s::Real)
    t      = _hc_t_from_s(seg, Float64(s))
    dp     = _hc_dp(seg.a1, seg.a2, seg.a3, t)
    dp_v   = [dp[1], dp[2], dp[3]]
    ddp    = _hc_ddp(seg.a2, seg.a3, t)
    ddp_v  = [ddp[1], ddp[2], ddp[3]]
    dddp_v = [6*seg.a3[1], 6*seg.a3[2], 6*seg.a3[3]]
    c      = cross(dp_v, ddp_v)
    denom  = dot(c, c)
    denom < 1e-30 && return 0.0
    return dot(c, dddp_v) / denom
end

function position_local(seg::HermiteConnector, s::Real)
    t = _hc_t_from_s(seg, Float64(s))
    a0, a1, a2, a3 = seg.a0, seg.a1, seg.a2, seg.a3
    t2, t3 = t*t, t*t*t
    return [a0[i] + a1[i]*t + a2[i]*t2 + a3[i]*t3 for i in 1:3]
end

function tangent_local(seg::HermiteConnector, s::Real)
    t    = _hc_t_from_s(seg, Float64(s))
    dp   = _hc_dp(seg.a1, seg.a2, seg.a3, t)
    dp_v = [dp[1], dp[2], dp[3]]
    spd  = norm(dp_v)
    spd < 1e-15 && return [0.0, 0.0, 1.0]
    return dp_v ./ spd
end

function normal_local(seg::HermiteConnector, s::Real)
    t    = _hc_t_from_s(seg, Float64(s))
    dp   = _hc_dp(seg.a1, seg.a2, seg.a3, t)
    dp_v = [dp[1], dp[2], dp[3]]
    spd  = norm(dp_v)
    T    = spd > 1e-15 ? dp_v ./ spd : [0.0, 0.0, 1.0]
    ddp  = _hc_ddp(seg.a2, seg.a3, t)
    acc  = [ddp[1], ddp[2], ddp[3]] .- dot([ddp[1], ddp[2], ddp[3]], T) .* T
    an   = norm(acc)
    return an > 1e-15 ? acc ./ an : _perp_unit(T)
end

function binormal_local(seg::HermiteConnector, s::Real)
    return cross(tangent_local(seg, s), normal_local(seg, s))
end

function end_position_local(seg::HermiteConnector)
    a0, a1, a2, a3 = seg.a0, seg.a1, seg.a2, seg.a3
    return [a0[i] + a1[i] + a2[i] + a3[i] for i in 1:3]
end

function end_frame_local(seg::HermiteConnector)
    # t = 1.0 exactly at the endpoint; skip the arc-length table lookup.
    dp  = _hc_dp(seg.a1, seg.a2, seg.a3, 1.0)
    dp_v = [dp[1], dp[2], dp[3]]
    spd  = norm(dp_v)
    T    = spd > 1e-15 ? dp_v ./ spd : [0.0, 0.0, 1.0]
    ddp  = _hc_ddp(seg.a2, seg.a3, 1.0)
    acc  = [ddp[1], ddp[2], ddp[3]] .- dot([ddp[1], ddp[2], ddp[3]], T) .* T
    an   = norm(acc)
    N    = an > 1e-15 ? acc ./ an : _perp_unit(T)
    return (T, N, cross(T, N))
end

# Resolve JumpBy / JumpTo to a HermiteConnector at build() time.

function _resolve_at_placement(seg::JumpBy, pos::Vector{Float64}, frame_mat::Matrix{Float64})
    p1_local  = collect(seg.delta)
    chord     = norm(p1_local)
    t_hat_out = isnothing(seg.tangent_out) ?
        (chord > 1e-15 ? p1_local ./ chord : [0.0, 0.0, 1.0]) :
        normalize(collect(seg.tangent_out))
    return _build_hermite_connector(p1_local, t_hat_out;
                                    min_bend_radius = seg.min_bend_radius,
                                    meta = seg.meta)
end

function _resolve_at_placement(seg::JumpTo, pos::Vector{Float64}, frame_mat::Matrix{Float64})
    p1_local  = frame_mat' * (collect(seg.destination) .- pos)
    chord     = norm(p1_local)
    t_hat_out = if isnothing(seg.tangent_out)
        chord > 1e-15 ? p1_local ./ chord : [0.0, 0.0, 1.0]
    else
        normalize(frame_mat' * normalize(collect(seg.tangent_out)))
    end
    return _build_hermite_connector(p1_local, t_hat_out;
                                    min_bend_radius = seg.min_bend_radius,
                                    meta = seg.meta)
end

_resolve_at_placement(seg::AbstractPathSegment, ::AbstractVector, ::AbstractMatrix) = seg

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
# PathSpecBuilder (mutable authoring) → PathSpec (immutable spec)
# -----------------------------------------------------------------------

mutable struct PathSpecBuilder
    segments::Vector{AbstractPathSegment}
    s_start::Real
    PathSpecBuilder() = new(AbstractPathSegment[], 0.0)
end

struct PathSpec
    segments::Vector{AbstractPathSegment}
    s_start::Real
end

"""
    freeze(builder::PathSpecBuilder) → PathSpec

Snapshot a builder into an immutable `PathSpec`. The segment vector is
deep-copied so subsequent mutation of `builder` does not affect the spec.
"""
freeze(b::PathSpecBuilder) = PathSpec(deepcopy(b.segments), b.s_start)

function straight!(spec::PathSpecBuilder; length, meta = AbstractMeta[])
    push!(spec.segments, StraightSegment(length; meta))
    return spec
end

function bend!(spec::PathSpecBuilder; radius::Real, angle::Real, axis_angle::Real = 0.0,
               meta = AbstractMeta[])
    push!(spec.segments, BendSegment(radius, angle, axis_angle; meta))
    return spec
end

function helix!(spec::PathSpecBuilder; radius::Real, pitch::Real, turns::Real,
                axis_angle::Real = 0.0, meta = AbstractMeta[])
    push!(spec.segments, HelixSegment(radius, pitch, turns, axis_angle; meta))
    return spec
end

function catenary!(spec::PathSpecBuilder; a::Real, length::Real, axis_angle::Real = 0.0,
                   meta = AbstractMeta[])
    push!(spec.segments, CatenarySegment(a, length, axis_angle; meta))
    return spec
end

function jumpby!(spec::PathSpecBuilder; delta, tangent = nothing, min_bend_radius = nothing,
                 meta = AbstractMeta[])
    push!(spec.segments, JumpBy(delta; tangent_out = tangent, min_bend_radius, meta))
    return spec
end

function jumpto!(spec::PathSpecBuilder; destination, tangent = nothing, min_bend_radius = nothing,
                 meta = AbstractMeta[])
    push!(spec.segments, JumpTo(destination; tangent_out = tangent, min_bend_radius, meta))
    return spec
end

# -----------------------------------------------------------------------
# PathSpecCached  (immutable derived layout)
# -----------------------------------------------------------------------

struct PathSpecCached
    spec::PathSpec
    placed_segments::Vector{PlacedSegment}
    s_end::Real
    resolved_twists::Vector{ResolvedTwistRate}
end

# Convenience accessor: callers historically read `path.spec.s_start`. Provided
# as a function so external code can be retyped without reaching into
# `path.spec.s_start` everywhere.
s_start(path::PathSpecCached) = path.spec.s_start

# -----------------------------------------------------------------------
# build()
# -----------------------------------------------------------------------

"""
    build(builder_or_spec) → PathSpecCached

Compile a `PathSpecBuilder` or `PathSpec` into an immutable `PathSpecCached`.
"""
_safe_normalize(v::AbstractVector) = v ./ sqrt(sum(abs2, v))

build(b::PathSpecBuilder) = build(freeze(b))

function build(spec::PathSpec)
    pos     = zeros(3)
    N_frame = [1.0, 0.0, 0.0]
    B_frame = [0.0, 1.0, 0.0]
    T_frame = [0.0, 0.0, 1.0]

    isempty(spec.segments) && error("build: PathSpec contains no segments")

    s_eff  = Float64(spec.s_start)
    placed = PlacedSegment[]

    for seg_orig in spec.segments
        frame      = hcat(N_frame, B_frame, T_frame)   # columns: [N | B | T]
        seg_placed = _resolve_at_placement(seg_orig, pos, frame)
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

        s_eff += arc_length(seg_placed)
    end

    # Skip twist resolution when no segment carries a Twist. This also keeps
    # MCM-valued (Particles) paths working: _resolve_twists requires a Float64
    # s_end which is incompatible with Particles-typed cumulative arc length.
    has_twist = any(ps -> any(m -> m isa Twist, segment_meta(ps.segment)), placed)
    resolved_twists = has_twist ? _resolve_twists(placed, Float64(s_eff)) : ResolvedTwistRate[]
    return PathSpecCached(spec, placed, s_eff, resolved_twists)
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
            push!(anchors, (Float64(ps.s_offset_eff), twists_here[1]))
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

function _find_placed_segment(path::PathSpecCached, s)
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
    error("s = $s out of path bounds [$(path.spec.s_start), $(path.s_end)]")
end

function _local_to_global(ps::PlacedSegment, v_local::AbstractVector)
    return ps.frame * v_local
end

# -----------------------------------------------------------------------
# Differential geometry interface on Path
# -----------------------------------------------------------------------

arc_length(path::PathSpecCached) = path.s_end - path.spec.s_start

function arc_length(::PathSpecCached, s1, s2)
    @assert s2 >= s1 "arc_length: require s2 >= s1"
    return s2 - s1
end

function curvature(path::PathSpecCached, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return curvature(ps.segment, s_local)
end

function geometric_torsion(path::PathSpecCached, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return geometric_torsion(ps.segment, s_local)
end

"""
    material_twist(path, s)

Material twist rate (rad/m) at effective arc length `s`, summed over all
resolved twist runs that contain `s`. Runs are disjoint by construction; the
sum exists only as a robustness guard.
"""
function material_twist(path::PathSpecCached, s)
    τ = zero(s isa AbstractFloat ? s : Float64(s))
    for r in path.resolved_twists
        if r.s_eff_start <= s <= r.s_eff_end
            τ += r.rate isa Function ? r.rate(s - r.s_eff_start) : r.rate
        end
    end
    return τ
end

function position(path::PathSpecCached, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return ps.origin + _local_to_global(ps, position_local(ps.segment, s_local))
end

function tangent(path::PathSpecCached, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return _local_to_global(ps, tangent_local(ps.segment, s_local))
end

function normal(path::PathSpecCached, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return _local_to_global(ps, normal_local(ps.segment, s_local))
end

function binormal(path::PathSpecCached, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return _local_to_global(ps, binormal_local(ps.segment, s_local))
end

function frame(path::PathSpecCached, s::Real)
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

start_point(path::PathSpecCached)   = position(path, path.spec.s_start)
end_point(path::PathSpecCached)     = position(path, path.s_end)
start_tangent(path::PathSpecCached) = tangent(path, path.spec.s_start)
end_tangent(path::PathSpecCached)   = tangent(path, path.s_end)

# -----------------------------------------------------------------------
# Path measures
# -----------------------------------------------------------------------

path_length(path::PathSpecCached) = arc_length(path)

function cartesian_distance(path::PathSpecCached, s1::Real, s2::Real)
    return norm(position(path, s2) - position(path, s1))
end

function bounding_box(path::PathSpecCached; n::Int = 512)
    ss = range(path.spec.s_start, path.s_end; length = n)
    pts = [position(path, s) for s in ss]
    lo = minimum(reduce(hcat, pts); dims = 2) |> vec
    hi = maximum(reduce(hcat, pts); dims = 2) |> vec
    return (; lo, hi)
end

function total_turning_angle(path::PathSpecCached)
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

function total_torsion(path::PathSpecCached)
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
`s_start` to `s_end` (defaults: full path).  Both
endpoints must lie in `[path.spec.s_start, path.s_end]`.

Constant twist rates integrate analytically; `Function` rates use adaptive
Gauss–Kronrod (QuadGK) at the supplied `rtol`/`atol`.
"""
function total_material_twist(
    path::PathSpecCached;
    s_start::Real = path.spec.s_start,
    s_end::Real   = path.s_end,
    rtol::Real    = 1e-8,
    atol::Real    = 0.0,
)
    s_lo = Float64(s_start)
    s_hi = Float64(s_end)
    if s_lo > s_hi
        throw(ArgumentError(
            "total_material_twist: require s_start ≤ s_end; got s_start=$(s_lo), s_end=$(s_hi)"))
    end
    ps0 = Float64(path.spec.s_start)
    ps1 = Float64(path.s_end)
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
    path::PathSpecCached;
    s_start::Real = path.spec.s_start,
    s_end::Real   = path.s_end,
    rtol::Real    = 1e-8,
    atol::Real    = 0.0,
)
    s_lo = Float64(s_start)
    s_hi = Float64(s_end)
    if s_lo > s_hi
        throw(ArgumentError(
            "total_frame_rotation: require s_start ≤ s_end; got s_start=$(s_lo), s_end=$(s_hi)"))
    end
    ps0 = Float64(path.spec.s_start)
    ps1 = Float64(path.s_end)
    if !(ps0 - 1e-12 <= s_lo <= ps1 + 1e-12) || !(ps0 - 1e-12 <= s_hi <= ps1 + 1e-12)
        throw(ArgumentError(
            "total_frame_rotation: require path.s_start ≤ s ≤ path.s_end for both endpoints; " *
            "got [$(s_lo), $(s_hi)] m vs path domain [$(ps0), $(ps1)] m"))
    end
    s_lo == s_hi && return 0.0

    # Geometric torsion: integrate segment-by-segment, with closed-form
    # shortcuts where available.
    τ_total = 0.0
    rtolf = Float64(rtol)
    atolf = Float64(atol)
    for ps in path.placed_segments
        seg = ps.segment
        seg_s_start = Float64(ps.s_offset_eff)
        seg_s_end   = seg_s_start + Float64(arc_length(seg))
        a = max(seg_s_start, s_lo)
        b = min(seg_s_end,   s_hi)
        b <= a && continue
        a_loc = a - seg_s_start
        b_loc = b - seg_s_start
        if seg isa StraightSegment || seg isa BendSegment || seg isa CatenarySegment
            # τ_geom = 0
        elseif seg isa HelixSegment
            τ_total += geometric_torsion(seg, 0.0) * (b_loc - a_loc)
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
function writhe(path::PathSpecCached; n::Int = 256)
    ss = collect(range(path.spec.s_start, path.s_end; length = n))
    rs = [position(path, s) for s in ss]
    ts = [tangent(path, s)  for s in ss]
    ds = (path.s_end - path.spec.s_start) / (n - 1)

    Wr = 0.0
    for i in 1:n, j in 1:n
        i == j && continue
        r_ij = rs[i] - rs[j]
        d = norm(r_ij)
        d < 1e-14 * (path.s_end - path.spec.s_start) && continue
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

function _segment_total_angle(seg::HermiteConnector)
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
    path::PathSpecCached,
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
function sample_path(path::PathSpecCached, s1::Real, s2::Real; fidelity::Float64 = 1.0)
    @assert s2 > s1    "sample_path: require s2 > s1"
    @assert fidelity > 0.0 "sample_path: fidelity must be positive"

    s_lo = Float64(s1)
    s_hi = Float64(s2)

    # Collect arc-length stations segment by segment, sharing junction points.
    all_s = Float64[]
    for ps in path.placed_segments
        seg_s_start = ps.s_offset_eff
        seg_s_end   = ps.s_offset_eff + arc_length(ps.segment)

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

function path_segment_breakpoints(path::PathSpecCached)
    points = Real[path.spec.s_start]
    for ps in path.placed_segments
        push!(points, ps.s_offset_eff)
        push!(points, ps.s_offset_eff + arc_length(ps.segment))
    end
    push!(points, path.s_end)
    return normalize_breakpoints(points)
end

function path_twist_breakpoints(path::PathSpecCached)
    points = Real[path.spec.s_start, path.s_end]
    for r in path.resolved_twists
        push!(points, r.s_eff_start)
        push!(points, r.s_eff_end)
    end
    return normalize_breakpoints(points)
end

function breakpoints(path::PathSpecCached)
    return normalize_breakpoints(vcat(path_segment_breakpoints(path), path_twist_breakpoints(path)))
end

function sample(path::PathSpecCached, s_values)
    return [frame(path, s) for s in s_values]
end

function sample_uniform(path::PathSpecCached; n::Int = 256)
    ss = range(path.spec.s_start, path.s_end; length = n)
    return sample(path, ss)
end
