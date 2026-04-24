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
    total_material_twist(path; s_start, s_end, n_quad)
    total_frame_rotation(path; s_start, s_end, n_quad)
    writhe(path)
    sample(path, s_values)
    sample_uniform(path; n)
"""

using LinearAlgebra

# -----------------------------------------------------------------------
# Abstract segment type
# -----------------------------------------------------------------------

abstract type AbstractPathSegment end

# Required interface for each concrete segment (local arc-length s ∈ [0, arc_length(seg)]).
# Return element type T follows the segment's own type parameter (Float64 for
# deterministic segments, Particles for MCM-valued fields):
#   arc_length(seg)                   → T
#   nominal_arc_length(seg)           → T          (arc length at shrinkage = 1)
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

segment_shrinkage(seg::AbstractPathSegment) = seg.shrinkage
segment_nickname(seg::AbstractPathSegment)  = :nickname ∈ fieldnames(typeof(seg)) ? seg.nickname : nothing

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

struct StraightSegment{T} <: AbstractPathSegment
    length::T
    shrinkage::T
    nickname::Union{Nothing,String}
end

function StraightSegment(length; shrinkage = 1.0, nickname = nothing)
    L, α = promote(length, shrinkage)
    StraightSegment{typeof(L)}(L, α, isnothing(nickname) ? nothing : String(nickname))
end

arc_length(seg::StraightSegment)         = seg.length * seg.shrinkage
nominal_arc_length(seg::StraightSegment) = seg.length
curvature(seg::StraightSegment, _)       = zero(seg.length)
geometric_torsion(seg::StraightSegment, _) = zero(seg.length)

position_local(seg::StraightSegment, s)   = [zero(s), zero(s), s]
tangent_local(seg::StraightSegment, _)    = [zero(seg.length), zero(seg.length), one(seg.length)]
normal_local(seg::StraightSegment, _)     = [one(seg.length), zero(seg.length), zero(seg.length)]
binormal_local(seg::StraightSegment, _)   = [zero(seg.length), one(seg.length), zero(seg.length)]
end_position_local(seg::StraightSegment)  = [zero(seg.length), zero(seg.length), arc_length(seg)]
end_frame_local(seg::StraightSegment) = ([zero(seg.length), zero(seg.length), one(seg.length)],
                                          [one(seg.length), zero(seg.length), zero(seg.length)],
                                          [zero(seg.length), one(seg.length), zero(seg.length)])

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
struct BendSegment{T} <: AbstractPathSegment
    radius::T
    angle::T       # total angle swept (rad), preserved under shrinkage
    axis_angle::T  # orientation of inward normal in transverse plane (rad)
    shrinkage::T
    nickname::Union{Nothing,String}

    function BendSegment(radius, angle, axis_angle = 0.0;
                         shrinkage = 1.0, nickname = nothing)
        @assert radius > 0 "BendSegment: radius must be positive"
        r, a, x, α = promote(radius, angle, axis_angle, shrinkage)
        new{typeof(r)}(r, a, x, α, isnothing(nickname) ? nothing : String(nickname))
    end
end

arc_length(seg::BendSegment)         = seg.shrinkage * seg.radius * abs(seg.angle)
nominal_arc_length(seg::BendSegment) = seg.radius * abs(seg.angle)
curvature(seg::BendSegment, _)       = one(seg.radius) / (seg.shrinkage * seg.radius)
geometric_torsion(seg::BendSegment, _) = zero(seg.radius)

function position_local(seg::BendSegment, s)
    R   = seg.shrinkage * seg.radius
    θ   = s / R
    φ   = seg.axis_angle
    n̂   = [cos(φ), sin(φ), zero(φ)]
    return R * (1 - cos(θ)) * n̂ + [zero(R), zero(R), R * sin(θ)]
end

function tangent_local(seg::BendSegment, s)
    R = seg.shrinkage * seg.radius
    θ = s / R
    φ = seg.axis_angle
    return [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)]
end

function normal_local(seg::BendSegment, s)
    R = seg.shrinkage * seg.radius
    θ = s / R
    φ = seg.axis_angle
    return [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)]
end

function binormal_local(seg::BendSegment, _)
    φ = seg.axis_angle
    return [-sin(φ), cos(φ), zero(φ)]   # constant for circular arc (zero torsion)
end

function end_position_local(seg::BendSegment)
    R = seg.shrinkage * seg.radius
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
    CatenarySegment(a, length, axis_angle; shrinkage)

A catenary curve in the plane whose horizontal direction is at `axis_angle`
from the local N-axis.  `a` (m) is the catenary parameter (a = T₀/(ρg),
the ratio of horizontal tension to weight per unit length).  The curve
starts with tangent along the local z-axis (vertical catenary vertex) and
curves horizontally.

Under shrinkage α: effective parameter a → α·a, arc length → α·length,
curvature κ(s) = a_eff / (a_eff² + s²).  Angles are preserved.
"""
struct CatenarySegment{T} <: AbstractPathSegment
    a::T
    length::T
    axis_angle::T
    shrinkage::T
    nickname::Union{Nothing,String}

    function CatenarySegment(a, length, axis_angle = 0.0;
                             shrinkage = 1.0, nickname = nothing)
        @assert a > 0      "CatenarySegment: a must be positive"
        @assert length > 0 "CatenarySegment: length must be positive"
        av, L, x, α = promote(a, length, axis_angle, shrinkage)
        new{typeof(av)}(av, L, x, α, isnothing(nickname) ? nothing : String(nickname))
    end
end

arc_length(seg::CatenarySegment)         = seg.length * seg.shrinkage
nominal_arc_length(seg::CatenarySegment) = seg.length
geometric_torsion(seg::CatenarySegment, _) = zero(seg.a)

function curvature(seg::CatenarySegment, s)
    a = seg.a * seg.shrinkage
    return a / (a^2 + s^2)
end

function position_local(seg::CatenarySegment, s)
    a = seg.a * seg.shrinkage
    φ = seg.axis_angle
    n̂ = [cos(φ), sin(φ), zero(φ)]
    horiz = a * (sqrt(1 + (s / a)^2) - 1)
    vert  = a * asinh(s / a)
    return horiz * n̂ + [zero(a), zero(a), vert]
end

function tangent_local(seg::CatenarySegment, s)
    a = seg.a * seg.shrinkage
    φ = seg.axis_angle
    q = sqrt(1 + (s / a)^2)
    return [(s / a) / q * cos(φ), (s / a) / q * sin(φ), one(q) / q]
end

function normal_local(seg::CatenarySegment, s)
    # N = dT/ds / |dT/ds|, derived analytically: N = [n̂_horiz/q, -s/a/q] normalised
    a = seg.a * seg.shrinkage
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
    HelixSegment(radius, pitch, turns, axis_angle; shrinkage)

A helix whose entry tangent is ẑ (the incoming sliding-frame tangent), ensuring
continuity with the prior segment.  `axis_angle` (rad) selects which transverse
direction n̂ = cos(axis_angle)·x̂ + sin(axis_angle)·ŷ the helix curves toward.

The helix axis â = (h·ẑ + R·n̂) / ℓ' is tilted from the transverse plane by
arctan(h/R) toward ẑ, where h = pitch/(2π) and ℓ' = √(R²+h²).  This tilt is
the geometric consequence of demanding tangent(0) = ẑ: a zero-pitch helix
reduces to a circular arc (BendSegment) in the n̂ direction.

    κ(s) = R / (R² + h²)          (constant)
    τ_geom(s) = h / (R² + h²)     (constant)
    arc_length = turns · 2π · √(R² + h²) · shrinkage

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
    shrinkage::T
    nickname::Union{Nothing,String}

    function HelixSegment(radius, pitch, turns,
                          axis_angle = 0.0; shrinkage = 1.0, nickname = nothing)
        @assert radius > 0 "HelixSegment: radius must be positive"
        @assert turns  > 0 "HelixSegment: turns must be positive"
        r, p, n, x, α = promote(radius, pitch, turns, axis_angle, shrinkage)
        new{typeof(r)}(r, p, n, x, α, isnothing(nickname) ? nothing : String(nickname))
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

function curvature(seg::HelixSegment, _)
    R = seg.radius * seg.shrinkage
    h = _helix_h(seg) * seg.shrinkage
    return R / (R^2 + h^2)
end

function geometric_torsion(seg::HelixSegment, _)
    R = seg.radius * seg.shrinkage
    h = _helix_h(seg) * seg.shrinkage
    return h / (R^2 + h^2)
end

function _helix_basis(seg::HelixSegment)
    R  = seg.radius * seg.shrinkage
    h  = _helix_h(seg) * seg.shrinkage
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
    JumpBy(delta, tangent_out, shrinkage, min_bend_radius)

Connects the current position to current_position + shrinkage·delta using an
smooth function.  The incoming tangent is the current sliding frame
tangent.  `tangent_out` is the desired outgoing tangent direction; if nothing,
the connector computes an implicit tangent that minimizes curvature variation
(requires a nonlinear solve).

`min_bend_radius` (metres, nothing = unconstrained) sets a lower bound on the
radius of curvature of the Hermite connector.  The tangent handle length is
extended beyond the chord default when necessary to keep κ ≤ 1/R_min.
Currently a stub. 

"""
struct JumpBy <: AbstractPathSegment
    delta::NTuple{3, Float64}
    tangent_out::Union{Nothing, NTuple{3, Float64}}
    shrinkage::Float64
    min_bend_radius::Union{Nothing, Float64}
end

function JumpBy(delta; tangent_out = nothing, shrinkage::Real = 1.0,
                min_bend_radius = nothing)
    d = (Float64(delta[1]), Float64(delta[2]), Float64(delta[3]))
    t = isnothing(tangent_out) ? nothing :
        (Float64(tangent_out[1]), Float64(tangent_out[2]), Float64(tangent_out[3]))
    r = isnothing(min_bend_radius) ? nothing : Float64(min_bend_radius)
    JumpBy(d, t, Float64(shrinkage), r)
end

"""
    JumpTo(destination, tangent_out, shrinkage, min_bend_radius)

Connects the current position to the fixed lab-frame `destination` using smooth
function.  Shrinkage changes the arc length of the connector (the
connector absorbs the geometry change) but does not move the destination.

`min_bend_radius` (metres, nothing = unconstrained) sets a lower bound on the
radius of curvature of the Hermite connector.  The tangent handle length is
extended beyond the chord default when necessary to keep κ ≤ 1/R_min.
Currently a stub. 
"""
struct JumpTo <: AbstractPathSegment
    destination::NTuple{3, Float64}
    tangent_out::Union{Nothing, NTuple{3, Float64}}
    shrinkage::Float64
    min_bend_radius::Union{Nothing, Float64}
end

function JumpTo(destination; tangent_out = nothing, shrinkage::Real = 1.0,
                min_bend_radius = nothing)
    d = (Float64(destination[1]), Float64(destination[2]), Float64(destination[3]))
    t = isnothing(tangent_out) ? nothing :
        (Float64(tangent_out[1]), Float64(tangent_out[2]), Float64(tangent_out[3]))
    r = isnothing(min_bend_radius) ? nothing : Float64(min_bend_radius)
    JumpTo(d, t, Float64(shrinkage), r)
end

# JumpBy and JumpTo are context-dependent: geometry is resolved at build() time
# into a HermiteConnector.  Calling these methods on the raw structs is unsupported.
for T in (JumpBy, JumpTo)
    @eval begin
        arc_length(::$T)             = error($(string(T)) * ": call build() to resolve jump geometry")
        nominal_arc_length(::$T)     = error($(string(T)) * ": call build() to resolve jump geometry")
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
    arithmetic) are not lifted to `Particles`. A `PathSpec` that mixes an
    uncertain `BendSegment`/`HelixSegment`/etc. with a `jumpby!`/`jumpto!`
    will fail at `build()`.
"""
struct HermiteConnector <: AbstractPathSegment
    a0        :: NTuple{3,Float64}   # P(t) = a0 + a1 t + a2 t² + a3 t³, t ∈ [0,1]
    a1        :: NTuple{3,Float64}
    a2        :: NTuple{3,Float64}
    a3        :: NTuple{3,Float64}
    s_table   :: Vector{Float64}     # s_table[i] = arc-length at t = (i-1)/(n-1)
    shrinkage :: Float64             # always 1.0; stored for interface uniformity
end

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
                                   min_bend_radius::Union{Nothing,Float64} = nothing)
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
    return HermiteConnector(a0, a1, a2, a3, _hc_build_table(a1, a2, a3, n_table), 1.0)
end

arc_length(seg::HermiteConnector)         = seg.s_table[end]
nominal_arc_length(seg::HermiteConnector) = arc_length(seg)

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
    p1_local  = collect(seg.delta) .* seg.shrinkage
    chord     = norm(p1_local)
    t_hat_out = isnothing(seg.tangent_out) ?
        (chord > 1e-15 ? p1_local ./ chord : [0.0, 0.0, 1.0]) :
        normalize(collect(seg.tangent_out))
    return _build_hermite_connector(p1_local, t_hat_out;
                                    min_bend_radius = seg.min_bend_radius)
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
                                    min_bend_radius = seg.min_bend_radius)
end

_resolve_at_placement(seg::AbstractPathSegment, ::AbstractVector, ::AbstractMatrix) = seg

# -----------------------------------------------------------------------
# PathSpec  (mutable authoring struct)
# -----------------------------------------------------------------------

mutable struct PathSpec
    segments::Vector{AbstractPathSegment}
    twist_overlays::Vector{TwistOverlay}
    PathSpec() = new(AbstractPathSegment[], TwistOverlay[])
end

function straight!(spec::PathSpec; length::Real, shrinkage::Real = 1.0, nickname = nothing)
    push!(spec.segments, StraightSegment(length; shrinkage, nickname))
    return spec
end

function bend!(spec::PathSpec; radius::Real, angle::Real, axis_angle::Real = 0.0,
               shrinkage::Real = 1.0, nickname = nothing)
    push!(spec.segments, BendSegment(radius, angle, axis_angle; shrinkage, nickname))
    return spec
end

function helix!(spec::PathSpec; radius::Real, pitch::Real, turns::Real,
                axis_angle::Real = 0.0, shrinkage::Real = 1.0, nickname = nothing)
    push!(spec.segments, HelixSegment(radius, pitch, turns, axis_angle; shrinkage, nickname))
    return spec
end

function catenary!(spec::PathSpec; a::Real, length::Real, axis_angle::Real = 0.0,
                   shrinkage::Real = 1.0, nickname = nothing)
    push!(spec.segments, CatenarySegment(a, length, axis_angle; shrinkage, nickname))
    return spec
end

function jumpby!(spec::PathSpec; delta, tangent = nothing, shrinkage::Real = 1.0,
                 min_bend_radius = nothing)
    push!(spec.segments, JumpBy(delta; tangent_out = tangent, shrinkage, min_bend_radius))
    return spec
end

function jumpto!(spec::PathSpec; destination, tangent = nothing, shrinkage::Real = 1.0,
                 min_bend_radius = nothing)
    push!(spec.segments, JumpTo(destination; tangent_out = tangent, shrinkage, min_bend_radius))
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
    s_offset_eff::Real          # cumulative effective arc-length at segment start
    s_offset_nom::Real          # cumulative nominal arc-length at segment start
    origin::AbstractVector      # global start position (length 3)
    frame::AbstractMatrix       # 3×3, columns [N_global | B_global | T_global]
                                # transforms local vectors → global: v_g = frame * v_l
end

struct ResolvedTwistRate
    s_eff_start::Real
    s_eff_end::Real
    rate::Function                     # τ_mat(s_eff) in rad/m
end

struct ResolvedTwistOverlay
    rates::Vector{ResolvedTwistRate}
end

# -----------------------------------------------------------------------
# Path  (immutable built struct)
# -----------------------------------------------------------------------

struct Path
    s_start::Real
    s_end::Real
    placed_segments::Vector{PlacedSegment}
    resolved_overlays::Vector{ResolvedTwistOverlay}
end

# -----------------------------------------------------------------------
# build()
# -----------------------------------------------------------------------

function _apply_shrinkage_override(seg::S, override) where {S <: AbstractPathSegment}
    isnothing(override) && return seg
    # shrinkage and nickname are keyword arguments in the outer constructors.
    excluded = (:shrinkage, :nickname)
    flds = fieldnames(S)
    vals = [getfield(seg, f) for f in flds if f ∉ excluded]
    nick = :nickname ∈ flds ? getfield(seg, :nickname) : nothing
    # Use the UnionAll (unparameterized) name so the outer keyword constructor dispatches.
    return S.name.wrapper(vals...; shrinkage = override, nickname = nick)
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
_safe_normalize(v::AbstractVector) = v ./ sqrt(sum(abs2, v))

function build(spec::PathSpec; shrinkage = nothing)
    pos     = zeros(3)
    N_frame = [1.0, 0.0, 0.0]
    B_frame = [0.0, 1.0, 0.0]
    T_frame = [0.0, 0.0, 1.0]

    isempty(spec.segments) && error("build: PathSpec contains no segments")

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

        frame      = hcat(N_frame, B_frame, T_frame)   # columns: [N | B | T]
        seg_placed = _resolve_at_placement(seg, pos, frame)
        push!(placed, PlacedSegment(seg_placed, s_eff, s_nom, copy(pos), copy(frame)))

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
        s_nom += nominal_arc_length(seg_placed)
    end

    resolved = [_resolve_overlay(ov, placed) for ov in spec.twist_overlays]

    return Path(0.0, s_eff, placed, resolved)
end

# -----------------------------------------------------------------------
# Segment lookup helpers
# -----------------------------------------------------------------------

function _find_placed_segment(path::Path, s)
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
    error("s = $s out of path bounds [$(path.s_start), $(path.s_end)]")
end

function _local_to_global(ps::PlacedSegment, v_local::AbstractVector)
    return ps.frame * v_local
end

# -----------------------------------------------------------------------
# Differential geometry interface on Path
# -----------------------------------------------------------------------

arc_length(path::Path) = path.s_end - path.s_start

function arc_length(::Path, s1, s2)
    @assert s2 >= s1 "arc_length: require s2 >= s1"
    return s2 - s1
end

function curvature(path::Path, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return curvature(ps.segment, s_local)
end

function geometric_torsion(path::Path, s::Real)
    ps, s_local = _find_placed_segment(path, s)
    return geometric_torsion(ps.segment, s_local)
end

"""
    material_twist

Material twist rate at effective arc-length `s`.
"""
function material_twist(path::Path, s)
    τ = zero(s)
    for ov in path.resolved_overlays
        for r in ov.rates
            if r.s_eff_start <= s <= r.s_eff_end
                τ += r.rate(s)
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

"""
    total_material_twist(path; s_start=path.s_start, s_end=path.s_end, n_quad=128) → Float64

Integrated material twist ``∫ τ_{\\mathrm{mat}}(s) \\, ds`` over effective arc length from
`s_start` to `s_end` (defaults: full path). Require `s_start ≤ s_end`; if `s_start > s_end`,
an `ArgumentError` is thrown (no reordering). If `s_start == s_end`, returns `0.0`. Both
endpoints must lie in `[path.s_start, path.s_end]`. Quadrature matches the previous rule on
each overlap of `[s_start, s_end]` with resolved twist overlays.
"""
function total_material_twist(
    path::Path;
    s_start::Real = path.s_start,
    s_end::Real = path.s_end,
    n_quad::Int = 128,
)
    s_lo = min(Float64(s_start), Float64(s_end))
    s_hi = max(Float64(s_start), Float64(s_end))
    ps0 = path.s_start
    ps1 = path.s_end
    if !(ps0 <= s_lo <= ps1) || !(ps0 <= s_hi <= ps1)
        throw(ArgumentError(
            "total_material_twist: require path.s_start ≤ s ≤ path.s_end for both endpoints; " *
            "got [$(s_lo), $(s_hi)] m vs path domain [$(ps0), $(ps1)] m",
        ))
    end
    s_hi <= s_lo && return 0.0

    total = 0.0
    for ov in path.resolved_overlays
        for r in ov.rates
            a = max(r.s_eff_start, s_lo)
            b = min(r.s_eff_end, s_hi)
            b <= a && continue
            ss = range(a, b; length = n_quad + 1)
            h = (b - a) / n_quad
            total += h * (sum(r.rate(s) for s in ss) - r.rate(a)/2 - r.rate(b)/2)
        end
    end
    return total
end

"""
    total_frame_rotation(path; s_start=path.s_start, s_end=path.s_end, n_quad=128) → Float64

Total rotation of the polarization reference frame over effective arc length from `s_start`
to `s_end`, integrating both contributions:

    dψ/ds = τ_geom(s) + Ω_material(s)

where `τ_geom` is the geometric torsion of the centerline (nonzero for helices; zero for
straight segments and circular bends) and `Ω_material` is the applied material twist rate
from `TwistOverlay`s.  Returns the integral in radians.

See `path-geometry.md` for a discussion of how this differs from `total_torsion` and
`total_material_twist` individually.
"""
function total_frame_rotation(
    path::Path;
    s_start::Real = path.s_start,
    s_end::Real   = path.s_end,
    n_quad::Int   = 128,
)
    s_lo = min(Float64(s_start), Float64(s_end))
    s_hi = max(Float64(s_start), Float64(s_end))
    ps0 = path.s_start
    ps1 = path.s_end
    if !(ps0 <= s_lo <= ps1) || !(ps0 <= s_hi <= ps1)
        throw(ArgumentError(
            "total_frame_rotation: require path.s_start ≤ s ≤ path.s_end for both endpoints; " *
            "got [$(s_lo), $(s_hi)] m vs path domain [$(ps0), $(ps1)] m",
        ))
    end
    s_hi <= s_lo && return 0.0

    # Geometric torsion: integrate segment-by-segment over the requested window.
    τ_total = 0.0
    s_seg_start = ps0
    for ps in path.placed_segments
        seg = ps.segment
        s_seg_end = s_seg_start + arc_length(seg)
        a = max(s_seg_start, s_lo)
        b = min(s_seg_end,   s_hi)
        if b > a
            # Convert global s to local s within the segment.
            a_loc = a - s_seg_start
            b_loc = b - s_seg_start
            if seg isa StraightSegment || seg isa BendSegment
                # τ_geom = 0 exactly; skip.
            elseif seg isa HelixSegment
                τ_total += geometric_torsion(seg, 0.0) * (b_loc - a_loc)
            else
                ss = range(a_loc, b_loc; length = n_quad + 1)
                h  = (b_loc - a_loc) / n_quad
                τ_total += h * (sum(geometric_torsion(seg, s) for s in ss) -
                                geometric_torsion(seg, a_loc)/2 - geometric_torsion(seg, b_loc)/2)
            end
        end
        s_seg_start = s_seg_end
    end

    # Material twist: reuse total_material_twist logic over the same window.
    Ω_total = total_material_twist(path; s_start = s_lo, s_end = s_hi, n_quad)

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

"""
    Sample

One evaluated point on a `Path`: arc-length coordinate `s` plus all Frenet frame quantities.

Fields: `s`, `position`, `tangent`, `normal`, `binormal`, `curvature`,
`geometric_torsion`, `material_twist`.
"""
struct Sample
    s                 :: Float64
    position          :: Vector{Float64}
    tangent           :: Vector{Float64}
    normal            :: Vector{Float64}
    binormal          :: Vector{Float64}
    curvature         :: Float64
    geometric_torsion :: Float64
    material_twist    :: Float64
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
    L_eff = arc_length(seg)
    a_eff = seg.a * seg.shrinkage
    return atan(L_eff / a_eff)
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
    path::Path,
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
    geom_angle  = _segment_total_angle(seg) * frac
    geom_budget = max(2, ceil(Int, fidelity * geom_angle / (2π) * 32))

    # Twist budget: exact integral of material twist rate over [a, b]
    twist_angle  = abs(total_material_twist(path; s_start = a, s_end = b, n_quad = 32))
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
function sample_path(path::Path, s1::Real, s2::Real; fidelity::Float64 = 1.0)
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

function sample(path::Path, s_values)
    return [frame(path, s) for s in s_values]
end

function sample_uniform(path::Path; n::Int = 256)
    ss = range(path.s_start, path.s_end; length = n)
    return sample(path, ss)
end

