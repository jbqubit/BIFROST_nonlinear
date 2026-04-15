"""
Fiber path, source specification, and segment-based authoring.

This file has two layers.

Authoring layer:
- `FiberSpec`, `BendSegment`, and `TwistSegment` provide a readable way to
  build a fiber segment by segment on absolute `s`-intervals.
- `FiberSpec` also holds a `FiberCrossSection` (step-index core/cladding model)
  so the authored path is tied to a specific transverse fiber type.
- `bend!`, `twist!`, and `build` compile that authored description into the
  lower-level runtime representation used by the solver in path-integral.jl.

Runtime/source layer:
- `AbstractBirefringenceSource` and concrete subtypes such as `BendSource` and
  `TwistSource` encode the local geometric quantities of curvature and twist
  rate.
- Each source answers three distinct questions:
  (1) What is the underlying local geometric quantity?
  (2) What matrix does that produce in `K(s)`?
  (3) What matrix does that produce in `Kω(s)`?

Supporting functionality includes:
- coverage validation and breakpoint handling
- `Fiber` assembly from low-level sources or from `FiberSpec`
- source-level and fiber-level `K(s)` and `Kω(s)` assembly
- bend/twist diagnostic helpers used by plotting
"""

if !isdefined(Main, :FiberCrossSection)
    include("fiber-cross-section.jl")
end

# ----------------------------
# Example Use
# ----------------------------

"""
xs = FiberCrossSection(
    GermaniaSilicaGlass(0.036),
    GermaniaSilicaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like"
)

spec = FiberSpec(0.0, 20.0; cross_section = xs)

twist!(spec, 0.0, 5.0; rate = X15)

bend!(spec, 0.0, 2.0;
    angle = π / 2,
    axis = 0.0
)

bend!(spec, 2.0, 3.0;
    angle = 0.0,
    axis = 0.0
)

twist!(spec, 5.0, 7.0; rate = X57)

bend!(spec, 3.0, 4.0;
    angle = π / 4,
    axis = π / 6
)

bend!(spec, 4.0, 7.0;
    angle = 0.0,
    axis = π / 6
)

fiber = build(spec)
"""


# ----------------------------
# Piecewise scalar profile
# ----------------------------

struct PiecewiseProfile{F}
    breaks::Vector{Float64}   # [s0, s1, ..., sN]
    pieces::Vector{F}         # length N-1  
    function PiecewiseProfile(breaks::Vector{Float64}, pieces::Vector{F}) where {F}
        @assert length(breaks) >= 2
        @assert length(pieces) == length(breaks) - 1
        @assert issorted(breaks)
        new{F}(breaks, pieces)
    end
end

function (p::PiecewiseProfile)(s::Real)
    # this syntax: take an object p of type PiecewiseProfile, allow it to be called like 
    # a function with one argument s::Real
    @assert p.breaks[1] <= s <= p.breaks[end] "s out of bounds"
    i = searchsortedlast(p.breaks, s)
    i = min(i, length(p.pieces))   # handles s == last breakpoint
    return p.pieces[i](s)
end

# ----------------------------
# Fiber sources and assembly
# ----------------------------

abstract type AbstractBirefringenceSource end
# this syntax: declares abstract type, no instances can be created of it

struct SourceResponseCallback{F}
    f::F
end

function coverage_segments_from_breaks(breaks::Vector{Float64})
    return [(breaks[i], breaks[i + 1]) for i in 1:length(breaks)-1]
end

function normalize_breakpoints(breakpoints::Vector{Float64})
    return sort(unique(copy(breakpoints)))
end

struct BendSource{FR,FTH,FRESP} <: AbstractBirefringenceSource
    Rb::FR
    theta_b::FTH
    response::FRESP
    coverage::Vector{Tuple{Float64, Float64}}
    breakpoints::Vector{Float64}
end

function BendSource(
    Rb,
    theta_b,
    response::SourceResponseCallback;
    coverage::Union{Nothing, Vector{Tuple{Float64, Float64}}} = nothing,
    breakpoints::Vector{Float64} = Float64[]
)
    cov = isnothing(coverage) ? coverage_segments_from_breaks(breakpoints) : copy(coverage)
    return BendSource{typeof(Rb), typeof(theta_b), typeof(response.f)}(
        Rb,
        theta_b,
        response.f,
        cov,
        normalize_breakpoints(breakpoints)
    )
end

pair_response_callback(strength_K::Function, strength_Kω::Function) =
    SourceResponseCallback((mode, x) -> (; Δβ = strength_K(x), dω = strength_Kω(x)))

function BendSource(
    Rb,
    theta_b,
    bend_strength::Function,
    bend_strength_omega::Function = k2 -> 0.0;
    coverage::Union{Nothing, Vector{Tuple{Float64, Float64}}} = nothing,
    breakpoints::Vector{Float64} = Float64[]
)
    return BendSource(
        Rb,
        theta_b,
        pair_response_callback(bend_strength, bend_strength_omega);
        coverage = coverage,
        breakpoints = breakpoints
    )
end

function bend_components(src::BendSource, s::Real)
    R = src.Rb(s)
    if isinf(R)
        return (kx = 0.0, ky = 0.0, k2 = 0.0)
    end

    invR = 1.0 / R
    c = cos(src.theta_b(s))
    sn = sin(src.theta_b(s))
    kx = invR * c
    ky = invR * sn
    return (kx = kx, ky = ky, k2 = kx * kx + ky * ky)
end

struct TwistSource{FDTW,FRESP} <: AbstractBirefringenceSource
    dtwist::FDTW
    response::FRESP
    coverage::Vector{Tuple{Float64, Float64}}
    breakpoints::Vector{Float64}
end

function TwistSource(
    dtwist,
    response::SourceResponseCallback;
    coverage::Union{Nothing, Vector{Tuple{Float64, Float64}}} = nothing,
    breakpoints::Vector{Float64} = Float64[]
)
    cov = isnothing(coverage) ? coverage_segments_from_breaks(breakpoints) : copy(coverage)
    return TwistSource{typeof(dtwist), typeof(response.f)}(
        dtwist,
        response.f,
        cov,
        normalize_breakpoints(breakpoints)
    )
end

function TwistSource(
    dtwist,
    twist_strength::Function,
    twist_strength_omega::Function = tau -> 0.0;
    coverage::Union{Nothing, Vector{Tuple{Float64, Float64}}} = nothing,
    breakpoints::Vector{Float64} = Float64[]
)
    return TwistSource(
        dtwist,
        pair_response_callback(twist_strength, twist_strength_omega);
        coverage = coverage,
        breakpoints = breakpoints
    )
end

source_breakpoints(src::AbstractBirefringenceSource) = src.breakpoints
source_coverage(src::AbstractBirefringenceSource) = src.coverage

function source_name(src::AbstractBirefringenceSource)
    return string(nameof(typeof(src)))
end

"""
    uncovered_intervals()

Given a set of coverage intervals and a target interval [s_start, s_end], 
return the sub-intervals of [s_start, s_end] that are not covered by any 
of the coverage intervals. Intervals are considered uncovered if they are 
not covered by any coverage interval within an absolute tolerance `atol` to 
account for floating-point precision issues.

We hope that return value is empty, which would indicate that the coverage 
intervals fully cover the target interval.  
"""
function uncovered_intervals(
    coverage::Vector{Tuple{Float64, Float64}},
    s_start::Float64,
    s_end::Float64;
    atol::Float64 = 1e-12
)
    @assert s_end > s_start "Require s_end > s_start"
    if isempty(coverage)
        return [(s_start, s_end)]
    end

    segs = sort(copy(coverage), by = first)
    uncovered = Tuple{Float64, Float64}[]
    cursor = s_start

    for (a_raw, b_raw) in segs
        a = max(a_raw, s_start)
        b = min(b_raw, s_end)
        if b <= a + atol
            continue
        end
        if a > cursor + atol
            push!(uncovered, (cursor, a))
        end
        cursor = max(cursor, b)
        if cursor >= s_end - atol
            cursor = s_end
            break
        end
    end

    if cursor < s_end - atol
        push!(uncovered, (cursor, s_end))
    end

    return uncovered
end

function format_intervals(intervals::Vector{Tuple{Float64, Float64}})
    return join(["[$(a), $(b)]" for (a, b) in intervals], ", ")
end

function validate_source_coverage(src::AbstractBirefringenceSource, s_start::Float64, s_end::Float64)
    gaps = uncovered_intervals(source_coverage(src), s_start, s_end)
    isempty(gaps) || error("Source $(source_name(src)) does not cover fiber domain [$s_start, $s_end]; uncovered intervals: $(format_intervals(gaps))")
    return nothing
end

struct Fiber
    s_start::Float64
    s_end::Float64
    sources::Vector{AbstractBirefringenceSource}
    function Fiber(s_start::Real, s_end::Real, sources::Vector{<:AbstractBirefringenceSource})
        s0 = Float64(s_start)
        s1 = Float64(s_end)
        @assert s1 > s0 "Fiber requires s_end > s_start"
        srcs = AbstractBirefringenceSource[s for s in sources]
        @assert !isempty(srcs) "Fiber requires at least one birefringence source"
        for src in srcs
            validate_source_coverage(src, s0, s1)
        end
        new(s0, s1, srcs)
    end
end

# ----------------------------
# Fiber authoring API
# ----------------------------

default_bend_strength_K(k2) = k2
default_bend_strength_Kω(k2) = 0.0
default_twist_strength_K(tau) = tau
default_twist_strength_Kω(tau) = 0.0

struct BendSegment
    s_start::Float64
    s_end::Float64
    angle::Float64
    axis::Float64
    strength_K::Function
    strength_Kω::Function
end

function BendSegment(
    s_start::Real,
    s_end::Real;
    angle::Real,
    axis::Real,
    strength_K::Function = default_bend_strength_K,
    strength_Kω::Function = default_bend_strength_Kω
)
    return BendSegment(
        Float64(s_start),
        Float64(s_end),
        Float64(angle),
        Float64(axis),
        strength_K,
        strength_Kω
    )
end

struct TwistSegment
    s_start::Float64
    s_end::Float64
    rate::Float64
    strength_K::Function
    strength_Kω::Function
end

function TwistSegment(
    s_start::Real,
    s_end::Real;
    rate::Real,
    strength_K::Function = default_twist_strength_K,
    strength_Kω::Function = default_twist_strength_Kω
)
    return TwistSegment(
        Float64(s_start),
        Float64(s_end),
        Float64(rate),
        strength_K,
        strength_Kω
    )
end

mutable struct FiberSpec
    cross_section::FiberCrossSection
    s_start::Float64
    s_end::Float64
    bend_segments::Vector{BendSegment}
    twist_segments::Vector{TwistSegment}
    default_bend_strength_K::Function
    default_bend_strength_Kω::Function
    default_twist_strength_K::Function
    default_twist_strength_Kω::Function
end

function FiberSpec(
    s_start::Real,
    s_end::Real;
    cross_section::FiberCrossSection,
    default_bend_strength_K::Function = default_bend_strength_K,
    default_bend_strength_Kω::Function = default_bend_strength_Kω,
    default_twist_strength_K::Function = default_twist_strength_K,
    default_twist_strength_Kω::Function = default_twist_strength_Kω
)
    return FiberSpec(
        cross_section,
        Float64(s_start),
        Float64(s_end),
        BendSegment[],
        TwistSegment[],
        default_bend_strength_K,
        default_bend_strength_Kω,
        default_twist_strength_K,
        default_twist_strength_Kω
    )
end

function bend!(
    spec::FiberSpec,
    s_start::Real,
    s_end::Real;
    angle::Real,
    axis::Real,
    strength_K::Function = spec.default_bend_strength_K,
    strength_Kω::Function = spec.default_bend_strength_Kω
)
    push!(spec.bend_segments, BendSegment(s_start, s_end; angle = angle, axis = axis, strength_K = strength_K, strength_Kω = strength_Kω))
    return spec
end

function twist!(
    spec::FiberSpec,
    s_start::Real,
    s_end::Real;
    rate::Real,
    strength_K::Function = spec.default_twist_strength_K,
    strength_Kω::Function = spec.default_twist_strength_Kω
)
    push!(spec.twist_segments, TwistSegment(s_start, s_end; rate = rate, strength_K = strength_K, strength_Kω = strength_Kω))
    return spec
end

function sorted_bend_segments(spec::FiberSpec)
    return sort(copy(spec.bend_segments), by = seg -> (seg.s_start, seg.s_end))
end

function sorted_twist_segments(spec::FiberSpec)
    return sort(copy(spec.twist_segments), by = seg -> (seg.s_start, seg.s_end))
end

function validate_segment_layout(segments, s_start::Float64, s_end::Float64, label::AbstractString)
    prev_end = s_start
    first_segment = true
    for seg in segments
        @assert s_start <= seg.s_start < seg.s_end <= s_end "$label segment [$((seg.s_start)), $((seg.s_end))] lies outside fiber domain [$s_start, $s_end]"
        if !first_segment
            @assert prev_end <= seg.s_start "$label segments may not overlap"
        end
        prev_end = seg.s_end
        first_segment = false
    end
    return nothing
end

function segment_profile_breaks(domain_start::Float64, domain_end::Float64, seg_start::Float64, seg_end::Float64)
    breaks = Float64[domain_start]
    if seg_start > domain_start
        push!(breaks, seg_start)
    end
    if seg_end < domain_end
        push!(breaks, seg_end)
    end
    push!(breaks, domain_end)
    return breaks
end

function compile_bend_source(spec::FiberSpec, seg::BendSegment)
    breaks = segment_profile_breaks(spec.s_start, spec.s_end, seg.s_start, seg.s_end)
    angle = seg.angle
    axis = angle >= 0.0 ? seg.axis : seg.axis + π
    active_radius = iszero(angle) ? Inf : (seg.s_end - seg.s_start) / abs(angle)
    pieces_Rb = Function[]
    pieces_theta = Function[]

    for i in 1:length(breaks)-1
        mid = 0.5 * (breaks[i] + breaks[i + 1])
        if seg.s_start <= mid <= seg.s_end
            push!(pieces_Rb, _ -> active_radius)
            push!(pieces_theta, _ -> axis)
        else
            push!(pieces_Rb, _ -> Inf)
            push!(pieces_theta, _ -> axis)
        end
    end

    return BendSource(
        PiecewiseProfile(copy(breaks), copy(pieces_Rb)),
        PiecewiseProfile(copy(breaks), copy(pieces_theta)),
        seg.strength_K,
        seg.strength_Kω;
        breakpoints = copy(breaks)
    )
end

function compile_twist_source(spec::FiberSpec, seg::TwistSegment)
    breaks = segment_profile_breaks(spec.s_start, spec.s_end, seg.s_start, seg.s_end)
    pieces_dtwist = Function[]

    for i in 1:length(breaks)-1
        mid = 0.5 * (breaks[i] + breaks[i + 1])
        if seg.s_start <= mid <= seg.s_end
            push!(pieces_dtwist, _ -> seg.rate)
        else
            push!(pieces_dtwist, _ -> 0.0)
        end
    end

    return TwistSource(
        PiecewiseProfile(copy(breaks), copy(pieces_dtwist)),
        seg.strength_K,
        seg.strength_Kω;
        breakpoints = copy(breaks)
    )
end

function zero_bend_source(spec::FiberSpec)
    return BendSource(
        _ -> Inf,
        _ -> 0.0,
        spec.default_bend_strength_K,
        spec.default_bend_strength_Kω;
        coverage = [(spec.s_start, spec.s_end)],
        breakpoints = [spec.s_start, spec.s_end]
    )
end

function build(spec::FiberSpec)
    bends = sorted_bend_segments(spec)
    twists = sorted_twist_segments(spec)
    validate_segment_layout(bends, spec.s_start, spec.s_end, "Bend")
    validate_segment_layout(twists, spec.s_start, spec.s_end, "Twist")

    sources = AbstractBirefringenceSource[]
    isempty(bends) && isempty(twists) && push!(sources, zero_bend_source(spec))

    for seg in bends
        push!(sources, compile_bend_source(spec, seg))
    end
    for seg in twists
        push!(sources, compile_twist_source(spec, seg))
    end

    return Fiber(spec.s_start, spec.s_end, sources)
end

# ----------------------------
# Generator K(s) and Curvature Kω(s)
# ----------------------------

zero_generator() = zeros(ComplexF64, 2, 2)

function generator_K_contribution(src::BendSource, s::Real)
    curv = bend_components(src, s)
    if curv.k2 == 0.0
        return zero_generator()
    end

    Δβb = src.response(:value, curv.k2).Δβ
    c2φ = (curv.kx * curv.kx - curv.ky * curv.ky) / curv.k2
    s2φ = (2 * curv.kx * curv.ky) / curv.k2
    return ComplexF64[
         0.5im * Δβb * c2φ             0.5im * Δβb * s2φ
         0.5im * Δβb * s2φ            -0.5im * Δβb * c2φ
    ]
end

function generator_Kω_contribution(src::BendSource, s::Real)
    curv = bend_components(src, s)
    if curv.k2 == 0.0
        return zero_generator()
    end

    Δβbω = src.response(:spectral, curv.k2).dω
    c2φ = (curv.kx * curv.kx - curv.ky * curv.ky) / curv.k2
    s2φ = (2 * curv.kx * curv.ky) / curv.k2
    return ComplexF64[
         0.5im * Δβbω * c2φ             0.5im * Δβbω * s2φ
         0.5im * Δβbω * s2φ            -0.5im * Δβbω * c2φ
    ]
end

function generator_K_contribution(src::TwistSource, s::Real)
    tau = src.dtwist(s)
    Δβt = src.response(:value, tau).Δβ
    return ComplexF64[
         0.0 + 0.0im   -0.5 * Δβt
         0.5 * Δβt      0.0 + 0.0im
    ]
end

function generator_Kω_contribution(src::TwistSource, s::Real)
    tau = src.dtwist(s)
    Δβtω = src.response(:spectral, tau).dω
    return ComplexF64[
         0.0 + 0.0im   -0.5 * Δβtω
         0.5 * Δβtω     0.0 + 0.0im
    ]
end

# ----------------------------
# Fiber construction from segments
# ----------------------------

function fiber_breakpoints(f::Fiber)
    points = Float64[f.s_start, f.s_end]
    for src in f.sources
        append!(points, [p for p in source_breakpoints(src) if f.s_start <= p <= f.s_end])
    end
    return sort(unique(points))
end

function generator_K(f::Fiber)
    return function (s::Real)
        K = zero_generator()
        for src in f.sources
            K .+= generator_K_contribution(src, s)
        end
        return K
    end
end

function generator_Kω(f::Fiber)
    return function (s::Real)
        Kω = zero_generator()
        for src in f.sources
            Kω .+= generator_Kω_contribution(src, s)
        end
        return Kω
    end
end

# ----------------------------
# Fiber diagnostics for plotting
# ----------------------------

function bend_geometry(f::Fiber, s::Real)
    kx = 0.0
    ky = 0.0
    for src in f.sources
        if src isa BendSource
            curv = bend_components(src, s)
            kx += curv.kx
            ky += curv.ky
        end
    end

    k2 = kx * kx + ky * ky
    if k2 == 0.0
        return (Rb = Inf, theta_b = 0.0, kx = 0.0, ky = 0.0, k2 = 0.0)
    end

    return (Rb = inv(sqrt(k2)), theta_b = atan(ky, kx), kx = kx, ky = ky, k2 = k2)
end

function twist_rate(f::Fiber, s::Real)
    τ = 0.0
    for src in f.sources
        if src isa TwistSource
            τ += src.dtwist(s)
        end
    end
    return τ
end
