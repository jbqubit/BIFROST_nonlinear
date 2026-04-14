"""
Fiber path and source specification.

This file defines:
- PiecewiseProfile
- birefringence source types and coverage validation
- Fiber assembly
- source-level and fiber-level K(s) / Kω(s) assembly
- breakpoint computation
- bend/twist diagnostic helpers used by plotting
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
    @assert p.breaks[1] <= s <= p.breaks[end] "s out of bounds"
    i = searchsortedlast(p.breaks, s)
    i = min(i, length(p.pieces))   # handles s == last breakpoint
    return p.pieces[i](s)
end

# ----------------------------
# Fiber sources and assembly
# ----------------------------

abstract type AbstractBirefringenceSource end

function coverage_segments_from_breaks(breaks::Vector{Float64})
    return [(breaks[i], breaks[i + 1]) for i in 1:length(breaks)-1]
end

function normalize_breakpoints(breakpoints::Vector{Float64})
    return sort(unique(copy(breakpoints)))
end

struct BendSource{FR,FTH,FB,FBW} <: AbstractBirefringenceSource
    Rb::FR
    theta_b::FTH
    bend_strength::FB
    bend_strength_omega::FBW
    coverage::Vector{Tuple{Float64, Float64}}
    breakpoints::Vector{Float64}
end

function BendSource(
    Rb,
    theta_b,
    bend_strength,
    bend_strength_omega = k2 -> 0.0;
    coverage::Union{Nothing, Vector{Tuple{Float64, Float64}}} = nothing,
    breakpoints::Vector{Float64} = Float64[]
)
    cov = isnothing(coverage) ? coverage_segments_from_breaks(breakpoints) : copy(coverage)
    return BendSource{typeof(Rb), typeof(theta_b), typeof(bend_strength), typeof(bend_strength_omega)}(
        Rb,
        theta_b,
        bend_strength,
        bend_strength_omega,
        cov,
        normalize_breakpoints(breakpoints)
    )
end

struct TwistSource{FDTW,FT,FTW} <: AbstractBirefringenceSource
    dtwist::FDTW
    twist_strength::FT
    twist_strength_omega::FTW
    coverage::Vector{Tuple{Float64, Float64}}
    breakpoints::Vector{Float64}
end

function TwistSource(
    dtwist,
    twist_strength,
    twist_strength_omega = tau -> 0.0;
    coverage::Union{Nothing, Vector{Tuple{Float64, Float64}}} = nothing,
    breakpoints::Vector{Float64} = Float64[]
)
    cov = isnothing(coverage) ? coverage_segments_from_breaks(breakpoints) : copy(coverage)
    return TwistSource{typeof(dtwist), typeof(twist_strength), typeof(twist_strength_omega)}(
        dtwist,
        twist_strength,
        twist_strength_omega,
        cov,
        normalize_breakpoints(breakpoints)
    )
end

source_breakpoints(src::AbstractBirefringenceSource) = src.breakpoints
source_coverage(src::AbstractBirefringenceSource) = src.coverage

function source_name(src::AbstractBirefringenceSource)
    return string(nameof(typeof(src)))
end

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
# Generator K(s) and Kω(s)
# ----------------------------

zero_generator() = zeros(ComplexF64, 2, 2)

function curvature_components(src::BendSource, s::Real)
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

function generator_contribution(src::BendSource, s::Real)
    curv = curvature_components(src, s)
    if curv.k2 == 0.0
        return zero_generator()
    end

    Δβb = src.bend_strength(curv.k2)
    c2φ = (curv.kx * curv.kx - curv.ky * curv.ky) / curv.k2
    s2φ = (2 * curv.kx * curv.ky) / curv.k2
    return ComplexF64[
         0.5im * Δβb * c2φ             0.5im * Δβb * s2φ
         0.5im * Δβb * s2φ            -0.5im * Δβb * c2φ
    ]
end

function generator_omega_contribution(src::BendSource, s::Real)
    curv = curvature_components(src, s)
    if curv.k2 == 0.0
        return zero_generator()
    end

    Δβbω = src.bend_strength_omega(curv.k2)
    c2φ = (curv.kx * curv.kx - curv.ky * curv.ky) / curv.k2
    s2φ = (2 * curv.kx * curv.ky) / curv.k2
    return ComplexF64[
         0.5im * Δβbω * c2φ             0.5im * Δβbω * s2φ
         0.5im * Δβbω * s2φ            -0.5im * Δβbω * c2φ
    ]
end

function generator_contribution(src::TwistSource, s::Real)
    tau = src.dtwist(s)
    Δβt = src.twist_strength(tau)
    return ComplexF64[
         0.0 + 0.0im   -0.5 * Δβt
         0.5 * Δβt      0.0 + 0.0im
    ]
end

function generator_omega_contribution(src::TwistSource, s::Real)
    tau = src.dtwist(s)
    Δβtω = src.twist_strength_omega(tau)
    return ComplexF64[
         0.0 + 0.0im   -0.5 * Δβtω
         0.5 * Δβtω     0.0 + 0.0im
    ]
end

function fiber_breakpoints(f::Fiber)
    points = Float64[f.s_start, f.s_end]
    for src in f.sources
        append!(points, [p for p in source_breakpoints(src) if f.s_start <= p <= f.s_end])
    end
    return sort(unique(points))
end

function make_generator(f::Fiber)
    return function (s::Real)
        K = zero_generator()
        for src in f.sources
            K .+= generator_contribution(src, s)
        end
        return K
    end
end

function make_generator_omega(f::Fiber)
    return function (s::Real)
        Kω = zero_generator()
        for src in f.sources
            Kω .+= generator_omega_contribution(src, s)
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
            curv = curvature_components(src, s)
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
