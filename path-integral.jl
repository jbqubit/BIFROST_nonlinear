using LinearAlgebra


"""
Code Overview: 
- PiecewiseProfile and FiberInput define the fiber model.
- make_generator turns that model into the local Jones dynamics (generator).
- the numerical propagation stack consists of
    - exp_jones_generator
    - exp_midpoint_step
    - phase_insensitive_error
    - propagate_interval!
    - propagate_piecewise

path-plot.jl sits on top of that stack; it provides visual diagnostics 

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
# Fiber input
# ----------------------------

struct FiberInput{FR,FTH,FDTW,FB,FT}
    Rb::FR                 # s -> bend radius [m]
    theta_b::FTH          # s -> bend angle [rad]
    dtwist::FDTW          # s -> d(twist)/ds [rad/m]
    bend_strength::FB     # k2 -> Δβ_b [1/m]
    twist_strength::FT    # tau -> Δβ_t [1/m], tau in rad/m
end

include("path-plot.jl")

# ----------------------------
# Generator K(s)
# ----------------------------

function make_generator(f::FiberInput)
    return function (s::Real)
        R = f.Rb(s)

        if isinf(R)
            kx = 0.0
            ky = 0.0
            k2 = 0.0
        else
            invR = 1.0 / R
            c = cos(f.theta_b(s))
            sn = sin(f.theta_b(s))
            kx = invR * c
            ky = invR * sn
            k2 = kx*kx + ky*ky
        end

        tau = f.dtwist(s)   # rad/m

        Δβt = f.twist_strength(tau)

        if k2 == 0.0
            Δβb = 0.0
            c2φ = 1.0
            s2φ = 0.0
        else
            Δβb = f.bend_strength(k2)
            # double-angle quantities from curvature vector directly
            c2φ = (kx*kx - ky*ky) / k2
            s2φ = (2*kx*ky) / k2
        end

        return ComplexF64[
             0.5im*Δβb*c2φ             0.5im*Δβb*s2φ - 0.5*Δβt
             0.5im*Δβb*s2φ + 0.5*Δβt   -0.5im*Δβb*c2φ
        ]
    end
end

# ----------------------------
# One exponential-midpoint step
# ----------------------------

"""
    sinhc()

sinhc(μ) = sinh(μ)/μ, with Taylor expansion around μ=0 to avoid numerical issues
(catastrophic cancellation).
"""
function sinhc(μ::ComplexF64)
    if abs(μ) < 1e-8
        μ2 = μ * μ
        return 1 + μ2 / 6 + μ2 * μ2 / 120 + μ2 * μ2 * μ2 / 5040
    end
    return sinh(μ) / μ
end

"""
    exp_jones_generator(A::Matrix{ComplexF64})

Closed-form exponential for a `2×2` Jones generator.

For an exactly traceless matrix `A`, Cayley-Hamilton gives `A^2 = -det(A) I`, so

`exp(A) = cosh(μ) I + sinh(μ)/μ * A`, where `μ^2 = -det(A)`.

This implementation also removes any tiny numerical trace drift by factoring out
`exp(tr(A)/2)` and applying the closed form to the centered traceless part.
"""
function exp_jones_generator(A::Matrix{ComplexF64})
    @assert size(A) == (2, 2) "exp_jones_generator expects a 2x2 matrix"

    halftr = 0.5 * (A[1, 1] + A[2, 2])

    b11 = A[1, 1] - halftr
    b12 = A[1, 2]
    b21 = A[2, 1]
    b22 = A[2, 2] - halftr

    μ2 = -(b11 * b22 - b12 * b21)
    μ = sqrt(μ2)
    pref = exp(halftr)
    c = cosh(μ)
    s = sinhc(μ)

    return ComplexF64[
        pref * (c + s * b11)  pref * (s * b12)
        pref * (s * b21)      pref * (c + s * b22)
    ]
end

function exp_midpoint_step(K, s::Float64, h::Float64, J::Matrix{ComplexF64})
    M = K(s + 0.5h)
    return exp_jones_generator(h * M) * J
end

# ----------------------------
# Matrix error metric
# Global phase insensitive
# ----------------------------

function phase_insensitive_error(A::Matrix{ComplexF64}, B::Matrix{ComplexF64})
    # Remove best common phase from B relative to A
    α = tr(A' * B)
    if abs(α) == 0
        return opnorm(A - B)
    end
    ϕ = α / abs(α)
    return opnorm(A - ϕ * B)
end

# ----------------------------
# Adaptive step-doubling on one smooth interval
# ----------------------------

struct PropagatorStats
    accepted_steps::Int
    rejected_steps::Int
end

function propagate_interval!(
    K,
    s0::Float64,
    s1::Float64,
    J0::Matrix{ComplexF64};
    rtol::Float64 = 1e-9,
    atol::Float64 = 1e-12,
    h_init::Float64 = (s1 - s0) / 100,
    h_min::Float64 = 1e-12,
    h_max::Float64 = s1 - s0,
    safety::Float64 = 0.9,
    growth_max::Float64 = 2.0,
    shrink_min::Float64 = 0.2
)
    s = s0
    J = copy(J0)
    h = min(h_init, s1 - s0)

    accepted = 0
    rejected = 0

    while s < s1
        h = min(h, s1 - s)
        @assert h >= h_min "Step size underflow near s=$s"

        # one full step
        J_full = exp_midpoint_step(K, s, h, J)

        # two half steps
        J_half = exp_midpoint_step(K, s, 0.5h, J)
        J_twohalf = exp_midpoint_step(K, s + 0.5h, 0.5h, J_half)

        err_abs = phase_insensitive_error(J_full, J_twohalf)
        scale = max(opnorm(J_twohalf), opnorm(J_full))
        tol = atol + rtol * scale

        if err_abs <= tol
            # accept better solution
            s += h
            J = J_twohalf
            accepted += 1

            # midpoint method is 2nd order, step-doubling estimate ~ O(h^3)
            if err_abs == 0.0
                fac = growth_max
            else
                fac = safety * (tol / err_abs)^(1/3)
                fac = clamp(fac, 1.0, growth_max)
            end
            h = min(h * fac, h_max)
        else
            rejected += 1
            fac = safety * (tol / err_abs)^(1/3)
            fac = clamp(fac, shrink_min, 0.5)
            h *= fac
            @assert h >= h_min "Step size fell below h_min near s=$s"
        end
    end

    return J, PropagatorStats(accepted, rejected)
end

# ----------------------------
# Domain decomposition with optional jump matrices
# ----------------------------

"""
    propagate_piecewise(K, breaks; jumps=Dict(), kwargs...)

Propagate dJ/ds = K(s) J from breaks[1] to breaks[end].

Arguments
---------
- K: callable generator
- breaks: sorted vector [s0, s1, ..., sN]
- jumps: optional Dict{Float64, Matrix{ComplexF64}} giving lumped jump matrices
         applied immediately AFTER reaching that breakpoint.

Returns
-------
- J_final
- stats_per_interval
"""
function propagate_piecewise(
    K,
    breaks::Vector{Float64};
    jumps::Dict{Float64, Matrix{ComplexF64}} = Dict{Float64, Matrix{ComplexF64}}(),
    J0::Matrix{ComplexF64} = Matrix{ComplexF64}(I, 2, 2),
    kwargs...
)
    @assert issorted(breaks)
    J = copy(J0)
    stats = PropagatorStats[]

    for i in 1:length(breaks)-1
        sL = breaks[i]
        sR = breaks[i+1]

        J, st = propagate_interval!(K, sL, sR, J; kwargs...)
        push!(stats, st)

        if haskey(jumps, sR)
            J = jumps[sR] * J
        end
    end

    return J, stats
end
