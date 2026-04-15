"""
Local optical properties of an ideal step-index fiber cross section.

This file intentionally models only quantities that are meaningful for a single
transverse slice of fiber of infinitesimal length. It excludes any property
that depends on fiber length, path through space, accumulated phase, or
concatenation of segments.

The baseline object is a circular step-index fiber cross section described by
core/cladding materials and diameters. Perturbations such as core ellipticity,
bending, axial tension, and twist are handled by separate functions with
explicit arguments rather than stored on the type.

NOTE: All the dω derivatives computed in this file are computed by ChatGPT-5.4 and have
not been checked by a human.

Example
-------
fiber = FiberCrossSection(
    GermaniaSilicaGlass(0.036),
    GermaniaSilicaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like"
)

λ = 1550e-9
T = 297.15

v = normalized_frequency(fiber, λ, T)
β = propagation_constant(fiber, λ, T)
Aeff = effective_mode_area(fiber, λ, T)
Δβ_bend = total_bending_birefringence(fiber, λ, T; bend_radius_m = 0.03, axial_tension_N = 0.5)
"""

if !isdefined(Main, :GermaniaSilicaGlass)
    include("material-properties.jl")
end

const STEP_INDEX_GLASS = Union{GermaniaSilicaGlass, FluorinatedSilicaGlass}
const SPEED_OF_LIGHT_M_PER_S = 299_792_458.0
const LP11_CUTOFF_V = 2.405
const MARCUSE_V_MIN = 1.2
const MARCUSE_V_MAX = 2.4

struct FiberCrossSection{T<:Real}
    manufacturer::Union{Nothing, String}
    model_number::Union{Nothing, String}
    core_material::STEP_INDEX_GLASS
    cladding_material::STEP_INDEX_GLASS
    core_diameter_m::T
    cladding_diameter_m::T
end

struct BirefringenceResponse{T}
    Δβ::T
    dω::T
end

function FiberCrossSection(
    core_material::STEP_INDEX_GLASS,
    cladding_material::STEP_INDEX_GLASS,
    core_diameter_m::Real,
    cladding_diameter_m::Real;
    manufacturer::Union{Nothing, AbstractString} = nothing,
    model_number::Union{Nothing, AbstractString} = nothing
)
    core_diameter = validate_positive_length(core_diameter_m, "core diameter")
    cladding_diameter = validate_positive_length(cladding_diameter_m, "cladding diameter")
    core_diameter < cladding_diameter || throw(ArgumentError(
        "core diameter must be smaller than cladding diameter"
    ))

    T = promote_type(typeof(core_diameter), typeof(cladding_diameter))
    return FiberCrossSection{T}(
        isnothing(manufacturer) ? nothing : String(manufacturer),
        isnothing(model_number) ? nothing : String(model_number),
        core_material,
        cladding_material,
        core_diameter,
        cladding_diameter
    )
end

function validate_positive_length(value::Real, name::AbstractString)
    x = float(value)
    if !(isfinite(x) && x > zero(x))
        throw(ArgumentError("$(name) must be a finite positive value in meters"))
    end
    return x
end

function validate_nonnegative(value::Real, name::AbstractString)
    x = float(value)
    if !(isfinite(x) && x >= zero(x))
        throw(ArgumentError("$(name) must be a finite nonnegative value"))
    end
    return x
end

function validate_positive_finite(value::Real, name::AbstractString)
    x = float(value)
    if !(isfinite(x) && x > zero(x))
        throw(ArgumentError("$(name) must be a finite positive value"))
    end
    return x
end

function validate_bend_radius(bend_radius_m::Real)
    R = float(bend_radius_m)
    if isinf(R) && R > zero(R)
        return R
    end
    if !(isfinite(R) && R > zero(R))
        throw(ArgumentError(
            "bend_radius_m must be a finite positive value in meters or Inf"
        ))
    end
    return R
end

function validate_axis_ratio(axis_ratio::Real)
    ε = float(axis_ratio)
    if !(isfinite(ε) && ε > zero(ε))
        throw(ArgumentError("axis_ratio must be a finite positive value"))
    end
    return ε
end

core_radius(fiber::FiberCrossSection) = fiber.core_diameter_m / 2
cladding_radius(fiber::FiberCrossSection) = fiber.cladding_diameter_m / 2

function core_refractive_index(style::SpectralStyle, fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    return refractive_index(style, fiber.core_material, λ_m, T_K)
end

core_refractive_index(fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    core_refractive_index(ValueOnly(), fiber, λ_m, T_K)

function cladding_refractive_index(style::SpectralStyle, fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    return refractive_index(style, fiber.cladding_material, λ_m, T_K)
end

cladding_refractive_index(fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    cladding_refractive_index(ValueOnly(), fiber, λ_m, T_K)

function waveguide_factor(V)
    α = one(V) + sqrt(one(V) + one(V))
    t = (4 + V^4)^(one(V) / 4)
    return α * V / (one(V) + t)
end

function waveguide_factor_prime(V)
    α = one(V) + sqrt(one(V) + one(V))
    t = (4 + V^4)^(one(V) / 4)
    dt_dV = V^3 / (4 + V^4)^(3 * one(V) / 4)
    den = one(V) + t
    return α * (den - V * dt_dV) / den^2
end

function modal_prefactor(V)
    α = one(V) + sqrt(one(V) + one(V))
    den = one(V) + (4 + V^4)^(one(V) / 4)
    return one(V) - α^2 / den^2
end

function modal_prefactor_prime(V)
    α = one(V) + sqrt(one(V) + one(V))
    t = (4 + V^4)^(one(V) / 4)
    dt_dV = V^3 / (4 + V^4)^(3 * one(V) / 4)
    den = one(V) + t
    return 2 * α^2 * dt_dV / den^3
end

function mode_terms(::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    λ = validate_model_wavelength(λ_m)
    n_core = core_refractive_index(fiber, λ, T_K)
    n_clad = cladding_refractive_index(fiber, λ, T_K)
    n_core > n_clad || throw(ArgumentError(
        "guided-mode calculations require n_core > n_cladding; got " *
        "n_core=$(n_core), n_cladding=$(n_clad)"
    ))

    r_core = core_radius(fiber)
    r_clad = cladding_radius(fiber)
    k0 = 2π / λ
    dk0_dω = one(λ) / SPEED_OF_LIGHT_M_PER_S
    na = sqrt(n_core^2 - n_clad^2)
    V = r_core * k0 * na
    g = waveguide_factor(V)
    q = modal_prefactor(V)
    β = sqrt((n_core^2) * k0^2 - g^2 / r_core^2)
    z = zero(β)

    return (
        core_radius = r_core,
        cladding_radius = r_clad,
        k0 = k0,
        dk0_dω = dk0_dω,
        n_core = n_core,
        dn_core_dω = z,
        n_clad = n_clad,
        dn_clad_dω = z,
        na = na,
        dna_dω = z,
        V = V,
        dV_dω = z,
        waveguide_factor = g,
        dwaveguide_factor_dω = z,
        modal_prefactor = q,
        dmodal_prefactor_dω = z,
        β = β,
        dβ_dω = z
    )
end

function mode_terms(::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    λ = validate_model_wavelength(λ_m)
    n_core_resp = core_refractive_index(WithDerivative(), fiber, λ, T_K)
    n_clad_resp = cladding_refractive_index(WithDerivative(), fiber, λ, T_K)
    n_core = n_core_resp.value
    n_clad = n_clad_resp.value
    n_core > n_clad || throw(ArgumentError(
        "guided-mode calculations require n_core > n_cladding; got " *
        "n_core=$(n_core), n_cladding=$(n_clad)"
    ))

    r_core = core_radius(fiber)
    r_clad = cladding_radius(fiber)
    k0 = 2π / λ
    dk0_dω = one(λ) / SPEED_OF_LIGHT_M_PER_S
    na = sqrt(n_core^2 - n_clad^2)
    dna_dω = (n_core * n_core_resp.dω - n_clad * n_clad_resp.dω) / na
    V = r_core * k0 * na
    dV_dω = r_core * (dk0_dω * na + k0 * dna_dω)
    g = waveguide_factor(V)
    dg_dω = waveguide_factor_prime(V) * dV_dω
    q = modal_prefactor(V)
    dq_dω = modal_prefactor_prime(V) * dV_dω
    β = sqrt((n_core^2) * k0^2 - g^2 / r_core^2)
    dβ_inner_dω = 2 * n_core * n_core_resp.dω * k0^2 +
                  2 * n_core^2 * k0 * dk0_dω -
                  2 * g * dg_dω / r_core^2
    dβ_dω = dβ_inner_dω / (2 * β)

    return (
        core_radius = r_core,
        cladding_radius = r_clad,
        k0 = k0,
        dk0_dω = dk0_dω,
        n_core = n_core,
        dn_core_dω = n_core_resp.dω,
        n_clad = n_clad,
        dn_clad_dω = n_clad_resp.dω,
        na = na,
        dna_dω = dna_dω,
        V = V,
        dV_dω = dV_dω,
        waveguide_factor = g,
        dwaveguide_factor_dω = dg_dω,
        modal_prefactor = q,
        dmodal_prefactor_dω = dq_dω,
        β = β,
        dβ_dω = dβ_dω
    )
end

function guided_refractive_indices(style::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    terms = mode_terms(style, fiber, λ_m, T_K)
    return terms.n_core, terms.n_clad
end

function guided_refractive_indices(style::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    terms = mode_terms(style, fiber, λ_m, T_K)
    return (
        SpectralResponse(terms.n_core, terms.dn_core_dω),
        SpectralResponse(terms.n_clad, terms.dn_clad_dω)
    )
end

guided_refractive_indices(fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    guided_refractive_indices(ValueOnly(), fiber, λ_m, T_K)

function relative_index_difference(style::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    n_core, n_clad = guided_refractive_indices(style, fiber, λ_m, T_K)
    return (n_core - n_clad) / n_clad
end

relative_index_difference(fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    relative_index_difference(ValueOnly(), fiber, λ_m, T_K)

numerical_aperture(style::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    mode_terms(style, fiber, λ_m, T_K).na

function numerical_aperture(style::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    terms = mode_terms(style, fiber, λ_m, T_K)
    return SpectralResponse(terms.na, terms.dna_dω)
end

numerical_aperture(fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    numerical_aperture(ValueOnly(), fiber, λ_m, T_K)

normalized_frequency(style::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    mode_terms(style, fiber, λ_m, T_K).V

function normalized_frequency(style::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    terms = mode_terms(style, fiber, λ_m, T_K)
    return SpectralResponse(terms.V, terms.dV_dω)
end

normalized_frequency(fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    normalized_frequency(ValueOnly(), fiber, λ_m, T_K)

propagation_constant(style::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    mode_terms(style, fiber, λ_m, T_K).β

function propagation_constant(style::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    terms = mode_terms(style, fiber, λ_m, T_K)
    return SpectralResponse(terms.β, terms.dβ_dω)
end

propagation_constant(fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    propagation_constant(ValueOnly(), fiber, λ_m, T_K)

function effective_mode_index(style::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    terms = mode_terms(style, fiber, λ_m, T_K)
    return terms.β / terms.k0
end

function effective_mode_index(style::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    terms = mode_terms(style, fiber, λ_m, T_K)
    value = terms.β / terms.k0
    dω = (terms.dβ_dω * terms.k0 - terms.β * terms.dk0_dω) / terms.k0^2
    return SpectralResponse(value, dω)
end

effective_mode_index(fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    effective_mode_index(ValueOnly(), fiber, λ_m, T_K)

function effective_group_index(
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    dλ::Real = 0.1e-9
)
    λ = validate_model_wavelength(λ_m)
    step = validate_positive_finite(dλ, "dλ")
    n_center = effective_mode_index(fiber, λ, T_K)
    n_minus = effective_mode_index(fiber, λ - step, T_K)
    n_plus = effective_mode_index(fiber, λ + step, T_K)
    dndλ = (n_plus - n_minus) / (2 * step)
    return n_center - λ * dndλ
end

function effective_mode_area(fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    v = normalized_frequency(fiber, λ_m, T_K)
    if !(MARCUSE_V_MIN <= v <= MARCUSE_V_MAX)
        @warn "Marcuse effective-area approximation is calibrated for 1.2 <= V <= 2.4; got V=$(v)"
    end
    w_over_r = 0.65 + 1.619 / v^1.5 + 2.879 / v^6
    w = w_over_r * core_radius(fiber)
    return π * w^2
end

function chromatic_dispersion_parameter(
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    dλ::Real = 0.1e-9
)
    λ = validate_model_wavelength(λ_m)
    step = validate_positive_finite(dλ, "dλ")
    n_center = effective_mode_index(fiber, λ, T_K)
    n_minus = effective_mode_index(fiber, λ - step, T_K)
    n_plus = effective_mode_index(fiber, λ + step, T_K)
    return -λ / SPEED_OF_LIGHT_M_PER_S * (n_plus - 2 * n_center + n_minus) / step^2 * 1e6
end

function group_velocity_dispersion_parameter(
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    dλ::Real = 0.1e-9
)
    D_SI = chromatic_dispersion_parameter(fiber, λ_m, T_K; dλ = dλ) * 1e-6
    λ = validate_model_wavelength(λ_m)
    return -(λ^2 / (2π * SPEED_OF_LIGHT_M_PER_S)) * D_SI
end

core_nonlinear_refractive_index(fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    nonlinear_refractive_index(fiber.core_material, λ_m, T_K)

function nonlinear_coefficient(fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    n2 = core_nonlinear_refractive_index(fiber, λ_m, T_K)
    Aeff = effective_mode_area(fiber, λ_m, T_K)
    λ = validate_model_wavelength(λ_m)
    return (2π / λ) * n2 / Aeff
end

is_single_mode(fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    normalized_frequency(fiber, λ_m, T_K) < LP11_CUTOFF_V

function cutoff_wavelength(
    fiber::FiberCrossSection,
    T_K::Real;
    λ_min::Real = MIN_VALID_WAVELENGTH_M,
    λ_max::Real = MAX_VALID_WAVELENGTH_M,
    atol::Real = 1e-12,
    maxiter::Integer = 200
)
    λ_lo = validate_positive_length(λ_min, "λ_min")
    λ_hi = validate_positive_length(λ_max, "λ_max")
    λ_lo < λ_hi || throw(ArgumentError("λ_min must be smaller than λ_max"))
    tolerance = validate_positive_finite(atol, "atol")

    f_lo = normalized_frequency(fiber, λ_lo, T_K) - LP11_CUTOFF_V
    f_hi = normalized_frequency(fiber, λ_hi, T_K) - LP11_CUTOFF_V

    if f_lo == 0.0
        return λ_lo
    elseif f_hi == 0.0
        return λ_hi
    elseif signbit(f_lo) == signbit(f_hi)
        throw(ArgumentError(
            "cutoff wavelength is not bracketed in [$(λ_lo), $(λ_hi)] m"
        ))
    end

    a = λ_lo
    b = λ_hi
    fa = f_lo

    for _ in 1:maxiter
        mid = (a + b) / 2
        fm = normalized_frequency(fiber, mid, T_K) - LP11_CUTOFF_V
        if abs(fm) <= tolerance || (b - a) / 2 <= tolerance
            return mid
        elseif signbit(fm) == signbit(fa)
            a = mid
            fa = fm
        else
            b = mid
        end
    end

    return (a + b) / 2
end

function eccentricity_squared(axis_ratio::Real; signed::Bool = false)
    ε = validate_axis_ratio(axis_ratio)
    if ε >= one(ε)
        return one(ε) - inv(ε)^2
    end
    value = one(ε) - ε^2
    return signed ? -value : value
end

function core_noncircularity_dω(style::SpectralStyle, fiber::FiberCrossSection, λ_m::Real, T_K::Real; axis_ratio::Real)
    terms = mode_terms(style, fiber, λ_m, T_K)
    ε = validate_axis_ratio(axis_ratio)
    ε == one(ε) && return BirefringenceResponse(zero(terms.β), zero(terms.β))
    χ = one(terms.n_core) - terms.n_clad^2 / terms.n_core^2
    dχ_dω = -2 * terms.n_clad * terms.dn_clad_dω / terms.n_core^2 +
            2 * terms.n_clad^2 * terms.dn_core_dω / terms.n_core^3
    h = 4 * log(terms.V)^3 / (terms.V^3 * (one(terms.V) + log(terms.V)))
    h_prime = h / terms.V * (3 / log(terms.V) - 3 - inv(one(terms.V) + log(terms.V)))
    prefactor = eccentricity_squared(ε; signed = true) / terms.core_radius
    Δβ = prefactor * χ^(3 / 2) * h
    dω = prefactor * ((3 / 2) * sqrt(χ) * dχ_dω * h + χ^(3 / 2) * h_prime * terms.dV_dω)
    return BirefringenceResponse(Δβ, dω)
end

function asymmetric_thermal_stress_dω(
    style::SpectralStyle,
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    axis_ratio::Real
)
    terms = mode_terms(style, fiber, λ_m, T_K)
    ε = validate_axis_ratio(axis_ratio)
    ε == one(ε) && return BirefringenceResponse(zero(terms.β), zero(terms.β))
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    α_core = cte(fiber.core_material, T_K)
    α_clad = cte(fiber.cladding_material, T_K)
    T_soft = softening_temperature(fiber.core_material, T_K)
    ν = poisson_ratio(fiber.core_material, T_K)
    const_factor = 0.5 * (p11 - p12) * (α_clad - α_core) * abs(T_soft - T_K) / (1 - ν^2) * ((ε - 1) / (ε + 1))
    Δβ = terms.k0 * terms.modal_prefactor * terms.n_core^3 * const_factor
    dω = const_factor * (
        terms.dk0_dω * terms.modal_prefactor * terms.n_core^3 +
        terms.k0 * terms.dmodal_prefactor_dω * terms.n_core^3 +
        terms.k0 * terms.modal_prefactor * 3 * terms.n_core^2 * terms.dn_core_dω
    )
    return BirefringenceResponse(Δβ, dω)
end

function bending_dω(
    style::SpectralStyle,
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    bend_radius_m::Real
)
    terms = mode_terms(style, fiber, λ_m, T_K)
    R = validate_bend_radius(bend_radius_m)
    isinf(R) && return BirefringenceResponse(zero(terms.β), zero(terms.β))
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    ν = poisson_ratio(fiber.core_material, T_K)
    geom = 0.5 * (terms.cladding_radius^2 / R^2)
    const_factor = (p11 - p12) * (1 + ν) * geom / 2
    Δβ = terms.k0 * terms.n_core^3 * const_factor
    dω = const_factor * (terms.dk0_dω * terms.n_core^3 + terms.k0 * 3 * terms.n_core^2 * terms.dn_core_dω)
    return BirefringenceResponse(Δβ, dω)
end

function axial_tension_dω(
    style::SpectralStyle,
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    bend_radius_m::Real,
    axial_tension_N::Real
)
    terms = mode_terms(style, fiber, λ_m, T_K)
    R = validate_bend_radius(bend_radius_m)
    tf = validate_nonnegative(axial_tension_N, "axial_tension_N")
    (isinf(R) || tf == zero(tf)) && return BirefringenceResponse(zero(terms.β), zero(terms.β))

    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    ν = poisson_ratio(fiber.core_material, T_K)
    E = youngs_modulus(fiber.core_material, T_K)
    geom = ((2 - 3 * ν) / (1 - ν)) * (terms.cladding_radius / R) * (tf / (π * terms.cladding_radius^2 * E))
    const_factor = (p11 - p12) * (1 + ν) * geom / 2
    Δβ = terms.k0 * terms.n_core^3 * const_factor
    dω = const_factor * (terms.dk0_dω * terms.n_core^3 + terms.k0 * 3 * terms.n_core^2 * terms.dn_core_dω)
    return BirefringenceResponse(Δβ, dω)
end

function total_bending_dω(
    style::SpectralStyle,
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    bend_radius_m::Real,
    axial_tension_N::Real = 0.0
)
    bend = bending_dω(style, fiber, λ_m, T_K; bend_radius_m = bend_radius_m)
    tension = axial_tension_dω(
        style,
        fiber,
        λ_m,
        T_K;
        bend_radius_m = bend_radius_m,
        axial_tension_N = axial_tension_N
    )
    return BirefringenceResponse(bend.Δβ + tension.Δβ, bend.dω + tension.dω)
end

function twisting_dω(
    style::SpectralStyle,
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    twist_rate_rad_per_m::Real
)
    terms = mode_terms(style, fiber, λ_m, T_K)
    tr = float(twist_rate_rad_per_m)
    isfinite(tr) || throw(ArgumentError("twist_rate_rad_per_m must be finite"))
    tr == zero(tr) && return BirefringenceResponse(zero(terms.β), zero(terms.β))
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    coeff = (p11 - p12) / 2
    Δβ = (one(terms.n_core) + coeff * terms.n_core^2) * tr
    dω = 2 * coeff * terms.n_core * terms.dn_core_dω * tr
    return BirefringenceResponse(Δβ, dω)
end

function total_birefringence_dω(
    style::SpectralStyle,
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    axis_ratio::Real = 1.0,
    bend_radius_m::Real = Inf,
    axial_tension_N::Real = 0.0,
    twist_rate_rad_per_m::Real = 0.0
)
    b_cnc = core_noncircularity_dω(style, fiber, λ_m, T_K; axis_ratio = axis_ratio)
    b_ats = asymmetric_thermal_stress_dω(style, fiber, λ_m, T_K; axis_ratio = axis_ratio)
    b_bend = total_bending_dω(
        style,
        fiber,
        λ_m,
        T_K;
        bend_radius_m = bend_radius_m,
        axial_tension_N = axial_tension_N
    )
    b_twist = twisting_dω(style, fiber, λ_m, T_K; twist_rate_rad_per_m = twist_rate_rad_per_m)
    return BirefringenceResponse(
        b_cnc.Δβ + b_ats.Δβ + b_bend.Δβ + b_twist.Δβ,
        b_cnc.dω + b_ats.dω + b_bend.dω + b_twist.dω
    )
end

core_noncircularity_birefringence(::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real; axis_ratio::Real) =
    core_noncircularity_dω(ValueOnly(), fiber, λ_m, T_K; axis_ratio = axis_ratio).Δβ

core_noncircularity_birefringence(::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real; axis_ratio::Real) =
    core_noncircularity_dω(WithDerivative(), fiber, λ_m, T_K; axis_ratio = axis_ratio)

core_noncircularity_birefringence(fiber::FiberCrossSection, λ_m::Real, T_K::Real; axis_ratio::Real) =
    core_noncircularity_birefringence(ValueOnly(), fiber, λ_m, T_K; axis_ratio = axis_ratio)

asymmetric_thermal_stress_birefringence(::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real; axis_ratio::Real) =
    asymmetric_thermal_stress_dω(ValueOnly(), fiber, λ_m, T_K; axis_ratio = axis_ratio).Δβ

asymmetric_thermal_stress_birefringence(::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real; axis_ratio::Real) =
    asymmetric_thermal_stress_dω(WithDerivative(), fiber, λ_m, T_K; axis_ratio = axis_ratio)

asymmetric_thermal_stress_birefringence(fiber::FiberCrossSection, λ_m::Real, T_K::Real; axis_ratio::Real) =
    asymmetric_thermal_stress_birefringence(ValueOnly(), fiber, λ_m, T_K; axis_ratio = axis_ratio)

bending_birefringence(::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real; bend_radius_m::Real) =
    bending_dω(ValueOnly(), fiber, λ_m, T_K; bend_radius_m = bend_radius_m).Δβ

bending_birefringence(::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real; bend_radius_m::Real) =
    bending_dω(WithDerivative(), fiber, λ_m, T_K; bend_radius_m = bend_radius_m)

bending_birefringence(fiber::FiberCrossSection, λ_m::Real, T_K::Real; bend_radius_m::Real) =
    bending_birefringence(ValueOnly(), fiber, λ_m, T_K; bend_radius_m = bend_radius_m)

axial_tension_birefringence(::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real; bend_radius_m::Real, axial_tension_N::Real) =
    axial_tension_dω(ValueOnly(), fiber, λ_m, T_K; bend_radius_m = bend_radius_m, axial_tension_N = axial_tension_N).Δβ

axial_tension_birefringence(::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real; bend_radius_m::Real, axial_tension_N::Real) =
    axial_tension_dω(WithDerivative(), fiber, λ_m, T_K; bend_radius_m = bend_radius_m, axial_tension_N = axial_tension_N)

axial_tension_birefringence(fiber::FiberCrossSection, λ_m::Real, T_K::Real; bend_radius_m::Real, axial_tension_N::Real) =
    axial_tension_birefringence(ValueOnly(), fiber, λ_m, T_K; bend_radius_m = bend_radius_m, axial_tension_N = axial_tension_N)

total_bending_birefringence(::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real; bend_radius_m::Real, axial_tension_N::Real = 0.0) =
    total_bending_dω(ValueOnly(), fiber, λ_m, T_K; bend_radius_m = bend_radius_m, axial_tension_N = axial_tension_N).Δβ

total_bending_birefringence(::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real; bend_radius_m::Real, axial_tension_N::Real = 0.0) =
    total_bending_dω(WithDerivative(), fiber, λ_m, T_K; bend_radius_m = bend_radius_m, axial_tension_N = axial_tension_N)

total_bending_birefringence(fiber::FiberCrossSection, λ_m::Real, T_K::Real; bend_radius_m::Real, axial_tension_N::Real = 0.0) =
    total_bending_birefringence(ValueOnly(), fiber, λ_m, T_K; bend_radius_m = bend_radius_m, axial_tension_N = axial_tension_N)

twisting_birefringence(::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real; twist_rate_rad_per_m::Real) =
    twisting_dω(ValueOnly(), fiber, λ_m, T_K; twist_rate_rad_per_m = twist_rate_rad_per_m).Δβ

twisting_birefringence(::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real; twist_rate_rad_per_m::Real) =
    twisting_dω(WithDerivative(), fiber, λ_m, T_K; twist_rate_rad_per_m = twist_rate_rad_per_m)

twisting_birefringence(fiber::FiberCrossSection, λ_m::Real, T_K::Real; twist_rate_rad_per_m::Real) =
    twisting_birefringence(ValueOnly(), fiber, λ_m, T_K; twist_rate_rad_per_m = twist_rate_rad_per_m)

function total_birefringence(::ValueOnly, fiber::FiberCrossSection, λ_m::Real, T_K::Real; axis_ratio::Real = 1.0, bend_radius_m::Real = Inf, axial_tension_N::Real = 0.0, twist_rate_rad_per_m::Real = 0.0)
    return total_birefringence_dω(
        ValueOnly(),
        fiber,
        λ_m,
        T_K;
        axis_ratio = axis_ratio,
        bend_radius_m = bend_radius_m,
        axial_tension_N = axial_tension_N,
        twist_rate_rad_per_m = twist_rate_rad_per_m
    ).Δβ
end

function total_birefringence(::WithDerivative, fiber::FiberCrossSection, λ_m::Real, T_K::Real; axis_ratio::Real = 1.0, bend_radius_m::Real = Inf, axial_tension_N::Real = 0.0, twist_rate_rad_per_m::Real = 0.0)
    return total_birefringence_dω(
        WithDerivative(),
        fiber,
        λ_m,
        T_K;
        axis_ratio = axis_ratio,
        bend_radius_m = bend_radius_m,
        axial_tension_N = axial_tension_N,
        twist_rate_rad_per_m = twist_rate_rad_per_m
    )
end

function total_birefringence(fiber::FiberCrossSection, λ_m::Real, T_K::Real; axis_ratio::Real = 1.0, bend_radius_m::Real = Inf, axial_tension_N::Real = 0.0, twist_rate_rad_per_m::Real = 0.0)
    return total_birefringence(
        ValueOnly(),
        fiber,
        λ_m,
        T_K;
        axis_ratio = axis_ratio,
        bend_radius_m = bend_radius_m,
        axial_tension_N = axial_tension_N,
        twist_rate_rad_per_m = twist_rate_rad_per_m
    )
end
