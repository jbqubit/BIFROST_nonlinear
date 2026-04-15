"""
Local optical properties of an ideal step-index fiber cross section.

This file intentionally models only quantities that are meaningful for a single
transverse slice of fiber of infinitesimal length. It excludes any property that depends on fiber
length, path through space, accumulated phase, or concatenation of segments.

The baseline object is a circular step-index fiber cross section described by
core/cladding materials and diameters. Perturbations such as core ellipticity,
bending, axial tension, and twist are handled by separate functions with
explicit arguments rather than stored on the type.

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

struct FiberCrossSection
    manufacturer::Union{Nothing, String}
    model_number::Union{Nothing, String}
    core_material::STEP_INDEX_GLASS
    cladding_material::STEP_INDEX_GLASS
    core_diameter_m::Float64
    cladding_diameter_m::Float64
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

    return FiberCrossSection(
        isnothing(manufacturer) ? nothing : String(manufacturer),
        isnothing(model_number) ? nothing : String(model_number),
        core_material,
        cladding_material,
        core_diameter,
        cladding_diameter
    )
end

function validate_positive_length(value::Real, name::AbstractString)
    x = Float64(value)
    if !(isfinite(x) && x > 0.0)
        throw(ArgumentError("$(name) must be a finite positive value in meters"))
    end
    return x
end

function validate_nonnegative(value::Real, name::AbstractString)
    x = Float64(value)
    if !(isfinite(x) && x >= 0.0)
        throw(ArgumentError("$(name) must be a finite nonnegative value"))
    end
    return x
end

function validate_positive_finite(value::Real, name::AbstractString)
    x = Float64(value)
    if !(isfinite(x) && x > 0.0)
        throw(ArgumentError("$(name) must be a finite positive value"))
    end
    return x
end

function validate_bend_radius(bend_radius_m::Real)
    R = Float64(bend_radius_m)
    if isinf(R) && R > 0.0
        return R
    end
    if !(isfinite(R) && R > 0.0)
        throw(ArgumentError(
            "bend_radius_m must be a finite positive value in meters or Inf"
        ))
    end
    return R
end

function validate_axis_ratio(axis_ratio::Real)
    ε = Float64(axis_ratio)
    if !(isfinite(ε) && ε > 0.0)
        throw(ArgumentError("axis_ratio must be a finite positive value"))
    end
    return ε
end

core_radius(fiber::FiberCrossSection) = fiber.core_diameter_m / 2
cladding_radius(fiber::FiberCrossSection) = fiber.cladding_diameter_m / 2

function core_refractive_index(fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    return refractive_index(fiber.core_material, λ_m, T_K)
end

function cladding_refractive_index(fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    return refractive_index(fiber.cladding_material, λ_m, T_K)
end

function guided_refractive_indices(fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    n_core = core_refractive_index(fiber, λ_m, T_K)
    n_clad = cladding_refractive_index(fiber, λ_m, T_K)
    n_core > n_clad || throw(ArgumentError(
        "guided-mode calculations require n_core > n_cladding; got " *
        "n_core=$(n_core), n_cladding=$(n_clad)"
    ))
    return n_core, n_clad
end

function relative_index_difference(fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    n_core, n_clad = guided_refractive_indices(fiber, λ_m, T_K)
    return (n_core - n_clad) / n_clad
end

function numerical_aperture(fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    n_core, n_clad = guided_refractive_indices(fiber, λ_m, T_K)
    return sqrt(n_core^2 - n_clad^2)
end

function normalized_frequency(fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    n_core, n_clad = guided_refractive_indices(fiber, λ_m, T_K)
    return core_radius(fiber) * (2π / Float64(λ_m)) * sqrt(n_core^2 - n_clad^2)
end

function propagation_constant(fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    n_core, n_clad = guided_refractive_indices(fiber, λ_m, T_K)
    r_core = core_radius(fiber)
    v = normalized_frequency(fiber, λ_m, T_K)
    k0 = 2π / Float64(λ_m)
    waveguide_term = ((1 + sqrt(2)) * v) / (1 + (4 + v^4)^(1 / 4))
    return sqrt((n_core^2) * k0^2 - (waveguide_term^2) / r_core^2)
end

effective_mode_index(fiber::FiberCrossSection, λ_m::Real, T_K::Real) =
    propagation_constant(fiber, λ_m, T_K) / (2π / Float64(λ_m))

function effective_group_index(
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    dλ::Real = 0.1e-9
)
    λ = Float64(λ_m)
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
    λ = Float64(λ_m)
    step = validate_positive_finite(dλ, "dλ")
    n_center = effective_mode_index(fiber, λ, T_K)
    n_minus = effective_mode_index(fiber, λ - step, T_K)
    n_plus = effective_mode_index(fiber, λ + step, T_K)
    return -λ / SPEED_OF_LIGHT_M_PER_S * (n_plus - 2n_center + n_minus) / step^2 * 1e6
end

function group_velocity_dispersion_parameter(
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    dλ::Real = 0.1e-9
)
    D_SI = chromatic_dispersion_parameter(fiber, λ_m, T_K; dλ = dλ) * 1e-6
    return -(Float64(λ_m)^2 / (2π * SPEED_OF_LIGHT_M_PER_S)) * D_SI
end

function core_nonlinear_refractive_index(
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real
)
    return nonlinear_refractive_index(fiber.core_material, λ_m, T_K)
end

function nonlinear_coefficient(fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    n2 = core_nonlinear_refractive_index(fiber, λ_m, T_K)
    Aeff = effective_mode_area(fiber, λ_m, T_K)
    return (2π / Float64(λ_m)) * n2 / Aeff
end

function is_single_mode(fiber::FiberCrossSection, λ_m::Real, T_K::Real)
    return normalized_frequency(fiber, λ_m, T_K) < LP11_CUTOFF_V
end

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
    if ε >= 1.0
        return 1 - inv(ε)^2
    end
    value = 1 - ε^2
    return signed ? -value : value
end

function core_noncircularity_birefringence(
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    axis_ratio::Real
)
    validate_axis_ratio(axis_ratio) == 1.0 && return 0.0
    n_core, n_clad = guided_refractive_indices(fiber, λ_m, T_K)
    v = normalized_frequency(fiber, λ_m, T_K)
    e2 = eccentricity_squared(axis_ratio; signed = true)
    return (e2 * (1 - n_clad^2 / n_core^2)^(3 / 2)) / core_radius(fiber) *
           (4 / v^3) * log(v)^3 / (1 + log(v))
end

function asymmetric_thermal_stress_birefringence(
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    axis_ratio::Real
)
    ε = validate_axis_ratio(axis_ratio)
    ε == 1.0 && return 0.0
    n_core, _ = guided_refractive_indices(fiber, λ_m, T_K)
    β = propagation_constant(fiber, λ_m, T_K)
    v = normalized_frequency(fiber, λ_m, T_K)
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    α_core = cte(fiber.core_material, T_K)
    α_clad = cte(fiber.cladding_material, T_K)
    T_soft = softening_temperature(fiber.core_material, T_K)
    ν = poisson_ratio(fiber.core_material, T_K)
    k0 = 2π / Float64(λ_m)
    prefactor = k0 * (1 - core_radius(fiber)^2 * (n_core^2 * k0^2 - β^2) / v^2)
    stress_factor = 0.5 * n_core^3 * (p11 - p12) * (α_clad - α_core) *
                    abs(T_soft - Float64(T_K)) / (1 - ν^2)
    return prefactor * stress_factor * ((ε - 1) / (ε + 1))
end

function bending_birefringence(
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    bend_radius_m::Real
)
    R = validate_bend_radius(bend_radius_m)
    isinf(R) && return 0.0
    n_core, _ = guided_refractive_indices(fiber, λ_m, T_K)
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    ν = poisson_ratio(fiber.core_material, T_K)
    r = cladding_radius(fiber)
    k0 = 2π / Float64(λ_m)
    return k0 * (n_core^3 / 2) * (p11 - p12) * (1 + ν) * 0.5 * (r^2 / R^2)
end

function axial_tension_birefringence(
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    bend_radius_m::Real,
    axial_tension_N::Real
)
    R = validate_bend_radius(bend_radius_m)
    tf = validate_nonnegative(axial_tension_N, "axial_tension_N")
    (isinf(R) || tf == 0.0) && return 0.0

    n_core, _ = guided_refractive_indices(fiber, λ_m, T_K)
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    ν = poisson_ratio(fiber.core_material, T_K)
    E = youngs_modulus(fiber.core_material, T_K)
    r = cladding_radius(fiber)
    k0 = 2π / Float64(λ_m)
    tension_term = ((2 - 3ν) / (1 - ν)) * (r / R) * (tf / (π * r^2 * E))
    return k0 * (n_core^3 / 2) * (p11 - p12) * (1 + ν) * tension_term
end

function total_bending_birefringence(
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    bend_radius_m::Real,
    axial_tension_N::Real = 0.0
)
    return bending_birefringence(fiber, λ_m, T_K; bend_radius_m = bend_radius_m) +
           axial_tension_birefringence(
               fiber,
               λ_m,
               T_K;
               bend_radius_m = bend_radius_m,
               axial_tension_N = axial_tension_N
           )
end

function twisting_birefringence(
    fiber::FiberCrossSection,
    λ_m::Real,
    T_K::Real;
    twist_rate_rad_per_m::Real
)
    tr = Float64(twist_rate_rad_per_m)
    isfinite(tr) || throw(ArgumentError("twist_rate_rad_per_m must be finite"))
    tr == 0.0 && return 0.0
    n_core, _ = guided_refractive_indices(fiber, λ_m, T_K)
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    return (1 + (n_core^2 / 2) * (p11 - p12)) * tr
end


