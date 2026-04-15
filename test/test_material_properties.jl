using Test

if !isdefined(Main, :refractive_index)
    include("../material-properties.jl")
end

const SILICA_B_COEFFS = (
    (1.10127, -4.94251e-5, 5.27414e-7, -1.59700e-9, 1.75949e-12),
    (1.78752e-5, 4.76391e-5, -4.49019e-7, 1.44546e-9, -1.57223e-12),
    (7.93552e-1, -1.27815e-3, 1.84595e-5, -9.20275e-8, 1.48829e-10)
)

const SILICA_C_COEFFS = (
    (-8.906e-2, 9.0873e-6, -6.53638e-8, 7.77072e-11, 6.84605e-14),
    (2.97562e-1, -8.59578e-4, 6.59069e-6, -1.09482e-8, 7.85145e-13),
    (9.34454, -70.9788e-3, 1.01968e-4, -5.07660e-7, 8.21348e-10)
)

function reference_polynomial(coeffs, T_kelvin)
    return sum(coeffs[i] * T_kelvin^(i - 1) for i in eachindex(coeffs))
end

function reference_sellmeier_coefficients(T_celsius)
    T_kelvin = T_celsius + 273.15
    return ntuple(i -> (
        reference_polynomial(SILICA_B_COEFFS[i], T_kelvin),
        reference_polynomial(SILICA_C_COEFFS[i], T_kelvin)
    ), 3)
end

function reference_refractive_index(λ_meters, T_celsius)
    λ_um = λ_meters * 1e6
    total = 1.0
    for (B, C) in reference_sellmeier_coefficients(T_celsius)
        total += B * λ_um^2 / (λ_um^2 - C^2)
    end
    return sqrt(total)
end

@testset "Material property helpers" begin
    poly = TemperaturePolynomial((1.0, 2.0, 3.0, 4.0, 5.0))
    @test poly(2.0) ≈ 129.0
    @test temperature_kelvin(25.0) ≈ 298.15
    @test wavelength_microns(1.55e-6) ≈ 1.55
    @test_throws ArgumentError wavelength_microns(0.0)
    @test_throws ArgumentError wavelength_microns(-1.0)
end

@testset "SiO2 Sellmeier coefficients" begin
    silica = SiO2()
    terms = sellmeier_terms(silica)
    @test length(terms) == 3
    @test terms isa NTuple{3, SellmeierTerm}

    for T_celsius in (-40.0, 24.0, 120.0)
        actual = sellmeier_coefficients(silica, T_celsius)
        expected = reference_sellmeier_coefficients(T_celsius)
        @test length(actual) == 3
        for i in 1:3
            @test actual[i][1] ≈ expected[i][1] atol = 1e-12 rtol = 1e-12
            @test actual[i][2] ≈ expected[i][2] atol = 1e-12 rtol = 1e-12
        end
    end
end

@testset "SiO2 refractive index" begin
    silica = SiO2()

    cases = [
        (1.064e-6, 24.0),
        (1.31e-6, 0.0),
        (1.55e-6, 24.0),
        (1.625e-6, 80.0)
    ]

    for (λ_meters, T_celsius) in cases
        @test refractive_index(silica, λ_meters, T_celsius) ≈ reference_refractive_index(λ_meters, T_celsius) atol = 1e-12 rtol = 1e-12
    end
end
