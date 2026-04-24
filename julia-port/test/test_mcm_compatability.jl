using Test
using MonteCarloMeasurements

if !isdefined(Main, :refractive_index)
    include("../material-properties.jl")
end

# Each section below exercises MCM (Particles) compatibility for one source file.
# Sections are grouped by the file under test and delineated by a banner. Add new
# sections in the same format as coverage expands to fiber-cross-section.jl,
# path-geometry.jl, etc.

# ═══════════════════════════════════════════════════════════════════════════════
# ▓▓▓  material-properties.jl  ▓▓▓
# ═══════════════════════════════════════════════════════════════════════════════

@testset "MCM :: material-properties.jl" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        λ = 1550e-9
        T_nom = 293.0
        T = T_nom ± 2.0

        # T-GUARDRAIL: Particles flow through each material's refractive_index
        for mat in (PURE_SILICA, PURE_GERMANIA,
                    GermaniaSilicaGlass(0.036), FluorinatedSilicaGlass(0.01))
            n = refractive_index(mat, λ, T)
            @test n isa Particles
            # pmean differs from nominal at O(δ²·∂²n/∂T²); loose tolerance.
            @test pmean(n) ≈ refractive_index(mat, λ, T_nom) rtol=1e-4
        end

        # T-GUARDRAIL: WithDerivative SpectralResponse lifts both fields
        for mat in (PURE_SILICA, PURE_GERMANIA,
                    GermaniaSilicaGlass(0.036), FluorinatedSilicaGlass(0.01))
            resp = refractive_index(WithDerivative(), mat, λ, T)
            @test resp.value isa Particles
            @test resp.dω isa Particles
            ref = refractive_index(WithDerivative(), mat, λ, T_nom)
            @test pmean(resp.value) ≈ ref.value rtol=1e-4
            @test pmean(resp.dω) ≈ ref.dω atol=1e-20 rtol=1e-3
        end

        # T-PHYSICS: mean of n(T_nom ± δ) matches n(T_nom) to O(δ²)
        δ = 2.0
        T_small = T_nom ± δ
        n_mean = pmean(refractive_index(PURE_SILICA, λ, T_small))
        n_ref = refractive_index(PURE_SILICA, λ, T_nom)
        @test abs(n_mean - n_ref) < 1e-3  # second-order curvature in δ=2 K

        # T-GUARDRAIL: out-of-range T (mean outside validity window) errors
        @test_throws ArgumentError refractive_index(PURE_SILICA, λ, 500.0 ± 1.0)

        # T-GUARDRAIL: T-independent material constants stay Float64 even under Particles T
        @test cte(PURE_SILICA, T) isa Float64
        @test softening_temperature(PURE_SILICA, T) isa Float64
        @test poisson_ratio(PURE_SILICA, T) isa Float64
        @test youngs_modulus(PURE_SILICA, T) isa Float64
        @test photoelastic_constants(PURE_SILICA, T) isa Tuple{Float64, Float64}
        @test nonlinear_refractive_index(PURE_SILICA, λ, T) isa Float64
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# ▓▓▓  fiber-cross-section.jl  ▓▓▓
# ═══════════════════════════════════════════════════════════════════════════════

if !isdefined(Main, :FiberCrossSection)
    include("../fiber-cross-section.jl")
end

@testset "MCM :: fiber-cross-section.jl" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        λ = 1550e-9
        T_nom = 293.0
        T = T_nom ± 2.0
        fiber = FiberCrossSection(
            GermaniaSilicaGlass(0.036),
            GermaniaSilicaGlass(0.0),
            8.2e-6,
            125e-6,
        )

        # ---- T-dependent mode quantities lift through Particles T
        β = propagation_constant(fiber, λ, T)
        @test β isa Particles
        @test pmean(β) ≈ propagation_constant(fiber, λ, T_nom) rtol=1e-4

        V = normalized_frequency(fiber, λ, T)
        @test V isa Particles
        @test pmean(V) ≈ normalized_frequency(fiber, λ, T_nom) rtol=1e-4

        @test effective_mode_area(fiber, λ, T) isa Particles
        @test nonlinear_coefficient(fiber, λ, T) isa Particles

        # WithDerivative carries Particles in both fields
        β_resp = propagation_constant(WithDerivative(), fiber, λ, T)
        @test β_resp.value isa Particles
        @test β_resp.dω isa Particles

        # ---- bending_birefringence with uncertain bend_radius_m
        R = 0.03 ± 0.003
        Δβ_bend = bending_birefringence(fiber, λ, T_nom; bend_radius_m = R)
        @test Δβ_bend isa Particles
        @test pmean(Δβ_bend) ≈
              bending_birefringence(fiber, λ, T_nom; bend_radius_m = 0.03) rtol=5e-2

        # ---- axial_tension_birefringence with uncertain axial_tension_N
        Δβ_tension = axial_tension_birefringence(fiber, λ, T_nom;
                                                  bend_radius_m = 0.03,
                                                  axial_tension_N = 0.5 ± 0.05)
        @test Δβ_tension isa Particles

        # ---- twisting_birefringence with uncertain twist_rate
        Δβ_twist = twisting_birefringence(fiber, λ, T_nom;
                                           twist_rate_rad_per_m = 10.0 ± 1.0)
        @test Δβ_twist isa Particles
        @test pmean(Δβ_twist) ≈
              twisting_birefringence(fiber, λ, T_nom; twist_rate_rad_per_m = 10.0) rtol=1e-3

        # ---- core_noncircularity and asymmetric_thermal_stress with uncertain axis_ratio
        ε = 1.02 ± 0.005
        Δβ_nc = core_noncircularity_birefringence(fiber, λ, T_nom; axis_ratio = ε)
        @test Δβ_nc isa Particles
        Δβ_ats = asymmetric_thermal_stress_birefringence(fiber, λ, T_nom; axis_ratio = ε)
        @test Δβ_ats isa Particles

        # ---- Combined: uncertain T AND uncertain bend radius (both lift together)
        Δβ_both = bending_birefringence(fiber, λ, T; bend_radius_m = R)
        @test Δβ_both isa Particles

        # ---- Guardrails
        # Invalid uncertain bend radius (mean negative) → error.
        @test_throws ArgumentError bending_birefringence(fiber, λ, T_nom;
                                                         bend_radius_m = -0.01 ± 0.001)
        # cutoff_wavelength requires scalar T (documented pmean-only).
        @test cutoff_wavelength(fiber, T_nom) isa Float64
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# ▓▓▓  path-integral.jl  ▓▓▓
# ═══════════════════════════════════════════════════════════════════════════════

if !isdefined(Main, :exp_jones_generator)
    include("../fiber-path.jl")
    include("../path-integral.jl")
end

@testset "MCM :: path-integral.jl (Frechet 4×4 exp)" begin
    # Confirm the Particles-compatible 4×4 block exp agrees with the generic exp
    # on a deterministic Float64 instance, and lifts through Particles.
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        M_det = ComplexF64[0.1+0.2im  0.3; -0.2  -0.05-0.1im]
        Mω_det = ComplexF64[0.01  0.02+0.01im; -0.01  0.005]
        A = zeros(ComplexF64, 4, 4)
        A[1:2, 1:2] = M_det;  A[1:2, 3:4] = Mω_det;  A[3:4, 3:4] = M_det
        E_ref = exp(A)[1:2, 1:2]
        F_ref = exp(A)[1:2, 3:4]
        E, F = exp_block_upper_triangular_2x2(M_det, Mω_det)
        @test isapprox(E, E_ref; atol = 1e-10)
        @test isapprox(F, F_ref; atol = 1e-10)

        # Same call with Particles-valued entries lifts cleanly.
        ε = 1.0 ± 0.01
        M_p = M_det .* ε
        Mω_p = Mω_det .* ε
        E_p, F_p = exp_block_upper_triangular_2x2(M_p, Mω_p)
        @test eltype(E_p) <: Complex{<:Particles}
        @test eltype(F_p) <: Complex{<:Particles}
        # Mean should track the deterministic result.
        @test isapprox(pmean.(real.(E_p)), real.(E_ref); atol = 1e-3)
        @test isapprox(pmean.(imag.(E_p)), imag.(E_ref); atol = 1e-3)
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end

@testset "MCM :: path-integral.jl (single-interval propagation)" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        # T-PHYSICS: pure twist at constant rate τ over length L produces
        # J = rotation by τ·L. Under uncertain τ, the mean output should track
        # the rotation for the nominal τ.
        τ_nom = 0.5
        τ = τ_nom ± 0.05
        L = 2.0
        # K(s) for a pure twist (no bend) with material-twist rate τ is a 2×2
        # skew-Hermitian with entries [0 -τ; τ 0]. Build a direct generator:
        zT = zero(τ) + 0im
        K = s -> [zT  -(τ + 0im); (τ + 0im)  zT]
        J0 = Matrix{Complex{Particles{Float64, 2000}}}(I, 2, 2)
        # Promote the identity to the Particles eltype so eltype propagation is uniform.
        J0_det = Matrix{ComplexF64}(I, 2, 2)
        J_out, stats = propagate_interval!(K, 0.0, L, J0_det; rtol = 1e-8, atol = 1e-10)
        @test eltype(J_out) <: Complex{<:Particles}
        # Expected rotation by θ = τ·L:  [cos θ  -sin θ; sin θ  cos θ]
        θ_nom = τ_nom * L
        expected = ComplexF64[cos(θ_nom)  -sin(θ_nom); sin(θ_nom)  cos(θ_nom)]
        for i in 1:2, j in 1:2
            @test pmean(real(J_out[i, j])) ≈ real(expected[i, j]) rtol=5e-2
            @test pmean(imag(J_out[i, j])) ≈ imag(expected[i, j]) rtol=5e-2 atol=1e-3
        end
        # Step controller should have made a reasonable number of accepted steps.
        @test stats.accepted_steps > 0
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end

@testset "MCM :: path-integral.jl (sensitivity propagation)" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        # Sensitivity propagation with an uncertain twist rate.
        τ_nom = 0.5
        τ = τ_nom ± 0.05
        L = 1.5
        zT = zero(τ) + 0im
        K = s -> [zT  -(τ + 0im); (τ + 0im)  zT]
        Kω = s -> zeros(ComplexF64, 2, 2)  # no frequency dispersion for pure twist
        J0 = Matrix{ComplexF64}(I, 2, 2)
        J_out, G_out, _ = propagate_interval_sensitivity!(K, Kω, 0.0, L, J0;
                                                          rtol = 1e-8, atol = 1e-10)
        @test eltype(J_out) <: Complex{<:Particles}
        @test eltype(G_out) <: Complex{<:Particles}
        # With Kω = 0 the group-delay operator G should be identically zero.
        for i in 1:2, j in 1:2
            @test pmean(abs(G_out[i, j])) < 1e-8
        end

        # 2×2 closed-form DGD: when G = 0, PMD generator H = 0, so DGD = 0.
        dgd = output_dgd_2x2(J_out, G_out)
        @test pmean(dgd) ≈ 0.0 atol=1e-6
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# ▓▓▓  path-geometry.jl  ▓▓▓
# ═══════════════════════════════════════════════════════════════════════════════

if !isdefined(Main, :BendSegment)
    include("../path-geometry.jl")
end

@testset "MCM :: path-geometry.jl (segment-level)" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        # ---- BendSegment with uncertain radius and axis_angle
        R_nom = 0.05
        R = R_nom ± 0.005
        bend = BendSegment(R, π/2, 0.1 ± 0.01)
        @test arc_length(bend) isa Particles
        @test pmean(arc_length(bend)) ≈ R_nom * (π/2) rtol=1e-3
        @test curvature(bend, 0.01) isa Particles
        # Jensen: E[1/R] > 1/E[R]; second-order bias scales as Var[R]/R³ ≈ 2% here.
        @test pmean(curvature(bend, 0.01)) ≈ 1 / R_nom rtol=5e-2
        @test tangent_local(bend, 0.01) isa Vector{<:Particles}

        # ---- BendSegment with uncertain shrinkage
        bend2 = BendSegment(0.05, π/2; shrinkage = 1.0 ± 0.01)
        @test arc_length(bend2) isa Particles
        @test pmean(arc_length(bend2)) ≈ 0.05 * (π/2) rtol=1e-3

        # ---- HelixSegment with uncertain radius and pitch
        helix = HelixSegment(0.03 ± 0.003, 0.01 ± 0.001, 2.0)
        @test arc_length(helix) isa Particles
        @test curvature(helix, 0.0) isa Particles
        @test geometric_torsion(helix, 0.0) isa Particles

        # ---- CatenarySegment with uncertain axis_angle
        cat_seg = CatenarySegment(0.1, 0.05, 0.2 ± 0.02)
        @test position_local(cat_seg, 0.025) isa Vector{<:Particles}

        # ---- StraightSegment with uncertain shrinkage
        straight = StraightSegment(0.1; shrinkage = 1.0 ± 0.005)
        @test arc_length(straight) isa Particles
        @test pmean(arc_length(straight)) ≈ 0.1 rtol=1e-3

        # ---- Float64 segments still produce Float64 output (no regression)
        bend_det = BendSegment(0.05, π/2)
        @test arc_length(bend_det) isa Float64
        @test curvature(bend_det, 0.01) isa Float64
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end

@testset "MCM :: path-geometry.jl (Path-level)" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        # ---- Single uncertain bend, queried through the Path interface
        R_nom = 0.05
        spec = PathSpec()
        bend!(spec; radius = R_nom ± 0.005, angle = π/2)
        path = build(spec)

        @test arc_length(path) isa Particles
        @test pmean(arc_length(path)) ≈ R_nom * (π/2) rtol=1e-3

        # Query at a local s well inside the segment (s = 0.01 < arc_length min)
        @test curvature(path, 0.01) isa Particles
        @test position(path, 0.01) isa Vector{<:Particles}
        @test tangent(path, 0.01) isa Vector{<:Particles}

        # ---- Mixed certain/uncertain: straight (Float64) + bend (Particles)
        spec2 = PathSpec()
        straight!(spec2; length = 0.02)
        bend!(spec2; radius = 0.05 ± 0.005, angle = π/2)
        path2 = build(spec2)
        @test length(path2.placed_segments) == 2
        # Query in the Float64 straight segment
        @test curvature(path2, 0.01) == 0.0
        # Query in the uncertain bend
        @test curvature(path2, 0.03) isa Particles

        # ---- Build-time shrinkage override with Particles
        spec3 = PathSpec()
        bend!(spec3; radius = 0.05, angle = π/2)
        path3 = build(spec3; shrinkage = 1.0 ± 0.01)
        @test arc_length(path3) isa Particles
        @test pmean(arc_length(path3)) ≈ 0.05 * (π/2) rtol=1e-3

        # ---- T-PHYSICS: straight fiber → zero curvature everywhere (no ensemble bias)
        spec4 = PathSpec()
        straight!(spec4; length = 0.1, shrinkage = 1.0 ± 0.005)
        path4 = build(spec4)
        for s in (0.01, 0.05, 0.09)
            # curvature should be identically zero (not Particles), but our ruleset promotes it.
            c = curvature(path4, s)
            @test pmean(c) ≈ 0.0 atol=1e-12
        end
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end
