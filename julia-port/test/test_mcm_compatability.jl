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

# (pending Stage 1 Part B)

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
