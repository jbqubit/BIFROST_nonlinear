using Test
using LinearAlgebra

if !isdefined(Main, :PathSpec)
    include("../path-geometry.jl")
end
if !isdefined(Main, :shrink)
    include("../fiber-path-shrinkage.jl")
end

# -----------------------------------------------------------------------
# Uniform shrinkage
# -----------------------------------------------------------------------

@testset "shrink — uniform scales path length" begin
    # T-PHYSICS: uniform scaling α multiplies total arc length by α
    spec = PathSpec()
    straight!(spec; length = 1.0)
    bend!(spec; radius = 0.1, angle = π / 2)
    path = build(spec)

    path_scaled = shrink(path, 1.05)
    @test path_length(path_scaled) ≈ path_length(path) * 1.05 atol = 1e-10
end

@testset "shrink — preserves joint tangent continuity" begin
    # T-GUARDRAIL
    spec = PathSpec()
    straight!(spec; length = 0.5)
    bend!(spec; radius = 0.1, angle = π / 3)
    path = shrink(build(spec), 0.98)

    ps = path.placed_segments
    for i in 1:(length(ps) - 1)
        s_joint = ps[i + 1].s_offset_eff
        T_before = tangent(path, s_joint - 1e-9)
        T_after  = tangent(path, s_joint + 1e-9)
        @test norm(T_before - T_after) < 1e-6
    end
end

@testset "shrink — BendSegment preserves angle, scales radius" begin
    # T-PHYSICS: shrink α scales R_eff = α·R; swept angle preserved.
    spec = PathSpec()
    bend!(spec; radius = 0.1, angle = π / 3, axis_angle = 0.0)
    path1 = build(spec)
    path2 = shrink(path1, 1.1)

    seg1 = path1.placed_segments[1].segment
    seg2 = path2.placed_segments[1].segment
    @test seg2.radius ≈ 1.1 * seg1.radius
    @test seg2.angle == seg1.angle
    @test arc_length(seg2) ≈ 1.1 * arc_length(seg1)
end

@testset "shrink — CatenarySegment scales arc length and parameter a" begin
    # T-PHYSICS
    spec = PathSpec()
    catenary!(spec; a = 0.2, length = 1.0)
    path = build(spec)
    path_s = shrink(path, 0.95)

    seg = path.placed_segments[1].segment
    seg_s = path_s.placed_segments[1].segment
    @test seg_s.a ≈ 0.95 * seg.a
    @test seg_s.length ≈ 0.95 * seg.length
    @test arc_length(seg_s) ≈ 0.95 * arc_length(seg)
end

@testset "shrink — HelixSegment scales arc length but preserves turns" begin
    spec = PathSpec()
    helix!(spec; radius = 0.03, pitch = 0.01, turns = 2.0)
    path = build(spec)
    path_s = shrink(path, 0.9)

    seg = path.placed_segments[1].segment
    seg_s = path_s.placed_segments[1].segment
    @test seg_s.turns == seg.turns
    @test seg_s.radius ≈ 0.9 * seg.radius
    @test seg_s.pitch ≈ 0.9 * seg.pitch
    @test arc_length(seg_s) ≈ 0.9 * arc_length(seg)
end

@testset "shrink — JumpBy endpoint scales" begin
    # T-PHYSICS: uniform α scales all positions
    spec = PathSpec()
    jumpby!(spec; delta = (0.0, 0.0, 1.0))
    path = build(spec)
    path_s = shrink(path, 0.8)
    @test end_point(path_s) ≈ [0.0, 0.0, 0.8] atol = 1e-10
end

# -----------------------------------------------------------------------
# Per-index shrinkage
# -----------------------------------------------------------------------

@testset "shrink — per-index scaling applies independently" begin
    spec = PathSpec()
    straight!(spec; length = 1.0)
    straight!(spec; length = 1.0)
    path = build(spec)

    path_s = shrink(path, Dict(1 => 1.0, 2 => 0.5))
    @test path_length(path_s) ≈ 1.5 atol = 1e-10

    # Segment 1 unchanged; segment 2 halved
    @test arc_length(path_s.placed_segments[1].segment) ≈ 1.0
    @test arc_length(path_s.placed_segments[2].segment) ≈ 0.5
end

@testset "shrink — per-index missing keys default to 1.0" begin
    spec = PathSpec()
    straight!(spec; length = 1.0)
    straight!(spec; length = 2.0)
    path = build(spec)

    path_s = shrink(path, Dict(2 => 0.5))  # segment 1 absent → α = 1
    @test arc_length(path_s.placed_segments[1].segment) ≈ 1.0
    @test arc_length(path_s.placed_segments[2].segment) ≈ 1.0
end

# -----------------------------------------------------------------------
# Twist overlay remapping under shrinkage
# -----------------------------------------------------------------------

@testset "shrink — constant twist rate preserved; total scales with length" begin
    # T-PHYSICS (chosen semantic): the material twist rate τ (rad/m) is a
    # property of the fiber itself, not a baked-in turn count. Under a shrinkage
    # that only scales *geometric* arc length, τ stays constant and the total
    # accumulated twist ∫τ ds scales by α.
    spec = PathSpec()
    straight!(spec; length = 2.0)
    twist!(spec; s_start = 0.0, length = 2.0, rate = 1.5)
    path = build(spec)
    α = 0.75
    path_s = shrink(path, α)

    total_orig   = total_material_twist(path)
    total_shrunk = total_material_twist(path_s)
    @test total_shrunk ≈ α * total_orig atol = 1e-10

    # Rate at any shrunk-arc-length point equals the original rate.
    @test material_twist(path_s, 0.5)  ≈ 1.5 atol = 1e-12
    @test material_twist(path_s, 1.0)  ≈ 1.5 atol = 1e-12
end

# -----------------------------------------------------------------------
# MCM: Particles flow through shrink
# -----------------------------------------------------------------------

if !isdefined(Main, :Particles)
    using MonteCarloMeasurements
end

@testset "shrink — MCM Particles α scales arc length" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        spec = PathSpec()
        bend!(spec; radius = 0.05, angle = π / 2)
        path = build(spec)
        α = 1.0 ± 0.01
        path_s = shrink(path, α)

        @test arc_length(path_s) isa Particles
        @test pmean(arc_length(path_s)) ≈ 0.05 * (π / 2) rtol = 1e-3
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end
