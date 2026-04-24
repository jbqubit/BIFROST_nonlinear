using Test
using LinearAlgebra

if !isdefined(Main, :PathSpec)
    include("../path-geometry.jl")
end
if !isdefined(Main, :Fiber)
    include("../fiber-path.jl")
end
if !isdefined(Main, :modify)
    include("../fiber-path-modify.jl")
end

# -----------------------------------------------------------------------
# Test setup: a small helper to build a Fiber with a pure-silica cladding
# (so α_lin = SILICA_CTE) and a reference temperature.
# -----------------------------------------------------------------------

const _MODIFY_TEST_XS = FiberCrossSection(
    GermaniaSilicaGlass(0.036),
    GermaniaSilicaGlass(0.0),     # pure silica cladding → α_lin = SILICA_CTE
    8.2e-6,
    125e-6,
)

const _MODIFY_TEST_T_REF = 297.15
const _MODIFY_ALPHA_LIN  = cte(_MODIFY_TEST_XS.cladding_material, _MODIFY_TEST_T_REF)

# Helper: ΔT that produces a given total-length α = 1 + α_lin·ΔT.
_ΔT_for(α) = (α - 1) / _MODIFY_ALPHA_LIN
_mcm(α)    = AbstractMeta[MCMadd(:T_K, _ΔT_for(α))]

_modify_fiber_path(path) = modify(Fiber(path;
    cross_section = _MODIFY_TEST_XS, T_ref_K = _MODIFY_TEST_T_REF))

# -----------------------------------------------------------------------
# Uniform thermal shrinkage via :T_K
# -----------------------------------------------------------------------

@testset "modify — :T_K uniform scales path length" begin
    # T-PHYSICS: uniform α on every segment multiplies total arc length by α
    spec = PathSpec()
    straight!(spec; length = 1.0, meta = _mcm(1.05))
    bend!(spec; radius = 0.1, angle = π / 2, meta = _mcm(1.05))
    path = build(spec)

    path_scaled = _modify_fiber_path(path)
    @test path_length(path_scaled) ≈ path_length(path) * 1.05 atol = 1e-10
end

@testset "modify — :T_K preserves joint tangent continuity" begin
    # T-GUARDRAIL
    spec = PathSpec()
    straight!(spec; length = 0.5, meta = _mcm(0.98))
    bend!(spec; radius = 0.1, angle = π / 3, meta = _mcm(0.98))
    path = _modify_fiber_path(build(spec))

    ps = path.placed_segments
    for i in 1:(length(ps) - 1)
        s_joint = ps[i + 1].s_offset_eff
        T_before = tangent(path, s_joint - 1e-9)
        T_after  = tangent(path, s_joint + 1e-9)
        @test norm(T_before - T_after) < 1e-6
    end
end

@testset "modify — :T_K on BendSegment preserves angle, scales radius" begin
    # T-PHYSICS: α scales R_eff = α·R; swept angle preserved.
    spec = PathSpec()
    bend!(spec; radius = 0.1, angle = π / 3, axis_angle = 0.0, meta = _mcm(1.1))
    path1 = build(spec)
    path2 = _modify_fiber_path(path1)

    seg1 = path1.placed_segments[1].segment
    seg2 = path2.placed_segments[1].segment
    @test seg2.radius ≈ 1.1 * seg1.radius
    @test seg2.angle == seg1.angle
    @test arc_length(seg2) ≈ 1.1 * arc_length(seg1)
end

@testset "modify — :T_K on CatenarySegment scales arc length and parameter a" begin
    # T-PHYSICS
    spec = PathSpec()
    catenary!(spec; a = 0.2, length = 1.0, meta = _mcm(0.95))
    path = build(spec)
    path_s = _modify_fiber_path(path)

    seg = path.placed_segments[1].segment
    seg_s = path_s.placed_segments[1].segment
    @test seg_s.a ≈ 0.95 * seg.a
    @test seg_s.length ≈ 0.95 * seg.length
    @test arc_length(seg_s) ≈ 0.95 * arc_length(seg)
end

@testset "modify — :T_K on HelixSegment scales arc length but preserves turns" begin
    spec = PathSpec()
    helix!(spec; radius = 0.03, pitch = 0.01, turns = 2.0, meta = _mcm(0.9))
    path = build(spec)
    path_s = _modify_fiber_path(path)

    seg = path.placed_segments[1].segment
    seg_s = path_s.placed_segments[1].segment
    @test seg_s.turns == seg.turns
    @test seg_s.radius ≈ 0.9 * seg.radius
    @test seg_s.pitch ≈ 0.9 * seg.pitch
    @test arc_length(seg_s) ≈ 0.9 * arc_length(seg)
end

@testset "modify — :T_K on JumpBy endpoint scales" begin
    # T-PHYSICS: uniform α scales all positions
    spec = PathSpec()
    jumpby!(spec; delta = (0.0, 0.0, 1.0), meta = _mcm(0.8))
    path = build(spec)
    path_s = _modify_fiber_path(path)
    @test end_point(path_s) ≈ [0.0, 0.0, 0.8] atol = 1e-10
end

# -----------------------------------------------------------------------
# Per-segment independence (meta on some segments, not others)
# -----------------------------------------------------------------------

@testset "modify — per-segment scaling applies independently" begin
    spec = PathSpec()
    straight!(spec; length = 1.0, meta = _mcm(1.0))   # α = 1 via ΔT = 0
    straight!(spec; length = 1.0, meta = _mcm(0.5))
    path = build(spec)

    path_s = _modify_fiber_path(path)
    @test path_length(path_s) ≈ 1.5 atol = 1e-10
    @test arc_length(path_s.placed_segments[1].segment) ≈ 1.0
    @test arc_length(path_s.placed_segments[2].segment) ≈ 0.5
end

@testset "modify — segments without MCM annotations default to α = 1.0" begin
    spec = PathSpec()
    straight!(spec; length = 1.0)                     # no meta → α = 1
    straight!(spec; length = 2.0, meta = _mcm(0.5))
    path = build(spec)

    path_s = _modify_fiber_path(path)
    @test arc_length(path_s.placed_segments[1].segment) ≈ 1.0
    @test arc_length(path_s.placed_segments[2].segment) ≈ 1.0
end

@testset "modify — non-matching MCM symbols are ignored" begin
    spec = PathSpec()
    straight!(spec; length = 1.0,
              meta = AbstractMeta[MCMadd(:not_a_field, 1e6),
                                  Nickname("labelled")])
    path = build(spec)

    path_s = _modify_fiber_path(path)
    @test arc_length(path_s.placed_segments[1].segment) ≈ 1.0
end

# -----------------------------------------------------------------------
# Direct field-level MCMs (new in fiber-path-modify.jl)
# -----------------------------------------------------------------------

@testset "modify — direct field :length on StraightSegment" begin
    spec = PathSpec()
    straight!(spec; length = 1.0, meta = AbstractMeta[MCMadd(:length, 0.05)])
    path = build(spec)
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.length ≈ 1.05 atol = 1e-12
end

@testset "modify — direct field :radius on BendSegment" begin
    spec = PathSpec()
    bend!(spec; radius = 0.05, angle = π / 2,
          meta = AbstractMeta[MCMadd(:radius, 0.01)])
    path = build(spec)
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.radius ≈ 0.06 atol = 1e-12
    @test seg.angle ≈ π / 2
end

@testset "modify — direct field :angle on BendSegment" begin
    spec = PathSpec()
    bend!(spec; radius = 0.05, angle = π / 2,
          meta = AbstractMeta[MCMadd(:angle, π / 12)])
    path = build(spec)
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.angle ≈ π / 2 + π / 12 atol = 1e-12
    @test seg.radius ≈ 0.05
end

@testset "modify — direct field :pitch on HelixSegment (MCMmul)" begin
    spec = PathSpec()
    helix!(spec; radius = 0.03, pitch = 0.01, turns = 2.0,
           meta = AbstractMeta[MCMmul(:pitch, 1.1)])
    path = build(spec)
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.pitch ≈ 0.01 * 1.1 atol = 1e-14
    @test seg.radius ≈ 0.03
end

@testset "modify — direct field :a on CatenarySegment" begin
    spec = PathSpec()
    catenary!(spec; a = 0.2, length = 1.0,
              meta = AbstractMeta[MCMadd(:a, 0.002)])
    path = build(spec)
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.a ≈ 0.202 atol = 1e-12
    @test seg.length ≈ 1.0
end

@testset "modify — combined :T_K and direct :length on StraightSegment" begin
    # Order: T_K first (multiplicative length scaling via CTE), then the direct
    # :length additive offset on top.
    spec = PathSpec()
    straight!(spec; length = 1.0,
              meta = AbstractMeta[MCMadd(:T_K, 10.0),
                                  MCMadd(:length, 0.001)])
    path = build(spec)
    seg = _modify_fiber_path(path).placed_segments[1].segment
    expected = 1.0 * (1 + _MODIFY_ALPHA_LIN * 10.0) + 0.001
    @test seg.length ≈ expected atol = 1e-14
end

# -----------------------------------------------------------------------
# Twist overlay remapping
# -----------------------------------------------------------------------

@testset "modify — :T_K preserves constant twist rate; total scales with length" begin
    # T-PHYSICS: under a geometric shrinkage, τ (rad/m) stays constant and
    # the total ∫τ ds scales by α.
    α = 0.75
    spec = PathSpec()
    straight!(spec; length = 2.0, meta = _mcm(α))
    twist!(spec; s_start = 0.0, length = 2.0, rate = 1.5)
    path = build(spec)
    path_s = _modify_fiber_path(path)

    total_orig   = total_material_twist(path)
    total_shrunk = total_material_twist(path_s)
    @test total_shrunk ≈ α * total_orig atol = 1e-10

    # Rate at any shrunk-arc-length point equals the original rate.
    @test material_twist(path_s, 0.5)  ≈ 1.5 atol = 1e-12
    @test material_twist(path_s, 1.0)  ≈ 1.5 atol = 1e-12
end

# -----------------------------------------------------------------------
# MCMadd: Particles flow through modify via a Particles-valued ΔT
# -----------------------------------------------------------------------

if !isdefined(Main, :Particles)
    using MonteCarloMeasurements
end

@testset "modify — MCMadd Particles ΔT scales arc length" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        # Particles-valued ΔT with zero mean → nominal arc length preserved,
        # output lifts to Particles.
        ΔT = 0.0 ± (0.01 / _MODIFY_ALPHA_LIN)    # σ_α = 0.01 around α = 1
        spec = PathSpec()
        bend!(spec; radius = 0.05, angle = π / 2,
              meta = AbstractMeta[MCMadd(:T_K, ΔT)])
        path = build(spec)
        path_s = _modify_fiber_path(path)

        @test arc_length(path_s) isa Particles
        @test pmean(arc_length(path_s)) ≈ 0.05 * (π / 2) rtol = 1e-3
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end
