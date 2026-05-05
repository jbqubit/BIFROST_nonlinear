# Pass 3 unit tests for the migrated `fiber-path.jl` and
# `fiber-path-modify.jl` API. Standalone; not yet wired into runtests.jl
# (other tests in the test/ directory are still on the old API and will
# be migrated in Pass 4).

using Test
using LinearAlgebra

if !isdefined(Main, :Fiber)
    include("../fiber-path.jl")
end
if !isdefined(Main, :modify)
    include("../fiber-path-modify.jl")
end

const _XS = FiberCrossSection(
    GermaniaSilicaGlass(0.036), GermaniaSilicaGlass(0.0),
    8.2e-6, 125e-6,
)
const _T_REF = 297.15
const _α_LIN = cte(_XS.cladding_material, _T_REF)

# ----------------------------
# Helpers
# ----------------------------

# Single-Subpath builder with a straight-line interior; seal at (0,0,L).
function _trivial_subpath(L::Float64; meta = AbstractMeta[])
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = L, meta = meta)
    jumpto!(sb; point = (0.0, 0.0, L), incoming_tangent = (0.0, 0.0, 1.0))
    return Subpath(sb)
end

# ----------------------------
# 3a. Fiber construction
# ----------------------------

@testset "Fiber — construction from SubpathBuilt" begin
    # T-GUARDRAIL
    b = build(_trivial_subpath(2.0))
    f = Fiber(b; cross_section = _XS, T_ref_K = _T_REF)
    @test f.s_start == 0.0
    @test isapprox(f.s_end, arc_length(b); atol = 1e-12)
    @test f.path === b
end

@testset "Fiber — construction from PathBuilt" begin
    # T-GUARDRAIL
    sub1 = _trivial_subpath(1.0)
    sb2 = SubpathBuilder()
    start!(sb2; point = (0.0, 0.0, 1.0), outgoing_tangent = (0.0, 0.0, 1.0))
    straight!(sb2; length = 0.5)
    jumpto!(sb2; point = (0.0, 0.0, 1.5), incoming_tangent = (0.0, 0.0, 1.0))
    p = build([sub1, Subpath(sb2)])
    f = Fiber(p; cross_section = _XS, T_ref_K = _T_REF)
    @test f.s_start == 0.0
    @test isapprox(f.s_end, arc_length(p); atol = 1e-12)
    @test f.path === p
end

@testset "Fiber — generator closures are callable" begin
    # T-GUARDRAIL
    b = build(_trivial_subpath(0.5))
    f = Fiber(b; cross_section = _XS, T_ref_K = _T_REF)
    K = generator_K(f, 1550e-9)
    Kω = generator_Kω(f, 1550e-9)
    M  = K(0.25)
    Mω = Kω(0.25)
    @test size(M) == (2, 2)
    @test size(Mω) == (2, 2)
    @test eltype(M) == ComplexF64
end

# ----------------------------
# 3b. modify() — thermal scaling on a straight segment
# ----------------------------

@testset "modify — :T_K scales StraightSegment length by (1 + α·ΔT)" begin
    # T-PHYSICS
    L = 0.5
    ΔT = 100.0
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = L, meta = [MCMadd(:T_K, ΔT)])
    jumpto!(sb; point = (0.0, 0.0, L), incoming_tangent = (0.0, 0.0, 1.0))
    f  = Fiber(build(sb); cross_section = _XS, T_ref_K = _T_REF)
    f2 = modify(f)
    expected = L * (1 + _α_LIN * ΔT)
    @test isapprox(f2.path.placed_segments[1].segment.length, expected;
                   atol = 1e-12)
    # The original Subpath length was unchanged.
    @test f.path.placed_segments[1].segment.length == L
end

# ----------------------------
# 3c. modify() — bend radius MCMadd
# ----------------------------

@testset "modify — MCMadd(:radius) on BendSegment shifts radius only" begin
    # T-PHYSICS
    R0   = 0.10
    Δr   = 0.02
    θ    = π / 3
    sb = SubpathBuilder(); start!(sb)
    bend!(sb; radius = R0, angle = θ, axis_angle = 0.0,
          meta = [MCMadd(:radius, Δr)])
    # Pin endpoint at the natural exit of an *unperturbed* bend; the modify
    # pipeline uses the same authored jumpto_point regardless.
    jumpto!(sb; point = (R0 * (1 - cos(θ)), 0.0, R0 * sin(θ)),
            incoming_tangent = (sin(θ), 0.0, cos(θ)))
    f  = Fiber(build(sb); cross_section = _XS, T_ref_K = _T_REF)
    f2 = modify(f)
    seg = f2.path.placed_segments[1].segment
    @test isapprox(seg.radius, R0 + Δr; atol = 1e-12)
    @test seg.angle == θ            # untouched
    @test seg.axis_angle == 0.0      # untouched
end

# ----------------------------
# 3d. modify() — conserve_path_length invariant
# ----------------------------

@testset "modify — conserve_path_length preserves total Subpath length" begin
    # T-PHYSICS: with conserve_path_length=true, the Subpath's total arc
    # length under modify() matches the baseline (within the connector
    # solver's tolerance), even when an interior segment is thermally
    # expanded.
    L  = 0.5
    ΔT = 100.0
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = L, meta = [MCMadd(:T_K, ΔT)])
    jumpto!(sb; point = (0.1, 0.0, L), incoming_tangent = (1.0, 0.0, 0.0),
            conserve_path_length = true)
    b = build(sb)
    L_baseline = Float64(_qc_nominalize(arc_length(b)))
    f  = Fiber(b; cross_section = _XS, T_ref_K = _T_REF)
    f2 = modify(f)
    L_modified = Float64(_qc_nominalize(arc_length(f2.path)))

    # The connector's arc-length search tolerates ~1e-6 rel_tol, so use a
    # tolerance that reflects that.
    @test isapprox(L_modified, L_baseline; atol = 1e-5 * max(1.0, L_baseline))

    # Sanity: the interior straight did expand thermally.
    expected_interior = L * (1 + _α_LIN * ΔT)
    @test isapprox(f2.path.placed_segments[1].segment.length, expected_interior;
                   atol = 1e-12)
end

@testset "modify — conserve_path_length=false grows total under :T_K" begin
    # T-GUARDRAIL: the constraint is opt-in; without it, the total Subpath
    # arc length grows when an interior segment expands thermally.
    L  = 0.5
    ΔT = 100.0
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = L, meta = [MCMadd(:T_K, ΔT)])
    jumpto!(sb; point = (0.1, 0.0, L), incoming_tangent = (1.0, 0.0, 0.0))
    b = build(sb)
    L_baseline = Float64(_qc_nominalize(arc_length(b)))
    f2 = modify(Fiber(b; cross_section = _XS, T_ref_K = _T_REF))
    L_modified = Float64(_qc_nominalize(arc_length(f2.path)))
    # The straight grew by (α·ΔT)·L; the connector solve at the unchanged
    # endpoint adjusts but doesn't compensate. Net change should be > 0.
    @test L_modified > L_baseline + 1e-9
end

# ----------------------------
# 3e. modify() — PathBuilt across Subpaths
# ----------------------------

@testset "modify — PathBuilt rebuilds each Subpath independently" begin
    # T-GUARDRAIL: a PathBuilt of two Subpaths with :T_K only on Subpath 2
    # leaves Subpath 1 unchanged after modify().
    sb1 = SubpathBuilder(); start!(sb1)
    straight!(sb1; length = 0.4)
    jumpto!(sb1; point = (0.0, 0.0, 0.4), incoming_tangent = (0.0, 0.0, 1.0))

    sb2 = SubpathBuilder()
    start!(sb2; point = (0.0, 0.0, 0.4), outgoing_tangent = (0.0, 0.0, 1.0))
    straight!(sb2; length = 0.6, meta = [MCMadd(:T_K, 50.0)])
    jumpto!(sb2; point = (0.0, 0.0, 1.0), incoming_tangent = (0.0, 0.0, 1.0))

    p  = build([Subpath(sb1), Subpath(sb2)])
    f  = Fiber(p; cross_section = _XS, T_ref_K = _T_REF)
    f2 = modify(f)
    @test f2.path isa PathBuilt
    @test length(f2.path.subpaths) == 2

    # Subpath 1 unchanged.
    seg1_orig = f.path.subpaths[1].placed_segments[1].segment
    seg1_mod  = f2.path.subpaths[1].placed_segments[1].segment
    @test seg1_orig.length == seg1_mod.length

    # Subpath 2's straight expanded.
    expected = 0.6 * (1 + _α_LIN * 50.0)
    seg2_mod = f2.path.subpaths[2].placed_segments[1].segment
    @test isapprox(seg2_mod.length, expected; atol = 1e-12)
end

# ----------------------------
# 3f. modify() — JumpBy interior segment
# ----------------------------

@testset "modify — :T_K on JumpBy scales the resolved connector by τ" begin
    # T-PHYSICS: JumpBy's resolved QuinticConnector arc length scales by
    # (1 + α_lin·ΔT) under :T_K.
    ΔT = 200.0
    sb = SubpathBuilder(); start!(sb)
    jumpby!(sb; delta = (0.0, 0.0, 0.4), meta = [MCMadd(:T_K, ΔT)])
    jumpto!(sb; point = (0.0, 0.0, 0.4), incoming_tangent = (0.0, 0.0, 1.0))
    b = build(sb)
    L_jumpby_baseline = Float64(_qc_nominalize(
        arc_length(b.placed_segments[1].segment)))

    f  = Fiber(b; cross_section = _XS, T_ref_K = _T_REF)
    f2 = modify(f)
    L_jumpby_modified = Float64(_qc_nominalize(
        arc_length(f2.path.placed_segments[1].segment)))

    expected = L_jumpby_baseline * (1 + _α_LIN * ΔT)
    @test isapprox(L_jumpby_modified, expected; rtol = 1e-12)
end
