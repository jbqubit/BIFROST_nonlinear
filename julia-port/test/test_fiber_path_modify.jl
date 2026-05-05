using Test
using LinearAlgebra

if !isdefined(Main, :SubpathBuilt)
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
_mcm(α)    = [MCMadd(:T_K, _ΔT_for(α))]

# Trial-build a SubpathBuilder, read the natural exit position/tangent, and
# seal the builder with a jumpto! at that point. Used by tests whose paths
# don't naturally have a global-anchored endpoint.
function _seal_natural!(sb::SubpathBuilder)
    @assert isnothing(sb.jumpto_point) "_seal_natural!: builder already sealed"
    tmp = deepcopy(sb)
    jumpto!(tmp; point = (1e9, 1e9, 1e9))
    b_tmp = build(Subpath(tmp))
    s_end_int = Float64(_qc_nominalize(b_tmp.jumpto_placed.s_offset_eff))
    if s_end_int <= 0.0
        natural_pos = collect(sb.start_point::NTuple{3, Float64})
        natural_tan = collect(sb.start_outgoing_tangent::NTuple{3, Float64})
    else
        natural_pos = collect(position(b_tmp, s_end_int))
        natural_tan = collect(tangent(b_tmp, s_end_int))
    end
    jumpto!(sb;
        point = (natural_pos[1], natural_pos[2], natural_pos[3]),
        incoming_tangent = (natural_tan[1], natural_tan[2], natural_tan[3]),
    )
    return sb
end

# Build a sealed SubpathBuilt from a do-block that authors interior
# segments. Seals at the natural exit if no jumpto! was called.
function _build_path(f::Function)
    sb = SubpathBuilder(); start!(sb)
    f(sb)
    if isnothing(sb.jumpto_point)
        _seal_natural!(sb)
    end
    return build(Subpath(sb))
end

_modify_fiber_path(path) = modify(Fiber(path;
    cross_section = _MODIFY_TEST_XS, T_ref_K = _MODIFY_TEST_T_REF)).path

# -----------------------------------------------------------------------
# Uniform thermal scaling via :T_K
# -----------------------------------------------------------------------

@testset "modify — :T_K uniform scales path length" begin
    # T-PHYSICS: uniform α on every segment multiplies total interior arc
    # length by α. The terminal connector (degenerate) is essentially
    # zero-length; check arc length of the interior segments only.
    α = 1.05
    path = _build_path() do sb
        straight!(sb; length = 1.0, meta = _mcm(α))
        bend!(sb; radius = 0.1, angle = π / 2, meta = _mcm(α))
    end
    interior_baseline = sum(arc_length(ps.segment) for ps in path.placed_segments)
    path_scaled = _modify_fiber_path(path)
    interior_scaled = sum(arc_length(ps.segment)
                          for ps in path_scaled.placed_segments)
    @test interior_scaled ≈ α * interior_baseline atol = 1e-9
end

@testset "modify — :T_K preserves joint tangent continuity" begin
    # T-GUARDRAIL
    path = _build_path() do sb
        straight!(sb; length = 0.5, meta = _mcm(0.98))
        bend!(sb; radius = 0.1, angle = π / 3, meta = _mcm(0.98))
    end
    path_s = _modify_fiber_path(path)
    ps = path_s.placed_segments
    for i in 1:(length(ps) - 1)
        s_joint = Float64(_qc_nominalize(ps[i + 1].s_offset_eff))
        T_before = tangent(path_s, s_joint - 1e-9)
        T_after  = tangent(path_s, s_joint + 1e-9)
        @test norm(T_before - T_after) < 1e-6
    end
end

@testset "modify — :T_K on BendSegment preserves angle, scales radius" begin
    # T-PHYSICS: α scales R_eff = α·R; swept angle preserved.
    path1 = _build_path() do sb
        bend!(sb; radius = 0.1, angle = π / 3, axis_angle = 0.0, meta = _mcm(1.1))
    end
    path2 = _modify_fiber_path(path1)
    seg1 = path1.placed_segments[1].segment
    seg2 = path2.placed_segments[1].segment
    @test seg2.radius ≈ 1.1 * seg1.radius
    @test seg2.angle == seg1.angle
    @test arc_length(seg2) ≈ 1.1 * arc_length(seg1)
end

@testset "modify — :T_K on CatenarySegment scales arc length and parameter a" begin
    # T-PHYSICS
    path = _build_path() do sb
        catenary!(sb; a = 0.2, length = 1.0, meta = _mcm(0.95))
    end
    path_s = _modify_fiber_path(path)
    seg = path.placed_segments[1].segment
    seg_s = path_s.placed_segments[1].segment
    @test seg_s.a ≈ 0.95 * seg.a
    @test seg_s.length ≈ 0.95 * seg.length
    @test arc_length(seg_s) ≈ 0.95 * arc_length(seg)
end

@testset "modify — :T_K on HelixSegment scales arc length but preserves turns" begin
    path = _build_path() do sb
        helix!(sb; radius = 0.03, pitch = 0.01, turns = 2.0, meta = _mcm(0.9))
    end
    path_s = _modify_fiber_path(path)
    seg = path.placed_segments[1].segment
    seg_s = path_s.placed_segments[1].segment
    @test seg_s.turns == seg.turns
    @test seg_s.radius ≈ 0.9 * seg.radius
    @test seg_s.pitch ≈ 0.9 * seg.pitch
    @test arc_length(seg_s) ≈ 0.9 * arc_length(seg)
end

@testset "modify — :T_K on JumpBy scales connector arc length" begin
    # T-PHYSICS: JumpBy's resolved QuinticConnector arc length scales by α.
    path = _build_path() do sb
        jumpby!(sb; delta = (0.0, 0.0, 1.0), meta = _mcm(0.8))
    end
    L_jumpby_baseline = arc_length(path.placed_segments[1].segment)
    path_s = _modify_fiber_path(path)
    L_jumpby_modified = arc_length(path_s.placed_segments[1].segment)
    @test L_jumpby_modified ≈ 0.8 * L_jumpby_baseline rtol = 1e-12
end

# -----------------------------------------------------------------------
# Per-segment independence (meta on some segments, not others)
# -----------------------------------------------------------------------

@testset "modify — per-segment scaling applies independently" begin
    path = _build_path() do sb
        straight!(sb; length = 1.0, meta = _mcm(1.0))   # α = 1 via ΔT = 0
        straight!(sb; length = 1.0, meta = _mcm(0.5))
    end
    path_s = _modify_fiber_path(path)
    @test arc_length(path_s.placed_segments[1].segment) ≈ 1.0
    @test arc_length(path_s.placed_segments[2].segment) ≈ 0.5
end

@testset "modify — segments without MCM annotations default to α = 1.0" begin
    path = _build_path() do sb
        straight!(sb; length = 1.0)                     # no meta → α = 1
        straight!(sb; length = 2.0, meta = _mcm(0.5))
    end
    path_s = _modify_fiber_path(path)
    @test arc_length(path_s.placed_segments[1].segment) ≈ 1.0
    @test arc_length(path_s.placed_segments[2].segment) ≈ 1.0
end

@testset "modify — non-matching MCM symbols are ignored" begin
    path = _build_path() do sb
        straight!(sb; length = 1.0,
                  meta = [MCMadd(:not_a_field, 1e6),
                          Nickname("labelled")])
    end
    path_s = _modify_fiber_path(path)
    @test arc_length(path_s.placed_segments[1].segment) ≈ 1.0
end

# -----------------------------------------------------------------------
# Direct field-level MCMs
# -----------------------------------------------------------------------

@testset "modify — direct field :length on StraightSegment" begin
    path = _build_path() do sb
        straight!(sb; length = 1.0, meta = [MCMadd(:length, 0.05)])
    end
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.length ≈ 1.05 atol = 1e-12
end

@testset "modify — direct field :radius on BendSegment" begin
    path = _build_path() do sb
        bend!(sb; radius = 0.05, angle = π / 2,
              meta = [MCMadd(:radius, 0.01)])
    end
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.radius ≈ 0.06 atol = 1e-12
    @test seg.angle ≈ π / 2
end

@testset "modify — direct field :angle on BendSegment" begin
    path = _build_path() do sb
        bend!(sb; radius = 0.05, angle = π / 2,
              meta = [MCMadd(:angle, π / 12)])
    end
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.angle ≈ π / 2 + π / 12 atol = 1e-12
    @test seg.radius ≈ 0.05
end

@testset "modify — direct field :pitch on HelixSegment (MCMmul)" begin
    path = _build_path() do sb
        helix!(sb; radius = 0.03, pitch = 0.01, turns = 2.0,
               meta = [MCMmul(:pitch, 1.1)])
    end
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.pitch ≈ 0.01 * 1.1 atol = 1e-14
    @test seg.radius ≈ 0.03
end

@testset "modify — direct field :a on CatenarySegment" begin
    path = _build_path() do sb
        catenary!(sb; a = 0.2, length = 1.0,
                  meta = [MCMadd(:a, 0.002)])
    end
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.a ≈ 0.202 atol = 1e-12
    @test seg.length ≈ 1.0
end

@testset "modify — combined :T_K and direct :length on StraightSegment" begin
    # T_K first (multiplicative length scaling via CTE), then the direct
    # :length additive offset on top.
    path = _build_path() do sb
        straight!(sb; length = 1.0,
                  meta = [MCMadd(:T_K, 10.0), MCMadd(:length, 0.001)])
    end
    seg = _modify_fiber_path(path).placed_segments[1].segment
    expected = 1.0 * (1 + _MODIFY_ALPHA_LIN * 10.0) + 0.001
    @test seg.length ≈ expected atol = 1e-14
end

# -----------------------------------------------------------------------
# Twist overlay remapping
# -----------------------------------------------------------------------

# TODO: twist refactor — pending per-segment-meta twist subsystem.
@testset "modify — :T_K preserves constant twist rate; total scales with length" begin
    @test_skip true
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
        # Particles-valued ΔT with zero mean → nominal interior arc length
        # preserved; output segment length lifts to Particles.
        ΔT = 0.0 ± (0.01 / _MODIFY_ALPHA_LIN)    # σ_α = 0.01 around α = 1
        path = _build_path() do sb
            bend!(sb; radius = 0.05, angle = π / 2,
                  meta = [MCMadd(:T_K, ΔT)])
        end
        path_s = _modify_fiber_path(path)

        seg_s = path_s.placed_segments[1].segment
        @test arc_length(seg_s) isa Particles
        @test pmean(arc_length(seg_s)) ≈ 0.05 * (π / 2) rtol = 1e-3
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end

# -----------------------------------------------------------------------
# QuinticConnector :T_K scaling
# -----------------------------------------------------------------------

@testset "modify — :T_K scales QuinticConnector arc length (JumpBy)" begin
    # A JumpBy with :T_K MCM meta must scale through to the resolved
    # QuinticConnector after build(). The geometric scale τ should multiply
    # the coefficient matrix and the arc-length table uniformly so the
    # JumpBy segment length scales by τ.
    path = _build_path() do sb
        jumpby!(sb; delta = (0.4, 0.0, 0.4),
                tangent = (1.0, 0.0, 0.0),
                meta = _mcm(0.5))
    end
    L_before = arc_length(path.placed_segments[1].segment)
    path_s = _modify_fiber_path(path)
    L_after = arc_length(path_s.placed_segments[1].segment)
    @test isapprox(L_after, 0.5 * L_before; rtol = 1e-10, atol = 1e-12)
end

@testset "modify — T-GUARDRAIL: upstream bend changes recompute connector K0" begin
    path_s = _modify_fiber_path(_build_path() do sb
        bend!(sb; radius = 1.0, angle = π / 3,
              meta = [MCMmul(:radius, 2.0)])
        jumpby!(sb; delta = (1.0, 0.0, 0.2))
    end)
    bend_seg = path_s.placed_segments[1].segment
    connector = path_s.placed_segments[2].segment

    @test curvature(bend_seg, arc_length(bend_seg)) ≈ 0.5 atol = 1e-12
    @test curvature(connector, 0.0) ≈ 0.5 atol = 1e-12
end

@testset "modify — T-GUARDRAIL: Twist anchors tolerate MCM-valued modified length" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        path = _build_path() do sb
            straight!(sb; length = 1.0,
                      meta = [
                          Twist(; rate = 2.0),
                          MCMadd(:length, 0.0 ± 0.01),
                      ])
        end
        path_s = _modify_fiber_path(path)
        seg_s = path_s.placed_segments[1].segment
        L_seg = arc_length(seg_s)
        @test L_seg isa Particles
        # Twist run extends from the segment's start over its full length.
        # Integrated material twist over the segment = 2.0 * L_seg.
        @test total_material_twist(path_s; s_start = 0.0,
                                   s_end = Float64(_qc_nominalize(L_seg))) ≈
              2.0 * pmean(L_seg) rtol = 1e-3
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end

# -----------------------------------------------------------------------
# conserve_path_length on a Subpath's terminal jumpto
# -----------------------------------------------------------------------
# Pass 1/3 introduced `conserve_path_length=true` as the new mechanism for
# pinning a Subpath's total arc length under meta-induced length changes.
# Pass 3's test_fiber_path_pass3.jl covers this in detail; the rest of the
# old `:T_K`-on-JumpTo testsets have been deleted.
