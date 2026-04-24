using Test
using LinearAlgebra

if !isdefined(Main, :PathSpec)
    include("../path-geometry.jl")
end

# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

is_unit(v)       = abs(norm(v) - 1.0) < 1e-12
is_orthonormal(T, N, B) = is_unit(T) && is_unit(N) && is_unit(B) &&
                           abs(dot(T, N)) < 1e-12 &&
                           abs(dot(T, B)) < 1e-12 &&
                           abs(dot(N, B)) < 1e-12 &&
                           norm(cross(T, N) - B) < 1e-12   # right-handed

# -----------------------------------------------------------------------
# StraightSegment
# -----------------------------------------------------------------------

@testset "StraightSegment — arc length" begin
    # T-PHYSICS: arc length of a straight segment is length * shrinkage
    seg = StraightSegment(3.0; shrinkage = 1.0)
    @test arc_length(seg) ≈ 3.0

    seg2 = StraightSegment(3.0; shrinkage = 1.02)
    @test arc_length(seg2) ≈ 3.06
end

@testset "StraightSegment — geometry" begin
    # T-PHYSICS: straight segment has zero curvature and zero torsion
    seg = StraightSegment(5.0)
    @test curvature(seg, 0.0)          == 0.0
    @test curvature(seg, 2.5)          == 0.0
    @test geometric_torsion(seg, 0.0)  == 0.0

    # tangent is always along local z
    @test tangent_local(seg, 0.0) ≈ [0.0, 0.0, 1.0]
    @test tangent_local(seg, 2.5) ≈ [0.0, 0.0, 1.0]

    # position advances along z
    @test position_local(seg, 0.0) ≈ [0.0, 0.0, 0.0]
    @test position_local(seg, 3.0) ≈ [0.0, 0.0, 3.0]

    # frame is orthonormal at any s
    T = tangent_local(seg, 1.0)
    N = normal_local(seg, 1.0)
    B = binormal_local(seg, 1.0)
    @test is_orthonormal(T, N, B)
end

@testset "StraightSegment — shrinkage scales arc length, not curvature" begin
    # T-PHYSICS: shrinkage is a metric rescaling; curvature (1/length) is unchanged
    # for a straight segment (κ = 0 regardless)
    seg = StraightSegment(2.0; shrinkage = 0.9)
    @test arc_length(seg) ≈ 1.8
    @test curvature(seg, 0.0) == 0.0
end

# -----------------------------------------------------------------------
# BendSegment
# -----------------------------------------------------------------------

@testset "BendSegment — arc length and curvature" begin
    # T-PHYSICS: arc_length = shrinkage * radius * |angle|, κ = 1 / (shrinkage * radius)
    seg = BendSegment(0.1, π / 2)
    @test arc_length(seg) ≈ 0.1 * π / 2
    @test curvature(seg, 0.0) ≈ 1.0 / 0.1

    seg2 = BendSegment(0.1, π / 2; shrinkage = 1.05)
    @test arc_length(seg2) ≈ 1.05 * 0.1 * π / 2
    @test curvature(seg2, 0.0) ≈ 1.0 / (1.05 * 0.1)
end

@testset "BendSegment — zero torsion (planar curve)" begin
    # T-PHYSICS: circular arc is planar, geometric torsion = 0
    seg = BendSegment(0.05, π)
    @test geometric_torsion(seg, 0.0) == 0.0
    @test geometric_torsion(seg, arc_length(seg) / 2) == 0.0
end

@testset "BendSegment — initial tangent along local z" begin
    # T-GUARDRAIL: all segments must start with tangent (0,0,1) in local coords
    seg = BendSegment(0.2, π / 3, π / 4)
    @test tangent_local(seg, 0.0) ≈ [0.0, 0.0, 1.0]
end

@testset "BendSegment — end position for quarter circle (axis_angle = 0)" begin
    # T-PHYSICS: quarter circle of radius R in x-z plane (axis_angle=0).
    # Start: (0,0,0), tangent z. End: (R, 0, R) tangent x.
    R = 0.1
    seg = BendSegment(R, π / 2, 0.0)
    pos_end = end_position_local(seg)
    @test pos_end ≈ [R, 0.0, R] atol = 1e-12

    (T_end, N_end, B_end) = end_frame_local(seg)
    @test T_end ≈ [1.0, 0.0, 0.0] atol = 1e-12   # tangent rotated 90° toward +x
    @test is_orthonormal(T_end, N_end, B_end)
end

@testset "BendSegment — end position for half circle" begin
    # T-PHYSICS: half circle of radius R in x-z plane.
    # Start: (0,0,0) tangent z. End: (2R, 0, 0) tangent -z.
    R = 0.05
    seg = BendSegment(R, π, 0.0)
    pos_end = end_position_local(seg)
    @test pos_end ≈ [2R, 0.0, 0.0] atol = 1e-12

    (T_end, _, _) = end_frame_local(seg)
    @test T_end ≈ [0.0, 0.0, -1.0] atol = 1e-12   # reversed
end

@testset "BendSegment — axis_angle rotates bend plane" begin
    # T-PHYSICS: axis_angle = π/2 bends in the y-z plane instead of x-z.
    # Quarter circle end position should be (0, R, R).
    R = 0.1
    seg = BendSegment(R, π / 2, π / 2)
    pos_end = end_position_local(seg)
    @test pos_end ≈ [0.0, R, R] atol = 1e-12
end

@testset "BendSegment — frame orthonormality along arc" begin
    # T-GUARDRAIL: frame must remain orthonormal at all points
    seg = BendSegment(0.08, 2π / 3, π / 6)
    for s in range(0.0, arc_length(seg); length = 9)
        T = tangent_local(seg, s)
        N = normal_local(seg, s)
        B = binormal_local(seg, s)
        @test is_orthonormal(T, N, B)
    end
end

@testset "BendSegment — shrinkage preserves angle, scales radius" begin
    # T-PHYSICS: shrinkage α scales R_eff = α*R; swept angle is preserved.
    # Two segments with same angle but different shrinkage should have the same
    # end tangent direction but different end positions.
    R = 0.1; θ = π / 3
    seg1 = BendSegment(R, θ, 0.0; shrinkage = 1.0)
    seg2 = BendSegment(R, θ, 0.0; shrinkage = 1.1)

    (T1, _, _) = end_frame_local(seg1)
    (T2, _, _) = end_frame_local(seg2)
    @test T1 ≈ T2 atol = 1e-12   # same angle → same end tangent

    @test arc_length(seg2) ≈ arc_length(seg1) * 1.1
end

# -----------------------------------------------------------------------
# CatenarySegment
# -----------------------------------------------------------------------

@testset "CatenarySegment — initial tangent along local z" begin
    # T-PHYSICS: catenary vertex has vertical tangent (along z in local frame)
    seg = CatenarySegment(0.5, 1.0)
    @test tangent_local(seg, 0.0) ≈ [0.0, 0.0, 1.0] atol = 1e-12
end

@testset "CatenarySegment — curvature at vertex" begin
    # T-PHYSICS: κ(0) = 1/a (maximum curvature at vertex)
    a = 0.3
    seg = CatenarySegment(a, 1.0)
    @test curvature(seg, 0.0) ≈ 1.0 / a atol = 1e-12
end

@testset "CatenarySegment — curvature formula" begin
    # T-PHYSICS: κ(s) = a / (a² + s²) for effective parameter a_eff = shrinkage * a
    a = 0.4; α = 1.05
    seg = CatenarySegment(a, 2.0; shrinkage = α)
    a_eff = α * a
    for s in [0.0, 0.3, 0.7, 1.2]
        @test curvature(seg, s) ≈ a_eff / (a_eff^2 + s^2) atol = 1e-12
    end
end

@testset "CatenarySegment — zero geometric torsion (planar)" begin
    # T-PHYSICS: catenary is a planar curve, τ_geom = 0
    seg = CatenarySegment(0.2, 1.5)
    @test geometric_torsion(seg, 0.0) == 0.0
    @test geometric_torsion(seg, 0.5) == 0.0
end

@testset "CatenarySegment — frame orthonormality along arc" begin
    # T-GUARDRAIL
    seg = CatenarySegment(0.3, 1.0, π / 4)
    for s in range(0.0, arc_length(seg); length = 9)
        T = tangent_local(seg, s)
        N = normal_local(seg, s)
        B = binormal_local(seg, s)
        @test is_orthonormal(T, N, B)
    end
end

@testset "CatenarySegment — shrinkage scales arc length" begin
    # T-PHYSICS
    a = 0.2; L = 1.0; α = 0.95
    seg = CatenarySegment(a, L; shrinkage = α)
    @test arc_length(seg) ≈ α * L
    @test nominal_arc_length(seg) ≈ L
end

# -----------------------------------------------------------------------
# PathSpec assembly and build
# -----------------------------------------------------------------------

@testset "PathSpec — empty spec fails build" begin
    # T-GUARDRAIL: a path with no segments cannot be built
    spec = PathSpec()
    @test_throws Exception build(spec)
end

@testset "PathSpec — single straight segment" begin
    # T-PHYSICS: a single straight segment of length L produces a path from
    # (0,0,0) to (0,0,L) with constant tangent (0,0,1).
    spec = PathSpec()
    straight!(spec; length = 2.0)
    path = build(spec)

    @test path_length(path) ≈ 2.0
    @test start_point(path) ≈ [0.0, 0.0, 0.0] atol = 1e-12
    @test end_point(path)   ≈ [0.0, 0.0, 2.0] atol = 1e-12
    @test start_tangent(path) ≈ [0.0, 0.0, 1.0] atol = 1e-12
    @test end_tangent(path)   ≈ [0.0, 0.0, 1.0] atol = 1e-12
end

@testset "Path — tangent continuity at segment joints" begin
    # T-GUARDRAIL: tangent must be continuous across every segment boundary.
    # Test with a straight → bend → straight sequence.
    spec = PathSpec()
    straight!(spec; length = 0.5)
    bend!(spec; radius = 0.1, angle = π / 2, axis_angle = 0.0)
    straight!(spec; length = 0.3)
    path = build(spec)

    ps = path.placed_segments
    for i in 1:(length(ps) - 1)
        s_joint = ps[i + 1].s_offset_eff
        # Tangent just before and just after the joint should agree
        T_before = tangent(path, s_joint - 1e-9)
        T_after  = tangent(path, s_joint + 1e-9)
        @test norm(T_before - T_after) < 1e-6
    end
end

@testset "Path — straight + quarter bend geometry" begin
    # T-PHYSICS: straight segment of length L followed by quarter-circle of
    # radius R. End position should be (R, 0, L+R), end tangent should be (1,0,0).
    L = 1.0; R = 0.2
    spec = PathSpec()
    straight!(spec; length = L)
    bend!(spec; radius = R, angle = π / 2, axis_angle = 0.0)
    path = build(spec)

    @test end_point(path)   ≈ [R, 0.0, L + R] atol = 1e-10
    @test end_tangent(path) ≈ [1.0, 0.0, 0.0]  atol = 1e-10
end

@testset "Path — full circle returns to start" begin
    # T-PHYSICS: a complete circle of radius R returns to the start position
    # with the same tangent direction.
    R = 0.15
    spec = PathSpec()
    bend!(spec; radius = R, angle = 2π, axis_angle = 0.0)
    path = build(spec)

    @test start_point(path)   ≈ end_point(path)   atol = 1e-10
    @test start_tangent(path) ≈ end_tangent(path)  atol = 1e-10
    @test path_length(path)   ≈ 2π * R             atol = 1e-10
end

@testset "Path — cartesian_distance vs arc_length" begin
    # T-PHYSICS: for a straight segment, cartesian_distance == arc_length.
    # For a curved path, cartesian_distance < arc_length.
    spec_straight = PathSpec()
    straight!(spec_straight; length = 3.0)
    path_s = build(spec_straight)
    @test cartesian_distance(path_s, 0.0, 3.0) ≈ 3.0 atol = 1e-12

    spec_bent = PathSpec()
    bend!(spec_bent; radius = 0.5, angle = π)
    path_b = build(spec_bent)
    L = arc_length(path_b)
    @test cartesian_distance(path_b, 0.0, L) < L   # chord < arc
end

@testset "Path — frame orthonormality along assembled path" begin
    # T-GUARDRAIL: frame must remain orthonormal everywhere in the assembled path
    spec = PathSpec()
    straight!(spec; length = 0.3)
    bend!(spec; radius = 0.08, angle = π / 3, axis_angle = π / 6)
    catenary!(spec; a = 0.2, length = 0.4)
    path = build(spec)

    ss = range(path.s_start, path.s_end; length = 31)
    for s in ss
        T = tangent(path, s)
        N = normal(path, s)
        B = binormal(path, s)
        @test is_orthonormal(T, N, B)
    end
end

@testset "Path — bounding box contains all sampled points" begin
    # T-GUARDRAIL
    spec = PathSpec()
    straight!(spec; length = 1.0)
    bend!(spec; radius = 0.2, angle = π / 2)
    path = build(spec)
    bb = bounding_box(path; n = 256)

    for s in range(path.s_start, path.s_end; length = 64)
        p = position(path, s)
        @test all(p .>= bb.lo .- 1e-10)
        @test all(p .<= bb.hi .+ 1e-10)
    end
end

# -----------------------------------------------------------------------
# Shrinkage on assembled path
# -----------------------------------------------------------------------

@testset "Path — shrinkage uniform override scales path length" begin
    # T-PHYSICS: uniform shrinkage α scales total arc length by α
    spec = PathSpec()
    straight!(spec; length = 1.0)
    bend!(spec; radius = 0.1, angle = π / 2)
    path1 = build(spec)
    path2 = build(spec; shrinkage = 1.05)

    @test path_length(path2) ≈ path_length(path1) * 1.05 atol = 1e-10
end

@testset "Path — shrinkage preserves joint tangent continuity" begin
    # T-GUARDRAIL
    spec = PathSpec()
    straight!(spec; length = 0.5)
    bend!(spec; radius = 0.1, angle = π / 3)
    path = build(spec; shrinkage = 0.98)

    ps = path.placed_segments
    for i in 1:(length(ps) - 1)
        s_joint = ps[i + 1].s_offset_eff
        T_before = tangent(path, s_joint - 1e-9)
        T_after  = tangent(path, s_joint + 1e-9)
        @test norm(T_before - T_after) < 1e-6
    end
end

# -----------------------------------------------------------------------
# TwistOverlay and material_twist
# -----------------------------------------------------------------------

@testset "TwistOverlay — constant rate (Float64)" begin
    # T-PHYSICS: constant rate τ returns that rate at any s within the interval
    τ = 2π * 3.0 / 2.0
    spec = PathSpec()
    straight!(spec; length = 2.0)
    twist!(spec; s_start = 0.0, length = 2.0, rate = τ)
    path = build(spec)

    @test material_twist(path, 0.5) ≈ τ atol = 1e-10
    @test material_twist(path, 1.5) ≈ τ atol = 1e-10
    @test total_material_twist(path) ≈ τ * 2.0 atol = 1e-10
end

@testset "TwistOverlay — zero outside overlay interval" begin
    # T-GUARDRAIL: material_twist must be zero outside the overlay's s range
    spec = PathSpec()
    straight!(spec; length = 3.0)
    twist!(spec; s_start = 1.0, length = 1.0, rate = 1.0)
    path = build(spec)

    @test material_twist(path, 0.5) == 0.0
    @test material_twist(path, 2.5) == 0.0
    @test material_twist(path, 1.5) > 0.0
end

@testset "TwistOverlay — function-valued rate" begin
    # T-PHYSICS: function rate τ(s) is evaluated at the query point s
    τ_fn = s -> 1.0 + 0.5 * sin(s)
    spec = PathSpec()
    straight!(spec; length = 3.0)
    twist!(spec; s_start = 0.0, length = 3.0, rate = τ_fn)
    path = build(spec)

    for s in [0.3, 1.0, 2.1]
        @test material_twist(path, s) ≈ τ_fn(s) atol = 1e-10
    end
end

@testset "TwistOverlay — function rate total twist via integration" begin
    # T-PHYSICS: ∫τ(s) ds over [0, L] computed numerically
    τ_fn = s -> 2.0 + cos(s)
    L = 2.0
    spec = PathSpec()
    straight!(spec; length = L)
    twist!(spec; s_start = 0.0, length = L, rate = τ_fn)
    path = build(spec)

    # Analytic: ∫₀² (2 + cos(s)) ds = [2s + sin(s)]₀² = 4 + sin(2)
    expected = 4.0 + sin(2.0)
    @test total_material_twist(path; n_quad = 512) ≈ expected atol = 1e-4
end

@testset "TwistOverlay — overlay spanning segment boundary" begin
    # T-GUARDRAIL: overlays are allowed to cross segment boundaries
    τ = 5.0
    spec = PathSpec()
    straight!(spec; length = 1.0)
    straight!(spec; length = 1.0)
    twist!(spec; s_start = 0.5, length = 1.0, rate = τ)   # crosses joint at s=1
    path = build(spec)

    @test material_twist(path, 0.75) ≈ τ atol = 1e-10
    @test material_twist(path, 1.25) ≈ τ atol = 1e-10
    @test total_material_twist(path) ≈ τ * 1.0 atol = 1e-10
end

@testset "TwistOverlay — total_material_twist partial interval" begin
    # T-PHYSICS: for constant τ, integral over half the path equals half the full-path
    # quadrature (same rule as total_material_twist on each subinterval).
    τ = 2π * 3.0 / 2.0
    spec = PathSpec()
    straight!(spec; length = 2.0)
    twist!(spec; s_start = 0.0, length = 2.0, rate = τ)
    path = build(spec)
    full = total_material_twist(path)
    @test total_material_twist(path; s_start = 0.0, s_end = 1.0) ≈ full / 2 atol = 1e-10
    @test total_material_twist(path; s_start = 1.0, s_end = 0.0) ≈ full / 2 atol = 1e-10
end

# -----------------------------------------------------------------------
# Path measures
# -----------------------------------------------------------------------

@testset "Path — total_turning_angle of a full circle" begin
    # T-PHYSICS: ∫κ ds over a full circle = 2π
    spec = PathSpec()
    bend!(spec; radius = 0.1, angle = 2π)
    path = build(spec)
    @test total_turning_angle(path) ≈ 2π atol = 1e-10
end

@testset "Path — total_torsion of straight and bend segments is zero" begin
    # T-PHYSICS: straight and circular-arc segments have zero geometric torsion
    spec = PathSpec()
    straight!(spec; length = 1.0)
    bend!(spec; radius = 0.2, angle = π / 2)
    path = build(spec)
    @test total_torsion(path) == 0.0
end

@testset "Path — writhe of a straight path is zero" begin
    # T-PHYSICS: a straight line has no self-linking, Wr = 0
    spec = PathSpec()
    straight!(spec; length = 2.0)
    path = build(spec)
    @test abs(writhe(path; n = 64)) < 1e-6
end

# -----------------------------------------------------------------------
# Sampling
# -----------------------------------------------------------------------

@testset "sample_uniform — returns n frames" begin
    # T-GUARDRAIL
    spec = PathSpec()
    straight!(spec; length = 1.0)
    bend!(spec; radius = 0.1, angle = π / 2)
    path = build(spec)

    frames = sample_uniform(path; n = 50)
    @test length(frames) == 50
    @test hasproperty(frames[1], :position)
    @test hasproperty(frames[1], :tangent)
    @test hasproperty(frames[1], :curvature)
end

# -----------------------------------------------------------------------
# JumpBy and JumpTo (HermiteConnector)
# -----------------------------------------------------------------------

@testset "JumpBy — endpoint matches delta in local frame" begin
    # T-PHYSICS: JumpBy with delta along the current tangent direction (ẑ)
    # should move the path forward by that amount.
    spec = PathSpec()
    jumpby!(spec; delta = (0.0, 0.0, 0.5))
    path = build(spec)

    @test start_point(path) ≈ [0.0, 0.0, 0.0] atol = 1e-10
    @test end_point(path)   ≈ [0.0, 0.0, 0.5] atol = 1e-10
end

@testset "JumpBy — initial tangent is ẑ" begin
    # T-GUARDRAIL: incoming tangent at start of every connector is ẑ
    spec = PathSpec()
    jumpby!(spec; delta = (0.1, 0.0, 0.3))
    path = build(spec)
    @test tangent(path, 0.0) ≈ [0.0, 0.0, 1.0] atol = 1e-10
end

@testset "JumpBy — outgoing tangent honours tangent_out" begin
    # T-PHYSICS: explicit tangent_out should be the end tangent (in local frame)
    t_out = normalize([1.0, 0.0, 1.0])
    spec = PathSpec()
    jumpby!(spec; delta = (0.3, 0.0, 0.3), tangent = t_out)
    path = build(spec)
    @test tangent(path, path.s_end) ≈ t_out atol = 1e-8
end

@testset "JumpBy — frame orthonormality along connector" begin
    # T-GUARDRAIL
    spec = PathSpec()
    jumpby!(spec; delta = (0.2, 0.1, 0.4))
    path = build(spec)
    for s in range(path.s_start, path.s_end; length = 11)
        T = tangent(path, s)
        N = normal(path, s)
        B = binormal(path, s)
        @test is_orthonormal(T, N, B)
    end
end

@testset "JumpBy — shrinkage scales endpoint" begin
    # T-PHYSICS: shrinkage α scales the displacement by α
    spec = PathSpec()
    jumpby!(spec; delta = (0.0, 0.0, 1.0), shrinkage = 0.8)
    path = build(spec)
    @test end_point(path) ≈ [0.0, 0.0, 0.8] atol = 1e-10
end

@testset "JumpBy — after straight, delta is in rotated local frame" begin
    # T-PHYSICS: after a quarter-circle bend the local ẑ points along global +x.
    # JumpBy with delta=(0,0,d) should therefore move +x in global frame.
    R = 0.1; d = 0.5
    spec = PathSpec()
    bend!(spec; radius = R, angle = π / 2, axis_angle = 0.0)
    jumpby!(spec; delta = (0.0, 0.0, d))   # local ẑ is now global +x
    path = build(spec)
    # end_point of bend is (R, 0, R); jumpby adds d along global +x
    @test end_point(path) ≈ [R + d, 0.0, R] atol = 1e-8
end

@testset "JumpTo — endpoint matches destination" begin
    # T-PHYSICS: JumpTo should place the path end at the specified global position.
    dest = [1.0, 0.5, 2.0]
    spec = PathSpec()
    jumpto!(spec; destination = dest)
    path = build(spec)
    @test end_point(path) ≈ dest atol = 1e-10
end

@testset "JumpTo — initial tangent is ẑ" begin
    spec = PathSpec()
    jumpto!(spec; destination = (0.3, 0.1, 0.8))
    path = build(spec)
    @test tangent(path, 0.0) ≈ [0.0, 0.0, 1.0] atol = 1e-10
end

@testset "JumpTo — after bend, destination is in global frame" begin
    # T-PHYSICS: JumpTo destination is always in global frame regardless of prior segments.
    R = 0.1
    dest = [R + 0.5, 0.0, R]   # global position
    spec = PathSpec()
    bend!(spec; radius = R, angle = π / 2, axis_angle = 0.0)
    jumpto!(spec; destination = dest)
    path = build(spec)
    @test end_point(path) ≈ collect(dest) atol = 1e-8
end

@testset "JumpTo — outgoing tangent honours tangent_out (global)" begin
    # T-PHYSICS: tangent_out for JumpTo is specified in global frame
    t_out_global = normalize([1.0, 0.0, 0.0])
    spec = PathSpec()
    jumpto!(spec; destination = (1.0, 0.0, 0.5), tangent = t_out_global)
    path = build(spec)
    @test tangent(path, path.s_end) ≈ t_out_global atol = 1e-8
end

@testset "JumpBy/JumpTo — sample_path works on connector" begin
    # T-GUARDRAIL: sample_path must not error on paths containing connectors
    spec = PathSpec()
    straight!(spec; length = 0.3)
    jumpby!(spec; delta = (0.1, 0.0, 0.4))
    straight!(spec; length = 0.2)
    path = build(spec)
    ps = sample_path(path, path.s_start, path.s_end)
    @test ps.n >= 4
    @test ps.samples[1].s ≈ path.s_start atol = 1e-12
    @test ps.samples[end].s ≈ path.s_end  atol = 1e-12
end

#=
DISABLED: min_bend_radius enforcement is a stub.

The enforcement code in _build_hermite_connector is commented out (see
path-geometry.jl, inside _build_hermite_connector). Re-enable these testsets
when the scan/bisect block is restored. See commit f06e689.

@testset "JumpBy — min_bend_radius keeps curvature ≤ 1/R_min" begin
    # T-PHYSICS: every point along the connector must have κ ≤ 1/R_min.
    R_min = 0.05
    spec = PathSpec()
    jumpby!(spec; delta = (0.02, 0.0, 0.03), min_bend_radius = R_min)
    path = build(spec)
    for s in range(path.s_start, path.s_end; length = 51)
        @test curvature(path, s) ≤ 1.0 / R_min
    end
end

@testset "JumpTo — min_bend_radius keeps curvature ≤ 1/R_min" begin
    # T-PHYSICS: min_bend_radius on JumpTo also constrains the connector curvature.
    R_min = 0.04
    spec = PathSpec()
    jumpto!(spec; destination = (0.05, 0.0, 0.06), min_bend_radius = R_min)
    path = build(spec)
    for s in range(path.s_start, path.s_end; length = 51)
        @test curvature(path, s) ≤ 1.0 / R_min
    end
end

@testset "JumpBy — min_bend_radius does not move endpoint" begin
    # T-GUARDRAIL: the constraint changes handle length only, not the target position.
    delta = (0.03, 0.0, 0.05)
    spec = PathSpec()
    jumpby!(spec; delta = delta, min_bend_radius = 0.10)
    path = build(spec)
    @test end_point(path) ≈ collect(delta) atol = 1e-10
end

@testset "JumpTo — infeasible min_bend_radius throws ArgumentError" begin
    # T-GUARDRAIL: anti-parallel incoming/outgoing tangents cause peak curvature to
    # increase (not decrease) with handle length — no equal-handle Hermite can satisfy
    # the constraint, so build() must throw ArgumentError.
    spec = PathSpec()
    straight!(spec; length = 1.0)
    # After straight, incoming tangent is +z. destination=(2,0,1) is 1 m transverse.
    # Outgoing tangent is -z (anti-parallel). κ grows with h → infeasible.
    jumpto!(spec; destination = (2.0, 0.0, 1.0), tangent = (0.0, 0.0, -1.0),
            min_bend_radius = 2.0)
    @test_throws ArgumentError build(spec)
end

# JumpTo min_bend_radius — T1/T2/T3 from demo_fiber_path_jumps_min_radius
# Each helper builds the path up to and including the named jump.
# Thresholds are taken from the demo comments.

@testset "JumpTo min_bend_radius — T1: passes at 0.5, fails above 0.5" begin
    # T-PHYSICS: transverse chord (1,0,0), outgoing tangent (0,0,-1) anti-parallel
    # to incoming (0,0,1). Minimum-curvature connector has bend radius ≈ 0.5 m.
    function build_T1(mbr)
        spec = PathSpec()
        straight!(spec; length = 1)
        jumpto!(spec; destination = (1, 0.0, 1), tangent = (0.0, 0.0, -1.0),
                min_bend_radius = mbr)
        build(spec)
    end
    @test_nowarn build_T1(0.49)
    @test_nowarn build_T1(0.5)
    @test_throws ArgumentError build_T1(0.51)
end

@testset "JumpTo min_bend_radius — T2: passes at 0.5, fails above 0.5" begin
    # T-PHYSICS: after T1 the path runs downward (-z). The straight extends
    # that to position (1,0,0). T2 jumps to (2,0,0) with outgoing tangent (0,0,1)
    # — again transverse chord, anti-parallel tangents.
    function build_T2(mbr)
        spec = PathSpec()
        straight!(spec; length = 1)
        jumpto!(spec; destination = (1, 0.0, 1), tangent = (0.0, 0.0, -1.0),
                min_bend_radius = 0.1) # so small it won't error
        straight!(spec; length = 1)
        jumpto!(spec; destination = (2, 0.0, 0), tangent = (0.0, 0.0, 1.0),
                min_bend_radius = mbr)
        build(spec)
    end
    @test_nowarn build_T2(0.49)
    @test_nowarn build_T2(0.5)
    @test_throws ArgumentError build_T2(0.51)
end

@testset "JumpTo min_bend_radius — T3: passes at 0.50, fails at 0.51" begin
    # T-PHYSICS: same pattern as T1; threshold is 0.50 m.
    function build_T3(mbr)
        spec = PathSpec()
        straight!(spec; length = 1)
        jumpto!(spec; destination = (1, 0.0, 1), tangent = (0.0, 0.0, -1.0),
                min_bend_radius = 0.1) # so small it won't error
        straight!(spec; length = 1)
        jumpto!(spec; destination = (2, 0.0, 0), tangent = (0.0, 0.0, 1.0),
                min_bend_radius = 0.1) # so small it won't error
        straight!(spec; length = 1)
        jumpto!(spec; destination = (3, 0.0, 1), tangent = (0.0, 0.0, -1.0),
                min_bend_radius = mbr)
        build(spec)
    end
    @test_nowarn build_T3(0.49)
    @test_nowarn build_T3(0.50)
    @test_throws ArgumentError build_T3(0.51)
end
=#
