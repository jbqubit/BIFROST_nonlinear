using Test
using LinearAlgebra

if !isdefined(Main, :sample_fiber_centerline)
    include("../fiber-path.jl")
    include("../path-integral.jl")
    include("../fiber-path-plot.jl")
end

# Seal a SubpathBuilder at the natural exit of its current interior segments
# via a trial-build. Used by tests that author paths whose natural exit is
# not analytically convenient.
function _test_seal_natural!(sb::SubpathBuilder)
    @assert isnothing(sb.jumpto_point) "_test_seal_natural!: builder already sealed"
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

const CENTERLINE_ATOL = 5e-8
const TANGENT_ATOL = 5e-13
const SAMPLE_COUNT = 8001

function test_cross_section()
    return FiberCrossSection(
        GermaniaSilicaGlass(0.036),
        GermaniaSilicaGlass(0.0),
        8.2e-6,
        125e-6
    )
end

function canonical_bend_segment(seg)
    if length(seg) == 2
        R, θ = seg
        α = 0.0
    else
        R, α, θ = seg
    end

    @assert R > 0.0 "Bend radius must be positive"
    αeff = θ >= 0 ? α : α + π
    return (R = Float64(R), α = Float64(αeff), θ = abs(Float64(θ)))
end

function skew_matrix(u::AbstractVector{<:Real})
    return Float64[
         0.0   -u[3]   u[2]
         u[3]   0.0   -u[1]
        -u[2]   u[1]   0.0
    ]
end

function rotate_vector_about_axis(v::Vector{Float64}, u::Vector{Float64}, θ::Float64)
    c = cos(θ)
    s = sin(θ)
    return c .* v + s .* cross(u, v) + (1 - c) .* dot(u, v) .* u
end

function exact_endpoint_and_tangent(segments)
    r = zeros(3)
    N = [1.0, 0.0, 0.0]
    B = [0.0, 1.0, 0.0]
    T = [0.0, 0.0, 1.0]

    for raw_seg in segments
        seg = canonical_bend_segment(raw_seg)
        if seg.θ == 0.0
            continue
        end

        n_local = [cos(seg.α), sin(seg.α), 0.0]
        local_frame = hcat(N, B, T)
        local_step = seg.R * (1 - cos(seg.θ)) * n_local + [0.0, 0.0, seg.R * sin(seg.θ)]
        r .+= local_frame * local_step

        T_local = [sin(seg.θ) * cos(seg.α), sin(seg.θ) * sin(seg.α), cos(seg.θ)]
        N_local = [cos(seg.θ) * cos(seg.α), cos(seg.θ) * sin(seg.α), -sin(seg.θ)]
        B_local = [-sin(seg.α), cos(seg.α), 0.0]

        T = local_frame * T_local
        N = local_frame * N_local
        B = local_frame * B_local
        T ./= norm(T)
        N ./= norm(N)
        B ./= norm(B)
    end

    return r, T
end

function build_centerline_test_fiber(segments)
    @assert !isempty(segments) "Use build_straight_test_fiber for the straight case"
    sb = SubpathBuilder(); start!(sb)
    for raw_seg in segments
        seg = canonical_bend_segment(raw_seg)
        bend!(sb; radius = seg.R, angle = seg.θ, axis_angle = seg.α)
    end
    _test_seal_natural!(sb)
    path = build(Subpath(sb))
    return Fiber(path; cross_section = test_cross_section())
end

function build_straight_test_fiber(length::Real)
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = length)
    _test_seal_natural!(sb)
    path = build(Subpath(sb))
    return Fiber(path; cross_section = test_cross_section())
end

function sample_final_centerline_state(fiber::Fiber; n::Int = SAMPLE_COUNT)
    path = sample_fiber_centerline(fiber, fiber.s_start, fiber.s_end; n = n)
    endpoint = [path.x[end], path.y[end], path.zc[end]]
    tangent = [path.tx[end], path.ty[end], path.tz[end]]
    return path, endpoint, tangent
end

const EXPLICIT_PATH_CASES = [
    (
        name = "single_quarter_circle_xz_plane",
        segments = [(1.0, 0.0, π / 2)],
        xyz = [1.0, 0.0, 1.0],
        tangent = [1.0, 0.0, 0.0],
        planar = true
    ),
    (
        name = "single_half_circle_xz_plane",
        segments = [(2.0, 0.0, π)],
        xyz = [4.0, 0.0, 0.0],
        tangent = [0.0, 0.0, -1.0],
        planar = true
    ),
    (
        name = "two_quarters_make_half_circle",
        segments = [(1.0, 0.0, π / 2), (1.0, 0.0, π / 2)],
        xyz = [2.0, 0.0, 0.0],
        tangent = [0.0, 0.0, -1.0],
        planar = true
    ),
    (
        name = "three_sixty_degree_bends_make_half_circle",
        segments = [(1.0, 0.0, π / 3), (1.0, 0.0, π / 3), (1.0, 0.0, π / 3)],
        xyz = [2.0, 0.0, 0.0],
        tangent = [0.0, 0.0, -1.0],
        planar = true
    ),
    (
        name = "single_quarter_circle_yz_plane",
        segments = [(1.0, π / 2, π / 2)],
        xyz = [0.0, 1.0, 1.0],
        tangent = [0.0, 1.0, 0.0],
        planar = false
    ),
    (
        name = "single_oblique_quarter_circle",
        segments = [(2.0, π / 4, π / 2)],
        xyz = [sqrt(2), sqrt(2), 2.0],
        tangent = [inv(sqrt(2)), inv(sqrt(2)), 0.0],
        planar = false
    ),
    (
        name = "orthogonal_quarters_lab_frame",
        segments = [(1.0, 0.0, π / 2), (1.0, π / 2, π / 2)],
        xyz = [2.0, 1.0, 1.0],
        tangent = [0.0, 1.0, 0.0],
        planar = false
    ),
    (
        name = "orthogonal_quarters_swapped_order",
        segments = [(1.0, π / 2, π / 2), (1.0, 0.0, π / 2)],
        xyz = [0.0, 2.0, 0.0],
        tangent = [0.0, 0.0, -1.0],
        planar = false
    )
]

@testset "fiber-path centerline geometry" begin
    @testset "Straight section" begin
        fiber = build_straight_test_fiber(2.5)
        path, endpoint, tangent = sample_final_centerline_state(fiber)

        @test endpoint ≈ [0.0, 0.0, 2.5] atol = CENTERLINE_ATOL
        @test tangent ≈ [0.0, 0.0, 1.0] atol = TANGENT_ATOL
        @test maximum(abs.(path.x)) <= CENTERLINE_ATOL
        @test maximum(abs.(path.y)) <= CENTERLINE_ATOL
    end

    @testset "Explicit endpoint cases" begin
        for case in EXPLICIT_PATH_CASES
            @testset "$(case.name)" begin
                exact_xyz, exact_tangent = exact_endpoint_and_tangent(case.segments)
                fiber = build_centerline_test_fiber(case.segments)
                path, sampled_xyz, sampled_tangent = sample_final_centerline_state(fiber)

                @test exact_xyz ≈ case.xyz atol = 1e-12
                @test exact_tangent ≈ case.tangent atol = 1e-12
                @test sampled_xyz ≈ exact_xyz atol = CENTERLINE_ATOL
                @test sampled_tangent ≈ exact_tangent atol = TANGENT_ATOL
                @test norm(sampled_tangent) ≈ 1.0 atol = TANGENT_ATOL
                @test path.s[end] ≈ fiber.s_end atol = eps(Float64)

                if case.planar
                    @test maximum(abs.(path.y)) <= CENTERLINE_ATOL
                end
            end
        end
    end

    @testset "Segmentation invariance" begin
        reference = [(1.0, 0.0, π / 2)]
        split_two = [(1.0, 0.0, π / 4), (1.0, 0.0, π / 4)]
        split_four = [(1.0, 0.0, π / 8), (1.0, 0.0, π / 8), (1.0, 0.0, π / 8), (1.0, 0.0, π / 8)]

        ref_xyz_exact, ref_t_exact = exact_endpoint_and_tangent(reference)
        ref_path, ref_xyz, ref_t = sample_final_centerline_state(build_centerline_test_fiber(reference))
        two_path, two_xyz, two_t = sample_final_centerline_state(build_centerline_test_fiber(split_two))
        four_path, four_xyz, four_t = sample_final_centerline_state(build_centerline_test_fiber(split_four))

        @test ref_xyz_exact ≈ [1.0, 0.0, 1.0] atol = 1e-12
        @test ref_t_exact ≈ [1.0, 0.0, 0.0] atol = 1e-12

        @test ref_xyz ≈ ref_xyz_exact atol = CENTERLINE_ATOL
        @test two_xyz ≈ ref_xyz_exact atol = CENTERLINE_ATOL
        @test four_xyz ≈ ref_xyz_exact atol = CENTERLINE_ATOL
        @test ref_t ≈ ref_t_exact atol = TANGENT_ATOL
        @test two_t ≈ ref_t_exact atol = TANGENT_ATOL
        @test four_t ≈ ref_t_exact atol = TANGENT_ATOL

        @test ref_xyz ≈ two_xyz atol = CENTERLINE_ATOL
        @test ref_xyz ≈ four_xyz atol = CENTERLINE_ATOL
        @test two_xyz ≈ four_xyz atol = CENTERLINE_ATOL
        @test maximum(abs.(ref_path.y)) <= CENTERLINE_ATOL
        @test maximum(abs.(two_path.y)) <= CENTERLINE_ATOL
        @test maximum(abs.(four_path.y)) <= CENTERLINE_ATOL
    end

    @testset "Hard exact oracle case" begin
        segments = [(1.0, π / 4, π / 3), (1.0, -π / 4, π / 3)]
        exact_xyz, exact_tangent = exact_endpoint_and_tangent(segments)
        fiber = build_centerline_test_fiber(segments)
        _, sampled_xyz, sampled_tangent = sample_final_centerline_state(fiber)

        @test sampled_xyz ≈ exact_xyz atol = CENTERLINE_ATOL
        @test sampled_tangent ≈ exact_tangent atol = TANGENT_ATOL
    end

    @testset "Order matters in 3D" begin
        case_a = [(1.0, 0.0, π / 2), (1.0, π / 2, π / 2)]
        case_b = [(1.0, π / 2, π / 2), (1.0, 0.0, π / 2)]

        _, xyz_a, _ = sample_final_centerline_state(build_centerline_test_fiber(case_a))
        _, xyz_b, _ = sample_final_centerline_state(build_centerline_test_fiber(case_b))

        @test norm(xyz_a - xyz_b) > 1.0
    end
end

@testset "Fiber reference temperature" begin
    xs = test_cross_section()

    function straight_fiber(; T_ref_K = DEFAULT_T_REF_K, length = 10.0)
        sb = SubpathBuilder(); start!(sb)
        straight!(sb; length = length)
        _test_seal_natural!(sb)
        return Fiber(build(Subpath(sb)); cross_section = xs, T_ref_K = T_ref_K)
    end

    @testset "default T_ref_K is 297.15" begin
        fiber = straight_fiber()
        @test fiber.T_ref_K ≈ 297.15
    end

    @testset "T_ref_K round-trips" begin
        fiber = straight_fiber(T_ref_K = 310.0)
        @test fiber.T_ref_K == 310.0
    end
end

@testset "Path-backed fiber assembly" begin
    xs = test_cross_section()

    # TODO: spinning refactor — the spinning! call and spinning breakpoints in this test
    # are pending the per-segment-meta spinning subsystem.
    @testset "T-GUARDRAIL: domain, coverage, and breakpoint derivation" begin
        sb = SubpathBuilder(); start!(sb)
        straight!(sb; length = 1.0)
        bend!(sb; radius = 0.2, angle = π / 2)
        _test_seal_natural!(sb)
        path = build(Subpath(sb))
        fiber = Fiber(path; cross_section = xs)

        @test fiber.s_start == 0.0
        @test fiber.s_end ≈ Float64(_qc_nominalize(arc_length(path)))
        @test fiber.path === path
        @test fiber.cross_section === xs

        # Breakpoints: 0, end-of-first-straight, end-of-bend, plus the
        # terminal connector start (which equals end-of-bend in the
        # natural-exit seal).
        bps = fiber_breakpoints(fiber)
        @test 0.0 in bps
        @test any(b -> isapprox(b, 1.0; atol = 1e-9), bps)
        @test any(b -> isapprox(b, 1.0 + 0.2 * (π / 2); atol = 1e-9), bps)
        @test breakpoints(fiber.path) == bps
    end

    @testset "T-PHYSICS: straight path gives zero generators" begin
        sb = SubpathBuilder(); start!(sb)
        straight!(sb; length = 0.8)
        _test_seal_natural!(sb)
        fiber = Fiber(build(Subpath(sb)); cross_section = xs)

        @test generator_K(fiber, 1550e-9)(0.4) ≈ zeros(ComplexF64, 2, 2) atol = 1e-14
        @test generator_Kω(fiber, 1550e-9)(0.4) ≈ zeros(ComplexF64, 2, 2) atol = 1e-14
    end

    @testset "T-PHYSICS: circular bend uses bending_birefringence" begin
        λ = 1550e-9
        T = 297.15
        R = 0.04
        sb = SubpathBuilder(); start!(sb)
        bend!(sb; radius = R, angle = π / 3)
        _test_seal_natural!(sb)
        fiber = Fiber(build(Subpath(sb)); cross_section = xs, T_ref_K = T)
        # Query at the midpoint of the bend (which is interior, not the
        # connector). The bend's arc length is R*π/3.
        K = generator_K(fiber, λ)(0.5 * R * π / 3)
        Δβ = bending_birefringence(xs, λ, T; bend_radius_m = R)

        @test K[1, 1] ≈ 0.5im * Δβ atol = 1e-12
        @test K[2, 2] ≈ -0.5im * Δβ atol = 1e-12
        @test K[1, 2] ≈ 0.0 atol = 1e-12
        @test K[2, 1] ≈ 0.0 atol = 1e-12
    end

    # TODO: spinning refactor — pending per-segment-meta spinning subsystem.
    @testset "T-PHYSICS: spinning overlay uses twisting_birefringence" begin
        @test_skip true
    end

end
