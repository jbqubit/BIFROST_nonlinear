using Test
using LinearAlgebra

include("../path-demo.jl")

@testset "Fiber source API" begin
    bend = BendSource(
        s -> 2.0,
        s -> 0.0,
        k2 -> 3.0,
        k2 -> 5.0;
        coverage = [(0.0, 1.0)],
        breakpoints = [0.0, 0.4, 1.0]
    )
    twist = TwistSource(
        s -> 7.0,
        tau -> 11.0,
        tau -> 13.0;
        coverage = [(0.0, 1.0)],
        breakpoints = [0.0, 0.6, 1.0]
    )

    fiber = Fiber(0.0, 1.0, AbstractBirefringenceSource[bend, twist])

    @test fiber_breakpoints(fiber) == [0.0, 0.4, 0.6, 1.0]

    s = 0.25
    K_expected = generator_contribution(bend, s) + generator_contribution(twist, s)
    Kω_expected = generator_omega_contribution(bend, s) + generator_omega_contribution(twist, s)
    @test make_generator(fiber)(s) ≈ K_expected
    @test make_generator_omega(fiber)(s) ≈ Kω_expected

    bad_twist = TwistSource(
        s -> 0.0,
        tau -> tau,
        tau -> 0.0;
        coverage = [(0.0, 0.4), (0.5, 1.0)],
        breakpoints = [0.0, 0.4, 0.5, 1.0]
    )
    @test_throws Exception Fiber(0.0, 1.0, AbstractBirefringenceSource[bend, bad_twist])

    zero_bend = BendSource(
        s -> Inf,
        s -> 0.0,
        k2 -> 1.0,
        k2 -> 0.0;
        coverage = [(0.0, 1.0)],
        breakpoints = [0.0, 1.0]
    )
    zero_twist = TwistSource(
        s -> 0.0,
        tau -> tau,
        tau -> 0.0;
        coverage = [(0.0, 1.0)],
        breakpoints = [0.0, 1.0]
    )
    zero_fiber = Fiber(0.0, 1.0, AbstractBirefringenceSource[zero_bend, zero_twist])
    @test make_generator(zero_fiber)(0.5) ≈ zeros(ComplexF64, 2, 2)

    J_zero, G_zero, _ = propagate_fiber_sensitivity(zero_fiber; h_init = 0.1)
    @test output_dgd(J_zero, G_zero) ≈ 0.0

    custom = Fiber(0.0, 2.0, AbstractBirefringenceSource[
        BendSource(s -> Inf, s -> 0.0, k2 -> 0.0, k2 -> 0.0; coverage = [(0.0, 2.0)], breakpoints = [0.0, 2.0]),
        TwistSource(s -> 0.0, tau -> 0.0, tau -> 0.0; coverage = [(0.0, 2.0)], breakpoints = [0.0, 2.0])
    ])
    K = s -> zeros(ComplexF64, 2, 2)
    Kω = s -> 0.5im .* ComplexF64[1 0; 0 -1]
    J_test, G_test, _ = propagate_piecewise_sensitivity(K, Kω, fiber_breakpoints(custom); h_init = 0.1, rtol = 1e-11, atol = 1e-13)
    @test output_dgd(J_test, G_test) ≈ 2.0 atol = 1e-8

    demo = demofiber1()
    J_demo, stats_demo = propagate_fiber(demo.fiber; rtol = 1e-10, atol = 1e-12, h_init = 1e-2)
    @test size(J_demo) == (2, 2)
    @test !isempty(stats_demo)

    J_demo_sens, G_demo_sens, stats_demo_sens = propagate_fiber_sensitivity(demo.fiber; rtol = 1e-10, atol = 1e-12, h_init = 1e-2)
    @test size(J_demo_sens) == (2, 2)
    @test size(G_demo_sens) == (2, 2)
    @test !isempty(stats_demo_sens)
    @test isfinite(output_dgd(J_demo_sens, G_demo_sens))

    plot_path = write_fiber_input_plot3d(
        demo.fiber,
        demo.fiber.s_start,
        demo.fiber.s_end;
        n = 25,
        output = "/tmp/test_path_integral_sources_plot.html"
    )
    @test isfile(plot_path)
end
