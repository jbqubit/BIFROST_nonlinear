if !isdefined(Main, :propagate_fiber)
    include("path-integral.jl")
end

const DEMO_FIBER_CROSS_SECTION = FiberCrossSection(
    GermaniaSilicaGlass(0.036),
    GermaniaSilicaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like demo"
)

"""
    demofiber1()
The shape of the fiber in this file is that of polarization paddles.
The first loop is L/4, the second is L/2, and the third is L/4, where 
L is the length of fiber in each loop. The angles of the loops are 0, π/2, 
and 0, respectively. The input state is vertical polarization.
"""
function demofiber1()
    lead_in = 4.0
    spacer = 8.5
    lead_out = 4.0
    loop_turns = [1.0, 1.0, 1.0]
    target_retardance = [π / 2, π, π / 2]
    loop_radii = 2π .* loop_turns ./ target_retardance
    loop_angles = [0.0, π / 2, 0.0]

    total_length = lead_in + sum(2π .* loop_radii .* loop_turns) + 2 * spacer + lead_out
    spec = FiberSpec(0.0, total_length; cross_section = DEMO_FIBER_CROSS_SECTION)

    s = 0.0
    bend!(spec, s, s + lead_in; angle = 0.0, axis = loop_angles[1])
    s += lead_in

    loop1_angle = 2π * loop_turns[1]
    loop1_length = loop_radii[1] * loop1_angle
    bend!(spec, s, s + loop1_length; angle = loop1_angle, axis = loop_angles[1])
    s += loop1_length

    bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[1])
    twist!(spec, s, s + spacer; rate = (loop_angles[2] - loop_angles[1]) / spacer)
    s += spacer

    loop2_angle = 2π * loop_turns[2]
    loop2_length = loop_radii[2] * loop2_angle
    bend!(spec, s, s + loop2_length; angle = loop2_angle, axis = loop_angles[2])
    s += loop2_length

    bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[2])
    twist!(spec, s, s + spacer; rate = (loop_angles[3] - loop_angles[2]) / spacer)
    s += spacer

    loop3_angle = 2π * loop_turns[3]
    loop3_length = loop_radii[3] * loop3_angle
    bend!(spec, s, s + loop3_length; angle = loop3_angle, axis = loop_angles[3])
    s += loop3_length

    bend!(spec, s, s + lead_out; angle = 0.0, axis = loop_angles[3])

    return (
        fiber = build(spec),
        title = "Fiber polarization control paddles",
        input_state = ComplexF64[0.0 + 0.0im, 1.0 + 0.0im],
        output = "fiber-path-plot.html",
        n = 4001
    )
end


"""
    demofiber2()
This fiber is a more complex version of the first demo, with four loops instead of three, 
and with the loop angles changing by π/3 between each loop. The retardances of the loops 
are π/2, π, π/2, and π/3. The input state is diagonal polarization.  
"""
function demofiber2()
    lead_in = 3.0
    spacer = 6.0
    lead_out = 3.0
    loop_turns = [0.75, 1.0, 0.75, 0.5]
    target_retardance = [π / 2, π, π / 2, π / 3]
    loop_radii = 2π .* loop_turns ./ target_retardance
    loop_angles = [0.0, π / 3, 2π / 3, π / 6]

    total_length = lead_in + sum(2π .* loop_radii .* loop_turns) + 3 * spacer + lead_out
    spec = FiberSpec(0.0, total_length; cross_section = DEMO_FIBER_CROSS_SECTION)

    s = 0.0
    bend!(spec, s, s + lead_in; angle = 0.0, axis = loop_angles[1])
    s += lead_in

    loop1_angle = 2π * loop_turns[1]
    loop1_length = loop_radii[1] * loop1_angle
    bend!(spec, s, s + loop1_length; angle = loop1_angle, axis = loop_angles[1])
    s += loop1_length

    bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[1])
    twist!(spec, s, s + spacer; rate = (loop_angles[2] - loop_angles[1]) / spacer)
    s += spacer

    loop2_angle = 2π * loop_turns[2]
    loop2_length = loop_radii[2] * loop2_angle
    bend!(spec, s, s + loop2_length; angle = loop2_angle, axis = loop_angles[2])
    s += loop2_length

    bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[2])
    twist!(spec, s, s + spacer; rate = (loop_angles[3] - loop_angles[2]) / spacer)
    s += spacer

    loop3_angle = 2π * loop_turns[3]
    loop3_length = loop_radii[3] * loop3_angle
    bend!(spec, s, s + loop3_length; angle = loop3_angle, axis = loop_angles[3])
    s += loop3_length

    bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[3])
    twist!(spec, s, s + spacer; rate = (loop_angles[4] - loop_angles[3]) / spacer)
    s += spacer

    loop4_angle = 2π * loop_turns[4]
    loop4_length = loop_radii[4] * loop4_angle
    bend!(spec, s, s + loop4_length; angle = loop4_angle, axis = loop_angles[4])
    s += loop4_length

    bend!(spec, s, s + lead_out; angle = 0.0, axis = loop_angles[4])

    return (
        fiber = build(spec),
        title = "Four-loop twisted fiber controller",
        input_state = ComplexF64[inv(sqrt(2)) + 0.0im, 0.0 + inv(sqrt(2)) * im],
        output = "fiber-path-plot.html",
        n = 4001
    )
end

function run_demo(demo = demofiber1(); output::Union{Nothing,AbstractString} = nothing)
    fiber = demo.fiber
    title = demo.title
    input_state = demo.input_state
    plot_output = isnothing(output) ? demo.output : output
    n = demo.n

    J_final, stats = propagate_fiber(
        fiber;
        rtol = 1e-10,
        atol = 1e-12,
        h_init = 1e-2
    )

    println("Final Jones matrix:")
    println(J_final)

    println("\nInterval stats:")
    for (i, st) in enumerate(stats)
        println("interval $i: accepted=$(st.accepted_steps), rejected=$(st.rejected_steps)")
    end

    plot_path = write_fiber_input_plot3d(
        fiber,
        fiber.s_start,
        fiber.s_end;
        n = n,
        output = plot_output,
        title = title,
        input_state = input_state
    )

    println("\nWrote interactive 3D plot to:")
    println(plot_path)

    return (; demo..., J_final, stats, plot_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_demo()
end
