include("path-integral.jl")

"""
# Execution Flow #
path-demo.jl includes path-integral.jl, which includes path-plot.jl.

run_demo() gets a demo bundle from demofiber1() or demofiber2().

It computes the final Jones matrix with 
    make_generator -> propagate_piecewise -> 
    propagate_interval! -> exp_midpoint_step -> exp_jones_generator.

It then plans the visualization with write_fiber_input_plot3d() via subroutines:
- geometry via sample_fiber_centerline()
- polarization via sample_polarization_trajectory()
It then assembles the main 3D fiber plot consisting of
- moving polarization-circle cursor overlay
- Poincare sphere inset
The output is an HTML file

"""


function build_piecewise_demo_fiber(
    lead_in,
    spacer,
    lead_out,
    loop_turns,
    loop_radii,
    loop_angles;
    title::AbstractString,
    input_state::Vector{ComplexF64} = ComplexF64[0.0 + 0.0im, 1.0 + 0.0im],
    output::AbstractString = "fiberinput_demo.html",
    n::Int = 4001
)
    @assert length(loop_turns) == length(loop_radii) == length(loop_angles)

    breaks = Float64[0.0]
    pieces_Rb = Function[]
    pieces_theta = Function[]
    pieces_dtwist = Function[]

    function push_constant_segment!(Δs, Rval, θval, τval)
        push!(pieces_Rb, _ -> Rval)
        push!(pieces_theta, _ -> θval)
        push!(pieces_dtwist, _ -> τval)
        push!(breaks, breaks[end] + Δs)
    end

    function push_twist_segment!(Δs, θ0, θ1)
        sL = breaks[end]
        dθds = (θ1 - θ0) / Δs
        push!(pieces_Rb, _ -> Inf)
        push!(pieces_theta, s -> θ0 + dθds * (s - sL))
        push!(pieces_dtwist, _ -> dθds)
        push!(breaks, sL + Δs)
    end

    push_constant_segment!(lead_in, Inf, loop_angles[1], 0.0)
    for i in eachindex(loop_turns)
        push_constant_segment!(loop_turns[i] * 2π * loop_radii[i], loop_radii[i], loop_angles[i], 0.0)
        if i < length(loop_turns)
            push_twist_segment!(spacer, loop_angles[i], loop_angles[i + 1])
        end
    end
    push_constant_segment!(lead_out, Inf, loop_angles[end], 0.0)

    fiber = FiberInput(
        PiecewiseProfile(copy(breaks), copy(pieces_Rb)),
        PiecewiseProfile(copy(breaks), copy(pieces_theta)),
        PiecewiseProfile(copy(breaks), copy(pieces_dtwist)),
        k2  -> k2,
        tau -> tau
    )

    return (; fiber, breaks, title, input_state, output, n)
end

function demofiber1()
    lead_in = 4.0
    spacer = 8.5
    lead_out = 4.0
    loop_turns = [1.0, 1.0, 1.0]
    target_retardance = [π / 2, π, π / 2]
    loop_radii = 2π .* loop_turns ./ target_retardance
    loop_angles = [0.0, π / 2, 0.0]
    return build_piecewise_demo_fiber(
        lead_in,
        spacer,
        lead_out,
        loop_turns,
        loop_radii,
        loop_angles;
        title = "Fiber polarization control paddles",
        output = "fiberinput_paddles_3d.html"
    )
end

function demofiber2()
    lead_in = 3.0
    spacer = 6.0
    lead_out = 3.0
    loop_turns = [0.75, 1.0, 0.75, 0.5]
    target_retardance = [π / 2, π, π / 2, π / 3]
    loop_radii = 2π .* loop_turns ./ target_retardance
    loop_angles = [0.0, π / 3, 2π / 3, π / 6]
    return build_piecewise_demo_fiber(
        lead_in,
        spacer,
        lead_out,
        loop_turns,
        loop_radii,
        loop_angles;
        title = "Four-loop twisted fiber controller",
        input_state = ComplexF64[inv(sqrt(2)) + 0.0im, 0.0 + inv(sqrt(2)) * im],
        output = "fiberinput_fourloop_3d.html"
    )
end

function run_demo(demo = demofiber1(); output::Union{Nothing,AbstractString} = nothing)
    fiber = demo.fiber
    breaks = demo.breaks
    title = demo.title
    input_state = demo.input_state
    plot_output = isnothing(output) ? demo.output : output
    n = demo.n

    K = make_generator(fiber)
    J_final, stats = propagate_piecewise(
        K,
        copy(breaks);
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
        first(breaks),
        last(breaks);
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
