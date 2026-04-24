using MonteCarloMeasurements

include("material-properties.jl")
include("fiber-cross-section.jl")
include("path-geometry.jl")
include("path-geometry-plot.jl")
include("path-integral.jl")
include("fiber-path-plot.jl")


const DEMO_FIBER_CROSS_SECTION = FiberCrossSection(
    GermaniaSilicaGlass(0.036),
    GermaniaSilicaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like demo"
)
const DEMO_λ_M = 1550e-9
const DEMO_T_K = 297.15

# """
#     demofiber1()
# The shape of the fiber in this file is that of polarization paddles.
# The first loop is L/4, the second is L/2, and the third is L/4, where 
# L is the length of fiber in each loop. The angles of the loops are 0, π/2, 
# and 0, respectively. The input state is vertical polarization.
# """
# function demofiber1()
#     lead_in = 4.0
#     spacer = 8.5
#     lead_out = 4.0
#     loop_turns = [1.0, 1.0, 1.0]
#     target_retardance = [π / 2, π, π / 2]
#     loop_radii = 2π .* loop_turns ./ target_retardance
#     loop_angles = [0.0, π / 2, 0.0]

#     total_length = lead_in + sum(2π .* loop_radii .* loop_turns) + 2 * spacer + lead_out
#     spec = FiberSpec(0.0, total_length; cross_section = DEMO_FIBER_CROSS_SECTION, λ_m = DEMO_λ_M)

#     s = 0.0
#     bend!(spec, s, s + lead_in; angle = 0.0, axis = loop_angles[1], T_K = DEMO_T_K)
#     s += lead_in

#     loop1_angle = 2π * loop_turns[1]
#     loop1_length = loop_radii[1] * loop1_angle
#     bend!(spec, s, s + loop1_length; angle = loop1_angle, axis = loop_angles[1], T_K = DEMO_T_K)
#     s += loop1_length

#     bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[1], T_K = DEMO_T_K)
#     twist!(spec, s, s + spacer; rate = (loop_angles[2] - loop_angles[1]) / spacer, T_K = DEMO_T_K)
#     s += spacer

#     loop2_angle = 2π * loop_turns[2]
#     loop2_length = loop_radii[2] * loop2_angle
#     bend!(spec, s, s + loop2_length; angle = loop2_angle, axis = loop_angles[2], T_K = DEMO_T_K)
#     s += loop2_length

#     bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[2], T_K = DEMO_T_K)
#     twist!(spec, s, s + spacer; rate = (loop_angles[3] - loop_angles[2]) / spacer, T_K = DEMO_T_K)
#     s += spacer

#     loop3_angle = 2π * loop_turns[3]
#     loop3_length = loop_radii[3] * loop3_angle
#     bend!(spec, s, s + loop3_length; angle = loop3_angle, axis = loop_angles[3], T_K = DEMO_T_K)
#     s += loop3_length

#     bend!(spec, s, s + lead_out; angle = 0.0, axis = loop_angles[3], T_K = DEMO_T_K)

#     return (
#         fiber = build(spec),
#         title = "Fiber polarization control paddles",
#         input_state = ComplexF64[0.0 + 0.0im, 1.0 + 0.0im],
#         output = "fiber-path-plot.html",
#         n = 4001
#     )
# end


# """
#     demofiber2()
# This fiber is a more complex version of the first demo, with four loops instead of three, 
# and with the loop angles changing by π/3 between each loop. The retardances of the loops 
# are π/2, π, π/2, and π/3. The input state is diagonal polarization.  
# """
# function demofiber2()
#     lead_in = 3.0
#     spacer = 6.0
#     lead_out = 3.0
#     loop_turns = [0.75, 1.0, 0.75, 0.5]
#     target_retardance = [π / 2, π, π / 2, π / 3]
#     loop_radii = 2π .* loop_turns ./ target_retardance
#     loop_angles = [0.0, π / 3, 2π / 3, π / 6]

#     total_length = lead_in + sum(2π .* loop_radii .* loop_turns) + 3 * spacer + lead_out
#     spec = FiberSpec(0.0, total_length; cross_section = DEMO_FIBER_CROSS_SECTION, λ_m = DEMO_λ_M)

#     s = 0.0
#     bend!(spec, s, s + lead_in; angle = 0.0, axis = loop_angles[1], T_K = DEMO_T_K)
#     s += lead_in

#     loop1_angle = 2π * loop_turns[1]
#     loop1_length = loop_radii[1] * loop1_angle
#     bend!(spec, s, s + loop1_length; angle = loop1_angle, axis = loop_angles[1], T_K = DEMO_T_K)
#     s += loop1_length

#     bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[1], T_K = DEMO_T_K)
#     twist!(spec, s, s + spacer; rate = (loop_angles[2] - loop_angles[1]) / spacer, T_K = DEMO_T_K)
#     s += spacer

#     loop2_angle = 2π * loop_turns[2]
#     loop2_length = loop_radii[2] * loop2_angle
#     bend!(spec, s, s + loop2_length; angle = loop2_angle, axis = loop_angles[2], T_K = DEMO_T_K)
#     s += loop2_length

#     bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[2], T_K = DEMO_T_K)
#     twist!(spec, s, s + spacer; rate = (loop_angles[3] - loop_angles[2]) / spacer, T_K = DEMO_T_K)
#     s += spacer

#     loop3_angle = 2π * loop_turns[3]
#     loop3_length = loop_radii[3] * loop3_angle
#     bend!(spec, s, s + loop3_length; angle = loop3_angle, axis = loop_angles[3], T_K = DEMO_T_K)
#     s += loop3_length

#     bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[3], T_K = DEMO_T_K)
#     twist!(spec, s, s + spacer; rate = (loop_angles[4] - loop_angles[3]) / spacer, T_K = DEMO_T_K)
#     s += spacer

#     loop4_angle = 2π * loop_turns[4]
#     loop4_length = loop_radii[4] * loop4_angle
#     bend!(spec, s, s + loop4_length; angle = loop4_angle, axis = loop_angles[4], T_K = DEMO_T_K)
#     s += loop4_length

#     bend!(spec, s, s + lead_out; angle = 0.0, axis = loop_angles[4], T_K = DEMO_T_K)

#     return (
#         fiber = build(spec),
#         title = "Four-loop twisted fiber controller",
#         input_state = ComplexF64[inv(sqrt(2)) + 0.0im, 0.0 + inv(sqrt(2)) * im],
#         output = "fiber-path-plot.html",
#         n = 4001
#     )
# end

# function run_demo(demo ; output::Union{Nothing,AbstractString} = nothing)
#     fiber = demo.fiber
#     title = demo.title
#     input_state = demo.input_state
#     plot_output = isnothing(output) ? demo.output : output
#     n = demo.n

#     J_final, stats = propagate_fiber(
#         fiber;
#         rtol = 1e-10,
#         atol = 1e-12,
#         h_init = 1e-2
#     )

#     println("Final Jones matrix:")
#     println(J_final)

#     println("\nInterval stats:")
#     for (i, st) in enumerate(stats)
#         println("interval $i: accepted=$(st.accepted_steps), rejected=$(st.rejected_steps)")
#     end

#     plot_path = write_fiber_input_plot3d(
#         fiber,
#         fiber.s_start,
#         fiber.s_end;
#         n = n,
#         output = plot_output,
#         title = title,
#         input_state = input_state
#     )

#     println("\nWrote interactive 3D plot to:")
#     println(plot_path)

#     return (; demo..., J_final, stats, plot_path)
# end

"""
    demo_fiber_path(; output, n, title)

Demonstrate the analytic 3D path layer in [`path-geometry.jl`](path-geometry.jl) (via the
`PathGeometry` module in [`path-geometry-plot.jl`](path-geometry-plot.jl)): a short centerline
with circular bends, a catenary segment, and a material twist overlay on a straight spacer.
Writes an interactive Plotly HTML file.
"""
function demo_fiber_path(;
    output::AbstractString = joinpath(@__DIR__, "..", "output", "path-geometry-demo.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Fiber path geometry: bends, catenary, twist",
)
    PG = PathGeometry
    spec = PG.PathSpec()
    L_lead = 0.10
    R1 = 0.05
    θ1 = π / 2
    L_spacer = 0.12
    PG.straight!(spec; length = L_lead)
    PG.bend!(spec; radius = R1, angle = θ1)
    PG.straight!(spec; length = L_spacer)
    PG.catenary!(spec; a = 0.03, length = 0.10, axis_angle = 0.0)
    PG.bend!(spec; radius = 0.06, angle = π / 3)
    PG.straight!(spec; length = 0.08)
    s_twist = L_lead + R1 * abs(θ1)
    PG.twist!(spec; s_start = s_twist, length = L_spacer, rate = 35.0)
    path = PG.build(spec)
    println("Arc length (effective): ", PG.path_length(path))
    println("Writhe: ", PG.writhe(path))
    plot_path = write_path_geometry_plot3d(
        path,
        path.s_start,
        path.s_end;
        fidelity = fidelity,
        output = output,
        title = title,
    )
    println("Wrote path geometry plot to:")
    println(plot_path)
    return (; path, plot_path)
end

"""
    demo_fiber_path_segment_labels(; output, fidelity, title)

Illustrate optional **segment nicknames**: `straight!`, `bend!`, `catenary!`, and `helix!` are
called with `nickname` strings. [`write_path_geometry_plot3d`](path-geometry-plot.jl) renders each
name as 3D text near the segment midpoint (offset along the principal normal, in the osculating
plane). Uses only straight / bend / catenary / helix segments so the demo does not depend on jump
connectors or `min_bend_radius`.
"""
function demo_fiber_path_segment_labels(;
    output::AbstractString = joinpath(@__DIR__, "..", "output", "path-geometry-segment-labels-demo.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Path geometry: segment nicknames",
)
    PG = PathGeometry
    spec = PG.PathSpec()
    PG.straight!(spec; length = 0.08, nickname = "lead-in")
    PG.bend!(spec; radius = 0.06, angle = π / 2, nickname = "90° bend")
    PG.straight!(spec; length = 0.06, nickname = "spacer")
    PG.catenary!(spec; a = 0.04, length = 0.08, axis_angle = 0.0, nickname = "sag")
    PG.helix!(spec; radius = 0.025, pitch = 0.015, turns = 1.2, axis_angle = 0.0, nickname = "twist section")
    PG.straight!(spec; length = 0.06, nickname = "lead-out")
    path = PG.build(spec)
    println("Arc length (effective): ", PG.path_length(path), " m")
    plot_path = write_path_geometry_plot3d(
        path,
        path.s_start,
        path.s_end;
        fidelity = fidelity,
        output = output,
        title = title,
    )
    println("Wrote segment-label demo to:")
    println(plot_path)
    return (; path, plot_path)
end

"""
    demo_fiber_path_helix(; output, n, title)

Demonstrate `HelixSegment` with three separate single-segment paths, each a helix
with a different `axis_angle` (0, π/3, 2π/3).  Each helix has the same radius and
pitch, so the shape is the same but the winding plane rotates around the fiber axis.
Writes an interactive Plotly HTML file via `write_path_geometry_plot3d`.
"""
function demo_fiber_path_helix(;
    output::AbstractString = joinpath(@__DIR__, "..", "output", "path-geometry-helix-demo.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "HelixSegment: three axis angles",
)
    PG = PathGeometry
    paths = map([0.0, π/3, 2π/3]) do axis_angle
        spec = PG.PathSpec()
        PG.straight!(spec; length = 0.05)
        PG.helix!(spec; radius = 0.03, pitch = 0.02, turns = 2.0, axis_angle)
        PG.straight!(spec; length = 0.05)
        PG.build(spec)
    end
    for (i, path) in enumerate(paths)
        println("Helix $(i) (axis_angle=$(round((i-1)*π/3, digits=3)) rad): " *
                "arc_length=$(round(PG.path_length(path), digits=4)) m")
    end
    base, ext = splitext(output)
    plot_paths = map(enumerate(paths)) do (i, path)
        angle_label = ["0", "pi_3", "2pi_3"][i]
        out = "$(base)_$(angle_label)$(ext)"
        write_path_geometry_plot3d(
            path, path.s_start, path.s_end;
            fidelity = fidelity,
            output = out,
            title = "$(title) — axis_angle=$(["0", "π/3", "2π/3"][i])",
        )
    end
    for p in plot_paths
        println("Wrote helix demo plot to: ", p)
    end
    return (; paths, plot_paths)
end

"""
    demo_fiber_path_jumps_min_radius(; output, fidelity, title)

Demonstrate `jumpby!` and `jumpto!` with focus on the min_bend_radius parameter.
"""
function demo_fiber_path_jumps_min_radius(;
    output::AbstractString = joinpath(@__DIR__, "..", "output", "path-geometry-jumps-demo.html"),
    fidelity::Float64 = 4.0,
    title::AbstractString = "JumpBy and JumpTo: Hermite connectors, min_bend_radius",
)
    PG = PathGeometry
    spec = PG.PathSpec()

    PG.straight!(spec; length = 1)
    PG.jumpto!(spec; destination = (1, 0.0, 1), tangent = (0.0, 0.0, -1.0),
            min_bend_radius = 0.4) # T1: fails only if min_bend_radius >0.5
    PG.straight!(spec; length = 1)
    PG.jumpto!(spec; destination = (2, 0.0, 0), tangent = (0.0, 0.0, 1.0),
            min_bend_radius = 0.1) # T2: fails only if min_bend_radius >0.5
    PG.straight!(spec; length = 1)
    PG.jumpto!(spec; destination = (3, 0.0, 1), tangent = (0.0, 0.0, -1.0),
            min_bend_radius = 0.2) # T3: fails only if min_bend_radius >0.50
    PG.straight!(spec; length = 1)
    PG.jumpby!(spec; delta = (-1, 0.0, 0), tangent = (0.0, 0.0, -1.0),
            min_bend_radius = 0.1)
    PG.straight!(spec; length = 1)

    path = PG.build(spec)

    plot_path = write_path_geometry_plot3d(
        path, path.s_start, path.s_end;
        fidelity = fidelity, output   = output, title    = title )
    return (; path, plot_path)
end

"""
    demo_fiber_path_jumps2(; output, fidelity, title)

Demonstrate a more complex routing scenario with multiple `jumpby!` and `jumpto!`
calls interspersed with straights, a helix, and a circular bend.

Path layout:
  1. Short straight lead-in.
  2. **JumpBy** with an explicit `tangent_out` — the Hermite connector prescribes
     the exit direction (pointing diagonally upward), not the chord direction.
  3. Single-turn helix section.
  4. Short straight spacer.
  5. **JumpBy** without `tangent_out` — exit tangent follows the chord.
  6. Quarter-circle bend.
  7. **JumpTo** to a fixed lab-frame waypoint with no prescribed tangent.
  8. Short straight.
  9. **JumpTo** to a second waypoint with a prescribed outgoing tangent (global +z).
  10. Straight lead-out.
"""
function demo_fiber_path_jumps2(;
    output::AbstractString = joinpath(@__DIR__, "..", "output", "path-geometry-jumps2-demo.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Multi-jump routing: JumpBy (with/without tangent_out) + JumpTo",
)
    PG = PathGeometry
    spec = PG.PathSpec()

    PG.straight!(spec; length = 0.05)

    # Prescribed tangent_out + min_bend_radius: handle length grows to keep κ ≤ 1/R.
    PG.jumpby!(spec; delta = (0.03, 0.0, 0.04), tangent = (0.0, 0.0, 1.0),
               min_bend_radius = 0.02)

    PG.helix!(spec; radius = 0.02, pitch = 0.01, turns = 1.0, axis_angle = 0.0)
    PG.straight!(spec; length = 0.04)

    # No tangent_out → exit tangent follows chord direction.
    PG.jumpby!(spec; delta = (0.05, 0.0, 0.03))

    PG.bend!(spec; radius = 0.04, angle = π / 2, axis_angle = 0.0)
    PG.jumpto!(spec; destination = (0.22, 0.0, 0.20))
    PG.straight!(spec; length = 0.03)
    PG.jumpto!(spec; destination = (0.22, 0.0, 0.28), tangent = (0.0, 0.0, 1.0),
               min_bend_radius = 0.1)
    PG.straight!(spec; length = 0.05)

    path = PG.build(spec)
    println("Arc length : ", round(PG.path_length(path) * 100; digits = 2), " cm")
    println("Start      : ", round.(PG.start_point(path) .* 100; digits = 1), " cm")
    println("End        : ", round.(PG.end_point(path) .* 100; digits = 1), " cm")
    println("Writhe     : ", round(PG.writhe(path); digits = 6))

    plot_path = write_path_geometry_plot3d(
        path, path.s_start, path.s_end;
        fidelity = fidelity,
        output   = output,
        title    = title,
    )
    println("Wrote jumps2 demo to: ", plot_path)
    return (; path, plot_path)
end

"""
    demo_fiber_paddles_mcm(; μ_T = 297.15, σ_T = 2.0, N = 200, ...)

Four-paddle polarization controller (Λ/4 – Λ/2 – Λ/2 – Λ/4) with paddle 2's
temperature modeled as a normal random variable via MonteCarloMeasurements.
Propagates a horizontal input state through the fiber in a single
`propagate_fiber` call — the Particles-valued `T_K` on paddle 2 lifts through
`bending_birefringence` and the Jones integrator so the output state carries
an ensemble of `N` samples.

Emits two HTML artefacts:

1. `paddles-mcm-poincare.html` — Poincare sphere with the output-state point
   cloud (one marker per temperature sample, colored by T₂) and the
   ensemble-mean state vector drawn as an arrow.
2. `paddles-mcm-fiber.html` — 3D fiber centerline rendered via
   [`write_fiber_input_plot3d`](fiber-path-plot.jl) at the nominal temperature
   (the MCM ensemble collapses to its mean for the centerline view — the
   geometry is temperature-independent at this level; only the endpoint
   polarization scatters).
"""
function demo_fiber_paddles_mcm(;
    μ_T::Float64 = 297.15,
    σ_T::Float64 = 10.0,
    N::Int = 200,
    poincare_output::AbstractString = joinpath(@__DIR__, "..", "output", "paddles-mcm-poincare.html"),
    fiber_output::AbstractString = joinpath(@__DIR__, "..", "output", "paddles-mcm-fiber.html"),
)
    # Four-paddle controller: (turns, radius_m, axis_angle_rad).
    # Geometry matches the Thorlabs FPC030 small-diameter polarization controller
    # at 1550 nm with SMF-28-class fiber.
    paddles = [
        (turns = 1.0, radius = 0.030, axis = 0.0),
        (turns = 2.0, radius = 0.030, axis =  π/4),
        (turns = 2.0, radius = 0.030, axis = -π/4),
        (turns = 1.0, radius = 0.030, axis = 0.0),
    ]

    # Paddle 2's temperature as Particles. Explicit N-sample ensemble from a
    # normal distribution centered at μ_T with standard deviation σ_T.
    T_paddle2 = Particles(N, MonteCarloMeasurements.Distributions.Normal(μ_T, σ_T))
    paddle2_temperatures = T_paddle2.particles  # Vector{Float64}, length N

    # Assemble fiber arc-length domain [0, L_total] with one low-level custom
    # bend segment per paddle.
    s = 0.0
    segment_ranges = Tuple{Float64,Float64}[]
    path_spec = PathSpec()
    for p in paddles
        L = 2π * p.radius * p.turns
        push!(segment_ranges, (s, s + L))
        bend!(path_spec; radius = p.radius, angle = 2π * p.turns, axis_angle = p.axis)
        s += L
    end
    L_total = s

    path = build(path_spec)
    function T_profile(s)
        for (i, (s0, s1)) in enumerate(segment_ranges)
            if s0 <= s <= s1
                return i == 2 ? T_paddle2 : μ_T
            end
        end
        return μ_T
    end
    fiber = Fiber(path; cross_section = DEMO_FIBER_CROSS_SECTION)

    # Propagate once with in-band MCM. unsafe_comparisons lets branching inside
    # the adaptive integrator operate on Particles-valued scalars.
    MonteCarloMeasurements.unsafe_comparisons(true)
    J = try
        J_, _stats = propagate_fiber(fiber; λ_m = DEMO_λ_M, T_K = T_profile, rtol = 1e-9, atol = 1e-12, h_init = 0.1)
        J_
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end

    H_STATE = ComplexF64[1.0 + 0.0im, 0.0 + 0.0im]
    ψ_out = J * H_STATE                                  # 2-vector of Complex{Particles}

    # Extract the per-sample Jones vectors.
    samples_re = [real(ψ_out[k]).particles for k in 1:2]
    samples_im = [imag(ψ_out[k]).particles for k in 1:2]
    ψ_samples = [ComplexF64[samples_re[1][j] + im * samples_im[1][j],
                            samples_re[2][j] + im * samples_im[2][j]] for j in 1:N]
    stokes_samples = [stokes_from_jones(ψ_samples[j]) for j in 1:N]

    # Mean output state — normalize after averaging; Stokes of the mean is
    # built directly from pmean so phase is consistent.
    ψ_mean_unnorm = ComplexF64[
        pmean(real(ψ_out[1])) + im * pmean(imag(ψ_out[1])),
        pmean(real(ψ_out[2])) + im * pmean(imag(ψ_out[2])),
    ]
    ψ_mean = ψ_mean_unnorm / norm(ψ_mean_unnorm)
    rep_mean = poincare_vector_representation(ψ_mean)

    _write_poincare_cloud_plot(
        rep_mean, stokes_samples, paddle2_temperatures;
        output = poincare_output,
        title = "Four-paddle controller — Poincare sphere MCM cloud (σ_T2 = $(σ_T) K, N = $(N))",
    )

    # Build a scalar-temperature fiber for the 3D centerline view.
    fiber_mean = Fiber(path; cross_section = DEMO_FIBER_CROSS_SECTION)
    plot_fiber = write_fiber_input_plot3d(
        fiber_mean, 0.0, L_total;
        λ_m = DEMO_λ_M,
        T_K = μ_T,
        n = 801,
        output = fiber_output,
        title = "Four-paddle controller — fiber centerline (nominal T)",
        input_state = H_STATE,
    )

    println("Wrote Poincare MCM cloud to:    ", poincare_output)
    println("Wrote fiber centerline plot to: ", plot_fiber)

    return (
        fiber = fiber,
        J_mean_sample = ψ_mean,
        ψ_samples = ψ_samples,
        stokes_samples = stokes_samples,
        paddle2_temperatures = paddle2_temperatures,
        plot_poincare = poincare_output,
        plot_fiber = plot_fiber,
    )
end

# Minimal Plotly HTML writer for a Poincare sphere with a per-sample scatter
# cloud. Uses `render_poincare_sphere` for the sphere/equator scaffolding and
# adds the point-cloud trace on top.
function _write_poincare_cloud_plot(
    rep_mean::NamedTuple,
    stokes_samples::Vector,
    sample_colors::Vector{Float64};
    output::AbstractString,
    title::AbstractString,
)
    sphere = render_poincare_sphere(rep_mean)
    xs = [s[1] for s in stokes_samples]
    ys = [s[2] for s in stokes_samples]
    zs = [s[3] for s in stokes_samples]

    js_num(x) = isnan(x) ? "NaN" : string(Float64(x))
    js_arr(xs) = "[" * join(js_num.(xs), ",") * "]"
    js_surf(xss) = "[" * join([js_arr(r) for r in xss], ",") * "]"
    js_strarr(xs) = "[" * join(["\"$(x)\"" for x in xs], ",") * "]"

    html = """
<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>$title</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>html,body{margin:0;padding:0;width:100%;height:100%;font-family:sans-serif;}#plot{width:100%;height:100%;}</style>
</head><body><div id="plot"></div>
<script>
const sphU = $(js_surf(sphere.surface.x));
const sphV = $(js_surf(sphere.surface.y));
const sphW = $(js_surf(sphere.surface.z));
const eqXYx = $(js_arr(sphere.circles.xy.x));
const eqXYy = $(js_arr(sphere.circles.xy.y));
const eqXYz = $(js_arr(sphere.circles.xy.z));
const eqXZx = $(js_arr(sphere.circles.xz.x));
const eqXZy = $(js_arr(sphere.circles.xz.y));
const eqXZz = $(js_arr(sphere.circles.xz.z));
const eqYZx = $(js_arr(sphere.circles.yz.x));
const eqYZy = $(js_arr(sphere.circles.yz.y));
const eqYZz = $(js_arr(sphere.circles.yz.z));
const vecX = $(js_arr(sphere.vector.x));
const vecY = $(js_arr(sphere.vector.y));
const vecZ = $(js_arr(sphere.vector.z));
const labX = $(js_arr(sphere.labels.x));
const labY = $(js_arr(sphere.labels.y));
const labZ = $(js_arr(sphere.labels.z));
const labT = $(js_strarr(sphere.labels.text));
const cloudX = $(js_arr(xs));
const cloudY = $(js_arr(ys));
const cloudZ = $(js_arr(zs));
const cloudC = $(js_arr(sample_colors));

const traces = [
  {type:"surface", x:sphU, y:sphV, z:sphW, opacity:0.12, showscale:false,
   colorscale:[[0,"#dddddd"],[1,"#dddddd"]], hoverinfo:"skip", contours:{x:{highlight:false},y:{highlight:false},z:{highlight:false}}},
  {type:"scatter3d", mode:"lines", x:eqXYx, y:eqXYy, z:eqXYz, line:{width:2,color:"#888"}, hoverinfo:"skip", showlegend:false},
  {type:"scatter3d", mode:"lines", x:eqXZx, y:eqXZy, z:eqXZz, line:{width:2,color:"#888"}, hoverinfo:"skip", showlegend:false},
  {type:"scatter3d", mode:"lines", x:eqYZx, y:eqYZy, z:eqYZz, line:{width:2,color:"#888"}, hoverinfo:"skip", showlegend:false},
  {type:"scatter3d", mode:"markers", x:cloudX, y:cloudY, z:cloudZ,
   marker:{size:3, color:cloudC, colorscale:"Viridis", showscale:true, colorbar:{title:"T₂ (K)"}, opacity:0.85},
   name:"samples"},
  {type:"scatter3d", mode:"lines", x:vecX, y:vecY, z:vecZ, line:{width:6,color:"#c00"}, name:"mean state"},
  {type:"scatter3d", mode:"markers", x:[vecX[1]], y:[vecY[1]], z:[vecZ[1]],
   marker:{size:6, color:"#c00"}, showlegend:false},
  {type:"scatter3d", mode:"text", x:labX, y:labY, z:labZ, text:labT, textfont:{size:14}, hoverinfo:"skip", showlegend:false}
];
const layout = {
  title: "$title",
  scene: {aspectmode:"cube",
          xaxis:{title:"S1", range:[-1.3,1.3]},
          yaxis:{title:"S2", range:[-1.3,1.3]},
          zaxis:{title:"S3", range:[-1.3,1.3]}},
  margin:{l:0,r:0,t:40,b:0}
};
Plotly.newPlot("plot", traces, layout, {responsive:true});
</script></body></html>
"""
    open(output, "w") do io
        write(io, html)
    end
    return output
end

const DEMO_INDEX = [
    (
        fn   = demo_fiber_path,
        kwargs = NamedTuple(),
        desc = "Mixed-segment path: straight leads, circular bends, a catenary sag, and a " *
               "material twist overlay — illustrates the segment assembly API and " *
               "the Frenet–Serret sliding frame.",
    ),
    (
        fn   = demo_fiber_path_segment_labels,
        kwargs = NamedTuple(),
        desc = "Same style of path as the mixed-segment demo, but each segment has a `nickname` " *
               "string; the HTML plot shows those names as 3D labels offset in the osculating plane.",
    ),
    (
        fn   = demo_fiber_path_helix,
        kwargs = NamedTuple(),
        desc = "Three HelixSegment paths sharing the same radius and pitch but with " *
               "axis_angle ∈ {0, π/3, 2π/3} — shows how the winding plane rotates " *
               "while the entry tangent stays aligned with the incoming fiber direction.",
    ),
    (
        fn   = demo_fiber_path_jumps_min_radius,
        kwargs = NamedTuple(),
        desc = "Demonstrate `jumpby!` and `jumpto!` with focus on the min_bend_radius parameter." 
    ),
    (
        fn   = demo_fiber_path_jumps2,
        kwargs = NamedTuple(),
        desc = "Multi-jump routing: two JumpBy calls (one with, one without tangent_out) " *
               "and two JumpTo calls (one with, one without tangent_out) interleaved with " *
               "a helix, a circular bend, and straight spacers — exercises every " *
               "combination of jump type × tangent prescription.",
    ),
    (
        fn   = demo_fiber_paddles_mcm,
        kwargs = NamedTuple(),
        desc = "Four-paddle polarization controller (Λ/4 – Λ/2 – Λ/2 – Λ/4) with paddle 2's " *
               "temperature drawn from a normal distribution. Uses MonteCarloMeasurements to " *
               "propagate the ensemble through the Jones integrator in one shot, then plots " *
               "the output-state point cloud on the Poincaré sphere and the fiber centerline " *
               "with the mean-temperature polarization trajectory.",
    ),
]

"""
    demo_all(; index_output)

Run every demo in `DEMO_INDEX` and write an `index.html` that links to each output file
with a short description of what it illustrates.
"""
function demo_all(; index_output::AbstractString = joinpath(@__DIR__, "..", "output", "index.html"))
    entries = Tuple{String, String, String}[]   # (title, rel_path, desc)

    for d in DEMO_INDEX
        result = d.fn(; d.kwargs...)
        # Collect all HTML output paths from the result (handles single or vector).
        paths = result isa NamedTuple ? values(result) : (result,)
        html_paths = String[]
        for v in paths
            if v isa AbstractString && endswith(v, ".html")
                push!(html_paths, v)
            elseif v isa AbstractVector
                for item in v
                    item isa AbstractString && endswith(item, ".html") && push!(html_paths, item)
                end
            end
        end
        for p in html_paths
            push!(entries, (basename(p), p, d.desc))
        end
    end

    open(index_output, "w") do io
        println(io, """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>BIFROST path-geometry demos</title>
  <style>
    body { font-family: sans-serif; max-width: 800px; margin: 2em auto; color: #222; }
    h1   { font-size: 1.5em; border-bottom: 1px solid #ccc; padding-bottom: 0.3em; }
    ul   { padding-left: 1.2em; }
    li   { margin: 1em 0; }
    a    { font-weight: bold; color: #1a6; }
    p.desc { margin: 0.3em 0 0 0; color: #555; font-size: 0.95em; }
  </style>
</head>
<body>
  <h1>BIFROST path-geometry demos</h1>
  <ul>""")
        for (title, path, desc) in entries
            println(io, "    <li>")
            println(io, "      <a href=\"$(path)\">$(title)</a>")
            println(io, "      <p class=\"desc\">$(desc)</p>")
            println(io, "    </li>")
        end
        println(io, """  </ul>
</body>
</html>""")
    end

    println("Wrote demo index to: ", index_output)
    return index_output
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_all()
end
