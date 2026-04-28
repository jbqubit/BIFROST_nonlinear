using MonteCarloMeasurements

include("material-properties.jl")
include("fiber-cross-section.jl")
include("path-geometry.jl")
include("path-geometry-plot.jl")
include("path-integral.jl")
include("fiber-path-plot.jl")
include("fiber-path-modify.jl")


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
    spec = PG.PathSpecBuilder()
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
    # TODO: twist refactor —     PG.twist!(spec; s_start = s_twist, length = L_spacer, rate = 35.0)
    path = PG.build(spec)
    println("Arc length (effective): ", PG.path_length(path))
    println("Writhe: ", PG.writhe(path))
    plot_path = write_path_geometry_plot3d(
        path,
        path.spec.s_start,
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
    _nick(s) = PG.AbstractMeta[PG.Nickname(s)]
    spec = PG.PathSpecBuilder()
    PG.straight!(spec; length = 0.08, meta = _nick("lead-in"))
    PG.bend!(spec; radius = 0.06, angle = π / 2, meta = _nick("90° bend"))
    PG.straight!(spec; length = 0.06, meta = _nick("spacer"))
    PG.catenary!(spec; a = 0.04, length = 0.08, axis_angle = 0.0, meta = _nick("sag"))
    PG.helix!(spec; radius = 0.025, pitch = 0.015, turns = 1.2, axis_angle = 0.0, meta = _nick("twist section"))
    PG.straight!(spec; length = 0.06, meta = _nick("lead-out"))
    path = PG.build(spec)
    println("Arc length (effective): ", PG.path_length(path), " m")
    plot_path = write_path_geometry_plot3d(
        path,
        path.spec.s_start,
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
        spec = PG.PathSpecBuilder()
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
            path, path.spec.s_start, path.s_end;
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
    spec = PG.PathSpecBuilder()

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
        path, path.spec.s_start, path.s_end;
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
    spec = PG.PathSpecBuilder()

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
        path, path.spec.s_start, path.s_end;
        fidelity = fidelity,
        output   = output,
        title    = title,
    )
    println("Wrote jumps2 demo to: ", plot_path)
    return (; path, plot_path)
end

"""
    demo_fiber_path_modify(; output_dir)

Exercise [`modify`](fiber-path-modify.jl) on a common inverted-U baseline path
(straight up → bend π → straight down) by attaching MCM annotations to a
single target segment and comparing the modified paths against the
unperturbed baseline. Produces HTML files (one per experiment):

1. `modify-straight-length.html` — shortens the first straight via
   `MCMadd(:length, −0.4)`.
2. `modify-bend-radius.html` — halves and doubles the bend radius via
   `MCMmul(:radius, −0.5)` / `MCMmul(:radius, 1.0)`.
3. `modify-bend-angle.html` — reduces and extends the bend angle via
   `MCMadd(:angle, −π/2)` / `MCMadd(:angle, π/4)`.

In every figure the unmodified segments of a variant are drawn in green and
the modified segment in red; variants are offset along *x* so the baseline
and each perturbation appear side-by-side in a single 3D scene.
"""
function demo_fiber_path_modify(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "output"),
)
    L, R = 1.0, 0.5

    row1 = _modify_row_html(
        joinpath(output_dir, "modify-straight-length.html"),
        "modify(:length) on the first straight segment";
        L = L, R = R,
        variants = [
            ("baseline",                      1, AbstractMeta[]),
            ("MCMadd(:length, -0.4)",         1, AbstractMeta[MCMadd(:length, -0.4)]),
        ],
    )

    row2 = _modify_row_html(
        joinpath(output_dir, "modify-bend-radius.html"),
        "modify(:radius) on the bend segment";
        L = L, R = R,
        variants = [
            ("baseline",                      2, AbstractMeta[]),
            ("MCMadd(:radius, -0.25)",        2, AbstractMeta[MCMadd(:radius, -0.25)]),
            ("MCMadd(:radius, +0.50)",        2, AbstractMeta[MCMadd(:radius,  0.50)]),
        ],
    )

    row3 = _modify_row_html(
        joinpath(output_dir, "modify-bend-angle.html"),
        "modify(:angle) on the bend segment";
        L = L, R = R,
        variants = [
            ("baseline",                      2, AbstractMeta[]),
            ("MCMadd(:angle, -π/2)",          2, AbstractMeta[MCMadd(:angle, -π/2)]),
            ("MCMadd(:angle, +π/4)",          2, AbstractMeta[MCMadd(:angle,  π/4)]),
            ("MCMadd(:angle, +π)",          2, AbstractMeta[MCMadd(:angle,  π)]),
        ],
    )

    row4 = _modify_row_html(
        joinpath(output_dir, "modify-straight-length-mul.html"),
        "modify(:length) on the first straight segment — MCMmul";
        L = L, R = R,
        variants = [
            ("baseline",                      1, AbstractMeta[]),
            ("MCMmul(:length, -0.4)",         1, AbstractMeta[MCMmul(:length, -0.4)]),
            ("MCMmul(:length, +0.5)",         1, AbstractMeta[MCMmul(:length,  0.5)]),
        ],
    )

    row5 = _modify_row_html(
        joinpath(output_dir, "modify-bend-radius-mul.html"),
        "modify(:radius) on the bend segment — MCMmul";
        L = L, R = R,
        variants = [
            ("baseline",                      2, AbstractMeta[]),
            ("MCMmul(:radius, 0.5)",          2, AbstractMeta[MCMmul(:radius, 0.5)]),
#            ("MCMmul(:radius, -0.5)",         2, AbstractMeta[MCMmul(:radius, -0.5)]), # neg radius not supported, throws @assertion error
            ("MCMmul(:radius, 2.0)",          2, AbstractMeta[MCMmul(:radius, 2.0)]),
        ],
    )

    row6 = _modify_row_html(
        joinpath(output_dir, "modify-bend-angle-mul.html"),
        "modify(:angle) on the bend segment — MCMmul";
        L = L, R = R,
        variants = [
            ("baseline",                      2, AbstractMeta[]),
            ("MCMmul(:angle, 0.5)",           2, AbstractMeta[MCMmul(:angle, 0.5)]),
            ("MCMmul(:angle, 1.25)",          2, AbstractMeta[MCMmul(:angle, 1.25)]),
        ],
    )

    # -------------------------------------------------------------------
    # Rows 7–12: helix-baseline (straight · bend · helix · straight).
    # Rows 7–9 use MCMadd on the helix's scalar fields; rows 10–12 are the
    # MCMmul counterparts. Target index 3 = helix segment.
    # -------------------------------------------------------------------

    row7 = _modify_row_html(
        joinpath(output_dir, "modify-helix-radius.html"),
        "modify(:radius) on the helix segment — MCMadd";
        L = L, R = R,
        variant_spacing = 3.5,
        builder = _build_modify_variant_helix,
        variants = [
            ("baseline",                      3, AbstractMeta[]),
            ("MCMadd(:radius, -0.05)",        3, AbstractMeta[MCMadd(:radius, -0.05)]),
            ("MCMadd(:radius, +0.10)",        3, AbstractMeta[MCMadd(:radius,  0.10)]),
        ],
    )

    row8 = _modify_row_html(
        joinpath(output_dir, "modify-helix-pitch.html"),
        "modify(:pitch) on the helix segment — MCMadd";
        L = L, R = R,
        variant_spacing = 3.5,
        builder = _build_modify_variant_helix,
        variants = [
            ("baseline",                      3, AbstractMeta[]),
            ("MCMadd(:pitch, -0.10)",         3, AbstractMeta[MCMadd(:pitch, -0.10)]),
            ("MCMadd(:pitch, +0.20)",         3, AbstractMeta[MCMadd(:pitch,  0.20)]),
        ],
    )

    row9 = _modify_row_html(
        joinpath(output_dir, "modify-helix-turns.html"),
        "modify(:turns) on the helix segment — MCMadd";
        L = L, R = R,
        variant_spacing = 3.5,
        builder = _build_modify_variant_helix,
        variants = [
            ("baseline",                      3, AbstractMeta[]),
            ("MCMadd(:turns, -0.5)",          3, AbstractMeta[MCMadd(:turns, -0.5)]),
            ("MCMadd(:turns, +0.5)",          3, AbstractMeta[MCMadd(:turns,  0.5)]),
        ],
    )

    row10 = _modify_row_html(
        joinpath(output_dir, "modify-helix-radius-mul.html"),
        "modify(:radius) on the helix segment — MCMmul";
        L = L, R = R,
        variant_spacing = 3.5,
        builder = _build_modify_variant_helix,
        variants = [
            ("baseline",                      3, AbstractMeta[]),
            ("MCMmul(:radius, 0.5)",          3, AbstractMeta[MCMmul(:radius, 0.5)]),
            ("MCMmul(:radius, 2.0)",          3, AbstractMeta[MCMmul(:radius, 2.0)]),
        ],
    )

    row11 = _modify_row_html(
        joinpath(output_dir, "modify-helix-pitch-mul.html"),
        "modify(:pitch) on the helix segment — MCMmul";
        L = L, R = R,
        variant_spacing = 3.5,
        builder = _build_modify_variant_helix,
        variants = [
            ("baseline",                      3, AbstractMeta[]),
            ("MCMmul(:pitch, 0.5)",           3, AbstractMeta[MCMmul(:pitch, 0.5)]),
            ("MCMmul(:pitch, 2.0)",           3, AbstractMeta[MCMmul(:pitch, 2.0)]),
        ],
    )

    row12 = _modify_row_html(
        joinpath(output_dir, "modify-helix-turns-mul.html"),
        "modify(:turns) on the helix segment — MCMmul";
        L = L, R = R,
        variant_spacing = 3.5,
        builder = _build_modify_variant_helix,
        variants = [
            ("baseline",                      3, AbstractMeta[]),
            ("MCMmul(:turns, 0.67)",          3, AbstractMeta[MCMmul(:turns, 0.67)]),
            ("MCMmul(:turns, 1.5)",           3, AbstractMeta[MCMmul(:turns, 1.5)]),
        ],
    )

    return (
        plot_length            = row1,
        plot_radius            = row2,
        plot_angle             = row3,
        plot_length_mul        = row4,
        plot_radius_mul        = row5,
        plot_angle_mul         = row6,
        plot_helix_radius      = row7,
        plot_helix_pitch       = row8,
        plot_helix_turns       = row9,
        plot_helix_radius_mul  = row10,
        plot_helix_pitch_mul   = row11,
        plot_helix_turns_mul   = row12,
    )
end

# Build the baseline 3-segment inverted-U and attach `target_meta` to the
# segment indexed by `target_idx` (1 = first straight, 2 = bend, 3 = second
# straight). Returns the modified path and the target segment index.
function _build_modify_variant(L::Float64, R::Float64,
                               target_idx::Int, target_meta::Vector{AbstractMeta})
    spec = PathSpecBuilder()
    straight!(spec; length = L, meta = target_idx == 1 ? target_meta : AbstractMeta[])
    bend!(spec; radius = R, angle = π, axis_angle = 0.0,
          meta = target_idx == 2 ? target_meta : AbstractMeta[])
    straight!(spec; length = L, meta = target_idx == 3 ? target_meta : AbstractMeta[])
    path  = build(spec)
    fiber = Fiber(path; cross_section = DEMO_FIBER_CROSS_SECTION, T_ref_K = DEMO_T_K)
    return modify(fiber)
end

# Baseline helix parameters for the 4-segment variant (inverted-U plus a helix
# between the bend and the second straight). `axis_angle = 0` tilts the helix
# axis into the x–z plane so the coil projects recognisably onto the viewer's
# plane.
const _MODIFY_HELIX_RADIUS = 0.15
const _MODIFY_HELIX_PITCH  = 0.25
const _MODIFY_HELIX_TURNS  = 1.5

# Build a 4-segment baseline (straight · bend π · helix · straight) and attach
# `target_meta` to the segment indexed by `target_idx`:
#   1 = first straight, 2 = bend, 3 = helix, 4 = second straight.
function _build_modify_variant_helix(L::Float64, R::Float64,
                                     target_idx::Int, target_meta::Vector{AbstractMeta})
    spec = PathSpecBuilder()
    straight!(spec; length = L, meta = target_idx == 1 ? target_meta : AbstractMeta[])
    bend!(spec; radius = R, angle = π, axis_angle = 0.0,
          meta = target_idx == 2 ? target_meta : AbstractMeta[])
    helix!(spec; radius = _MODIFY_HELIX_RADIUS,
                 pitch  = _MODIFY_HELIX_PITCH,
                 turns  = _MODIFY_HELIX_TURNS,
                 axis_angle = 0.0,
                 meta = target_idx == 3 ? target_meta : AbstractMeta[])
    straight!(spec; length = L, meta = target_idx == 4 ? target_meta : AbstractMeta[])
    path  = build(spec)
    fiber = Fiber(path; cross_section = DEMO_FIBER_CROSS_SECTION, T_ref_K = DEMO_T_K)
    return modify(fiber)
end

# Sample a single placed segment's centerline in the global frame.
function _sample_segment_xyz(path::PathSpecCached, seg_index::Int; n::Int = 128)
    ps = path.placed_segments[seg_index]
    s0 = ps.s_offset_eff
    s1 = s0 + arc_length(ps.segment)
    ss = range(s0, s1; length = n)
    xs = Float64[]; ys = Float64[]; zs = Float64[]
    for s in ss
        p = position(path, s)
        push!(xs, Float64(p[1]))
        push!(ys, Float64(p[2]))
        push!(zs, Float64(p[3]))
    end
    return (x = xs, y = ys, z = zs)
end

# Render one experiment row: each variant's 3-segment path laid out side-by-
# side along x, with the target segment drawn in red and the others in green.
function _modify_row_html(output::AbstractString, title::AbstractString;
                          L::Float64, R::Float64,
                          variants::Vector,
                          variant_spacing::Float64 = 2.5,
                          builder = _build_modify_variant)
    js_num(x)   = isnan(x) ? "NaN" : string(Float64(x))
    js_arr(xs)  = "[" * join(js_num.(xs), ",") * "]"
    js_strarr(xs) = "[" * join(("\"" * replace(string(x), "\"" => "\\\"") * "\"" for x in xs), ",") * "]"

    trace_strs  = String[]
    label_xs    = Float64[]; label_ys = Float64[]; label_zs = Float64[]
    label_texts = String[]
    start_xs    = Float64[]; start_ys = Float64[]; start_zs = Float64[]

    for (k, (label, target_idx, target_meta)) in enumerate(variants)
        path = builder(L, R, target_idx, target_meta)
        dx = (k - 1) * variant_spacing
        n_segs = length(path.placed_segments)
        for i in 1:n_segs
            s     = _sample_segment_xyz(path, i)
            xs    = s.x .+ dx
            color = (i == target_idx && !isempty(target_meta)) ? "#e55" : "#6c6"
            name  = "$(label) — seg $(i)"
            push!(trace_strs, """
{
  type: 'scatter3d',
  mode: 'lines',
  x: $(js_arr(xs)),
  y: $(js_arr(s.y)),
  z: $(js_arr(s.z)),
  line: {color: '$(color)', width: 6},
  name: $(repr(name)),
  showlegend: $(i == 1)
}""")
        end
        # Green marker at the path start (origin of this variant).
        p0 = position(path, path.spec.s_start)
        push!(start_xs, Float64(p0[1]) + dx)
        push!(start_ys, Float64(p0[2]))
        push!(start_zs, Float64(p0[3]))

        # Label above each variant at the top of the inverted U.
        push!(label_xs, dx + R)
        push!(label_ys, 0.0)
        push!(label_zs, L + R + 0.15)
        push!(label_texts, label)
    end

    start_trace = """
{
  type: 'scatter3d',
  mode: 'markers',
  x: $(js_arr(start_xs)),
  y: $(js_arr(start_ys)),
  z: $(js_arr(start_zs)),
  marker: {color: '#6c6', size: 6, symbol: 'circle'},
  name: 'start',
  showlegend: false
}"""

    labels_trace = """
{
  type: 'scatter3d',
  mode: 'text',
  x: $(js_arr(label_xs)),
  y: $(js_arr(label_ys)),
  z: $(js_arr(label_zs)),
  text: $(js_strarr(label_texts)),
  textfont: {color: '#ddd', size: 13},
  showlegend: false
}"""

    html = """<!DOCTYPE html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><title>$(title)</title>
<script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>
<style>html,body{margin:0;padding:0;width:100%;height:100%;background:#111;color:#eee;font-family:sans-serif;}#plot{width:100%;height:100%;}</style>
</head><body><div id=\"plot\"></div>
<script>
const traces = [$(join(trace_strs, ",\n")),
$(start_trace)$(isempty(label_texts) ? "" : ",\n")$(labels_trace)];
const layout = {
  title: {text: $(repr(title)), font: {color: '#eee'}},
  paper_bgcolor: '#111',
  plot_bgcolor:  '#111',
  scene: {
    bgcolor: '#111',
    xaxis: {gridcolor: '#333', color: '#aaa'},
    yaxis: {gridcolor: '#333', color: '#aaa'},
    zaxis: {gridcolor: '#333', color: '#aaa'},
    aspectmode: 'data',
    camera: {eye: {x: 0.0, y: -2.5, z: 0.5}}
  },
  margin: {l:0, r:0, t:40, b:0},
  legend: {font: {color: '#eee'}}
};
Plotly.newPlot('plot', traces, layout);
</script>
</body></html>
"""
    open(output, "w") do io
        write(io, html)
    end
    return output
end



"""
    demo_fiber_paddles_mcm(; μ_T = 297.15, σ_T = 2.0, N = 200, ...)

Four-paddle polarization controller (Λ/4 – Λ/2 – Λ/2 – Λ/4) with the fiber
reference temperature modeled as a normal random variable via
MonteCarloMeasurements. Propagates a horizontal input state through the fiber
in a single `propagate_fiber` call — the Particles-valued `T_ref_K` lifts
through `bending_birefringence` and the Jones integrator so the output state
carries an ensemble of `N` samples.

Note: per-segment temperature variation (e.g., only paddle 2 uncertain) will
return once path-geometry gains segment-level metadata (MetaList/MCM). For now
the whole fiber shares one Particles-valued T.

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

    # Fiber reference temperature as Particles. Explicit N-sample ensemble from
    # a normal distribution centered at μ_T with standard deviation σ_T.
    T_ensemble = Particles(N, MonteCarloMeasurements.Distributions.Normal(μ_T, σ_T))
    paddle2_temperatures = T_ensemble.particles  # Vector{Float64}, length N

    # Assemble fiber arc-length domain [0, L_total] with one low-level custom
    # bend segment per paddle.
    s = 0.0
    path_spec = PathSpecBuilder()
    for p in paddles
        L = 2π * p.radius * p.turns
        bend!(path_spec; radius = p.radius, angle = 2π * p.turns, axis_angle = p.axis)
        s += L
    end
    L_total = s

    path = build(path_spec)
    fiber = Fiber(path; cross_section = DEMO_FIBER_CROSS_SECTION, T_ref_K = T_ensemble)

    # Propagate once with in-band MCM. unsafe_comparisons lets branching inside
    # the adaptive integrator operate on Particles-valued scalars.
    MonteCarloMeasurements.unsafe_comparisons(true)
    J = try
        J_, _stats = propagate_fiber(fiber; λ_m = DEMO_λ_M, rtol = 1e-9, atol = 1e-12, h_init = 0.1)
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

"""
    demo_seven_segment_mcm_temperature(; ΔT_K, N, ...)

End-to-end MCM stack test: a 7-segment fiber whose middle segment (1000 m
straight, with a constant material-twist overlay) carries a per-segment
`MCMadd(:T_K, ΔT)` annotation drawn from `Uniform(-ΔT_K, +ΔT_K)`.

This demo's purpose is to surface problems in the full meta → modify →
propagate stack. As of writing it has surfaced three:

1. (FIXED) `path_segment_breakpoints` / `path_twist_breakpoints` in
   `path-geometry.jl` previously coerced every breakpoint to `Float64`,
   crashing as soon as a perturbed segment's `s_offset_eff` or
   `arc_length` carried `Particles`. Now both functions accept `<:Real`.

2. (OPEN) `_find_segment_index` in `fiber-path-modify.jl` returns the
   first segment whose `seg_old_end >= s_eff_start - 1e-12`, so a twist
   overlay whose `s_start` exactly equals a segment boundary is mapped
   onto the *earlier* segment (with αi = 1) and the temperature signal
   on the twist window is silently dropped. Workaround in this demo:
   start the overlay 0.5 m into the middle segment.

3. (OPEN) `propagate_piecewise` reduces interval bounds to `Float64` via
   `scalar_reduce` before calling `propagate_interval!`. For a constant-K
   interval, this discards the interval-length Particles factor entirely
   — a uniform-rate twist over an uncertain-length straight produces
   zero output spread even though the optical path length is uncertain.
   Surfacing the ΔT signal on this demo's geometry would require either
   a non-constant K(s) on the middle segment (e.g. a bend) or a change
   to how length-uncertainty composes with constant generators.
"""
function demo_seven_segment_mcm_temperature(;
    ΔT_K::Float64 = 5.0,
    N::Int = 200,
    twist_rate::Float64 = 0.05,        # rad / m on the middle segment
    poincare_output::AbstractString =
        joinpath(@__DIR__, "..", "output", "seven-seg-mcm-temperature-poincare.html"),
    title::AbstractString = "Seven-segment fiber — middle-segment ΔT MCM (±$(ΔT_K) K, N=$(N))",
)
    # TODO: twist refactor — this demo depends on a material-twist overlay
    # to surface segment-length uncertainty as polarization spread. Pending
    # the per-segment-meta twist subsystem, the demo is skipped.
    @info "demo_seven_segment_mcm_temperature: skipped pending twist refactor"
    return (; poincare = "")
    L_lead   = 10.0
    L_spacer =  5.0
    L_mid    = 1000.0
    R_bend   =  5.0
    angle    = π / 12

    ΔT_ensemble = Particles(N, MonteCarloMeasurements.Distributions.Uniform(-ΔT_K, ΔT_K))
    ΔT_samples  = ΔT_ensemble.particles

    spec = PathSpecBuilder()
    straight!(spec; length = L_lead,   meta = AbstractMeta[Nickname("lead-in")])
    bend!(spec; radius = R_bend, angle = angle, axis_angle = 0.0,
                meta = AbstractMeta[Nickname("bend-1")])
    straight!(spec; length = L_spacer, meta = AbstractMeta[Nickname("spacer-1")])
    # Segment 4: the long middle straight with per-segment ΔT MCM.
    straight!(spec; length = L_mid,
              meta = AbstractMeta[Nickname("middle"),
                                  MCMadd(:T_K, ΔT_ensemble)])
    straight!(spec; length = L_spacer, meta = AbstractMeta[Nickname("spacer-2")])
    bend!(spec; radius = R_bend, angle = -angle, axis_angle = 0.0,
                meta = AbstractMeta[Nickname("bend-2")])
    straight!(spec; length = L_lead,   meta = AbstractMeta[Nickname("lead-out")])

    # Constant material twist sitting strictly inside the middle segment so
    # `_find_segment_index` in fiber-path-modify.jl picks segment 4 (it
    # otherwise returns the earlier segment for boundary-touching s values,
    # which leaves αi = 1.0 and discards the temperature signal on the
    # twist overlay).
    s_mid_start = L_lead + R_bend * angle + L_spacer + 0.5
    # TODO: twist refactor —     twist!(spec; s_start = s_mid_start, length = L_mid - 1.0, rate = twist_rate)

    path = build(spec)
    fiber_baseline = Fiber(path;
        cross_section = DEMO_FIBER_CROSS_SECTION, T_ref_K = DEMO_T_K)

    # Apply per-segment :T_K MCM by thermally rescaling segment geometries.
    MonteCarloMeasurements.unsafe_comparisons(true)
    result = try
        path_mod  = modify(fiber_baseline)
        fiber_mod = Fiber(path_mod;
            cross_section = DEMO_FIBER_CROSS_SECTION, T_ref_K = DEMO_T_K)
        J_, stats_ = propagate_fiber(fiber_mod;
            λ_m  = DEMO_λ_M,
            rtol = 1e-7,
            atol = 1e-10,
            h_init = 0.5)
        (J = J_, stats = stats_, fiber = fiber_mod)
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
    J = result.J

    H_STATE = ComplexF64[1.0 + 0.0im, 0.0 + 0.0im]
    ψ_out   = J * H_STATE

    samples_re = [real(ψ_out[k]).particles for k in 1:2]
    samples_im = [imag(ψ_out[k]).particles for k in 1:2]
    ψ_samples = [ComplexF64[samples_re[1][j] + im * samples_im[1][j],
                            samples_re[2][j] + im * samples_im[2][j]] for j in 1:N]
    stokes_samples = [stokes_from_jones(ψ_samples[j]) for j in 1:N]

    ψ_mean_unnorm = ComplexF64[
        pmean(real(ψ_out[1])) + im * pmean(imag(ψ_out[1])),
        pmean(real(ψ_out[2])) + im * pmean(imag(ψ_out[2])),
    ]
    ψ_mean   = ψ_mean_unnorm / norm(ψ_mean_unnorm)
    rep_mean = poincare_vector_representation(ψ_mean)

    _write_poincare_cloud_plot(
        rep_mean, stokes_samples, ΔT_samples;
        output = poincare_output,
        title  = title,
    )

    println("Wrote Poincaré cloud to: ", poincare_output)

    return (
        fiber          = result.fiber,
        ΔT_samples     = ΔT_samples,
        ψ_samples      = ψ_samples,
        stokes_samples = stokes_samples,
        plot_poincare  = poincare_output,
    )
end

"""
    demo_fiber_helix_mcm_twist(; output, fidelity, title)

Illustrate application of MCM twist to fiber segment with a helical twist at the center. MCM is used
to wobble the twist in the fiber before the helix to mimic vibration/air flow in the lab. 
"""
function demo_fiber_helix_mcm_twist(;
    output::AbstractString = joinpath(@__DIR__, "..", "output", "demo_fiber_helix_mcm_twist.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "demo_fiber_helix_mcm_twist()",
)
    PG = PathGeometry
    _nick(s) = PG.AbstractMeta[PG.Nickname(s)]
    spec = PG.PathSpecBuilder()
    PG.straight!(spec; length = 1, meta = _nick("lead-in"))
    # TODO: twist refactor —     PG.twist!(spec; s_start = 0, length = 1, rate = 1 ± 0.1) # MCM twist with mean 1 turn/unit length and ±0.1 uncertainty
    PG.helix!(spec; radius = 0.5, pitch = 0.05, turns = 4, axis_angle = 0.0, meta = _nick("helix"))
    PG.straight!(spec; length = 1, meta = _nick("lead-out"))
    path = PG.build(spec)
    plot_path = write_path_geometry_plot3d(
        path,
        path.spec.s_start,
        path.s_end;
        fidelity = fidelity,
        output = output,
        title = title,
    )
    return (; path, plot_path)
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
    (
        fn   = demo_fiber_path_modify,
        kwargs = NamedTuple(),
        desc = "Exercises `modify` on a shared inverted-U baseline path. Three HTML outputs " *
               "show how `MCMadd(:length, …)`, `MCMmul(:radius, …)`, and `MCMadd(:angle, …)` " *
               "perturb a single target segment (red) against the unmodified baseline (green)."
    ),
    (
        fn   = demo_seven_segment_mcm_temperature,
        kwargs = NamedTuple(),
        desc = "Seven-segment fiber (lead · bend · spacer · 1000 m middle straight w/ twist · " *
               "spacer · bend · lead). Middle segment carries `MCMadd(:T_K, Uniform(-5,5))` " *
               "via meta; `modify` rescales its length, then `propagate_fiber` lifts the " *
               "Particles ensemble through the Jones integrator. End-to-end MCM stack test."
    ),
    (
        fn   = demo_fiber_helix_mcm_twist,
        kwargs = NamedTuple(),
        desc = "Illustrate application of MCM twist to fiber segment with a helical twist at the " *
            "center. MCM is used to wobble the twist in the fiber before the helix to mimic" * 
            "vibration/air flow in the lab."
    )
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
