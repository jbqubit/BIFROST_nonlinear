# =====================================================================
# demo-path-geometry.jl
#
# Geometry-only demos used as a human-in-the-loop visual debug tool for
# `path-geometry*.jl`. Loads only the geometry layer (no fiber, no
# cross-section, no propagation). Each demo function authors a Subpath
# (or PathBuilt) using the new SubpathBuilder API, builds it, and writes
# an interactive Plotly HTML to `julia-port/output/`.
#
# Run all demos:
#
#     include("julia-port/demo-path-geometry.jl")
#     demo_path_geometry_all()
# =====================================================================

include("path-geometry.jl")
include("path-geometry-plot.jl")

const _OUTPUT_DIR = joinpath(@__DIR__, "..", "output")

_path_html(name::AbstractString) = joinpath(_OUTPUT_DIR, name)

# Convenience wrapper: SubpathBuilder built object → plot full domain.
function _plot_full(b::PathGeometry.SubpathBuilt; output, title, fidelity)
    s_hi = Float64(PathGeometry._qc_nominalize(PathGeometry.s_end(b)))
    return write_path_geometry_plot3d(b, 0.0, s_hi;
        output = output, title = title, fidelity = fidelity)
end

function _plot_full(p::PathGeometry.PathBuilt; output, title, fidelity)
    s_hi = PathGeometry.s_end(p)
    return write_path_geometry_plot3d(p, 0.0, s_hi;
        output = output, title = title, fidelity = fidelity)
end

# =====================================================================
# 2.1 Simple multi-segment Subpath (mirrors old demo_path_geometry).
# =====================================================================

function demo_path_geometry_simple(;
    output::AbstractString = _path_html("path-geometry.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Fiber path geometry: bends, catenary",
)
    PG = PathGeometry
    sb = PG.SubpathBuilder()
    PG.start!(sb)
    PG.straight!(sb; length = 0.10, meta = [PG.Nickname("Straight")])
    PG.bend!(sb; radius = 0.05, angle = π / 2, meta = [PG.Nickname("Bend")])
    PG.straight!(sb; length = 0.12, meta = [PG.Nickname("Straight")])
    PG.catenary!(sb; a = 0.03, length = 0.10, axis_angle = 0.0,
                 meta = [PG.Nickname("Catenary")])
    PG.bend!(sb; radius = 0.06, angle = π / 3, meta = [PG.Nickname("Bend")])
    PG.straight!(sb; length = 0.08, meta = [PG.Nickname("Straight")])
    # Seal at a point along the natural exit direction. The connector
    # collapses to a near-degenerate arc and contributes negligible length.
    sub = PG.Subpath(_seal_natural(sb))
    b = PG.build(sub)
    println("Arc length: ", PG.path_length(b))
    println("Writhe:     ", PG.writhe(b; n = 128))
    plot_path = _plot_full(b; output, title, fidelity)
    println("Wrote ", plot_path)
    return (; path = b, plot_path)
end

# Internal: build the SubpathBuilder once with a placeholder jumpto, peek
# at the natural exit pos/tangent, then re-author with that pos/tangent.
# This is a small price for not requiring each demo to compute exit state
# analytically. Returns a *new sealed* SubpathBuilder.
function _seal_natural(sb::PathGeometry.SubpathBuilder;
                       extra::Float64 = 0.0)
    PG = PathGeometry
    @assert isnothing(sb.jumpto_point) "_seal_natural: builder already sealed"
    # Trial seal at an arbitrary point. We just want to access sb.segments
    # placed-state via a temporary build. Use a far-away point so the
    # connector solver is well-conditioned.
    tmp = deepcopy(sb)
    PG.jumpto!(tmp; point = (1e9, 1e9, 1e9))
    b_tmp = PG.build(PG.Subpath(tmp))
    # Read pos and tangent at the END of the last interior segment
    # (= start of the terminal connector).
    s_end_interior = Float64(PG._qc_nominalize(b_tmp.jumpto_placed.s_offset_eff))
    if s_end_interior <= 0.0
        # No interior segments: use start state.
        natural_pos = collect(sb.start_point::NTuple{3, Float64})
        natural_tan = collect(sb.start_outgoing_tangent::NTuple{3, Float64})
    else
        natural_pos = collect(PG.position(b_tmp, s_end_interior))
        natural_tan = collect(PG.tangent(b_tmp, s_end_interior))
    end
    # Apply optional extra offset along the natural tangent.
    if extra != 0.0
        natural_pos = natural_pos .+ extra .* natural_tan
    end
    PG.jumpto!(sb;
        point = (natural_pos[1], natural_pos[2], natural_pos[3]),
        incoming_tangent = (natural_tan[1], natural_tan[2], natural_tan[3]),
    )
    return sb
end

# =====================================================================
# 2.2 Segment labels (mirrors old demo_path_geometry_segment_labels).
# =====================================================================

function demo_path_geometry_segment_labels(;
    output::AbstractString = _path_html("path-geometry-segment-labels.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Path geometry: segment nicknames",
)
    PG = PathGeometry
    sb = PG.SubpathBuilder()
    PG.start!(sb)
    PG.straight!(sb; length = 0.08, meta = [PG.Nickname("lead-in")])
    PG.bend!(sb; radius = 0.06, angle = π / 2, meta = [PG.Nickname("90° bend")])
    PG.straight!(sb; length = 0.06, meta = [PG.Nickname("spacer")])
    PG.catenary!(sb; a = 0.04, length = 0.08, axis_angle = 0.0,
                 meta = [PG.Nickname("sag")])
    PG.helix!(sb; radius = 0.025, pitch = 0.015, turns = 1.2,
              axis_angle = 0.0, meta = [PG.Nickname("twist section")])
    PG.straight!(sb; length = 0.06, meta = [PG.Nickname("lead-out")])
    b = PG.build(PG.Subpath(_seal_natural(sb)))
    println("Arc length: ", PG.path_length(b), " m")
    plot_path = _plot_full(b; output, title, fidelity)
    println("Wrote ", plot_path)
    return (; path = b, plot_path)
end

# =====================================================================
# 2.3 Helix demos for axis_angle = 0, π/3, 2π/3.
# =====================================================================

function _demo_helix(axis_angle::Float64; output, fidelity, title)
    PG = PathGeometry
    sb = PG.SubpathBuilder()
    PG.start!(sb)
    PG.straight!(sb; length = 0.05, meta = [PG.Nickname("Straight")])
    PG.helix!(sb; radius = 0.03, pitch = 0.02, turns = 2.0,
              axis_angle = axis_angle, meta = [PG.Nickname("Helix")])
    PG.straight!(sb; length = 0.05, meta = [PG.Nickname("Straight")])
    b = PG.build(PG.Subpath(_seal_natural(sb)))
    println("Helix axis_angle=", axis_angle, ": arc_length=",
            round(PG.path_length(b); digits = 4), " m")
    plot_path = _plot_full(b; output, title, fidelity)
    println("Wrote ", plot_path)
    return (; path = b, plot_path)
end

demo_path_geometry_helix_0(;
    output::AbstractString = _path_html("path-geometry-helix-0.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "HelixSegment: axis_angle = 0",
) = _demo_helix(0.0; output, fidelity, title)

demo_path_geometry_helix_pi_3(;
    output::AbstractString = _path_html("path-geometry-helix-pi-3.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "HelixSegment: axis_angle = π/3",
) = _demo_helix(π / 3; output, fidelity, title)

demo_path_geometry_helix_2pi_3(;
    output::AbstractString = _path_html("path-geometry-helix-2pi-3.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "HelixSegment: axis_angle = 2π/3",
) = _demo_helix(2π / 3; output, fidelity, title)

# =====================================================================
# 2.4 jumps_min_radius — paddle pattern realized as a PathBuilt.
#
# Old demo: alternating straight / jumpto / straight / jumpto / ...
# In the new architecture each `jumpto!` seals a Subpath, so the
# original 4-segment 4-jump pattern becomes 5 Subpaths.
# =====================================================================

function demo_path_geometry_jumps_min_radius(;
    output::AbstractString = _path_html("path-geometry-jumps-min-radius.html"),
    fidelity::Float64 = 4.0,
    title::AbstractString = "JumpBy/JumpTo paddle: PathBuilt of 5 Subpaths",
)
    PG = PathGeometry
    # Subpath 1: straight up to (0,0,1), seal with jumpto landing at (1,0,1)
    # with incoming tangent (0,0,-1) (heading down at landing). This gives
    # a transverse chord; min_bend_radius=0.4 keeps the connector smooth.
    sb1 = PG.SubpathBuilder()
    PG.start!(sb1)
    PG.straight!(sb1; length = 1.0, meta = [PG.Nickname("Sub1 straight")])
    PG.jumpto!(sb1; point = (1.0, 0.0, 1.0),
               incoming_tangent = (0.0, 0.0, -1.0),
               min_bend_radius = 0.4)

    # Subpath 2: starts at (1,0,1) heading -z, straight to (1,0,0), seals
    # to (2,0,0) with incoming tangent (0,0,1).
    sb2 = PG.SubpathBuilder()
    PG.start!(sb2; point = (1.0, 0.0, 1.0),
                  outgoing_tangent = (0.0, 0.0, -1.0))
    PG.straight!(sb2; length = 1.0, meta = [PG.Nickname("Sub2 straight")])
    PG.jumpto!(sb2; point = (2.0, 0.0, 0.0),
               incoming_tangent = (0.0, 0.0, 1.0),
               min_bend_radius = 0.1)

    # Subpath 3: starts at (2,0,0) heading +z, straight to (2,0,1),
    # seals to (3,0,1) with incoming tangent (0,0,-1).
    sb3 = PG.SubpathBuilder()
    PG.start!(sb3; point = (2.0, 0.0, 0.0),
                  outgoing_tangent = (0.0, 0.0, 1.0))
    PG.straight!(sb3; length = 1.0, meta = [PG.Nickname("Sub3 straight")])
    PG.jumpto!(sb3; point = (3.0, 0.0, 1.0),
               incoming_tangent = (0.0, 0.0, -1.0),
               min_bend_radius = 0.05)

    # Subpath 4: starts at (3,0,1) heading -z, straight + interior JumpBy.
    # JumpBy delta is in the local frame; after the straight, local +z is
    # global -z. So delta=(-1,0,0) (local) means -1 m along global +x.
    # We'll seal this Subpath at the interior JumpBy's natural endpoint.
    sb4 = PG.SubpathBuilder()
    PG.start!(sb4; point = (3.0, 0.0, 1.0),
                  outgoing_tangent = (0.0, 0.0, -1.0))
    PG.straight!(sb4; length = 1.0, meta = [PG.Nickname("Sub4 straight")])
    PG.jumpby!(sb4; delta = (-1.0, 0.0, 0.0),
               tangent = (0.0, 0.0, -1.0),
               min_bend_radius = 0.1,
               meta = [PG.Nickname("Sub4 JumpBy")])
    sb4 = _seal_natural(sb4)

    # Subpath 5: continues straight from sb4's natural endpoint.
    sb5 = PG.SubpathBuilder()
    sub4 = PG.Subpath(sb4)
    b4_tmp = PG.build(sub4)
    end_pos = collect(PG.end_point(b4_tmp))
    end_tan = collect(PG.end_tangent(b4_tmp))
    PG.start!(sb5; point = (end_pos[1], end_pos[2], end_pos[3]),
                  outgoing_tangent = (end_tan[1], end_tan[2], end_tan[3]))
    PG.straight!(sb5; length = 1.0, meta = [PG.Nickname("Sub5 straight")])
    sb5 = _seal_natural(sb5)

    p = PG.build([PG.Subpath(sb1), PG.Subpath(sb2), PG.Subpath(sb3),
                  sub4, PG.Subpath(sb5)])
    println("PathBuilt arc length: ", PG.path_length(p), " m")
    plot_path = _plot_full(p; output, title, fidelity)
    println("Wrote ", plot_path)
    return (; path = p, plot_path)
end

# =====================================================================
# 2.5 NEW: Multi-Subpath demo showing PathBuilt assembly.
# =====================================================================

function demo_path_geometry_pathbuilt(;
    output::AbstractString = _path_html("path-geometry-pathbuilt.html"),
    fidelity::Float64 = 2.0,
    title::AbstractString = "PathBuilt: three Subpaths (straight, bend, helix)",
)
    PG = PathGeometry

    # Subpath 1: straight, sealed at (0,0,0.2) with tangent +z.
    sb1 = PG.SubpathBuilder(meta = [PG.Nickname("Subpath 1: straight")])
    PG.start!(sb1)
    PG.straight!(sb1; length = 0.2, meta = [PG.Nickname("Straight")])
    PG.jumpto!(sb1; point = (0.0, 0.0, 0.2),
               incoming_tangent = (0.0, 0.0, 1.0))

    # Subpath 2: starts at (0,0,0.2) tangent +z, quarter bend (axis_angle=0),
    # so end position (R, 0, 0.2 + R) and end tangent +x.
    R2 = 0.05
    sb2 = PG.SubpathBuilder(meta = [PG.Nickname("Subpath 2: bend")])
    PG.start!(sb2; point = (0.0, 0.0, 0.2),
                  outgoing_tangent = (0.0, 0.0, 1.0))
    PG.bend!(sb2; radius = R2, angle = π / 2,
             meta = [PG.Nickname("90° bend")])
    PG.jumpto!(sb2; point = (R2, 0.0, 0.2 + R2),
               incoming_tangent = (1.0, 0.0, 0.0))

    # Subpath 3: starts at (R2, 0, 0.2 + R2) tangent +x, helix in transverse
    # plane (axis_angle=0). Seal at the helix's natural exit.
    sb3 = PG.SubpathBuilder(meta = [PG.Nickname("Subpath 3: helix")])
    PG.start!(sb3; point = (R2, 0.0, 0.2 + R2),
                  outgoing_tangent = (1.0, 0.0, 0.0))
    PG.helix!(sb3; radius = 0.025, pitch = 0.02, turns = 1.5,
              axis_angle = 0.0, meta = [PG.Nickname("Helix")])
    sb3 = _seal_natural(sb3)

    p = PG.build([PG.Subpath(sb1), PG.Subpath(sb2), PG.Subpath(sb3)])
    println("PathBuilt: ", length(p.subpaths), " subpaths")
    println("PathBuilt arc length: ", PG.path_length(p), " m")
    plot_path = _plot_full(p; output, title, fidelity)
    println("Wrote ", plot_path)
    return (; path = p, plot_path)
end

# =====================================================================
# 2.6 Aggregator + index HTML
# =====================================================================

# Registry of demos. Each entry: (group, fn, desc). Used by
# `demo_path_geometry_all` to run every demo and emit an index HTML that
# links to each output with a short description.
const _DEMO_PATH_GEOMETRY_INDEX = [
    (group = "Subpath",
     fn    = demo_path_geometry_simple,
     desc  = "Single Subpath: straight + bend + straight + catenary + bend + straight, sealed at the natural exit."),
    (group = "Subpath",
     fn    = demo_path_geometry_segment_labels,
     desc  = "Same shape as 'simple' but every segment carries a Nickname so labels render in the plot."),
    (group = "Subpath — helix",
     fn    = demo_path_geometry_helix_0,
     desc  = "HelixSegment with axis_angle = 0."),
    (group = "Subpath — helix",
     fn    = demo_path_geometry_helix_pi_3,
     desc  = "HelixSegment with axis_angle = π/3."),
    (group = "Subpath — helix",
     fn    = demo_path_geometry_helix_2pi_3,
     desc  = "HelixSegment with axis_angle = 2π/3."),
    (group = "PathBuilt",
     fn    = demo_path_geometry_jumps_min_radius,
     desc  = "Paddle pattern: 5 Subpaths joined at jumpto endpoints, each with its own min_bend_radius. Includes one interior JumpBy segment."),
    (group = "PathBuilt",
     fn    = demo_path_geometry_pathbuilt,
     desc  = "Three Subpaths (straight, quarter-bend, helix) stitched into a PathBuilt. Demonstrates the conformity check at each Subpath boundary."),
]

function demo_path_geometry_all(;
    index_output::AbstractString = _path_html("demo-path-geometry-index.html"),
)
    isdir(_OUTPUT_DIR) || mkpath(_OUTPUT_DIR)

    # entries: (group, title, path, desc)
    entries = Tuple{String, String, String, String}[]
    for d in _DEMO_PATH_GEOMETRY_INDEX
        println("[ demo ] $(nameof(d.fn))")
        result = d.fn()
        push!(entries, (d.group, basename(result.plot_path),
                        result.plot_path, d.desc))
    end

    # Group ordering preserves insertion order.
    seen_groups = String[]
    for (g, _, _, _) in entries
        g in seen_groups || push!(seen_groups, g)
    end

    open(index_output, "w") do io
        println(io, """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>BIFROST path-geometry-only demos</title>
  <style>
    body { font-family: sans-serif; max-width: 800px; margin: 2em auto; background: #111; color: #ddd; }
    h1   { font-size: 1.5em; border-bottom: 1px solid #444; padding-bottom: 0.3em; }
    h2   { font-size: 1.15em; margin-top: 1.8em; color: #4db87a; }
    ul   { padding-left: 1.2em; }
    li   { margin: 1em 0; }
    a    { font-weight: bold; color: #4db87a; }
    p.desc { margin: 0.3em 0 0 0; color: #999; font-size: 0.95em; }
    p.intro { color: #aaa; font-size: 0.95em; }
  </style>
</head>
<body>
  <h1>BIFROST path-geometry-only demos</h1>
  <p class="intro">
    These demos exercise <code>path-geometry*.jl</code> end-to-end without
    the fiber, cross-section, or propagation layers. Each demo authors a
    <code>Subpath</code> (or a <code>PathBuilt</code> of multiple Subpaths)
    via the <code>SubpathBuilder</code> API and writes an interactive
    Plotly HTML.
  </p>""")
        for g in seen_groups
            println(io, "  <h2>$(g)</h2>")
            println(io, "  <ul>")
            for (eg, title, path, desc) in entries
                eg == g || continue
                href = basename(path)
                println(io, "    <li>")
                println(io, "      <a href=\"$(href)\">$(title)</a>")
                println(io, "      <p class=\"desc\">$(desc)</p>")
                println(io, "    </li>")
            end
            println(io, "  </ul>")
        end
        println(io, """</body>
</html>""")
    end

    println()
    println("Wrote demo index to: ", index_output)
    return index_output
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_path_geometry_all()
end
