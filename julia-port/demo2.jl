# =====================================================================
# demo2.jl — JumpBy / JumpTo visual experiments
# =====================================================================
#
# Companion file to demo.jl, dedicated to illustrating `JumpBy` and
# `JumpTo`. Two flavours:
#
#   * 2D scenes (inline SVG, no external library):
#     demo_fiber_path_jumps_2d() and demo_fiber_path_jumps_min_radius_2d()
#     — the path centerline projected onto the x–z plane. Useful as a
#     quick-look that needs no JavaScript.
#
#   * 3D scenes (Plotly): demo_fiber_path_jumps() — the same scenes,
#     with the resolved QuinticConnector drawn in red and the surrounding
#     fixed segments drawn in green.
#
# `demo2_all()` runs both groups (2D first, 3D second) and writes
# `output/index2.html` with the two groups under separate headings.
#
# This file expects to be `include`d after demo.jl is in scope (it
# reuses `_sample_segment_xyz` and the path-builder API).

if !isdefined(Main, :_sample_segment_xyz)
    include(joinpath(@__DIR__, "demo.jl"))
end

# ---------------------------------------------------------------------
# Shared builders (used by both 3D and 2D demos)
# ---------------------------------------------------------------------

function _build_jumpby_variant(L::Float64, R::Float64, target_idx::Int,
                                jump_kwargs::NamedTuple;
                                priors::Vector = [(:straight, (length = L,))],
                                tails::Vector  = [(:straight, (length = L,))])
    spec = PathSpecBuilder()
    for (kind, kw) in priors
        kind === :straight ? straight!(spec; kw...) :
        kind === :bend     ? bend!(spec; kw...)     :
        kind === :helix    ? helix!(spec; kw...)    :
        error("unknown prior kind: $kind")
    end
    jumpby!(spec; jump_kwargs...)
    for (kind, kw) in tails
        kind === :straight ? straight!(spec; kw...) :
        kind === :bend     ? bend!(spec; kw...)     :
        kind === :helix    ? helix!(spec; kw...)    :
        error("unknown tail kind: $kind")
    end
    return (build(spec), target_idx)
end

function _build_jumpto_variant(L::Float64, R::Float64, target_idx::Int,
                                jump_kwargs::NamedTuple;
                                priors::Vector = [(:straight, (length = L,))],
                                tails::Vector  = [(:straight, (length = L,))])
    spec = PathSpecBuilder()
    for (kind, kw) in priors
        kind === :straight ? straight!(spec; kw...) :
        kind === :bend     ? bend!(spec; kw...)     :
        kind === :helix    ? helix!(spec; kw...)    :
        error("unknown prior kind: $kind")
    end
    jumpto!(spec; jump_kwargs...)
    for (kind, kw) in tails
        kind === :straight ? straight!(spec; kw...) :
        kind === :bend     ? bend!(spec; kw...)     :
        kind === :helix    ? helix!(spec; kw...)    :
        error("unknown tail kind: $kind")
    end
    return (build(spec), target_idx)
end

function _build_jump_composite(segments::Vector)
    spec = PathSpecBuilder()
    for (kind, kw) in segments
        if kind === :straight; straight!(spec; kw...)
        elseif kind === :bend; bend!(spec; kw...)
        elseif kind === :helix; helix!(spec; kw...)
        elseif kind === :jumpby; jumpby!(spec; kw...)
        elseif kind === :jumpto; jumpto!(spec; kw...)
        else error("unknown segment kind: $kind") end
    end
    return build(spec)
end

# =====================================================================
# 2D GROUP — inline-SVG scenes
# =====================================================================
#
# Renderer + scenes for the 2D demos (x–z projection). No JavaScript,
# no external dependencies. Scene functions live below the renderer.
# ---------------------------------------------------------------------
# 2D inline-SVG renderer
# ---------------------------------------------------------------------
#
# All paths in this demo set lie in the y = 0 plane (we only ever use
# axis_angle = 0 and never set y-components in delta/destination), so
# the natural 2D projection is (x, z). The SVG is emitted directly —
# no JS, no external dependency.

# Sample every placed segment of `path`. Returns a Vector of
# (xs, zs, is_red) tuples — one per segment.
function _sample_path_2d(path, red_indices::Vector{Int})
    rows = Tuple{Vector{Float64}, Vector{Float64}, Bool}[]
    for i in 1:length(path.placed_segments)
        s = _sample_segment_xyz(path, i)
        push!(rows, (s.x, s.z, i in red_indices))
    end
    return rows
end

# Render a single 2D scene to `output` as a standalone SVG-in-HTML file.
# `variants` is a list of (label, build_fn) where build_fn returns
# `(path, red_indices)`. Variants are offset along x by `variant_spacing`.
# Setting variant_spacing=0 overlays them.
function _jump_row_svg(output::AbstractString, title::AbstractString;
                       variants::Vector,
                       variant_spacing::Float64 = 2.5,
                       width::Int = 1100, height::Int = 520,
                       margin::Int = 60)
    # Collect all sampled curves across all variants, applying x-offset.
    all_curves = Vector{Tuple{Vector{Float64}, Vector{Float64}, Bool, Int}}()
    starts     = Vector{Tuple{Float64, Float64}}()
    labels     = Vector{Tuple{Float64, Float64, String}}()

    xmin =  Inf; xmax = -Inf
    zmin =  Inf; zmax = -Inf

    for (k, (label, build_fn)) in enumerate(variants)
        path, red_spec = build_fn()
        red_indices = red_spec isa Integer ? Int[red_spec] :
                      Int.(collect(red_spec))
        dx   = (k - 1) * variant_spacing
        rows = _sample_path_2d(path, red_indices)
        for (xs, zs, is_red) in rows
            xs2 = xs .+ dx
            push!(all_curves, (xs2, zs, is_red, k))
            xmin = min(xmin, minimum(xs2)); xmax = max(xmax, maximum(xs2))
            zmin = min(zmin, minimum(zs));  zmax = max(zmax, maximum(zs))
        end
        p0 = position(path, path.spec.s_start)
        push!(starts, (Float64(p0[1]) + dx, Float64(p0[3])))
        # Label sits just above each variant.
        x_label = dx
        push!(labels, (x_label, zmax, label))
    end

    # Padding so the curves don't touch the axes.
    pad_x = 0.05 * max(xmax - xmin, 1e-6)
    pad_z = 0.10 * max(zmax - zmin, 1e-6)
    xmin -= pad_x; xmax += pad_x
    zmin -= pad_z; zmax += pad_z + 0.15  # extra room for labels

    # Map data (x, z) to SVG (px, py). SVG y grows downward, so flip z.
    plot_w = width  - 2*margin
    plot_h = height - 2*margin
    sx = plot_w / (xmax - xmin)
    sz = plot_h / (zmax - zmin)
    s  = min(sx, sz)  # equal aspect
    # Center the data in the plot region.
    ox = margin + 0.5 * (plot_w - s * (xmax - xmin)) - s * xmin
    oz = margin + 0.5 * (plot_h - s * (zmax - zmin)) + s * zmax
    px(x) = ox + s * x
    pz(z) = oz - s * z

    # Build the SVG body.
    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$(width)" height="$(height)" viewBox="0 0 $(width) $(height)">""")
    println(io, """<rect width="100%" height="100%" fill="#111"/>""")

    # Axis box.
    bx0 = px(xmin); bz0 = pz(zmax); bw = s*(xmax-xmin); bh = s*(zmax-zmin)
    println(io, """<rect x="$(round(bx0;digits=2))" y="$(round(bz0;digits=2))" width="$(round(bw;digits=2))" height="$(round(bh;digits=2))" fill="none" stroke="#333" stroke-width="1"/>""")

    # Light grid at integer coordinates (data units of metres).
    for ix in ceil(Int, xmin):floor(Int, xmax)
        x = px(Float64(ix))
        println(io, """<line x1="$(round(x;digits=2))" y1="$(round(bz0;digits=2))" x2="$(round(x;digits=2))" y2="$(round(bz0+bh;digits=2))" stroke="#222" stroke-width="0.5"/>""")
    end
    for iz in ceil(Int, zmin):floor(Int, zmax)
        z = pz(Float64(iz))
        println(io, """<line x1="$(round(bx0;digits=2))" y1="$(round(z;digits=2))" x2="$(round(bx0+bw;digits=2))" y2="$(round(z;digits=2))" stroke="#222" stroke-width="0.5"/>""")
    end

    # Axis labels.
    println(io, """<text x="$(round(bx0+bw/2;digits=2))" y="$(height-15)" fill="#aaa" font-family="sans-serif" font-size="13" text-anchor="middle">x (m)</text>""")
    println(io, """<text x="18" y="$(round(bz0+bh/2;digits=2))" fill="#aaa" font-family="sans-serif" font-size="13" text-anchor="middle" transform="rotate(-90 18,$(round(bz0+bh/2;digits=2)))">z (m)</text>""")
    println(io, """<text x="$(round(width/2;digits=2))" y="22" fill="#eee" font-family="sans-serif" font-size="15" text-anchor="middle">$(title)</text>""")

    # Curves.
    for (xs, zs, is_red, _k) in all_curves
        color = is_red ? "#e55" : "#6c6"
        # Build a polyline points string.
        parts = String[]
        for i in eachindex(xs)
            push!(parts, string(round(px(xs[i]); digits=2), ",", round(pz(zs[i]); digits=2)))
        end
        println(io, """<polyline points="$(join(parts, ' '))" fill="none" stroke="$(color)" stroke-width="2.5" stroke-linejoin="round"/>""")
    end

    # Start markers (green dots).
    for (sxv, szv) in starts
        println(io, """<circle cx="$(round(px(sxv);digits=2))" cy="$(round(pz(szv);digits=2))" r="4" fill="#6c6" stroke="#9d9" stroke-width="1"/>""")
    end

    # White circles at the end of every atomic segment (segment boundaries).
    for (xs, zs, _is_red, _k) in all_curves
        ex = xs[end]; ez = zs[end]
        println(io, """<circle cx="$(round(px(ex);digits=2))" cy="$(round(pz(ez);digits=2))" r="3" fill="#fff" stroke="#000" stroke-width="0.6"/>""")
    end

    # Variant labels (top of each column).
    for (lx, lz, txt) in labels
        println(io, """<text x="$(round(px(lx);digits=2))" y="$(round(pz(lz)-6;digits=2))" fill="#ddd" font-family="sans-serif" font-size="12" text-anchor="middle">$(txt)</text>""")
    end

    println(io, "</svg>")
    svg_body = String(take!(io))

    html = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>$(title)</title>
<style>html,body{margin:0;padding:0;background:#111;color:#eee;font-family:sans-serif;}
.wrap{display:flex;align-items:center;justify-content:center;min-height:100vh;}
svg{max-width:100%;height:auto;}</style>
</head><body><div class="wrap">$(svg_body)</div></body></html>
"""
    open(output, "w") do io
        write(io, html)
    end
    return output
end

# Helper: build a path from an explicit list of builder calls (one per
# atomic segment) and return it in the (path, red_indices) shape that
# `_jump_row_svg` expects. Use it like:
#
#     _path() do spec
#         straight!(spec; length = 1.0)
#         jumpby!(spec;   delta  = (0.2, 0.0, 0.5))
#         straight!(spec; length = 1.0)
#     end
#
# `red` is optional and lists which atomic segments to colour red; with
# the default `Int[]` everything draws in the single non-red colour.
function _path(f::Function; red::Vector{Int} = Int[])
    spec = PathSpecBuilder()
    f(spec)
    return (build(spec), red)
end

# Build a path one segment at a time, stopping at the first segment whose
# `build()` throws. Returns `(partial_path_or_nothing, n_built, error_or_nothing)`.
# Used by `min_bend_radius` demos to display a partial path with the
# infeasible segment missing.
#
# `f(spec)` should call atomic builders (straight!, bend!, jumpby!, etc.).
# Each call must append exactly one segment to `spec.segments`; we trial-
# build after every append and roll back the last segment if it fails.
function _build_with_failure(f::Function)
    probe_spec = PathSpecBuilder()
    f(probe_spec)
    declared = copy(probe_spec.segments)

    spec = PathSpecBuilder()
    last_good = nothing
    for (i, seg) in enumerate(declared)
        push!(spec.segments, seg)
        try
            last_good = build(spec)
        catch e
            pop!(spec.segments)  # discard the failing segment
            partial = isempty(spec.segments) ? nothing : build(spec)
            return (partial, i - 1, e)
        end
    end
    return (last_good, length(declared), nothing)
end

"""
    demo_fiber_path_jumps_2d(; output_dir = …)

2D (x, z projection) demos for `JumpBy` / `JumpTo`, rendered as
inline SVG with no external library. Each variant is constructed
inline with explicit calls to `straight!`, `bend!`, `jumpby!`,
`jumpto!`, etc., so the atomic segment layout and numerical
arguments are visible in the call site.

White circles mark the end of every atomic segment.
"""
function demo_fiber_path_jumps_2d(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "output"),
)
    s1 = _jump_row_svg(
        joinpath(output_dir, "jumpby-2d-delta-transverse.html"),
        "2D — JumpBy(delta=(d, 0, 0.5)): transverse offset sweep";
        variants = [
            ("d = 0.0", () -> _path() do spec
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta  = (0.0, 0.0, 0.5), tangent=(0, 0.0, 1.0))
                straight!(spec; length = 1.0)
            end),
            ("d = 0.2", () -> _path() do spec
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta  = (0.2, 0.0, 0.5), tangent=(0, 0.0, 1.0))
                straight!(spec; length = 1.0)
            end),
            ("d = 0.5", () -> _path() do spec
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta  = (0.5, 0.0, 0.5), tangent=(0, 0.0, 1.0))
                straight!(spec; length = 1.0)
            end),
            ("d = 0.8", () -> _path() do spec
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta  = (0.8, 0.0, 0.5), tangent=(0, 0.0, 1.0), min_bend_radius=0.3)
                straight!(spec; length = 1.0)
            end),
        ],
        variant_spacing = 2.0,
    )

    s2 = _jump_row_svg(
        joinpath(output_dir, "jumpby-2d-tangent-out.html"),
        "2D — JumpBy: outgoing tangent (local frame)";
        variants = [
            ("t = (+1,0,1)/√2", () -> _path() do spec
                straight!(spec; length  = 1.0)
                jumpby!(spec;   delta   = (0.4, 0.0, 0.4),
                                tangent = (1/sqrt(2), 0.0, 1/sqrt(2)))
                straight!(spec; length  = 1.0)
            end),
            ("t = (0,0,1)", () -> _path() do spec
                straight!(spec; length  = 1.0)
                jumpby!(spec;   delta   = (0.4, 0.0, 0.4),
                                tangent = (0.0, 0.0, 1.0))
                straight!(spec; length  = 1.0)
            end),
            ("t = (-1,0,1)/√2", () -> _path() do spec
                straight!(spec; length  = 1.0)
                jumpby!(spec;   delta   = (0.4, 0.0, 0.4),
                                tangent = (-1/sqrt(2), 0.0, 1/sqrt(2)))
                straight!(spec; length  = 1.0)
            end),
        ],
        variant_spacing = 2.0,
    )

    s3 = _jump_row_svg(
        joinpath(output_dir, "jumpto-2d-destination.html"),
        "2D — JumpTo(destination=(x, 0, 1.5)): transverse sweep";
        variants = [
            ("x = 0.0", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.0, 0.0, 1.5))
                straight!(spec; length      = 1.0)
            end),
            ("x = 0.3", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.3, 0.0, 1.5))
                straight!(spec; length      = 1.0)
            end),
            ("x = 0.6", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.6, 0.0, 1.5))
                straight!(spec; length      = 1.0)
            end),
            ("x = 1.0", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (1.0, 0.0, 1.5))
                straight!(spec; length      = 1.0)
            end),
        ],
        variant_spacing = 2.5,
    )

    s4 = _jump_row_svg(
        joinpath(output_dir, "jumpto-2d-tangent-global.html"),
        "2D — JumpTo: outgoing tangent (GLOBAL frame)";
        variants = [
            ("t = (+1,0,0)", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.5, 0.0, 1.5),
                                tangent     = (1.0, 0.0, 0.0))
                straight!(spec; length      = 1.0)
            end),
            ("t = (0,0,+1)", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.5, 0.0, 1.5),
                                tangent     = (0.0, 0.0, 1.0))
                straight!(spec; length      = 1.0)
            end),
            ("t = (-1,0,1)/√2", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.5, 0.0, 1.5),
                                tangent     = (-1/sqrt(2), 0.0, 1/sqrt(2)))
                straight!(spec; length      = 1.0)
            end),
        ],
        variant_spacing = 2.5,
    )

    s5 = _jump_row_svg(
        joinpath(output_dir, "jumpto-2d-routing.html"),
        "2D — JumpTo routing: T1/T2/T3 composite (anti-parallel tangents)";
        variants = [
            ("composite", () -> _path(red = [2, 4, 6]) do spec
                straight!(spec; length          = 1.0)
                jumpto!(spec;   destination     = (1.0, 0.0, 1.0),
                                tangent         = (0.0, 0.0, -1.0),
                                min_bend_radius = 0.1)
                straight!(spec; length          = 1.0)
                jumpto!(spec;   destination     = (2.0, 0.0, 0.0),
                                tangent         = (0.0, 0.0, 1.0),
                                min_bend_radius = 0.1)
                straight!(spec; length          = 1.0)
                jumpto!(spec;   destination     = (3.0, 0.0, 1.0),
                                tangent         = (0.0, 0.0, -1.0),
                                min_bend_radius = 0.1)
            end),
        ],
        variant_spacing = 0.0,
        height = 420,
    )

    return (
        svg_jumpby_delta_transverse = s1,
        svg_jumpby_tangent_out      = s2,
        svg_jumpto_destination      = s3,
        svg_jumpto_tangent_global   = s4,
        svg_jumpto_routing          = s5,
    )
end

"""
    demo_fiber_path_jumps_min_radius_2d(; output_dir = …)

Sweep `min_bend_radius` on a `JumpTo` whose chord is transverse and whose
outgoing tangent is anti-parallel to the incoming one — the canonical
infeasibility case from the unit tests (threshold ≈ 0.5 m).

Small values succeed; values past the threshold intentionally fail. The
`ArgumentError` is trapped per-variant and the SVG renders the partial
path with the failed jump (and any subsequent segments) omitted; the
variant label shows `(infeasible)`.
"""
function demo_fiber_path_jumps_min_radius_2d(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "output"),
)
    # Build one variant for the given min_bend_radius value. The geometry
    # is fixed (`straight · JumpTo · straight`); only `mbr` varies. If the
    # build fails we render whatever was successfully built.
    function variant(mbr::Float64)
        partial, n_built, err = _build_with_failure() do spec
            straight!(spec; length = 1.0)
            jumpto!(spec;   destination     = (1.0, 0.0, 1.0),
                            tangent         = (0.0, 0.0, -1.0),
                            min_bend_radius = mbr)
            straight!(spec; length = 1.0)
        end
        # Always mark the JumpTo (segment 2) red when present so the
        # connector is visually distinguishable from the straights.
        red = (partial !== nothing && length(partial.placed_segments) >= 2) ?
              [2] : Int[]
        label = err === nothing ?
                "mbr = $(mbr)" :
                "mbr = $(mbr) (infeasible — $(n_built)/3 built)"
        return (label, partial, red)
    end

    # Five values: three feasible, two infeasible. Threshold is 0.5 m.
    raw = [variant(m) for m in (0.10, 0.30, 0.49, 0.51, 0.70)]

    # Drop any variant that produced no path at all (shouldn't happen here,
    # since the lead-in straight always succeeds, but guard anyway).
    raw = filter(t -> t[2] !== nothing, raw)

    variants = [(label, () -> (path, red)) for (label, path, red) in raw]

    return (svg_jumpto_min_radius = _jump_row_svg(
        joinpath(output_dir, "jumpto-2d-min-radius.html"),
        "2D — JumpTo min_bend_radius sweep (transverse chord, anti-parallel tangents; threshold ≈ 0.5 m)";
        variants = variants,
        variant_spacing = 2.0,
        height = 460,
    ),)
end

# =====================================================================
# 3D GROUP — Plotly scenes
# =====================================================================
#
# Renderer + scenes for the 3D demos. Output: standalone HTML files
# loading plotly.js from CDN. Scene functions live below the renderer.
# ---------------------------------------------------------------------
# 3D Plotly renderer
# ---------------------------------------------------------------------

function _jump_row_html(output::AbstractString, title::AbstractString;
                        variants::Vector,
                        variant_spacing::Float64 = 2.5,
                        z_label_offset::Float64 = 1.6)
    js_num(x)   = isnan(x) ? "NaN" : string(Float64(x))
    js_arr(xs)  = "[" * join(js_num.(xs), ",") * "]"
    js_strarr(xs) = "[" * join(("\"" * replace(string(x), "\"" => "\\\"") * "\"" for x in xs), ",") * "]"

    trace_strs  = String[]
    label_xs    = Float64[]; label_ys = Float64[]; label_zs = Float64[]
    label_texts = String[]
    start_xs    = Float64[]; start_ys = Float64[]; start_zs = Float64[]

    for (k, (label, build_fn)) in enumerate(variants)
        path, red_spec = build_fn()
        red_indices = red_spec isa Integer ? Int[red_spec] : Int.(collect(red_spec))
        dx = (k - 1) * variant_spacing
        n_segs = length(path.placed_segments)
        for i in 1:n_segs
            s     = _sample_segment_xyz(path, i)
            xs    = s.x .+ dx
            color = i in red_indices ? "#e55" : "#6c6"
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
        p0 = position(path, path.spec.s_start)
        push!(start_xs, Float64(p0[1]) + dx)
        push!(start_ys, Float64(p0[2]))
        push!(start_zs, Float64(p0[3]))

        push!(label_xs, dx)
        push!(label_ys, 0.0)
        push!(label_zs, z_label_offset)
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
    demo_fiber_path_jumps(; output_dir = …)

3D Plotly demos for `JumpBy` / `JumpTo`. The connector segment is
drawn in red; surrounding fixed segments are green. Variants are
offset along *x* for side-by-side comparison.
"""
function demo_fiber_path_jumps(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "output"),
)
    L = 1.0
    R = 0.5

    jb(kw) = () -> _build_jumpby_variant(L, R, 2, kw)

    row1 = _jump_row_html(
        joinpath(output_dir, "jumpby-delta-axial.html"),
        "JumpBy(delta=(0,0,d)) — axial sweep";
        variants = [
            ("d = 0.3",  jb((delta = (0.0, 0.0, 0.3),))),
            ("d = 0.6",  jb((delta = (0.0, 0.0, 0.6),))),
            ("d = 1.0",  jb((delta = (0.0, 0.0, 1.0),))),
        ],
        z_label_offset = 2.2,
    )

    row2 = _jump_row_html(
        joinpath(output_dir, "jumpby-delta-transverse.html"),
        "JumpBy(delta=(d,0,0.5)) — transverse sweep";
        variants = [
            ("d = 0.0",  jb((delta = (0.0, 0.0, 0.5),))),
            ("d = 0.2",  jb((delta = (0.2, 0.0, 0.5),))),
            ("d = 0.5",  jb((delta = (0.5, 0.0, 0.5),))),
        ],
        z_label_offset = 1.8,
    )

    row3 = _jump_row_html(
        joinpath(output_dir, "jumpby-tangent-out.html"),
        "JumpBy — outgoing tangent (local frame)";
        variants = [
            ("tangent = (+1,0,1)/√2", jb((delta = (0.4, 0.0, 0.4),
                                            tangent = (1/sqrt(2), 0.0, 1/sqrt(2))))),
            ("tangent = (0,0,1)",     jb((delta = (0.4, 0.0, 0.4),
                                            tangent = (0.0, 0.0, 1.0)))),
            ("tangent = (-1,0,1)/√2", jb((delta = (0.4, 0.0, 0.4),
                                            tangent = (-1/sqrt(2), 0.0, 1/sqrt(2))))),
        ],
        z_label_offset = 1.8,
    )

    row4 = _jump_row_html(
        joinpath(output_dir, "jumpby-curvature-out.html"),
        "JumpBy — outgoing curvature (G2 knob, local frame)";
        variants = [
            ("κ_out = 0",         jb((delta = (0.5, 0.0, 0.5),
                                       tangent = (1.0, 0.0, 1.0) ./ sqrt(2),
                                       curvature_out = (0.0, 0.0, 0.0)))),
            ("κ_out = (0,+2,0)",  jb((delta = (0.5, 0.0, 0.5),
                                       tangent = (1.0, 0.0, 1.0) ./ sqrt(2),
                                       curvature_out = (0.0,  2.0, 0.0)))),
            ("κ_out = (0,-2,0)",  jb((delta = (0.5, 0.0, 0.5),
                                       tangent = (1.0, 0.0, 1.0) ./ sqrt(2),
                                       curvature_out = (0.0, -2.0, 0.0)))),
        ],
        z_label_offset = 1.8,
    )

    jt(kw) = () -> _build_jumpto_variant(L, R, 2, kw)

    row5 = _jump_row_html(
        joinpath(output_dir, "jumpto-destination.html"),
        "JumpTo(destination=(0,0,z)) — axial sweep";
        variants = [
            ("z = 1.3",  jt((destination = (0.0, 0.0, 1.3),))),
            ("z = 1.6",  jt((destination = (0.0, 0.0, 1.6),))),
            ("z = 2.0",  jt((destination = (0.0, 0.0, 2.0),))),
        ],
        z_label_offset = 2.6,
    )

    row6 = _jump_row_html(
        joinpath(output_dir, "jumpto-destination-transverse.html"),
        "JumpTo(destination=(x,0,1.5)) — transverse sweep";
        variants = [
            ("x = 0.0",  jt((destination = (0.0, 0.0, 1.5),))),
            ("x = 0.3",  jt((destination = (0.3, 0.0, 1.5),))),
            ("x = 0.6",  jt((destination = (0.6, 0.0, 1.5),))),
        ],
        z_label_offset = 2.4,
    )

    row7 = _jump_row_html(
        joinpath(output_dir, "jumpto-tangent-global.html"),
        "JumpTo — outgoing tangent (GLOBAL frame)";
        variants = [
            ("tangent = (+1,0,0)", jt((destination = (0.5, 0.0, 1.5),
                                         tangent = (1.0, 0.0, 0.0)))),
            ("tangent = (0,0,+1)", jt((destination = (0.5, 0.0, 1.5),
                                         tangent = (0.0, 0.0, 1.0)))),
            ("tangent = (-1,0,1)/√2", jt((destination = (0.5, 0.0, 1.5),
                                         tangent = (-1/sqrt(2), 0.0, 1/sqrt(2))))),
        ],
        z_label_offset = 2.4,
    )

    row8 = _jump_row_html(
        joinpath(output_dir, "jumpto-curvature-global.html"),
        "JumpTo — outgoing curvature (GLOBAL frame)";
        variants = [
            ("κ_out = 0",          jt((destination = (0.5, 0.0, 1.5),
                                        tangent = (0.0, 0.0, 1.0),
                                        curvature_out = (0.0, 0.0, 0.0)))),
            ("κ_out = (10,0,0)",   jt((destination = (0.5, 0.0, 1.5),
                                        tangent = (0.0, 0.0, 1.0),
                                        curvature_out = (10.0, 0.0, 0.0)))),
            ("κ_out = (-10,0,0)",  jt((destination = (0.5, 0.0, 1.5),
                                        tangent = (0.0, 0.0, 1.0),
                                        curvature_out = (-10.0, 0.0, 0.0)))),
        ],
        variant_spacing = 1.5,
        z_label_offset = 2.4,
    )

    bend_priors = [(:straight, (length = L,)), (:bend, (radius = R, angle = π/2, axis_angle = 0.0))]
    jb_post_bend(kw) = () -> _build_jumpby_variant(L, R, 3, kw;
                                                    priors = bend_priors,
                                                    tails  = [(:straight, (length = L,))])

    row9 = _jump_row_html(
        joinpath(output_dir, "jumpby-after-bend.html"),
        "JumpBy after 90° bend — delta in ROTATED local frame";
        variants = [
            ("delta = (0,0,0.5)",   jb_post_bend((delta = (0.0, 0.0, 0.5),))),
            ("delta = (0.3,0,0.5)", jb_post_bend((delta = (0.3, 0.0, 0.5),))),
            ("delta = (-0.3,0,0.5)",jb_post_bend((delta = (-0.3, 0.0, 0.5),))),
        ],
        variant_spacing = 3.5,
        z_label_offset = 2.0,
    )

    jt_post_bend(kw) = () -> _build_jumpto_variant(L, R, 3, kw;
                                                    priors = bend_priors,
                                                    tails  = [(:straight, (length = L,))])

    row10 = _jump_row_html(
        joinpath(output_dir, "jumpto-after-bend.html"),
        "JumpTo after 90° bend — destination in GLOBAL frame";
        variants = [
            ("dest = (1.0, 0, 1.5)", jt_post_bend((destination = (1.0, 0.0, 1.5),))),
            ("dest = (1.3, 0, 1.5)", jt_post_bend((destination = (1.3, 0.0, 1.5),))),
            ("dest = (0.5, 0, 2.0)", jt_post_bend((destination = (0.5, 0.0, 2.0),))),
        ],
        variant_spacing = 3.5,
        z_label_offset = 2.4,
    )

    function jb_g2(R_b)
        () -> _build_jumpby_variant(L, R, 2, (delta = (0.3, 0.0, 0.3),);
                                     priors = [(:bend, (radius = R_b, angle = π/4, axis_angle = 0.0))],
                                     tails  = [(:straight, (length = L,))])
    end

    row11 = _jump_row_html(
        joinpath(output_dir, "jumpby-g2-inheritance.html"),
        "JumpBy after bend — G2 inheritance of incoming κ from prior bend";
        variants = [
            ("R_bend = 0.30 (κ_in ≈ 3.33)", jb_g2(0.30)),
            ("R_bend = 0.50 (κ_in = 2.00)", jb_g2(0.50)),
            ("R_bend = 0.80 (κ_in = 1.25)", jb_g2(0.80)),
        ],
        variant_spacing = 3.0,
        z_label_offset = 1.6,
    )

    function build_routing()
        path = _build_jump_composite([
            (:straight, (length = 1.0,)),
            (:jumpto,   (destination = (1.0, 0.0, 1.0),
                         tangent = (0.0, 0.0, -1.0),
                         min_bend_radius = 0.1)),
            (:straight, (length = 1.0,)),
            (:jumpto,   (destination = (2.0, 0.0, 0.0),
                         tangent = (0.0, 0.0, 1.0),
                         min_bend_radius = 0.1)),
            (:straight, (length = 1.0,)),
            (:jumpto,   (destination = (3.0, 0.0, 1.0),
                         tangent = (0.0, 0.0, -1.0),
                         min_bend_radius = 0.1)),
        ])
        return (path, [2, 4, 6])
    end

    row12 = _jump_row_html(
        joinpath(output_dir, "jumpto-routing.html"),
        "JumpTo routing — T1/T2/T3 composite (anti-parallel tangents)";
        variants = [("composite", build_routing)],
        variant_spacing = 0.0,
        z_label_offset = 1.6,
    )

    return (
        plot_jumpby_delta_axial      = row1,
        plot_jumpby_delta_transverse = row2,
        plot_jumpby_tangent_out      = row3,
        plot_jumpby_curvature_out    = row4,
        plot_jumpto_destination      = row5,
        plot_jumpto_dest_transverse  = row6,
        plot_jumpto_tangent_global   = row7,
        plot_jumpto_curvature_global = row8,
        plot_jumpby_after_bend       = row9,
        plot_jumpto_after_bend       = row10,
        plot_jumpby_g2_inheritance   = row11,
        plot_jumpto_routing          = row12,
    )
end

# ---------------------------------------------------------------------
# Index page (index2.html)
# ---------------------------------------------------------------------

const DEMO2_INDEX = [
    (
        group = "2D",
        fn   = demo_fiber_path_jumps_2d,
        kwargs = (;),
        desc = "2D (x–z projection) of the same JumpBy / JumpTo experiments, rendered as inline SVG (no JavaScript). Faster to scan than the 3D scenes.",
    ),
    (
        group = "2D",
        fn   = demo_fiber_path_jumps_min_radius_2d,
        kwargs = (;),
        desc = "2D sweep of the JumpTo min_bend_radius parameter through its feasibility threshold (~0.5 m). Infeasible variants are trapped at build time and rendered with the missing segment omitted; the variant label is annotated with '(infeasible)'.",
    ),
    (
        group = "3D",
        fn   = demo_fiber_path_jumps,
        kwargs = (;),
        desc = "3D Plotly visualizations of JumpBy / JumpTo. Connector segment in red; fixed segments in green. Each scene sweeps one parameter (delta, destination, tangent_out, curvature_out) with variants offset along x for side-by-side comparison.",
    ),
]

"""
    demo2_all(; index_output)

Run every demo in `DEMO2_INDEX` and write `index2.html` linking each
output file with a short description.
"""
function demo2_all(; index_output::AbstractString = joinpath(@__DIR__, "..", "output", "index2.html"))
    # Each entry: (group, title, path, desc).
    entries = Tuple{String, String, String, String}[]

    for d in DEMO2_INDEX
        result = d.fn(; d.kwargs...)
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
            push!(entries, (d.group, basename(p), p, d.desc))
        end
    end

    open(index_output, "w") do io
        println(io, """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>BIFROST JumpBy / JumpTo demos</title>
  <style>
    body { font-family: sans-serif; max-width: 800px; margin: 2em auto; color: #222; }
    h1   { font-size: 1.5em; border-bottom: 1px solid #ccc; padding-bottom: 0.3em; }
    h2   { font-size: 1.15em; margin-top: 1.8em; color: #1a6; }
    ul   { padding-left: 1.2em; }
    li   { margin: 1em 0; }
    a    { font-weight: bold; color: #1a6; }
    p.desc { margin: 0.3em 0 0 0; color: #555; font-size: 0.95em; }
  </style>
</head>
<body>
  <h1>BIFROST JumpBy / JumpTo demos</h1>""")

        # Render each group under its own <h2>, preserving the order in
        # which groups first appear in DEMO2_INDEX (2D first, then 3D).
        seen_groups = String[]
        for (g, _, _, _) in entries
            g in seen_groups || push!(seen_groups, g)
        end
        group_titles = Dict(
            "2D" => "2D scenes (inline SVG)",
            "3D" => "3D scenes (Plotly)",
        )
        for g in seen_groups
            heading = get(group_titles, g, g)
            println(io, "  <h2>$(heading)</h2>")
            println(io, "  <ul>")
            for (eg, title, path, desc) in entries
                eg == g || continue
                println(io, "    <li>")
                println(io, "      <a href=\"$(path)\">$(title)</a>")
                println(io, "      <p class=\"desc\">$(desc)</p>")
                println(io, "    </li>")
            end
            println(io, "  </ul>")
        end
        println(io, """</body>
</html>""")
    end

    println("Wrote demo2 index to: ", index_output)
    return index_output
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo2_all()
end
