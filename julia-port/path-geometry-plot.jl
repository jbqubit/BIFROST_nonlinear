"""
path-geometry-plot.jl

Interactive Plotly HTML for a [`Path`](path-geometry.jl): 3D centerline, a draggable cursor
driven by horizontal mouse position, a movable transverse square (normal–binormal plane),
short tangent/normal/binormal axes at the cursor, and an inset showing the unit disk in
local Frenet **n̂**–**b̂** coordinates.

This file depends only on `path-geometry.jl` (via the nested `PathGeometry` module) and the
Plotly CDN. It does not load or reference other `julia-port/` sources.

# Usage

    include("path-geometry-plot.jl")

    spec = PathGeometry.PathSpec()
    PathGeometry.straight!(spec; length = 0.2)
    PathGeometry.bend!(spec; radius = 0.4, angle = π / 2)
    path = PathGeometry.build(spec)
    write_path_geometry_plot3d(path, path.s_start, path.s_end; title = "Demo", output = "path.html")
"""

using LinearAlgebra

"""
    PathGeometry

Submodule wrapping `path-geometry.jl` (space-curve `Path` / `PathSpec`). Nesting keeps the
analytic path API in one namespace and avoids redefining those names in `Main` when this
file is `include`d alongside unrelated code.
"""
module PathGeometry
using LinearAlgebra
include(joinpath(@__DIR__, "path-geometry.jl"))
end

# ---------------------------------------------------------------------------
# Sampling (uses `frame` from PathGeometry — analytic Frenet data on `Path`)
# ---------------------------------------------------------------------------

"""
    sample_path_for_plot(path::PathGeometry.Path, s1, s2; n = 801)

Dense samples of position, Frenet frame, curvature, and torsions on `[s1, s2]`,
for plotting and for embedding in generated HTML.
"""
function sample_path_for_plot(path::PathGeometry.Path, s1::Real, s2::Real; n::Int = 801)
    @assert s2 > s1 "sample_path_for_plot: require s2 > s1"
    @assert n >= 2 "sample_path_for_plot: need at least two sample points"

    ss = collect(range(Float64(s1), Float64(s2); length = n))
    x = zeros(Float64, n)
    y = zeros(Float64, n)
    z = zeros(Float64, n)
    tx = zeros(Float64, n)
    ty = zeros(Float64, n)
    tz = zeros(Float64, n)
    nx = zeros(Float64, n)
    ny = zeros(Float64, n)
    nz = zeros(Float64, n)
    bx = zeros(Float64, n)
    by = zeros(Float64, n)
    bz = zeros(Float64, n)
    kappa = zeros(Float64, n)
    tau_geom = zeros(Float64, n)
    tau_mat = zeros(Float64, n)

    for i in eachindex(ss)
        fr = PathGeometry.frame(path, ss[i])
        p = fr.position
        x[i], y[i], z[i] = p[1], p[2], p[3]
        T = fr.tangent
        N = fr.normal
        B = fr.binormal
        tx[i], ty[i], tz[i] = T[1], T[2], T[3]
        nx[i], ny[i], nz[i] = N[1], N[2], N[3]
        bx[i], by[i], bz[i] = B[1], B[2], B[3]
        kappa[i] = fr.curvature
        tau_geom[i] = fr.geometric_torsion
        tau_mat[i] = fr.material_twist
    end

    return (;
        s = ss,
        x,
        y,
        z,
        tx,
        ty,
        tz,
        nx,
        ny,
        nz,
        bx,
        by,
        bz,
        kappa,
        tau_geom,
        tau_mat
    )
end

# ---------------------------------------------------------------------------
# Plotly HTML helpers (small JSON serializers for embedded numeric arrays)
# ---------------------------------------------------------------------------

function _js_real(x::Real)
    xf = Float64(x)
    if isnan(xf)
        return "NaN"
    elseif xf == Inf
        return "Infinity"
    elseif xf == -Inf
        return "-Infinity"
    else
        return string(xf)
    end
end

_js_array(xs::AbstractVector{<:Real}) = "[" * join(_js_real.(xs), ", ") * "]"

function _js_string_array(xs::AbstractVector{<:AbstractString})
    escaped = replace.(xs, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n")
    return "[" * join(["\"" * x * "\"" for x in escaped], ", ") * "]"
end

"""
    write_path_geometry_plot3d(path::PathGeometry.Path, s1, s2; n, output, title, plane_extent_frac, axis_extent_frac, twist_n_quad)

Write a standalone Plotly HTML file. Horizontal mouse position (without a mouse button pressed)
scrubs the arc-length parameter; the transverse plane and frame axes update accordingly.

Keywords `plane_extent_frac` and `axis_extent_frac` scale the half-width of the movable square
and the length of the **T**, **N**, **B** segments relative to the axis-aligned bounding box
diagonal of the path (before padding).

The main 3D scene axis ranges are fixed at generation time from
[`PathGeometry.bounding_box`](path-geometry.jl), expanded to keep the transverse square and
**T**/**N**/**B** axes in view. After each scrub update, the same ranges are re-applied with
`Plotly.relayout` so Plotly does not re-fit the scene to the moving mesh.

Open-circle markers are drawn at effective arc-length joins between authored segments (interior
segment boundaries) that fall within `[s1, s2]`.

A red arrow in the local **N̂**–**B̂** plane at the cursor points along
cos(Φ) N̂ + sin(Φ) B̂, where Φ is [`PathGeometry.total_material_twist`](path-geometry.jl) from `s1`
to the cursor arc length (same effective-``s`` as the plot), using `twist_n_quad` for overlay
quadrature.

Keyword `twist_n_quad` is passed to `total_material_twist` when building Φ at each sample index
(default 128).
"""
function write_path_geometry_plot3d(
    path::PathGeometry.Path,
    s1::Real,
    s2::Real;
    n::Int = 801,
    output::AbstractString = "path_geometry_3d.html",
    title::AbstractString = "Path geometry",
    plane_extent_frac::Float64 = 0.08,
    axis_extent_frac::Float64 = 0.06,
    twist_n_quad::Int = 128,
)
    samples = sample_path_for_plot(path, s1, s2; n = n)
    xs = samples.x
    ys = samples.y
    zs = samples.z

    bb = PathGeometry.bounding_box(path; n = max(n, 512))
    lo = Vector{Float64}(bb.lo)
    hi = Vector{Float64}(bb.hi)
    lo = min.(lo, [minimum(xs), minimum(ys), minimum(zs)])
    hi = max.(hi, [maximum(xs), maximum(ys), maximum(zs)])
    diag = norm(hi - lo)
    diag = diag > 0 ? diag : 1.0
    plane_half = Float64(plane_extent_frac) * diag
    axis_len = max(Float64(axis_extent_frac) * diag, 1e-9)
    pad = sqrt(2) * plane_half + axis_len
    lo_plot = lo .- pad
    hi_plot = hi .+ pad
    x_range_js = _js_array([lo_plot[1], hi_plot[1]])
    y_range_js = _js_array([lo_plot[2], hi_plot[2]])
    z_range_js = _js_array([lo_plot[3], hi_plot[3]])

    s1f = Float64(s1)
    s2f = Float64(s2)
    seg_bx = Float64[]
    seg_by = Float64[]
    seg_bz = Float64[]
    seg_bound_hover = String[]
    placed = path.placed_segments
    if length(placed) >= 2
        for i in 2:length(placed)
            sj = placed[i].s_offset_eff
            if sj < s1f || sj > s2f
                continue
            end
            p = PathGeometry.position(path, sj)
            push!(seg_bx, p[1])
            push!(seg_by, p[2])
            push!(seg_bz, p[3])
            push!(
                seg_bound_hover,
                "Segment boundary<br>s = $(sj) m<br>x, y, z = $(p[1]), $(p[2]), $(p[3])"
            )
        end
    end

    title_html = replace(replace(title, "&" => "&amp;"), "<" => "&lt;")

    # Unit circle in local (n̂, b̂) coordinates for the inset
    nθ = 121
    θs = collect(range(0.0, 2π; length = nθ))
    inset_x = cos.(θs)
    inset_y = sin.(θs)

    inset_tick_x = Float64[]
    inset_tick_y = Float64[]
    for ϕ in (0.0, 0.5π, π, 1.5π)
        push!(inset_tick_x, 0.88 * cos(ϕ), cos(ϕ), NaN)
        push!(inset_tick_y, 0.88 * sin(ϕ), sin(ϕ), NaN)
    end

    s_samples = Vector{Float64}(samples.s)
    integrated_tau_mat = Vector{Float64}(undef, length(s_samples))
    for i in eachindex(s_samples)
        integrated_tau_mat[i] = PathGeometry.total_material_twist(
            path;
            s_start = s1f,
            s_end = s_samples[i],
            n_quad = twist_n_quad,
        )
    end

    html = """
    <!--
      Main 3D plot (legend names match Plotly traces):
      - path: centerline polyline; faint gray line; small markers along the curve colored by
        arc length s (Turbo colorscale).
      - segment joins: open circles at effective arc-length boundaries between authored
        segments (within the plotted s interval).
      - start / end: filled markers at the first and last sample points of the plotted
        interval.
      - cursor: filled marker at the scrub position along the path.
      - normal–binormal plane: semi-transparent square in the local plane spanned by N̂ and B̂
        at the cursor (moves with scrub).
      - T̂: orange segment, unit tangent at the cursor.
      - N̂: blue segment, principal normal at the cursor.
      - B̂: green segment, binormal T̂×N̂ at the cursor.
      - ∫τ_mat: red arrow in the N̂–B̂ plane at the cursor; Φ = total_material_twist(path; s_start
        = plot start, s_end = cursor) (same length scale as T̂/N̂/B̂ axes).
      Inset plot: unit circle and axes in local transverse (n̂, b̂) coordinates at the same
      scrub position.
    -->
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="utf-8" />
      <meta name="viewport" content="width=device-width, initial-scale=1" />
      <title>$title_html</title>
      <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
      <style>
        html, body {
          width: 100%;
          height: 100%;
          margin: 0;
          font-family: sans-serif;
        }
        #viewer {
          position: relative;
          width: 100%;
          height: 100%;
          overflow: hidden;
        }
        #plot {
          width: 100%;
          height: 100%;
        }
        #inset {
          position: absolute;
          top: 16px;
          right: 16px;
          width: min(28vw, 340px);
          height: min(28vw, 340px);
          min-width: 220px;
          min-height: 220px;
          border: 1px solid rgba(0, 0, 0, 0.15);
          background: rgba(255, 255, 255, 0.9);
          box-shadow: 0 8px 24px rgba(0, 0, 0, 0.12);
        }
        #status {
          position: absolute;
          top: calc(16px + min(28vw, 340px) + 10px);
          right: 16px;
          width: min(28vw, 340px);
          min-width: 220px;
          padding: 10px 12px;
          box-sizing: border-box;
          border: 1px solid rgba(0, 0, 0, 0.15);
          background: rgba(255, 255, 255, 0.92);
          box-shadow: 0 8px 24px rgba(0, 0, 0, 0.12);
          font-size: 13px;
          line-height: 1.4;
          white-space: pre-line;
        }
        #scrub-help {
          position: absolute;
          left: 16px;
          bottom: 18px;
          padding: 8px 10px;
          font-size: 13px;
          background: rgba(255, 255, 255, 0.86);
          border: 1px solid rgba(0, 0, 0, 0.12);
        }
      </style>
    </head>
    <body>
      <div id="viewer">
        <div id="plot"></div>
        <div id="inset"></div>
        <div id="status"></div>
        <div id="scrub-help">Move the mouse left-to-right over the main view to move the plane along the path (arc length).</div>
      </div>
      <script>
        const xs = $(_js_array(samples.x));
        const ys = $(_js_array(samples.y));
        const zs = $(_js_array(samples.z));
        const ss = $(_js_array(samples.s));
        const tx = $(_js_array(samples.tx));
        const ty = $(_js_array(samples.ty));
        const tz = $(_js_array(samples.tz));
        const nx = $(_js_array(samples.nx));
        const ny = $(_js_array(samples.ny));
        const nz = $(_js_array(samples.nz));
        const bx = $(_js_array(samples.bx));
        const by = $(_js_array(samples.by));
        const bz = $(_js_array(samples.bz));
        const kappa = $(_js_array(samples.kappa));
        const tauGeom = $(_js_array(samples.tau_geom));
        const tauMat = $(_js_array(samples.tau_mat));
        const integratedTwist = $(_js_array(integrated_tau_mat));
        const planeHalf = $(_js_real(plane_half));
        const axisLen = $(_js_real(axis_len));

        const insetCircleX = $(_js_array(inset_x));
        const insetCircleY = $(_js_array(inset_y));
        const insetTickX = $(_js_array(inset_tick_x));
        const insetTickY = $(_js_array(inset_tick_y));

        const xRange = $x_range_js;
        const yRange = $y_range_js;
        const zRange = $z_range_js;

        const segBoundX = $(_js_array(seg_bx));
        const segBoundY = $(_js_array(seg_by));
        const segBoundZ = $(_js_array(seg_bz));
        const segBoundHover = $(_js_string_array(seg_bound_hover));

        const sceneAxisLayout = {
          "scene.xaxis.autorange": false,
          "scene.xaxis.range": xRange,
          "scene.yaxis.autorange": false,
          "scene.yaxis.range": yRange,
          "scene.zaxis.autorange": false,
          "scene.zaxis.range": zRange
        };

        function planeMesh(i) {
          const r = [xs[i], ys[i], zs[i]];
          const N = [nx[i], ny[i], nz[i]];
          const B = [bx[i], by[i], bz[i]];
          const h = planeHalf;
          const corners = [
            [r[0] + h * ( N[0] + B[0]), r[1] + h * ( N[1] + B[1]), r[2] + h * ( N[2] + B[2])],
            [r[0] + h * (-N[0] + B[0]), r[1] + h * (-N[1] + B[1]), r[2] + h * (-N[2] + B[2])],
            [r[0] + h * (-N[0] - B[0]), r[1] + h * (-N[1] - B[1]), r[2] + h * (-N[2] - B[2])],
            [r[0] + h * ( N[0] - B[0]), r[1] + h * ( N[1] - B[1]), r[2] + h * ( N[2] - B[2])]
          ];
          const x = corners.map(c => c[0]);
          const y = corners.map(c => c[1]);
          const z = corners.map(c => c[2]);
          const iTri = [0, 0];
          const jTri = [1, 2];
          const kTri = [2, 3];
          return { x, y, z, i: iTri, j: jTri, k: kTri };
        }

        function axisSegments(i) {
          const r = [xs[i], ys[i], zs[i]];
          const T = [tx[i], ty[i], tz[i]];
          const N = [nx[i], ny[i], nz[i]];
          const B = [bx[i], by[i], bz[i]];
          const L = axisLen;
          const seg = (v) => [
            [r[0], r[0] + L * v[0]],
            [r[1], r[1] + L * v[1]],
            [r[2], r[2] + L * v[2]]
          ];
          const Ta = seg(T);
          const Na = seg(N);
          const Ba = seg(B);
          return {
            T: { x: Ta[0], y: Ta[1], z: Ta[2] },
            N: { x: Na[0], y: Na[1], z: Na[2] },
            B: { x: Ba[0], y: Ba[1], z: Ba[2] }
          };
        }

        function twistArrowInNBPlane(i) {
          const phi = integratedTwist[i];
          const r = [xs[i], ys[i], zs[i]];
          const N = [nx[i], ny[i], nz[i]];
          const B = [bx[i], by[i], bz[i]];
          const ux = Math.cos(phi) * N[0] + Math.sin(phi) * B[0];
          const uy = Math.cos(phi) * N[1] + Math.sin(phi) * B[1];
          const uz = Math.cos(phi) * N[2] + Math.sin(phi) * B[2];
          const L = axisLen;
          return {
            x: [r[0], r[0] + L * ux],
            y: [r[1], r[1] + L * uy],
            z: [r[2], r[2] + L * uz]
          };
        }

        const pathTrace = {
          type: "scatter3d",
          mode: "lines+markers",
          x: xs,
          y: ys,
          z: zs,
          hoverinfo: "skip",
          line: {
            width: 4,
            color: "rgba(30, 30, 30, 0.25)"
          },
          marker: {
            size: 3.5,
            color: ss,
            colorscale: "Turbo",
            cmin: ss[0],
            cmax: ss[ss.length - 1],
            showscale: false
          },
          name: "path"
        };

        const segmentBoundaryTrace = {
          type: "scatter3d",
          mode: "markers",
          x: segBoundX,
          y: segBoundY,
          z: segBoundZ,
          text: segBoundHover,
          hoverinfo: "text",
          marker: {
            size: 8,
            symbol: "circle-open",
            color: "#7f7f7f",
            line: { width: 2, color: "#444444" }
          },
          name: "segment joins",
          visible: segBoundX.length > 0,
          showlegend: true
        };

        const startTrace = {
          type: "scatter3d",
          mode: "markers",
          x: [$(samples.x[1])],
          y: [$(samples.y[1])],
          z: [$(samples.z[1])],
          hovertemplate: "Start<br>x=%{x} m<br>y=%{y} m<br>z=%{z} m<extra></extra>",
          marker: { size: 7, color: "#2ca02c", symbol: "circle" },
          name: "start"
        };

        const endTrace = {
          type: "scatter3d",
          mode: "markers",
          x: [$(samples.x[end])],
          y: [$(samples.y[end])],
          z: [$(samples.z[end])],
          hovertemplate: "End<br>x=%{x} m<br>y=%{y} m<br>z=%{z} m<extra></extra>",
          marker: { size: 7, color: "#d62728", symbol: "circle" },
          name: "end"
        };

        const cursorTrace = {
          type: "scatter3d",
          mode: "markers",
          x: [xs[0]],
          y: [ys[0]],
          z: [zs[0]],
          hoverinfo: "skip",
          marker: { size: 6, color: "#111111", symbol: "circle" },
          name: "cursor"
        };

        const pm0 = planeMesh(0);
        const ax0 = axisSegments(0);

        const planeTrace = {
          type: "mesh3d",
          x: pm0.x,
          y: pm0.y,
          z: pm0.z,
          i: pm0.i,
          j: pm0.j,
          k: pm0.k,
          opacity: 0.38,
          color: "#6ba3ff",
          flatshading: true,
          hoverinfo: "skip",
          name: "normal–binormal plane"
        };

        const traceT = {
          type: "scatter3d",
          mode: "lines",
          x: ax0.T.x,
          y: ax0.T.y,
          z: ax0.T.z,
          line: { width: 5, color: "#ff7f0e" },
          hoverinfo: "skip",
          name: "T̂"
        };
        const traceN = {
          type: "scatter3d",
          mode: "lines",
          x: ax0.N.x,
          y: ax0.N.y,
          z: ax0.N.z,
          line: { width: 5, color: "#1f77b4" },
          hoverinfo: "skip",
          name: "N̂"
        };
        const traceB = {
          type: "scatter3d",
          mode: "lines",
          x: ax0.B.x,
          y: ax0.B.y,
          z: ax0.B.z,
          line: { width: 5, color: "#2ca02c" },
          hoverinfo: "skip",
          name: "B̂"
        };

        const tw0 = twistArrowInNBPlane(0);
        const traceTwistInt = {
          type: "scatter3d",
          mode: "lines+markers",
          x: tw0.x,
          y: tw0.y,
          z: tw0.z,
          line: { width: 5, color: "#e31a1c" },
          marker: {
            size: [2, 9],
            color: "#e31a1c",
            symbol: ["circle", "diamond"],
            line: { width: 1, color: "#7f0000" }
          },
          hoverinfo: "skip",
          name: "∫τ_mat ds"
        };

        const layout = {
          title: $(repr(title)),
          uirevision: "path-geometry-plot",
          showlegend: true,
          legend: { x: 0.02, y: 0.98 },
          margin: { l: 0, r: 0, b: 0, t: 48 },
          scene: {
            xaxis: { title: "x (m)", range: xRange, autorange: false, showspikes: false },
            yaxis: { title: "y (m)", range: yRange, autorange: false, showspikes: false },
            zaxis: { title: "z (m)", range: zRange, autorange: false, showspikes: false },
            camera: { eye: { x: 1.6, y: 1.6, z: 0.9 } },
            aspectmode: "data"
          }
        };

        Plotly.newPlot("plot", [pathTrace, segmentBoundaryTrace, startTrace, endTrace, cursorTrace, planeTrace, traceT, traceN, traceB, traceTwistInt], layout, {
          responsive: true,
          scrollZoom: true,
          displaylogo: false
        });

        const insetCircleTrace = {
          type: "scatter",
          mode: "lines",
          x: insetCircleX,
          y: insetCircleY,
          line: { width: 2, color: "rgba(30, 30, 30, 0.45)" },
          hoverinfo: "skip",
          name: "unit circle"
        };
        const insetAxesTrace = {
          type: "scatter",
          mode: "lines",
          x: [-1.25, 1.25, null, 0, 0],
          y: [0, 0, null, -1.25, 1.25],
          line: { width: 1.5, color: "rgba(0, 0, 0, 0.35)" },
          hoverinfo: "skip",
          showlegend: false
        };
        const insetTicksTrace = {
          type: "scatter",
          mode: "lines",
          x: insetTickX,
          y: insetTickY,
          line: { width: 2, color: "#111111" },
          hoverinfo: "skip",
          showlegend: false
        };
        const insetLabels = {
          type: "scatter",
          mode: "text",
          x: [1.15, 0.08, 0.08],
          y: [0.08, 1.15, -1.2],
          text: ["n̂", "b̂", "Transverse (n̂, b̂) coordinates"],
          textposition: "middle center",
          textfont: { size: 14, color: "#111111", family: "sans-serif" },
          hoverinfo: "skip",
          showlegend: false
        };

        const insetLayout = {
          margin: { l: 36, r: 18, b: 36, t: 28 },
          title: { text: "Local Frenet section", font: { size: 14 } },
          xaxis: {
            title: "component along N̂",
            range: [-1.35, 1.35],
            scaleanchor: "y",
            scaleratio: 1,
            zeroline: false
          },
          yaxis: {
            title: "component along B̂",
            range: [-1.35, 1.35],
            zeroline: false
          },
          showlegend: false
        };

        Plotly.newPlot("inset", [insetCircleTrace, insetAxesTrace, insetTicksTrace, insetLabels], insetLayout, {
          responsive: true,
          displaylogo: false,
          staticPlot: true
        });

        const statusBox = document.getElementById("status");
        let activeIndex = 0;

        function formatStatus(index) {
          const deg = integratedTwist[index] * (180 / Math.PI);
          return [
            "Arc length s = " + ss[index].toFixed(5) + " m",
            "x, y, z = " + xs[index].toFixed(4) + ", " + ys[index].toFixed(4) + ", " + zs[index].toFixed(4) + " m",
            "κ = " + kappa[index].toExponential(4) + " 1/m",
            "τ_geom = " + tauGeom[index].toExponential(4) + " rad/m",
            "τ_mat  = " + tauMat[index].toExponential(4) + " rad/m",
            "∫τ_mat ds = " + integratedTwist[index].toFixed(5) + " rad (" + deg.toFixed(2) + "°)"
          ].join("\\n");
        }

        function updateCursor(index) {
          activeIndex = Math.max(0, Math.min(xs.length - 1, index));
          const pm = planeMesh(activeIndex);
          const ax = axisSegments(activeIndex);
          statusBox.textContent = formatStatus(activeIndex);
          const chain = (p, fn) => (p && typeof p.then === "function" ? p.then(fn) : fn());
          let p = Plotly.restyle("plot", {
            x: [[xs[activeIndex]]],
            y: [[ys[activeIndex]]],
            z: [[zs[activeIndex]]]
          }, [4]);
          p = chain(p, () => Plotly.restyle("plot", {
            x: [pm.x],
            y: [pm.y],
            z: [pm.z],
            i: [pm.i],
            j: [pm.j],
            k: [pm.k]
          }, [5]));
          p = chain(p, () => Plotly.restyle("plot", {
            x: [ax.T.x],
            y: [ax.T.y],
            z: [ax.T.z]
          }, [6]));
          p = chain(p, () => Plotly.restyle("plot", {
            x: [ax.N.x],
            y: [ax.N.y],
            z: [ax.N.z]
          }, [7]));
          p = chain(p, () => Plotly.restyle("plot", {
            x: [ax.B.x],
            y: [ax.B.y],
            z: [ax.B.z]
          }, [8]));
          p = chain(p, () => {
            const tw = twistArrowInNBPlane(activeIndex);
            return Plotly.restyle("plot", {
              x: [tw.x],
              y: [tw.y],
              z: [tw.z]
            }, [9]);
          });
          chain(p, () => Plotly.relayout("plot", sceneAxisLayout));
        }

        const viewer = document.getElementById("viewer");
        viewer.addEventListener("mousemove", event => {
          if (event.buttons !== 0) return;
          const rect = viewer.getBoundingClientRect();
          const t = Math.max(0, Math.min(1, (event.clientX - rect.left) / rect.width));
          const idx = Math.round(t * (xs.length - 1));
          if (idx !== activeIndex) updateCursor(idx);
        });

        viewer.addEventListener("touchmove", event => {
          if (event.touches.length === 0) return;
          const rect = viewer.getBoundingClientRect();
          const t = Math.max(0, Math.min(1, (event.touches[0].clientX - rect.left) / rect.width));
          const idx = Math.round(t * (xs.length - 1));
          if (idx !== activeIndex) updateCursor(idx);
        }, { passive: true });

        updateCursor(0);
      </script>
    </body>
    </html>
    """

    open(output, "w") do io
        write(io, html)
    end
    return output
end
