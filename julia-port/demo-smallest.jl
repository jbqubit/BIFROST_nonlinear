include("material-properties.jl")
include("fiber-cross-section.jl")
include("path-geometry.jl")
include("path-integral.jl")

xs = FiberCrossSection(
    GermaniaSilicaGlass(0.036),
    GermaniaSilicaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like",
)

spec = PathSpecBuilder()
straight!(spec; length = 0.5, meta = [Nickname("lead-in")])
bend!(spec; radius = 0.05, angle = pi / 2, meta = [Nickname("90 deg bend")])
straight!(spec; length = 0.5, meta = [Nickname("lead-out")])

path = build(spec)
fiber = Fiber(path; cross_section = xs, T_ref_K = 297.15)

J, stats = propagate_fiber(fiber; λ_m = 1550e-9, rtol = 1e-9, verbose = false)

println("J =")
display(J)
println("intervals = ", length(stats))
