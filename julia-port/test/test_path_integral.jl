using Test
using LinearAlgebra

if !isdefined(Main, :propagate_piecewise)
    include("../path-integral.jl")
end

const PI_MATRIX_ATOL = 1e-9
const SIGMA_X = ComplexF64[0 1; 1 0]
const SIGMA_Y = ComplexF64[0 -im; im 0]
const SIGMA_Z = ComplexF64[1 0; 0 -1]
const DGD_ATOL = 1e-9
const MATRIX_ATOL = 1e-9

pauli_axis(axis::AbstractVector{<:Real}) = axis[1] .* SIGMA_X + axis[2] .* SIGMA_Y + axis[3] .* SIGMA_Z

function piecewise_constant_matrix(breaks::Vector{Float64}, values::Vector{Matrix{ComplexF64}})
    @assert issorted(breaks)
    @assert length(values) == length(breaks) - 1
    return function (s::Real)
        idx = min(searchsortedlast(breaks, Float64(s)), length(values))
        return values[idx]
    end
end

function propagate_dgd_case(K, Komega, breaks; rtol = 1e-11, atol = 1e-13, h_init = 0.05)
    J, G, stats = propagate_piecewise_Kω(K, Komega, breaks; rtol = rtol, atol = atol, h_init = h_init)
    return J, G, stats
end

accepted_steps(stats) = sum(st.accepted_steps for st in stats)

# ─── T-PHYSICS: noncommuting piecewise-constant generators ──────────────────────────────────────
#
# For piecewise-constant K with a breakpoint, path-ordering gives
#
#   J = exp(K₂ Δs₂) · exp(K₁ Δs₁)
#
# where K₁ and K₂ are on different Pauli axes, so they do NOT commute.
# The expected result is the explicit matrix product of the two closed-form exponentials.
# This is the minimal noncommuting benchmark: if the propagator ignored ordering it would
# compute exp((K₁+K₂)/2 · L) instead, which is wrong when [K₁,K₂] ≠ 0.

@testset "Path-ordering: noncommuting piecewise-constant generators" begin
    # T-PHYSICS
    # K₁ = β₁·i·σ_x on [0, L₁],  K₂ = β₂·i·σ_z on [L₁, L₁+L₂]
    # J_expected = exp(β₂·i·σ_z·L₂) · exp(β₁·i·σ_x·L₁)   (right-to-left composition)
    L1 = 0.7
    L2 = 1.3
    β1 = 0.6
    β2 = 0.4
    breaks = [0.0, L1, L1 + L2]

    K = function (s::Real)
        s < L1 ? β1 * im .* SIGMA_X : β2 * im .* SIGMA_Z
    end

    J_num, _ = propagate_piecewise(K, breaks; rtol = 1e-11, atol = 1e-13)

    J_seg1 = exp_jones_generator(L1 * β1 * im .* SIGMA_X)
    J_seg2 = exp_jones_generator(L2 * β2 * im .* SIGMA_Z)
    J_expected = J_seg2 * J_seg1

    @test J_num ≈ J_expected atol = PI_MATRIX_ATOL

    # Confirm that the wrong (commuting) result is actually different, so this test has
    # power to detect an ordering bug.
    J_wrong = exp_jones_generator(L1 * β1 * im .* SIGMA_X + L2 * β2 * im .* SIGMA_Z)
    @test opnorm(J_expected - J_wrong) > 1e-3
end

@testset "Path-ordering: three noncommuting segments" begin
    # T-PHYSICS
    # Three segments on three different Pauli axes — each pair is noncommuting.
    # J_expected = exp(K₃ L₃) · exp(K₂ L₂) · exp(K₁ L₁)
    L1, L2, L3 = 0.5, 0.8, 0.6
    β1, β2, β3 = 0.5, -0.3, 0.7
    breaks = [0.0, L1, L1 + L2, L1 + L2 + L3]

    function K3(s::Real)
        if s < L1
            return β1 * im .* SIGMA_X
        elseif s < L1 + L2
            return β2 * im .* SIGMA_Y
        else
            return β3 * im .* SIGMA_Z
        end
    end

    J_num, _ = propagate_piecewise(K3, breaks; rtol = 1e-11, atol = 1e-13)

    J1 = exp_jones_generator(L1 * β1 * im .* SIGMA_X)
    J2 = exp_jones_generator(L2 * β2 * im .* SIGMA_Y)
    J3 = exp_jones_generator(L3 * β3 * im .* SIGMA_Z)
    J_expected = J3 * J2 * J1

    @test J_num ≈ J_expected atol = PI_MATRIX_ATOL
end

@testset "Noncommuting smooth K: convergence to analytic result" begin
    # T-PHYSICS
    # K(s) = α·i·σ_x·cos(s) + β·i·σ_z·sin(s)  on [0, π].
    # σ_x and σ_z do not commute, so no closed form exists; instead we check that
    # the adaptive integrator converges: coarse and fine tolerances agree to fine-tol precision.
    α = 0.8
    β = 0.5
    breaks = [0.0, Float64(π)]

    K_smooth = s -> α * im * cos(s) .* SIGMA_X + β * im * sin(s) .* SIGMA_Z

    J_fine, stats_fine = propagate_piecewise(K_smooth, breaks; rtol = 1e-11, atol = 1e-13)
    J_coarse, stats_coarse = propagate_piecewise(K_smooth, breaks; rtol = 1e-6, atol = 1e-8)

    # Fine result is the reference; coarse should agree to its own tolerance (~1e-6 rtol on a unit matrix).
    @test opnorm(J_fine - J_coarse) <= 1e-4

    # Fine integrator should be unitary (SU(2) invariant).
    @test opnorm(J_fine' * J_fine - Matrix{ComplexF64}(I, 2, 2)) <= PI_MATRIX_ATOL

    # Sanity: the result is not the identity — K is non-trivial.
    @test opnorm(J_fine - Matrix{ComplexF64}(I, 2, 2)) > 1e-2
end

@testset "Noncommuting: unitary preservation" begin
    # T-GUARDRAIL
    # For any noncommuting piecewise K, the propagated J must remain unitary to solver tolerance.
    breaks = [0.0, 0.5, 1.2, 2.0]
    mixed_axis = (SIGMA_X .+ SIGMA_Y .+ SIGMA_Z) ./ sqrt(3)
    K_mixed = function (s::Real)
        if s < 0.5
            return 0.6im .* SIGMA_X
        elseif s < 1.2
            return -0.4im .* SIGMA_Z
        else
            return 0.3im .* mixed_axis
        end
    end

    J, _ = propagate_piecewise(K_mixed, breaks; rtol = 1e-11, atol = 1e-13)

    @test opnorm(J' * J - Matrix{ComplexF64}(I, 2, 2)) <= PI_MATRIX_ATOL
    @test abs(det(J) - 1.0) <= PI_MATRIX_ATOL
end

@testset "Noncommuting: reverse path gives inverse Jones matrix" begin
    # T-PHYSICS
    # If K_rev(s) = -K(L - s), then propagating K_rev over [0, L] gives J⁻¹.
    # This tests that path-ordering is globally consistent: reversing the path
    # undoes the forward propagation.
    L = 1.8
    breaks_fwd = [0.0, 0.6, 1.2, L]
    β_vals = [0.5, -0.3, 0.7]
    axes = [SIGMA_X, SIGMA_Z, SIGMA_Y]

    K_fwd = function (s::Real)
        if s < 0.6
            return β_vals[1] * im .* axes[1]
        elseif s < 1.2
            return β_vals[2] * im .* axes[2]
        else
            return β_vals[3] * im .* axes[3]
        end
    end

    # K_rev(s) = -K_fwd(L - s)
    K_rev = s -> -K_fwd(L - s)
    breaks_rev = [0.0, L - 1.2, L - 0.6, L]

    J_fwd, _ = propagate_piecewise(K_fwd, breaks_fwd; rtol = 1e-11, atol = 1e-13)
    J_rev, _ = propagate_piecewise(K_rev, breaks_rev; rtol = 1e-11, atol = 1e-13)

    @test J_rev * J_fwd ≈ Matrix{ComplexF64}(I, 2, 2) atol = PI_MATRIX_ATOL
end

@testset "DGD exact cases" begin
    @testset "Zero Komega gives zero DGD" begin
        K = s -> 0.3im .* SIGMA_Z
        Komega = s -> zeros(ComplexF64, 2, 2)
        J, G, stats = propagate_dgd_case(K, Komega, [0.0, 2.0])

        @test accepted_steps(stats) > 0
        @test opnorm(J - Matrix{ComplexF64}(I, 2, 2)) > 1e-3
        @test opnorm(G) <= MATRIX_ATOL
        @test output_dgd(J, G) <= DGD_ATOL
    end

    @testset "Closed form along sigma_z" begin
        length_scale = 2.0
        beta = 0.4
        beta_prime = 0.7
        K = s -> beta * im .* SIGMA_Z
        Komega = s -> beta_prime * im .* SIGMA_Z
        J, G, _ = propagate_dgd_case(K, Komega, [0.0, length_scale])

        expected_pmd = length_scale * beta_prime .* SIGMA_Z
        @test pmd_generator(J, G) ≈ expected_pmd atol = MATRIX_ATOL
        @test output_dgd(J, G) ≈ 2 * abs(length_scale * beta_prime) atol = DGD_ATOL
    end

    @testset "Closed form along arbitrary axis" begin
        length_scale = 2.0
        beta = 0.4
        beta_prime = 0.7
        axis = [1.0, 1.0, 1.0]
        axis ./= norm(axis)
        axis_matrix = pauli_axis(axis)
        K = s -> beta * im .* axis_matrix
        Komega = s -> beta_prime * im .* axis_matrix
        J, G, _ = propagate_dgd_case(K, Komega, [0.0, length_scale])

        expected_pmd = length_scale * beta_prime .* axis_matrix
        @test pmd_generator(J, G) ≈ expected_pmd atol = MATRIX_ATOL
        @test output_dgd(J, G) ≈ 2 * abs(length_scale * beta_prime) atol = DGD_ATOL
    end

    @testset "Piecewise cancellation on one axis" begin
        breaks = [0.0, 1.0, 2.0]
        K = s -> 0.5im .* SIGMA_Z
        Komega = piecewise_constant_matrix(
            breaks,
            [
                0.8im .* SIGMA_Z,
                -0.8im .* SIGMA_Z
            ]
        )
        J, G, stats = propagate_dgd_case(K, Komega, breaks)

        @test length(stats) == 2
        @test pmd_generator(J, G) ≈ zeros(ComplexF64, 2, 2) atol = MATRIX_ATOL
        @test output_dgd(J, G) <= DGD_ATOL
    end

    @testset "Piecewise commuting accumulation" begin
        breaks = [0.0, 0.5, 1.25, 2.0]
        beta_vals = [0.3, 0.4, -0.2]
        beta_prime_vals = [0.25, -0.5, 1.25]
        K = piecewise_constant_matrix(breaks, [b * im .* SIGMA_Z for b in beta_vals])
        Komega = piecewise_constant_matrix(breaks, [bp * im .* SIGMA_Z for bp in beta_prime_vals])
        J, G, _ = propagate_dgd_case(K, Komega, breaks)

        expected_coeff =
            (breaks[2] - breaks[1]) * beta_prime_vals[1] +
            (breaks[3] - breaks[2]) * beta_prime_vals[2] +
            (breaks[4] - breaks[3]) * beta_prime_vals[3]

        @test pmd_generator(J, G) ≈ expected_coeff .* SIGMA_Z atol = MATRIX_ATOL
        @test output_dgd(J, G) ≈ 2 * abs(expected_coeff) atol = DGD_ATOL
    end
end

@testset "DGD invariants and regression" begin
    @testset "Constant basis rotation preserves DGD" begin
        breaks = [0.0, 0.5, 1.25, 2.0]
        K = piecewise_constant_matrix(
            breaks,
            [
                0.3im .* SIGMA_Z,
                0.4im .* SIGMA_Z,
                -0.2im .* SIGMA_Z
            ]
        )
        Komega = piecewise_constant_matrix(
            breaks,
            [
                0.25im .* SIGMA_Z,
                -0.5im .* SIGMA_Z,
                1.25im .* SIGMA_Z
            ]
        )
        J_ref, G_ref, _ = propagate_dgd_case(K, Komega, breaks)
        H_ref = pmd_generator(J_ref, G_ref)

        rotation = exp_jones_generator(0.37im .* SIGMA_Y)
        K_rot = s -> rotation * K(s) * rotation'
        Komega_rot = s -> rotation * Komega(s) * rotation'
        J_rot, G_rot, _ = propagate_dgd_case(K_rot, Komega_rot, breaks)
        H_rot = pmd_generator(J_rot, G_rot)

        @test H_rot ≈ rotation * H_ref * rotation' atol = MATRIX_ATOL
        @test output_dgd(J_rot, G_rot) ≈ output_dgd(J_ref, G_ref) atol = DGD_ATOL
    end

    @testset "Phase scaling leaves PMD and DGD unchanged" begin
        J, G, _ = propagate_dgd_case(s -> 0.4im .* SIGMA_Z, s -> 0.7im .* SIGMA_Z, [0.0, 2.0])
        phase = cis(0.73)

        @test pmd_generator(phase .* J, phase .* G) ≈ pmd_generator(J, G) atol = MATRIX_ATOL
        @test output_dgd(phase .* J, phase .* G) ≈ output_dgd(J, G) atol = DGD_ATOL
    end

    @testset "Noncommuting K with zero Komega still gives zero DGD" begin
        breaks = [0.0, 0.7, 1.4, 2.0]
        mixed_axis = (SIGMA_X .+ SIGMA_Z) ./ sqrt(2)
        K = piecewise_constant_matrix(
            breaks,
            [
                0.6im .* SIGMA_X,
                -0.2im .* SIGMA_Z,
                0.4im .* mixed_axis
            ]
        )
        Komega = s -> zeros(ComplexF64, 2, 2)
        J, G, _ = propagate_dgd_case(K, Komega, breaks)

        @test opnorm(J - Matrix{ComplexF64}(I, 2, 2)) > 1e-3
        @test opnorm(G) <= MATRIX_ATOL
        @test output_dgd(J, G) <= DGD_ATOL
    end

    @testset "Direct sensitivity agrees with finite difference" begin
        length_scale = 1.75
        breaks = [0.0, length_scale]
        a0, a1 = 0.31, 0.07
        b0, b1 = -0.22, 0.05
        c0, c1 = 0.18, -0.03

        function K_param(s, omega)
            return im * ((a0 + a1 * omega) * sin(0.7 * s)) .* SIGMA_X +
                   im * ((b0 + b1 * omega) * cos(0.4 * s)) .* SIGMA_Y +
                   im * (c0 + c1 * omega) .* SIGMA_Z
        end

        K_of_omega(omega) = s -> K_param(s, omega)
        Komega = s -> im * (a1 * sin(0.7 * s)) .* SIGMA_X +
                      im * (b1 * cos(0.4 * s)) .* SIGMA_Y +
                      im * c1 .* SIGMA_Z

        J_direct, G_direct, _ = propagate_dgd_case(
            K_of_omega(0.0),
            Komega,
            breaks;
            rtol = 1e-10,
            atol = 1e-12,
            h_init = 0.02
        )

        delta_omega = 1e-5
        J_plus, _ = propagate_piecewise(K_of_omega(delta_omega), breaks; rtol = 1e-10, atol = 1e-12, h_init = 0.02)
        J_minus, _ = propagate_piecewise(K_of_omega(-delta_omega), breaks; rtol = 1e-10, atol = 1e-12, h_init = 0.02)
        G_fd = (J_plus - J_minus) ./ (2 * delta_omega)

        @test opnorm(G_direct - G_fd) <= 1e-8
        @test opnorm(pmd_generator(J_direct, G_direct) - pmd_generator(J_direct, G_fd)) <= 1e-8
        @test abs(output_dgd(J_direct, G_direct) - output_dgd(J_direct, G_fd)) <= 1e-8
    end
end

@testset "DGD fiber integration" begin
    bend_breaks = [0.0, 1.0, 3.0]
    twist_breaks = [0.0, 2.0, 3.0]
    explicit_union = [0.0, 1.0, 2.0, 3.0]

    bend = BendSource(
        PiecewiseProfile(bend_breaks, [s -> 1.0, s -> 2.0]),
        PiecewiseProfile(bend_breaks, [s -> 0.0, s -> pi / 4]),
        k2 -> 0.9 * k2,
        k2 -> -0.35 * k2;
        breakpoints = copy(bend_breaks)
    )
    twist = TwistSource(
        PiecewiseProfile(twist_breaks, [s -> 0.2, s -> -0.5]),
        tau -> 0.6 * tau,
        tau -> -0.8 * tau;
        breakpoints = copy(twist_breaks)
    )

    fiber = Fiber(0.0, 3.0, AbstractBirefringenceSource[bend, twist])

    @testset "Automatic breakpoint union matches explicit partition" begin
        J_auto, G_auto, stats_auto = propagate_fiber_Kω(fiber; rtol = 1e-11, atol = 1e-13, h_init = 0.05)
        J_explicit, G_explicit, stats_explicit = propagate_piecewise_Kω(
            generator_K(fiber),
            generator_Kω(fiber),
            explicit_union;
            rtol = 1e-11,
            atol = 1e-13,
            h_init = 0.05
        )

        @test fiber_breakpoints(fiber) == explicit_union
        @test length(stats_auto) == length(explicit_union) - 1
        @test length(stats_explicit) == length(explicit_union) - 1
        @test J_auto ≈ J_explicit atol = MATRIX_ATOL
        @test G_auto ≈ G_explicit atol = MATRIX_ATOL
        @test output_dgd(J_auto, G_auto) ≈ output_dgd(J_explicit, G_explicit) atol = DGD_ATOL
    end
end
