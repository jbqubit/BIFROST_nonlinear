# Julia Unit-Test Guide

This note explains the physical and numerical logic behind the Julia test suite in [`runtests.jl`](runtests.jl). It is written for a reader who wants to understand what each test means physically, what identity it is checking, and which Julia subroutines are being exercised.
It is intended to be self-sufficient and not to rely on any separate design-note file remaining in the repository.

The test entrypoint runs:

- `test_fiber_path.jl`
- `test_path_integral.jl` (Jones propagation, path-ordering checks, and DGD / Kω)
- `test_material_properties.jl`
- `fiber-cross-section.jl` (in-repo include)
- `test_paddle_transfer.jl`
- `test_fiber_cross_section.jl`

## Governing equations

The Jones propagator is defined by

$$
\frac{dJ}{ds} = K(s)\,J, \qquad J(0)=I.
$$

The DGD machinery introduces

$$
G(s,\omega) \equiv \partial_\omega J(s,\omega),
$$

and propagates the coupled system

$$
\frac{dJ}{ds}=KJ, \qquad \frac{dG}{ds}=K_\omega J + KG,
$$

with

$$
J(0)=I, \qquad G(0)=0.
$$

At the fiber output, the PMD generator used by the code is

$$
H_{\mathrm{PMD}} = -i\,J^{-1}G,
$$

implemented by `pmd_generator(J, G)`, and the scalar DGD reported by the code is

$$
\mathrm{DGD} = \lambda_{\max}(H_{\mathrm{PMD}})-\lambda_{\min}(H_{\mathrm{PMD}}),
$$

implemented by `output_dgd(J, G)`.

## `test_fiber_path.jl`

This file is the geometry-first test suite for the fiber path layer. It is discussed first because it validates the centerline construction that sits underneath the visualization and higher-level fiber infrastructure. These tests do not probe Jones propagation or DGD directly. They probe the geometric reconstruction returned by `sample_fiber_centerline`.

The guiding idea is to use one explicit $3$-dimensional rigid-body oracle for both planar and nonplanar bend sequences. A planar test is then just the special case in which the bend azimuth remains fixed and the resulting path stays in a coordinate plane.

### Geometry convention used by the tests

For each bend segment, the test-side oracle uses:

- bend radius $R>0$
- bend azimuth $\alpha$
- signed bend sweep $\theta$

The current implementation in `sample_fiber_centerline` interprets the bend azimuth as a lab-frame curvature direction, with signed bend handled by reversing that direction:

$$
\alpha_{\mathrm{eff}} =
\begin{cases}
\alpha, & \theta \ge 0, \\
\alpha + \pi, & \theta < 0.
\end{cases}
$$

and the propagated arc length is

$$
\Delta s = R\,|\theta|.
$$

The local curvature vector is therefore

$$
\kappa = \left(\frac{\cos \alpha_{\mathrm{eff}}}{R},
\frac{\sin \alpha_{\mathrm{eff}}}{R},
0\right),
$$

and the instantaneous rotation axis for the tangent is

$$
u = \frac{\Omega}{\|\Omega\|}
=
\left(-\sin \alpha_{\mathrm{eff}},\,
\cos \alpha_{\mathrm{eff}},\,
0\right),
$$

because the code advances the tangent using

$$
\Omega = (-\kappa_y,\kappa_x,0).
$$

For a constant-curvature segment, the exact endpoint update used by the oracle is

$$
r_{\mathrm{next}} = r + A(R,\theta,u)\,T,
$$

where $T$ is the incoming tangent and

$$
A(R,\theta,u)
=
R\theta I
+ R(1-\cos\theta)[u]_\times
+ R(\theta-\sin\theta)[u]_\times^2.
$$

The outgoing tangent is

$$
T_{\mathrm{next}} = \mathcal{R}(u,\theta)\,T,
$$

with $\mathcal{R}(u,\theta)$ the Rodrigues rotation about axis $u$ by angle $\theta$.

This is implemented in the test-side subroutines:

- `canonical_bend_segment`
- `skew_matrix`
- `rotate_vector_about_axis`
- `exact_endpoint_and_tangent`

The actual infrastructure under test is built with:

- `build_centerline_test_fiber`
- `build_straight_test_fiber`
- `sample_final_centerline_state`

and ultimately exercises:

- `BendSource`
- `TwistSource`
- `Fiber`
- `bend_geometry`
- `twist_rate`
- `sample_fiber_centerline`

### Test set `fiber-path centerline geometry`

#### `Straight section`

This is the zero-curvature baseline. The test builds a straight fiber of length $L=2.5$ and checks

$$
r(L) = (0,0,L),
\qquad
T(L) = (0,0,1).
$$

It also checks that the sampled path remains on the $z$ axis:

$$
x(s)=0,\qquad y(s)=0.
$$

Why this matters:
it catches bad initialization, accidental offsets, and any spurious curvature in the centerline reconstruction.

#### `Explicit endpoint cases`

These are exact benchmark trajectories listed in `EXPLICIT_PATH_CASES`.
For each case, the test compares:

1. the closed-form endpoint and tangent from `exact_endpoint_and_tangent`
2. the expected analytic endpoint/tangent hard-coded into the test table
3. the numerically sampled endpoint and tangent from `sample_fiber_centerline`

The cases include:

- `single_quarter_circle_xz_plane`
  $$
  (R,\alpha,\theta)=(1,0,\pi/2)
  \Rightarrow
  r=(1,0,1),\quad T=(1,0,0)
  $$
- `single_half_circle_xz_plane`
  $$
  (R,\alpha,\theta)=(2,0,\pi)
  \Rightarrow
  r=(4,0,0),\quad T=(0,0,-1)
  $$
- `two_quarters_make_half_circle`
  two consecutive quarter circles giving the same final geometry as a half circle
- `three_sixty_degree_bends_make_half_circle`
  the same half-circle endpoint realized by three $60^\circ$ segments
- `single_quarter_circle_yz_plane`
  $$
  (R,\alpha,\theta)=(1,\pi/2,\pi/2)
  \Rightarrow
  r=(0,1,1),\quad T=(0,1,0)
  $$
- `single_oblique_quarter_circle`
  $$
  (R,\alpha,\theta)=(2,\pi/4,\pi/2)
  \Rightarrow
  r=(\sqrt{2},\sqrt{2},2),\quad
  T=\left(\frac{1}{\sqrt 2},\frac{1}{\sqrt 2},0\right)
  $$
- `orthogonal_quarters_lab_frame`
  a genuinely $3$-dimensional composition under the current lab-frame azimuth convention
- `orthogonal_quarters_swapped_order`
  the swapped-order companion to the previous case

Why this matters:
these are the anchor tests for sign convention, arc-length convention, and the distinction between planar and fully $3$-dimensional curvature sequences.

For the planar cases, the suite also checks

$$
\max_s |y(s)| \approx 0,
$$

so the path stays in the expected plane.

#### `Segmentation invariance`

This test compares three discretizations of the same total bend:

$$
(1,0,\pi/2),
$$

versus

$$
(1,0,\pi/4),(1,0,\pi/4),
$$

versus

$$
4 \times (1,0,\pi/8).
$$

All three should have the same endpoint and final tangent:

$$
r=(1,0,1),\qquad T=(1,0,0).
$$

Why this matters:
the geometry layer should treat a sequence of smaller bends as a faithful segmentation of the same smooth circular path, not as a different physical object.

#### `Hard exact oracle case`

This is the less obvious two-segment case

$$
(1,\pi/4,\pi/3),\qquad (1,-\pi/4,\pi/3),
$$

for which the test does not rely on visual intuition. It relies on the exact rigid-body oracle and checks that the sampler reproduces both the final endpoint and the final tangent.

Why this matters:
it probes mixed azimuths in a way that is hard to fake by accidental symmetry.

#### `Order matters in 3D`

This compares

$$
(1,0,\pi/2),\ (1,\pi/2,\pi/2)
$$

against

$$
(1,\pi/2,\pi/2),\ (1,0,\pi/2).
$$

The test asserts that the two endpoints are substantially different.

Why this matters:
geometry composition is noncommutative. Changing the order of rotated bend segments should change the final centerline.

## `test_path_integral_sources.jl`

This file checks the basic source-based fiber abstraction and a few calibration points for the propagator.

### Test set `Fiber source API`

#### Breakpoint union

The test constructs a bend source and a twist source with different discontinuity sets and checks that

$$
\mathrm{fiber\_breakpoints}(F)
=
\{s_{\min},s_{\max}\}\cup
\mathrm{breakpoints}(\mathrm{bend})\cup
\mathrm{breakpoints}(\mathrm{twist}).
$$

Why this matters:
the adaptive propagator should never step across a true discontinuity in either $K(s)$ or $K_\omega(s)$.

Subroutines exercised:
`Fiber`, `fiber_breakpoints`

#### Additivity of source contributions

At a sample point $s=0.25$, the test checks

$$
K(s)=K_{\mathrm{bend}}(s)+K_{\mathrm{twist}}(s),
$$

and likewise

$$
K_\omega(s)=K_{\omega,\mathrm{bend}}(s)+K_{\omega,\mathrm{twist}}(s),
$$

via the subroutines `generator_K`, `generator_Kω`, `generator_K_contribution`, and `generator_Kω_contribution`.

Why this matters:
the entire fiber architecture assumes linear superposition of source contributions at the generator level.

Subroutines exercised:
`generator_K`, `generator_Kω`, `generator_K_contribution`, `generator_Kω_contribution`

#### Coverage-gap rejection

The `bad_twist` construction leaves part of the interval uncovered, and the `Fiber` constructor is required to throw.

The model assumption is:
if a source is inactive on some interval, it should still be present and return zero there. It should not simply disappear from the domain.

Mathematically, if the fiber domain is $[s_0,s_1]$, then each source must be defined for all

$$
s\in[s_0,s_1].
$$

Subroutines exercised:
`Fiber`, `validate_source_coverage`

#### Zero generator gives zero DGD

The `zero_fiber` case sets both bend and twist contributions to zero. Then

$$
K(s)=0,\qquad K_\omega(s)=0,
$$

so the exact solution is

$$
J(s)=I,\qquad G(s)=0,
$$

and therefore

$$
\mathrm{DGD}=0.
$$

Subroutines exercised:
`propagate_fiber_sensitivity`, `output_dgd`

#### Constant synthetic `K_\omega` calibration

The `custom` case uses `propagate_piecewise_sensitivity` directly with

$$
K(s)=0,
\qquad
K_\omega(s)=\frac{i}{2}\sigma_3,
\qquad
L=2.
$$

Since $J(L)=I$, one gets

$$
G(L)=\int_0^L K_\omega\,ds = L\frac{i}{2}\sigma_3 = i\sigma_3,
$$

and therefore

$$
H_{\mathrm{PMD}}
=
-iJ^{-1}G
=
-iG
=
\sigma_3.
$$

The eigenvalues are $\pm 1$, so the DGD is

$$
\mathrm{DGD}=1-(-1)=2.
$$

This is the simplest exact normalization test for `pmd_generator` and `output_dgd`.

Subroutines exercised:
`propagate_piecewise_sensitivity`, `pmd_generator`, `output_dgd`

#### Demo-fiber smoke tests

The `demofiber1()` checks verify that:

- `propagate_fiber` returns a valid $2\times 2$ Jones matrix and interval statistics
- `propagate_fiber_sensitivity` returns valid $(J,G)$ output and finite DGD
- `write_fiber_input_plot3d` produces an output file

This is not a sharp analytic test. It is a full-stack smoke test for the user-facing path.

Subroutines exercised:
`demofiber1`, `propagate_fiber`, `propagate_fiber_sensitivity`, `write_fiber_input_plot3d`

## `test_paddle_transfer.jl`

This file tests the polarization-transfer-function part of the code using fiber paddles interpreted as ideal retarders.

### Test-side helper subroutines

The file introduces:

- `build_paddle_test_fiber(paddles)`
- `propagate_test_state(fiber, input_state)`
- `state_phase_error(actual, expected)`

The helper `build_paddle_test_fiber` converts a list of paddle specifications

$$
(R,\;N,\;\theta)
$$

into a `Fiber` object. The calibration used in the tests is the same one stated in `unit-test-idea.md`:

$$
\delta(R,N)=\frac{\pi}{2}\,N\,\frac{30\ \mathrm{mm}}{R}.
$$

The ideal retarder model is

$$
J(\theta,\delta)
=
R(-\theta)
\begin{pmatrix}
e^{i\delta/2} & 0 \\
0 & e^{-i\delta/2}
\end{pmatrix}
R(\theta),
$$

with

$$
R(\theta)=
\begin{pmatrix}
\cos\theta & -\sin\theta \\
\sin\theta & \cos\theta
\end{pmatrix}.
$$

The helper `state_phase_error` compares Jones states modulo global phase. That is physically necessary because

$$
\psi \sim e^{i\phi}\psi
$$

represents the same polarization state.

### Test set `Paddle transfer cases`

Each named case `1P-*`, `2P-*`, `3P-*`, or `5P-*` checks a known retarder composition. All cases verify:

- `propagate_fiber` returns a $2\times 2$ Jones matrix
- each interval was actually propagated
- norm is preserved
- the output state matches the expected state up to global phase

The core propagation path under test is:

$$
\texttt{exp\_jones\_generator}
\rightarrow
\texttt{exp\_midpoint\_step}
\rightarrow
\texttt{propagate\_interval!}
\rightarrow
\texttt{propagate\_fiber}.
$$

#### One-paddle tests

The `1P-*` cases are the cleanest convention tests.

- `1P-1`: half-wave plate at $45^\circ$, so
  $$
  H \rightarrow V
  $$
- `1P-2`: quarter-wave plate at $45^\circ$, so
  $$
  H \rightarrow R
  $$
- `1P-3`: half-wave plate aligned with the $H/V$ axes, so
  $$
  D \rightarrow A
  $$
- `1P-4`: a fractional-retardance case producing a nontrivial elliptic state

These are mostly tests of sign, angle, and retardance convention.

#### Two-paddle tests

The `2P-*` cases probe composition:

- `2P-1`: two quarter-wave plates on the same axis add to a half-wave plate
- `2P-2`: quarter-wave plus half-wave on different axes checks noncommuting composition
- `2P-3`: analytic two-stage mapping of $D$
- `2P-4`: circular input passing through a mixed half-wave sequence and returning to the same circular state

The key physical point is that retarder compositions generally do not commute:

$$
J_2J_1 \neq J_1J_2.
$$

So these tests are sensitive to multiplication order.

#### Three-paddle tests

The `3P-*` cases extend the same idea to longer compositions:

- `3P-1`: QWP-HWP-QWP identity on $H$
- `3P-2`: analogous identity on $D$
- `3P-3`, `3P-4`: harder arbitrary compositions with known outputs

Identity-style cases are especially effective for exposing angle-sign mistakes because a small convention error destroys the exact return-to-input behavior.

#### Five-paddle tests

The `5P-*` cases are long-composition regression tests:

- `5P-1`: grouped same-axis retardances add to an effective half-wave plate
- `5P-2`: obvious identity sub-blocks multiply to the identity
- `5P-3`, `5P-4`, `5P-5`: longer arbitrary analytic compositions

These are useful because small sign or ordering defects can accumulate and only become obvious in longer sequences.

## DGD tests (`test_path_integral.jl`)

These test sets live in `test_path_integral.jl` alongside the path-integral tests. They exercise the DGD channel directly and are the most analytic part of that file.

### Test-side helper subroutines

The test file defines:

- `piecewise_constant_matrix(breaks, values)`
- `propagate_dgd_case(K, Komega, breaks)`
- `accepted_steps(stats)`

The point of these helpers is to test the DGD solver directly with analytically prescribed matrices, rather than only through bend and twist physics.

### Test set `DGD exact cases`

These are the highest-value DGD tests because they have closed-form answers.

#### `Zero Komega gives zero DGD`

Here $K(s)$ is nontrivial but

$$
K_\omega(s)=0
\qquad
\forall s.
$$

Then the sensitivity equation reduces to

$$
\frac{dG}{ds}=KG,\qquad G(0)=0,
$$

whose exact solution is

$$
G(s)=0.
$$

So the test correctly expects:

- $J$ is nontrivial
- $G=0$
- $\mathrm{DGD}=0$

This is the cleanest check that DGD is not spuriously generated by ordinary Jones evolution.

Subroutines exercised:
`propagate_piecewise_sensitivity`, `output_dgd`

#### `Closed form along sigma_z`

This is the canonical commuting benchmark for the DGD solver:

$$
K(\omega)=i\,\beta(\omega)\,\sigma_3,
\qquad
K_\omega(\omega)=i\,\beta'(\omega)\,\sigma_3.
$$

Since everything commutes,

$$
J(L)=e^{\,i\beta L\sigma_3},
$$

and

$$
G(L)=\partial_\omega J(L)=iL\beta' \sigma_3\,J(L).
$$

Therefore

$$
H_{\mathrm{PMD}}
=
-iJ^{-1}G
=
L\beta'\sigma_3,
$$

so the DGD is the eigenvalue separation

$$
\mathrm{DGD}=2|L\beta'|.
$$

This test fixes the code's normalization convention for both `pmd_generator` and `output_dgd`.

Subroutines exercised:
`propagate_piecewise_sensitivity`, `pmd_generator`, `output_dgd`

#### `Closed form along arbitrary axis`

This repeats the previous test with

$$
K(\omega)=i\,\beta(\omega)\,\hat n\cdot\vec{\sigma},
\qquad
K_\omega(\omega)=i\,\beta'(\omega)\,\hat n\cdot\vec{\sigma},
$$

for a fixed unit vector

$$
\hat n=\frac{1}{\sqrt 3}(1,1,1).
$$

The exact scalar DGD is still

$$
\mathrm{DGD}=2|L\beta'|.
$$

The purpose is to show that DGD is not tied to the special $H/V$ basis.

Subroutines exercised:
`propagate_piecewise_sensitivity`, `pmd_generator`, `output_dgd`

#### `Piecewise cancellation on one axis`

Two equal-length segments are placed on the same axis with equal and opposite sensitivity:

$$
K_{1,\omega}=ia\sigma_3,\qquad K_{2,\omega}=-ia\sigma_3.
$$

Then the exact accumulated PMD generator is

$$
H_{\mathrm{PMD}}
=
(L_1a-L_2a)\sigma_3.
$$

With $L_1=L_2$, this is

$$
H_{\mathrm{PMD}}=0,
\qquad
\mathrm{DGD}=0.
$$

This is more stringent than the $K_\omega=0$ case because the sensitivity channel is active, and only correct piecewise accumulation makes the cancellation happen.

Subroutines exercised:
`piecewise_constant_matrix`, `propagate_piecewise_sensitivity`, `output_dgd`

#### `Piecewise commuting accumulation`

For commuting segments on a common axis,

$$
K_j=i\beta_j\sigma_3,
\qquad
K_{\omega,j}=ia_j\sigma_3,
$$

the exact PMD generator is

$$
H_{\mathrm{PMD}}
=
\left(\sum_j L_j a_j\right)\sigma_3,
$$

and therefore

$$
\mathrm{DGD}
=
2\left|\sum_j L_j a_j\right|.
$$

This is the sharpest breakpoint-accumulation test in the DGD suite.

Subroutines exercised:
`piecewise_constant_matrix`, `propagate_piecewise_sensitivity`, `pmd_generator`, `output_dgd`

### Test set `DGD invariants and regression`

These tests do not always rely on a unique closed form, but they express exact invariants or controlled compatibility checks.

#### `Constant basis rotation preserves DGD`

If

$$
K'(s)=QK(s)Q^{-1},
\qquad
K'_\omega(s)=QK_\omega(s)Q^{-1},
$$

for constant unitary $Q$, then

$$
J'(L)=QJ(L)Q^{-1},
\qquad
G'(L)=QG(L)Q^{-1},
$$

and therefore

$$
H'_{\mathrm{PMD}}
=
QH_{\mathrm{PMD}}Q^{-1}.
$$

So the eigenvalue spectrum, and hence the DGD, is unchanged.

This is the natural basis-invariance test for any PMD/DGD implementation.

Subroutines exercised:
`exp_jones_generator`, `propagate_piecewise_sensitivity`, `pmd_generator`, `output_dgd`

#### `Phase scaling leaves PMD and DGD unchanged`

If one multiplies the final outputs by a common scalar phase,

$$
J\rightarrow e^{i\phi}J,
\qquad
G\rightarrow e^{i\phi}G,
$$

then

$$
J^{-1}G \rightarrow J^{-1}G.
$$

So both `pmd_generator(J,G)` and `output_dgd(J,G)` must be invariant.

This is a compact exact check that the final observables are phase-insensitive.

Subroutines exercised:
`pmd_generator`, `output_dgd`

#### `Noncommuting K with zero Komega still gives zero DGD`

This case uses piecewise changes of birefringence axis, so

$$
[K(s_1),K(s_2)]\neq 0
$$

in general, while still enforcing

$$
K_\omega(s)=0.
$$

The exact prediction remains

$$
G(s)=0,
\qquad
\mathrm{DGD}=0.
$$

This is important because it proves that noncommuting Jones evolution alone does not pollute the DGD channel.

Subroutines exercised:
`propagate_piecewise_sensitivity`, `output_dgd`

#### `Direct sensitivity agrees with finite difference`

This is the compatibility test between the direct sensitivity solve and the older finite-difference idea.

The direct method computes

$$
G_{\mathrm{direct}}(L)
$$

from the coupled system. The finite-difference reference computes

$$
G_{\mathrm{FD}}(L)
\approx
\frac{J(L,\omega+\Delta\omega)-J(L,\omega-\Delta\omega)}{2\Delta\omega},
$$

for a small but not excessively small $\Delta\omega$.

The test then compares:

- $G_{\mathrm{direct}}$ against $G_{\mathrm{FD}}$
- the PMD generators built from each
- the resulting DGD values

This is a regression test, not the fundamental truth oracle. The direct sensitivity method is the one the code is designed to trust; finite differencing is used here as a controlled cross-check in a benign regime.

Subroutines exercised:
`propagate_piecewise_sensitivity`, `propagate_piecewise`, `pmd_generator`, `output_dgd`

### Test set `DGD fiber integration`

#### `Automatic breakpoint union matches explicit partition`

This test returns to the actual `Fiber` abstraction with a bend source and a twist source carrying different breakpoint sets. It compares

$$
\texttt{propagate\_fiber\_sensitivity(fiber)}
$$

against

$$
\texttt{propagate\_piecewise\_sensitivity(}
\texttt{make\_generator(fiber),}
\texttt{make\_generator\_omega(fiber),}
\texttt{explicit\_union)}
$$

on the manually supplied union of all discontinuities.

The expected identity is:

$$
\{ \text{automatic breakpoints} \}
=
\{ \text{explicit union} \},
$$

and therefore the two propagations should give identical $(J,G)$ and identical DGD.

This is the most direct test of the fiber-level convenience wrapper for DGD propagation.

Subroutines exercised:
`BendSource`, `TwistSource`, `Fiber`, `fiber_breakpoints`, `generator_K`, `generator_Kω`, `propagate_fiber_sensitivity`, `propagate_piecewise_sensitivity`

## What these tests establish, and what they do not

The strongest tests in this suite are the ones with exact analytic answers or exact invariance statements. They establish:

- correct source additivity
- correct breakpoint bookkeeping
- correct Jones-matrix multiplication order
- correct normalization of `pmd_generator` and `output_dgd`
- absence of spurious DGD when $K_\omega=0$
- agreement of the direct sensitivity solver with a finite-difference reference in a regime where the latter is trustworthy

They do not, by themselves, prove that every bend or twist constitutive law in the model is the final law for every real fiber experiment. A unit test can show that the solver reproduces a chosen model correctly. It cannot, by itself, prove that the chosen model is the unique correct physical model of the laboratory system.
