**path-integral approach to birefringence propagation**

J. Britton, 4/13/2026 

---

These julia files describe a different approach to calculating the end-to-end polarization 
transfer function (PTF) of an optical fiber link. It's based on a Magnus-type Lie-group integrator which automatically subdivides the calculation to keep errors bounded.

# current approach
The approach in BIFROST paper calculates
$$J_{\text{total}}=\prod_{i}J_{i}$$
with the matrix order matching the order light encounters along the length of the fiber. There is no consideration of the non-commutative behavior of eg linear birefringence and fiber twist. Therefore it's hard to put a bound on the PTF error. 

# a new approach using generator $K(s)$
This approach is refactored around a generator $K(s)$, not around pre-sliced segments. The generator is defined by

$$\frac{dJ}{ds}=K(s)\,J,\qquad J(0)=I. \tag{1}$$

The design is
- The physical inputs describing the fiber must be available as functions so that the fiber is described parametrically (not just on a coarse sampled grid). 
- Each birefringence mechanism returns a local generator contribution $K_m(s)$.
- The solver assembles $K(s)=\sum_m K_m(s)$.
- The physical model is organized as a `Fiber` made from typed birefringence sources, rather than one monolithic fiber-input struct.
- Each source carries its own parametric description and declares the points where its behavior can change discontinuously.
- The fiber computes the global breakpoint set automatically by taking the union of the breakpoints declared by its sources.
- The propagator advances $J$ adaptively over smooth intervals. 
- The interval size is selected dynamically to keep error below a specified threshold.
    - smooth fiber → large steps,
    - rapidly changing bend/twist → small steps,
    - discontinuities → explicit interval boundaries.
    - Net: The computational effort is concentrated where the physics actually varies.


## simple analytic case
If the length of a fiber segment is short, its Jones matrix is close to identity $J\approx I+K(s)\,\Delta s,$ where $K(s)$ is the local generator of polarization evolution.

$$ J(L)\approx\prod_{i=1}^{N}\bigl(I+K(s_{i})\Delta s\bigr)$$

If the matrices commute you can integrate analytically. $$\prod_{i}\bigl(I+K(s_{i})\Delta s\bigr)\to\exp\!\left(\int_{0}^{L}K(s)\,ds\right)$$
 

## numerical solving

The most natural framework is a Magnus-type Lie-group integrator. These methods are built for non-autonomous linear matrix ODEs and propagate by products of exponentials, such as Jones matrices. The simplest useful version is the exponential midpoint rule 

$$J_{n+1}=e^{h\,K(s_{n}+h/2)}J_{n}$$

where $h$ is a small step. 

> **INSIGHT** The use of this exponential is motivated by a simple geometric idea. The solution of a constant-coefficient ODE over one interval is exactly an exponential. By using exponentials step by step, the method respects the multiplicative structure of Jones propagation rather than approximating it as an additive update. If the fiber is lossless and you factor out common phase, $K$ lives in the Lie algebra of $SU(2)$, so exponential-based propagation preserves the physical structure much better than a generic entrywise ODE solver.

Use a posteriori error estimator for adaptive step-size control, such as adaptive step-doubling. Adaptive solvers do this automatically (eg RK45). For exponential midpoint, the simplest is:
- one full step of size $h$,
- two half steps of size $h/2$,
- compare the results,
- accept/reject and update $h$.

There is some additional detail to the algorithm. 

- Because the error behaves roughly like $h^3$, the next step size is chosen with a cubic-root controller:
    
$$h_{\text{new}} \sim h \left(\frac{\text{tol}}{\text{err}}\right)^{1/3}$$

- In the interest of accuracy the adaptive mechanism for choosing $h$ is not permitted to cross discontinuities in $K$. For example, if the bend radius jumps from one value to another, the algorithm integrates on each side and then combines. 

In the current code this breakpoint logic is part of the fiber-specification layer rather than the propagator itself. Each birefringence source declares its own breakpoints, and the assembled fiber computes

$$\text{breakpoints}(F)=\{s_{\mathrm{start}},s_{\mathrm{end}}\}\cup\bigcup_m \text{breakpoints}(m).$$

The propagator then integrates separately on each interval between consecutive breakpoints.

> **INSIGHT** It turns out there is a closed-form exponential for a $2×2$ Jones generator. More generally, for an exactly traceless matrix $A$ the Cayley-Hamilton theorem gives $A^2 = -\det(A) I$, so
> $$\exp(A) = \cosh(μ) I + \sinh(μ)/μ * A$$
> where  $μ^2 = -\det(A)$. This is exact, not an approximation. See Appendix A.1.

> **INSIGHT** Jones matrices are physically meaningful only up to an overall global phase in many applications. Two matrices differing only by a scalar phase factor often represent the same observable polarization action. So the step-doubling comparison should not use a raw matrix difference alone. Instead, align the phases first:
> $$\text{err} = \|A-e^{i\phi}B\|,$$
>where $e^{i\phi}$ is chosen to best match $A$ and $B$.
>
>This prevents the adaptive controller from wasting effort on irrelevant common-phase differences. It's implemented in `phase_insensitive_error(A, B)`.


# Integrate $\partial_\omega J$ directly for DGD

For a polarization maintaining fiber the differential group delay (DGD) is a measure of the difference in transit time for light launched into the fast axis and light launched in the slow axis. It can be shown that DGD can be written as 

$$ G(s, \omega) \equiv  \frac{\partial J(s, \omega)}{\partial \omega}$$

The BIFROST paper’s DGD formulation is Eq (19) (BIFROST `calcDGD()`) 

$$\partial_\omega J \approx \frac{J(\omega + \Delta\omega) - J(\omega)}{\Delta\omega}.$$

This relies on an additional numerical parameter $\Delta \omega$ and the finite difference is subject to round-off error. 

An alternate approach removes this avoidable source of numerical error and puts DGD under the same adaptive tolerance as the main propagator. Recall the definition of the propagator where we make the frequency dependence explicit. 

$$\frac{dJ}{ds}=K(s,\omega)\,J,\qquad J(0,\omega)=I. $$

Differentiate with respect to $\omega$. We get a coupled system

$$\frac{dJ}{ds} = KJ,\quad J(0) = I,\tag{2}$$

$$\frac{dG}{ds} = K_\omega J + KG,\quad G(0) = 0. \tag{3}$$

Then form $J^{-1}G$ directly at the end. 




## Code overview

The code is now split into two conceptual layers.

### Fiber specification layer

`fiber-path.jl` defines how a fiber is described before any propagation happens.

- `PiecewiseProfile` is still the basic helper for piecewise-parametric scalar functions of $s$.
- `BendSource` and `TwistSource` are concrete birefringence-source types.
- `Fiber` is an assembled object with a domain $[s_{\mathrm{start}},s_{\mathrm{end}}]$ and a list of sources.
- Each source provides:
  - a local contribution to $K(s)$,
  - a local contribution to $K_\omega(s)$,
  - a declared breakpoint set,
  - explicit coverage intervals in $s$.

The explicit coverage is important. A source is not allowed to be undefined on part of the fiber and silently treated as zero there. Instead, every source must cover the full fiber domain. If a source is physically inactive on some interval, it must still be defined there and return a zero contribution. The `Fiber` constructor checks this and throws an error if any source has a gap in its coverage.

At the assembled-fiber level,

$$K(s)=\sum_m K_m(s), \qquad K_\omega(s)=\sum_m K_{\omega,m}(s).$$

The fiber also computes the global breakpoint set automatically from its sources.

### Propagation layer

`path-integral.jl` sits on top of that specification layer and implements the actual propagators.

- `make_generator(f)` assembles the total local generator $K(s)$ from the sources stored in a `Fiber`.
- `make_generator_omega(f)` assembles the total frequency derivative $K_\omega(s)$ from the same sources.
- `exp_jones_generator(A)` computes $\exp(A)$ for a $2\times2$ generator using the closed-form Cayley-Hamilton formula.
- `exp_midpoint_step(K, s, h, J)` is one exponential-midpoint step for the Jones propagator.
- `propagate_interval!()` is the adaptive controller on one smooth interval. It takes one full step and two half steps, compares them with `phase_insensitive_error`, accepts or rejects the step, and updates $h$ with the cubic-root rule $(\text{tol}/\text{err})^{1/3}$.
- `propagate_piecewise()` performs the same propagation on a prescribed breakpoint partition.
- `propagate_fiber()` is the fiber-level convenience wrapper: it computes the breakpoints from the `Fiber` automatically and then calls the piecewise propagator.

### DGD extension
The original code was built around a very specialized solver for one 2x2 matrix ODE.
Once you add the pair of equations (2)-(3), this becomes a coupled block system. It requires extending the machinery for $J$ to now include both $J$ and $G$. This is implemented as
- `exp_sensitivity_midpoint_step()` propagates the coupled ($J$, $G$) system over one midpoint step
- `propagate_interval_sensitivity!()` adaptive step-doubling over smooth intervals
- `propagate_piecewise_sensitivity()` adaptive step-doubling over breakpoints
- `propagate_fiber_sensitivity()` computes the breakpoints from the `Fiber`, assembles both $K(s)$ and $K_\omega(s)$, and propagates the coupled system
- `pmd_generator(J, G)` forms -im * J^{-1}G
- `output_dgd(J, G)` extracts the DGD 

At the API level, the important change is that $K_\omega(s)$ is now part of the source abstraction. Each source has a method for its $K(s)$ contribution and a second method for its $K_\omega(s)$ contribution. Right now some source types may legitimately return the zero matrix for $K_\omega$ when their frequency dependence has not yet been modeled, but the interface is there and the DGD propagator always has a well-defined assembled $K_\omega(s)$ to use.

# Example path-demo.jl
`path-demo.jl` is just a thin driver on top of that stack. It builds a `Fiber` from a bend source and a twist source, then calls `propagate_fiber()` to get the final Jones matrix and per-interval stats. The demo no longer carries a separate manually maintained `breaks` array through the propagation API; the breakpoints are derived automatically from the fiber's sources.


## Simplified physics encoded in the generator K(s)

*Bending → linear birefringence*

For a bend of radius $R(s)$, the fiber experiences stress-induced birefringence with magnitude proportional to $1/R(s)^2$ in the simplest bending model. In the paper, the bend-induced birefringence is given by Eq. (9).  ￼The bend axis has an orientation in the transverse plane. If the bend angle is $\theta_b(s)$ in turns, the corresponding physical axis angle is

$$\phi_b(s)=2\pi\,\theta_b(s).$$

Linear birefringence aligned with this axis gives a generator proportional to a rotated Pauli-matrix combination. The double-angle dependence appears because Jones matrices represent polarization axes modulo $180^\circ$, not $360^\circ$.

*Twisting → circular birefringence*

Twisting creates circular birefringence proportional to the local twist rate $\tau(s)$. In the paper, the left-right circular propagation-constant difference is linear in $\tau$ [Eq. (11)], and the corresponding Jones matrix has the circular-birefringence form shown in Eq. (16).  ￼

If $\text{twist}(s)$ is the accumulated twist angle in turns, then the solver really needs
$\tau(s)=2\pi \frac{d}{ds}\text{twist}(s)$,
in rad/m.

**Combined generator**

The total local generator is the sum

$$K(s)=K_{\text{bend}}(s)+K_{\text{twist}}(s).$$

Because the bending and twisting terms generally point along different Pauli directions, they do not commute:
$[K_{\text{bend}}(s),K_{\text{twist}}(s)]\neq 0$. That noncommutation is why one cannot usually collapse the whole problem into a single scalar integral.

**Fiber representation**

The adaptive solver evaluates the fiber at arbitrary points such as $s+h/2$, $s+h/4$, and so on. Therefore the physical inputs must be available as functions, not just a coarse sampled grid.

The clean representation is now source-based rather than one large input struct.

- A `BendSource` owns the bend-specific parametrization and strength law.
- A `TwistSource` owns the twist-specific parametrization and strength law.
- A `Fiber` assembles those sources into one physical fiber model over a specified interval in $s$.

For the present bend/twist model, the bend source is naturally described by

	•	$R_b(s)$: bend radius
	•	$\theta_b(s)$: bend-axis orientation

and the twist source by

	•	$d(\text{twist})/ds$: twist gradient

These are then converted internally into solver-friendly quantities:

$$\kappa_x(s)=\frac{\cos(2\pi\theta_b(s))}{R_b(s)},\qquad
\kappa_y(s)=\frac{\sin(2\pi\theta_b(s))}{R_b(s)},$$

and

$$\tau(s)=2\pi\,\frac{d}{ds}\text{twist}(s).$$

This formulation is numerically better because it avoids angle singularities when the bend goes to zero and makes straight fiber simply $\kappa_x=\kappa_y=0$.

The same source abstraction also carries the data needed for DGD. In addition to its contribution to $K(s)$, each source also contributes its piece of $K_\omega(s)$. That keeps the PMD/DGD machinery aligned with the ordinary Jones propagation machinery: the fiber assembles both objects from the same list of sources, and the sensitivity propagator uses the same breakpoint structure and the same adaptive error-control strategy.


# APPENDICES

## A1 Cayley–Hamilton

For any $2\times2$ matrix $A$, Cayley–Hamilton gives
$$
A^2-(\operatorname{tr}A)A+(\det A)I=0.
$$
If $A$ is exactly traceless, $\operatorname{tr}A=0$, so
$$
A^2=-(\det A)I.
$$
Define
$$
\mu^2=-\det A.
$$
Then
$$
A^2=\mu^2 I,\qquad A^3=\mu^2 A,\qquad A^4=\mu^4 I,\ \dots
$$
So the exponential series splits into even and odd powers:
$$
e^A=\sum_{n=0}^\infty \frac{A^n}{n!}
=\sum_{k=0}^\infty \frac{A^{2k}}{(2k)!}
+\sum_{k=0}^\infty \frac{A^{2k+1}}{(2k+1)!}.
$$
Using $A^{2k}=\mu^{2k}I$ and $A^{2k+1}=\mu^{2k}A$,
$$
e^A=
\left(\sum_{k=0}^\infty \frac{\mu^{2k}}{(2k)!}\right)I
+\left(\sum_{k=0}^\infty \frac{\mu^{2k}}{(2k+1)!}\right)A.
$$
Recognizing the series,
$$
\sum_{k=0}^\infty \frac{\mu^{2k}}{(2k)!}=\cosh\mu,
\qquad
\sum_{k=0}^\infty \frac{\mu^{2k}}{(2k+1)!}=\frac{\sinh\mu}{\mu},
$$
we get
$$
e^A=\cosh(\mu),I+\frac{\sinh(\mu)}{\mu},A,
\qquad \mu^2=-\det A.
$$

At $\mu=0$, interpret $\sinh(\mu)/\mu\to1$, so
$$
e^A=I+A.
$$
