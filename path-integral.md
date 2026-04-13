These julia files describe a different approach to calculating the end-to-end polarization transfer function (PTF) of an optical fiber link. It's based on a Magnus-type Lie-group integrator which automatically subdivides the calculation to keep errors bounded.

# current approach
The approach in BIFROST paper calculates
$$J_{\text{total}}=\prod_{i}J_{i}$$
with the matrix order matching the order light encounters along the length of the fiber. There is no consideration of the non-commutative behavior of eg linear birefringence and fiber twist. Therefore it's hard to put a bound on the PTF error. 

# a new approach
This approach is refactored around a generator $K(z)$, not around pre-sliced segments. The generator is defined by

$$\frac{dJ}{dz}=K(z)\,J,\qquad J(0)=I.$$

The design is
- The physical inputs describing the fiber must be available as functions so that the fiber is described parametrically (not just on a coarse sampled grid). 
- Each birefringence mechanism returns a local generator contribution $K_m(z,\lambda)$.
- The solver assembles $K(z)=\sum_m K_m(z)$.
- The propagator advances $J$ adaptively over smooth intervals. 
- The interval size is selected dynamically to keep error below a specified threshold.
    - smooth fiber → large steps,
    - rapidly changing bend/twist → small steps,
    - discontinuities → explicit interval boundaries.
    - Net: The computational effort is concentrated where the physics actually varies.


## simple analytic case
If the length of a fiber segment is short, its Jones matrix is close to identity $J\approx I+K(z)\,\Delta z,$ where $K(z)$ is the local generator of polarization evolution.

$$ J(L)\approx\prod_{i=1}^{N}\bigl(I+K(z_{i})\Delta z\bigr)$$

If the matrices commute you can integrate analytically. $$\prod_{i}\bigl(I+K(z_{i})\Delta z\bigr)\to\exp\!\left(\int_{0}^{L}K(z)\,dz\right)$$
 

## numerical solving

The most natural framework is a Magnus-type Lie-group integrator. These methods are built for non-autonomous linear matrix ODEs and propagate by products of exponentials, such as Jones matrices. The simplest useful version is the exponential midpoint rule 

$$J_{n+1}=e^{h\,K(z_{n}+h/2)}J_{n}$$

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

- [ ] TODO Add explanation of the codebase with reference to the functions in `path-integral.jl`. 

> **INSIGHT** It turns out there is a closed-form exponential for a $2×2$ Jones generator. More generally, for an exactly traceless matrix $A$ the Cayley-Hamilton theorem gives $A^2 = -\det(A) I$, so
> $$\exp(A) = \cosh(μ) I + \sinh(μ)/μ * A$$
> where  $μ^2 = -\det(A)$. This is exact, not an approximation. See Appendix A.1.

- [ ] TODO: An observation by one of the GPTs is the following. Is this implemented in the codbase? What does it even mean? 
 > A small-$|h\mu|$ series fallback avoids catastrophic cancellation:
 > $$\frac{\sinh x}{x}\approx 1+\frac{x^2}{6}+\frac{x^4}{120}$$


> **INSIGHT** Jones matrices are physically meaningful only up to an overall global phase in many applications. Two matrices differing only by a scalar phase factor often represent the same observable polarization action. So the step-doubling comparison should not use a raw matrix difference alone. Instead, align the phases first:
> $$\text{err} = \|A-e^{i\phi}B\|,$$
>where $e^{i\phi}$ is chosen to best match $A$ and $B$.
>
>This prevents the adaptive controller from wasting effort on irrelevant common-phase differences.


## Integrate $\partial_\omega J$ directly for DGD

TODO: Explain this better. 

The paper’s DGD formulation uses $J^{-1}\partial_\omega J$, but both the paper and current code obtain that by finite-differencing wavelength. Instead, propagate the sensitivity matrix $G=\partial_\omega J$ alongside $J$:

$$\frac{dG}{dz}=K_\omega J + K G,\qquad G(0)=0.$$

Then form $J^{-1}G$ directly at the end. That removes the extra finite-difference parameter in wavelength and puts the DGD error under the same adaptive z-tolerance as the propagation itself.  

# Example path-demo.jl
Here's what's in the path-demo.jl file. 

## Simplified physics encoded in the generator K(z)

*Bending → linear birefringence*

For a bend of radius $R(z)$, the fiber experiences stress-induced birefringence with magnitude proportional to $1/R(z)^2$ in the simplest bending model. In the paper, the bend-induced birefringence is given by Eq. (9).  ￼The bend axis has an orientation in the transverse plane. If the bend angle is $\theta_b(z)$ in turns, the corresponding physical axis angle is

$$\phi_b(z)=2\pi\,\theta_b(z).$$

Linear birefringence aligned with this axis gives a generator proportional to a rotated Pauli-matrix combination. The double-angle dependence appears because Jones matrices represent polarization axes modulo $180^\circ$, not $360^\circ$.

*Twisting → circular birefringence*

Twisting creates circular birefringence proportional to the local twist rate $\tau(z)$. In the paper, the left-right circular propagation-constant difference is linear in $\tau$ [Eq. (11)], and the corresponding Jones matrix has the circular-birefringence form shown in Eq. (16).  ￼

If $\text{twist}(z)$ is the accumulated twist angle in turns, then the solver really needs
$\tau(z)=2\pi \frac{d}{dz}\text{twist}(z)$,
in rad/m.

**Combined generator**

The total local generator is the sum

$$K(z)=K_{\text{bend}}(z)+K_{\text{twist}}(z).$$

Because the bending and twisting terms generally point along different Pauli directions, they do not commute:
$[K_{\text{bend}}(z),K_{\text{twist}}(z)]\neq 0$. That noncommutation is why one cannot usually collapse the whole problem into a single scalar integral.

**Fiber representation**

The adaptive solver evaluates the fiber at arbitrary points such as $z+h/2$, $z+h/4$, and so on. Therefore the physical inputs must be available as functions, not just a coarse sampled grid.

The clean representation is:
	•	$R_b(z)$: bend radius
	•	$\theta_b(z)$: bend-axis orientation
	•	$d(\text{twist})/dz$: twist gradient

These are then converted internally into more solver-friendly quantities:

$$\kappa_x(z)=\frac{\cos(2\pi\theta_b(z))}{R_b(z)},\qquad
\kappa_y(z)=\frac{\sin(2\pi\theta_b(z))}{R_b(z)},$$

and

$$\tau(z)=2\pi\,\frac{d}{dz}\text{twist}(z).$$

This formulation is numerically better because it avoids angle singularities when the bend goes to zero and makes straight fiber simply $\kappa_x=\kappa_y=0$.


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