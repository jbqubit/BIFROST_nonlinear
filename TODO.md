This list of tasks helps keep humans organized but also provides agents with an idea of
upcoming features. Agents should not start work on these without explicit authorization.

- [ ] add struct in path-integral.jl that reflects simulation parameters (eg `rtol`,
      `atol`, step controls)

- [ ] Resume work on branch benchmark-integrate. 

- [ ] Finish ensuring MCM compatability across the rest of the codebase.

- [ ] Add unit tests for fiber paddles that reproduce Thorlabs website. 

- [x] Update AGENT.md to provide guidance on how to structure unit tests. Focus on
      clearly deliniating high-value physics motivated tests, test involving validation
      data, tests that accompany select simulation test runs and finally tests that help
      keep agents from wandering

- [x] Is there a julia convention for building a 3D path one bit at a time? I guess this
      will likey involve the Frenet-Serret frame. I'm open to the idea that there's
      already a convention for this. In that case I'd use multiple dispatch to express
      the detailed intent for my fiber based application.

- [x] The method uncovered_intervals in path-geometry.jl seems redundant since fiber segments 
  are built piecewise. My construction it's not possible for there to be gaps. Don't do this.... they are needed. 

- [ ] `HermiteConnector`, `JumpBy`, and `JumpTo` are presently
  incompatible with MCM. 

- [ ] implement for T in (JumpBy, JumpTo)
  - [ ] properly implement sample_path for these

- [ ] Remove the "Local Frenet section" subplot in path-geometry.plot.jl.

- [ ] In light of the recent addition of ARTHITECTURE.md and AGENT.md and README.md is
      there any refactoring that should take place? Are there any consequential
      contradictions between these new files and the existing organization, logic,
      represented physics or physics modeling in the codebase?

- [ ] Move the portion of fiber-path-plot.jl that relates to the 3D geometry of the fiber
      to to path-geometry-plot.jl. The 3D plotting should be informed by the geometries
      and geometric calculations in path-geometry.jl. It should include a movable plane
      that is moved by the mouse along the path length from start to finish. At each
      intermediate point the fernet frame coordinates should be graphically illustrated
      in a cross- section plot. Nothing in path-geometry-plot.jl should relate to optical
      fiber. Create some simple example paths in demo.jl that illustrate
      path-geometry-plot.jl.

  - [ ] test the refactoring from Sunday night related to above

- [ ] Move the portion of fiber-path-plot.jl related to plotting the Poincare sphere and
      move it into a new file called poincare-sphere-plot.jl.

- [ ] Restructure fiber-path-plot.jl to draw upon on path-geometry-plot.jl and
      poincare-sphere-plot.jl while still retaining its overall goals. That is,
      illustrating the transformation of an input polarization state as a function of
      distance s along the length of an optical fiber.


- [ ] Create a path-geometry.md that documents how it works. Add specific 
  illustratings for important features of path-geometry.jl. Each is described in
  .md and accompanied by code-generated .png that illustrate the points. The
  file defining these illustrations is called path-geometry-illustrated.jl. 
    - [ ] Show that the reference frame for adding additional segments does not 
    rotate with fiber twist.
    - [ ] Illustrate how the axis_angle option is defined and how it
    changes segment addition.
    - [ ] Show how the orientation of the prior segment influences the orientation of a helix and the helix exit path.
    - [ ] Devise examples that illustrate how segments respond to shrinkage and
    contrast it with 



- [ ] Verify the `Twist` `is_continuous = true` carry-over semantics are what we
      want: when a `Twist` meta has `is_continuous = true`, the resolver in
      path-geometry.jl computes its `phi_0` as `prev.phi_0 + ∫_0^{prev_run_length}
      prev.rate(s_local) ds_local` — i.e. the absolute phase at the start of the
      new run equals the accumulated phase at the end of the prior run. Confirm
      this is the intended physical meaning before any downstream consumer
      (polarization propagator, plotting, etc.) starts relying on `phi_0`.

#########################
#### path-geometry ######
#########################

- [x] I want to make some modifications that focus on path-geometry.jl. Let's worry about the downstream consequences of these changes later. 

      Currently TwistOverlay is specified by s_start and length. I want to make changes so that the start and end of each twist is defined with respect to segment boundaries. There are several ways  I can think of implementing this. 

      OPTION 1 :: use meta
      Create a new struct Twist <: AbstractMeta with members
            rate::Function
            \phi_0::Real
            is_continuous::Bool
      In this approach Twist meta is associated with the segment where a particular twist rate commences and continues until the end of the Path or until another segment has an associated Twist meta. The bool is_continuous specifies if the twist phase remains continuous with the prior Twist specification. If is_continuous is True then \phi_0 shouldn't be specified. If is_continuous is False than \phi_0 must be specified as this is the starting phase. The rate::Function must accept \phi as a parameter. 

      OPTION 2 :: Use zero-length Segment
      Create a new struct Twist <: AbstractPathSegment with members
            rate::Function
            \phi_0::Real
            is_continuous::Bool
      In this approach, a Twist segment is inserted into Path placed_segment as a zero-length segment that demarks the start of a particular twist specified by its member data. The same twist rate  continues until the end of the Path or until another Twist segment is added to the Path. The bool is_continuous specifies if the twist phase remains continuous with the prior Twist specification. If is_continuous is True then \phi_0 shouldn't be specified. If is_continuous is False than \phi_0 must be specified as this is the starting phase. The rate::Function must accept \phi as a parameter. 

      One detail common to both approaches is that if is_continuous is True the length of prior segments is important in calculating the phase continuity at the boundary between old and new twist rates. 

      While it's outside the context of the current refactoring it's important to note that properties of individual path segments can be modified using meta

      Please help me think which is the right approach of if there is another that's better.  

- [x] TODO001
      One of the physical aspects not yet properly handled by fiber-path-modify.jl is fiber path
      length due to temperature. Everything works for a series of segments with variable temperature (their length changes) the fiber is rebuilt piecewise segment by segment after the meta is applied. To summarize the desired behavior: if fiber path length L is subject to thermal expansion by factor k the new path is k*L (see eg _modify_segment(seg::StraightSegment)).

      Things break down as soon as add a JumpTo segment is added. Recall that JumpTo 
      specifies a destination which is a global invariant. It's connection to the prior segment is
      resolved at build() time to a QuinticConnector. We also see that the QuinticConnector responds to temperature changes via  
            _modify_segment(seg::QuinticConnector, T_ref, α_lin)
      There's a contradiction here. We want the honor the accumulated pathlength variation due to
      temperature meta prior to JumpTo. This suggests the need for some added constraint that the path solution for JumpTo is also length constrained. See Claude chat "JumpTo
      and meta ::T_K".


- [x] Add new section to the end of demo2.jl titled "meta and JumpTo interplay"


- [x] Apply the principles in AGENT.md to refactor the following in demo2.jl
  - demo_jumpto_2d_min_radius
  - for all 3D demos remove crap like _build_jumpby_variant()
  - 8afdc32951ad6ef1711dfd4ce81844d65504dee3

  - [x] Why is such a lumbering use required for the following. PG.AbstractMeta[PG.Nickname(s)]
404d48ecdabccd3ed92bb3b8b264906b7863c604

  - [ ] Explore what's saved in claude /memory

  - [x] Create a new demo in demo2.jl that illustrates the use of MCM. It should be end-to-end
  showing fiber creation, simulation and summary of the output PTF. The segment
  structure should look like the following.
      straight 5 m, helix of diameter 10 m with 100 turns (pitch 0.05 m), straight 5 m, 
      single helix of diameter 2 m with 4 turns (pitch 0.05m), straight 5 m
      Add temperature meta for the first helix and make it an MCM Particle. Choose 
      the temperature to be 30C +/- 5 C. Illustrate the output by showing degree of linear polarization
      as a function of temperature in an output plot (plus other ways to illustrate the
      temperature dependence of the PTF)

- [ ] One of the questions about this julia port of BIFROST is its speed. Let's come up with some
speed benchmarks and add them to demo3.jl. Extend the two existing MCM demos. Benchmark the compile time and run time of each. Compare with and without StaticParticles for each case. 

- [ ] The current API for fiber meta excludes some of the potential sources
(and influences on) birefringence in fiber-cross-section.jl. Namly 
      - core noncircularity
      - assymetric thermal stress
      - axial tension
      - are there any others? 
Make a plan for adding them as meta (if you agree that's where they belong). Take note that some effects depend on segment temperature.  

- [ ] The current implementation is light on exercise of the DGD path. Make sure that it's working. Are there adequate unit tests? Are there any unit tests where the answers are obvious from simple physics principles? 

- [ ] I want to shift to running all development in docker container that is on a linux box.
In particular I want claude to run locally on that machine either via CLI or with ssh login from my MacBook. I want things sandboxed so that I can give full permission to claude. I'm OK with
giving both claude and the container access to the internet. Use a mounted filesystem
to hold the source code and julia depot. The container should be based on something from nvidia
that supports cuda. Help me draft a docker file and related infrastructure to make this a 
plan come to life. 

- [ ] Systematically look for places where MCM is using minimum, maximum or mean values (or any other coerced values) to handle branches. I want to be cognizant of them. 

- [x] Adjust the color scheme of index*.html to be dark. 

- [ ] Improve api...  MCMadd(:T_K, ΔT) is the existing general mechanism for “segment-local perturbation of named quantity.”

- [x] at the top of index*.html add quick links to all the index files: index1.html, index2.html, etc (all in a single compact line)

- [x] Add some sort of CLI indicator of how far along a particular path path-integral.jl has progressed

- [ ] Consider adding a new layer of abstraction fiber-path-integral.jl that encapsulates the process of taking an extant fiber and performing integration operations. The FiberPathIntegral struct would contain computation specific details like rtol, atol, what type of MCM options should be used and the like. 

- [ ] Add periodic tasks to AGENT.md. Like #2 below. 

- [ ] Fix demo example modify-jumpto-anchor-thermal-2d.html to show path with meta temperature. 

- [ ] VS Code doesn't provide popup hints for path-geometry.jl jumpby!

########## ITEMS FOR TALK ##############
- What does the agent actually have loaded into context at the start of every session? See
chats labeled `AGENTS.md battle`.
- Systematically make passes over README.md, inline examples, inline method documentation to make sure it is up to date. This doesn't seem to happen automatically. 