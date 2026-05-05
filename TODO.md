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

- [ ] Is this what we want? Piecewise bend! loops don't accumulate geometric twist in total_material_twist — but they shouldn't, because a BendSegment has geometric_torsion = 0 (circular arcs have zero torsion). A helix does have nonzero geometric_torsion, but that's captured in geometric_torsion(seg, s), not in total_material_twist.

- [ ] Update demo.jl to use physical birefringences if possible.

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


- [ ] TODO fix the MCM demo in demo.jl 4/28 task

- [ ] Verify the `Twist` `is_continuous = true` carry-over semantics are what we
      want: when a `Twist` meta has `is_continuous = true`, the resolver in
      path-geometry.jl computes its `phi_0` as `prev.phi_0 + ∫_0^{prev_run_length}
      prev.rate(s_local) ds_local` — i.e. the absolute phase at the start of the
      new run equals the accumulated phase at the end of the prior run. Confirm
      this is the intended physical meaning before any downstream consumer
      (polarization propagator, plotting, etc.) starts relying on `phi_0`.

- [ ] Add a `resolved::Bool` flag (or similar mechanism) to `AbstractMeta`
      so the system can certify that all meta on a `SubpathBuilt` /
      `PathBuilt` has been fully processed. Some meta is interpreted by code
      internal to `path-geometry*.jl` (e.g. `Twist`); other meta is
      interpreted by external code (e.g. `MCMadd`/`MCMmul` in
      `fiber-path-modify.jl`). A `resolved` flag would let downstream
      consumers (or a `complete(::SubpathBuilt)` / `complete(::PathBuilt)`
      check) verify that no meta has been silently ignored before a built
      structure is treated as authoritative.

#########################
#### path-geometry ######
#########################

- [ ] I want to make some modifications that focus on path-geometry.jl. Let's worry about the downstream consequences of these changes later. 

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

##################################
### true Subpath-local storage ###
##################################

- [ ] TODO20260505 PlacedSegment.origin and PlacedSegment.frame are stored in 
  global coordinates, evolving forward as segments are added.   
  The "first and last point anchored globally" part of your 
  vision is correct, and "coordinate frame evolving per segment"
   is correct — the frame on segment N is the global frame at
  the start of segment N. But the storage is global, not
  Subpath-local. The path just happens to start at global
  (0,0,0) because there's no Subpath start_point concept yet.

  What the new plan currently says                              
   
  Same as today, except pos starts at collect(sub.start_point)  
  (still global, just anchored at the Subpath's start_point 
  instead of (0,0,0)). PlacedSegments still hold global         
  origin/frame.                                             

  What "true Subpath-local storage" would mean

  PlacedSegment.origin would be the offset from start_point     
  (Subpath-local), and PlacedSegment.frame would be a rotation
  matrix whose columns are local N/B/T axes (with local +z =    
  start_outgoing_tangent). Global queries would compute     
  start_point + R_start · ps.origin + R_start · ps.frame ·
  v_segment_local, where R_start is the rotation taking
  Subpath-local axes to global.

  This is genuinely stronger independence: a SubpathBuilt could 
  be "relocated" by mutating its Subpath's start_point without
  rebuilding, and the internal data is invariant to where in    
  space the Subpath sits. The cost is more careful build() math
  (rotating jumpto_point into local before solving the
  connector) and a query layer that does R_start · ...
  everywhere.