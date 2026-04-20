This list of tasks helps keep humans organized but also
provides agents with an idea of upcoming features. Agents
should not start work on these without explicit authorization. 

- [ ] add struct in path-integral.jl that reflects simulation parameters 
  (eg `rtol`, `atol`, step controls)
  
- [x] Update AGENT.md to provide guidance on how to structure unit
  tests. Focus on clearly deliniating high-value physics motivated tests, 
  test involving validation data, tests that accompany select simulation 
  test runs and finally tests that help keep agents from wandering 

- [x] Is there a julia convention for building a 3D path one bit at a time? 
  I guess this will likey involve the Frenet-Serret frame. I'm open to the idea 
  that there's already a convention for this. In that case I'd use multiple 
  dispatch to express the detailed intent for my fiber based application.

- [ ] In light of the recent addition of ARTHITECTURE.md and AGENT.md 
  and README.md is there any refactoring that should take place? Are there any
  consequential contradictions between these new files and the existing 
  organization, logic, represented physics or physics modeling in the codebase?

- [ ] Move the portion of fiber-path-plot.jl that relates to the 3D
  geometry of the fiber to to path-geometry-plot.jl. The 3D plotting 
  should be informed by the geometries and geometric calculations in 
  path-geometry.jl. It should include a movable plane that is moved by the 
  mouse along the path length from start to finish. At each intermediate point 
  the fernet frame coordinates should be graphically illustrated in a cross-
  section plot. Nothing in path-geometry-plot.jl should relate to optical
  fiber. Create some simple example paths in demo.jl that illustrate 
  path-geometry-plot.jl. 

  - [ ] test the refactoring from Sunday night related to above

- [ ] Move the portion of fiber-path-plot.jl related to plotting the 
  Poincare sphere and move it into a new file called poincare-sphere-plot.jl.

- [ ] Restructure fiber-path-plot.jl to draw upon on path-geometry-plot.jl
  and poincare-sphere-plot.jl while still retaining its overall goals. That is, illustrating the transformation of an input polarization state as a function of distance s along the length of an optical fiber. 


