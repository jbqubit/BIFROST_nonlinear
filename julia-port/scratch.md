# scratch 

Several key points should be reiterated at the beginning of the plan. Confirm if you agree.  Let's focus just on modifications to path-geometry*.jl right now and then return in annother pass to build the full refactoring plan. 

Subpaths are constructed using SubpathBuilder. Each Subpath specification begins with start! and ends with jumpto!. It contains only user-supplied information (nothing calculated). Each Subpath is fully independent of all others. 

Each SubpathCached houses a Subpath and contains the results of the following. Starting at the start, for each segment in the Subpath, apply meta and extend the segment forward. When the jumpto is encountered do the curve-fitting needed to extend the path to the endpoint.  Each SubpathCached is fully independent of all others.

PathCached is a container that holds an ordered collection of SubpathCached elements. It's a lightweight organizational element that should feel optional to the enduser. The rest of the system should be made to ingest SubpathCached as the primary input (or bags of them in PathCached). Somehow I want it to be natural for many SubpathCached to be computed in parallel from a bunch of Subpaths... the result gets organized by PathCached. 

I don't think we need the notion of a Path anymore. Do you agree?

Would it make more linguistic sense to call SubpathCached SubpathBuilt and PathCached PathBuilt? This reflects the SubpathBuilder naming. 

As regards your questions...

1. I accept that every demo/test gets a ceremonial jumpto! It also gets a SubpathBuilder, Subpath and SubpathBuilt.

2. I don't understand what's going on with regard to _modified_rebuild. Give me some background so we can discuss further. 

3. The intent is that all other parts of the codebase should be focused on ingesting SubpathBuilt (or bags of them in PathBuilt). All downstream consumers should be focused on using SubpathBuilt as the primary input. They should look inside SubpathBuilt if they need to know anything about start, end, breakpoints, etc.

4. PathBuilt should be a high level bag that holds many SubpathBuilt. 

5. Twist is_continuous requires some glue logic at the interface between subsequent Subpaths. 
I want for a twist specification to extend no longer than a single Subpath. That is, each Subpath needs its own twist (if any). This is part of making each Subpath fully independent of all others. 

6. - The thermal-expansion-on-JumpTo behavior described in the existing
  docstring is being deleted, not refactored.  Agreed. 

7. Eliminate meta for the on the terminal's connector. That just adds too much complexity. 

8. OK

9. Good point. The place where I want things to be parallel is inside propagate_fiber(). The point of making all the Subpath independent is in support of this parallization. This needs to be made clearer in the motivation for this refactor. 

10. JumpBy should be demoted to just-another segment. Don't treat it special. 

11. It may in time. This could be used for fiber plotting in the future. 

12. I agree

13. Yes.

1 Phase A... I agree

2. Phase B ... I agree

3. Phase C .... 





-----

Regarding _modified_rebuild...

"A Fiber is built from PathSpecBuilt." No, A Fiber should be built from Subpath parts. Internal to Fiber is a PathBuilt that contains placed SubpathBuilt parts. Ideally we view SubpathBuilt and PathBuilt as internal with respect to the Fiber implementation. 
- We are doing away with PathSpec and PathSpecCached. 
- We can discuss this more when we get to the refactoring pass.

Regarding the new Plan for A.
- Go ahead and also modify test_path_geometry* to account for changes in syntax. Work hard to ensure that the underlying logical checks are not trivially modified to pass.
- "Once built, a Subpath's geometry can be queried in isolation." Yes, attempts to query geometry of path related aspects of Subpath should fail (as its not yet built).
- Another build(spb::SubpathBuilt) does whatever stitching is needed on the way to populating PathBuilt.
- Although PathBuilt is intended to be a bag it should still be a distinct struct. 
----

Rgarding, build(spb::SubpathBuilt)... It may need to do more to construct PathBuilt. This is the first layer of the stack that know about the overall ordering of the Subpath elements. All the methods that For example, stitching together twist adjacent SubpathBuilt. Make note of that and let's come back to it. Perhaps this can be delegated back to the PathBuilt level? 





