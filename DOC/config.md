* For symplectic Euler, `npoiper2=256` should be set in particular if all orbits
  (including strongly passing) are traced. This is especially the case when they can collide into this region if `sw_coll=.True.`. If one sets a limit via
  e.g. `contr_pp=-1d0` to trace only mildly passing orbits, or even switch
  them off via `notrace_passing=.True.`, one may resort to 128, or even 64
  in extreme cases with short tracing time. However, also for these cases,
  `npoiper2=256` appears to be a safer choice.
