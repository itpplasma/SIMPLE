function znew = neo_orb_timestep(z)
% Fortran routine timestep in neo_orb.f90
    s = libpointer('doublePtr', z(1));
    th = libpointer('doublePtr', z(2));
    ph = libpointer('doublePtr', z(3));
    lam = libpointer('doublePtr', z(4));
    ierr = libpointer('int32Ptr', 0);
    calllib('libneo_orb', 'neo_orb_global_MOD_timestep', s, th, ph, lam, ierr);
    znew(1) = s.Value;
    znew(2) = th.Value;
    znew(3) = ph.Value;
    znew(4) = lam.Value;
end
