function neo_orb_init_field(ns, ntp, multharm, integmode)
% Fortran routine init_field in neo_orb.f90
% arguments describe accuracy of splining VMEC fields
% currently file 'wout.nc' is read all the time
% integmode ... -1: RK VMEC, 0: RK CAN, 1: Symplectic Euler
    calllib('libneo_orb', 'neo_orb_MOD_init_field', ...
        libpointer('int32Ptr', ns), ...
        libpointer('int32Ptr', ntp), ...
        libpointer('int32Ptr', multharm), ...
        libpointer('int32Ptr', integmode));
end