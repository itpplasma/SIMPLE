function neo_orb_init_params(Z_charge, A_mass, E0, dtau, dtaumax, relerr)
% Fortran routine init_params in neo_orb.f90
% Initialized particle parameters for ions (electrons not supported yet)
% Z_charge  ... Charge number of ions
% A_mass    ... Atomic mass units of ions 
% E0        ... Kinetic/thermal energy in electron volts
% dtau      ... Normalized time step, tau = sqrt(2*E0/m)*t
%   For Runge-Kutta (RK) time-steps are adaptively refined with results 
%   supplied in steps of dtau. Symplectic integrators have fixed step dtau.
% dtaumax   ... Maximum allowed time step for Runge-Kutta (default: dtau)
% relerr    ... Relative error for Runge-Kutta (default: 1e-8)

    p_Z_charge = libpointer('int32Ptr', Z_charge);
    p_A_mass = libpointer('int32Ptr', A_mass);
    p_E0 = libpointer('doublePtr', E0);
    p_dtau = libpointer('doublePtr', dtau);
    p_dtaumax = libpointer('doublePtr', dtaumax);
    p_relerr = libpointer('doublePtr', relerr);

    calllib('libneo_orb', 'neo_orb_global_MOD_init_params', p_Z_charge, ...
        p_A_mass, p_E0, p_dtau, p_dtaumax, p_relerr);
end
