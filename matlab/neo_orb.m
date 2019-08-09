addpath('/home/calbert/build/NEO-ORB')

if not(libisloaded('libneo_orb'))
    loadlibrary('libneo_orb', 'neo_orb.h')
end

libfunctions('libneo_orb', '-full')

ns = libpointer('int32Ptr', 5);
ntp = libpointer('int32Ptr', 5);
multharm = libpointer('int32Ptr', 3);
calllib('libneo_orb', 'neo_orb_global_MOD_init_field', ns, ntp, multharm);

Z_charge = libpointer('int32Ptr', 2);
m_mass = libpointer('int32Ptr', 4);
E_kin = libpointer('doublePtr', 3.5e6);
dtau = libpointer('doublePtr', 100);
dtaumax = libpointer('doublePtr', 100);
integmode = libpointer('int32Ptr', -1);
relerr = libpointer('doublePtr', 1e-8);

calllib('libneo_orb', 'neo_orb_global_MOD_init_params', Z_charge, m_mass, ...
    E_kin, dtau, dtaumax, integmode, relerr);
%%
s = libpointer('doublePtr', 0.5);
th = libpointer('doublePtr', 0);
ph = libpointer('doublePtr', 0);
lam = libpointer('doublePtr', 0.1);
ierr = libpointer('int32Ptr', 0);

nt = 10000;
z = zeros(nt, 4);

for k = 1:nt
  calllib('libneo_orb', 'neo_orb_global_MOD_timestep', s, th, ph, lam, ierr);
  z(k, 1) = s.Value;
  z(k, 2) = th.Value;
  z(k, 3) = ph.Value;
  z(k, 4) = lam.Value;
end

plot(z(:,1).*cos(z(:,2)), z(:,1).*sin(z(:,2)))

%%
unloadlibrary('libneo_orb')
