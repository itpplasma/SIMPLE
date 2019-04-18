%% Initialize field
addpath('/proj/plasma/CODE/NEO-ORB/matlab');
neo_orb_init('/proj/plasma/CODE/NEO-ORB');
neo_orb_init_field(5, 5, 3, -1);

%% Initialize parameters
Z_charge = 2; A_mass = 4;  E0 = 3.5e6;  % fast fusion alphas (deuterons)
dtau = 100;  % fixed timestep for evaluation, RK45 adapts inside
neo_orb_init_params(Z_charge, A_mass, E0, dtau, dtau, 1e-8)


%% Test computing of field
% Create grid in poloidal plane in VMEC coordinates at fixed phv=ph0
[s, thv] = meshgrid(linspace(0.1, 1.0, 30), linspace(0.0, 2*pi, 30));
ph0 = 0.0;

bmodfield = zeros(size(s));

% Evaluate field on grid
for k = 1:numel(s)
  [bmod, sqrtg, bder, hcovar, hctrvar, hcurl] = neo_orb_magfie_vmec(...
                                                  [s(k), thv(k), ph0]);
  bmodfield(k) = bmod;
end

figure()  % Plot in "pseudo"-cylindrical coordinates in poloidal plane
contourf(s.*cos(thv), s.*sin(thv), bmodfield)
hold on

%% Test orbit integration
nt = 1000;
z = zeros(nt, 4);
z(1,:) = [0.35, 0.333, 0.97, 0.1];

for k = 1:nt-1
  z(k+1,:) = neo_orb_timestep(z(k,:));
end
plot(z(:,1).*cos(z(:,2)), z(:,1).*sin(z(:,2)), 'k')
hold off

figure() % Plot orbit in "pseudo"-Cartesian coordinates in 3D
R0 = 2.0;
plot3((R0+z(:,1).*cos(z(:,2))).*cos(z(:,3)), ...
      (R0+z(:,1).*cos(z(:,2))).*sin(z(:,3)), ...
      z(:,1).*sin(z(:,2)), 'k')
xlim([-3.0, 3.0])
ylim([-3.0, 3.0])
zlim([-1.0, 1.0])

%% Johanna comparisons
u=2;v=0.1;
hcompare=zeros(3,99);hhcompare=hcompare;s=linspace(0,1,99);
sqrtg=zeros(1,99); bmod=sqrtg;
for i=1:99
    [bmod(i), sqrtg(i), ~, hcovar,hcontrvar,~]=neo_orb_magfie_vmec([(s(i)), u, v]);
    hcompare(:,i)=hcovar(:,1);
    hhcompare(:,i)=hcontrvar(:,1)*bmod(i);
end

figure
plot(s,sqrtg)
xlabel('s/ fluxsurface label')
ylabel('g/')

figure
plot(s,bmod)
xlabel('s/ fluxsurface label')
ylabel('b/')

figure
plot(s,hhcompare(1,:),s,hhcompare(2,:),s,hhcompare(3,:));
xlabel('s/ fluxsurface label')

figure
plot(s,hcompare(1,:),s,hcompare(2,:),s,hcompare(3,:));
xlabel('s/ fluxsurface label')
%% Clean up
neo_orb_cleanup();