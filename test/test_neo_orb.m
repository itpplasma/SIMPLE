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
figure
subplot(2,2,1)
plot([1:1:length(z)],z(:,1));
title('s');
subplot(2,2,2)
plot([1:1:length(z)],z(:,2));
title('u');
subplot(2,2,3)
plot([1:1:length(z)],z(:,3));
title('v');
subplot(2,2,4)
plot([1:1:length(z)],z(:,4));
title('v_{||}/ v_{thermal}');

bmu=zeros(length(z),1);
for i=1:length(z)
    [bmu(i), ~, ~, ~, ~, ~] = neo_orb_magfie_vmec([z(i,1), z(i,2), z(i,3)]);
end

figure
plot([1:1:length(z)],(ones(size(bmu))-z(:,4).^2)./bmu)
%% Clean up
neo_orb_cleanup();