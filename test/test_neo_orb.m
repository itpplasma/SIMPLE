%% Initialize field
addpath('/proj/plasma/CODE/NEO-ORB/matlab');
neo_orb_init('/proj/plasma/CODE/NEO-ORB');
neo_orb_init_field(5, 5, 3, -1);

%% Initialize parameters
Z_charge = 2; A_mass = 4;  E0 = 3.5e6;  % fast fusion alphas (deuterons)
dtau = 100;  % fixed timestep for evaluation, RK45 adapts inside
neo_orb_init_params(Z_charge, A_mass, E0, dtau, dtau, 1e-8)


%% Initialize conditions
ss=0.35;uu=0.33;vv=0.97;vpp=0.1;


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

%% Test orbit integration due to neo_orb.f90 in SRC 
syst=now;
nt = 1000;
z = zeros(nt, 4);
z(1,:) = [ss, uu, vv, vpp];

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
neosyst=now;
disp(['Neo Integration time',num2str((neosyst-syst)*24*3600)]);
%% Johanna Orbit Integration
%%%include fastBspline for all functions
addpath('/proj/plasma/CODE/NEO-ORB/test/fastBSpline');
%CompileMexFiles() %has to be done once; works for linux
%information for windows: https://es.mathworks.com/matlabcentral/fileexchange/32509-fast-b-spline-class?focused=e49de530-d908-d96f-5bf5-4e178065c9e5&tab=function

%%%loads vmec file in the subdirectory playvmec 
%(can be commented in second run -already loaded)
javaaddpath('/proj/plasma/CODE/NEO-ORB/test/netcdfAll-4.6.13.jar');
data=read_vmec('/proj/plasma/CODE/NEO-ORB/test/wout.nc');

%%%integration parameters
RelTol=1e-8; AbsTol=1e-12;

%%%evaluation functions of VMEC quantities
J=@(s,u,v) ceval(s,u,v,data.gmnc,data.xm,data.xn,0,0);
B=@(s,u,v,derivindex) ceval(s,u,v,data.bmnc,data.xm,data.xn,derivindex,0);
B_s=@(s,u,v,derivindex) -seval(s,u,v,data.bsubsmns,data.xm,data.xn,derivindex,1); 
B_u=@(s,u,v,derivindex) -ceval(s,u,v,data.bsubumnc,data.xm,data.xn,derivindex,0);
B_v=@(s,u,v,derivindex) -ceval(s,u,v,data.bsubvmnc,data.xm,data.xn,derivindex,0);
Bu=@(s,u,v,derivindex) -ceval(s,u,v,data.bsupumnc,data.xm,data.xn,derivindex,0);
Bv=@(s,u,v,derivindex) -ceval(s,u,v,data.bsupvmnc,data.xm,data.xn,derivindex,0);

%%%%physical constants
m=6.644657230*1e-8; %kg %normalised to e (*1e19)
q=2;
e=1.6021766208; %C %normalised(*1e19)

%%%evaluation of mu and v_||
%v2...thermal velocity squared
%rho_c...gyration radius with whole thermal velocity
%mu_c...normalised magnetic moment (/thermal velocity ^2 /m )
%u0_c...normalised initial conditions with v_paral/v as u0_c(1)
%t1...changed integration time t=t* thermal velocity
v2=2*T*e/m;
rho_c=m*sqrt(v2)/(q*e); 
mu_c=(1-vpp.^2)./(B(ss,uu,vv,0)*2);
u0_c=[vpp;ss;uu;vv];

%solving differential equation
%guidingcenter_c...littlejohn normalised system of equations
%integrationtime dtau*nt cm = dtau*nt/100 m
opts = odeset('RelTol',RelTol,'AbsTol',AbsTol);
[t,u] = ode45( @(t,u) guidingcenter_c(t,u,data,J,B,B_s,B_u,B_v,Bu,Bv,mu_c,rho_c),[0,dtau*nt/100],u0_c,opts);
disp(['VMEC Integration time',num2str((now-neosyst)*24*3600)]);
%% Johanna Comparison

tneo=[0:dtau:(nt-1)*dtau]./100;
figure
plot(tneo,z(:,1),t,u(:,2));
xlabel('t^n / m')
ylabel('s')
%saveas(gcf,['../../../../../user/ganglb_j/masterarbeit/latextemplate/figures/scompareoriginal',num2str(ceil(vpp*10)),'.eps'],'epsc');

figure
plot(tneo,z(:,2),t,u(:,3));
xlabel('t^n / m')
ylabel('u')
h1=legend('Fortran Routine','VMEC coefficient evaluation');
set(h1,'Location','north')
%saveas(gcf,['../../../../../user/ganglb_j/masterarbeit/latextemplate/figures/ucompareoriginal',num2str(ceil(vpp*10)),'.eps'],'epsc');

figure
plot(tneo,z(:,3),t,u(:,4));
xlabel('t^n / m')
ylabel('v')
%saveas(gcf,['../../../../../user/ganglb_j/masterarbeit/latextemplate/figures/vcompareoriginal',num2str(ceil(vpp*10)),'.eps'],'epsc');

figure
plot(tneo,z(:,4),t,u(:,1));
xlabel('t^n / m')
ylabel('{v_{||}}^n')
%saveas(gcf,['../../../../../user/ganglb_j/masterarbeit/latextemplate/figures/vpcompareoriginal',num2str(ceil(vpp*10)),'.eps'],'epsc');
%% Clean up
neo_orb_cleanup();