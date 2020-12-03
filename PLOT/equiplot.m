%%
javaaddpath('/home/calbert/src/matlabVMECv3.00/netcdfAll-4.2.jar');

% %%
% data = read_vmec('../data/QI/wout_23_1900_fix_bdry.nc');
% VMECplot(data);
% 
% %%
% data = read_vmec('../data/QH/wout_qh_8_7.nc');
% VMECplot(data);
% 
% %%
%data = read_vmec('../data/QA/wout_ants1.nc');
%VMECplot(data);

%data = read_vmec('../data/QA/wout_publishedQA2019.nc');
%VMECplot(data);
%data = read_vmec('wout.nc');

%data = read_vmec('/home/calbert/net/cobra/run/VMEC/QA/12x12/wout_qa_12x12.nc');
data = read_vmec('/home/calbert/net/cobra/run/VMEC/QA/12x12_128/wout_qa_12x12.nc');


%%

%%
fsind = 64;  % s = 0.25

theta=0:2*pi/(359):2*pi;    % 36 points define theta [0,6.28]
phi=0:2*pi/(359):2*pi;     % 361 points define phi [0,6.28]
r=cfunct(theta,phi,data.rmnc,data.xm,data.xn);  % Cos transform
z=sfunct(theta,phi,data.zmns,data.xm,data.xn);  % Sin transform

%% see https://de.mathworks.com/matlabcentral/answers/
%      127342-access-a-single-element-of-an-anonymous-function-that-returns-an-array
indexat = @(expr, index) expr(index);
xdot = @(t, x) [
    indexat(cfunct(x(1),x(2),data.bsupumnc,data.xmnyq,-data.xnnyq),fsind);
    indexat(cfunct(x(1),x(2),data.bsupvmnc,data.xmnyq,-data.xnnyq),fsind)    
    ];

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t,x] = ode45(xdot, linspace(0,35.534,1000), [pi; 0], options);

nx = size(x,1);
bfield = zeros(nx,1); 
for k = 1:nx
    fie = cfunct(x(k,1),x(k,2),data.bmnc,data.xmnyq,-data.xnnyq);
    bfield(k) = fie(fsind);
end

refdata = dlmread('../data/bfield_s025_sophia.csv');
tnorm = t/max(t)*5;

figure()
plot(refdata(:,1), refdata(:,2))
hold on
plot(tnorm, bfield/4.9)
hold off
ylim([0.905, 1.1])
legend('Henneberg2019', 'wout publishedQA2019.nc')
legend('Location','north')

%%

theta=0:2*pi/(143):2*pi;    % 36 points define theta [0,6.28]
phi=0:2*pi/(143):pi;     % 361 points define phi [0,6.28]
r=cfunct(theta,phi,data.rmnc,data.xm,data.xn);  % Cos transform
z=sfunct(theta,phi,data.zmns,data.xm,data.xn);  % Sin transform

bu=cfunct(theta,phi,data.bsupumnc,data.xmnyq,-data.xnnyq);
bv=cfunct(theta,phi,data.bsupvmnc,data.xmnyq,-data.xnnyq);
b=cfunct(theta,phi,data.bmnc,data.xmnyq,-data.xnnyq);

drdu=sfunct(theta,phi,data.rumns,data.xm,data.xn);
drdv=sfunct(theta,phi,data.rvmns,data.xm,data.xn);
dzdu=cfunct(theta,phi,data.zumnc,data.xm,data.xn);
dzdv=cfunct(theta,phi,data.zvmnc,data.xm,data.xn);
brho=bu.*drdu+bv.*drdv;    % Make B_R
bphi=bv.*r;                % Make B_phi
bz=bu.*dzdu+bv.*dzdv;      % Make B_Z

figure()
refcont = dlmread('../data/bfield_s025_sophia_contour.csv');
refxy = dlmread('../data/bfield_s025_sophia_xy.csv');
[PH, TH] = meshgrid(phi, theta);
contour(PH*180/pi, TH*180/pi, squeeze(b(fsind,:,:)), 18)
hold on
plot(refxy(:,1), refxy(:,2))
%plot(x(:,2)*180/pi, x(:,1)*180/pi)
plot(refcont(:,1), refcont(:,2),'r.')
hold off
xlim([0, 180])

colormap('jet')
