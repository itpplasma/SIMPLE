%% Initialize field
parentpath='/afs/itp.tugraz.at/proj/plasma/CODE/NEO-ORB';
addpath([parentpath,'/matlab']);
neo_orb_init(parentpath);
neo_orb_init_field(5, 5, 3, -1);

%% Initialize parameters
Z_charge = 2; A_mass = 4;  E0 = 3.5e6;  % fast fusion alphas (deuterons)
dtau = 100;  % fixed timestep for evaluation, RK45 adapts inside
neo_orb_init_params(Z_charge, A_mass, E0, dtau, dtau, 1e-8)

%% include fastBspline
%enter path of directory in variable parentpath

addpath([parentpath,'/test/fastBSpline']);
%CompileMexFiles() %has to be done once; works for linux
%information for windows: https://es.mathworks.com/matlabcentral/fileexchange/32509-fast-b-spline-class?focused=e49de530-d908-d96f-5bf5-4e178065c9e5&tab=function

%% loads vmec file
%(can be commented in second run -already loaded)
javaaddpath([parentpath, '/test/netcdfAll-4.6.13.jar']);
data=read_vmec([parentpath, '/test/wout.nc']);


%% evaluation functions of VMEC quantities
J=@(s,u,v) ceval(s,u,v,data.gmnc,data.xm,data.xn,0,0);
B=@(s,u,v,derivindex) ceval(s,u,v,data.bmnc,data.xm,data.xn,derivindex,0);
B_s=@(s,u,v,derivindex) seval(s,u,v,data.bsubsmns,data.xm,data.xn,derivindex,1);
B_u=@(s,u,v,derivindex) ceval(s,u,v,data.bsubumnc,data.xm,data.xn,derivindex,0);
B_v=@(s,u,v,derivindex) ceval(s,u,v,data.bsubvmnc,data.xm,data.xn,derivindex,0);
Bu=@(s,u,v,derivindex) ceval(s,u,v,data.bsupumnc,data.xm,data.xn,derivindex,0);
Bv=@(s,u,v,derivindex) ceval(s,u,v,data.bsupvmnc,data.xm,data.xn,derivindex,0);
Rcc=@(s,u,v,derivindex)  ceval(s,u,v,data.rmnc,data.xm,data.xn,derivindex,1);
Zcc=@(s,u,v,derivindex)  seval(s,u,v,data.zmns,data.xm,data.xn,derivindex,1);
lambda=@(s,u,v,derivindex) seval(s,u,v,data.lmns,data.xm,data.xn,derivindex,0);

%% evaluate field on s grid

x=linspace(0.1,1,99);
u=2;v=0.1;
hcompare=zeros(3,99);hhcompare=hcompare;
sqrtg=zeros(1,99); bmod=sqrtg;
gi=zeros(size(x)); bi=zeros(size(x));
Bvi=zeros(size(x));Bui=zeros(size(x));B_si=zeros(size(x));B_vi=zeros(size(x));B_ui=zeros(size(x));

for i=1:length(x)
    [bmod(i), sqrtg(i), ~, hcovar,hcontrvar,~]=neo_orb_magfie_vmec([(x(i)), u, v]);
    hcompare(:,i)=hcovar(:,1);
    hhcompare(:,i)=hcontrvar(:,1);

    gi(i)=J(x(i),u,v); %gii(i)=Rcc(x(i),u,v,0)*(Rcc(x(i),u,v,1)*Zcc(x(i),u,v,2)-Rcc(x(i),u,v,2)*Zcc(x(i),u,v,1));
    bi(i)=B(x(i),u,v,0);
    Bvi(i)=Bv((x(i)),u,v,0); %Bvii(i)=data.phi(end)*(1+lambda(x(i),u,v,2))/J(x(i),u,v);
    Bui(i)=Bu((x(i)),u,v,0);
    B_si(i)=B_s((x(i)),u,v,0)/bi(i);
    B_ui(i)=B_u((x(i)),u,v,0)/bi(i);
    B_vi(i)=B_v((x(i)),u,v,0)/bi(i);
end

facb=mean(bmod*1e-4/bi);
facbu=mean(hhcompare(2,:)*1e4-Bui);
facbv=mean(hhcompare(3,:)*1e4-Bvi);


figure 
subplot(2,1,1)
plot(x,sqrtg)
xlabel('s/ fluxsurface label')
ylabel('g/cm^3')
title('Fortran Routine')
subplot(2,1,2)
plot(x,gi)
xlabel('s/ fluxsurface label')
ylabel('g/ m^3')
title('VMEC coefficient evaluation')

figure 
subplot(2,1,1)
plot(x,bmod)
xlabel('s/ fluxsurface label')
ylabel('|B|/G')
title('Fortran Routine')
subplot(2,1,2)
plot(x,bi)
xlabel('s/ fluxsurface label')
ylabel('|B|/T')
title('VMEC coefficient evaluation')

figure 
subplot(1,2,1)
plot(x,hhcompare(2,:),x,hhcompare(3,:));
xlabel('s/ fluxsurface label')
ylabel('contravariant unit vector')
title('Fortran Routine')
subplot(1,2,2)
plot(x,Bui,x,Bvi);
xlabel('s/ fluxsurface label')
ylabel('contravariant unit vector')
title('VMEC coefficient evaluation')
legend('u','v')

figure 
subplot(1,2,1)
plot(x,hcompare(1,:),x,hcompare(2,:),x,hcompare(3,:));
xlabel('s/ fluxsurface label')
ylabel('covariant unit vector')
title('Fortran Routine')
subplot(1,2,2)
plot(x,B_si,x,B_ui,x,B_vi);
xlabel('s/ fluxsurface label')
xlabel('s/ fluxsurface label')
ylabel('covariant unit vector')
title('VMEC coefficient evaluation')
legend('s','u','v')