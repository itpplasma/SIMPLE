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
B_s=@(s,u,v,derivindex) seval(s,u,v,data.bsubsmns,data.xm,data.xn,derivindex,0);
B_u=@(s,u,v,derivindex) ceval(s,u,v,data.bsubumnc,data.xm,data.xn,derivindex,0);
B_v=@(s,u,v,derivindex) ceval(s,u,v,data.bsubvmnc,data.xm,data.xn,derivindex,0);
Bu=@(s,u,v,derivindex) ceval(s,u,v,data.bsupumnc,data.xm,data.xn,derivindex,0);
Bv=@(s,u,v,derivindex) ceval(s,u,v,data.bsupvmnc,data.xm,data.xn,derivindex,0);
Rcc=@(s,u,v,derivindex)  ceval(s,u,v,data.rmnc,data.xm,data.xn,derivindex,1)*1e2;
Zcc=@(s,u,v,derivindex)  seval(s,u,v,data.zmns,data.xm,data.xn,derivindex,1)*1e2;
lambda=@(s,u,v,derivindex) seval(s,u,v,data.lmns,data.xm,data.xn,derivindex,0);


%% evaluate different fields

x=linspace(0.2,0.8,64);
u=2;v=0.1;
hcompare=zeros(3,length(x));hhcompare=hcompare;
sqrtg=zeros(1,length(x)); bmod=sqrtg;
Gi=zeros(size(x)); Bi=zeros(size(x));
Bvi=zeros(size(x));Bui=zeros(size(x));B_si=zeros(size(x));B_vi=zeros(size(x));B_ui=zeros(size(x));
bvii=zeros(size(x));buii=zeros(size(x));bvi=zeros(size(x));bui=zeros(size(x));
b_sii=zeros(size(x));b_uii=zeros(size(x));b_vii=zeros(size(x));
b_si=zeros(size(x));b_ui=zeros(size(x));b_vi=zeros(size(x));
bi=zeros(size(x));

for i=1:length(x)
    [bmod(i), sqrtg(i), ~, hcovar,hcontrvar,~]=neo_orb_magfie_vmec([(x(i)), u, v]);
    hcompare(:,i)=hcovar(:,1);
    hhcompare(:,i)=hcontrvar(:,1);
    ggg=metric(x(i),u,v,Rcc,Zcc,lambda);
    
    Gi(i)=J(x(i),u,v); 
    Bi(i)=B(x(i),u,v,0);
    Bvi(i)=Bv((x(i)),u,v,0); 
    Bui(i)=Bu((x(i)),u,v,0);
    B_si(i)=B_s((x(i)),u,v,0);
    B_ui(i)=B_u((x(i)),u,v,0);
    B_vi(i)=B_v((x(i)),u,v,0);
    
    bvii(i)=1e8*(1+lambda(x(i),u,v,2))/(sqrtg(i));
    buii(i)=bvii(i)*iotaeval(x(i),data.iotaf,0);
    bui(i)=(buii(i)-lambda(x(i),u,v,3)*bvii(i))/(1+lambda(x(i),u,v,2));
    bvi(i)=bvii(i);
    b_sii(i)=ggg(1,2)*buii(i)+ggg(1,3)*bvii(i);
    b_uii(i)=ggg(2,2)*buii(i)+ggg(2,3)*bvii(i);
    b_vii(i)=ggg(3,2)*buii(i)+ggg(3,3)*bvii(i);
    b_si(i)=b_sii(i)+lambda(x(i),u,v,1)*b_uii(i);
    b_ui(i)=(1+lambda(x(i),u,v,2))*b_uii(i);
    b_vi(i)=b_vii(i)+lambda(x(i),u,v,3)*b_uii(i);
    bi(i)=sqrt(b_uii(i)*buii(i)+b_vii(i)*bvii(i));
    
end

%plot jacobian
figure 
subplot(2,1,1)
plot(x,sqrtg)
xlabel('s/ fluxsurface label')
ylabel('g/cm^3')
title('Fortran Routine')
subplot(2,1,2)
plot(x,Gi)
xlabel('s/ fluxsurface label')
ylabel('g/ m^3')
title('VMEC coefficient evaluation')

%plot magnetic field magnitude
figure 
subplot(2,1,1)
plot(x,bi)
xlabel('s/ fluxsurface label')
ylabel('|B|/G')
title('Fortran Routine')
subplot(2,1,2)
plot(x,Bi)
xlabel('s/ fluxsurface label')
ylabel('|B|/T')
title('VMEC coefficient evaluation')

%plot contravariant vector
h=gobjects(2,1);
figure 
h(1)=subplot(2,1,1);
h1=plot(x,bui,x,bvi);
%hold on
%plot(x,hhcompare(3,:).*bmod);
ylabel('contravariant unit vector')
title('Fortran Routine')
h(2)=subplot(2,1,2);
plot(x,Bui,x,Bvi);
xlabel('s/ fluxsurface label')
title('VMEC coefficient evaluation')
lgd=legend([h1],'u','v');
set(h(1),'Position',[0.13,0.65,0.7750,0.3]);
set(h(2),'Position',[0.13,0.25,0.7750,0.3]);
set(lgd,'Position',[0.13,0.05,0.7750,0.05]);

%plot covariant vector
figure 
h(1)=subplot(2,1,1);
h1=plot(x,b_si,x,b_ui,x,b_vi);
title('Fortran Routine')
h(2)=subplot(2,1,2);
plot(x,B_si,x,B_ui,x,B_vi);
xlabel('s/ fluxsurface label')
ylabel('covariant unit vector')
title('VMEC coefficient evaluation')
lgd=legend([h1],'s','u','v');
set(h(1),'Position',[0.13,0.65,0.7750,0.3]);
set(h(2),'Position',[0.13,0.25,0.7750,0.3]);
set(lgd,'Position',[0.13,0.05,0.7750,0.05]);
