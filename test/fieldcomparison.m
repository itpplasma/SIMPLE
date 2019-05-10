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
B=@(s,u,v,derivindex) 2*pi*ceval(s,u,v,data.bmnc,data.xm,data.xn,derivindex,0);
B_s=@(s,u,v,derivindex) -2*pi*seval(s,u,v,data.bsubsmns,data.xm,data.xn,derivindex,0);
B_u=@(s,u,v,derivindex) -2*pi*ceval(s,u,v,data.bsubumnc,data.xm,data.xn,derivindex,0);
B_v=@(s,u,v,derivindex) -2*pi*ceval(s,u,v,data.bsubvmnc,data.xm,data.xn,derivindex,0);
Bu=@(s,u,v,derivindex) -2*pi*ceval(s,u,v,data.bsupumnc,data.xm,data.xn,derivindex,0);
Bv=@(s,u,v,derivindex) -2*pi*ceval(s,u,v,data.bsupvmnc,data.xm,data.xn,derivindex,0);
Rcc=@(s,u,v,derivindex)  2*pi*ceval(s,u,v,data.rmnc,data.xm,data.xn,derivindex,1)*1e2;
Zcc=@(s,u,v,derivindex)  -2*pi*seval(s,u,v,data.zmns,data.xm,data.xn,derivindex,1)*1e2;
lambda=@(s,u,v,derivindex) -2*pi*seval(s,u,v,data.lmns,data.xm,data.xn,derivindex,0);


%% evaluate different fields

x=linspace(0.2,0.8,64);
u=0.3333;v=0.97;
hcompare=zeros(3,length(x));hhcompare=hcompare;hbder=hcompare;hhcurl=hcompare;
sqrtg=zeros(1,length(x)); bmod=sqrtg;

Gi=zeros(size(x)); Bi=zeros(size(x));
Bbi=zeros(3,length(x));B_i=zeros(3,length(x));
Bderi=zeros(3,length(x));Curlbi=Bderi;

bvii=zeros(size(x));buii=zeros(size(x));bvi=zeros(size(x));bui=zeros(size(x));
b_sii=zeros(size(x));b_uii=zeros(size(x));b_vii=zeros(size(x));
b_si=zeros(size(x));b_ui=zeros(size(x));b_vi=zeros(size(x));
bi=zeros(size(x)); 


for i=1:length(x)
    [bmod(i), sqrtg(i), bder, hcovar,hcontrvar,hcurl]=neo_orb_magfie_vmec([(x(i)), u, v]);
    hcompare(:,i)=hcovar(:,1);
    hhcompare(:,i)=hcontrvar(:,1);
    hbder(:,i)=bder(:,1);
    hhcurl(:,i)=hcurl(:,1);
    ggg=metric(x(i),u,v,Rcc,Zcc,lambda);
    
    Gi(i)=J(x(i),u,v); 
    Bi(i)=B(x(i),u,v,0);
    Bbi(:,i)=[0;Bu((x(i)),u,v,0);Bv((x(i)),u,v,0)];
    B_i(:,i)=[B_s((x(i)),u,v,0);B_u((x(i)),u,v,0);B_v((x(i)),u,v,0)];
    Bderi(:,i)=[B(x(i),u,v,1);B(x(i),u,v,2);B(x(i),u,v,3)];
    Curlbi(:,i)= -cross(Bderi(:,i),B_i(:,i))./(Bi(i)^2*Gi(i))+...
            [B_v(x(i),u,v,2)-B_u(x(i),u,v,3);B_s(x(i),u,v,3)-B_v(x(i),u,v,1);B_u(x(i),u,v,1)-B_s(x(i),u,v,2)]./(Bi(i)*Gi(i));
   
    
    bvii(i)=data.phi(end)*1e8*(1+lambda(x(i),u,v,2))/(sqrtg(i));
    buii(i)=bvii(i)*iotaeval(x(i),data.iotaf,0);
    bui(i)=(buii(i)-lambda(x(i),u,v,3)*bvii(i))/(1+lambda(x(i),u,v,2));
    bvi(i)=bvii(i);
    b_sii(i)=ggg(1,2)*buii(i)+ggg(1,3)*bvii(i);
    b_uii(i)=ggg(2,2)*buii(i)+ggg(2,3)*bvii(i);
    b_vii(i)=ggg(3,2)*buii(i)+ggg(3,3)*bvii(i);
    b_si(i)=b_sii(i)+lambda(x(i),u,v,1)*b_uii(i);
    b_ui(i)=(1+lambda(x(i),u,v,2))*b_uii(i);
    b_vi(i)=b_vii(i)+lambda(x(i),u,v,3)*b_uii(i);
    bi(i)=sqrt(b_ui(i)*bui(i)+b_vi(i)*bvi(i));
    
    
end

hbmodu=zeros(1,1000); hbmodv=hbmodu; Bmodu=hbmodu; Bmodv=hbmodv;
uu=linspace(0,2*pi,1000);s=0.35;vv=linspace(0,2*pi/5,1000);
for i=1:length(uu)
    Bmodu(i)=B(s,uu(i),v,0);
    Bmodv(i)=B(s,u,vv(i),0);
    [hbmodu(i), ~, ~, ~,~,~]=neo_orb_magfie_vmec([s, uu(i), v]);
    [hbmodv(i), ~, ~, ~,~,~]=neo_orb_magfie_vmec([s, u, vv(i)]);
end
    

%plot jacobian
figure 
subplot(2,1,1)
plot(x,sqrtg)
xlabel('s / fluxsurface label')
ylabel('g / cm^3')
title('Fortran Routine')
subplot(2,1,2)
plot(x,Gi)
xlabel('s / fluxsurface label')
ylabel('g / m^3')
title('VMEC coefficient evaluation')

%plot magnetic field magnitude
figure 
subplot(2,1,1)
plot(x,bmod)
xlabel('s / fluxsurface label')
ylabel('|B| / G')
title('Fortran Routine')
subplot(2,1,2)
plot(x,Bi)
xlabel('s / fluxsurface label')
ylabel('|B| / T')
title('VMEC coefficient evaluation')

%plot contravariant vector
h=gobjects(2,1);
figure 
h(1)=subplot(2,1,1);
h1=plot(x,hhcompare(2,:),x,hhcompare(3,:));
ylabel('contravariant / cm^{-1}')
title('Fortran Routine')
h(2)=subplot(2,1,2);
plot(x,Bbi(2,:)./Bi,x,Bbi(3,:)./Bi);
xlabel('s / fluxsurface label')
ylabel('contravariant / m^{-1}')
title('VMEC coefficient evaluation')
lgd=legend([h1],'u','v');
set(h(1),'Position',[0.13,0.65,0.7750,0.3]);
set(h(2),'Position',[0.13,0.25,0.7750,0.3]);
set(lgd,'Position',[0.13,0.05,0.7750,0.05]);

%plot covariant vector
figure 
h(1)=subplot(2,1,1);
h1=plot(x,hcompare(1,:),x,hcompare(2,:),x,hcompare(3,:));
ylabel('covariant / cm')
title('Fortran Routine')
h(2)=subplot(2,1,2);
plot(x,B_i(1,:)./Bi,x,B_i(2,:)./Bi,x,B_i(3,:)./Bi);
xlabel('s/ fluxsurface label')
ylabel('covariant / m')
title('VMEC coefficient evaluation')
lgd=legend([h1],'s','u','v');
set(h(1),'Position',[0.13,0.65,0.7750,0.3]);
set(h(2),'Position',[0.13,0.25,0.7750,0.3]);
set(lgd,'Position',[0.13,0.05,0.7750,0.05]);

%plot derivative of B
h=gobjects(2,1);
figure 
h(1)=subplot(2,1,1);
h1=plot(x,hbder(1,:),x,hbder(2,:),x,hbder(3,:));
ylabel('\nabla B / |B|')
title('Fortran Routine')
h(2)=subplot(2,1,2);
plot(x,Bderi(1,:)./Bi,x,Bderi(2,:)./Bi,x,Bderi(3,:)./Bi);
xlabel('s / fluxsurface label')
ylabel('\nabla B / |B|')
title('VMEC coefficient evaluation')
lgd=legend([h1],'s','u','v');
set(h(1),'Position',[0.13,0.65,0.7750,0.3]);
set(h(2),'Position',[0.13,0.25,0.7750,0.3]);
set(lgd,'Position',[0.13,0.05,0.7750,0.05]);

%plot curl of b
h=gobjects(2,1);
figure 
h(1)=subplot(2,1,1);
h1=plot(x,hhcurl(1,:),x,hhcurl(2,:),x,hhcurl(3,:));
ylabel('\nabla x b / cm^-2')
title('Fortran Routine')
h(2)=subplot(2,1,2);
plot(x,Curlbi(1,:),x,Curlbi(2,:),x,Curlbi(3,:));
xlabel('s / fluxsurface label')
ylabel('\nabla x b / m^-2')
title('VMEC coefficient evaluation')
lgd=legend([h1],'s','u','v');
set(h(1),'Position',[0.13,0.65,0.7750,0.3]);
set(h(2),'Position',[0.13,0.25,0.7750,0.3]);
set(lgd,'Position',[0.13,0.05,0.7750,0.05]);


