function g=metric(s,u,v,Rcc,Zcc,lambda)


gV=zeros(3,3);
gV(1,1)= (Rcc(s,u,v,1)^2+Zcc(s,u,v,1)^2);
gV(1,2)=Rcc(s,u,v,1)*Rcc(s,u,v,2)+Zcc(s,u,v,1)*Zcc(s,u,v,2);
gV(1,3)=Rcc(s,u,v,1)*Rcc(s,u,v,3)+Zcc(s,u,v,1)*Zcc(s,u,v,3);
gV(2,1)=gV(1,2);
gV(2,2)=Rcc(s,u,v,2)^2+Zcc(s,u,v,2)^2;
gV(2,3)=Rcc(s,u,v,2)*Rcc(s,u,v,3)+Zcc(s,u,v,2)*Zcc(s,u,v,3);
gV(3,1)=gV(1,3);
gV(3,2)=gV(2,3);
gV(3,3)=Rcc(s,u,v,0)^2+Rcc(s,u,v,3)^2+Zcc(s,u,v,3)^2;

gmat=zeros(3,3);
cmat(1,2:3)=0.d0;
cmat(3,1:2)=0.d0;
cmat(1,1)=1.d0;
cmat(3,3)=1.d0;
cmat(2,1)=-lambda(s,u,v,1)/(1+lambda(s,u,v,2));
cmat(2,2)=1/(1+lambda(s,u,v,2));
cmat(2,3)=-lambda(s,u,v,3)/(1+lambda(s,u,v,2));

g=transpose(cmat)*(gV*cmat);
end