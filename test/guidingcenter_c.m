
%littlejohn normalised system of equations
function uprime=guidingcenter1(t,u,data,g,B,B_s,B_u,B_v,Bu,Bv,mu,rho_c)

%magnetic field quantities
mag_g=g(u(2),u(3),u(4));
mag_B=B(u(2),u(3),u(4),0);
nabla_B=[B(u(2),u(3),u(4),1);B(u(2),u(3),u(4),2);B(u(2),u(3),u(4),3)];
co_B=[B_s(u(2),u(3),u(4),0);B_u(u(2),u(3),u(4),0);B_v(u(2),u(3),u(4),0)];
contra_B=[0;Bu(u(2),u(3),u(4),0);Bv(u(2),u(3),u(4),0)];
curl_B=[B_v(u(2),u(3),u(4),2)-B_u(u(2),u(3),u(4),3);B_s(u(2),u(3),u(4),3)-B_v(u(2),u(3),u(4),1);B_u(u(2),u(3),u(4),1)-B_s(u(2),u(3),u(4),2)];

%Rotation of unit magnetic field B/abs_B
curl_b= -cross(nabla_B,co_B)./(mag_B^2*mag_g) +curl_B./(mag_B*mag_g);

%newly defined magnetic field 
b_star=(curl_b.*u(1)*rho_c+contra_B)./mag_B;

%parallelvalue of newly defined magnetic field 
b_star_b=dot(b_star,co_B)/mag_B;

%set up equation for concerning B-field
%%uprime(1)...parallel velocity U; uprime(2-4)...s,u,v- component of
%%driftvelocity vd in VMEC coordinates
uprime=zeros(4,1);
%uprime(1)=-((1-u(1).^2)./(2*b_star_b*mag_B))*dot(b_star,nabla_B);
uprime(1)=-(mu./b_star_b)*dot(b_star,nabla_B);
uprime(2:4)=(1./b_star_b)*(u(1)*b_star+mu*rho_c*cross(co_B,nabla_B)./(mag_B^2*mag_g));
end