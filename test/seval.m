%%%%%get values by setting togeter the given fourierseries
%%%%%f(s,u,v)=sum(R_mn (s) * sin(m*theta-n*phi))
%s...flux surface staying on
%u...vmec poloidal angle
%v... vmec toroidal angle
%coef...Fouriercoefficients of VMEC
%xm...slow running m for fast running n (1-11)
%xn...fast running n for slow running n (-60 - 60)
%derivindex....0 no derivative, 1 derivative in s, 2 derivative in u, 3 derivative in v
%mesh =0 for half mesh; =1 for full mesh

function val=seval(s,u,v,coef,xm,xn,derivindex,mesh)
    sdim=size(coef,2); 
    x=linspace(0,1,sdim);
    deltax=x(2)-x(1);
    if mesh==0
        sdim=sdim-1;
        x=linspace(0+deltax*0.5,1-deltax*0.5,sdim);
    end
    knots = [x(1)-2*deltax,x(1)-2*deltax,x(1)-deltax,x,x(end)+deltax,x(end)+2*deltax,x(end)+2*deltax];
    
    switch(derivindex)
        case 0
            valvec=sin(xm.*u+v.*xn)*coef; 
            if mesh==0, valvec=valvec(2:end);end
            k=(valvec(end)-valvec(1))/(x(end)-x(1)); d=valvec(1)-k*x(1);
            y = valvec-k*x-d*ones(size(x));
            sp = fastBSpline.lsqspline(knots,3,x,y);
            val=sp.evalAt(s)+k*s+d*ones(size(s));
        case 1
            valvec=sin(xm.*u+v.*xn)*coef; 
            if mesh==0, valvec=valvec(2:end);end
            k=(valvec(end)-valvec(1))/(x(end)-x(1)); d=valvec(1)-k*x(1);
            y = valvec-k*x-d*ones(size(x));
            sp = fastBSpline.lsqspline(knots,3,x,y);
            sp=sp.dx;
            val=sp.evalAt(s)+k*ones(size(s));
        case 2
            valvec=(cos(xm.*u+v.*xn).*xm)*coef; 
            if mesh==0, valvec=valvec(2:end);end
            k=(valvec(end)-valvec(1))/(x(end)-x(1)); d=valvec(1)-k*x(1);
            y = valvec-k*x-d*ones(size(x));
            sp = fastBSpline.lsqspline(knots,3,x,y);
            val=sp.evalAt(s)+k*s+d*ones(size(s));
        case 3
            valvec=(cos(xm.*u+v.*xn).*xn)*coef; 
            if mesh==0, valvec=valvec(2:end);end
            k=(valvec(end)-valvec(1))/(x(end)-x(1)); d=valvec(1)-k*x(1);
            y = valvec-k*x-d*ones(size(x));
            sp = fastBSpline.lsqspline(knots,3,x,y);
            val=sp.evalAt(s)+k*s+d*ones(size(s));
        case 11
            valvec=sin(xm.*u+v.*xn)*coef; 
            if mesh==0, valvec=valvec(2:end);end
            k=(valvec(end)-valvec(1))/(x(end)-x(1)); d=valvec(1)-k*x(1);
            y = valvec-k*x-d*ones(size(x));
            sp = fastBSpline.lsqspline(knots,3,x,y);
            sp=sp.dx.dx;
            val=sp.evalAt(s);
        case {12,21}
            valvec=(cos(xm.*u+v.*xn).*xm)*coef; 
            if mesh==0, valvec=valvec(2:end);end
            k=(valvec(end)-valvec(1))/(x(end)-x(1)); d=valvec(1)-k*x(1);
            y = valvec-k*x-d*ones(size(x));
            sp = fastBSpline.lsqspline(knots,3,x,y);
            sp=sp.dx;
            val=sp.evalAt(s)+k*ones(size(s));
        case {13,31}
            valvec=(cos(xm.*u+v.*xn).*xn)*coef; 
            if mesh==0, valvec=valvec(2:end);end
            k=(valvec(end)-valvec(1))/(x(end)-x(1)); d=valvec(1)-k*x(1);
            y = valvec-k*x-d*ones(size(x));
            sp = fastBSpline.lsqspline(knots,3,x,y);
            sp=sp.dx;
            val=sp.evalAt(s)+k*ones(size(s));
        case 22
            valvec=-(sin(xm.*u+v.*xn).*(xm.^2))*coef; 
            if mesh==0, valvec=valvec(2:end);end
            k=(valvec(end)-valvec(1))/(x(end)-x(1)); d=valvec(1)-k*x(1);
            y = valvec-k*x-d*ones(size(x));
            sp = fastBSpline.lsqspline(knots,3,x,y);
            val=sp.evalAt(s)+k*s+d*ones(size(s));
        case {23,32}
            valvec=-(sin(xm.*u+v.*xn).*(xm.*xn))*coef; 
            if mesh==0, valvec=valvec(2:end);end
            k=(valvec(end)-valvec(1))/(x(end)-x(1)); d=valvec(1)-k*x(1);
            y = valvec-k*x-d*ones(size(x));
            sp = fastBSpline.lsqspline(knots,3,x,y);
            val=sp.evalAt(s)+k*s+d*ones(size(s));
        case 33
            valvec=-(sin(xm.*u+v.*xn).*(xn.^2))*coef; 
            if mesh==0, valvec=valvec(2:end);end
            k=(valvec(end)-valvec(1))/(x(end)-x(1)); d=valvec(1)-k*x(1);
            y = valvec-k*x-d*ones(size(x));
            sp = fastBSpline.lsqspline(knots,3,x,y);
            val=sp.evalAt(s)+k*s+d*ones(size(s));
     end
end