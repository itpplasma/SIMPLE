function val=iotaeval(s,data,derivindex)
    x=linspace(0,1,length(data));
    deltax=x(2)-x(1);
    k=(data(end)-data(1))/(x(end)-x(1)); d=data(1)-k*x(1);
    y = data-k*x-d*ones(size(x));
    knots = [x(1)-2*deltax,x(1)-deltax,x(1)-deltax,x,x(end)+deltax,x(end)+deltax,x(end)+2*deltax];
    switch derivindex
        case 0
            sp = fastBSpline.lsqspline(knots,3,x,y);
            val=sp.evalAt(s)+k*s+d*ones(size(s));
        case 1
            sp = fastBSpline.lsqspline(knots,3,x,y);
            sp=sp.dx;
            val=sp.evalAt(s)+k*ones(size(s));
        case 2
            sp = fastBSpline.lsqspline(knots,3,x,y);
            sp=sp.dx.dx;
            val=sp.evalAt(s);
    end
end