#plot "euler1_quasi.out"   u 5 ps .001 
#replot "euler2_quasi.out" u 5 ps .001

#plot "gauss2_quasi.out" u 5 ps .001  
#replot "midpoint_quasi.out" u 5 ps .001  
#replot   "verlet_quasi.out" u 5 ps .001

#plot "kahan8q.out" u 5 ps .001
#replot "gauss8q.out"   u 5 ps .001

plot "gauss2.out" u 5 ps .001
#replot "verlet.out" u 5 ps .001
replot "lobatto4q.out"   u 5 ps .001
