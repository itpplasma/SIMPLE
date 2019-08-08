#plot "euler1_quasi.out"   u 5 ps .001 
#replot "euler2_quasi.out" u 5 ps .001

#plot "gauss2_quasi.out" u 5 ps .001  
#replot "midpoint_quasi.out" u 5 ps .001  
#replot   "verlet_quasi.out" u 5 ps .001

plot "mclachlan4_quasi.out" u 5 ps .001
replot "blanes4_quasi.out" u 5 ps .001
#replot "order4_quasi.out" u 5 ps .001
replot "gauss4_quasi.out"   u 5 ps .001

