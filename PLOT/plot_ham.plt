plot "euler1.out"   u 5 ps .001
replot "midpoint.out" u 5 ps .001
replot "mclachlan4.out" u 5 ps .001

# plot "euler1.out" u (1+$1*cos($2)):($1*sin($2)) ps .01
# replot "verlet.out" u (1+$1*cos($2)):($1*sin($2)) ps .01
