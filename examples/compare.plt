orbit_num = ARG1
phi0 = 3.2038182147970
R0 = 5.0

# Print orbit_num as overall title
set title 'Orbit '.orbit_num

set size square

# Plot the data using the variable
splot 'canflux/fort.'.orbit_num u ((R0+$2*cos($3))*cos($4)):((R0+$2*cos($3))*sin($4)):($2*sin($3)) w l lt rgb 'dark-green' \
    title 'Canonical', \
    'boozer/fort.'.orbit_num u ((R0+$2*cos($3))*cos($4+phi0)):((R0+$2*cos($3))*sin($4+phi0)):($2*sin($3)) w l  lt rgb 'black'\
    title 'Boozer', \
    'meiss/fort.'.orbit_num u ((R0+$2*cos($3))*cos($4)):((R0+$2*cos($3))*sin($4)):($2*sin($3)) w l lt rgb 'blue'  \
    title 'Meiss', \
    'albert/fort.'.orbit_num u ((R0+$2*cos($3))*cos($4)):((R0+$2*cos($3))*sin($4)):($2*sin($3)) w l lt rgb 'red' \
    title 'Albert'

pause -1
