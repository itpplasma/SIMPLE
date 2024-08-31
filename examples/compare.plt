orbit_num = ARG1
R0 = 5.0

# Print orbit_num as overall title
set title 'Orbit '.orbit_num

set size square

# Plot the data using the variable
splot 'vmec/fort.'.orbit_num u ((R0+$2*cos($3))*cos($4)):((R0+$2*cos($3))*sin($4)):($2*sin($3)) w l lt rgb 'dark-red' lw 2 \
    title 'VMEC', \
    'canflux/fort.'.orbit_num u ((R0+$2*cos($3))*cos($4)):((R0+$2*cos($3))*sin($4)):($2*sin($3)) w l lt rgb 'dark-green' lw 1.5 \
    title 'Canonical', \
    'boozer/fort.'.orbit_num u ((R0+$2*cos($3))*cos($4)):((R0+$2*cos($3))*sin($4)):($2*sin($3)) w l  lt rgb 'black' lw 1\
    title 'Boozer', \
    'meiss/fort.'.orbit_num u ((R0+$2*cos($3))*cos($4)):((R0+$2*cos($3))*sin($4)):($2*sin($3)) w l lt rgb 'blue' lw 2\
    title 'Meiss', \
    'albert/fort.'.orbit_num u ((R0+$2*cos($3))*cos($4)):((R0+$2*cos($3))*sin($4)):($2*sin($3)) w l lt rgb 'red' lw 1\
    title 'Albert'

pause -1
