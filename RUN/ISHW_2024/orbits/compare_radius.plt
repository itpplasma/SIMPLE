orbit_num = ARG1

# Print orbit_num as overall title
set title 'Orbit '.orbit_num

# Plot the data using the variable
plot 'vmec/fort.'.orbit_num u 1:2 w l lt rgb 'dark-red'  \
    title 'VMEC', \
    'canflux/fort.'.orbit_num u 1:2 w l  lt rgb 'dark-green'  \
    title 'Canonical', \
    'boozer/fort.'.orbit_num u 1:2 w l  lt rgb 'black' \
    title 'Boozer', \
    'meiss/fort.'.orbit_num u 1:2 w l lt rgb 'blue' \
    title 'Meiss', \
    'albert/fort.'.orbit_num u 1:2 w l lt rgb 'red' \
    title 'Albert'

pause -1
