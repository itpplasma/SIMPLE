orbit_num = ARG1

# Print orbit_num as overall title
set title 'Orbit '.orbit_num

# Plot the data using the variable
plot 'canflux/fort.'.orbit_num u 2 lt rgb 'dark-green' w l \
    title 'Canonical', \
    'boozer/fort.'.orbit_num u 2 lt rgb 'black' w l \
    title 'Boozer', \
    'meiss/fort.'.orbit_num u 2 w l lt rgb 'blue' lw 2 \
    title 'Meiss', \
    'albert/fort.'.orbit_num u 2 w l lt rgb 'red' lw 1 \
    title 'Albert'

pause -1
