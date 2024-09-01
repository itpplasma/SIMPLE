set title 'Confined fraction'

# Plot the data using the variable
plot 'vmec/confined_fraction.dat' u 1:($2+$3) w l lt rgb 'dark-red'  \
    title 'VMEC', \
    'canflux/confined_fraction.dat' u 1:($2+$3) w l lt rgb 'dark-green'  \
    title 'Canonical', \
    'meiss/confined_fraction.dat' u 1:($2+$3) w l lt rgb 'blue' \
    title 'Meiss'

pause -1
