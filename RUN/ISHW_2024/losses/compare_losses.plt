set title 'Confined fraction'

# Plot the data using the variable
plot 'vmec/confined_fraction.dat' u 1:($2+$3) w lp lt rgb 'dark-red'  \
    title 'VMEC', \
    'canflux/confined_fraction.dat' u 1:($2+$3) w lp lt rgb 'dark-green'  \
    title 'Canonical', \
    'meiss/confined_fraction.dat' u 1:($2+$3) w l lt rgb 'blue' \
    title 'Meiss Canonical', \
    'albert/confined_fraction.dat' u 1:($2+$3) w l lt rgb 'red' \
    title 'Meiss RK4/5', \
    'meiss512/confined_fraction.dat' u 1:($2+$3) w l lt rgb 'green' \
    title 'Meiss more accurate'
pause -1
