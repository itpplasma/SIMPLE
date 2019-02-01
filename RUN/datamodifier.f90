implicit none
double precision:: a1,a2,a3,a4,a5,a6
integer:: i

OPEN(1, FILE='orbit_vmec.out', recl=1023)
OPEN(2,FILE='orbit_vmec_plot.out',recl=1023)
do i=1,100
    READ(1,*) a1,a2,a3,a4,a5,a6
    WRITE (2, *) a2*cos(a3)*cos(a4),a2*cos(a3)*sin(a4)
enddo
CLOSE(1)

end
