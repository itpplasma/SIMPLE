FC       = gfortran

#DEBUGFLAG += -g -ggdb -C -p -pg -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -fbounds-check
DEBUGFLAG += -g -ggdb -C -p -pg -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check
#OPTS= -J OBJS -O
OPTS= -J OBJS $(DEBUGFLAG) -cpp -fopenmp # comment out -fopenmp for single-core version
NCINC = -I/usr/include
NCLIB = -L/usr/lib -lnetcdff -lnetcdf

OBJS =  OBJS/canonical_coordinates_mod.o \
	OBJS/nctools_module.o \
	OBJS/odeint_allroutines.o \
	OBJS/magfie.o \
	OBJS/chamb_m.o \
	OBJS/sub_alpha_lifetime_can.o \
	OBJS/vmecinm_m.o \
	OBJS/spline_vmec_data.o \
	OBJS/spl_three_to_five.o \
	OBJS/new_vmec_allocation_stuff.o \
	OBJS/get_canonical_coordinates.o \
	OBJS/zzg.o \
	OBJS/field_can.o \
	OBJS/orbit_symplectic.o \
	OBJS/binsrc.o \
	OBJS/common.o \
	OBJS/minpack.o \
	OBJS/alpha_lifetime_c.o

alpha_lifetime_vmec_c.x: $(OBJS) Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -o alpha_lifetime_vmec_c.x $(OBJS) $(NCLIB) -llapack
OBJS/canonical_coordinates_mod.o: SRC/canonical_coordinates_mod.f90 Alpha_lifetime_vmec_c.mk 
	$(FC) $(OPTS) -c SRC/canonical_coordinates_mod.f90
	mv canonical_coordinates_mod.o OBJS/
OBJS/nctools_module.o: SRC/nctools_module.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) $(NCINC) -c SRC/nctools_module.f90 
	mv nctools_module.o OBJS/
OBJS/odeint_allroutines.o: SRC/odeint_allroutines.f Alpha_lifetime_vmec_c.mk 
	$(FC) $(OPTS) -c SRC/odeint_allroutines.f
	mv odeint_allroutines.o OBJS/
OBJS/magfie.o: SRC/magfie.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/magfie.f90
	mv magfie.o OBJS/
OBJS/chamb_m.o: SRC/chamb_m.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90 
	$(FC) $(OPTS) -c SRC/chamb_m.f90
	mv chamb_m.o OBJS/
OBJS/sub_alpha_lifetime_can.o: SRC/sub_alpha_lifetime_can.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/sub_alpha_lifetime_can.f90
	mv sub_alpha_lifetime_can.o OBJS/
OBJS/vmecinm_m.o: SRC/vmecinm_m.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/vmecinm_m.f90
	mv vmecinm_m.o OBJS/
OBJS/spline_vmec_data.o: SRC/spline_vmec_data.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/spline_vmec_data.f90
	mv spline_vmec_data.o OBJS/
OBJS/spl_three_to_five.o: SRC/spl_three_to_five.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/spl_three_to_five.f90
	mv spl_three_to_five.o OBJS/
OBJS/new_vmec_allocation_stuff.o: SRC/new_vmec_allocation_stuff.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/new_vmec_allocation_stuff.f90
	mv new_vmec_allocation_stuff.o OBJS/
OBJS/get_canonical_coordinates.o: SRC/get_canonical_coordinates.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/get_canonical_coordinates.f90
	mv get_canonical_coordinates.o OBJS/
OBJS/zzg.o: SRC/zzg.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/zzg.f90
	mv zzg.o OBJS/
OBJS/field_can.o: SRC/field_can.f90 Canonical_coordinates.mk
	$(FC) $(OPTS) -c SRC/field_can.f90
	mv field_can.o OBJS/
OBJS/orbit_symplectic.o: SRC/orbit_symplectic.f90 Canonical_coordinates.mk SRC/field_can.f90
	$(FC) $(OPTS) -c SRC/orbit_symplectic.f90
	mv orbit_symplectic.o OBJS/
OBJS/binsrc.o: SRC/binsrc.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/binsrc.f90
	mv binsrc.o OBJS/
OBJS/common.o: SRC/common.f90 Canonical_coordinates.mk
	$(FC) $(OPTS) -c SRC/common.f90
	mv common.o OBJS/
OBJS/minpack.o: SRC/minpack.f90 Canonical_coordinates.mk
	$(FC) $(OPTS) -c SRC/minpack.f90
	mv minpack.o OBJS/
OBJS/alpha_lifetime_c.o: SRC/alpha_lifetime_c.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
	$(FC) $(OPTS) -c SRC/alpha_lifetime_c.f90
	mv alpha_lifetime_c.o OBJS/
#OBJS/.o: SRC/.f90 Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
#	$(FC) $(OPTS) -c SRC/.f90
#	mv .o OBJS/
#OBJS/.o: SRC/.f Alpha_lifetime_vmec_c.mk SRC/canonical_coordinates_mod.f90
#	$(FC) $(OPTS) -c SRC/.f
#	mv .o OBJS/
