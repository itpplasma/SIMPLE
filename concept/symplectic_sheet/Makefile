# Standalone prototype; not part of the SIMPLE CMake build.
FC = gfortran
FFLAGS ?= -O2 -march=native -fopenmp -Wall -Wextra -fimplicit-none -std=f2018
ROOT := concept/symplectic_sheet
BUILD := $(ROOT)/build

MODULES := sheet_field.f90 fo_sheet.f90 sheet_setup.f90 rk_reference.f90 gc_hybrid.f90
OBJECTS := $(addprefix $(BUILD)/,$(MODULES:%.f90=%.o))
PROGRAMS := study_mu_scattering study_hybrid

all: $(PROGRAMS:%=$(BUILD)/%.x)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: $(ROOT)/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -c $< -o $@

$(BUILD)/fo_sheet.o: $(BUILD)/sheet_field.o
$(BUILD)/sheet_setup.o: $(BUILD)/sheet_field.o
$(BUILD)/rk_reference.o: $(BUILD)/sheet_field.o
$(BUILD)/gc_hybrid.o: $(BUILD)/sheet_field.o

$(BUILD)/%.x: $(ROOT)/%.f90 $(OBJECTS)
	$(FC) $(FFLAGS) -I$(BUILD) $< $(OBJECTS) -o $@


clean:
	rm -rf $(BUILD)

.PHONY: all clean
