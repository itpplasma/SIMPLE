# Standalone prototype; not part of the SIMPLE CMake build.
FC = gfortran
FFLAGS ?= -O2 -march=native -fopenmp -Wall -Wextra -fimplicit-none -std=f2018
BUILD := build

MODULES := sheet_field.f90 fo_sheet.f90 sheet_setup.f90 rk_reference.f90
OBJECTS := $(MODULES:%.f90=$(BUILD)/%.o)
PROGRAMS := test_fo_symplectic test_fo_stationary study_mu_scattering

all: $(PROGRAMS:%=$(BUILD)/%.x)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: %.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -c $< -o $@

$(BUILD)/fo_sheet.o: $(BUILD)/sheet_field.o
$(BUILD)/sheet_setup.o: $(BUILD)/sheet_field.o
$(BUILD)/rk_reference.o: $(BUILD)/sheet_field.o

$(BUILD)/%.x: %.f90 $(OBJECTS)
	$(FC) $(FFLAGS) -I$(BUILD) $< $(OBJECTS) -o $@

test: all
	$(BUILD)/test_fo_symplectic.x
	$(BUILD)/test_fo_stationary.x

clean:
	rm -rf $(BUILD)

.PHONY: all test clean
