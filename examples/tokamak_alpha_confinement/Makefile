RUN := ../../build/simple.x

.PHONY: all run clean

all: run

run:
	@echo "Running alpha particle confinement test..."
	@$(RUN) simple.in
	@echo "Expected: Zero particles lost in axisymmetric tokamak"

clean:
	@rm -f fort.* *.dat *.nc
	@echo "Cleaned output files"
