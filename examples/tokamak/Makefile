EQDSK := EQDSK_I.geqdsk
URL   := https://crppwww.epfl.ch/~sauter/benchmark/EQDSK_I
RUN   := ../../build/simple.x

.PHONY: all run clean

all: $(EQDSK)

$(EQDSK):
	@printf 'Downloading %s...\n' $(EQDSK)
	@curl -L --fail --show-error --output $(EQDSK).tmp $(URL)
	@mv $(EQDSK).tmp $(EQDSK)
	@printf 'Saved %s\n' $(EQDSK)

run: $(EQDSK)
	@$(RUN) simple.in

clean:
	@rm -f $(EQDSK)
