# Thin top-level targets. The simulation itself still builds via src/compile.sh.

.PHONY: lassi characterization characterization-update clean-characterization unit-lattice

lassi:
	cd src && bash compile.sh

characterization: lassi
	LASSI_BIN="$(CURDIR)/src/lassi" bash tests/characterization/run_characterization.sh

characterization-update: lassi
	LASSI_BIN="$(CURDIR)/src/lassi" bash tests/characterization/run_characterization.sh --update-baselines

clean-characterization:
	rm -rf tests/characterization/.work

# Pure lattice-index / min-image distance helpers (Phase 5.1); no decks.
unit-lattice:
	gcc -std=gnu11 -O0 -Wall -Wextra -I. \
		tests/unit/test_lattice_index.c src/structure_lattice.c src/globals.c \
		-o tests/unit/test_lattice_index -lm -fcommon
	./tests/unit/test_lattice_index

