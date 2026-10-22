# Thin top-level targets. The simulation itself still builds via src/compile.sh.

.PHONY: lassi characterization characterization-update clean-characterization unit-lattice unit-legal unit-energy unit-energy-aniso

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

# Pure linker-distance / bead-list membership helpers (Phase 5.2); no decks.
unit-legal:
	gcc -std=gnu11 -O0 -Wall -Wextra -I. \
		tests/unit/test_linker_legal.c src/structure_legal.c src/structure_lattice.c \
		src/structure.c src/globals.c \
		-o tests/unit/test_linker_legal -lm -fcommon
	./tests/unit/test_linker_legal

# Pure isotropic OVLP/CONT/fSol kernels (Phase 5.3); no decks.
unit-energy:
	gcc -std=gnu11 -O0 -Wall -Wextra -I. \
		tests/unit/test_energy_iso.c src/energy_iso.c src/globals.c \
		-o tests/unit/test_energy_iso -lm -fcommon
	./tests/unit/test_energy_iso

# Pure anisotropic SC-SC kernels (Phase 5.4); no decks.
unit-energy-aniso:
	gcc -std=gnu11 -O0 -Wall -Wextra -I. \
		tests/unit/test_energy_aniso.c src/energy_aniso.c src/structure_legal.c \
		src/structure_lattice.c src/structure.c src/globals.c \
		-o tests/unit/test_energy_aniso -lm -fcommon
	./tests/unit/test_energy_aniso

