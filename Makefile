# Thin top-level targets. The simulation itself still builds via src/compile.sh.

.PHONY: lassi characterization characterization-update clean-characterization

lassi:
	cd src && bash compile.sh

characterization: lassi
	LASSI_BIN="$(CURDIR)/src/lassi" bash tests/characterization/run_characterization.sh

characterization-update: lassi
	LASSI_BIN="$(CURDIR)/src/lassi" bash tests/characterization/run_characterization.sh --update-baselines

clean-characterization:
	rm -rf tests/characterization/.work
