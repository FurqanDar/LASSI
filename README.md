<img src="https://github.com/Pappulab/LASSI/blob/master/docs/_imgs/LaSSI_Logo.png" alt="LaSSI" width="350"/>

LASSI - LAttice Simulation Engine For Stickers And Spacer Interactions
==========================================================================

LASSI is an engine designed for simulating the phase transitions
of polymers with specific anisotropic interactions.

For information about the underlying physics and an applied example, please see this 
[paper](https://doi.org/10.1371/journal.pcbi.1007028).

LASSI was written by Furqan Dar and Jeong-Mo Choi in the [Pappu Lab](http://pappulab.wustl.edu/).

Please address all questions to [Furqan](http://pappulab.wustl.edu/people.html#grads).

## Characterization / regression harness

Pinned-seed smoke decks under `tests/characterization/`. One command from the
repo root:

```bash
make characterization
```

Exits non-zero (CI-style) if trajectory / energy / MC-move outputs drift from
checked-in baselines. See
[`tests/characterization/README.md`](tests/characterization/README.md).
