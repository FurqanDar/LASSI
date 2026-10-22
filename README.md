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

A short fixed-seed smoke harness lives under `tests/characterization/`. From the
repo root, after a normal build:

```bash
make characterization
```

This runs two minimal decks and fails if trajectory / energy / MC-move outputs
differ bit-for-bit from checked-in baselines. See
[`tests/characterization/README.md`](tests/characterization/README.md) for
details, determinism notes, and how to refresh baselines.
