# Characterization harness

Short, fixed-seed smoke runs that catch unintended changes to LASSI’s Monte Carlo
trajectory. This is a characterization / regression check, not a physics test suite.

## Quick start

From the repository root:

```bash
make characterization
```

Or:

```bash
bash tests/characterization/run_characterization.sh
```

The script builds `src/lassi` (unless `LASSI_BIN` points at an existing binary),
runs the smoke decks, and **fails if any checked-in baseline file differs
bit-for-bit**.

## Decks

| Deck | Based on | Seed | Size | Steps |
|------|----------|------|------|-------|
| `decks/smoke_single` | Example A (single diblock) | `RANDOM_SEED 42` | 1×20 beads, box 40 | PREEQ 500 + N_STEPS 2000, 1 cycle |
| `decks/smoke_multi` | Example C (separated types) | `RANDOM_SEED 12345` | 2+2×10 beads, box 40 | same short schedule |

Compared artifacts (per deck), under `baselines/<deck>/`:

- `*_trj.lammpstrj`, `*_topo.lammpstrj`, `*_0_restart.lammpstrj`, `*_EQ_restart.lammpstrj`
- `*_0_energy.dat`, `*_energy.dat`, `*_0_mcmove.dat`, `*_mcmove.dat`

Stdout / `run.log` are **not** compared (they include wall-clock timings).

## Determinism

With a non-zero `RANDOM_SEED`, these smoke decks produce **bit-identical**
outputs on repeated runs using the same `src/compile.sh` binary on the same
platform (verified on Linux glibc with the stock `rand()`/`srand()` path).

Notes:

- `RANDOM_SEED 0` means “seed from `time(NULL)`” and is **not** reproducible.
- `rand()`/`srand()` behavior can differ across C libraries (e.g. glibc vs
  macOS/BSD). Baselines checked in here were generated on Linux with
  `src/compile.sh`. Cross-OS drift is an environment issue, not a simulation
  physics change.
- No simulation-code change is required for these smoke cases to be
  deterministic given a fixed seed.

## Updating baselines

Only after an **intentional** behavior change (or a deliberate reference refresh):

```bash
make characterization-update
```

Then review the diff under `tests/characterization/baselines/` and commit it
with the change that caused it.
