# Characterization harness

Short, fully pinned (seed + flags) smoke runs that catch unintended changes to
LASSI’s Monte Carlo trajectory.

## One command

From the repository root:

```bash
make characterization
```

That builds via `src/compile.sh`, runs the smoke decks, and **exits non-zero
(CI-style) if any baseline artifact differs bit-for-bit**.

## What is pinned

Each deck’s `param.key` sets:

- A **non-zero** `RANDOM_SEED` (seed `0` would call `time(NULL)` and is rejected
  by the harness).
- Explicit values for temperature / annealing, bias, all MC move frequencies
  (including `MV_PR_*`), reporting frequencies, `REPORT_CONFIG_MODE`, and
  `ANALYSIS_CLUSTER_MODE` — no reliance on ambient keyfile defaults.
- Analysis writers that are unused here are set to frequency `0`.

The harness itself also pins `LC_ALL=C`, `LANG=C`, and `TZ=UTC`, and checks the
run log for `RNG Seed = <expected>`.

Build flags are those in `src/compile.sh`
(`gcc -march=x86-64 -O3 -std=gnu11 ... -lm -fcommon`).

## Decks

| Deck | Based on | Seed | Size | Steps |
|------|----------|------|------|-------|
| `decks/smoke_single` | Example A (single diblock) | `42` | 1×20 beads, box 40 | PREEQ 500 + N_STEPS 2000, 1 cycle |
| `decks/smoke_multi` | Example C (separated types) | `12345` | 2+2×10 beads, box 40 | same short schedule |

Compared artifacts (per deck), under `baselines/<deck>/`:

- `*_trj.lammpstrj`, `*_topo.lammpstrj`, `*_0_restart.lammpstrj`, `*_EQ_restart.lammpstrj`
- `*_0_energy.dat`, `*_energy.dat`, `*_0_mcmove.dat`, `*_mcmove.dat`

Stdout / `run.log` are **not** compared (they include wall-clock timings).

## Determinism

With the pinned non-zero seeds and flags above, these decks produce
**bit-identical** outputs on repeated runs using the same `src/compile.sh`
binary on the same platform (Linux glibc `rand()`/`srand()`).

`rand()` can differ across C libraries (glibc vs macOS/BSD). Baselines here were
generated on Linux with `src/compile.sh`.

## Updating baselines

Only after an **intentional** reference change:

```bash
make characterization-update
```

Then review and commit the diff under `tests/characterization/baselines/`.
