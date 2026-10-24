# Performance bench

Small, deterministic wall-clock bench for **performance regression awareness**.
Not part of required CI; use locally (`make bench` / `make bench-compare`) or the
optional `workflow_dispatch` job.

## One command

```bash
make bench              # build Release-ish, run deck N times, print median
make bench-compare      # same + fail if median > baseline × (1 + threshold)
make bench-update       # same + rewrite tests/bench/baseline.json
```

Defaults: `BENCH_RUNS=3`, `BENCH_THRESHOLD_PCT=15`. Override as env vars.

## Deck: `decks/bench_diblock`

| Field | Value |
|-------|-------|
| Based on | Example B (`B_DiBlockSystem`) |
| System | 500 diblock chains × 20 beads (A10–B10), box 200 |
| Density | ~1.25×10⁻³ monomers / site |
| Seed | `20241024` (non-zero; harness rejects `0`) |
| Schedule | `PREEQ_STEPS` 1 000 000 + `N_STEPS` 10 000 000, 1 cycle, no annealing |
| Moves | Same frequencies as characterization multi / Example B production set |
| I/O | Trajectory / energy / analysis writers off (`*_FREQ` 0); sparse log |

Larger than characterization smoke (few thousand steps, few beads), still short:
about **tens of seconds per run** on a typical cloud VM (≈1–2 minutes for 3 runs).
Wall time is hardware-bound; the checked-in baseline is a reference point, not a
guarantee across machines.

## Build flags

Same as characterization: `src/compile.sh`
(`gcc -march=x86-64 -O3 -std=gnu11 ... -lm -fcommon`).

## Baseline

[`baseline.json`](baseline.json) stores median wall seconds plus metadata (host,
compiler, build flags, seed, step counts, how recorded). Refresh after an
intentional perf-relevant change or when re-homing the reference machine:

```bash
make bench-update
```

## Compare / threshold

```bash
make bench-compare
# or
BENCH_THRESHOLD_PCT=20 make bench-compare
```

Exits **non-zero** if current median is slower than
`baseline.median_seconds × (1 + threshold/100)`. Start loose (15%) so optional CI
does not flake on noisy VMs.

## Optional: `perf` cycles

This harness records **wall-clock** only. Local cycle counts are optional, e.g.:

```bash
perf stat -e cycles,instructions make bench
```

Not checked in; not used by `--compare`.

## Optional CI

[`.github/workflows/bench.yml`](../../.github/workflows/bench.yml) is
**manual only** (`workflow_dispatch`). It runs `make bench-compare` and is not
required on push/PR. Characterization remains the behavior gate.
