<img src="https://github.com/Pappulab/LASSI/blob/master/docs/_imgs/LaSSI_Logo.png" alt="LaSSI" width="350"/>

LASSI - LAttice Simulation Engine For Stickers And Spacer Interactions
==========================================================================

LASSI is an engine designed for simulating the phase transitions
of polymers with specific anisotropic interactions.

For information about the underlying physics and an applied example, please see this 
[paper](https://doi.org/10.1371/journal.pcbi.1007028).

LASSI was written by Furqan Dar and Jeong-Mo Choi in the [Pappu Lab](http://pappulab.wustl.edu/).

Please address all questions to [Furqan](http://pappulab.wustl.edu/people.html#grads).

## Source layout

After Phase 1 file splits, modules under `src/` are grouped by concern
(moves, structure, I/O). See [`docs/src-layout.md`](docs/src-layout.md).

## Code formatting

C/H sources under `src/` and `tests/` are formatted with
[clang-format](https://clang.llvm.org/docs/ClangFormat.html) using the root
[`.clang-format`](.clang-format) (LLVM-based). CI runs a dry-run check on every
push and pull request to `master` / `main` (see
[`.github/workflows/format.yml`](.github/workflows/format.yml)). Local
pre-commit hooks are optional later; not required today.

```bash
# Reformat project C/H (same scope as CI)
find src tests -type f \( -name '*.c' -o -name '*.h' \) -print0 | xargs -0 clang-format -i

# Check without writing (CI equivalent)
find src tests -type f \( -name '*.c' -o -name '*.h' \) -print0 | xargs -0 clang-format --dry-run --Werror
```

## Characterization / regression harness

Pinned-seed smoke decks under `tests/characterization/`. One command from the
repo root:

```bash
make characterization
```

Exits non-zero (CI-style) if trajectory / energy / MC-move outputs drift from
checked-in baselines. GitHub Actions runs `make characterization` on every push
and pull request to `master` / `main`. See
[`tests/characterization/README.md`](tests/characterization/README.md).
