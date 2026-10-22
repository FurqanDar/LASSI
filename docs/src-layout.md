# Source layout (after Phase 1 splits)

Brief map of where C modules live under `src/`. Public APIs stay in the
existing headers (`mcmove.h`, `structure.h`, `print.h`, …). No behavior change
is implied by this split.

## Entry / config

| Files | Role |
|-------|------|
| `main.c` | Program entry and run loop |
| `parsekey.c` / `parsekey.h` | Keyfile parsing |
| `initialize.c` / `initialize.h` | System setup |
| `global.h` | Shared globals and constants |

## Monte Carlo moves (`mcmove.h`)

| Files | Role |
|-------|------|
| `mcmove.c` | Thin `MC_Step` / `MC_Step_Therm` over `MV_*` runner tables |
| `mcmove_moves.c` | Production `Move_*` |
| `mcmove_moves_equil.c` | Thermalization `Move_*_Equil` |
| `mcmove_ops.c` | `OP_*` lattice/bead helpers, Rosenbluth, MH |
| `mcmove_rot.c` | `Rot_*`, rot-state checks, shuffle |

## Lattice / structure (`structure.h`)

| Files | Role |
|-------|------|
| `structure_lattice.c` | Lattice index, min-image distances, `LatPos_*`, `Vec3n_*` (pure cores + `g_system` wrappers) |
| `structure_legal.c` | Linker-constraint / bead-list legality (`LinkerConstraint_DistOK`, `Check_Linker*`, `Check_BeadID_InList`) |
| `structure.c` | Topology helpers, `LatticeUtil_*`, bead angle helpers |
| `structure_sanity.c` | Runtime sanity orchestration (uses legal/lattice helpers) |
| `structure_neighbor.c` | `NeighborSearch_*` |
| `structure_lists.c` | Bead/chain/list helper ops |

## Analysis / I/O

| Files | Role |
|-------|------|
| `cluster.c` / `cluster.h` | Cluster finding and analysis |
| `cluster_mc.h` | Move-facing cluster queries (see architecture.md) |
| `analysis.c` | RDF, radial density, COM, gyration (decls still in `structure.h`) |
| `print.c` | `DataPrinting_*` / `DataAnalysis_*` orchestration |
| `print_screen.c` | Screen / log / sanity-fail messages |
| `print_fileio.c` | Trajectories, energy/MC files, restarts, `CopyData_*` |

## Other

| Files | Role |
|-------|------|
| `energy.c` / `energy.h` | Energy evaluation (system/chain totals, anisotropic, Iso For* move deltas) |
| `energy_iso.c` | Isotropic OVLP/CONT/fSol pair kernels + empty-site solvation (pure FromCoupling + wrappers) |

Build: `src/compile.sh` (also wired in top-level `Makefile` / `CMakeLists.txt`).
Regression gate: `make characterization`.

## Dependency direction

See [`docs/architecture.md`](architecture.md) for the allowed
`structure ← energy ← moves ← analysis/I/O ← main/parse` include edges.
