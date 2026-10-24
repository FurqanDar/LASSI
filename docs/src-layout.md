# Source layout (after Phase 1 splits)

Brief map of where C modules live under `src/`. Public APIs stay in the
existing headers (`mcmove.h`, `structure.h`, `print.h`, …). No behavior change
is implied by this split.

## Entry / config

| Files | Role |
|-------|------|
| `main.c` | Program entry and run loop |
| `parsekey.c` / `parsekey.h` | Keyfile parsing → `SimulationConfig` (`g_config`); `Config_ApplyToRuntime` shims into domain globals |
| `initialize.c` / `initialize.h` | System setup; AtStart / initial-conditions entrypoints take `const SimulationConfig*` |
| `global.h` / `globals.c` | Shared globals, domain structs, `SimulationConfig` (scratch vs durable tags — see below) |

## Monte Carlo moves (`mcmove.h` / `mcmove_internal.h`)

| Files | Role |
|-------|------|
| `mcmove.h` | Public API: `MC_Step` / `MC_Step_Therm` only |
| `mcmove_internal.h` | Cross-TU `Move_*` / `OP_*` / rot helpers (mcmove*.c only) |
| `mcmove.c` | Thin `MC_Step` / `MC_Step_Therm` over `MV_*` runner tables |
| `mcmove_moves.c` | Production `Move_*` |
| `mcmove_moves_equil.c` | Thermalization `Move_*_Equil` |
| `mcmove_ops.c` | `OP_*` lattice/bead helpers, Rosenbluth, MH |
| `mcmove_rot.c` | `Rot_*` (file-local), rot-state checks, shuffle |

## Lattice / structure (`structure.h` / `structure_internal.h`)

| Files | Role |
|-------|------|
| `structure.h` | Public lattice / legal / neighbor / list / analysis-index API |
| `structure_internal.h` | Family-only cross-TU helpers (structure*.c only) |
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
| `analysis.c` | RDF, radial density, COM, gyration (public entrypoints still in `structure.h`) |
| `print.c` | `DataPrinting_*` / `DataAnalysis_*` orchestration |
| `print_screen.c` | Screen / log / sanity-fail messages |
| `print_fileio.c` | Trajectories, energy/MC files, restarts, `CopyData_*` |

## Other

| Files | Role |
|-------|------|
| `energy.h` | Public API: system total, bias/topo, Iso For* move deltas, pure FromCoupling cores |
| `energy_internal.h` | Family-only iso/aniso wrappers (energy*.c only) |
| `energy.c` | Totals, bias/topo, Iso For* move deltas; file-local neigh-list / isotropic helpers |
| `energy_iso.c` | Isotropic OVLP/CONT/fSol pair kernels + empty-site solvation |
| `energy_aniso.c` | Anisotropic SC-SC bond kernels |

Build: `src/compile.sh` (also wired in top-level `Makefile` / `CMakeLists.txt`).
Regression gate: `make characterization`.

## Scratch vs durable

See [`docs/scratch-vs-durable.md`](scratch-vs-durable.md) (Phase 7.1 inventory;
Phase 7.2 named nests: `scratch` / `cycle` / `totals`).

## Dependency direction

See [`docs/architecture.md`](architecture.md) for the allowed
`structure ← energy ← moves ← analysis/I/O ← main/parse` include edges.
