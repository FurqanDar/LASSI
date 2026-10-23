# Scratch vs durable state (Phase 7.1 inventory + 7.2 nests)

Ownership/lifetime tags for process-wide buffers. Phase 7.2 nests scratch /
cycle / totals under named members; **no intentional behavior change**.

## Definitions

| Kind | Meaning |
|------|---------|
| **Durable** | Meaningful across many MC steps / cycles: topology, lattice occupancy, energy tables, config-derived control, cycle accumulators (`ldaTOT*`). |
| **Scratch** | Ephemeral working storage: contents valid only for the current move, neighbor search, or analysis pass. Allocation may still live for the whole run. |
| **Cycle buffer** | Cleared or rolled at cycle boundaries; holds in-progress histograms before fold into `ldaTOT*`. |

## Nest layout (Phase 7.2)

| Domain | Nest member | Typedef | Kind |
|--------|-------------|---------|------|
| `g_system` | `.scratch` | `SystemScratch` | Per-move / per-search scratch |
| `g_energy` | `.scratch` | `EnergyScratch` | Rotational / Boltzmann scratch |
| `g_mc` | `.cycle` | `McCycleScratch` | Cycle accept matrix |
| `g_reports` | `.cycle` | `ReportsCycleScratch` | Cycle traj accumulation |
| `g_analysis` | `.cycle` | `AnalysisCycleScratch` | Cycle working sets / histograms |
| `g_analysis` | `.totals` | `AnalysisTotals` | Durable cross-cycle `ldaTOT*` |

Call sites use nested access, e.g. `g_system.scratch.old_bead`,
`g_analysis.cycle.naList`, `g_analysis.totals.ldaTOTRDF_Arr`.

## By domain

### `SystemState` (`g_system`)

| Fields | Kind | Notes |
|--------|------|-------|
| `bead_info`, `chain_info`, `topo_info`, `linker_len`, `naTotLattice` | Durable | Configuration / lattice occupancy |
| counts, box, mode flags, capability flags, linker floats | Durable | |
| `naLocalArr` | Durable (setup) | Fixed neighbor offsets filled once at init |
| `scratch.old_bead` | Scratch contents | Trial backup for reject/restore; allocated for the run |
| `scratch.naOld/NewOvlpNeighs`, `scratch.naOld/NewContNeighs`, `scratch.naTempR` | Scratch | Neighbor-search / trial coordinates |
| `scratch.naRot_IndArr` | Scratch contents | Rotation shuffle indices |

### `EnergyTables` (`g_energy`)

| Fields | Kind | Notes |
|--------|------|-------|
| `faEnergy`, `faEnRad`, `nBeadTypes` | Durable | Parsed tables |
| `faCurrEn` | Durable runtime | Current system energy vector |
| `scratch.naRotTrial`, `scratch.ldaBoltzFac*`, `scratch.ldaBoltzNorm` | Scratch | Rotational trial / Boltzmann temps |

### `McControl` (`g_mc`)

| Fields | Kind | Notes |
|--------|------|-------|
| temps, freqs, step counts, seed, `faKT_Cycle` | Durable control | From config / init |
| `fCuTemp`, `fSquishRad_Sq` | Durable runtime | |
| `cycle.naMCAccepMat` | Cycle buffer | Zeroed in `Reset_McCycleScratch` |

### `Reports` (`g_reports`)

| Fields | Kind | Notes |
|--------|------|-------|
| paths, freqs, `nTrajMode` | Durable | From config |
| `cycle.naTOTTRAJ_Arr`, `cycle.nTrajCurFrame` | Cycle buffer | Experimental traj accumulation |

### `AnalysisState` (`g_analysis`)

| Fields | Kind | Notes |
|--------|------|-------|
| `nClusteringMode`, size/comp metadata | Durable | |
| `cycle.naClusterMatrix`, `cycle.naList`, `cycle.laClusHistList`, `cycle.naChainCheckList` | Scratch / cycle | Cluster pass working sets |
| `cycle.ldaRDF_Arr`, `cycle.ldaRadDen_Arr`, `cycle.ldaMOLCLUS_Arr`, `cycle.faGyrTensor` | Cycle buffer | Current-cycle histograms |
| `totals.ldaTOTRDF_Arr`, `totals.ldaTOTRadDen_Arr`, `totals.ldaTOTCLUS_Arr`, `totals.ldaTOTMOLCLUS_Arr`, `totals.ldaTOTRg_Arr` | Durable accumulators | Across cycles |
| `cycle.nTot*Counter`, `cycle.nLargestClusterRightNow` | Cycle counters | Cleared by `Reset_Counters` |

### `SimulationConfig` (`g_config`)

Durable **inputs** only (keyfile). Not runtime scratch.

## Reset hooks (cycle boundary)

Called from the run loop via `Reset_Global_Arrays` between cycles:

| Helper | Clears |
|--------|--------|
| `Reset_AnalysisCycleScratch` | `g_analysis.cycle` working sets / histograms |
| `Reset_Counters` | `g_analysis.cycle` counters + `g_reports.cycle.nTrajCurFrame` |
| `Reset_McCycleScratch` | `g_mc.cycle.naMCAccepMat` |
| `Reset_Global_Arrays` | Orchestrates the three above (same order/effect as pre-7.2) |

No process-exit teardown of nest allocations (still owned for the whole run).
