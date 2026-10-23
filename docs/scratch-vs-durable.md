# Scratch vs durable state (Phase 7.1 inventory)

Ownership/lifetime tags for process-wide buffers. **No behavior change** in
7.1 — see section comments (and light anonymous groups) in `src/global.h`.
Phase 7.2+ may nest further and add explicit reset/teardown between cycles.

## Definitions

| Kind | Meaning |
|------|---------|
| **Durable** | Meaningful across many MC steps / cycles: topology, lattice occupancy, energy tables, config-derived control, cycle accumulators (`ldaTOT*`). |
| **Scratch** | Ephemeral working storage: contents valid only for the current move, neighbor search, or analysis pass. Allocation may still live for the whole run. |
| **Cycle buffer** | Cleared or rolled at cycle boundaries (`Reset_Global_Arrays` / `Reset_Counters`); holds in-progress histograms before fold into `ldaTOT*`. |

## By domain

### `SystemState` (`g_system`)

| Fields | Kind | Notes |
|--------|------|-------|
| `bead_info`, `chain_info`, `topo_info`, `linker_len`, `naTotLattice` | Durable | Configuration / lattice occupancy |
| counts, box, mode flags, capability flags, linker floats | Durable | |
| `naLocalArr` | Durable (setup) | Fixed neighbor offsets filled once at init |
| `old_bead` | Scratch contents | Trial backup for reject/restore; allocated for the run |
| `naOld/NewOvlpNeighs`, `naOld/NewContNeighs`, `naTempR` | Scratch | Neighbor-search / trial coordinates |
| `naRot_IndArr` | Scratch contents | Rotation shuffle indices |

### `EnergyTables` (`g_energy`)

| Fields | Kind | Notes |
|--------|------|-------|
| `faEnergy`, `faEnRad`, `nBeadTypes` | Durable | Parsed tables |
| `faCurrEn` | Durable runtime | Current system energy vector |
| `naRotTrial`, `ldaBoltzFac*`, `ldaBoltzNorm` | Scratch | Rotational trial / Boltzmann temps |

### `McControl` (`g_mc`)

| Fields | Kind | Notes |
|--------|------|-------|
| temps, freqs, step counts, seed, `faKT_Cycle` | Durable control | From config / init |
| `fCuTemp`, `fSquishRad_Sq` | Durable runtime | |
| `naMCAccepMat` | Cycle buffer | Zeroed in `Reset_Global_Arrays` |

### `Reports` (`g_reports`)

| Fields | Kind | Notes |
|--------|------|-------|
| paths, freqs, `nTrajMode` | Durable | From config |
| `naTOTTRAJ_Arr`, `nTrajCurFrame` | Cycle buffer | Experimental traj accumulation |

### `AnalysisState` (`g_analysis`)

| Fields | Kind | Notes |
|--------|------|-------|
| `nClusteringMode`, size/comp metadata | Durable | |
| `naClusterMatrix`, `naList`, `laClusHistList`, `naChainCheckList` | Scratch / cycle | Cluster pass working sets |
| `ldaRDF_Arr`, `ldaRadDen_Arr`, `ldaMOLCLUS_Arr`, `faGyrTensor` | Cycle buffer | Current-cycle histograms |
| `ldaTOTRDF_Arr`, `ldaTOTRadDen_Arr`, `ldaTOTCLUS_Arr`, `ldaTOTMOLCLUS_Arr`, `ldaTOTRg_Arr` | Durable accumulators | Across cycles |
| `nTot*Counter`, `nLargestClusterRightNow` | Cycle counters | Cleared by `Reset_Counters` |

### `SimulationConfig` (`g_config`)

Durable **inputs** only (keyfile). Not runtime scratch.

## Reset hooks

- `Reset_Counters` — analysis/report cycle counters (scratch/cycle).
- `Reset_Global_Arrays` — analysis cycle buffers + MC accept matrix; then `Reset_Counters`.
