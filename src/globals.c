#include "global.h"

/*
 * Process-wide domain instances for Phase 2.1.
 * Zero-initialized like the former common-block globals; call sites still use
 * *_glb shims from global.h until later migration PRs.
 */
SystemState g_system;
EnergyTables g_energy;
McControl g_mc;
Reports g_reports;
AnalysisState g_analysis;
SimulationConfig g_config;
