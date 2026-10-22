#ifndef _GLOBAL_H_ // include guard
#define _GLOBAL_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_AA         10
#define MAX_CHAINTYPES 10
#define MAX_CHAINLEN   1000
#define MAX_BONDS      5
#define MAX_VALENCY    500

// energy parameters
#define E_TOT   0 // index zero must be assigned for total energy
#define E_OVLP  1
#define E_CONT  2
#define E_SC_SC 3
#define E_F_SOL 4
#define E_T_IND 5
#define E_STIFF 6
#define E_BIAS  7
#define MAX_E   8 // just for counting; must be the last number of this list

// MC move parameters
#define MV_NULL       0 // index zero indicates null move
#define MV_STROT      1
#define MV_LOCAL      2
#define MV_COLOCAL    3
#define MV_MTLOCAL    4
#define MV_SNAKE      5
#define MV_TRANS      6
#define MV_SMCLSTR    7
#define MV_CLSTR      8
#define MV_PIVOT      9
#define MV_BRROT      10
#define MV_DBPVT      11
#define MV_PR_SMCLSTR 12
#define MV_PR_CLSTR   13
#define MAX_MV        14 // just for counting; must be the last number of this list

// report parameters
#define REPORT_NULL    0 // index zero
#define REPORT_LOG     1
#define REPORT_ENERGY  2
#define REPORT_CONFIG  3
#define REPORT_MCMOVE  4
#define REPORT_NETWORK 5
#define REPORT_RDFTOT  6
#define REPORT_COMDEN  7
#define MAX_REPORT     8 // just for counting; must be the last number of this list

// auxiliary definitions
#define POS_X   0
#define POS_Y   1
#define POS_Z   2
#define POS_MAX 3

#define BEAD_CHAINID 3
#define BEAD_TYPE    4
#define BEAD_FACE    5
#define BEADINFO_MAX 6

#define CHAIN_TYPE    0
#define CHAIN_LENGTH  1
#define CHAIN_START   2
#define CHAININFO_MAX 3

#define MAX_ROTSTATES       27
#define CLUS_CONTACT_NEIGHS 27

typedef int lInt;
typedef long lLong;
typedef double lDub;
typedef long double lLDub;

#define LARGEST_RADIUS 3
#define LINKER_RSCALE  1.74f

/*
 * Phase 2.1: domain structs with a single process-wide instance each.
 * Compatibility shims (#define …_glb) keep existing call sites working;
 * later 2.2+ PRs migrate call sites to struct fields / accessors.
 */

/// Lattice configuration, topology flags, neighbor scratch, linker params.
typedef struct SystemState
{
    // configurations and structural info
    lInt** bead_info;
    lInt** old_bead;
    lInt** chain_info;
    lInt** topo_info;
    lInt** linker_len;

    int tot_beads;
    int tot_chains;
    int tot_chain_types;
    lInt nAnnealing_Mode, nTemp_inv;
    lInt nBiasPotential_Mode, RotBias_Mode;
    lInt nBiasPotential_KeepON, nBiasPotential_CoupledToTemp;

    // system setup
    lInt naBoxSize[POS_MAX];
    lInt naLocalArr[MAX_ROTSTATES - 1][3];
    lInt naRot_IndArr[MAX_ROTSTATES - 1];
    char bReadConf;
    char bSystemHasTopo;
    char bSystemHasCont;
    char bSystemHasOvlp;
    char bSystemHasFSol;
    char bSystemHasSCSC;

    // topology / interaction capability flags
    lInt nBeadTypeIsSticker[MAX_AA];
    lInt nChainTypeIsLinear[MAX_CHAINTYPES];
    lInt nBeadTypeCanOvlp[MAX_AA];
    lInt nBeadTypeCanCont[MAX_AA];
    lInt nBeadTypeCanFSol[MAX_AA];
    lInt nBeadTypeCanTInd[MAX_AA];

    float fLinkerLength;
    float fLinkerSprCon;
    float fLinkerEqLen;

    // neighbor-search scratch
    lInt* naOldOvlpNeighs;
    lInt* naNewOvlpNeighs;
    lInt* naOldContNeighs;
    lInt* naNewContNeighs;

    lInt naTempR[POS_MAX];

    // lattice
    lInt* naTotLattice;
} SystemState;

/// Energy matrices, Boltzmann factors, current energy vector.
typedef struct EnergyTables
{
    lInt nBeadTypes;
    float faEnergy[MAX_AA][MAX_AA][MAX_E];
    float faEnRad[MAX_AA][MAX_AA][MAX_E];
    lInt naRotTrial[MAX_VALENCY][MAX_ROTSTATES];
    lLDub ldaBoltzFac[MAX_ROTSTATES - 1];
    lLDub ldaBoltzNorm[MAX_VALENCY];
    lLDub ldaBoltzFacNorm[MAX_AA][MAX_AA];
    lLDub ldLogOfSmallestPossibleProb;
    float faCurrEn[MAX_E];
} EnergyTables;

/// MC temperatures, frequencies, step counts, acceptance matrix, RNG seed.
typedef struct McControl
{
    float fKT, fPreKT, fCuTemp, fDeltaTemp, fMC_TempRate, fSquishRad, fSquishRad_Sq, fSquish_Stiff;
    float* faKT_Cycle;
    lLong nMCStepsPerCycle, nMCStepsForTherm;
    float faMCFreq[MAX_MV];
    lInt nTotCycleNum;
    lInt nRNG_Seed;
    lLong naMCAccepMat[2][MAX_MV];
} McControl;

/// Report filenames / prefixes and trajectory buffers.
typedef struct Reports
{
    char strReportPrefix[400];
    char strFileEnergy[512];
    char strFileTraj[512];
    char strFileMCMove[512];
    char strFileSysProp[512];
    char strRestartFile[512];
    lLong naReportFreqs[MAX_REPORT];
    lLong nTrajMode;
    lInt* naTOTTRAJ_Arr;
    lLong nTraj_FramesPerCycle;
    lInt nTrajCurFrame;
} Reports;

/// Cluster, RDF / radial density, and gyration analysis state.
typedef struct AnalysisState
{
    // cluster analysis
    lInt nClusteringMode;
    lInt** naClusterMatrix;
    lInt* naList;
    lLong* laClusHistList;
    lInt* naChainCheckList;
    lInt nTotClusCounter;
    lLDub** ldaTOTCLUS_Arr;
    lLDub* ldaMOLCLUS_Arr;
    lLDub* ldaTOTMOLCLUS_Arr;
    lInt nLargestClusterRightNow;
    lInt nLimitedClusterSize;

    // radial densities and PDFs
    lLDub* ldaTOTRDF_Arr;
    lLDub* ldaTOTRadDen_Arr;
    lLDub* ldaRDF_Arr;
    lLDub* ldaRadDen_Arr;
    lInt nRDF_TotComps;
    lInt nTotRDFCounter;
    lInt nRDF_TotBins;
    lInt nRadDen_TotComps;
    lInt nRadDen_CompShift;
    lInt nTotRadDenCounter;

    // gyration tensor
    float faGyrTensor[7];
    float faSysGyrRad;
    lLDub** ldaTOTRg_Arr;
    lInt nTotGyrRadCounter;
} AnalysisState;

/* Single process-wide instances (defined in globals.c). */
extern SystemState g_system;
extern EnergyTables g_energy;
extern McControl g_mc;
extern Reports g_reports;
extern AnalysisState g_analysis;

/*
 * Compatibility shims: existing *_glb names alias into the domain structs.
 * Do not take the address of a shim macro expecting a distinct global object
 * beyond the underlying struct field (same storage as before, now nested).
 */
/* SystemState *_glb shims removed after Phase 2.2 call-site migration to g_system.* */

/* EnergyTables *_glb shims removed after Phase 2 call-site migration to g_energy.* */

/* McControl *_glb shims removed after Phase 2.2 call-site migration to g_mc.* */

/* Reports *_glb shims removed after Phase 2.5 call-site migration to g_reports.* */

#define nClusteringMode_glb                g_analysis.nClusteringMode
#define naClusterMatrix_glb                g_analysis.naClusterMatrix
#define naList_glb                         g_analysis.naList
#define laClusHistList_glb                 g_analysis.laClusHistList
#define naChainCheckList_glb               g_analysis.naChainCheckList
#define nTotClusCounter_glb                g_analysis.nTotClusCounter
#define ldaTOTCLUS_Arr_glb                 g_analysis.ldaTOTCLUS_Arr
#define ldaMOLCLUS_Arr_glb                 g_analysis.ldaMOLCLUS_Arr
#define ldaTOTMOLCLUS_Arr_glb              g_analysis.ldaTOTMOLCLUS_Arr
#define nLargestClusterRightNow_glb        g_analysis.nLargestClusterRightNow
#define nLimitedClusterSize_glb            g_analysis.nLimitedClusterSize
#define ldaTOTRDF_Arr_glb                  g_analysis.ldaTOTRDF_Arr
#define ldaTOTRadDen_Arr_glb               g_analysis.ldaTOTRadDen_Arr
#define ldaRDF_Arr_glb                     g_analysis.ldaRDF_Arr
#define ldaRadDen_Arr_glb                  g_analysis.ldaRadDen_Arr
#define nRDF_TotComps_glb                  g_analysis.nRDF_TotComps
#define nTotRDFCounter_glb                 g_analysis.nTotRDFCounter
#define nRDF_TotBins_glb                   g_analysis.nRDF_TotBins
#define nRadDen_TotComps_glb               g_analysis.nRadDen_TotComps
#define nRadDen_CompShift_glb              g_analysis.nRadDen_CompShift
#define nTotRadDenCounter_glb              g_analysis.nTotRadDenCounter
#define faGyrTensor_glb                    g_analysis.faGyrTensor
#define faSysGyrRad_glb                    g_analysis.faSysGyrRad
#define ldaTOTRg_Arr_glb                   g_analysis.ldaTOTRg_Arr
#define nTotGyrRadCounter_glb              g_analysis.nTotGyrRadCounter

#endif // _GLOBAL_H_
