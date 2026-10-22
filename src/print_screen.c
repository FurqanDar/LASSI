#include "print.h"
#include "initialize.h"
#include "cluster.h"
#include "energy.h"
#include "structure.h"


/// PrintToScreen_EnergyMatrix - fancy function to print symmetric matrices with labels
/// \param strTitle
/// \param nSeqEn
/// \param fArray
/// \param param
void PrintToScreen_EnergyMatrix(char* strTitle, int nSeqEn, float fArray[MAX_AA][MAX_AA][MAX_E], int param)
{
    const int nLen   = nSeqEn;
    const int outLen = (nLen + 3) * 5;
    int i, j;

    char sSectionHead[512];
    memset(sSectionHead, '-', outLen);
    sSectionHead[outLen] = '\0';

    printf("%s\n", sSectionHead);

    printf("| %-5s ", strTitle);
    for (i = 0; i < nLen; i++)
        {
            printf("%4d  ", i);
        }
    printf("%3s\n", "|");

    for (i = 0; i < nLen; i++)
        {
            printf("| %-5d ", i);
            for (j = 0; j < nLen; j++)
                {
                    printf("%5.2f ", fArray[i][j][param]);
                }
            printf("%3s\n", "|");
        }

    printf("%s\n", sSectionHead);
}


/// Print_LogToScreen - prints current status of the run. Overall move
/// acceptance ratios, and energies. \param nGen \param run_it
void Print_LogToScreen(long nGen, const int run_it)
{
    ScreenIO_Print_Log_FullRun(nGen, run_it);
}


/// PrintToScreen_AllEnergyMatrices
void PrintToScreen_AllEnergyMatrices(void)
{

    const char lBrace[] = "<======      ";
    const char rBrace[] = "      ======>";

    printf("%s Energy Matrices %s\n", lBrace, rBrace);

    PrintToScreen_EnergyMatrix("OVLP ", nBeadTypes_glb, faEnergy_glb, E_OVLP);

    PrintToScreen_EnergyMatrix("CONT ", nBeadTypes_glb, faEnergy_glb, E_CONT);

    PrintToScreen_EnergyMatrix("SC_SC", nBeadTypes_glb, faEnergy_glb, E_SC_SC);

    PrintToScreen_EnergyMatrix("FSOL ", nBeadTypes_glb, faEnergy_glb, E_F_SOL);

    PrintToScreen_EnergyMatrix("STIFF", nBeadTypes_glb, faEnergy_glb, E_STIFF);
}


/// PrintToScreen_MCMoveFreqs
void PrintToScreen_MCMoveFreqs(void)
{
    const char* MoveName[MAX_MV];
    MoveName[MV_PIVOT]      = "Pivot           ";
    MoveName[MV_DBPVT]      = "Double Pivot    ";
    MoveName[MV_CLSTR]      = "Larger Cluster  ";
    MoveName[MV_SMCLSTR]    = "Smaller Cluster ";
    MoveName[MV_STROT]      = "Face Change     ";
    MoveName[MV_LOCAL]      = "Local           ";
    MoveName[MV_COLOCAL]    = "Co-local        ";
    MoveName[MV_MTLOCAL]    = "Shake           ";
    MoveName[MV_BRROT]      = "Rotate Branched ";
    MoveName[MV_SNAKE]      = "Slithering Snake";
    MoveName[MV_TRANS]      = "Translation     ";
    MoveName[MV_PR_SMCLSTR] = "Pr. Smal Cluster";
    MoveName[MV_PR_CLSTR]   = "Pr. Cluster     ";

    float freqMin = 1e10f;
    int i;
    for (i = MV_NULL + 1; i < MAX_MV; i++)
        {
            if (freqMin >= g_mc.faMCFreq[i] && g_mc.faMCFreq[i] != 0.)
                {
                    freqMin = g_mc.faMCFreq[i];
                }
        }

    for (i = 0; i < 30; i++)
        {
            printf("-");
        }
    printf("\n");
    printf("| MC Move Frequencies%9s\n", "|");
    for (i = MV_NULL + 1; i < MAX_MV; i++)
        {
            printf("| %-16s %5d %5s\n", MoveName[i], (int) ceilf(g_mc.faMCFreq[i] / freqMin), "|");
        }
    for (i = 0; i < 30; i++)
        {
            printf("-");
        }
}


void ScreenIOUtil_PrintTemperatures(void)
{
    printf("Thermalizing Temperature       = %.2f\n", g_mc.fPreKT);

    if (g_system.nTemp_inv == 1)
        {
            puts("Temperature Inverted           = YES");
            printf("MC Temperatures: (First, Last) = (%.2f, %.2f)\n", g_mc.fKT,
                   1.f / (g_mc.fKT + (float) (g_mc.nTotCycleNum - 1) * g_mc.fDeltaTemp));
        }
    else
        {
            puts("Temperature Inverted           = NO");
            printf("MC Temperatures: (First, Last) = (%.2f, %.2f)\n", g_mc.fKT,
                   g_mc.fKT + (float) (g_mc.nTotCycleNum - 1) * g_mc.fDeltaTemp);
        }

    printf("Temperature Annealing Mode     = %d\n", g_system.nAnnealing_Mode);
}


/// ScreenIO_Print_KeyFile - print the keyfile that was read in to the screen
void ScreenIO_Print_KeyFile(void)
{ // should be output-dependent (stdout, stderr, other files)

    int i;
    const char lBrace[] = "<======      ";
    const char rBrace[] = "      ======>";
    printf("%s System Settings %s\n", lBrace, rBrace);
    printf("Number of Bead Types = %d\n", nBeadTypes_glb);
    printf("Number of Beads      = %ld\n", (size_t) g_system.tot_beads);
    printf("Number of Chains     = %ld\n", (size_t) g_system.tot_chains);
    printf("Number of Components = %ld\n", (size_t) g_system.tot_chain_types);
    printf("Lattice Dimensions   = %3d %3d %3d\n", g_system.naBoxSize[0], g_system.naBoxSize[1], g_system.naBoxSize[2]);
    printf("Monomer Density      = %1.2e\n",
           (float) g_system.tot_beads / (float) g_system.naBoxSize[0] / (float) g_system.naBoxSize[1] / (float) g_system.naBoxSize[2]);
    printf("\n");

    PrintToScreen_AllEnergyMatrices();

    printf("\n");

    printf("%s Linker Info %s\n", lBrace, rBrace);
    printf("Linker length             = %.2f\n", g_system.fLinkerLength);
    printf("Linker spring constant    = %.2f\n", g_system.fLinkerSprCon);
    printf("Linker equilibrium length = %.2f\n", g_system.fLinkerEqLen);
    printf("\n");

    printf("%s MC Setup %s\n", lBrace, rBrace);
    ScreenIOUtil_PrintTemperatures();

    printf("Indent Mode                    = %d\n", g_system.nBiasPotential_Mode);
    //    printf("Rotational Bias Mode           = %d\n", g_system.RotBias_Mode);
    printf("Number of MC Cycles            = %d\n", g_mc.nTotCycleNum);
    printf("Number of MC Steps/Cycle       = %e\n", (float) g_mc.nMCStepsPerCycle);
    printf("Number of Thermalizing Steps   = %e\n", (float) g_mc.nMCStepsForTherm);
    printf("RNG Seed                       = %d\n", g_mc.nRNG_Seed);
    printf("Clustering Mode                = %d\n", nClusteringMode_glb);

    PrintToScreen_MCMoveFreqs();
    printf("\n");
}


/// ScreenIO_Print_SystemEnergy
void ScreenIO_Print_SystemEnergy(void)
{

    int i;

    char sSectionHead[32];
    memset(sSectionHead, '-', 19);
    sSectionHead[19] = '\0';

    printf("%s\n", sSectionHead);
    puts("Energies          |");
    printf("Tot  : %-10.2e |\n", faCurrEn_glb[E_TOT]);
    printf("Ovlp : %-10.2e |\n", faCurrEn_glb[E_OVLP]);
    printf("Cont : %-10.2e |\n", faCurrEn_glb[E_CONT]);
    printf("Aniso: %-10.2e |\n", faCurrEn_glb[E_SC_SC]);
    printf("FSol : %-10.2e |\n", faCurrEn_glb[E_F_SOL]);
    printf("Stiff: %-10.2e |\n", faCurrEn_glb[E_STIFF]);
    printf("Bias : %-10.2e |\n", faCurrEn_glb[E_BIAS]);
    printf("%s\n", sSectionHead);
}


/// void ScreenIO_Print_AcceptanceRatios
void ScreenIO_Print_AcceptanceRatios(void)
{

    const char* MoveName[MAX_MV];
    MoveName[MV_PIVOT]      = "Pivot       ";
    MoveName[MV_DBPVT]      = "Double Pivot";
    MoveName[MV_CLSTR]      = "La Cluster  ";
    MoveName[MV_SMCLSTR]    = "Sm Cluster  ";
    MoveName[MV_STROT]      = "Face Change ";
    MoveName[MV_LOCAL]      = "Local       ";
    MoveName[MV_COLOCAL]    = "Co-local    ";
    MoveName[MV_MTLOCAL]    = "Multi Local ";
    MoveName[MV_BRROT]      = "Rot. Br.    ";
    MoveName[MV_SNAKE]      = "Sli. Snake  ";
    MoveName[MV_TRANS]      = "Translation ";
    MoveName[MV_PR_SMCLSTR] = "Pr. Sm. Cls.";
    MoveName[MV_PR_CLSTR]   = "Pr. Cls.    ";

    int i, j;
    float fAccRatio = 0.f;
    lLong nMoveSum  = 0L;

    char sSectionHead[32];
    memset(sSectionHead, '-', 22);
    sSectionHead[22] = '\0';

    printf("%s\n", sSectionHead);
    puts("Acceptance Ratios    |");
    for (i = 1; i < MAX_MV; i++)
        {
            nMoveSum = g_mc.naMCAccepMat[0][i] + g_mc.naMCAccepMat[1][i];
            if (nMoveSum)
                {
                    fAccRatio = 100.f * (float) g_mc.naMCAccepMat[1][i] / (float) nMoveSum;
                    printf("%-12s: %-7.2f|\n", MoveName[i], fAccRatio);
                }
            else
                {
                    printf("%-12s: %-7s|\n", MoveName[i], "NA");
                }
        }
    printf("%s\n", sSectionHead);
}


/// ScreenIO_Print_Log_Thermalization - print the log to the screen.
void ScreenIO_Print_Log_Thermalization(const long nGen)
{
    puts("****************************************");
    printf("Run Cycle: Thermalization\n");
    printf("Step     : %8.3e\n", (float) nGen);
    printf("MC Temp  : %8.3e\n", g_mc.fCuTemp);

    if (nTotClusCounter_glb > 0)
        {
            printf("Perc Phi : %8.3f\n",
                   ((float) nLargestClusterRightNow_glb) / ((float) nTotClusCounter_glb) / (float) g_system.tot_chains);
        }

    ScreenIO_Print_SystemEnergy();
    ScreenIO_Print_AcceptanceRatios();
    puts("****************************************");
}


/// ScreenIO_Print_Log_FullRun - print the log to the screen.
void ScreenIO_Print_Log_FullRun(const long nGen, const int run_cycle)
{
    puts("****************************************");
    printf("Run Cycle: %d\n", run_cycle);
    printf("Step     : %8.3e\n", (float) nGen);
    printf("MC Temp  : %8.3e\n", g_mc.fCuTemp);

    if (nTotClusCounter_glb > 0)
        {
            printf("Perc Phi : %8.3f\n",
                   ((float) nLargestClusterRightNow_glb) / ((float) nTotClusCounter_glb) / (float) g_system.tot_chains);
        }

    ScreenIO_Print_SystemEnergy();
    ScreenIO_Print_AcceptanceRatios();
    puts("****************************************");
}


/// ScreenIO_PrintStdErr_BeadPosition. Prints the beadID and it's position to STDERR.
/// \param beadID
void ScreenIO_PrintStdErr_BeadPosition(const int beadID)
{
    const int x = g_system.bead_info[beadID][POS_X];
    const int y = g_system.bead_info[beadID][POS_Y];
    const int z = g_system.bead_info[beadID][POS_Z];
    fprintf(stderr, "bead: %d is a position: <%d, %d, %d>\n", beadID, x, y, z);
}


/// ScreenIO_Print_SanityCheckFailurePreamble. Prints the preamble to the error message for sanity check failures.
/// \param nGen Which MC Step this crash has occurred on.
/// \param run_cycle Which run_cycle this crash has occurred on.
void ScreenIO_Print_SanityCheckFailurePreamble(const long nGen, const int run_cycle)
{
    fputs("\n-------------------------------------------------------------------------------", stderr);
    fputs("\nERROR! Sanity check failed! Crashing!\n", stderr);
    fprintf(stderr, "Crash occurred at (run_cycle: %d; mc_step: %ld)\n", run_cycle, nGen);
    fputs("Snapshot of system saved to crash_snapshot.txt\n", stderr);
    fputs("CRASH TYPE: ", stderr);
}


/// ScreenIO_Print_SanityFail_BeadPosAndLattPos. Prints the crash information for when the bead position does not
/// correspond to the lattice representation. Prints the bead position, and what the lattice thinks is at that position.
/// \param badBead The bead for which to print the information.
void ScreenIO_Print_SanityFail_BeadPosAndLattPos(const int badBead)
{

    const int x = g_system.bead_info[badBead][POS_X];
    const int y = g_system.bead_info[badBead][POS_Y];
    const int z = g_system.bead_info[badBead][POS_Z];
    fprintf(stderr, "bead: %d should be at: <%d, %d, %d>\n", badBead, x, y, z);

    const int latInd = Lat_Ind_OfBead(badBead);
    const int latVal = g_system.naTotLattice[latInd];

    fprintf(stderr, "while\nlattice: <%d, %d, %d> (ind: %d) has: %d\n", x, y, z, latInd, latVal);
}


/// ScreenIO_Print_SanityFail_MolecularStructure. Prints the crash information for when the linkers for badBead are
/// broken. Prints out the chain badBead is in, all the bonded partners for badBead, where the beads are at, and the
/// distances between the beads.
/// \param badBead The bead for which to print the information.
void ScreenIO_Print_SanityFail_MolecularStructure(const int badBead)
{

    const int badChain = g_system.bead_info[badBead][BEAD_CHAINID];

    fprintf(stderr, "chain: %d has broken structure\n", badChain);
    fprintf(stderr, "bead:  %d has broken linkers\n", badBead);
    ScreenIO_PrintStdErr_BeadPosition(badBead);
    fprintf(stderr, "\nbead:  %d has partners:\n", badBead);
    int idx      = 0;
    int bondPart = g_system.topo_info[badBead][idx];
    while (bondPart != -1)
        {
            ScreenIO_PrintStdErr_BeadPosition(bondPart);
            idx++;
            bondPart = g_system.topo_info[badBead][idx];
        }
    fputs("\n\nThe beads have the following distances and linker constraints:\n", stderr);

    idx             = 0;
    bondPart        = g_system.topo_info[badBead][idx];
    float bdDist   = 0.;
    float linkCons = 0.;

    while (bondPart != -1)
        {

            bdDist   = Dist_BeadToBead(badBead, bondPart);
            linkCons = (float) g_system.linker_len[badBead][idx] * LINKER_RSCALE;

            fprintf(stderr, "beads: (%d, %d) have dist: %5.3f with linker-constraint: %5.3f\n", badBead, bondPart,
                    bdDist, linkCons);
            idx++;
            bondPart = g_system.topo_info[badBead][idx];
        }
}


/// ScreenIO_Print_SanityFail_SelfBond. Prints the crash information for when the a bead is self-bonded. Prints out the
/// id of the bead.
/// \param badBead The bead for which to print the information.
void ScreenIO_Print_SanityFail_SelfBond(const int badBead)
{
    fprintf(stderr, "bead: %d is bonded to itself!\n", badBead);
}


/// ScreenIO_Print_SanityFail_BeadBondSymmetry. Prints the crash information for when a bond is not symmetric. Prints
/// out the id of the bead that failed, which bead badBead is bonded to, and which bead the beadPartner is bonded to.
/// We also print the positions of the beads.
/// \param badBead The bead for which to print the information.
void ScreenIO_Print_SanityFail_BeadBondSymmetry(const int badBead)
{
    const int beadPart = g_system.bead_info[badBead][BEAD_FACE];

    fprintf(stderr, "bead: %d is bonded to bead: %d\n", badBead, beadPart);

    const int partBond = g_system.bead_info[beadPart][BEAD_FACE];

    fprintf(stderr, "while\nbead: %d is bonded to bead: %d\n", beadPart, partBond);

    ScreenIO_PrintStdErr_BeadPosition(badBead);
    ScreenIO_PrintStdErr_BeadPosition(beadPart);
}


/// ScreenIO_Print_SanityFail_BeadBondDistance. Prints the crash information for when the distance between bonded beads
/// is too large. Prints the relevant beads' ids, positions and the distance between them.
/// \param badBead The bead for which to print the information.
void ScreenIO_Print_SanityFail_BeadBondDistance(const int badBead)
{
    const int beadPart = g_system.bead_info[badBead][BEAD_FACE];

    fprintf(stderr, "bead: %d is bonded to bead: %d\n", badBead, beadPart);

    const float bDist = Dist_BeadToBead(badBead, beadPart);
    fprintf(stderr, "while\ndistance: %5.3f\n", bDist);

    ScreenIO_PrintStdErr_BeadPosition(badBead);
    ScreenIO_PrintStdErr_BeadPosition(beadPart);
}

