#include "structure.h"

/* print.h is I/O-layer; forward-declare only the crash/sanity printers we call
 * so structure does not take a hard include edge upward into I/O (Phase 3.1). */
void ScreenIO_Print_SanityCheckFailurePreamble(long nGen, int run_cycle);
void ScreenIO_Print_SanityFail_BeadPosAndLattPos(int badBead);
void ScreenIO_Print_SanityFail_MolecularStructure(int badBead);
void ScreenIO_Print_SanityFail_SelfBond(int badBead);
void ScreenIO_Print_SanityFail_BeadBondSymmetry(int badBead);
void ScreenIO_Print_SanityFail_BeadBondDistance(int badBead);
void FileIO_PrintCrashSnapshot(void);

/// CheckSystemUtil_BeadPosAndLattPosOK. Loop over all beads and make sure that the bead's internal position
/// matches the beadID that is placed at that location.
/// \return beadID of the bead that failed, -1 if all good.
int CheckSystemUtil_BeadPosAndLattPosOK(void)
{
    int beadID;
    int latt_pos_ind;
    int latt_bead;

    for (beadID = 0; beadID < g_system.tot_beads; beadID++)
        {
            latt_pos_ind = Lat_Ind_OfBead(beadID);
            latt_bead    = g_system.naTotLattice[latt_pos_ind];
            if (latt_bead != beadID)
                {
                    return beadID;
                }
        }

    return -1;
}


/// Check_LinkerConstraintForBead - we loop over all possible covalently bonded beads for beadID and check if the
/// linkers are OK.
/// \param beadID
/// \return beadID of the bead-bond-partner that failed, -1 if all good.
int Check_LinkerConstraintForBead(const int beadID)
{
    int idx; // Iterator to loop over bond Partners
    int bondPartner;
    idx         = 0;
    bondPartner = g_system.topo_info[beadID][idx]; // Initializing the two.
    while (idx < MAX_BONDS && g_system.topo_info[beadID][idx] != -1)
        { // Keep going till we run out of partners
            bondPartner = g_system.topo_info[beadID][idx];
            if (Dist_BeadToBead(beadID, bondPartner) > LINKER_RSCALE * (float) g_system.linker_len[beadID][idx])
                {
                    return bondPartner; // This means that we have broken one of the linkers.
                }
            idx++;
        }
    return -1; // This means that all linker constraints are satisfied.
}


/// CheckSystemUtil_MolecularStructuresOK. Loop over all beads and check if the linker constraints of all the
/// beads are satisfied correctly.
/// \return
int CheckSystemUtil_MolecularStructuresOK(void)
{
    int beadID;

    for (beadID = 0; beadID < g_system.tot_chains; beadID++)
        {
            if (Check_LinkerConstraintForBead(beadID) >= 0)
                {
                    return beadID;
                }
        }

    return -1;
}


/// CheckSystemUtil_BeadBondsSymmetricOK. Loops over all beads and makes sure that the bead I am bonded to is also
/// bonded to me. \return beadID of the bead that failed, -1 if all good.
int CheckSystemUtil_BeadBondsSymmetricOK(void)
{
    int beadID;

    for (beadID = 0; beadID < g_system.tot_beads; beadID++)
        {
            int bondPartner = g_system.bead_info[beadID][BEAD_FACE];
            if (bondPartner != -1)
                {
                    if (g_system.bead_info[bondPartner][BEAD_FACE] != beadID)
                        {
                            return beadID;
                        }
                }
        }

    return -1;
}


/// CheckSystemUtil_BeadBondsDistanceOK. Loops over all beads and checks that the distance between bonded beads is not
/// greater than the allowed threshold.
/// \return beadID of the bead that failed, -1 if all good.
int CheckSystemUtil_BeadBondsDistanceOK(void)
{
    int beadID;
    //    float bDist;
    for (beadID = 0; beadID < g_system.tot_beads; beadID++)
        {
            int bondPartner = g_system.bead_info[beadID][BEAD_FACE];
            if (bondPartner != -1)
                {
                    const float bDist = Dist_BeadToBead(beadID, bondPartner);
                    if (bDist > LINKER_RSCALE)
                        {
                            return beadID;
                        }
                }
        }

    return -1;
}


/// CheckSystemUtil_NoSelfBonds. Loops over all beads and makes sure that no bead is bonded to itself.
/// \return beadID of the bead that failed, -1 if all good.
int CheckSystemUtil_NoSelfBonds(void)
{
    int beadID;

    for (beadID = 0; beadID < g_system.tot_beads; beadID++)
        {
            int bondPartner = g_system.bead_info[beadID][BEAD_FACE];
            if (bondPartner == beadID)
                {
                    return beadID;
                }
        }

    return -1;
}


/// PerformRuntimeSanityCheck_BeadPosAndLattPos. Performs sanity check to make sure that the lattice has the bead
/// at the bead's location. If there is an error, we exit out of the program after printing out the relevant information
/// \param nGen Which MC Step this crash has occurred on.
/// \param run_cycle Which run_cycle this crash has occurred on.
void PerformRuntimeSanityCheck_BeadPosAndLattPos(const long nGen, const int run_cycle)
{
    int badBead = CheckSystemUtil_BeadPosAndLattPosOK();
    if (badBead != -1)
        {
            ScreenIO_Print_SanityCheckFailurePreamble(nGen, run_cycle);
            FileIO_PrintCrashSnapshot();
            fputs("Lattice positions and bead positions do not match!\n\nDetails of crash:\n", stderr);

            ScreenIO_Print_SanityFail_BeadPosAndLattPos(badBead);

            fputs("-------------------------------------------------------------------------------", stderr);
            exit(1);
        }
}


/// PerformRuntimeSanityCheck_MolecularStructure. Performs a sanity check to make sure all linkers are unbroken.
/// If there is an error, we exit out of the program after printing out the relevant information
/// \param nGen Which MC Step this crash has occurred on.
/// \param run_cycle Which run_cycle this crash has occurred on.
void PerformRuntimeSanityCheck_MolecularStructure(const long nGen, const int run_cycle)
{
    const int badBead = CheckSystemUtil_MolecularStructuresOK();
    if (badBead != -1)
        {
            ScreenIO_Print_SanityCheckFailurePreamble(nGen, run_cycle);
            FileIO_PrintCrashSnapshot();
            fputs("Molecular structure has been broken!\n\nDetails of crash:\n", stderr);

            ScreenIO_Print_SanityFail_MolecularStructure(badBead);

            fputs("-------------------------------------------------------------------------------", stderr);
            exit(1);
        }
}


/// PerformRuntimeSanityCheck_SelfBonds. Performs a sanity check to make sure no bead is self-bonded.
/// If there is an error, we print out the relevant information and crash the simulation.
/// \param nGen Which MC Step this crash has occurred on.
/// \param run_cycle Which run_cycle this crash has occurred on.
void PerformRuntimeSanityCheck_SelfBonds(const long nGen, const int run_cycle)
{
    const int badBead = CheckSystemUtil_NoSelfBonds();
    if (badBead != -1)
        {
            ScreenIO_Print_SanityCheckFailurePreamble(nGen, run_cycle);
            FileIO_PrintCrashSnapshot();
            fputs("A bead is self-bonded!\n\nDetails of crash:\n", stderr);

            ScreenIO_Print_SanityFail_SelfBond(badBead);

            fputs("-------------------------------------------------------------------------------", stderr);
            exit(1);
        }
}


/// PerformRuntimeSanityCheck_BondSymmetry. Performs a sanity check to make sure all bonded beads are bonded to each
/// other, or that the bonds are symmetric. If there is an error, we print out the relevant information and crash
/// the simulation.
/// \param nGen Which MC Step this crash has occurred on.
/// \param run_cycle Which run_cycle this crash has occurred on.
void PerformRuntimeSanityCheck_BondSymmetry(const long nGen, const int run_cycle)
{
    const int badBead = CheckSystemUtil_BeadBondsSymmetricOK();
    if (badBead != -1)
        {
            ScreenIO_Print_SanityCheckFailurePreamble(nGen, run_cycle);
            FileIO_PrintCrashSnapshot();
            fputs("Anisotropic bonds are not symmetric!\n\nDetails of crash:\n", stderr);

            ScreenIO_Print_SanityFail_BeadBondSymmetry(badBead);

            fputs("-------------------------------------------------------------------------------", stderr);
            exit(1);
        }
}


/// PerformRuntimeSanityCheck_BondDistance. Performs a sanity check to make sure that the distance between all bonded
/// beads is not greater than the global allowed threshold.
/// \param nGen Which MC Step this crash has occurred on.
/// \param run_cycle Which run_cycle this crash has occurred on.
void PerformRuntimeSanityCheck_BondDistance(const long nGen, const int run_cycle)
{
    const int badBead = CheckSystemUtil_BeadBondsDistanceOK();
    if (badBead != -1)
        {
            ScreenIO_Print_SanityCheckFailurePreamble(nGen, run_cycle);
            FileIO_PrintCrashSnapshot();
            fputs("Bond distance is too large!\n\nDetails of crash:\n", stderr);

            ScreenIO_Print_SanityFail_BeadBondDistance(badBead);

            fputs("-------------------------------------------------------------------------------", stderr);
            exit(1);
        }
}


/// PerformRuntimeSanityChecks. Performs the following sanity checks:
/// 1: PerformRuntimeSanityCheck_BeadPosAndLattPos
/// 2: PerformRuntimeSanityCheck_MolecularStructure
/// 3: PerformRuntimeSanityCheck_SelfBonds
/// 4: PerformRuntimeSanityCheck_BondSymmetry
/// 5: PerformRuntimeSanityCheck_BondDistance
/// \param nGen
/// \param run_cycle
void PerformRuntimeSanityChecks(const long nGen, const int run_cycle)
{
    /*
     * The following set of snippets cause specific failures.
     */
    //    //Causes lattice failure only
    //    g_system.naTotLattice[Lat_Ind_OfBead(10)]=5;
    //
    //    //Causes structure failure only
    //    g_system.bead_info[2][0] = 30;
    //    g_system.bead_info[3][0] = 0;
    //    g_system.naTotLattice[Lat_Ind_OfBead(2)]=2;
    //    g_system.naTotLattice[Lat_Ind_OfBead(3)]=3;
    //
    //    //Causes bond-symmetry failure only
    //    g_system.bead_info[4][BEAD_FACE] = 5;
    //    g_system.bead_info[5][BEAD_FACE] = -1;
    //
    //    //Causes self-bond failure only
    //    g_system.bead_info[6][BEAD_FACE] = 6;
    //
    //    //Causes bond-distance failure only
    //    g_system.bead_info[2][BEAD_FACE] = 3;
    //    g_system.bead_info[3][BEAD_FACE] = 2;

    PerformRuntimeSanityCheck_BeadPosAndLattPos(nGen, run_cycle);

    PerformRuntimeSanityCheck_MolecularStructure(nGen, run_cycle);

    PerformRuntimeSanityCheck_SelfBonds(nGen, run_cycle);

    PerformRuntimeSanityCheck_BondSymmetry(nGen, run_cycle);

    PerformRuntimeSanityCheck_BondDistance(nGen, run_cycle);
}


/// Check_LinkerConstraint - if I move beadID to tmpR, do I still satisfy the linker lengths for beadID?
/// For the proposed location, we loop over all bonded partners of beadID and check if the distance is within the
/// g_system.linker_len distance.
/// \param beadID
/// \param tmpR
/// \return 1 means all is good, 0 means bad.
int Check_LinkerConstraint(const int beadID, const int* const tmpR)
{
    // Check if the proposed new location for beadID is such that all the linkers are unbroken.
    int idx;         // Iterator to loop over bond Partners
    int bondPartner; // It is what it is.
    idx         = 0;
    bondPartner = g_system.topo_info[beadID][idx]; // Initializing the two.
    while (idx < MAX_BONDS && g_system.topo_info[beadID][idx] != -1)
        { // Keep going till we run out of partners
            bondPartner = g_system.topo_info[beadID][idx];
            if (Dist_PointToPoint(g_system.bead_info[bondPartner], tmpR) >
                LINKER_RSCALE * (float) g_system.linker_len[beadID][idx])
                {
                    return 0; // This means that we have broken one of the linkers.
                }
            idx++;
        }
    return 1; // This means that all linker constraints are satisfied.
}


/// Check_LinkerConstraints_ForBeadList - given that all the beads in beadList are at the new
/// positions already, as part of the MTLocal move, we iterate over every bead and check if
/// any linker constraints are broken. If so, return 0. If not, return 1 for success.
/// \param beadID
/// \param tmpR
/// \return
int Check_LinkerConstraints_ForBeadList(const int listSize, const int* beadList)
{

    int i, j, bead1, bead2;
    int dumBonds[MAX_BONDS + 1];
    int dumNum;
    float xDis;
    for (i = 0; i < listSize; i++)
        {
            bead1  = beadList[i];
            dumNum = OP_GetTopoBonds(bead1, dumBonds);
            for (j = 0; j < dumNum - 1; j++)
                { 
                    bead2 = dumBonds[j+1];
                    xDis  = Dist_BeadToBead(bead1, bead2);
                    if (xDis > LINKER_RSCALE * (float) g_system.linker_len[bead1][j])
                        {
                            return 0;
                        }
                }
        }

    return 1;
}


/// Check_BeadID_InList: Iterate over the list of beads to see if the bead is in the list.
/// If so, return 1. Otherwise return 0.
/// \param thisBeadID
/// \param listSize
/// \param beadList
/// \return
int Check_BeadID_InList(const int thisBeadID, const int listSize, const int beadList[MAX_BONDS + 1])
{
    int i;
    for (i = 0; i < listSize; i++)
        {
            if (thisBeadID == beadList[i])
                {
                    return 1;
                }
        }
    return 0;
}

