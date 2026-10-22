#include "mcmove.h"
#include "cluster.h"
#include "energy.h"
#include "structure.h"

/// Check_ChainDisp - checks if moving chainID by vec_disp[POS_MAX] will lead to
/// sterix clash.
/// \param chainID
/// \param vec_disp
/// \return 0 if clash, 1 if no clash.
int Check_ChainDisp(const int chainID, const int* vec_disp)
{ // Checks if chain can be displaced by vec_disp
    int i;
    int r_chck[POS_MAX];
    const int firstB = chain_info_glb[chainID][CHAIN_START];
    const int lastB  = firstB + chain_info_glb[chainID][CHAIN_LENGTH];
    for (i = firstB; i < lastB; i++)
        {
            LatPos_add_wPBC(r_chck, vec_disp, bead_info_glb[i]);
            if (naTotLattice_glb[Lat_Ind_FromVec(r_chck)] != -1)
                { // Steric clash
                    return 0;
                }
        }

    return 1;
}

/// OP_System_DispChain - displaced chainID by movR[POS_MAX], while remembering the
/// chain in old_beads. Also handles the lattice placement.
/// \param chainID
/// \param movR
void OP_System_DispChain(int chainID, const int* movR)
{
    int i, l;
    const int fB = chain_info_glb[chainID][CHAIN_START];
    const int lB = fB + chain_info_glb[chainID][CHAIN_LENGTH];
    int r_pos0[POS_MAX], r_posNew[POS_MAX];

    // Copy beads to old.
    OP_CopyBeadsToOld(fB, lB);

    for (i = fB; i < lB; i++)
        {
            LatPos_copy(r_pos0, bead_info_glb[i]);
            LatPos_add_wPBC(bead_info_glb[i], r_pos0, movR);
            LatPos_copy(r_posNew, bead_info_glb[i]);

            naTotLattice_glb[Lat_Ind_FromVec(r_pos0)]   = -1; // Removing from old place
            naTotLattice_glb[Lat_Ind_FromVec(r_posNew)] = i;  // Adding to new place
        }
}

/// OP_System_DispChain_ForTrans - displaced chainID by movR[POS_MAX], while
/// remembering the chain in old_beads. Also handles the lattice placement.
/// Move_Trans variant where I break all the physical bonds.
/// \param chainID
/// \param movR
void OP_System_DispChain_ForTrans(const int chainID, const int* movR)
{
    // Displaces current chain by movR and handles lattice
    // Specific for Move_Trans because it breaks old bonds!
    // Also remembers where everything was moved and saves into old_bead_glb

    int i, l;
    const int fB = chain_info_glb[chainID][CHAIN_START];
    const int lB = fB + chain_info_glb[chainID][CHAIN_LENGTH];
    int r_pos0[POS_MAX], r_posNew[POS_MAX];

    // Copy beads to old.
    OP_CopyBeadsToOld(fB, lB);

    for (i = fB; i < lB; i++)
        {
            LatPos_copy(r_pos0, bead_info_glb[i]);
            LatPos_add_wPBC(bead_info_glb[i], r_pos0, movR);
            LatPos_copy(r_posNew, bead_info_glb[i]);

            naTotLattice_glb[Lat_Ind_FromVec(r_pos0)]   = -1; // Removing from old place
            naTotLattice_glb[Lat_Ind_FromVec(r_posNew)] = i;  // Adding to new place
            OP_Beads_BreakBond(i);
        }
}

/// OP_System_RestoreChain - uses old_bead_glb to undo what OP_System_DispChain does.
/// \param chainID
void OP_System_RestoreChain(int chainID)
{ // Uses old_bead_glb to undo what OP_System_DispChain does.
    int i, l;
    const int fB = chain_info_glb[chainID][CHAIN_START];
    const int lB = fB + chain_info_glb[chainID][CHAIN_LENGTH];

    // Removing from 'new' place.
    for (i = fB; i < lB; ++i)
        {
            naTotLattice_glb[Lat_Ind_OfBead(i)] = -1;
        }

    // Remembering where we were.
    for (i = fB; i < lB; ++i)
        {
            OP_CopyBead(bead_info_glb[i], old_bead_glb[i]);
        }

    // Placing back.
    for (i = fB; i < lB; ++i)
        {
            naTotLattice_glb[Lat_Ind_OfBead(i)] = i;
        }
}

/// OP_CopyBead: Copies BEADINFO_MAX elements of orig_arr into copy_arr
/// \param copy_arr
/// \param orig_arr
void OP_CopyBead(int* restrict copy_arr, const int* restrict orig_arr)
{
    int j;
    for (j = 0; j < BEADINFO_MAX; j++)
        {
            copy_arr[j] = orig_arr[j];
        }
}

/// OP_CopyBeadsToOld: Copies all beads in [firstB, lastB) to old_beads. This is
/// sort of a wrapper for ease.
/// \param firstB
/// \param lastB
void OP_CopyBeadsToOld(const int firstB, const int lastB)
{
    int i;
    for (i = firstB; i < lastB; i++)
        {
            OP_CopyBead(old_bead_glb[i], bead_info_glb[i]);
        }
}

/// OP_RestoreBeadsFromOld: Copies all beads in [firstB, lastB) from old_beads.
/// This is sort of a wrapper for ease.
/// \param firstB
/// \param lastB
void OP_RestoreBeadsFromOld(const int firstB, const int lastB)
{
    int i, j;
    for (i = firstB; i < lastB; i++)
        {
            OP_CopyBead(bead_info_glb[i], old_bead_glb[i]);
        }
}

/// OP_System_RestoreChain_ForTrans - translation variant to restore chainID.
/// \param chainID
void OP_System_RestoreChain_ForTrans(int chainID)
{ // Uses old_bead_glb to undo what OP_System_DispChain_ForTrans does.
    // Note that removing and placing separately makes sure that no incorrect
    // bonds are formed, or unformed.
    int i, l;
    const int fB = chain_info_glb[chainID][CHAIN_START];
    const int lB = fB + chain_info_glb[chainID][CHAIN_LENGTH];

    // Breaking all current bonds
    for (i = fB; i < lB; ++i)
        {
            OP_Beads_BreakBond(i);
        }

    // Removing from 'new' place.
    for (i = fB; i < lB; ++i)
        {
            naTotLattice_glb[Lat_Ind_OfBead(i)] = -1;
        }

    // Remembering where we were.
    for (i = fB; i < lB; ++i)
        {
            OP_CopyBead(bead_info_glb[i], old_bead_glb[i]);
        }

    // Placing back.
    for (i = fB; i < lB; ++i)
        {
            naTotLattice_glb[Lat_Ind_OfBead(i)] = i;
        }

    // Restoring older bonds.
    for (i = fB; i < lB; ++i)
        {
            OP_Beads_RestoreBond(i);
        }
}

/// OP_RestoreChain_ForSnake - restores the chain: beadID's betweem fB and lB-1,
/// to before the snake move.
/// \param fB
/// \param lB
void OP_RestoreChain_ForSnake(const int fB, const int lB)
{
    int i, j;

    for (i = fB; i < lB; i++)
        { // Resetting the lattice
            if (bead_info_glb[i][BEAD_FACE] != -1)
                { // Need to break the newly proposed bond
                    bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_FACE] = -1;
                }
            naTotLattice_glb[Lat_Ind_OfBead(i)] = -1;
        }
    for (i = fB; i < lB; i++)
        {
            for (j = 0; j < BEADINFO_MAX; j++)
                {
                    bead_info_glb[i][j] = old_bead_glb[i][j]; // Restoring
                }
            if (bead_info_glb[i][BEAD_FACE] != -1)
                { // I was bonded so restore the bond
                    bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_FACE] = i;
                }
            naTotLattice_glb[Lat_Ind_OfBead(i)] = i;
        }
}

/// Check_MoveBeadTo - checks if newPos[POS_MAX] is empty on the lattice.
/// \param newPos
/// \return
int Check_MoveBeadTo(const int* newPos)
{ // Checks if I can move here

    if (naTotLattice_glb[Lat_Ind_FromVec(newPos)] != -1)
        {
            return 0; // There's something here already, brah
        }
    return 1;
}

/// OP_System_MoveBeadTo - move beadID to newPos[POS_MAX], while remembering the
/// bead properties, and handling the lattice.
/// \param beadID
/// \param newPos
void OP_System_MoveBeadTo(const int beadID, const int* newPos)
{ // Updates position to new one and handles lattice
    OP_CopyBead(old_bead_glb[beadID], bead_info_glb[beadID]);
    LatPos_copy(bead_info_glb[beadID], newPos);
    naTotLattice_glb[Lat_Ind_FromVec(old_bead_glb[beadID])] = -1;     // Removing from old place.
    naTotLattice_glb[Lat_Ind_FromVec(newPos)]               = beadID; // Placing in new place.
}

/// OP_System_MoveBeadTo_Inv - performs the inverse of OP_System_MoveBeadTo to
/// restore beadID using old_bead_glb.
/// \param beadID
void OP_System_MoveBeadTo_Inv(int beadID)
{                                                  // Undoes what OP_System_MoveBeadTo does
    naTotLattice_glb[Lat_Ind_OfBead(beadID)] = -1; // Removing from newly proposed place.
    OP_CopyBead(bead_info_glb[beadID], old_bead_glb[beadID]);
    naTotLattice_glb[Lat_Ind_OfBead(beadID)] = beadID; // Placing back where we used to be.
    OP_Beads_RestoreBond(beadID);
}

/// OP_MoveBeadTo_ForMTLocal - moves the bead over to newPos[POS_MAX], handles
/// the lattice and breaks bonds. This function is specifically for the
/// MultiLocal move.
/// \param beadID
/// \param newPos
void OP_MoveBeadTo_ForMTLocal(int beadID, const int* newPos)
{ // Updates position to new one and handles
  // lattice specifically for shake move
    int i;
    int tmpR2[POS_MAX];
    for (i = 0; i < POS_MAX; i++)
        {
            bead_info_glb[beadID][i] = newPos[i];
            tmpR2[i]                 = bead_info_glb[beadID][i];
        }
    naTotLattice_glb[Lat_Ind_FromVec(tmpR2)] = beadID;
    i                                        = bead_info_glb[beadID][BEAD_FACE];
    if (i != -1)
        {
            bead_info_glb[i][BEAD_FACE]      = -1;
            bead_info_glb[beadID][BEAD_FACE] = -1;
        }
}

void OP_System_MoveBeadsInListToPos(const int listSize, const int* beadList, const int (*newPos)[POS_MAX])
{

    OP_Lattice_EmptySitesForListOfBeads(listSize, beadList);
    OP_Beads_MoveBeadsInListToPos(listSize, beadList, newPos);
    OP_Lattice_PlaceBeadsInList(listSize, beadList);
}

void OP_Inv_MoveBeadsTo_InList(const int listSize, const int beadList[MAX_BONDS + 1],
                               const int newPos[MAX_BONDS + 1][POS_MAX])
{
    OP_Lattice_EmptySitesForListOfBeads(listSize, beadList);
    //    OP_Lattice_PlaceBeadsInList()
}

/// OP_SwapBeads - swaps all the properties, including bonding partners between
/// bead1 and bead2. Used in the DbPvt move.
/// \param bead1
/// \param bead2
void OP_SwapBeads(int bead1, int bead2)
{
    // Swaps ALL properties of the beads, and exchanged bonding partners
    int MyF1, MyF2;
    int i;
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    // First the coordinates
    for (i = 0; i < POS_MAX; i++)
        {
            tmpR[i]                 = bead_info_glb[bead1][i];
            tmpR2[i]                = bead_info_glb[bead2][i];
            bead_info_glb[bead1][i] = tmpR2[i];
            bead_info_glb[bead2][i] = tmpR[i];
        }

    // Onto bonds and partners
    MyF1 = bead_info_glb[bead1][BEAD_FACE]; // Bead 1's bonding partner
    MyF2 = bead_info_glb[bead2][BEAD_FACE]; // Bead 2's bonding partner.
    // If they are bonded to each other, don't do anything
    if (MyF1 != bead2)
        {
            bead_info_glb[bead1][BEAD_FACE] = MyF2;
            bead_info_glb[bead2][BEAD_FACE] = MyF1;
            if (MyF1 != -1)
                {                                           // Need to swap partners -- very 2019
                    bead_info_glb[MyF1][BEAD_FACE] = bead2; // It's bonded to bead2 now
                }
            if (MyF2 != -1)
                {
                    bead_info_glb[MyF2][BEAD_FACE] = bead1;
                }
        }
    // Swap them on the lattice
    naTotLattice_glb[Lat_Ind_FromVec(tmpR)]  = bead2;
    naTotLattice_glb[Lat_Ind_FromVec(tmpR2)] = bead1;
}

/// OP_GenMHValue - calculates the acceptance ratio for a given state.
/// In this implementation, it is assumed that all the inputs to the function
/// correspond to the logl() of the values (to keep the numbers small)
/// If the total value x, in log space, is larger than 0, it means that expl(x)
/// > 1. Therefore the move will be accepted regardless, so the function shall
/// return 2 (a value larger than 1). Similarly, if x < -21.5, expl(x)
/// ~1x10^{-10} which is smaller than 1/RAND_MAX; so it is essentially 0 -- thus
/// the function returns 0. In all other cases, the function returns expl(x),
/// which is then compared to a number between 0 and 1. Remember that
/// log(exp(a)) = a, so log(exp((old_en-new_en)/myTemp)) =
/// (old_en-new_en)/myTemp
/// \param fRos: Forward Rosenbluth weight.
/// \param bRos: Backwards Rosenbluth weight.
/// \param Delta_En: Difference in energy.
/// \param Cur_Temp: Current temperature.
/// \return The new proposed bonding partner.
lLDub OP_GenMHValue(lLDub fRos, lLDub bRos, lLDub Delta_En, lLDub Cur_Temp)
{
    lLDub MH_Value;

    MH_Value = fRos - bRos + Delta_En / Cur_Temp;
    if (MH_Value <= ldLogOfSmallestPossibleProb_glb)
        {
            MH_Value = 0.;
        }
    else if (MH_Value >= 0.)
        {
            MH_Value = 2.0;
        }
    else
        {
            MH_Value = expl(MH_Value);
        }
    return MH_Value;
}

///
/// \param firstB
/// \param lastB
/// \param r_posNew
void OP_System_Snake_SlitherFwd(const int firstB, const int lastB, const int* r_posNew)
{
    int i;
    int r_tmp[POS_MAX];
    naTotLattice_glb[Lat_Ind_OfBead(firstB)] = -1; // Only the firstB's location is empty
    for (i = firstB; i < lastB - 1; i++)
        {
            LatPos_copy(bead_info_glb[i], old_bead_glb[i + 1]); // Hopping over by one bead
            LatPos_copy(r_tmp, bead_info_glb[i]);
            naTotLattice_glb[Lat_Ind_FromVec(r_tmp)] = i;
        }
    // Moving the last bead, and it has to be done independently because lastB-1
    // -> r_posNew
    i = lastB - 1;
    LatPos_copy(bead_info_glb[i], r_posNew);
    naTotLattice_glb[Lat_Ind_FromVec(r_posNew)] = i;
}

void OP_System_Snake_SlitherBck(const int firstB, const int lastB, const int* r_posNew)
{
    int i;
    int r_tmp[POS_MAX];
    naTotLattice_glb[Lat_Ind_OfBead(lastB - 1)] = -1; // Only the lastB-1's location is empty
    for (i = firstB + 1; i < lastB; i++)
        {
            LatPos_copy(bead_info_glb[i], old_bead_glb[i - 1]); // Hopping back by one bead
            LatPos_copy(r_tmp, bead_info_glb[i]);
            naTotLattice_glb[Lat_Ind_FromVec(r_tmp)] = i;
        }

    // Moving the first bead, and it has to be done independently because firstB
    // -> r_posNew
    i = firstB;
    LatPos_copy(bead_info_glb[i], r_posNew);
    naTotLattice_glb[Lat_Ind_FromVec(r_posNew)] = i;
}

void OP_Beads_CopyBeadsInListToOld(const int listSize, const int* beadList)
{

    int i, thisBead;
    for (i = 0; i < listSize; i++)
        {
            thisBead = beadList[i];
            OP_CopyBead(old_bead_glb[thisBead], bead_info_glb[thisBead]);
        }
}

void OP_Beads_CopyBeadsInListFromOld(const int listSize, const int* beadList)
{

    int i, thisBead;
    for (i = 0; i < listSize; i++)
        {
            thisBead = beadList[i];
            OP_CopyBead(bead_info_glb[thisBead], old_bead_glb[thisBead]);
        }
}

void OP_Beads_CopyBeadsInListToPosList(const int listSize, const int* beadList, int (*beadPos)[POS_MAX])
{

    int i, thisBead;
    for (i = 0; i < listSize; i++)
        {
            thisBead = beadList[i];
            LatPos_copy(beadPos[i], bead_info_glb[thisBead]);
        }
}

void OP_Lattice_EmptySitesForListOfBeads(const int listSize, const int* beadList)
{
    int i, thisBead;
    for (i = 0; i < listSize; i++)
        {
            thisBead                                   = beadList[i];
            naTotLattice_glb[Lat_Ind_OfBead(thisBead)] = -1;
        }
}

void OP_Lattice_PlaceBeadsInList(const int listSize, const int* beadList)
{
    int i, thisBead;
    for (i = 0; i < listSize; i++)
        {
            thisBead                                   = beadList[i];
            naTotLattice_glb[Lat_Ind_OfBead(thisBead)] = thisBead;
        }
}

void OP_Lattice_EmptySitesForListOfPos(const int listSize, const int (*beadPos)[POS_MAX])
{
    int i;
    for (i = 0; i < listSize; i++)
        {
            naTotLattice_glb[Lat_Ind_FromVec(beadPos[i])] = -1;
        }
}

void OP_Beads_MoveBeadsInListToPos(const int listSize, const int* beadList, const int (*newPos)[POS_MAX])
{
    int i, tmpBead;
    for (i = 0; i < listSize; i++)
        {
            tmpBead = beadList[i];
            LatPos_copy(bead_info_glb[tmpBead], newPos[i]);
        }
}

void OP_Beads_BreakBondsInList(const int listSize, const int* beadList)
{
    int i, tmpID, bP;
    for (i = 0; i < listSize; i++)
        {
            tmpID = beadList[i];
            OP_Beads_BreakBond(tmpID);
        }
}

void OP_Beads_BreakBond(const int beadID)
{
    const int bP = bead_info_glb[beadID][BEAD_FACE];
    if (bP != -1)
        {
            bead_info_glb[bP][BEAD_FACE]     = -1;
            bead_info_glb[beadID][BEAD_FACE] = -1;
        }
}

void OP_Beads_RestoreBondsInList(const int listSize, const int* beadList)
{
    int i, tmpID, bP;
    for (i = 0; i < listSize; i++)
        {
            tmpID = beadList[i];
            OP_Beads_RestoreBond(tmpID);
        }
}

void OP_Beads_RestoreBond(const int beadID)
{
    const int bP = bead_info_glb[beadID][BEAD_FACE];
    if (bP != -1)
        {
            bead_info_glb[bP][BEAD_FACE] = beadID;
        }
}

void OP_Inv_MoveBeads_InList_ToPos(const int listSize, const int* beadList)
{
    int i, tmpBead;
    for (i = 0; i < listSize; i++)
        {
            tmpBead = beadList[i];
            LatPos_copy(bead_info_glb[tmpBead], old_bead_glb[tmpBead]);
        }
}

lLDub MC_RosenbluthSampling_ForLocal_AtOld(const int beadID, const int resi, long double* oldEn, const int neigh_num)
{
    int ros_num;
    if (nBeadTypeIsSticker_glb[resi])
        {
            // the anisotropic energy should not be added to the oldEn
            // because it is already counted in this biased smapling.
            // *oldEn  = *oldEn + Energy_Anisotropic(beadID);
            ros_num = Check_RotStates_wNeighList(beadID, resi, naOldOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);
            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForLocal_AtNew(const int beadID, const int resi, int* bead_part, long double* newEn,
                                           const int neigh_num)
{
    int ros_num;

    *bead_part = -1;
    if (nBeadTypeIsSticker_glb[resi])
        {
            // OP_ShuffleArray(neigh_num, naNewOvlpNeighs_glb);
            ros_num = Check_RotStates_wNeighList(beadID, resi, naNewOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);

            *bead_part = OP_PickRotState(ros_num);

            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForChains_AtOld(const int beadID, const int resi, long double* oldEn, const int neigh_num)
{
    int ros_num;
    if (nBeadTypeIsSticker_glb[resi])
        {
            // the anisotropic energy should not be added to the oldEn
            // because it is already counted in this biased smapling.
            // *oldEn  = *oldEn + Energy_Anisotropic_For_Chain(beadID); 
            ros_num = Check_RotStates_wNeighList_ForChains(beadID, resi, naOldOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);
            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForChains_AtNew(const int beadID, const int resi, int* bead_part, long double* newEn,
                                            const int neigh_num)
{
    int ros_num;
    const int cur_part = bead_info_glb[beadID][BEAD_FACE];
    *bead_part         = -1;
    if (nBeadTypeIsSticker_glb[resi])
        {
            // OP_ShuffleArray(neigh_num, naNewOvlpNeighs_glb);
            // ros_num = Check_RotStates_wNeighList_ForChains(beadID, resi, naNewOvlpNeighs_glb, neigh_num);
            ros_num = Check_RotStates_wNeighList_ForChains(beadID, resi, naNewOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);
            if (cur_part == -1)
                {
                    *bead_part = OP_PickRotState(ros_num);
                }
            return logl(ldaBoltzNorm_glb[0]); 
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForLists_AtOld(const int beadIdx, const int listSize, const int beadList[MAX_BONDS + 1],
                                           lLDub* oldEn, const int neigh_num)
{
    const int beadID = beadList[beadIdx];
    const int resi   = bead_info_glb[beadID][BEAD_TYPE];
    int ros_num;
    if (nBeadTypeIsSticker_glb[resi])
        {
            // the anisotropic energy should not be added to the oldEn
            // because it is already counted in this biased smapling.
            // *oldEn  = *oldEn + Energy_Anisotropic_For_List(beadID, listSize, beadList);
            // ros_num = Check_RotStates_wNeighList(beadID, resi, naOldOvlpNeighs_glb, neigh_num);
            ros_num = Check_RotStates_wNeighList_ForLists(beadIdx, listSize, beadList, resi, naOldOvlpNeighs_glb, neigh_num);           
            OP_NormalizeRotState(0, ros_num);
            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForLists_AtNew(const int beadIdx, const int listSize, const int beadList[MAX_BONDS + 1],
                                           int* bead_part, lLDub* newEn, const int neigh_num)
{
    *bead_part         = -1;
    const int beadID   = beadList[beadIdx];
    const int resi     = bead_info_glb[beadID][BEAD_TYPE];
    const int cur_part = bead_info_glb[beadID][BEAD_FACE];
    int ros_num;
    if (nBeadTypeIsSticker_glb[resi])
        {
            // the anisotropic energy should not be added to the oldEn
            // because it is already counted in this biased smapling.
            // *newEn  = *newEn + Energy_Anisotropic_For_List(beadID, listSize, beadList);
            // ros_num = Check_RotStates_wNeighList(beadID, resi, naNewOvlpNeighs_glb, neigh_num);
            ros_num = Check_RotStates_wNeighList_ForLists(beadIdx, listSize, beadList, resi, naNewOvlpNeighs_glb, neigh_num);           
            OP_NormalizeRotState(0, ros_num);
            if (cur_part == -1)
                {
                    *bead_part = OP_PickRotState(ros_num);
                }

            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForRange_AtOld(const int beadID, const int resi, const int smallestBead,
                                           const int largestBead, lLDub* oldEn, const int neigh_num)
{
    int ros_num;
    if (nBeadTypeIsSticker_glb[resi])
        {
            // the anisotropic energy should not be added to the oldEn
            // because it is already counted in this biased smapling.
            // *oldEn  = *oldEn + Energy_Anisotropic_For_Range(beadID, smallestBead, largestBead);
            // ros_num = Check_RotStates_wNeighList(beadID, resi, naOldOvlpNeighs_glb, neigh_num);
            ros_num = Check_RotState_wNeighList_ForRange(beadID, resi, smallestBead, largestBead, naOldOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);
            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForRange_AtNew(const int beadID, const int resi, const int smallestBead, 
                                           const int largestBead, int* bead_part, lLDub* newEn, const int neigh_num)
{
    int ros_num;
    const int cur_part = bead_info_glb[beadID][BEAD_FACE];
    *bead_part         = -1;
    if (nBeadTypeIsSticker_glb[resi])
        {
            // OP_ShuffleArray(neigh_num, naNewOvlpNeighs_glb);
            // ros_num = Check_RotStates_wNeighList(beadID, resi, naNewOvlpNeighs_glb, neigh_num);
            ros_num = Check_RotState_wNeighList_ForRange(beadID, resi, smallestBead, largestBead, naNewOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);
            if (cur_part == -1)
                {
                    *bead_part = OP_PickRotState(ros_num);
                }
            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}
