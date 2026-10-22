#include "mcmove_internal.h"
#include "energy.h"
#include "structure.h"

// All the _Equil variants of the moves are spatially the same as their non
// _Equil variants. The only difference is that the anisotropic interaction is
// ignored, and thus no bias is applied.

int Move_Local_Equil(int beadID, float MyTemp)
{ // Performs a local translation MC-move on beadID

    int bAccept = 0; // Used in MC steps
    int yTemp;
    int r_pos0[POS_MAX], r_posNew[POS_MAX],
        r_disp[POS_MAX]; // Vectors to stores coordinates.
    // printf("Beginning LOCAL\n");
    LatPos_copy(r_pos0, g_system.bead_info[beadID]);
    //    LatPos_gen_rand_wRad(r_disp, g_system.linker_len[beadID][0]);
    LatPos_gen_rand_wRad(r_disp, 2);
    LatPos_add_wPBC(r_posNew, r_pos0, r_disp);

    yTemp = Check_MoveBeadTo(r_posNew);

    if (yTemp == 1)
    { // This means we found an empty lattice site. So let's
      // check if the linkers are okay.
        yTemp = Check_LinkerConstraint(beadID, r_posNew);
    }

    if (yTemp != 1)
    {
        // This means that we have failed to find an appropriate spot for this
        // bead to be moved to.
        //  Therefore, the move is rejected!
        bAccept = 0;
        // printf("End LOCAL - No space\n");
        return bAccept;
    }
    // Have successfully found a good lattice spot.
    lLDub MCProb, oldEn, newEn; // For Metropolis Hastings
    oldEn = 0.;
    newEn = 0.;

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    const int resi = g_system.bead_info[beadID][BEAD_TYPE];

    oldEn = g_system.nBiasPotential_Mode == -1 ? 0.f : Energy_BiasingPotential(beadID);
    newEn = 0.;

    Energy_Iso_ForLocalEquil(
        beadID,
        resi,
        r_pos0,
        &oldEn,
        &newEn,
        &old_ovlp_num,
        &old_cont_num,
        g_system.naOldOvlpNeighs,
        g_system.naOldContNeighs);

    int bondList[MAX_BONDS + 1];
    int bondNum = g_system.bSystemHasTopo ? OP_GetTopoBonds(beadID, bondList) : 0;
    bondNum = BeadList_CanTopoAngle(bondNum, bondList);

    oldEn += bondNum ? Energy_Topo_Angle_ForList(bondNum, bondList) : 0.;

    OP_System_MoveBeadTo(beadID, r_posNew);

    newEn += g_system.nBiasPotential_Mode == -1 ? 0.f : Energy_BiasingPotential(beadID);

    Energy_Iso_ForLocalEquil(
        beadID,
        resi,
        r_posNew,
        &newEn,
        &oldEn,
        &new_ovlp_num,
        &new_cont_num,
        g_system.naNewOvlpNeighs,
        g_system.naNewContNeighs);

    newEn += bondNum ? Energy_Topo_Angle_ForList(bondNum, bondList) : 0.;

    MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub)MyTemp);
    if (MCProb < MHAcc)
    {                // Accept this state
        bAccept = 1; // Accepting!
        // printf("End LOCAL - Accept\n");
        return bAccept;
    }
    else
    {
        OP_System_MoveBeadTo_Inv(beadID);
        bAccept = 0;
        // printf("End LOCAL - Fail\n");
        return bAccept;
    }
}

int Move_Snake_Equil(int chainID, float MyTemp)
{ // Performs a slither MC-move on chainID

    int bAccept = 0; // Used in MC steps
    // Finding the bounds for looping over the molecule/chain
    const int firstB = g_system.chain_info[chainID][CHAIN_START];
    const int lastB = firstB + g_system.chain_info[chainID][CHAIN_LENGTH];
    if (lastB - firstB == 1)
    { // This means we have a monomer. Reject the
      // move, because Local or Trans
        // moves should be the ones that move monomers.
        bAccept = 0;
        return bAccept;
    }
    else
    {
        if (g_system.nChainTypeIsLinear[g_system.chain_info[chainID][CHAIN_TYPE]] != 1)
        { // If chain is not linear. Reject move because
          // slithering will not werk!
            bAccept = 0;
            return bAccept;
        }
    }
    // This chain is suitable to perform a slithering-snake move on.

    int i;     // Loop iterators
    int yTemp; // Random numbers to store things
    int r_posNew[POS_MAX], r_posTmp1[POS_MAX],
        r_posTmp2[POS_MAX];       // Vectors to store positions.
    int angBead_old, angBead_new; // BeadIDs for angular energies.

    lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX; // To decide if we slither forwards or backwards
    const char SnakeFwd = MCProb < 0.5f ? 1 : 0;    // 1 For forwards, 0 for backwards.
    if (SnakeFwd)
    { // Forwards slither, so lastB-1 (last bead) is anchor
        LatPos_gen_rand_wRad(r_posNew,
                             g_system.linker_len[lastB - 1][0]); // lastB-1 will be replaced by lastB-2
        LatPos_copy(r_posTmp1, r_posNew);
        LatPos_add_wPBC(r_posNew, g_system.bead_info[lastB - 1], r_posTmp1);
        yTemp = Check_MoveBeadTo(r_posNew); // 0: there is no space, 1: there is space
        angBead_old = firstB + 1;
        angBead_new = lastB - 2;
    }
    else
    { // Backwards slither, so firstB is anchor
        LatPos_gen_rand_wRad(r_posNew,
                             g_system.linker_len[firstB][0]); // firstB will be replaced by firstB+1
        LatPos_copy(r_posTmp1, r_posNew);
        LatPos_add_wPBC(r_posNew, g_system.bead_info[firstB], r_posTmp1);
        yTemp = Check_MoveBeadTo(r_posNew); // 0: there is no space, 1: there is space
        angBead_new = firstB + 1;
        angBead_old = lastB - 2;
    }

    if (yTemp == 0)
    { // Couldn't find a spot, so reject the damn move
        bAccept = 0;
        return bAccept;
    }
    // We should have a spot to move to! r_posNew has the location

    // Let's remember where this chain exists.
    OP_CopyBeadsToOld(firstB, lastB);

    lLDub newEn = 0.;
    lLDub oldEn = 0.;
    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = firstB; i < lastB; i++)
        {
            oldEn += Energy_BiasingPotential(i);
        }
    }

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    for (i = firstB; i < lastB; i++)
    {
        Energy_Iso_ForChainsEquil(
            i,
            &oldEn,
            &newEn,
            &old_ovlp_num,
            &old_cont_num,
            g_system.naOldOvlpNeighs,
            g_system.naOldContNeighs);
    }
    const int resi = g_system.bead_info[angBead_old][BEAD_TYPE];
    oldEn += g_energy.faEnergy[resi][resi][E_STIFF] ? Energy_Topo_Angle(angBead_old) : 0.;

    if (SnakeFwd)
    { // Slithering the chain forwards in ID-space
        OP_System_Snake_SlitherFwd(firstB, lastB, r_posNew);
    }
    else
    { // Slithering backwards in ID-space
        OP_System_Snake_SlitherBck(firstB, lastB, r_posNew);
    }

    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = firstB; i < lastB; i++)
        {
            newEn += Energy_BiasingPotential(i);
        }
    }

    lLDub FSum = 0.;
    for (i = firstB; i < lastB; i++)
    {
        Energy_Iso_ForChainsEquil(
            i,
            &newEn,
            &oldEn,
            &new_ovlp_num,
            &new_cont_num,
            g_system.naNewOvlpNeighs,
            g_system.naNewContNeighs);
    }

    newEn += g_energy.faEnergy[resi][resi][E_STIFF] ? Energy_Topo_Angle(angBead_new) : 0.;

    // Doing the Metropolis-Hastings thing
    MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub)MyTemp);
    if (MCProb < MHAcc)
    { // Accept the move. Remember that the bonds were assigned
      // above!Accept. Bonds have been handled before!
        bAccept = 1;
        return bAccept;
    }
    else
    {
        OP_RestoreChain_ForSnake(firstB, lastB);
        bAccept = 0;
        return bAccept;
    }
}

int Move_Trans_Equil(int chainID, float MyTemp)
{                               // Performs a translation move with orientational bias
    int bAccept = 0;            // Used in MC steps
    lLDub MCProb, oldEn, newEn; // For Metropolis Hastings
    oldEn = 0.;
    newEn = 0.;
    int i, j; // Loop iterators
    int firstB, lastB;
    int yTemp;           // Random numbers to store things
    int r_disp[POS_MAX]; // Vectors to store coordinates.
    // Finding the bounds for looping over the molecule/chain
    firstB = g_system.chain_info[chainID][CHAIN_START];
    lastB = firstB + g_system.chain_info[chainID][CHAIN_LENGTH];
    // Radii for translation moves. All moves are L/4 radius

    LatPos_gen_rand_wRad(r_disp, g_system.naBoxSize[0] / 2);

    yTemp = Check_ChainDisp(chainID, r_disp); // yTemp=0 means clash
    if (yTemp == 0)
    { // We have failed to find a good spot for this chain.
        bAccept = 0;
        // printf("Ending TRANS no space\n");
        return bAccept;
    }
    // We now have a chain which when moved does not overlap.

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    oldEn = 0.;
    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = firstB; i < lastB; i++)
        {
            oldEn += Energy_BiasingPotential(i);
        }
    }

    for (i = firstB; i < lastB; i++)
    {
        Energy_Iso_ForChainsEquil(
            i,
            &oldEn,
            &newEn,
            &old_ovlp_num,
            &old_cont_num,
            g_system.naOldOvlpNeighs,
            g_system.naOldContNeighs);
    }

    OP_System_DispChain_ForTrans(chainID,
                                 r_disp); // Moved the chain, broke bonds, and remembered stuff

    newEn = 0.;
    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = firstB; i < lastB; i++)
        {
            newEn += Energy_BiasingPotential(i);
        }
    }

    for (i = firstB; i < lastB; i++)
    {
        Energy_Iso_ForChainsEquil(
            i,
            &newEn,
            &oldEn,
            &new_ovlp_num,
            &new_cont_num,
            g_system.naNewOvlpNeighs,
            g_system.naNewContNeighs);
    }

    MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub)MyTemp);
    if (MCProb < MHAcc)
    { // Accept the move. Remember that the bonds were
      // assigned above!
        bAccept = 1;
        return bAccept;
    }
    else
    {
        OP_System_RestoreChain_ForTrans(chainID);
        bAccept = 0;
        return bAccept;
    }
}

int Move_MultiLocal_Equil(int beadID, float MyTemp)
{
    int beadsList[MAX_BONDS + 1];
    const int beadNum = OP_GetTopoBonds(beadID, beadsList);
    int bAccept;

    if (beadNum < 2)
    { // Bead must be bonded to at least 2 things.
        bAccept = 0;
        return bAccept;
    }

    int beads_pos0[MAX_BONDS + 1][POS_MAX];
    int beads_posNew[MAX_BONDS + 1][POS_MAX];
    int r_posTmp1[POS_MAX];
    int tmpBead;
    int i;
    int yTemp;

    OP_Beads_CopyBeadsInListToOld(beadNum, beadsList);
    OP_Beads_CopyBeadsInListToPosList(beadNum, beadsList, beads_pos0);
    OP_Lattice_EmptySitesForListOfPos(beadNum, beads_pos0);

    for (i = 0; i < beadNum; i++)
    {
        LatPos_gen_rand_wRad(r_posTmp1, 2);
        LatPos_add_wPBC(beads_posNew[i], r_posTmp1, beads_pos0[i]);
        yTemp = Check_MoveBeadTo(beads_posNew[i]);
        if (yTemp == 0)
        {
            break;
        }
        tmpBead = beadsList[i];
        g_system.naTotLattice[Lat_Ind_FromVec(beads_posNew[i])] = tmpBead;
    }

    OP_Lattice_EmptySitesForListOfPos(i, beads_posNew);

    if (yTemp == 1)
    { // No steric clash so check for topology constraint.
        OP_Beads_MoveBeadsInListToPos(beadNum, beadsList, beads_posNew);
        yTemp = Check_LinkerConstraints_ForBeadList(beadNum, beadsList);
        OP_Beads_MoveBeadsInListToPos(beadNum, beadsList, beads_pos0);
    }

    OP_Lattice_PlaceBeadsInList(beadNum, beadsList);

    if (yTemp == 0)
    {
        bAccept = 0;
        return bAccept;
    }

    lLDub oldEn = 0.;

    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            oldEn += Energy_BiasingPotential(tmpBead);
        }
    }

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    lLDub newEn = 0.;
    for (i = 0; i < beadNum; i++)
    {
        Energy_Iso_ForLists(
            i,
            beadNum,
            beadsList,
            beads_pos0,
            &oldEn,
            &newEn,
            &old_ovlp_num,
            &old_cont_num,
            g_system.naOldOvlpNeighs,
            g_system.naOldContNeighs);
    }

    int bondAngNum = 0;
    int tmpAngList[MAX_BONDS + 1];
    int bondAngList[(MAX_BONDS + 1) * (MAX_BONDS + 1)];
    if (g_system.bSystemHasTopo)
    {
        bondAngNum = BeadList_AppendBeads(0, bondAngList, beadsList, beadNum);
        for (i = 1; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            yTemp = OP_GetTopoBonds(tmpBead, tmpAngList);
            bondAngNum = BeadList_AppendBeads(bondAngNum, bondAngList, tmpAngList, yTemp);
        }
        qsort(bondAngList, bondAngNum, sizeof(int), UtilFunc_CompareInts);
        bondAngNum = ListOP_UniqueElementsOfSortedList_Int(bondAngNum, bondAngList);
    }
    else
    {
        bondAngNum = 0;
    }

    oldEn += bondAngNum ? Energy_Topo_Angle_ForList(bondAngNum, bondAngList) : 0.;

    OP_System_MoveBeadsInListToPos(beadNum, beadsList, beads_posNew);

    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            newEn += Energy_BiasingPotential(tmpBead);
        }
    }

    for (i = 0; i < beadNum; i++)
    {
        Energy_Iso_ForLists(
            i,
            beadNum,
            beadsList,
            beads_posNew,
            &newEn,
            &oldEn,
            &new_ovlp_num,
            &new_cont_num,
            g_system.naNewOvlpNeighs,
            g_system.naNewContNeighs);
    }

    newEn += bondAngNum ? Energy_Topo_Angle_ForList(bondAngNum, bondAngList) : 0.;

    lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub)MyTemp);

    if (MCProb < MHAcc)
    { // Accept the move. Remember that the bonds were
      // assigned above
        bAccept = 1;
        return bAccept;
    }
    else
    {
        OP_Lattice_EmptySitesForListOfBeads(beadNum, beadsList);
        OP_Beads_CopyBeadsInListFromOld(beadNum, beadsList);
        OP_Lattice_PlaceBeadsInList(beadNum, beadsList);
        bAccept = 0;
        return bAccept;
    }
}

int Move_Pivot_Equil(int chainID, float MyTemp)
{
    // Performs a pivot move on chainID
    /*
    Randomly pick a bead in this chain (anchorBead), and perform a symmetry
    operation on the beads after the anchorBead. Note that if anchorBead is the
    second-to-last, or last, bead, the move is rejected. Local and slither moves
    would be better for those. The protein must be linear and the move is
    rejected if the chain is branched. It is assumed that a linear chain has
    been passed!
    */
    int bAccept = 0;
    if (!g_system.nChainTypeIsLinear[g_system.chain_info[chainID][CHAIN_TYPE]])
    {
        return bAccept;
    }
    // Check if the chain is longer than 3 or if it is linear
    const int chainLength = g_system.chain_info[chainID][CHAIN_LENGTH];
    if (chainLength <= 3)
    {
        return bAccept;
    }

    int PivotM;
    PivotM = rand() % 10;
    if (PivotM == 0)
    { // Null move
        bAccept = 1;
        return bAccept;
    }

    // Get the beadID's for the first and last bead
    const int firstB = g_system.chain_info[chainID][CHAIN_START];
    const int lastB = firstB + chainLength;

    // Randomly select a bead that is neither the first nor last
    int anchorBead = chainLength - 2;
    anchorBead = 1 + (rand() % anchorBead);

    //-1 Means backwards, +1 means forwards. Always Pivot the smaller portion
    const int PivotDir = anchorBead > chainLength / 2 ? 1 : -1;
    anchorBead = firstB + anchorBead;
    // printf("Bead: %d\n", anchorBead);

    int i, j;
    int xTemp, yTemp;
    int anchorPos[POS_MAX];
    int beadsList[MAX_CHAINLEN];
    int beadNum = 0;

    if (PivotDir == 1)
    {
        for (i = anchorBead + 1; i < lastB; i++)
        {
            beadsList[beadNum++] = i;
        }
    }
    else
    {
        for (i = firstB; i < anchorBead; i++)
        {
            beadsList[beadNum++] = i;
        }
    }

    const int smBead = beadsList[0];
    const int lgBead = beadsList[beadNum - 1];

    LatPos_copy(anchorPos, g_system.bead_info[anchorBead]);

    int tmpBead;

    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        // Use anchorbead as origin, and check what happens after the rotation
        // to every bead
        OP_Rotation(PivotM, tmpBead, anchorPos);
        yTemp = Check_MoveBeadTo(g_system.naTempR);
        if (yTemp == 0)
        {
            break;
        }
    }
    if (yTemp == 0)
    {
        bAccept = 0;
        return bAccept;
    }

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            oldEn += Energy_BiasingPotential(tmpBead);
        }
    }

    int old_ovlp_num, old_cont_num;

    lLDub BSum = 0.;
    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        Energy_Iso_ForRangeEquil(
            tmpBead,
            smBead,
            lgBead,
            &oldEn,
            &newEn,
            &old_ovlp_num,
            &old_cont_num,
            g_system.naOldOvlpNeighs,
            g_system.naOldContNeighs);
    }

    const int resi = g_system.bead_info[anchorBead][BEAD_TYPE];
    oldEn += g_energy.faEnergy[resi][resi][E_STIFF] ? Energy_Topo_Angle(anchorBead) : 0.;

    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        OP_Rotation(PivotM, tmpBead, anchorPos);
        OP_System_MoveBeadTo(tmpBead, g_system.naTempR);
    }

    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            newEn += Energy_BiasingPotential(tmpBead);
        }
    }

    int new_ovlp_num, new_cont_num;

    lLDub FSum = 0.;
    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        Energy_Iso_ForRangeEquil(
            tmpBead,
            smBead,
            lgBead,
            &newEn,
            &oldEn,
            &new_ovlp_num,
            &new_cont_num,
            g_system.naNewOvlpNeighs,
            g_system.naNewContNeighs);
    }

    newEn += g_energy.faEnergy[resi][resi][E_STIFF] ? Energy_Topo_Angle(anchorBead) : 0.;

    lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(FSum, BSum, oldEn - newEn, (lLDub)MyTemp);
    if (MCProb < MHAcc)
    { // Accept the move. Remember that the bonds were
      // assigned above!
        bAccept = 1;
        return bAccept;
    }
    else
    { // Rejecting move
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            OP_System_MoveBeadTo_Inv(tmpBead);
        }
        bAccept = 0;
        return bAccept;
    }
}

int Move_BranchedRot_Equil(int chainID, float MyTemp)
{
    // Rotates a branched molecule about the branching, which is assumed to be
    // the firstB of chainID Performs a Move_Pivot() on molecule where the
    // rotation occurs around firstB
    /*
      Set the first bead as anchorBead, and perform a symmetry operation on the
      beads after the anchorBead (the whole molecule). Note that if the molecule
      is linear, the move is outright rejected. Again, it is assumed that the
      first bead in that molecule is the 'node'.
      */
    int bAccept = 0;
    // Reject if the molecule is linear
    if (g_system.nChainTypeIsLinear[g_system.chain_info[chainID][CHAIN_TYPE]])
    {
        return bAccept;
    }
    const int chainLength = g_system.chain_info[chainID][CHAIN_LENGTH];

    int PivotM;
    PivotM = rand() % 10;
    if (PivotM == 0)
    { // Null move
        bAccept = 1;
        return bAccept;
    }

    // Get the beadID's for the first and last bead
    const int firstB = g_system.chain_info[chainID][CHAIN_START];
    const int lastB = firstB + chainLength;

    // FirstB is the anchor.
    const int anchorBead = firstB;

    int i;
    int yTemp;
    int anchorPos[POS_MAX];
    int beadsList[MAX_CHAINLEN];
    int beadNum = 0;

    for (i = anchorBead + 1; i < lastB; i++)
    {
        beadsList[beadNum++] = i;
    }

    const int smBead = beadsList[0];
    const int lgBead = beadsList[beadNum - 1];

    LatPos_copy(anchorPos, g_system.bead_info[anchorBead]);

    int tmpBead;

    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        // Use anchorbead as origin, and check what happens after the rotation
        // to every bead
        OP_Rotation(PivotM, tmpBead, anchorPos);
        yTemp = Check_MoveBeadTo(g_system.naTempR);
        if (yTemp == 0)
        {
            break;
        }
    }
    if (yTemp == 0)
    {
        bAccept = 0;
        return bAccept;
    }

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            oldEn += Energy_BiasingPotential(tmpBead);
        }
    }

    int old_ovlp_num, old_cont_num;

    lLDub BSum = 0.;
    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        Energy_Iso_ForRangeEquil(
            tmpBead,
            smBead,
            lgBead,
            &oldEn,
            &newEn,
            &old_ovlp_num,
            &old_cont_num,
            g_system.naOldOvlpNeighs,
            g_system.naOldContNeighs);
    }

    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        OP_Rotation(PivotM, tmpBead, anchorPos);
        OP_System_MoveBeadTo(tmpBead, g_system.naTempR);
    }

    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            newEn += Energy_BiasingPotential(tmpBead);
        }
    }

    int new_ovlp_num, new_cont_num;

    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        Energy_Iso_ForRangeEquil(
            tmpBead,
            smBead,
            lgBead,
            &newEn,
            &oldEn,
            &new_ovlp_num,
            &new_cont_num,
            g_system.naNewOvlpNeighs,
            g_system.naNewContNeighs);
    }

    lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub)MyTemp);
    if (MCProb < MHAcc)
    { // Accept the move. Remember that the bonds were
      // assigned above!
        bAccept = 1;
        return bAccept;
    }
    else
    { // Rejecting move
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            OP_System_MoveBeadTo_Inv(tmpBead);
        }
        bAccept = 0;
        return bAccept;
    }
}

/// Check_ChainDisp - checks if moving chainID by vec_disp[POS_MAX] will lead to
/// sterix clash.
/// \param chainID
