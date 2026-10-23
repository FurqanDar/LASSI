#include "mcmove_internal.h"
#include "cluster_mc.h"
#include "energy.h"
#include "structure.h"

/*
 * The following functions are sub-routines that perform the various Monte-Carlo
 * (MC) moves.
 */

/// Move_Rot - searches the 3^3-1=26 possible sites for forming a bond, and
/// performs Metropolis-Hastings. Rejects move if beadID cannot forms bonds.
/// \param beadID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_Rot(int beadID, float MyTemp)
{
    // Performs a rotational MC-Move on beadID
    int bAccept;    // Used in MC steps
    int resi, resj; // To track bead types.
    // Firstly -- make sure that bead i can even have rotational states.
    resi = g_system.bead_info[beadID][BEAD_TYPE];
    // printf("Beginning ROT\n");
    if (g_system.nBeadTypeIsSticker[resi] == 0)
    { // Skip beads that cannot rotate!
        bAccept = 0;
        return bAccept;
    }

    lLDub MCProb, oldEn, newEn; // For Metropolis Hastings
    oldEn = 0.;
    newEn = 0.;
    int i; // General looping iterators
    int yTemp;
    int CandNum;

    if (g_system.bead_info[beadID][BEAD_FACE] != -1)
    {
        resj =
            g_system.bead_info[g_system.bead_info[beadID][BEAD_FACE]][BEAD_TYPE]; // This is who I'm currently bonded to
        oldEn = (lLDub)g_energy.faEnergy[resi][resj][E_SC_SC];
    }
    // OP_ShuffleRotIndecies();
    CandNum = Check_RotStatesOld(beadID, resi, MyTemp);
    OP_NormalizeRotState(0, CandNum);

    yTemp = OP_PickRotState(CandNum);

    if (yTemp != -1)
    { // There is a bead here
        resj = g_system.bead_info[yTemp][BEAD_TYPE];
        newEn = (lLDub)g_energy.faEnergy[resi][resj][E_SC_SC];
    }
    // See if we can accept this move
    MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    // There is no need to calculate the Rosenbluth weights for this case, because FW and BW are the
    // same.
    lLDub MHAcc = OP_GenMHValue(0., 0., 0., (lLDub)MyTemp);

    if (MCProb < MHAcc)
    { // Accept this state
        if (g_system.bead_info[beadID][BEAD_FACE] != -1)
        {                                                                              // Break old bond
            g_system.bead_info[g_system.bead_info[beadID][BEAD_FACE]][BEAD_FACE] = -1; // Breaking bond with old partner
        }
        g_system.bead_info[beadID][BEAD_FACE] = yTemp;
        if (yTemp != -1)
        { // New bond
            g_system.bead_info[yTemp][BEAD_FACE] = beadID;
        }
        bAccept = 1;
        return bAccept;
    }
    else
    {
        bAccept = 0;
        return bAccept;
    }
}

/// Move_Local - performs a biased local move on beadID by:
/// 1. Seeing if there is a free spot to move to in a <+-2,+-2,+-2> random
/// location. Rejects if no space found/exists.
/// 2. Calculating the Rosenbluth weight by searching the 3^3-1=26 possible
/// bonding locations, and calculate the rest of the energy.
/// 3. Move the bead. Calculate the Rosenbluth weight, and energy. Propose a
/// bond and perform Metropolis-Hastings
/// \param beadID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_Local(int beadID, float MyTemp)
{ // Performs a local translation MC-move on beadID

    int bAccept = 0;                                         // Used in MC steps
    int yTemp;                                               // Random numbers to store things
    int r_pos0[POS_MAX], r_posNew[POS_MAX], r_disp[POS_MAX]; // Vectors to stores coordinates.

    // Attempt to find an empty lattice point.
    LatPos_copy(r_pos0, g_system.bead_info[beadID]);
    LatPos_gen_rand_wRad(r_disp, 2);
    //    LatPos_gen_rand_wRad(r_disp, g_system.linker_len[beadID][0]);
    LatPos_add_wPBC(r_posNew, r_pos0, r_disp);

    // Checking to see validity of new point.
    yTemp = Check_MoveBeadTo(r_posNew);
    if (yTemp == 1)
    { // This means we found an empty lattice site. So let's
      // check if the linkers are okay.
        yTemp = Check_LinkerConstraint(beadID, r_posNew);
    }
    if (yTemp != 1)
    { // No space -- reject.
        bAccept = 0;
        return bAccept;
    }

    // Have successfully found a good lattice spot. Let's perform the usual
    // Metropolis-Hastings shenanigans.

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    const int resi = g_system.bead_info[beadID][BEAD_TYPE];

    lLDub oldEn = g_system.nBiasPotential_Mode == -1 ? 0. : Energy_BiasingPotential(beadID);
    lLDub newEn = 0.;

    Energy_Iso_ForLocal(
        beadID,
        resi,
        r_pos0,
        &oldEn,
        &newEn,
        &old_ovlp_num,
        &old_cont_num,
        g_system.scratch.naOldOvlpNeighs,
        g_system.scratch.naOldContNeighs);

    int bondList[MAX_BONDS + 1];
    int bondNum = g_system.bSystemHasTopo ? OP_GetTopoBonds(beadID, bondList) : 0;

    bondNum = BeadList_CanTopoAngle(bondNum, bondList);

    oldEn += bondNum ? Energy_Topo_Angle_ForList(bondNum, bondList) : 0.;

    lLDub BWRos = MC_RosenbluthSampling_ForLocal_AtOld(beadID, resi, &oldEn, old_ovlp_num);

    OP_System_MoveBeadTo(beadID, r_posNew);

    newEn += g_system.nBiasPotential_Mode == -1 ? 0.f : Energy_BiasingPotential(beadID);

    Energy_Iso_ForLocal(
        beadID,
        resi,
        r_posNew,
        &newEn,
        &oldEn,
        &new_ovlp_num,
        &new_cont_num,
        g_system.scratch.naNewOvlpNeighs,
        g_system.scratch.naNewContNeighs);

    newEn += bondNum ? Energy_Topo_Angle_ForList(bondNum, bondList) : 0.;

    lLDub FWRos = MC_RosenbluthSampling_ForLocal_AtNew(beadID, resi, &yTemp, &newEn, new_ovlp_num);

    lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(FWRos, BWRos, oldEn - newEn, (lLDub)MyTemp);
    // Accept or reject this state
    if (MCProb < MHAcc)
    { // Accept this state
        if (g_system.bead_info[beadID][BEAD_FACE] != -1)
        { // Breaking old bond
            g_system.bead_info[g_system.bead_info[beadID][BEAD_FACE]][BEAD_FACE] = -1;
        }
        g_system.bead_info[beadID][BEAD_FACE] = yTemp;
        if (yTemp != -1)
        { // Making new bond
            g_system.bead_info[yTemp][BEAD_FACE] = beadID;
        }
        bAccept = 1;
        return bAccept;
    }
    else
    {
        OP_System_MoveBeadTo_Inv(beadID);
        bAccept = 0;
        return bAccept;
    }
}

/// Move_Snake - performs a biased reptation move on chainID by:
/// 1. Randomly picking which end to reptate.
/// 2. Performing sort of a local move by finding an empty spot to move to,
/// while recording
///    the Rosenbluth weights, but for the whole chain. Calculate the total
///    energy of the chain.
/// Break all the physical bonds.
/// 3. Moving the chain along and recalculating the total weights. While
/// calculating the total weights, also randomly assign a physical bond from the
/// boltzmann distribution. Calculate the new energy. Note that rather than
/// \f\$prod_{i=1}^{N}W_{i}\f$, I calculate the \f$\log_{10}\f$ so I calculate
/// the sum of smallish numbers rather than the product of large numbers.
/// 4. Perform the Metropolis-Hastings step.
/// The move is automatically rejected for molecules that are not chains.
/// The move is NOT smart enough to detect that the given chain has all the same
/// linker lengths within the chain.
/// TODO: In initialization, add a sub-routine that checks if the molecules have
/// the same linker lengths throughout.
/// \param chainID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_Snake(int chainID, float MyTemp)
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
        {
            // If chain is not linear. Reject move because slithering will not
            // work!
            bAccept = 0;
            return bAccept;
        }
    }

    // This chain is suitable to perform a slithering-snake move on.

    int i; // Loop iterator
    int resi;
    int yTemp; // Random numbers to store things
    int r_posNew[POS_MAX], r_posTmp1[POS_MAX],
        r_posTmp2[POS_MAX]; // Vectors to store positions.
    int angBead_old, angBead_new;
    lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX; // To decide if we slither forwards or backwards
    const char SnakeFwd = MCProb < 0.5f ? 1 : 0;    // 1 for Fwds, and 0 for Bkwds.

    if (SnakeFwd)
    { // Forwards slither, so lastB-1 (last bead) is anchor
        LatPos_gen_rand_wRad(r_posTmp1,
                             g_system.linker_len[lastB - 1][0]); // lastB-1 will be replaced by lastB-2
        LatPos_copy(r_posTmp2, g_system.bead_info[lastB - 1]);
        LatPos_add_wPBC(r_posNew, r_posTmp1, r_posTmp2);
        yTemp = Check_MoveBeadTo(r_posNew); // 0: there is no space, 1: there is space
        angBead_old = firstB + 1;
        angBead_new = lastB - 2;
    }
    else
    { // Backwards slither, so firstB is anchor
        LatPos_gen_rand_wRad(r_posTmp1,
                             g_system.linker_len[firstB][0]); // firstB will be replaced by firstB+1
        LatPos_copy(r_posTmp2, g_system.bead_info[firstB]);
        LatPos_add_wPBC(r_posNew, r_posTmp1, r_posTmp2);
        yTemp = Check_MoveBeadTo(r_posNew); // 0: there is no space, 1: there is space
        angBead_new = firstB + 1;
        angBead_old = lastB - 2;
    }

    if (yTemp == 0)
    { // Couldn't find a spot, so reject the damn move
        bAccept = 0;
        return bAccept;
    }

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

    lLDub BSum = 0.;
    for (i = firstB; i < lastB; i++)
    {
        Energy_Iso_ForChains(
            i,
            &oldEn,
            &newEn,
            &old_ovlp_num,
            &old_cont_num,
            g_system.scratch.naOldOvlpNeighs,
            g_system.scratch.naOldContNeighs);
        resi = g_system.bead_info[i][BEAD_TYPE];
        BSum += MC_RosenbluthSampling_ForChains_AtOld(
            i,
            resi,
            &oldEn,
            old_ovlp_num); // Calculate the Rosenbluth weights. Summation in log
                           // space means multiplication in linear space.
    }

    resi = g_system.bead_info[angBead_old][BEAD_TYPE];
    oldEn += g_energy.faEnergy[resi][resi][E_STIFF] ? Energy_Topo_Angle(angBead_old) : 0.;

    if (SnakeFwd)
    { // Slithering the chain forwards in ID-space
        OP_System_Snake_SlitherFwd(firstB, lastB, r_posNew);
    }
    else
    { // Slithering backwards in ID-space
        OP_System_Snake_SlitherBck(firstB, lastB, r_posNew);
    }

    // Have slithered the chain whichever way
    for (i = firstB; i < lastB; i++)
    { // Break the old bonds!
        yTemp = g_system.bead_info[i][BEAD_FACE];
        if (yTemp != -1)
        { // Need to break this bond
            g_system.bead_info[yTemp][BEAD_FACE] = -1;
            g_system.bead_info[i][BEAD_FACE] = -1;
        }
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
        Energy_Iso_ForChains(
            i,
            &newEn,
            &oldEn,
            &new_ovlp_num,
            &new_cont_num,
            g_system.scratch.naNewOvlpNeighs,
            g_system.scratch.naNewContNeighs);
        resi = g_system.bead_info[i][BEAD_TYPE];
        FSum += MC_RosenbluthSampling_ForChains_AtNew(i, resi, &yTemp, &newEn, new_ovlp_num);
        if (yTemp != -1)
        { // An appropriate partner has been selected. Form
          // the bonds and add the energy
            g_system.bead_info[i][BEAD_FACE] = yTemp;
            g_system.bead_info[yTemp][BEAD_FACE] = i;
            // newEn += Energy_Anisotropic_For_Chain(i);
        }
    }

    resi = g_system.bead_info[angBead_new][BEAD_TYPE];
    newEn += g_energy.faEnergy[resi][resi][E_STIFF] ? Energy_Topo_Angle(angBead_new) : 0.;

    // Doing the Metropolis-Hastings thing
    MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(FSum, BSum, oldEn - newEn, (lLDub)MyTemp);

    if (MCProb < MHAcc)
    { // Accept this state. Bonds have already been formed!
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

/// Move_Trans - performs a biased translation of chainID by:
/// 1. Seeing if there is a spot to move the entire chain in a
/// <+-L/2,+-L/2,+-L/2> random location. Move fails if no spot available.
/// 2. Calculate the Rosenbluth weights of the whole chain like Move_Snake().
/// 3. Move the chain and recalculate the total weights.
/// 4. Perform the Metropolis-Hastings step.
/// \param chainID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_Trans(int chainID, float MyTemp)
{                        // Performs a translation move with orientational bias
    int bAccept = 0;     // Used in MC steps
    int yTemp;           // Random numbers to store things
    int r_disp[POS_MAX]; // Vectors to store coordinates.

    // All moves are L/2 radius
    LatPos_gen_rand_wRad(r_disp, g_system.naBoxSize[0] / 2);

    yTemp = Check_ChainDisp(chainID, r_disp); // yTemp=0 means clash
    if (yTemp == 0)
    { // We have failed to find a good spot for this chain.
        bAccept = 0;
        return bAccept;
    }

    // We now have a chain which when moved does not overlap.

    int i, j; // Loop iterators

    int resi, resj, firstB, lastB;
    // Finding the bounds for looping over the molecule/chain
    firstB = g_system.chain_info[chainID][CHAIN_START];
    lastB = firstB + g_system.chain_info[chainID][CHAIN_LENGTH];

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

    lLDub BSum = 0.;
    for (i = firstB; i < lastB; i++)
    {
        Energy_Iso_ForChains(
            i,
            &oldEn,
            &newEn,
            &old_ovlp_num,
            &old_cont_num,
            g_system.scratch.naOldOvlpNeighs,
            g_system.scratch.naOldContNeighs);
        resi = g_system.bead_info[i][BEAD_TYPE];
        BSum += MC_RosenbluthSampling_ForChains_AtOld(i, resi, &oldEn, old_ovlp_num);
    }

    OP_System_DispChain_ForTrans(chainID,
                                 r_disp); // Moved the chain, broke bonds, and remembered stuff.

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
        Energy_Iso_ForChains(
            i,
            &newEn,
            &oldEn,
            &new_ovlp_num,
            &new_cont_num,
            g_system.scratch.naNewOvlpNeighs,
            g_system.scratch.naNewContNeighs);
        resi = g_system.bead_info[i][BEAD_TYPE];
        FSum += MC_RosenbluthSampling_ForChains_AtNew(i, resi, &yTemp, &newEn, new_ovlp_num);
        if (yTemp != -1)
        { // An appropriate partner has been selected. Form
          // the bonds and add the energy
            g_system.bead_info[i][BEAD_FACE] = yTemp;
            g_system.bead_info[yTemp][BEAD_FACE] = i;
            // newEn += Energy_Anisotropic_For_Chain(i);
        }
    }

    lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(FSum, BSum, oldEn - newEn, (lLDub)MyTemp);
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

/// Move_Clus - translates the second largest cluster of the system by:
/// 1. Performs a total clustering analysis of the system and finds the second
/// largest cluster. If there are more than 1 second largest clusters, randomly
/// pick 1.
/// 2. Try to find an empty spot in a <+-L/2,+-L/2,+-L/2> random location.
/// Calculate old energy.
/// 3. Move the cluster over.
/// Note that since the definition of a cluster is based on the network of
/// existing physical bonds, the system energy, the clusters do not change post
/// translations.
/// 4. Calculate the new energy and perform the Metropolis-Hastings step.
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_Clus_Network(float MyTemp)
{
    // Attempts to move the second largest cluster

    int bAccept = 0; // Used in MC steps, assume that move fails initially.
    int* naClusList = (int*)malloc(sizeof(int) * (g_system.tot_chains + 1));

    const int ClusSize = ClusUtil_AnisoCluster_OfSystem_SecondLargest_ForMCMove(naClusList);

    if (ClusSize < 1)
    {
        bAccept = 0;
        free(naClusList);
        return bAccept;
    }

    int r_Disp[POS_MAX];
    int yTemp;
    int i, j;
    // Radii for translation moves. All moves are L/4 radius
    // I guess moving single chains around as well is not a bad idea

    LatPos_gen_rand_wRad(r_Disp, g_system.naBoxSize[0] / 2);
    for (i = 0; i < ClusSize; i++)
    {
        yTemp = Check_ChainDisp(naClusList[i], r_Disp); // Checking for steric clash
        if (yTemp == 0)
        {
            bAccept = 0;
            free(naClusList);

            // printf("End CLUS - No space\n");
            return bAccept;
        }
    }
    // This means that there is no steric clash when the cluster is moved.

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;
    int thisChain, firstB, lastB;

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (g_system.nBiasPotential_Mode != -1)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = naClusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                oldEn += Energy_BiasingPotential(i);
            }
        }
    }

    if (g_system.bSystemHasCont || g_system.bSystemHasOvlp || g_system.bSystemHasFSol)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = naClusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                Energy_Iso_ForChains(
                    i,
                    &oldEn,
                    &newEn,
                    &old_ovlp_num,
                    &old_cont_num,
                    g_system.scratch.naOldOvlpNeighs,
                    g_system.scratch.naOldContNeighs);
            }
        }
    }

    for (j = 0; j < ClusSize; j++)
    {
        OP_System_DispChain(naClusList[j], r_Disp); // Moving the cluster properly
    }

    if (g_system.nBiasPotential_Mode != -1)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = naClusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                newEn += Energy_BiasingPotential(i);
            }
        }
    }

    if (g_system.bSystemHasCont || g_system.bSystemHasOvlp || g_system.bSystemHasFSol)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = naClusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                Energy_Iso_ForChains(
                    i,
                    &newEn,
                    &oldEn,
                    &new_ovlp_num,
                    &new_cont_num,
                    g_system.scratch.naNewOvlpNeighs,
                    g_system.scratch.naNewContNeighs);
            }
        }
    }

    const lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    const lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub)MyTemp);
    if (MCProb < MHAcc)
    {                // Accept this state
        bAccept = 1; // Accept the move
        free(naClusList);
        return bAccept;
        // printf("End CLUS - Yes\n");
    }
    else
    {
        for (i = 0; i < ClusSize; i++)
        {
            OP_System_RestoreChain(naClusList[i]); // Placing  the cluster back properly
        }
        bAccept = 0;
        free(naClusList);
        return bAccept;
        // printf("End CLUS - Failed.\n");
    }
}

/// Move_SmallClus - translates this cluster only if it is smaller than 5 total
/// molecules by
/// 1. Performs a clustering analysis on chainID. Move fails if the cluster is
/// larger than 5.
/// 2. Try to find an empty spot in a <+-L/2,+-L/2,+-L/2> random location.
/// Calculate old energy.
/// 3. Move the cluster over.
/// Note that since the definition of a cluster is based on the network of
/// existing physical bonds, the system energy, the clusters do not change post
/// translations.
/// 4. Calculate the new energy and perform the Metropolis-Hastings step.
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_SmallClus_Network(int chainID, float MyTemp)
{
    // Performs a cluster move where a given chain and it's cluster are moved.
    // No new 'bonds' are made so the move is reversible....

    int bAccept = 0; // Used in MC steps, assume that move fails initially.

    //    const int ClusSize = Clus_Network_LimitedCluster(chainID); // Looking at everything that
    //    is connected to chainID

    char* oldHashTab = calloc(g_system.tot_chains, sizeof(char));
    int* clusList = calloc(g_system.tot_chains, sizeof(int));

    const int ClusSize =
        ClusUtil_AnisoCluster_OfChain_wMaxSize(chainID, oldHashTab, clusList, g_analysis.nLimitedClusterSize);
    if ((ClusSize < 2))
    {
        bAccept = 0;
        free(oldHashTab);
        free(clusList);
        return bAccept;
    }

    int r_Disp[POS_MAX];
    int yTemp;
    int i, j;
    // Radii for translation moves. All moves are L/2 radius
    // I guess moving single chains around as well is not a bad idea

    LatPos_gen_rand_wRad(r_Disp, g_system.naBoxSize[0] / 2);
    for (i = 0; i < ClusSize; i++)
    {
        // printf("%d\n", g_analysis.cycle.naList[i]);
        yTemp = Check_ChainDisp(clusList[i], r_Disp); // Checking for steric clash
        if (yTemp == 0)
        {
            bAccept = 0;
            // printf("End CLUS - No space\n");
            free(oldHashTab);
            free(clusList);
            return bAccept;
        }
    }
    // This means that there is no steric clash when the cluster is moved.

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;
    int thisChain, firstB, lastB;

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (g_system.nBiasPotential_Mode != -1)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = clusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                oldEn += Energy_BiasingPotential(i);
            }
        }
    }

    if (g_system.bSystemHasCont || g_system.bSystemHasOvlp || g_system.bSystemHasFSol)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = clusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                Energy_Iso_ForChains(
                    i,
                    &oldEn,
                    &newEn,
                    &old_ovlp_num,
                    &old_cont_num,
                    g_system.scratch.naOldOvlpNeighs,
                    g_system.scratch.naOldContNeighs);
            }
        }
    }

    for (j = 0; j < ClusSize; j++)
    {
        OP_System_DispChain(clusList[j], r_Disp); // Moving the cluster properly
    }

    if (g_system.nBiasPotential_Mode != -1)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = clusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                newEn += Energy_BiasingPotential(i);
            }
        }
    }

    if (g_system.bSystemHasCont || g_system.bSystemHasOvlp || g_system.bSystemHasFSol)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = clusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                Energy_Iso_ForChains(
                    i,
                    &newEn,
                    &oldEn,
                    &new_ovlp_num,
                    &new_cont_num,
                    g_system.scratch.naNewOvlpNeighs,
                    g_system.scratch.naNewContNeighs);
            }
        }
    }
    const lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    const lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub)MyTemp);
    if (MCProb < MHAcc)
    {                // Accept this state
        bAccept = 1; // Accept the move
        free(oldHashTab);
        free(clusList);
        return bAccept;
        // printf("End CLUS - Yes\n");
    }
    else
    {
        for (i = 0; i < ClusSize; i++)
        {
            OP_System_RestoreChain(clusList[i]); // Placing  the cluster back properly
        }
        bAccept = 0;
        free(oldHashTab);
        free(clusList);
        return bAccept;
        // printf("End CLUS - Failed.\n");
    }
}

/// Move_DbPvt - performs a Double-Pivot move on beadID by:
/// 1. Calculate the number of possible bridges in a <+-2,+-2,+-2> around
/// beadID. The move is rejected if the molecule is not a chain.
/// 2. Pick one of these bridges, and remember which beadID'+1 will be beadID+1
/// after the move.
/// 3. Calculate the number of bridges possible around beadID'+1.
/// 4. Calculate the ratio between the number of possible bridges, and perform
/// the Metropolis-Hastings step. Note that bridging only changes the local
/// connectivity of the chains, and thus not the energy.
/// \param beadID
/// \return 1 if accepted, 0 if rejected.
int Move_DbPvt(const int beadID, const float myTemp)
{ // Performs a double-pivot move.
    /* Molecule MUST be LINEAR
  The move requires selecting a random bead, which is beadID. Then, we'll search
  the lattice in +-min(l,3) sites around beadID. Let i be the position of beadID along
  it's chain. Let i' denote same position along another chain of the same type.
  We want dist(i,i'+1) < g_system.linker_len[i] && Dist_BeadToBead (i',i+1) <
  g_system.linker_len[i']. We'll count how many candidates there are and select one
  randomly. Then we just swap beads from i to then end of the chain. Do a
  Metropolis thing, and decide. In other words, i'+1 becomes i+1, i'+2 becomes
  i+2 until N, and i+1 become i'+1 and so on.
   */
    int bAccept = 0;                                               // Move acceptance and such LEL
    const int PChainID = g_system.bead_info[beadID][BEAD_CHAINID]; // The proposed chainID
    const int PType = g_system.chain_info[PChainID][CHAIN_TYPE];   // Type of chain.

    if (g_system.nChainTypeIsLinear[PType] == 0)
    {
        // Reject the move because the chain is not linear.
        return bAccept;
    }
    // Now make sure that the proposed bead is neither the start or end of a
    // chain.
    const int PStart = g_system.chain_info[PChainID][CHAIN_START];         // Start of this chain.
    const int PEnd = PStart + g_system.chain_info[PChainID][CHAIN_LENGTH]; // LastBead+1 for this chain.
    if (PStart == beadID || PEnd - 1 == beadID)
    {
        return bAccept;
    }

    int i, j; // Loop iterators
    int candList[345] = {0};
    int filtList[345] = {0};

    const int thisBead = beadID;
    const int thisBeadChainPos = thisBead - PStart; // 0-indexed chain position.
    const int nextBead = beadID + 1;
    const int nextBeadChainPos = thisBeadChainPos + 1;

    // The length of the box is min(LARGEST_RADIUS, g_system.linker_len[beadID][1])
    const int n_rad = LARGEST_RADIUS < g_system.linker_len[thisBead][1]
                          ? LARGEST_RADIUS
                          : g_system.linker_len[thisBead][1]; // The length of the box we'll be searching.

    const int r_pos0[POS_MAX] = {
        g_system.bead_info[thisBead][0],
        g_system.bead_info[thisBead][1],
        g_system.bead_info[thisBead][2]};
    int neigh_num = NeighborSearch_AroundPoint_wRad_IgnBead(thisBead, r_pos0, n_rad, candList);

    BeadListOP_GetChainIDs(neigh_num, candList, filtList);
    int other_chains_num = BeadListOP_InvFilter_wrt_SecondList(neigh_num, candList, filtList, PChainID);

    BeadListOP_GetChainTypes(other_chains_num, candList, filtList);
    int same_chain_types = BeadListOP_Filter_wrt_SecondList(other_chains_num, candList, filtList, PType);

    BeadListOP_GetIntraChainID(same_chain_types, candList, filtList);
    int corr_chain_pos = BeadListOP_Filter_wrt_SecondList(same_chain_types, candList, filtList, nextBeadChainPos);

    const int fwd_cand_num = BeadListOP_Filter_DbPvtLinkerConFwd(corr_chain_pos, candList, nextBead);

    if (fwd_cand_num == 0) // Rejecting if no forward moves
    {
        bAccept = 0;
        return bAccept;
    }

    // For detailed balance, we need the reverse move. So we look at i', and find suitable i+1
    // candidates.

    const int nextBead_P = candList[rand() % fwd_cand_num];
    const int thisBead_P = nextBead_P - 1;
    const int PChainID_P = g_system.bead_info[thisBead_P][BEAD_CHAINID];

    const int r_pos0_P[POS_MAX] = {
        g_system.bead_info[thisBead_P][0],
        g_system.bead_info[thisBead_P][1],
        g_system.bead_info[thisBead_P][2]};

    neigh_num = NeighborSearch_AroundPoint_wRad_IgnBead(thisBead_P, r_pos0_P, n_rad, candList);

    BeadListOP_GetChainIDs(neigh_num, candList, filtList);
    other_chains_num = BeadListOP_InvFilter_wrt_SecondList(neigh_num, candList, filtList, PChainID_P);

    BeadListOP_GetChainTypes(other_chains_num, candList, filtList);
    same_chain_types = BeadListOP_Filter_wrt_SecondList(other_chains_num, candList, filtList, PType);

    BeadListOP_GetIntraChainID(same_chain_types, candList, filtList);
    corr_chain_pos = BeadListOP_Filter_wrt_SecondList(same_chain_types, candList, filtList, nextBeadChainPos);

    const int bck_cand_num = BeadListOP_Filter_DbPvtLinkerConFwd(corr_chain_pos, candList, nextBead_P);

    if (bck_cand_num == 0)
    {          // We swap the chains.
        j = 0; // Tracks the beads
        for (i = nextBead; i < PEnd; i++)
        { // Swapping from beadID+1 onwards
            OP_SwapBeads(nextBead + j,
                         nextBead_P + j); // This is pretty dumb, but easier to read UGH
            j++;
        }
        bAccept = 1;
        return bAccept;
    }
    lLDub oldEn, newEn;

    const int resi = g_system.bead_info[thisBead][BEAD_TYPE];
    oldEn = 0.;
    if (g_energy.faEnergy[resi][resi][E_STIFF])
    {
        oldEn = Energy_Topo_Angle(thisBead) + Energy_Topo_Angle(nextBead);
        oldEn += Energy_Topo_Angle(thisBead_P) + Energy_Topo_Angle(nextBead_P);
    }

    // We swap the chains.
    j = 0; // Tracks the beads
    for (i = nextBead; i < PEnd; i++)
    { // Swapping from beadID+1 onwards
        OP_SwapBeads(nextBead + j,
                     nextBead_P + j); // This is pretty dumb, but easier to read UGH
        j++;
    }

    newEn = 0.;
    if (g_energy.faEnergy[resi][resi][E_STIFF])
    {
        newEn = Energy_Topo_Angle(thisBead) + Energy_Topo_Angle(nextBead);
        newEn += Energy_Topo_Angle(thisBead_P) + Energy_Topo_Angle(nextBead_P);
    }

    const lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    const lLDub MHAcc_W = (lLDub)(fwd_cand_num) / (lLDub)(bck_cand_num);
    if (MCProb < MHAcc_W * expl((oldEn - newEn) / myTemp))
    {
        bAccept = 1;
        return bAccept;
    }
    else
    { // Rejecting, so we swap back.
        // We swap the chains.
        j = 0; // Tracks the beads
        for (i = nextBead; i < PEnd; i++)
        { // Swapping from beadID+1 onwards
            OP_SwapBeads(nextBead + j,
                         nextBead_P + j); // This is pretty dumb, but easier to read UGH
            j++;
        }
        bAccept = 0;
        return bAccept;
    } //
}

/// Move_CoLocal - move thisBeadID and it's physical bond partner to a new
/// location. Move fails if thisBeadID has no partner, and if no space is found.
/// A standard unbiased translation of two beads where a Metropolis-Hastings
/// step is performed at the end.
/// \param thisBeadID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_CoLocal(int thisBeadID, float MyTemp)
{
    /*
  Translate a bead and its partner in tandem. If no partner, reject move.
  */

    int bAccept = 0; // Used in MC steps

    if (g_system.bead_info[thisBeadID][BEAD_FACE] == -1)
    {
        bAccept = 0;
        return bAccept;
    }
    int r_disp[POS_MAX]; // Random translation vector.
    int r_posNew1[POS_MAX], r_posNew2[POS_MAX];
    int r_pos1[POS_MAX], r_pos2[POS_MAX];
    const int otherBeadID = g_system.bead_info[thisBeadID][BEAD_FACE];
    int yTemp;

    LatPos_copy(r_pos1, g_system.bead_info[thisBeadID]);
    LatPos_copy(r_pos2, g_system.bead_info[otherBeadID]);

    LatPos_gen_rand_wRad(r_disp, 2);
    LatPos_add_wPBC(r_posNew1, r_disp, r_pos1);
    LatPos_add_wPBC(r_posNew2, r_disp, r_pos2);

    yTemp = Check_MoveBeadTo(r_posNew1) * Check_MoveBeadTo(r_posNew2);
    if (yTemp == 1)
    { // This means we found an empty lattice site. So let's
      // check if the linkers are okay.
        yTemp = Check_LinkerConstraint(thisBeadID, r_posNew1) * Check_LinkerConstraint(otherBeadID, r_posNew2);
    }
    if (yTemp == 0)
    { // Steric-clash or bad location for the beads.
        bAccept = 0;
        return bAccept;
    }

    lLDub newEn = 0.;
    lLDub oldEn = 0.;
    if (g_system.nBiasPotential_Mode != -1)
    {
        oldEn += Energy_BiasingPotential(thisBeadID);
        oldEn += Energy_BiasingPotential(otherBeadID);
    }

    int ovlp_num, cont_num;
    Energy_Iso_ForCoLocal(
        thisBeadID,
        otherBeadID,
        r_pos1,
        &oldEn,
        &newEn,
        &ovlp_num,
        &cont_num,
        g_system.scratch.naOldOvlpNeighs,
        g_system.scratch.naOldContNeighs);
    Energy_Iso_ForCoLocal(
        otherBeadID,
        thisBeadID,
        r_pos2,
        &oldEn,
        &newEn,
        &ovlp_num,
        &cont_num,
        g_system.scratch.naOldOvlpNeighs,
        g_system.scratch.naOldContNeighs);

    int bondList[2 * (MAX_BONDS + 1)];
    int bond_tmpList[MAX_BONDS + 1];
    int bondNum, bond_tmpNum;

    if (g_system.bSystemHasTopo)
    {
        bond_tmpNum = OP_GetTopoBonds(thisBeadID, bond_tmpList);
        bondNum = BeadList_AppendBeads(0, bondList, bond_tmpList, bond_tmpNum);
        bond_tmpNum = OP_GetTopoBonds(otherBeadID, bond_tmpList);
        bondNum = BeadList_AppendBeads(bondNum, bondList, bond_tmpList, bond_tmpNum);
        qsort(bondList, bondNum, sizeof(int), UtilFunc_CompareInts);
        bondNum = ListOP_UniqueElementsOfSortedList_Int(bondNum, bondList);
        bondNum = BeadList_CanTopoAngle(bondNum, bondList);
    }
    else
    {
        bondNum = 0;
    }

    oldEn += bondNum ? Energy_Topo_Angle_ForList(bondNum, bondList) : 0.;

    OP_System_MoveBeadTo(thisBeadID, r_posNew1);
    OP_System_MoveBeadTo(otherBeadID, r_posNew2);

    if (g_system.nBiasPotential_Mode != -1)
    {
        newEn += Energy_BiasingPotential(thisBeadID);
        newEn += Energy_BiasingPotential(otherBeadID);
    }

    Energy_Iso_ForCoLocal(
        thisBeadID,
        otherBeadID,
        r_posNew1,
        &newEn,
        &oldEn,
        &ovlp_num,
        &cont_num,
        g_system.scratch.naOldOvlpNeighs,
        g_system.scratch.naOldContNeighs);
    Energy_Iso_ForCoLocal(
        otherBeadID,
        thisBeadID,
        r_posNew2,
        &newEn,
        &oldEn,
        &ovlp_num,
        &cont_num,
        g_system.scratch.naOldOvlpNeighs,
        g_system.scratch.naOldContNeighs);

    newEn += bondNum ? Energy_Topo_Angle_ForList(bondNum, bondList) : 0.;

    lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub)MyTemp);
    if (MCProb < MHAcc)
    { // Accept this state
        bAccept = 1;
        // printf("Accepted!\n");
        return bAccept;
    }
    else
    {
        OP_System_MoveBeadTo_Inv(thisBeadID);
        OP_System_MoveBeadTo_Inv(otherBeadID);
        bAccept = 0;
        // printf("Rejected move; undoing!\n");
        return bAccept;
    }
}

/// Move_MultiLocal - performs a biased set of local moves where beadID and
/// every covalently bonded bead of beadID is moved by:
/// 1. Finding an empty spot for all the beads involved (in a <+-2,+-2,+-2>).
/// This is a little crude because I put the beads on the lattice into the new
/// proposed spot because the structure constraint sub-routine is designed as
/// such.
/// 2. If there is space for all of the beads. Calculate the energies and the
/// total Rosenbluth weights.
/// 3. Move the beads to new location and recalculate the energies and weights.
/// Perform a Metropolis-Hastings step.
/// \param beadID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_MultiLocal(int beadID, float MyTemp)
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

    lLDub BSum = 0.;
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
            g_system.scratch.naOldOvlpNeighs,
            g_system.scratch.naOldContNeighs);
        BSum += MC_RosenbluthSampling_ForLists_AtOld(i, beadNum, beadsList, &oldEn, old_ovlp_num);
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
    OP_Beads_BreakBondsInList(beadNum, beadsList);

    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            newEn += Energy_BiasingPotential(tmpBead);
        }
    }

    lLDub FSum = 0.;
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
            g_system.scratch.naNewOvlpNeighs,
            g_system.scratch.naNewContNeighs);
        FSum += MC_RosenbluthSampling_ForLists_AtNew(i, beadNum, beadsList, &yTemp, &newEn, new_ovlp_num);
        if (yTemp != -1)
        { // An appropriate partner has been selected. Form
          // the bonds and add the energy
            tmpBead = beadsList[i];
#if DEBUG_BUILD
            if (g_system.bead_info[tmpBead][BEAD_FACE] != -1)
            {
                printf("HOW?!\n");
                exit(1);
            }
#endif
            g_system.bead_info[tmpBead][BEAD_FACE] = yTemp;
            g_system.bead_info[yTemp][BEAD_FACE] = tmpBead;
            // newEn += Energy_Anisotropic_For_List(tmpBead, beadNum, beadsList);
        }
    }

    newEn += bondAngNum ? Energy_Topo_Angle_ForList(bondAngNum, bondAngList) : 0.;

    const lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    const lLDub MHAcc = OP_GenMHValue(FSum, BSum, oldEn - newEn, (lLDub)MyTemp);

    if (MCProb < MHAcc)
    { // Accept the move. Remember that the bonds were
      // assigned above
        bAccept = 1;
        return bAccept;
    }
    else
    {
        OP_Beads_BreakBondsInList(beadNum, beadsList);
        OP_Lattice_EmptySitesForListOfBeads(beadNum, beadsList);
        OP_Beads_CopyBeadsInListFromOld(beadNum, beadsList);
        OP_Beads_RestoreBondsInList(beadNum, beadsList);
        OP_Lattice_PlaceBeadsInList(beadNum, beadsList);
        bAccept = 0;
        return bAccept;
    }
}

/// Move_Pivot - performs a biased pivot move on chainID.
/// Move fails if chainID is branched. The process is:
/// 1. Pick a random bead within chainID that shall act as the pivot around
/// which a rotation operation shall be performed. The operation is performed on
/// the shorter part of the chain. If the selected bead is an end, the move is
/// rejected.
/// 2. Pick a rotation operation. The list of operation is a rotation by 90,
/// 180, 270 or 360 (null) degrees around the X, Y, or Z axis.
/// 3. Calculate the Rosenbluth weights and energies, and check if there is
/// steric clash after the proposed move. If there is a clash, reject the move.
/// 4. Rotate the beads, and recalculate the weights, propose new bonds and
/// recalculate the energies.
/// 5. Perform the Metropolis-Hastings step.
/// \param chainID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_Pivot(int chainID, float MyTemp)
{
    // Performs a pivot move on chainID
    /* Randomly pick a bead in this chain (anchorBead), and perform a symmetry
     * operation on the beads after the anchorBead. Note that if anchorBead is
     * the second-to-last, or last, bead, the move is rejected. Local and
     * slither moves would be better for those. The protein must be linear and
     * the move is rejected if the chain is branched. It is assumed that a
     * linear chain has been passed!
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
        yTemp = Check_MoveBeadTo(g_system.scratch.naTempR);
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
    int resi;

    lLDub BSum = 0.;
    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        Energy_Iso_ForRange(
            tmpBead,
            smBead,
            lgBead,
            &oldEn,
            &newEn,
            &old_ovlp_num,
            &old_cont_num,
            g_system.scratch.naOldOvlpNeighs,
            g_system.scratch.naOldContNeighs);
        resi = g_system.bead_info[tmpBead][BEAD_TYPE];
        BSum += MC_RosenbluthSampling_ForRange_AtOld(tmpBead, resi, smBead, lgBead, &oldEn, old_ovlp_num);
    }

    resi = g_system.bead_info[anchorBead][BEAD_TYPE];
    oldEn += g_energy.faEnergy[resi][resi][E_STIFF] ? Energy_Topo_Angle(anchorBead) : 0.;

    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        OP_Rotation(PivotM, tmpBead, anchorPos);
        OP_System_MoveBeadTo(tmpBead, g_system.scratch.naTempR);
    }

    OP_Beads_BreakBondsInList(beadNum, beadsList);

    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            newEn += Energy_BiasingPotential(tmpBead);
        }
    }

    int new_ovlp_num, new_cont_num;
    int resj;

    lLDub FSum = 0.;
    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        resi = g_system.bead_info[tmpBead][BEAD_TYPE];
        Energy_Iso_ForRange(
            tmpBead,
            smBead,
            lgBead,
            &newEn,
            &oldEn,
            &new_ovlp_num,
            &new_cont_num,
            g_system.scratch.naNewOvlpNeighs,
            g_system.scratch.naNewContNeighs);
        FSum += MC_RosenbluthSampling_ForRange_AtNew(tmpBead, resi, smBead, lgBead, &yTemp, &newEn, new_ovlp_num);
        if (yTemp != -1)
        {
            resj = g_system.bead_info[yTemp][BEAD_TYPE];
            g_system.bead_info[yTemp][BEAD_FACE] = tmpBead;
            g_system.bead_info[tmpBead][BEAD_FACE] = yTemp;
            // newEn += Energy_Anisotropic_For_Range(tmpBead, smBead, lgBead);
        }
    }

    resi = g_system.bead_info[anchorBead][BEAD_TYPE];
    newEn += g_energy.faEnergy[resi][resi][E_STIFF] ? Energy_Topo_Angle(anchorBead) : 0.;

    const lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(FSum, BSum, oldEn - newEn, (lLDub)MyTemp);
    if (MCProb < MHAcc)
    { // Accept the move. Remember that the bonds were
      // assigned above!
        bAccept = 1;
        return bAccept;
    }
    else
    { // Rejecting move
        OP_Beads_BreakBondsInList(beadNum, beadsList);
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            OP_System_MoveBeadTo_Inv(tmpBead);
        }
        bAccept = 0;
        return bAccept;
    }
}

/// Move_BranchedRot - performs the branched analog for the pivot move where the
/// central bead of the branched molecule is considered the anchor, and the
/// whole molecule is rotated with the biased sampling, and re-assigning of
/// bonds. The move is rejected if chainID is linear.
/// \param chainID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_BranchedRot(int chainID, float MyTemp)
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

    int i, j;
    int xTemp, yTemp;
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
        yTemp = Check_MoveBeadTo(g_system.scratch.naTempR);
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
    int resi;

    lLDub BSum = 0.;
    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        Energy_Iso_ForRange(
            tmpBead,
            smBead,
            lgBead,
            &oldEn,
            &newEn,
            &old_ovlp_num,
            &old_cont_num,
            g_system.scratch.naOldOvlpNeighs,
            g_system.scratch.naOldContNeighs);
        resi = g_system.bead_info[tmpBead][BEAD_TYPE];
        BSum += MC_RosenbluthSampling_ForRange_AtOld(tmpBead, resi, smBead, lgBead, &oldEn, old_ovlp_num);
    }

    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        OP_Rotation(PivotM, tmpBead, anchorPos);
        OP_System_MoveBeadTo(tmpBead, g_system.scratch.naTempR);
    }

    OP_Beads_BreakBondsInList(beadNum, beadsList);

    if (g_system.nBiasPotential_Mode != -1)
    {
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            newEn += Energy_BiasingPotential(tmpBead);
        }
    }

    int new_ovlp_num, new_cont_num;
    int resj;

    lLDub FSum = 0.;
    for (i = 0; i < beadNum; i++)
    {
        tmpBead = beadsList[i];
        resi = g_system.bead_info[tmpBead][BEAD_TYPE];
        Energy_Iso_ForRange(
            tmpBead,
            smBead,
            lgBead,
            &newEn,
            &oldEn,
            &new_ovlp_num,
            &new_cont_num,
            g_system.scratch.naNewOvlpNeighs,
            g_system.scratch.naNewContNeighs);
        FSum += MC_RosenbluthSampling_ForRange_AtNew(tmpBead, resi, smBead, lgBead, &yTemp, &newEn, new_ovlp_num);
        if (yTemp != -1)
        {
            resj = g_system.bead_info[yTemp][BEAD_TYPE];
            g_system.bead_info[yTemp][BEAD_FACE] = tmpBead;
            g_system.bead_info[tmpBead][BEAD_FACE] = yTemp;
            // newEn += Energy_Anisotropic_For_Range(tmpBead, smBead, lgBead);
        }
    }

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
        OP_Beads_BreakBondsInList(beadNum, beadsList);
        for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            OP_System_MoveBeadTo_Inv(tmpBead);
        }
        bAccept = 0;
        return bAccept;
    }
}

/// Move_SmallClus_Proximity
/// \param chainID
/// \param myTemp
/// \return
int Move_SmallClus_Proximity(const int chainID, const float myTemp)
{
    // Performs a cluster move where a given chain and it's cluster are moved.

    int bAccept = 0; // Used in MC steps, assume that move fails initially.

    int* oldClusList = calloc(g_system.tot_chains, sizeof(int));
    char* oldHashTab = calloc(g_system.tot_chains, sizeof(char));
    const int ClusSize =
        ClusUtil_OvlpCluster_OfChain_wMaxSize(chainID, oldHashTab, oldClusList, g_analysis.nLimitedClusterSize);

    if (ClusSize < 2)
    {
        bAccept = 0;
        free(oldClusList);
        free(oldHashTab);
        return bAccept;
    }
    //    printf("\n%d\n", ClusSize);

    // Radii for translation moves. All moves are L/2 radius
    int r_Disp[POS_MAX];
    int yTemp;
    int i, j;
    LatPos_gen_rand_wRad(r_Disp, g_system.naBoxSize[0] / 2);

    for (i = 0; i < ClusSize; i++)
    {
        yTemp = Check_ChainDisp(oldClusList[i], r_Disp); // Checking for steric clash
        if (yTemp == 0)
        {
            bAccept = 0;
            free(oldClusList);
            free(oldHashTab);
            return bAccept;
        }
    }
    // No  clash

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;
    int thisChain, firstB, lastB;

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (g_system.nBiasPotential_Mode != -1)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = oldClusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                oldEn += Energy_BiasingPotential(i);
            }
        }
    }

    // If no CONT interactions, then just checking the cluster is enough. If we have a different
    // cluster, we immediately reject. Therefore, the clusters are moved as rigid bodies and will
    // have the same OVLP energies before and after. If we have CONT interactions, we could have
    // different energies. Therefore, we should calculate the energies.
    if (g_system.bSystemHasCont)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = oldClusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                Energy_Iso_ForChains(
                    i,
                    &oldEn,
                    &newEn,
                    &old_ovlp_num,
                    &old_cont_num,
                    g_system.scratch.naOldOvlpNeighs,
                    g_system.scratch.naOldContNeighs);
            }
        }
    }

    for (i = 0; i < ClusSize; i++)
    {
        OP_System_DispChain(oldClusList[i], r_Disp); // Moving the cluster properly
    }
    // Recalculating cluster to see if we have the same cluster or not. If not, we reject.

    const int ClusCheckN = ClusUtil_OvlpCluster_OfChain_CheckForSame(oldHashTab, oldClusList, ClusSize);
    free(oldHashTab);

    if (ClusCheckN == -1)
    {
        for (i = 0; i < ClusSize; i++)
        {
            OP_System_RestoreChain(oldClusList[i]); // Placing  the cluster back properly
        }
        bAccept = 0;
        free(oldClusList);
        return bAccept;
    }

    if (g_system.nBiasPotential_Mode != -1)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = oldClusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                newEn += Energy_BiasingPotential(i);
            }
        }
    }

    if (g_system.bSystemHasCont)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = oldClusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                Energy_Iso_ForChains(
                    i,
                    &newEn,
                    &oldEn,
                    &new_ovlp_num,
                    &new_cont_num,
                    g_system.scratch.naNewOvlpNeighs,
                    g_system.scratch.naNewContNeighs);
            }
        }
    }

    const lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    const lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub)myTemp);
    if (MCProb < MHAcc)
    {                // Accept this state
        bAccept = 1; // Accept the move
        free(oldClusList);
        return bAccept;
        // printf("End CLUS - Yes\n");
    }
    else
    {
        for (i = 0; i < ClusSize; i++)
        {
            OP_System_RestoreChain(oldClusList[i]); // Placing  the cluster back properly
        }
        bAccept = 0;
        free(oldClusList);
        return bAccept;
        // printf("End CLUS - Failed.\n");
    }
}

/// Move_Clus_Proximity - Attempts to move the second largest cluster in the system
///  \param myTemp
///  \return
int Move_Clus_Proximity(const float myTemp)
{
    // Performs a cluster move where a given chain and it's cluster are moved.
    // No new 'bonds' are made so the move is reversible.... So we have to check if we get the same
    // cluster back

    int bAccept = 0; // Used in MC steps, assume that move fails initially.

    int* naClusList = (int*)malloc(sizeof(int) * (g_system.tot_chains + 1));
    const int ClusSize = ClusUtil_OvlpCluster_OfSystem_SecondLargest_ForMCMove(naClusList);
    if (ClusSize < 2)
    {
        bAccept = 0;
        free(naClusList);
        return bAccept;
    }

    // Radii for translation moves. All moves are L/2 radius
    int r_Disp[POS_MAX];
    int yTemp;
    int i, j;
    LatPos_gen_rand_wRad(r_Disp, g_system.naBoxSize[0] / 2);

    for (i = 0; i < ClusSize; i++)
    {
        yTemp = Check_ChainDisp(naClusList[i], r_Disp); // Checking for steric clash
        if (yTemp == 0)
        {
            bAccept = 0;
            free(naClusList);
            return bAccept;
        }
    }
    // No  clash

    char* oldHashTab = calloc(g_system.tot_chains, sizeof(char));

    for (i = 0; i < ClusSize; i++)
    {
        oldHashTab[naClusList[i]] = 1;
    }

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;
    int thisChain, firstB, lastB;

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (g_system.nBiasPotential_Mode != -1)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = naClusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                oldEn += Energy_BiasingPotential(i);
            }
        }
    }

    // If no CONT interactions, then just checking the cluster is enough. If we have a different
    // cluster, we immediately reject. Therefore, the clusters are moved as rigid bodies and will
    // have the same OVLP energies before and after. If we have CONT interactions, we could have
    // different energies. Therefore, we should calculate the energies.
    if (g_system.bSystemHasCont)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = naClusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                Energy_Iso_ForChains(
                    i,
                    &oldEn,
                    &newEn,
                    &old_ovlp_num,
                    &old_cont_num,
                    g_system.scratch.naOldOvlpNeighs,
                    g_system.scratch.naOldContNeighs);
            }
        }
    }

    for (i = 0; i < ClusSize; i++)
    {
        OP_System_DispChain(naClusList[i], r_Disp); // Moving the cluster properly
    }
    // Recalculating cluster to see if we have the same cluster or not. If not, we reject.
    const int ClusCheckN = ClusUtil_OvlpCluster_OfChain_CheckForSame(oldHashTab, naClusList, ClusSize);
    free(oldHashTab);

    if (ClusCheckN == -1)
    {
        for (i = 0; i < ClusSize; i++)
        {
            OP_System_RestoreChain(naClusList[i]); // Placing  the cluster back properly
        }
        bAccept = 0;
        free(naClusList);
        return bAccept;
    }

    if (g_system.nBiasPotential_Mode != -1)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = naClusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                newEn += Energy_BiasingPotential(i);
            }
        }
    }

    if (g_system.bSystemHasCont)
    {
        for (j = 0; j < ClusSize; j++)
        {
            thisChain = naClusList[j];
            firstB = g_system.chain_info[thisChain][CHAIN_START];
            lastB = firstB + g_system.chain_info[thisChain][CHAIN_LENGTH];
            for (i = firstB; i < lastB; i++)
            {
                Energy_Iso_ForChains(
                    i,
                    &newEn,
                    &oldEn,
                    &new_ovlp_num,
                    &new_cont_num,
                    g_system.scratch.naNewOvlpNeighs,
                    g_system.scratch.naNewContNeighs);
            }
        }
    }

    const lLDub MCProb = (lLDub)rand() / (lLDub)RAND_MAX;
    const lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub)myTemp);
    if (MCProb < MHAcc)
    {                // Accept this state
        bAccept = 1; // Accept the move
        free(naClusList);
        return bAccept;
        // printf("End CLUS - Yes\n");
    }
    else
    {
        for (i = 0; i < ClusSize; i++)
        {
            OP_System_RestoreChain(naClusList[i]); // Placing  the cluster back properly
        }
        bAccept = 0;
        free(naClusList);
        return bAccept;
        // printf("End CLUS - Failed.\n");
    }
}
