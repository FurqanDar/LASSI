#include "structure_internal.h"

/* File-local helpers (not in structure.h). */
static float BeadPos_CosThetaOfBeads(const int bead1, const int bead2);
static void LatticeUtil_GetOvlpLattIndecies(const int* restrict const r_pos0, int* restrict numList);
static void
LatticeUtil_GetLattVals_FromList(const int* restrict const nIndexList, int* restrict nLatValsList, const int listSize);
static int
LatticeUtil_GetBeadIDs_FromList(const int* restrict const nLatValsList, int* restrict nBeadsList, const int nListSize);

/// OP_GetTopoBonds - including beadID, get all the beads that are covalently bonded to beadID.
/// \param beadID
/// \param dum_list
/// \return Number_of_bonds. Can be used to loop over list which includes beadID.
int OP_GetTopoBonds(const int beadID, int* dum_list)
{
    int top_it = 0;
    int curID = beadID;

    while (curID != -1 && top_it < MAX_BONDS + 1)
    {
        dum_list[top_it] = curID;
        curID = g_system.topo_info[beadID][top_it++];
    }

    return top_it;
}

static float BeadPos_CosThetaOfBeads(const int bead1, const int bead2)
{
    const int r_pos1[POS_MAX] = {
        g_system.bead_info[bead1][0],
        g_system.bead_info[bead1][1],
        g_system.bead_info[bead1][2]};

    const int r_pos2[POS_MAX] = {
        g_system.bead_info[bead2][0],
        g_system.bead_info[bead2][1],
        g_system.bead_info[bead2][2]};

    return Vec3n_CosTheta(r_pos1, r_pos2);
}

/// LatticeUtil_GetOvlpLattIndecies: Given the starting location r_pos0, we store all of the
/// g_system.naTotLattice indecies corresponding to the +-1 neighbors. We store all the indecies in
/// numList. NOTE: This list _will_ include the point at r_pos0 as well. \param r_pos0 \param
/// numList - MUST be at least 27 = 3^3 elements long. Undefined behavior if not that long.
static void LatticeUtil_GetOvlpLattIndecies(const int* restrict const r_pos0, int* restrict numList)
{
    int i = 0;
    int r_disp[POS_MAX] = {0};
    int r_chck[POS_MAX] = {0};
    const int nRad = 1;
    int tmpBead;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
    {
        LatPos_add_wPBC_ofComp(r_chck, r_pos0, r_disp, POS_X);
        for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
        {
            LatPos_add_wPBC_ofComp(r_chck, r_pos0, r_disp, POS_Y);
            for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
            {
                LatPos_add_wPBC_ofComp(r_chck, r_pos0, r_disp, POS_Z);
                tmpBead = Lat_Ind_FromVec(r_chck);
                numList[i++] = tmpBead;
            }
        }
    }
}

/// LatticeUtil_GetLattVals_FromList - loop over nIndexList and store the g_system.naTotLattice
/// values in nLatValsList. Note that this will include -1 values as well, or empty sites. \param
/// nIndexList \param nLatValsList - Must be at least 27 elements long.
static void
LatticeUtil_GetLattVals_FromList(const int* restrict const nIndexList, int* restrict nLatValsList, const int listSize)
{
    int i;
    for (i = 0; i < listSize; ++i)
    {
        nLatValsList[i] = g_system.naTotLattice[nIndexList[i]];
    }
}

/// LatticeUtil_GetBeadIDs_FromList - Given the list of Lattice values, we get all the values that
/// are not -1. Or, the we get the values that correspond to non-empty lattice sites. \param
/// nLatValsList \param nBeadsList \return
static int
LatticeUtil_GetBeadIDs_FromList(const int* restrict const nLatValsList, int* restrict nBeadsList, const int nListSize)
{
    int nBeadsNum = 0;
    int tmpBead;
    int i;

    for (i = 0; i < nListSize; ++i)
    {
        tmpBead = nLatValsList[i];
        if (tmpBead != -1)
        {
            nBeadsList[nBeadsNum++] = tmpBead;
        }
    }
    return nBeadsNum;
}

/// LatticeUtil_GetNeighBeads_AtPos - given the starting position
/// \param r_pos0
/// \param nBeadsList
/// \return
int LatticeUtil_GetNeighBeads_AtPos(const int* restrict const r_pos0, int* restrict nBeadsList)
{
    int nLatInd[CLUS_CONTACT_NEIGHS];
    LatticeUtil_GetOvlpLattIndecies(r_pos0, nLatInd);

    int nLatVals[CLUS_CONTACT_NEIGHS];
    LatticeUtil_GetLattVals_FromList(nLatInd, nLatVals, CLUS_CONTACT_NEIGHS);

    return LatticeUtil_GetBeadIDs_FromList(nLatVals, nBeadsList, CLUS_CONTACT_NEIGHS);
}
