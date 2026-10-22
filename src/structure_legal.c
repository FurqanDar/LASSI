#include "structure.h"

/*
 * Move-facing legality / linker-constraint helpers.
 * Pure cores take explicit distance + linker length; wrappers bind g_system.
 */

/// LinkerConstraint_DistOK - whether Euclidean distance xDis is within the
/// scaled linker length (LINKER_RSCALE * linker_len). Pure: no globals.
/// \return 1 if OK, 0 if broken (same convention as Check_LinkerConstraint).
int LinkerConstraint_DistOK(const float xDis, const int linker_len)
{
    return (xDis > LINKER_RSCALE * (float)linker_len) ? 0 : 1;
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
    idx = 0;
    bondPartner = g_system.topo_info[beadID][idx]; // Initializing the two.
    while (idx < MAX_BONDS && g_system.topo_info[beadID][idx] != -1)
    { // Keep going till we run out of partners
        bondPartner = g_system.topo_info[beadID][idx];
        if (!LinkerConstraint_DistOK(
                Dist_PointToPoint(g_system.bead_info[bondPartner], tmpR),
                g_system.linker_len[beadID][idx]))
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
/// \param listSize
/// \param beadList
/// \return
int Check_LinkerConstraints_ForBeadList(const int listSize, const int* beadList)
{
    int i, j, bead1, bead2;
    int dumBonds[MAX_BONDS + 1];
    int dumNum;
    float xDis;
    for (i = 0; i < listSize; i++)
    {
        bead1 = beadList[i];
        dumNum = OP_GetTopoBonds(bead1, dumBonds);
        for (j = 0; j < dumNum - 1; j++)
        {
            bead2 = dumBonds[j + 1];
            xDis = Dist_BeadToBead(bead1, bead2);
            if (!LinkerConstraint_DistOK(xDis, g_system.linker_len[bead1][j]))
            {
                return 0;
            }
        }
    }

    return 1;
}

/// Check_BeadID_InList: Iterate over the list of beads to see if the bead is in the list.
/// Pure: no globals. If so, return 1. Otherwise return 0.
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
