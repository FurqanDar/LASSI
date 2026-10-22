#include "structure.h"

/// NeighborSearch_ForOvlp: Given the beadID, which should be at startVec, we calculate all the
/// neighbors within the +-1 cube. We return the number of neighbors. \param beadID \param startVec
/// \param neighList
/// \return Number of neighbors.
int NeighborSearch_ForOvlp(int const beadID, const int* startVec, int* neighList)
{
    int neigh_num = 0;

    int tmpBead;
    int r_disp[POS_MAX] = {0};
    int r_chck[POS_MAX] = {0};
    const int nRad = 1;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
    {
        LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_X);
        for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
        {
            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Y);
            for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
            {
                LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Z);
                tmpBead = g_system.naTotLattice[Lat_Ind_FromVec(r_chck)];
                if (tmpBead != -1 && tmpBead != beadID)
                {
                    neighList[neigh_num++] = tmpBead;
                }
            }
        }
    }
    neighList[neigh_num] = -1;
    return neigh_num;
}

/// NeighborSearch_ForOvlp: Given the beadID, which should be at startVec, we calculate all the
/// neighbors within the +-1 cube. We return the number of neighbors. \param beadID \param startVec
/// \param neighList
/// \return Number of neighbors.
int NeighborSearch_ForCont(int const beadID, const int* startVec, int* contList, int* ovlpList,
                           int* ovlpNum)
{
    int neigh_num = 0;

    int tmpBead;
    int r_disp[POS_MAX] = {0};
    int r_chck[POS_MAX] = {0};
    const int nRad = LARGEST_RADIUS;

    *ovlpNum = 0;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
    {
        LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_X);
        for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
        {
            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Y);
            for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
            {
                LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Z);
                tmpBead = g_system.naTotLattice[Lat_Ind_FromVec(r_chck)];
                if (tmpBead != -1 && tmpBead != beadID)
                {
                    if ((abs(r_disp[0]) <= 1) && (abs(r_disp[1]) <= 1) && (abs(r_disp[2]) <= 1))
                    {
                        ovlpList[*ovlpNum] = tmpBead;
                        *ovlpNum = *ovlpNum + 1;
                    }
                    contList[neigh_num++] = tmpBead;
                }
            }
        }
    }
    contList[neigh_num] = -1;
    ovlpList[*ovlpNum] = -1;
    return neigh_num;
}

/// NeighborSearch_EmptySitesAround_wRad: Given the starting point, and radius, we calculate how
/// many sites around this point are empty. \param startVec \param nRad \return Number of empty
/// lattice-sites around startVec
int NeighborSearch_EmptySitesAround_wRad(const int* startVec, const int nRad)
{
    int empty_num = 0;
    int tmpBead;
    int r_disp[POS_MAX] = {0};
    int r_search[POS_MAX] = {0};

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
    {
        for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
        {
            for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
            {
                LatPos_add_wPBC(r_search, startVec, r_disp);
                tmpBead = g_system.naTotLattice[Lat_Ind_FromVec(r_search)];
                if (tmpBead == -1)
                {
                    empty_num++;
                }
            }
        }
    }
    return empty_num;
}

/// NeighborSearch_SolSitesAround: number of solvation sites around this point. Finds the number of
/// empty lattice sites in the 26 closest sites. \param startVec \return Empty sites around this
/// point.
int NeighborSearch_SolSitesAround(const int* startVec)
{
    return NeighborSearch_EmptySitesAround_wRad(startVec, 1);
}

int NeighborSearch_AroundPoint_UptoIndex(const int beadID, const int* startVec, const int upTO,
                                         int* neighList)
{
    int i;
    int neigh_num = 0;
    int tmpBead;

    int r_search[POS_MAX] = {0};

    int r_all[100][POS_MAX] = {0};

    for (i = 0; i < upTO; ++i)
    {
        LatPos_add_wPBC(r_search, startVec, r_all[i]);
        tmpBead = g_system.naTotLattice[Lat_Ind_FromVec(r_search)];
        if (tmpBead != -1)
        {
            neighList[neigh_num++] = tmpBead;
        }
    }

    return neigh_num;
}

int NeighborSearch_AroundPoint_wRad_IgnBead(const int beadID, const int* startVec, const int nRad,
                                            int* neighList)
{
    int neigh_num = 0;
    int tmpBead;
    int r_disp[POS_MAX] = {0};
    int r_search[POS_MAX] = {0};

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
    {
        for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
        {
            for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
            {
                LatPos_add_wPBC(r_search, startVec, r_disp);
                tmpBead = g_system.naTotLattice[Lat_Ind_FromVec(r_search)];
                if (tmpBead != -1 && tmpBead != beadID)
                {
                    neighList[neigh_num++] = tmpBead;
                }
            }
        }
    }
    return neigh_num;
}

int NeighborSearch_AroundPoint_wRad_wDists(const int beadID, const int* startVec, const int nRad,
                                           int* neighList, float* distList)
{
    int neigh_num = 0;
    int tmpBead;
    int r_disp[POS_MAX] = {0};

    int r_search[POS_MAX] = {0};

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
    {
        for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
        {
            for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
            {
                LatPos_add_wPBC(r_search, startVec, r_disp);
                tmpBead = g_system.naTotLattice[Lat_Ind_FromVec(r_search)];
                if (tmpBead != -1 && tmpBead != beadID)
                {
                    neighList[neigh_num] = tmpBead;
                    distList[neigh_num] = Dist_Vec3n(r_disp);
                    neigh_num++;
                }
            }
        }
    }
    return neigh_num;
}

/// NeighborSearch_ForOvlp: Given the beadID, which should be at startVec, we calculate all the
/// neighbors within the +-1 cube. We return the number of neighbors. \param beadID \param startVec
/// \param neighList
/// \return Number of neighbors.
int NeighborSearch_ForCluster_Ovlp(int const beadID, const int* startVec, int* neighList)
{
    int neigh_num = 0;

    int tmpBead;
    int r_disp[POS_MAX] = {0};
    int r_chck[POS_MAX] = {0};
    const int nRad = 1;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
    {
        LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_X);
        for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
        {
            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Y);
            for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
            {
                LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Z);
                tmpBead = g_system.naTotLattice[Lat_Ind_FromVec(r_chck)];
                if (tmpBead != -1)
                {
                    neighList[neigh_num] = tmpBead;
                    neigh_num++;
                }
            }
        }
    }
    neighList[neigh_num] = -1;
    return neigh_num;
}
