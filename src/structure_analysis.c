#include "structure.h"
#include "cluster.h"
#include "energy.h"
#include "mcmove.h"
#include "print.h"


/// GyrTensor_GyrRad_Avg - calculates the total radius of gyration of the system, while not being smart about the
/// periodic boundaries. This is used as a proxy to detect phase separation, but is a relic of the old formalism. The
/// RDF should be used in general. Although this can be used without the need for a non-interacting prior.
void GyrTensor_GyrRad_Avg(void)
{
    /*
    Only calculates the diagonals of the gyration tensor, and calculates the sum of the
    diagonals. Remember that Rg^2 = Tr(GyrTen) so we only need to calculate the diagonals, and then then sum.
    I shall borrow most of the code from above, and so read GyrTensor_ClusterSpecific for what's happening here.
    */
    int i, j; // Loop indecies
    lDub tot_COM[POS_MAX] = {0.};

    Calc_SystemCenterOfMass(tot_COM);
    int dumArg[POS_MAX] = {0};
    for (i = 0; i < 7; i++)
        { // Initializing to 0
            faGyrTensor_glb[i] = 0.f;
        }

    for (i = 0; i < tot_beads_glb; i++)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    dumArg[j] = abs(bead_info_glb[i][j] - (int) tot_COM[j]);
                    dumArg[j] = dumArg[j] > naBoxSize_glb[j] / 2 ? naBoxSize_glb[j] - dumArg[j] : dumArg[j];
                    faGyrTensor_glb[j] += (float) (dumArg[j] * dumArg[j]);
                }
        }

    // Adding to the total faSysGyrRad_glb to be averaged at the end.
    faSysGyrRad_glb += sqrtf((faGyrTensor_glb[0] + faGyrTensor_glb[1] + faGyrTensor_glb[2]) / (float) tot_beads_glb);
    nTotGyrRadCounter_glb++; // Remembering that we have calculated the radius; for final averaging.
}


/// RDF_ComponentIndex - 1D index for the symmetric g_{ij} matrix given i and j.
/// \param i
/// \param j
/// \return The index of the array
/// The way it is set up, the indexing goes through the diagonal and then 0-1, 0-2, ... 0-N, 1-2, ... 1-N and so on
int RDF_ComponentIndex(const int i, const int j)
{
    if (i == j)
        {
            return i + 1;
        }
    else if (i < j)
        {
            return nBeadTypes_glb + j - (i * (3 + i - 2 * nBeadTypes_glb)) / 2;
        }
    else
        {
            return nBeadTypes_glb + i - (j * (3 + j - 2 * nBeadTypes_glb)) / 2;
        }
}


/// RDFArr_Index - 1D index for ldaTOTRDF_Arr_glb which is used to globally store the different RDFs
/// \param run_cycle
/// \param rdf_comp
/// \param x_pos
/// \return 1D index for the totalRDFArray
int RDFArr_Index(const int run_cycle, const int rdf_comp, const int x_pos)
{
    return x_pos + nRDF_TotBins_glb * (rdf_comp + nRDF_TotComps_glb * run_cycle);
}


int RadDen_ComponentIndex(const int i, const int j)
{
    if (i < 0)
        {
            return j;
        }
    else
        {
            return tot_chain_types_glb + j + tot_chain_types_glb * i;
        }
}


int RadDenArr_Index(const int run_cycle, const int rad_comp, const int x_pos)
{
    return x_pos + nRDF_TotBins_glb * (rad_comp + nRadDen_TotComps_glb * run_cycle);
}


int MolClusArr_Index(const int run_cycle, const int chain_type, const int clus_size)
{
    return clus_size + (int) tot_chains_glb * (chain_type + (int) tot_chain_types_glb * run_cycle);
}


/// RDF_ComponentWise_Avg - calculates the pair-distribution of the system where every bead acts as the center
/// of a radial histogram of pairs. Note that dr = 1/4 lattice units.
void RDF_ComponentWise_Avg(void)
{
    /*
    Calculates the RDF and adds it all up so that it can be averaged out at the end of the run.
    */
    float x; // For distance
    int i, j, k;
    int resi, resj;
    int myBin = 0;
    int array_pos;

    // Calculating where, and how many, pairs exist
    for (i = 0; i < tot_beads_glb; i++)
        {
            resi = bead_info_glb[i][BEAD_TYPE];
            for (j = i + 1; j < tot_beads_glb; j++)
                {
                    resj = bead_info_glb[j][BEAD_TYPE];
                    x    = Dist_BeadToBead(i, j);
                    // Note that Dist_BeadToBead(i,j) automatically ensures no distance is greater than (L/2)*sqrt(3)
                    myBin = (int) floor(4. * x);                      // I am assuming for now that dr=1/4
                    ldaRDF_Arr_glb[RDFArr_Index(0, 0, myBin)] += 2.0; // Adding a pair to that bin
                    array_pos = RDF_ComponentIndex(resi, resj);
                    ldaRDF_Arr_glb[RDFArr_Index(0, array_pos, myBin)] += 2.0;
                }
        }
    nTotRDFCounter_glb++;
}


///
/// \param tmpR
void Calc_SystemCenterOfMass(lDub* const tmpR)
{

    int i, j;                     // Iterators
    lDub tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lDub) naBoxSize_glb[j];
        }

    for (i = 0; i < tot_beads_glb; i++)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    dumArg = dumConst[j] * (lDub) bead_info_glb[i][j];
                    zeta[j] += sin(dumArg);
                    xi[j] += cos(dumArg);
                }
        }

    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lDub) tot_beads_glb;
                    zeta[j] /= (lDub) tot_beads_glb;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
            else
                {
                    tot_COM[j] = 0.;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            tmpR[j] = tot_COM[j];
        }
}


///
/// \param naCOM_out
/// \param cluster_size
/// \param naClusList_in
void Calc_CenterOfMass_OfCluster(lDub* const naCOM_out, const int cluster_size, const int* const restrict naClusList_in)
{
    // This version measures the COM of a cluster of size c
    //  cluster size, given the molecule ID's in naList_glb.
    // The COM from this is not necessarily the COM of the system as a whole.
    int thisMol, i, j, k;         // Iterators
    int fB, lB;                   // Keep track of first and last beads of a given molecule
    lDub tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored
    int bead_total_now    = 0;

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lLDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lLDub) naBoxSize_glb[j];
        }

    for (k = 0; k < cluster_size; k++)
        {
            thisMol = naClusList_in[k];
            fB      = chain_info_glb[thisMol][CHAIN_START];
            lB      = fB + chain_info_glb[thisMol][CHAIN_LENGTH];
            for (i = fB; i < lB; i++)
                {
                    bead_total_now++;
                    for (j = 0; j < POS_MAX; j++)
                        {
                            dumArg = dumConst[j] * (lDub) bead_info_glb[i][j];
                            zeta[j] += sin(dumArg);
                            xi[j] += cos(dumArg);
                        }
                }
        }
    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lLDub) bead_total_now;
                    zeta[j] /= (lLDub) bead_total_now;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            naCOM_out[j] = tot_COM[j];
        }
}


///
/// \param tmpR
/// \param thisType
void Calc_SystemCenterOfMass_OfMolType(lDub* tmpR, const int thisType)
{
    // This version measures the COM of only type thisType
    // The COM from this is not necessarily the COM of the system as a whole.
    int thisMol, i, j, k;         // Iterators
    int fB, lB;                   // Keep track of first and last beads of a given molecule
    lDub tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored
    int bead_total_now    = 0;

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lDub) naBoxSize_glb[j];
        }

    for (k = 0; k < tot_chains_glb; k++)
        {
            thisMol = chain_info_glb[k][CHAIN_TYPE];
            if (thisMol != thisType)
                {
                    continue;
                }
            fB = chain_info_glb[k][CHAIN_START];
            lB = fB + chain_info_glb[k][CHAIN_LENGTH];
            for (i = fB; i < lB; i++)
                {
                    bead_total_now++;
                    for (j = 0; j < POS_MAX; j++)
                        {
                            dumArg = dumConst[j] * (lDub) bead_info_glb[i][j];
                            zeta[j] += sin(dumArg);
                            xi[j] += cos(dumArg);
                        }
                }
        }
    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lDub) bead_total_now;
                    zeta[j] /= (lDub) bead_total_now;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
        }

    /*for (j=0; j<POS_MAX; j++){
        tmpR[j] = ((int) tot_COM[j]);
    }*/
    for (j = 0; j < POS_MAX; j++)
        {
            tmpR[j] = tot_COM[j];
        }
}


///
/// \param tmpR
/// \param thisType
void Calc_SystemCenterOfMass_WithoutMolType(lDub* tmpR, const int thisType)
{
    // This version measures the COM of everything except type thisType
    // The COM from this is not necessarily the COM of the system as a whole.
    int thisMol, i, j, k;         // Iterators
    int fB, lB;                   // Keep track of first and last beads of a given molecule
    lDub tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored
    int bead_total_now    = 0;

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lDub) naBoxSize_glb[j];
        }

    for (k = 0; k < tot_chains_glb; k++)
        {
            thisMol = chain_info_glb[k][CHAIN_TYPE];
            if (thisMol == thisType)
                {
                    continue;
                }
            fB = chain_info_glb[k][CHAIN_START];
            lB = fB + chain_info_glb[k][CHAIN_LENGTH];
            for (i = fB; i < lB; i++)
                {
                    bead_total_now++;
                    for (j = 0; j < POS_MAX; j++)
                        {
                            dumArg = dumConst[j] * (lDub) bead_info_glb[i][j];
                            zeta[j] += sin(dumArg);
                            xi[j] += cos(dumArg);
                        }
                }
        }
    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lDub) bead_total_now;
                    zeta[j] /= (lDub) bead_total_now;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            tmpR[j] = tot_COM[j];
        }
}


void RadDen_Avg_MolTypeWise_FromSysCen(void)
{

    int i;        // Iterator for loop
    int thisType; // Tracks the type of the chain
    lDub sysCOM[POS_MAX] = {0.};
    int myBin;
    float xDis = 0.; // Tracks the distance between the COM and the specific bead

    Calc_SystemCenterOfMass(sysCOM);
    for (i = 0; i < tot_beads_glb; i++)
        {
            thisType = bead_info_glb[i][BEAD_CHAINID];
            thisType = chain_info_glb[thisType][CHAIN_TYPE];
            xDis     = Dist_BeadToPoint_Double(i, sysCOM);
            myBin    = (int) floor(4. * xDis);
            ldaRadDen_Arr_glb[RadDenArr_Index(0, thisType, myBin)] += 1.0;
        }
    nTotRadDenCounter_glb++;
}


void RadDen_Avg_MolTypeWise_FromMolTypeCen_Old_CorrectVersion(void)
{

    int i, j;     // Iterators for loop
    int thisType; // Tracks the type of the chain
    int thisComp; // Tracks which component of ldRadDen
    lDub typeCOM[POS_MAX] = {0};
    int COM_int[POS_MAX]  = {0.};
    int myBin;
    float xDis   = 0.; // Tracks the distance between the COM and the specific bead
    int cur_type = 0;

    Calc_SystemCenterOfMass(typeCOM);
    for (j = 0; j < POS_MAX; j++)
        {
            COM_int[j] = (int) typeCOM[j];
        }
    for (i = 0; i < tot_beads_glb; i++)
        {
            thisType = bead_info_glb[i][BEAD_CHAINID];
            thisType = chain_info_glb[thisType][CHAIN_TYPE];
            thisComp = RadDen_ComponentIndex(-1, thisType);
            xDis     = Dist_BeadToPoint(i, COM_int);
            myBin    = (int) (4. * xDis);
            ldaRadDen_Arr_glb[RadDenArr_Index(0, thisComp, myBin)] += 1.0;
        }

    for (cur_type = 0; cur_type < tot_chain_types_glb; cur_type++)
        { // Go through each molType's center
            Calc_SystemCenterOfMass_OfMolType(typeCOM, cur_type);
            for (j = 0; j < POS_MAX; j++)
                {
                    COM_int[j] = (int) typeCOM[j];
                }
            for (i = 0; i < tot_beads_glb; i++)
                {
                    thisType = bead_info_glb[i][BEAD_CHAINID];
                    thisType = chain_info_glb[thisType][CHAIN_TYPE];
                    thisComp = RadDen_ComponentIndex(cur_type, thisType);
                    xDis     = Dist_BeadToPoint(i, COM_int);
                    myBin    = (int) (4. * xDis);
                    ldaRadDen_Arr_glb[RadDenArr_Index(0, thisComp, myBin)] += 1.0;
                }
        }

    /*int tmpBead;
    int cur_POS[POS_MAX] = {0.};
    int cur_DIS[POS_MAX] = {0.};
    int radRange;
    radRange = naBoxSize_glb[0];
    thisComp = RadDen_ComponentIndex(2,2);
    for(j=0; j<POS_MAX; j++){
        COM_int[j] = 0;
    }
    for (i = 0; i < radRange; i++) {
        //cur_POS[POS_X] = (COM_int[POS_X] + i) % naBoxSize_glb[POS_X];
        cur_DIS[POS_X] = i > naBoxSize_glb[0]/2 ? naBoxSize_glb[0] - i: i;
        cur_DIS[POS_X] *= cur_DIS[POS_X];
        for (j = 0; j < radRange; j++) {
            //cur_POS[POS_Y] = (COM_int[POS_Y] + j ) % naBoxSize_glb[POS_Y];
            cur_DIS[POS_Y] = j > naBoxSize_glb[0]/2 ? naBoxSize_glb[0] - j: j;
            cur_DIS[POS_Y] *= cur_DIS[POS_Y];
            for (k = 0; k < radRange; k++) {
                //cur_POS[POS_Z] = (COM_int[POS_Z] + k) % naBoxSize_glb[POS_Z];
                cur_DIS[POS_Z] = k > naBoxSize_glb[0]/2 ? naBoxSize_glb[0] - k: k;
                cur_DIS[POS_Z] *= cur_DIS[POS_Z];
                xDis = sqrtf((float)(cur_DIS[POS_X]+cur_DIS[POS_Y]+cur_DIS[POS_Z] ));
                myBin = (int) (4.*xDis);
                ldaRadDen_Arr_glb[RadDenArr_Index(0, thisComp, myBin)] += 1.0;
            }
        }
    }
    */

    nTotRadDenCounter_glb++;
}


void RadialDensityAnalysis_Perform_Analysis(void)
{

    RadDenHistUtil_ForSystem_FromLargestClusterOfMolTypes();
    RadDenHistUtil_ForSystem_FromCenterOfMassOfMolTypes();

    nTotRadDenCounter_glb++;
}


void RadDenHistUtil_ForSystem_FromCenterOfMassOfMolTypes(void)
{

    int i, j, k;
    int naCOM_r[POS_MAX] = {0};
    int nBin;
    float fDis = 0.f; // Tracks the distance between the COM and the specific bead
    int nCurrentType;
    int thisType, thisComp;
    lDub LaTypeCOM_r[POS_MAX];
    int fB, lB;

    for (nCurrentType = 0; nCurrentType < tot_chain_types_glb + 1; nCurrentType++)
        {
            if (nCurrentType == 0)
                {
                    Calc_SystemCenterOfMass(LaTypeCOM_r);
                }
            else
                {
                    Calc_SystemCenterOfMass_OfMolType(LaTypeCOM_r, nCurrentType - 1);
                }
            for (j = 0; j < POS_MAX; j++)
                {
                    naCOM_r[j] = (int) LaTypeCOM_r[j];
                }

            for (k = 0; k < tot_chains_glb; k++)
                {
                    fB       = chain_info_glb[k][CHAIN_START];
                    lB       = fB + chain_info_glb[k][CHAIN_LENGTH];
                    thisType = chain_info_glb[k][CHAIN_TYPE];
                    thisComp = nRadDen_CompShift_glb + RadDen_ComponentIndex(nCurrentType - 1, thisType);
                    for (i = fB; i < lB; i++)
                        {
                            fDis = Dist_BeadToPoint(i, naCOM_r);
                            nBin = (int) (4.f * fDis);
                            ldaRadDen_Arr_glb[RadDenArr_Index(0, thisComp, nBin)] += 1.0;
                        }
                }
        }
}


void RadDenHistUtil_ForSystem_FromLargestClusterOfMolTypes(void)
{

    int i, j, k;        // Iterators for loop
    int thisType;       // Tracks the type of the chain
    int nRadDenMolComp; // Tracks which component of ldRadDen
    int thisChain;
    lDub LaTypeCOM_r[POS_MAX] = {0.};
    int naCOM_r[POS_MAX]      = {0};
    int nRadBin;
    float fDis       = 0.f; // Tracks the distance between the COM and the specific bead
    int nCurrentType = 0;
    int fB, lB;
    int nMolWiseClusSize = 0;

    int* const naMolWiseClusIDs   = malloc((tot_chain_types_glb + 2) * sizeof(int));
    int* const naTmpClusChainList = malloc((tot_chains_glb + 1) * sizeof(int));
    int* const naFullClusList     = malloc((tot_chains_glb + 1) * sizeof(int));
    int* const naFullClusSizes    = calloc((tot_chains_glb + 1), sizeof(int));
    int* const naFullClusCumSizes = calloc((tot_chains_glb + 1), sizeof(int));

    const int nMolWiseClusNum = ClusUtil_OfSystem_MolWise_GetLargestClusters(naMolWiseClusIDs, naFullClusList,
                                                                             naFullClusSizes, naFullClusCumSizes);

    for (nCurrentType = 0; nCurrentType < tot_chain_types_glb + 1; nCurrentType++)
        {

            nMolWiseClusSize = ClusUtil_GetCluster_FromFullClusAndCumSizes(
                naMolWiseClusIDs[nCurrentType], naTmpClusChainList, naFullClusList, naFullClusCumSizes, naFullClusSizes,
                nMolWiseClusNum);

            Calc_CenterOfMass_OfCluster(LaTypeCOM_r, nMolWiseClusSize, naTmpClusChainList);

            for (j = 0; j < POS_MAX; j++)
                {
                    naCOM_r[j] = (int) LaTypeCOM_r[j];
                }

            for (k = 0; k < nMolWiseClusSize; k++)
                {
                    thisChain      = naTmpClusChainList[k];
                    fB             = chain_info_glb[thisChain][CHAIN_START];
                    lB             = fB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    thisType       = chain_info_glb[thisChain][CHAIN_TYPE];
                    nRadDenMolComp = RadDen_ComponentIndex(nCurrentType - 1, thisType);
                    for (i = fB; i < lB; i++)
                        {
                            fDis    = Dist_BeadToPoint(i, naCOM_r);
                            nRadBin = (int) (4.f * fDis);
                            ldaRadDen_Arr_glb[RadDenArr_Index(0, nRadDenMolComp, nRadBin)] += 1.0;
                        }
                }
        }

    free(naTmpClusChainList);
    free(naFullClusSizes);
    free(naFullClusList);
    free(naFullClusCumSizes);
    free(naMolWiseClusIDs);
}

