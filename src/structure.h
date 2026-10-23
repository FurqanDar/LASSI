#ifndef _STRUCTURE_H_ // include guard
#define _STRUCTURE_H_

#include "global.h"

/*
 * Public structure / lattice / list / analysis-index API.
 * Family-only decls live in structure_internal.h (structure*.c only).
 */

/* Lattice index + min-image distance (see structure_lattice.c).
 * Pure cores take explicit box; wrappers bind g_system. */

int Lat_Ind_FromCoords_Box(const int i, const int j, const int k, const int box_x, const int box_y);

int Lat_Ind_FromVec(const int* xArr);

int Lat_Ind_OfBead(const int beadID);

float Dist_PointToPoint_Box(const int* f1, const int* f2, const int* box);

float Dist_BeadToBead(const int bead1, const int bead2);

float Dist_BeadToPoint(const int bead1, const int* f1);

float Dist_BeadToPoint_Double(const int beadID, const lDub* f1);

float Dist_Vec3n(const int* f1);

int Dist_VecMagSq(const int* f1);

void LatPos_copy(int* copy_vec, const int* in_vec);

void LatPos_add_wPBC(int* outVec, const int* firVec, const int* secVec);

void BeadPos_sub_wPBC(int* outVec, const int* firVec, const int* secVec);

void LatPos_gen_rand_wRad(int outVec[POS_MAX], const int nRadius);

float Vec3n_CosTheta(const int* v1, const int* v2);

/* Legality / linker constraints (see structure_legal.c). */

int LinkerConstraint_DistOK(const float xDis, const int linker_len);

int Check_LinkerConstraint(const int beadID, const int* tmpR);

int Check_LinkerConstraints_ForBeadList(const int listSize, const int* beadList);

int Check_BeadID_InList(const int thisBeadID, const int listSize, const int beadList[MAX_BONDS + 1]);

void PerformRuntimeSanityChecks(long nGen, int run_cycle);

/* Analysis entrypoints / indices used outside analysis.c (impl in analysis.c). */

void RDF_ComponentWise_Avg(void);

int RDFArr_Index(const int run_cycle, const int rdf_comp, const int x_pos);

int RadDenArr_Index(const int run_cycle, const int rad_comp, const int x_pos);

int MolClusArr_Index(const int run_cycle, const int chain_type, const int clus_size);

void RadialDensityAnalysis_Perform_Analysis(void);

void GyrTensor_GyrRad_Avg(void);

/* Neighbor search (see structure_neighbor.c). */

int NeighborSearch_ForOvlp(int const beadID, const int* startVec, int* neighList);

int NeighborSearch_ForCont(int const beadID, const int* startVec, int* contList, int* ovlpList, int* ovlpNum);

int NeighborSearch_AroundPoint_wRad_IgnBead(const int beadID, const int* startVec, const int nRad, int* neighList);

int NeighborSearch_ForCluster_Ovlp(int const beadID, const int* startVec, int* neighList);

/* Topology / lattice utils (see structure.c). */

int OP_GetTopoBonds(const int beadID, int* dum_list);

int LatticeUtil_GetNeighBeads_AtPos(const int* restrict r_pos0, int* restrict nBeadsList);

/* Bead / chain / list helpers (see structure_lists.c). */

void BeadListOP_GetChainIDs(const int beadNum, const int* beadList, int* chainList);

void BeadListOP_GetChainTypes(const int beadNum, const int* beadList, int* chainList);

int BeadListOP_Filter_wrt_SecondList(const int beadNum, int* beadList, const int* propList, const int prop_val);

int BeadListOP_InvFilter_wrt_SecondList(const int beadNum, int* beadList, const int* propList, const int prop_val);

void BeadListOP_GetIntraChainID(const int beadNum, const int* beadList, int* chainList);

int BeadListOP_Filter_DbPvtLinkerConFwd(const int beadNum, int* beadList, const int thisBead);

int ListOP_UniqueElementsOfSortedList_Int(const int size, int* sorted_list);

int BeadList_AppendBeads(const int old_size, int* old_list, const int* app_list, const int app_size);

int UtilFunc_CompareInts(const void* a, const void* b);

int BeadList_CanTopoAngle(const int size, int* beadList);

int ChainListOP_AddUniqueChains_wSize(
    const int clusSize,
    int* clusList,
    const int newChainNums,
    const int* chainList,
    const int maxClus);

int ChainListOP_AddUniqueChains_wSize_Check(
    const int clusSize,
    int* clusList,
    const int newChainNums,
    const int* chainList,
    const int maxClus,
    const int* oldList);

void ChainListOP_GetChainTypes(int* const naTypesList, const int* const naChainIDs, const int nListSize);

int ListOP_GetMaxVal_Int(const int nListSize, const int* const naList);

int ListOP_GetRandomIndexForVal_Int(const int nVal, const int nListSize, const int* const naList);

int ListOP_Get2ndLargestVal_Int(const int nListSize, const int* const naInputList);

void ListOP_GenHistFromCounts_Int(
    int* const restrict naOutFreqHist,
    const int* const restrict naInputList,
    const int nSize);

#endif // _STRUCTURE_H_
