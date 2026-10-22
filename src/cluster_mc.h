#ifndef _CLUSTER_MC_H_ // include guard
#define _CLUSTER_MC_H_

#include "global.h"

/*
 * Move-facing cluster queries used by MC cluster moves.
 *
 * Implemented in cluster.c. Kept as a small neutral header so the moves layer
 * can call these helpers without #include "cluster.h" (analysis-layer API).
 * See docs/architecture.md.
 */

int ClusUtil_OvlpCluster_OfChain_wMaxSize(
    const int chainID,
    char* restrict caTotClusTable,
    int* restrict clusList,
    int const maxSize);

int ClusUtil_OvlpCluster_OfChain_CheckForSame(
    const char* const caTotClusTable,
    const int* const clusList,
    int const clusSize);

int ClusUtil_AnisoCluster_OfChain_wMaxSize(
    const int chainID,
    char* restrict caTotClusTable,
    int* restrict clusList,
    int const maxSize);

int ClusUtil_OvlpCluster_OfSystem_SecondLargest_ForMCMove(int* naOutClusList);

int ClusUtil_AnisoCluster_OfSystem_SecondLargest_ForMCMove(int* naOutClusList);

#endif // _CLUSTER_MC_H_
