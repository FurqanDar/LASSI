#include "structure.h"
#include "cluster.h"
#include "energy.h"
#include "mcmove.h"
#include "print.h"


/// BeadIDListOP_GetChainIDs: Given beadList of size beadNum, we record the chainIDs of every bead into chainList.
/// Make sure chainList is as big as beadList.
/// \param beadNum
/// \param beadList
/// \param chainList
void BeadListOP_GetChainIDs(const int beadNum, const int* restrict beadList, int* restrict chainList)
{
    int i, tmp_bead, dum_chain;

    for (i = 0; i < beadNum; ++i)
        {
            tmp_bead     = beadList[i];
            dum_chain    = g_system.bead_info[tmp_bead][BEAD_CHAINID];
            chainList[i] = dum_chain;
        }

    chainList[i] = -1;
}


/// BeadIDListOP_GetChainTypes: Given beadList of size beadNum, we record the chain_types of every bead into
/// chainList.
/// Make sure chainList is as big as beadList.
/// \param beadNum
/// \param beadList
/// \param chainList
void BeadListOP_GetChainTypes(const int beadNum, const int* beadList, int* chainList)
{
    int i, tmp_bead, dum_chain;

    for (i = 0; i < beadNum; ++i)
        {
            tmp_bead     = beadList[i];
            dum_chain    = g_system.bead_info[tmp_bead][BEAD_CHAINID];
            chainList[i] = g_system.chain_info[dum_chain][CHAIN_TYPE];
        }

    chainList[i] = -1;
}


/// BeadListOP_Filter_wrt_SecondList. We overwrite beadList such that whichever index propList[i] == prop_val, we
/// add to beadList, and ignore all others. Should act like a masking operation over beadList, where all masked
/// elements are removed.
/// \param beadNum
/// \param beadList
/// \param propList
/// \param prop_val
/// \return newSize - size of overwritten beadList.
int BeadListOP_Filter_wrt_SecondList(const int beadNum, int* beadList, const int* propList, const int prop_val)
{
    int newSize = 0;
    int i;
    for (i = 0; i < beadNum; i++)
        {
            if (propList[i] == prop_val)
                {
                    beadList[newSize] = beadList[i];
                    newSize++;
                }
        }
    beadList[newSize] = -1;
    return newSize;
}


/// BeadListOP_InvFilter_wrt_SecondList. We overwrite beadList such that whichever index propList[i] != prop_val, we
/// add to beadList, and ignore all others. Should act like a masking operation over beadList, where all masked
/// elements are removed.
/// Acts as the logical not for Filter_wrt_SecondList
/// \param beadNum
/// \param beadList
/// \param propList
/// \param prop_val
/// \return newSize - size of overwritten beadList.
int BeadListOP_InvFilter_wrt_SecondList(const int beadNum, int* beadList, const int* propList, const int prop_val)
{
    int newSize = 0;
    int i;
    for (i = 0; i < beadNum; i++)
        {
            if (propList[i] != prop_val)
                {
                    beadList[newSize] = beadList[i];
                    newSize++;
                }
        }
    beadList[newSize] = -1;
    return newSize;
}


/// BeadListOP_GetIntraChainID. For the beadList of size beadNum, we store the intra-chain ID of each bead.
/// tmpBeadID - chainStartID
/// \param beadNum
/// \param beadList
/// \param chainList
void BeadListOP_GetIntraChainID(const int beadNum, const int* beadList, int* chainList)
{
    int i, tmp_bead, dum_chain, dum_len;

    for (i = 0; i < beadNum; ++i)
        {
            tmp_bead     = beadList[i];
            dum_chain    = g_system.bead_info[tmp_bead][BEAD_CHAINID];
            dum_len      = g_system.chain_info[dum_chain][CHAIN_START];
            chainList[i] = tmp_bead - dum_len;
        }
    chainList[i] = -1;
}


/// BeadListOP_Filter_DbPvtLinkerConFwd. For every bead (tmpBead) in beadList, we calculate the
/// distance between thisBead and tmpBead-1. If the distance satisfies linker constraints,
/// we add into beadList.
/// This is for DbPvt move specifically. d(i+1, i') < linker[thisBead][0]
/// \param beadNum
/// \param beadList
/// \param thisBead
int BeadListOP_Filter_DbPvtLinkerConFwd(const int beadNum, int* beadList, const int thisBead)
{
    int newSize = 0;
    int i, tmpBead;
    const float linker_cons = (float) g_system.linker_len[thisBead][0] * LINKER_RSCALE;

    float xDis;

    for (i = 0; i < beadNum; ++i)
        {
            tmpBead = beadList[i];
            xDis    = Dist_BeadToBead(thisBead, tmpBead - 1);

            if (xDis < linker_cons)
                {
                    beadList[newSize] = tmpBead;
                    newSize++;
                }
        }

    return newSize;
}


/// BeadListOP_Filter_DbPvtLinkerConBck. For every bead (tmpBead) in beadList, we calculate the
/// distance between thisBead and tmpBead-1. If the distance satisfies linker constraints,
/// we add into beadList.
/// This is for DbPvt move specifically. d(i+1, i') < linker[thisBead][0]
/// \param beadNum
/// \param beadList
/// \param thisBead
int BeadListOP_Filter_DbPvtLinkerConBck(const int beadNum, int* beadList, const int thisBead)
{
    int newSize = 0;
    int i, tmpBead;
    const float linker_cons = (float) g_system.linker_len[thisBead][1] * LINKER_RSCALE;

    float xDis = 0.f;

    for (i = 0; i < beadNum; ++i)
        {
            tmpBead = beadList[i];
            xDis    = Dist_BeadToBead(thisBead, tmpBead);

            if (xDis < linker_cons)
                {
                    beadList[newSize] = tmpBead;
                    newSize++;
                }
        }

    return newSize;
}


/// ListOP_UniqueElementsOfSortedList_Int: Overwrite the provided _sorted_ list of integers with only the unique
/// elements, and return the size of the new array. The overwritten list is also sorted. The implementation only works
/// on sorted lists. \param dum_list: Sorted list of beadID's (or integers) \return number of unique elements.
int ListOP_UniqueElementsOfSortedList_Int(const int size, int* sorted_list)
{
    int* first  = sorted_list;
    int* second = first;
    int newSize = 1;
    while (++second != (sorted_list + size))
        {
            if (*second != *first)
                {
                    *(++first) = *second;
                    newSize++;
                }
        }
    return newSize;
}


/// BeadList_AppendBeads
/// Note that old_list should be at least as big as old_size + app_size. No bounds checking is done since
/// this dumb implementation has no idea about the size of the array.
/// \param old_size Number of elements in old_list.
/// \param old_list
/// \param app_list
/// \param app_size Number of elements of app_list to be appended to old_list.
/// \return old_size + app_size. Size of the new array.
int BeadList_AppendBeads(const int old_size, int* old_list, const int* app_list, const int app_size)
{
    int i;
    for (i = 0; i < app_size; i++)
        {
            old_list[i + old_size] = app_list[i];
        }

    return old_size + app_size;
}


/// UtilFunc_CompareInts: A simple function that is used in CSTDLIB qsort for integers.
/// \param a
/// \param b
/// \return
int UtilFunc_CompareInts(const void* a, const void* b)
{
    return (*(int*) a - *(int*) b);
}


/// BeadList_CanTopoAngle: For a list of beadIDs, we only keep the ones that have TopoAngle energies.
/// \param size
/// \param beadList
/// \return An overwritten list containing only beads that have a TopoAngle energy
int BeadList_CanTopoAngle(const int size, int* beadList)
{

    int i, tmpBead, resi;
    int newSize = 0;
    for (i = 0; i < size; i++)
        {
            tmpBead = beadList[i];
            resi    = g_system.bead_info[tmpBead][BEAD_TYPE];
            if (g_energy.faEnergy[resi][resi][E_STIFF])
                {
                    beadList[newSize++] = tmpBead;
                }
        }

    return newSize;
}


/// ChainListOP_AddUniqueChains_wSize: Add any chains that are not in clusList from chainList. This function assumes
/// that chainList contains only unique elements, and that clusList has the smallest chainID as the first element, and
/// the largest chainID the last element (at clusList[clusSize-1]). \param clusSize \param clusList \param newChainNums
/// \param chainList
int ChainListOP_AddUniqueChains_wSize(const int clusSize, int* clusList, const int newChainNums, const int* chainList,
                                      const int maxClus)
{
    int i;
    int curChain;
    int newClusSize = clusSize;

    for (i = 0; i < newChainNums; i++)
        {
            curChain    = chainList[i];
            newClusSize = ClusListOP_AddChainID(newClusSize, clusList, curChain);
            if (newClusSize >= maxClus)
                {
                    return newClusSize;
                }
        }

    return newClusSize;
}


/// ChainListOP_AddUniqueChains_wSize: Add any chains that are not in clusList from chainList. This function assumes
/// that chainList contains only unique elements, and that clusList has the smallest chainID as the first element, and
/// the largest chainID the last element (at clusList[clusSize-1]).
/// \param clusSize
/// \param clusList
/// \param newChainNums
/// \param chainList
int ChainListOP_AddUniqueChains_wSize_Check(const int clusSize, int* clusList, const int newChainNums,
                                            const int* chainList, const int maxClus, const int* oldList)
{
    int i;
    int curChain;
    int newClusSize = clusSize;

    for (i = 0; i < newChainNums; i++)
        {
            curChain    = chainList[i];
            newClusSize = ClusListOP_AddChainID(newClusSize, clusList, curChain);
            if (newClusSize >= maxClus)
                {
                    return newClusSize;
                }
            if (clusList[newClusSize - 1] != oldList[newClusSize - 1])
                {
                    return maxClus;
                }
        }

    return newClusSize;
}


/// ChainListOP_GetChainTypes - Given a list that contains chainIDs, we get the type of each of the chains.
/// \param naTypesList Array where the chain-types will be recorded.
/// \param naChainIDs Array that contains all the chainIDs.
/// \param nListSize Number of chains in the naChainIDs array
void ChainListOP_GetChainTypes(int* const naTypesList, const int* const naChainIDs, const int nListSize)
{
    int i;
    int tmpChain;
    for (i = 0; i < nListSize; i++)
        {
            tmpChain       = naChainIDs[i];
            naTypesList[i] = g_system.chain_info[tmpChain][CHAIN_TYPE];
        }
}


/// ClusListOP_AddChainID - Given a cluster-chainIDs list containing clusSize chains, we try to add the proposed
/// chainID. Assumes that clusList[0] is the smallest, and clusList[clusSize-1] is the largest.
/// - If the chain is smaller than the smallest, we add it to become the first element of the chain, because it is now
/// the smallest chainID. The largest chainID moves up one, and the old smallest is placed behind the largest.
/// - Else, if the chain is larger than the largest, we add it infront of the old largest.
/// - Else, we loop over the clusList to see if the proposed chainID is unique.
/// \param clusSize
/// \param clusList
/// \param chainID
/// \return
int ClusListOP_AddChainID(const int clusSize, int* clusList, const int chainID)
{

    int i;
    for (i = 0; i < clusSize; i++)
        {
            if (clusList[i] == chainID)
                {
                    return clusSize;
                }
        }

    clusList[clusSize] = chainID;

    return clusSize + 1;
}


/// ClusListOP_AddIfUniqueChainID - Loop over the clusList. If none of the chains match, we insert the chainID behind
/// last chain in clusList. So we move the last chain up by one.
/// \param clusSize
/// \param clusList
/// \param chainID
/// \return clusSize+1 if the chain was added, clusSize if not.
int ClusListOP_AddIfUniqueChainID(const int clusSize, int* clusList, const int chainID)
{
    int i;
    for (i = 0; i < clusSize; i++)
        {
            if (clusList[i] == chainID)
                {
                    return clusSize;
                }
        }
    clusList[clusSize]     = clusList[clusSize - 1];
    clusList[clusSize - 1] = chainID;
    return clusSize + 1;
}


/// ListOP_GetMaxVal_Int - Given an integer array naList, of size nListSize, we return the max val
/// \param nListSize
/// \param naList
/// \return Maximum value in the supplied array.
int ListOP_GetMaxVal_Int(const int nListSize, const int* const naList)
{
    int i;
    int maxVal = naList[0];

    for (i = 0; i < nListSize; i++)
        {
            if (naList[i] > maxVal)
                {
                    maxVal = naList[i];
                }
        }

    return maxVal;
}


///
/// \param nVal
/// \param naOutList
/// \param nListSize
/// \param naInputList
/// \return
int ListOP_IndicesForVal_Int(const int nVal, int* const naOutList, const int nListSize, const int* const naInputList)
{
    int i;
    int nNumberOfInds = 0;
    for (i = 0; i < nListSize; i++)
        {
            if (naInputList[i] == nVal)
                {
                    naOutList[nNumberOfInds] = i;
                    nNumberOfInds++;
                }
        }
    return nNumberOfInds;
}


///
/// \param nVal
/// \param nListSize
/// \param naList
/// \return
int ListOP_GetRandomIndexForVal_Int(const int nVal, const int nListSize, const int* const naList)
{
    int nOutIdx          = -1;
    int* const naIdxList = (int*) malloc(sizeof(int) * nListSize);

    const int nTotIds = ListOP_IndicesForVal_Int(nVal, naIdxList, nListSize, naList);

    if (nTotIds == 1)
        {
            nOutIdx = naIdxList[0];
        }
    else if (nTotIds > 1)
        {
            nOutIdx = naIdxList[rand() % nTotIds];
        }

    free(naIdxList);

    return nOutIdx;
}


/// ListOP_UniqueElementsOfList_Int: Overwrite the provided _unsorted_ list of integers with only the unique elements,
/// and return the nSize of the new array. The list is sorted first, and then ListOP_UniqueElementsOfSortedList_Int is
/// run on that newly sorted list. This is a wrapper-like function for ease. The overwritten list is also sorted. \param
/// naUnsortedList: Array of ints \return number of unique elements.
int ListOP_UniqueElementsOfList_Int(const int nSize, int* const naUnsortedList)
{
    qsort(naUnsortedList, nSize, sizeof(int), UtilFunc_CompareInts);
    return ListOP_UniqueElementsOfSortedList_Int(nSize, naUnsortedList);
}


/// ListOP_Get2ndLargestVal_Int - Given an integer array naInputList, of size nListSize, we return the second largest
/// value in the array. If the array only has one value, we just return that value back. We get the sorted unique
/// elements in the input array and return the second last. \param nListSize \param naInputList \return 2nd largest
/// value in the input array.
int ListOP_Get2ndLargestVal_Int(const int nListSize, const int* const naInputList)
{
    int* naTmpList = (int*) malloc(nListSize * sizeof(int));
    int i;
    for (i = 0; i < nListSize; i++)
        {
            naTmpList[i] = naInputList[i];
        }

    // Extract Unique-sizes (this is sorted)
    const int nUniqueVals = ListOP_UniqueElementsOfList_Int(nListSize, naTmpList);

    int thisVal = naTmpList[0];

    if (nUniqueVals > 1)
        {
            thisVal = naTmpList[nUniqueVals - 2];
        }

    free(naTmpList);

    return thisVal;
}


/// ListOP_GenHistFromCounts_Int - For every occurence in the naInputList, we increment that index by 1 in naOutFreqHist
/// naOutFreqHist[naInputList[i]]++ over all elements
/// \param naOutFreqHist
/// \param naInputList
/// \param nSize
void ListOP_GenHistFromCounts_Int(int* const restrict naOutFreqHist, const int* const restrict naInputList,
                                  const int nSize)
{
    int i;
    for (i = 0; i < nSize; i++)
        {
            naOutFreqHist[naInputList[i]]++;
        }
}

