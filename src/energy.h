#ifndef _ENERGY_H_ // include guard
#define _ENERGY_H_

#include "global.h"

/*
 * Public energy API (move deltas, system total, topo/bias) plus pure FromCoupling
 * cores kept public for unit tests. Family-only decls live in energy_internal.h.
 */

void Energy_Total_System(void);

float Energy_BiasingPotential(const int beadID);

float Energy_Topo_Angle(const int beadID);

float Energy_Topo_Angle_ForList(const int bondNum, const int* bondList);

/* Pure isotropic FromCoupling cores (see energy_iso.c); unit-tested. */

float Energy_Iso_Ovlp_FromCoupling(const float e_ovlp, const float xDis);

float Energy_Iso_Cont_FromCoupling(const float e_cont, const float xDis);

float Energy_Iso_fSol_FromCoupling(const float e_fsol);

float Energy_Iso_PairDelta_FromCoupling(const float e_ovlp, const float e_cont, const float xDis);

float Energy_Iso_EmptySiteSol_FromCoupling(const float e_fsol, const int empty_neigh_count);

/* Pure anisotropic FromCoupling cores (see energy_aniso.c); unit-tested. */

float Energy_Aniso_ScSc_FromCoupling(const float e_sc_sc);

float Energy_Aniso_ScSc_Half_FromCoupling(const float e_sc_sc);

float Energy_Aniso_ScSc_HalfDiv_FromCoupling(const float e_sc_sc);

/* Isotropic move-delta helpers used by mcmove (see energy.c). */

void Energy_Iso_ForLocal(
    const int beadID,
    const int resi,
    const int* r_pos0,
    long double* oldEn,
    long double* newEn,
    int* ovlp_num,
    int* cont_num,
    int* ovlp_neighs,
    int* cont_neighs);

void Energy_Iso_ForLocalEquil(
    const int beadID,
    const int resi,
    const int* r_pos0,
    long double* oldEn,
    long double* newEn,
    int* ovlp_num,
    int* cont_num,
    int* ovlp_neighs,
    int* cont_neighs);

void Energy_Iso_ForChains(
    const int beadID,
    long double* oldEn,
    long double* newEn,
    int* ovlp_num,
    int* cont_num,
    int* ovlp_neighs,
    int* cont_neighs);

void Energy_Iso_ForChainsEquil(
    const int beadID,
    long double* oldEn,
    long double* newEn,
    int* ovlp_num,
    int* cont_num,
    int* ovlp_neighs,
    int* cont_neighs);

void Energy_Iso_ForCoLocal(
    const int thisBeadID,
    const int otherBeadID,
    const int* r_pos0,
    long double* oldEn,
    long double* newEn,
    int* ovlp_num,
    int* cont_num,
    int* ovlp_neighs,
    int* cont_neighs);

void Energy_Iso_ForLists(
    const int beadIdx,
    int const listSize,
    const int beadList[MAX_BONDS + 1],
    const int beadPos[MAX_BONDS + 1][POS_MAX],
    long double* oldEn,
    long double* newEn,
    int* ovlp_num,
    int* cont_num,
    int* ovlp_neighs,
    int* cont_neighs);

void Energy_Iso_ForRange(
    const int beadID,
    const int smallestBead,
    const int largestBead,
    long double* oldEn,
    long double* newEn,
    int* ovlp_num,
    int* cont_num,
    int* ovlp_neighs,
    int* cont_neighs);

void Energy_Iso_ForRangeEquil(
    const int beadID,
    const int smallestBead,
    const int largestBead,
    long double* oldEn,
    long double* newEn,
    int* ovlp_num,
    int* cont_num,
    int* ovlp_neighs,
    int* cont_neighs);

#endif // _ENERGY_H_
