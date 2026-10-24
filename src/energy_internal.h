#ifndef _ENERGY_INTERNAL_H_ // include guard
#define _ENERGY_INTERNAL_H_

/*
 * Cross-TU declarations for the energy family (energy*.c only).
 * Not part of the public API; external callers should use energy.h.
 */

#include "energy.h"

/* Isotropic pair / solvation wrappers (energy_iso.c) — bind g_energy. */

float Energy_Iso_Ovlp(int const beadType1, int const beadType2, float const xDis);

float Energy_Iso_Cont(int const beadType1, int const beadType2, float const xDis);

float Energy_Iso_fSol(int const beadType);

float Energy_Iso_PairDelta(int const beadType1, int const beadType2, float const xDis);

float Energy_Iso_EmptySiteSol(const int beadType, const int empty_neigh_count);

/* Anisotropic SC-SC wrappers used by energy.c (energy_aniso.c). */

float Energy_Anisotropic(const int beadID);

float Energy_Anisotropic_Self(const int beadID);

float Energy_Anisotropic_For_Chain(const int beadID);

#endif // _ENERGY_INTERNAL_H_
