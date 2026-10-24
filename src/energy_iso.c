#include "energy_internal.h"

/*
 * Isotropic pair / solvation energy kernels (Phase 5.3).
 * Pure FromCoupling cores take table entries + distance; wrappers bind g_energy.
 */

/// Energy_Iso_Ovlp_FromCoupling - OVLP term from a coupling constant.
/// Pure: no globals. xDis is unused (OVLP is distance-independent in LASSI).
float Energy_Iso_Ovlp_FromCoupling(const float e_ovlp, const float xDis)
{
    (void)xDis;
    return e_ovlp;
}

/// Energy_Iso_Cont_FromCoupling - CONT term e_cont / xDis.
/// Pure: no globals.
float Energy_Iso_Cont_FromCoupling(const float e_cont, const float xDis)
{
    return e_cont / xDis;
}

/// Energy_Iso_fSol_FromCoupling - solvation coupling for one bead type.
/// Pure: no globals.
float Energy_Iso_fSol_FromCoupling(const float e_fsol)
{
    return e_fsol;
}

/// Energy_Iso_PairDelta_FromCoupling - OVLP + CONT for one neighbor pair.
/// Pure: no globals.
float Energy_Iso_PairDelta_FromCoupling(const float e_ovlp, const float e_cont, const float xDis)
{
    return Energy_Iso_Ovlp_FromCoupling(e_ovlp, xDis) + Energy_Iso_Cont_FromCoupling(e_cont, xDis);
}

/// Energy_Iso_EmptySiteSol_FromCoupling - solvation for empty neighbor slots
/// (empty_neigh_count * e_fsol). Pure: no globals.
float Energy_Iso_EmptySiteSol_FromCoupling(const float e_fsol, const int empty_neigh_count)
{
    return (float)empty_neigh_count * e_fsol;
}

/// Energy_Iso_Ovlp - calculate the OVLP interaction between two beads of types
/// beadType1 and beadType2, respectively. \param beadType1: Type of first bead.
/// \param beadType2: Type of second bead.
/// \param xDis: Distance between the two beads.
/// \return g_energy.faEnergy[beadType1][beadType2][E_OVLP].
float Energy_Iso_Ovlp(int const beadType1, int const beadType2, float const xDis)
{
    return Energy_Iso_Ovlp_FromCoupling(g_energy.faEnergy[beadType1][beadType2][E_OVLP], xDis);
}

/// Energy_Iso_Cont - calculate the CONT interaction between two beads of types
/// beadType1 and beadType2, respectively. \param beadType1: Type of first bead.
/// \param beadType2: Type of second bead.
/// \param xDis: Distance between the two beads.
/// \return g_energy.faEnergy[beadType1][beadType2][E_CONT]/xDis.
float Energy_Iso_Cont(int const beadType1, int const beadType2, float const xDis)
{
    return Energy_Iso_Cont_FromCoupling(g_energy.faEnergy[beadType1][beadType2][E_CONT], xDis);
}

/// Energy_Iso_fSol: Returns the solvation energy of this bead.
/// \param beadType
/// \return
float Energy_Iso_fSol(int const beadType)
{
    return Energy_Iso_fSol_FromCoupling(g_energy.faEnergy[beadType][beadType][E_F_SOL]);
}

/// Energy_Iso_PairDelta - OVLP + CONT between two bead types at distance xDis.
float Energy_Iso_PairDelta(int const beadType1, int const beadType2, float const xDis)
{
    return Energy_Iso_PairDelta_FromCoupling(
        g_energy.faEnergy[beadType1][beadType2][E_OVLP],
        g_energy.faEnergy[beadType1][beadType2][E_CONT],
        xDis);
}

/// Energy_Iso_EmptySiteSol - solvation contribution for empty_neigh_count vacant
/// neighbor sites around a bead of beadType (binds g_energy).
float Energy_Iso_EmptySiteSol(const int beadType, const int empty_neigh_count)
{
    return Energy_Iso_EmptySiteSol_FromCoupling(g_energy.faEnergy[beadType][beadType][E_F_SOL], empty_neigh_count);
}
