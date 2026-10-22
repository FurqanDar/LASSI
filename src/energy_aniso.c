#include "energy.h"
#include "structure.h"

/*
 * Anisotropic SC-SC bond energy helpers (Phase 5.4).
 * Pure FromCoupling cores take the table entry; wrappers bind g_system / g_energy.
 */

/// Energy_Aniso_ScSc_FromCoupling - SC-SC term from a coupling constant.
/// Pure: no globals.
float Energy_Aniso_ScSc_FromCoupling(const float e_sc_sc)
{
    return e_sc_sc;
}

/// Energy_Aniso_ScSc_Half_FromCoupling - SC-SC with float half (double-count).
/// Pure: no globals. Uses * 0.5f (same as For_Chain / For_List).
float Energy_Aniso_ScSc_Half_FromCoupling(const float e_sc_sc)
{
    return e_sc_sc * 0.5f;
}

/// Energy_Aniso_ScSc_HalfDiv_FromCoupling - SC-SC with / 2. (double-count).
/// Pure: no globals. Preserves Contiguous_Range / With_List float semantics.
float Energy_Aniso_ScSc_HalfDiv_FromCoupling(const float e_sc_sc)
{
    return e_sc_sc / 2.;
}

/// Energy_Anisotropic - returns the bond energy for beadID if it is bonded.
/// \param beadID
/// \return
float Energy_Anisotropic(const int beadID)
{ // Calculates the SC-SC energy of the bead in question.
    const int bP = g_system.bead_info[beadID][BEAD_FACE];
    const int resi = g_system.bead_info[beadID][BEAD_TYPE];
    int resj;
    if (bP != -1)
    {
        resj = g_system.bead_info[bP][BEAD_TYPE];
        return Energy_Aniso_ScSc_FromCoupling(g_energy.faEnergy[resi][resj][E_SC_SC]);
    }
    else
    {
        return 0.f;
    }
}

/// Energy_Anisotropic_Self - returns the bond energy for beadID only if it is
/// bonded to a bead in the same chain. Note that this does not account for
/// double counting \param beadID \return
float Energy_Anisotropic_Self(const int beadID)
{                      // Calculates the SC-SC energy of the bead in question.
    float totEn = 0.0; // Storing total overlap energy
    int bP = g_system.bead_info[beadID][BEAD_FACE];
    if (bP != -1)
    { // This bead does have a bond
        if (g_system.bead_info[beadID][BEAD_CHAINID] == g_system.bead_info[bP][BEAD_CHAINID])
        {
            totEn += Energy_Aniso_ScSc_FromCoupling(
                g_energy.faEnergy[g_system.bead_info[beadID][BEAD_TYPE]][g_system.bead_info[bP][BEAD_TYPE]][E_SC_SC]);
        }
    }
    return totEn;
}

/// Energy_Anisotropic_For_Chain - returns the bond energy for beadID if it is
/// bonded If it is bonded to another in the same chain, divide the energy by 2.
/// This should only be used when calculating the energy of an entire chain!
/// \param beadID
/// \return
float Energy_Anisotropic_For_Chain(const int beadID)
{                      // Calculates the SC-SC energy of the bead in question.
    float totEn = 0.f; // Storing total overlap energy
    const int bP = g_system.bead_info[beadID][BEAD_FACE];
    if (bP != -1)
    { // This bead does have a bond
        if (g_system.bead_info[beadID][BEAD_CHAINID] == g_system.bead_info[bP][BEAD_CHAINID])
        {
            totEn += Energy_Aniso_ScSc_Half_FromCoupling(
                g_energy.faEnergy[g_system.bead_info[beadID][BEAD_TYPE]][g_system.bead_info[bP][BEAD_TYPE]][E_SC_SC]);
        }
        else
        {
            totEn += Energy_Aniso_ScSc_FromCoupling(
                g_energy.faEnergy[g_system.bead_info[beadID][BEAD_TYPE]][g_system.bead_info[bP][BEAD_TYPE]][E_SC_SC]);
        }
    }
    return totEn;
}

/// Energy_Anisotropic_For_Range - returns the bond energy for beadID if it is
/// bonded. This function assumes that we get a beadID range [smallest_bead,
/// largest_bead], and that we will loop over every bead in this range. So if
/// the bead partner is within this range, we divide by 2 to account for double
/// counting. \param beadID \param smBead \param lgBead \return
float Energy_Anisotropic_For_Range(const int beadID, const int smBead, const int lgBead)
{
    float totEn = 0.f;
    const int bP = g_system.bead_info[beadID][BEAD_FACE];
    if (bP != -1)
    {
        if ((beadID <= lgBead) && (beadID >= smBead))
        {
            totEn += Energy_Aniso_ScSc_Half_FromCoupling(
                g_energy.faEnergy[g_system.bead_info[beadID][BEAD_TYPE]][g_system.bead_info[bP][BEAD_TYPE]][E_SC_SC]);
        }
        else
        {
            totEn += Energy_Aniso_ScSc_FromCoupling(
                g_energy.faEnergy[g_system.bead_info[beadID][BEAD_TYPE]][g_system.bead_info[bP][BEAD_TYPE]][E_SC_SC]);
        }
    }
    return totEn;
}

/// Energy_Anisotropic_Contiguous_Range - returns the bond energy for beadID if
/// it is bonded. This function assumes that we get a beadID range
/// [smallest_bead, largest_bead], and that we will loop over every bead in this
/// range. So if the bead partner is within this range, we divide by 2 to
/// account for double counting. \param beadID \return
float Energy_Anisotropic_Contiguous_Range(const int beadID, const int smallest_bead, const int largest_bead)
{
    float totEn = 0.0; // Storing total overlap energy
    int bP = g_system.bead_info[beadID][BEAD_FACE];
    if (bP != -1)
    { // This bead does have a bond
        if (g_system.bead_info[beadID][BEAD_CHAINID] == g_system.bead_info[bP][BEAD_CHAINID])
        { // Intra-molecular
            if (bP >= smallest_bead && bP <= largest_bead)
            { // Within subset
                totEn += Energy_Aniso_ScSc_HalfDiv_FromCoupling(
                    g_energy
                        .faEnergy[g_system.bead_info[beadID][BEAD_TYPE]][g_system.bead_info[bP][BEAD_TYPE]][E_SC_SC]);
            }
            else
            {
                totEn += Energy_Aniso_ScSc_FromCoupling(g_energy.faEnergy[g_system.bead_info[beadID][BEAD_TYPE]]
                                                                         [g_system.bead_info[bP][BEAD_TYPE]][E_SC_SC]);
            }
        }
        else
        {
            totEn += Energy_Aniso_ScSc_FromCoupling(
                g_energy.faEnergy[g_system.bead_info[beadID][BEAD_TYPE]][g_system.bead_info[bP][BEAD_TYPE]][E_SC_SC]);
        }
    }
    return totEn;
}

/// Energy_Anisotropic_With_List - returns the bond energy for beadID if it is
/// bonded. This function assumes that we get a list of beadID's of size
/// list_size, but that the last bead is at list_size-1. This also assumes that
/// we will loop over the entire list So if the bead partner is within this
/// list, we divide by 2 to account for double counting. \param beadID \return
float Energy_Anisotropic_With_List(const int beadID, const int* bead_list, const int list_size)
{
    float totEn = 0.0; // Storing total overlap energy
    int bP = g_system.bead_info[beadID][BEAD_FACE];
    int check_bead = 0;
    int i;
    if (bP != -1)
    { // This bead does have a bond
        for (i = 0; i < list_size; i++)
        { // Checking if we should worry about double counting
            if (bP == bead_list[i])
            {
                check_bead = 1;
                break;
            }
        }
        if (check_bead == 1)
        {
            totEn += Energy_Aniso_ScSc_HalfDiv_FromCoupling(
                g_energy.faEnergy[g_system.bead_info[beadID][BEAD_TYPE]][g_system.bead_info[bP][BEAD_TYPE]][E_SC_SC]);
        }
        else
        {
            totEn += Energy_Aniso_ScSc_FromCoupling(
                g_energy.faEnergy[g_system.bead_info[beadID][BEAD_TYPE]][g_system.bead_info[bP][BEAD_TYPE]][E_SC_SC]);
        }
    }
    return totEn;
}

float Energy_Anisotropic_For_List(const int beadID, const int listSize, const int beadList[MAX_BONDS + 1])
{
    const int bP = g_system.bead_info[beadID][BEAD_FACE];
    const int resi = g_system.bead_info[beadID][BEAD_TYPE];
    int resj;
    if (bP != -1)
    {
        resj = g_system.bead_info[bP][BEAD_TYPE];
        if (Check_BeadID_InList(bP, listSize, beadList))
        {
            return Energy_Aniso_ScSc_Half_FromCoupling(g_energy.faEnergy[resi][resj][E_SC_SC]);
        }
        else
        {
            return Energy_Aniso_ScSc_FromCoupling(g_energy.faEnergy[resi][resj][E_SC_SC]);
        }
    }
    else
    {
        return 0.f;
    }
}
