/*
 * Tiny pure-core checks for Phase 5.3 isotropic energy kernels.
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "src/energy.h"

int main(void)
{
    assert(fabsf(Energy_Iso_Ovlp_FromCoupling(1.5f, 9.0f) - 1.5f) < 1e-6f);
    assert(fabsf(Energy_Iso_Cont_FromCoupling(4.0f, 2.0f) - 2.0f) < 1e-6f);
    assert(fabsf(Energy_Iso_fSol_FromCoupling(-0.25f) - (-0.25f)) < 1e-6f);
    assert(fabsf(Energy_Iso_PairDelta_FromCoupling(1.0f, 4.0f, 2.0f) - 3.0f) < 1e-6f);
    assert(fabsf(Energy_Iso_EmptySiteSol_FromCoupling(0.5f, 26) - 13.0f) < 1e-6f);
    assert(fabsf(Energy_Iso_EmptySiteSol_FromCoupling(0.5f, 26 - 4) - 11.0f) < 1e-6f);

    puts("test_energy_iso: ok");
    return 0;
}
