/*
 * Tiny pure-core checks for Phase 5.4 anisotropic SC-SC kernels.
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "src/energy.h"

int main(void)
{
    assert(fabsf(Energy_Aniso_ScSc_FromCoupling(-2.5f) - (-2.5f)) < 1e-6f);
    assert(fabsf(Energy_Aniso_ScSc_Half_FromCoupling(4.0f) - 2.0f) < 1e-6f);
    assert(fabsf(Energy_Aniso_ScSc_HalfDiv_FromCoupling(4.0f) - 2.0f) < 1e-6f);
    assert(fabsf(Energy_Aniso_ScSc_Half_FromCoupling(-1.0f) - (-0.5f)) < 1e-6f);
    assert(fabsf(Energy_Aniso_ScSc_HalfDiv_FromCoupling(-1.0f) - (-0.5f)) < 1e-6f);

    puts("test_energy_aniso: ok");
    return 0;
}
