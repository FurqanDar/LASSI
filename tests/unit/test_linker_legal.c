/*
 * Tiny pure-core checks for Phase 5.2 linker legality helpers.
 */
#include <assert.h>
#include <stdio.h>

#include "src/structure.h"

int main(void)
{
    /* LINKER_RSCALE is 1.74f: dist 1.0 with len 1 is OK; 2.0 is broken. */
    assert(LinkerConstraint_DistOK(1.0f, 1) == 1);
    assert(LinkerConstraint_DistOK(1.74f, 1) == 1);
    assert(LinkerConstraint_DistOK(1.74001f, 1) == 0);
    assert(LinkerConstraint_DistOK(3.0f, 2) == 1); /* 3 <= 1.74*2 */
    assert(LinkerConstraint_DistOK(3.5f, 2) == 0); /* 3.5 > 3.48 */

    {
        const int list[MAX_BONDS + 1] = {10, 20, 30};
        assert(Check_BeadID_InList(20, 3, list) == 1);
        assert(Check_BeadID_InList(99, 3, list) == 0);
        assert(Check_BeadID_InList(10, 0, list) == 0);
    }

    puts("test_linker_legal: ok");
    return 0;
}
