/*
 * Tiny pure-core checks for Phase 5.1 lattice index / distance helpers.
 * No simulation globals required for these entry points.
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "src/structure.h"

int main(void)
{
    /* i + box_x * (j + box_y * k) */
    assert(Lat_Ind_FromCoords_Box(1, 2, 3, 10, 20) == 1 + 10 * (2 + 20 * 3));
    assert(Lat_Ind_FromCoords_Box(0, 0, 0, 5, 5) == 0);
    assert(Lat_Ind_FromCoords_Box(4, 4, 4, 5, 5) == 4 + 5 * (4 + 5 * 4));

    {
        const int a[POS_MAX] = {0, 0, 0};
        const int b[POS_MAX] = {3, 4, 0};
        const int box[POS_MAX] = {20, 20, 20};
        const float d = Dist_PointToPoint_Box(a, b, box);
        assert(fabsf(d - 5.f) < 1e-6f);
    }

    {
        /* Min-image: delta 18 on box 20 -> image distance 2 */
        const int a[POS_MAX] = {0, 0, 0};
        const int b[POS_MAX] = {18, 0, 0};
        const int box[POS_MAX] = {20, 20, 20};
        const float d = Dist_PointToPoint_Box(a, b, box);
        assert(fabsf(d - 2.f) < 1e-6f);
    }

    {
        const int v[POS_MAX] = {3, 4, 0};
        assert(Dist_VecMagSq(v) == 25);
        assert(fabsf(Dist_Vec3n(v) - 5.f) < 1e-6f);
    }

    puts("test_lattice_index: ok");
    return 0;
}
