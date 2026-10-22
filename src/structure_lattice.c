#include "structure.h"

/*
 * Lattice index, min-image distances, and lattice-position helpers.
 * Pure cores take explicit box extents; thin wrappers bind g_system.
 */

/// Lat_Ind_FromCoords_Box - 1D lattice index for (i,j,k) with explicit box extents.
/// Pure: no globals. Layout is i + box_x * (j + box_y * k).
int Lat_Ind_FromCoords_Box(const int i, const int j, const int k, const int box_x, const int box_y)
{
    return i + box_x * (j + box_y * k);
}

/// Lat_Ind_FromCoords - helper function to get the correct 1D index of this position
/// \param i
/// \param j
/// \param k
/// \return the 1D index location for (i,j,k) position
int Lat_Ind_FromCoords(const int i, const int j, const int k)
{
    return Lat_Ind_FromCoords_Box(i, j, k, g_system.naBoxSize[POS_X], g_system.naBoxSize[POS_Y]);
}

/// Lat_Ind_FromVec - returns the 1D index given the array xArr
/// \param xArr
/// \return
int Lat_Ind_FromVec(const int* const xArr)
{
    return Lat_Ind_FromCoords(xArr[POS_X], xArr[POS_Y], xArr[POS_Z]);
}

/// Lat_Ind_OfBead - returns the 1D index of this bead's location
/// \param beadID
/// \return
int Lat_Ind_OfBead(const int beadID)
{
    return Lat_Ind_FromVec(g_system.bead_info[beadID]);
}

// Note that the distance functions account for periodic boundaries

/// Dist_PointToPoint_Box - euclidean min-image distance between integer points with explicit box.
/// Pure: no globals. Same arithmetic as Dist_PointToPoint.
float Dist_PointToPoint_Box(const int* const restrict f1, const int* const restrict f2,
                            const int* const restrict box)
{
    int d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
    {
        d[i] = abs(f1[i] - f2[i]);
        d[i] = d[i] > box[i] / 2 ? box[i] - d[i] : d[i];
    }

    return sqrtf((float)(d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_PointTotPoint_Float - euclidean distance between the vectors (arrays) f1 and f2 where f1
/// and f2 are floats. \param f1 \param f2 \return
float Dist_PointToPoint_Float(const float* f1, const float* f2)
{
    float d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
    {
        d[i] = fabsf(f1[i] - f2[i]);
        d[i] =
            d[i] > (float)g_system.naBoxSize[i] / 2. ? (float)g_system.naBoxSize[i] - d[i] : d[i];
    }

    return sqrtf(d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]);
}

/// Dist_PointToPoint - euclidean distance between the vectors (arrays) f1 and f2 where f1 and f2
/// are integer points. \param f1 \param f2 \return
float Dist_PointToPoint(const int* const restrict f1, const int* const restrict f2)
{
    return Dist_PointToPoint_Box(f1, f2, g_system.naBoxSize);
}

/// Dist_BeadToPoint - euclidean distance between beadID and the vector f1.
/// \param beadID
/// \param f1
/// \return sqrt(dx^2 + dy^2 + dz^2)
float Dist_BeadToPoint(const int beadID, const int* const f1)
{
    return Dist_PointToPoint(g_system.bead_info[beadID], f1);
}

/// Dist_BeadToPoint_Double - euclidean distance between beadID and the double array f1.
/// \param beadID
/// \param f1
/// \return sqrt(dx^2 + dy^2 + dz^2)
float Dist_BeadToPoint_Double(const int beadID, const lDub* f1)
{
    lDub d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
    {
        d[i] = fabs((lDub)g_system.bead_info[beadID][i] - f1[i]);
        d[i] = d[i] > (lDub)g_system.naBoxSize[i] / 2. ? (lDub)g_system.naBoxSize[i] - d[i] : d[i];
    }
    return sqrtf((float)(d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_BeadToPoint_Float - euclidean distance between beadID and the float array f1.
/// \param beadID
/// \param f1
/// \return sqrt(dx^2 + dy^2 + dz^2)
float Dist_BeadToPoint_Float(const int beadID, const float* f1)
{
    float d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
    {
        d[i] = fabsf((float)g_system.bead_info[beadID][i] - f1[i]);
        d[i] =
            d[i] > (float)g_system.naBoxSize[i] / 2. ? (float)g_system.naBoxSize[i] - d[i] : d[i];
    }
    return sqrtf((d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_BeadToBead - euclidean distance between the two beads.
/// \param n1
/// \param n2
/// \return sqrt(dx^2 + dy^2 + dz^2)
float Dist_BeadToBead(const int n1, const int n2)
{
    return Dist_PointToPoint(g_system.bead_info[n1], g_system.bead_info[n2]);
}

/// Dist_Vec3n - non periodic boundary euclidean magnitude of vector
/// \param f1: The array where indicies 0,1 and 3 correspond to x y and z.
/// \return
float Dist_Vec3n(const int* const f1)
{ // Outputs the magnitude of the vector
    return sqrtf((float)(f1[0] * f1[0] + f1[1] * f1[1] + f1[2] * f1[2]));
}

/// Dist_VecMagSq - non periodic boundary euclidean magnitude of vector
/// \param f1: The array where indicies 0,1 and 3 correspond to x y and z.
/// \return Square of magnitude: x^2 + y^2 + z^2
int Dist_VecMagSq(const int* const f1)
{ // Outputs the magnitude of the vector
    return (f1[0] * f1[0] + f1[1] * f1[1] + f1[2] * f1[2]);
}

void Calculate_Distances_For_Radius(float* thisList, const int nRad)
{
    int r_disp[POS_MAX] = {0};

    int dum_iter = 0;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
    {
        for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
        {
            for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
            {
                thisList[dum_iter++] = Dist_Vec3n(r_disp);
            }
        }
    }
}

/// LatPos_copy: Copies the position elements to copy_vec.
/// \param copy_vec
/// \param in_vec
void LatPos_copy(int* copy_vec, const int* in_vec)
{
    copy_vec[0] = in_vec[0];
    copy_vec[1] = in_vec[1];
    copy_vec[2] = in_vec[2];
}

/// LatPos_gen_rand_wRad. Given a radius, we generate a random set of coordinates between -Radius
/// and Radius for each dimension. The result is stored in outVec \param outVec \param nRadius
void LatPos_gen_rand_wRad(int outVec[POS_MAX], const int nRadius)
{
    int j;
    int radUp = 2 * nRadius + 1;
    for (j = 0; j < POS_MAX; j++)
    {
        outVec[j] = rand() % radUp;
    }
    for (j = 0; j < POS_MAX; j++)
    {
        outVec[j] -= nRadius;
    }
}

/// LatPos_add_wPBC: Given firVec and secVec, we add the two vectors, and take care of PBCs, storing
/// in outVec. The implementation assumes that the sum of any coordinate can not be larger than
/// twice the lattice-size, because there is no explicit modulo operation. \param outVec \param
/// firVec \param secVec
void LatPos_add_wPBC(int* outVec, const int* firVec, const int* secVec)
{
    short j;
    for (j = 0; j < POS_MAX; j++)
    {
        outVec[j] = firVec[j] + secVec[j];
    }
    for (j = 0; j < POS_MAX; j++)
    {
        outVec[j] = outVec[j] < 0 ? outVec[j] + g_system.naBoxSize[j] : outVec[j];
    }
    for (j = 0; j < POS_MAX; j++)
    {
        outVec[j] =
            outVec[j] >= g_system.naBoxSize[j] ? outVec[j] - g_system.naBoxSize[j] : outVec[j];
    }
}

/// LatPos_add_wPBC_ofComp: Given firVec and secVec, we add the compNum component of the two
/// vectors, and take care of PBCs, storing in outVec. The implementation assumes that the sum of
/// any coordinate can not be larger than twice the lattice-size, because there is no explicit
/// modulo operation. \param outVec \param firVec \param secVec
void LatPos_add_wPBC_ofComp(int* outVec, const int* firVec, const int* secVec, const int compNum)
{
    outVec[compNum] = firVec[compNum] + secVec[compNum];
    outVec[compNum] =
        outVec[compNum] < 0 ? outVec[compNum] + g_system.naBoxSize[compNum] : outVec[compNum];
    outVec[compNum] = outVec[compNum] >= g_system.naBoxSize[compNum]
                          ? outVec[compNum] - g_system.naBoxSize[compNum]
                          : outVec[compNum];
}

/// LatPos_add_noPBC: Given the two arrays firVec and secVec, we store the sum of the two arrays in
/// outVec. \param outVec \param firVec \param secVec
void LatPos_add_noPBC(int* outVec, const int* firVec, const int* secVec)
{
    short j;
    for (j = 0; j < POS_MAX; j++)
    {
        outVec[j] = firVec[j] + secVec[j];
    }
}

/// BeadPos_sub_wPBC: Given firVec and secVec, we subtract second from first, and take care of PBCs,
/// storing in outVec. The implementation assumes that the sum of any coordinate can not be larger
/// than the lattice size because there is no explicit modulo operation. This version is
/// specifically for bead-positions, or for inter-bead calculations. \param outVec \param firVec
/// \param secVec
void BeadPos_sub_wPBC(int* outVec, const int* firVec, const int* secVec)
{
    short j;
    for (j = 0; j < POS_MAX; j++)
    {
        outVec[j] = firVec[j] - secVec[j];
    }
    for (j = 0; j < POS_MAX; j++)
    {
        outVec[j] = outVec[j] < -g_system.naBoxSize[j] / 2.f ? outVec[j] + g_system.naBoxSize[j]
                                                             : outVec[j]; // pht
    }
    for (j = 0; j < POS_MAX; j++)
    {
        outVec[j] = outVec[j] >= g_system.naBoxSize[j] / 2.f ? outVec[j] - g_system.naBoxSize[j]
                                                             : outVec[j]; // pht
    }
}

int Vec3n_DotProd(const int* vec_1, const int* vec_2)
{
    return vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2];
}

/// Vec3n_CosTheta - calculates the Cos(theta) between the vectors v1 and v2.
/// Cos(theta) = v_1 . v_2 / |v_1| / |v_2|
/// \param v1
/// \param v2
/// \return
float Vec3n_CosTheta(const int* v1, const int* v2)
{
    const int dotProd = Vec3n_DotProd(v1, v2);
    const int vec1MagSq = Vec3n_DotProd(v1, v1);
    const int vec2MagSq = Vec3n_DotProd(v2, v2);
    const float dumDenom = (float)vec1MagSq * (float)vec2MagSq;
    const float dumCosTheta = (float)dotProd / sqrtf(dumDenom);

    return dumCosTheta;
}

float Vec3n_AngleBetweenVecs(const int* vec_1, const int* vec_2)
{
    const float dumAng = Vec3n_CosTheta(vec_1, vec_2);
    return acosf(dumAng);
}
