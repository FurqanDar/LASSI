#ifndef _STRUCTURE_INTERNAL_H_ // include guard
#define _STRUCTURE_INTERNAL_H_

/*
 * Cross-TU declarations for the structure family (structure*.c only).
 * Not part of the public API; external callers should use structure.h.
 */

#include "structure.h"

float Dist_PointToPoint(const int* f1, const int* f2);

void LatPos_add_wPBC_ofComp(int* outVec, const int* firVec, const int* secVec, const int compNum);

#endif // _STRUCTURE_INTERNAL_H_
