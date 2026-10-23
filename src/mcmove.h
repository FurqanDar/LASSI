#ifndef _MCMOVE_H_ // include guard
#define _MCMOVE_H_

#include "global.h"

/*
 * Public MC step API. Internal Move_*, OP_*, and Rot_* helpers live in
 * mcmove_internal.h (mcmove translation units only).
 */

int MC_Step(float fMCTemp);

int MC_Step_Therm(float fMCTemp);

#endif // _MCMOVE_H_
