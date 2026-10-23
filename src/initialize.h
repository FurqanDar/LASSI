#ifndef _INITIALIZE_H_ // include guard
#define _INITIALIZE_H_

#include "global.h"

/* Phase 6.2: init entrypoints take const SimulationConfig* for keyfile-derived
 * inputs (box, MC/report freqs, restart path, …). Runtime topology (beads /
 * chains) still lives in g_system after structure parse. */

void Memory_Initialization_AtStart(const SimulationConfig* cfg);

void Memory_Allocate_NeighborLists(void);

void Memory_VerifyMalloc(void);

void Global_Array_Initialization_AtStart(const SimulationConfig* cfg);

void Reset_Counters(void);

void Reset_AnalysisCycleScratch(void);

void Reset_McCycleScratch(void);

void Reset_Global_Arrays(void);

void Initial_Conditions_Simple(const SimulationConfig* cfg);

void Initial_Conditions_FromFile(const SimulationConfig* cfg);

void Initial_Conditions_BreakBonds(void);

void BiasPotential_TurnOFF(void);

float Temperature_Function(const int mode, const long nGen);

void Calculate_Rot_Bias(const float CurrentTemp);

int* Create1DInt(const size_t xDim, const char* ArrName);

long* Create1DLong(const size_t xDim, const char* ArrName);

float* Create1DFloat(const size_t xDim, const char* ArrName);

long double* Create1DLongdouble(const size_t xDim, const char* ArrName);

int** Create2DInt(const size_t xDim, const size_t yDim, const char* ArrName);

long double** Create2DLongdouble(const size_t xDim, const size_t yDim, const char* ArrName);

#endif // _INITIALIZE_H_
