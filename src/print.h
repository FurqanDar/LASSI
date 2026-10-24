#ifndef _PRINT_H_ // include guard
#define _PRINT_H_

#include "global.h"

/*
 * Public print / I/O API used outside the print*.c family (mainly main.c).
 * Sanity-fail printers are also forward-declared from structure_sanity.c
 * (no #include "print.h" — Phase 3 dependency direction).
 * Cross-TU helpers for print*.c only live in print_internal.h.
 */

void ScreenIO_Print_KeyFile(void);

void ScreenIO_Print_SanityCheckFailurePreamble(long nGen, int run_cycle);

void ScreenIO_Print_SanityFail_BeadPosAndLattPos(int badBead);

void ScreenIO_Print_SanityFail_MolecularStructure(int badBead);

void ScreenIO_Print_SanityFail_SelfBond(int badBead);

void ScreenIO_Print_SanityFail_BeadBondSymmetry(int badBead);

void ScreenIO_Print_SanityFail_BeadBondDistance(int badBead);

void DataPrinting_Thermalization(const long nGen);

void DataPrinting_DuringRunCycles(const long nGen, const int run_it);

void DataAnalysis_DuringRunCycles(const long nGen, const int run_it);

void FileIO_CreateRunningDataFiles(void);

void FileIO_PrintCrashSnapshot(void);

void FileIOUtil_PreCycle_Init(const int run_it);

void FileIO_WriteRestart_ForThermalization(void);

void FileIO_PostCycle_WriteSystemRestart(const int run_it);

void FileIO_PostCycle_WriteCycleAvgData(const int run_it);

#endif // _PRINT_H_
