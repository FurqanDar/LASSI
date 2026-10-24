#ifndef _PRINT_INTERNAL_H_ // include guard
#define _PRINT_INTERNAL_H_

/*
 * Cross-TU declarations for the print family (print*.c only).
 * Not part of the public API; external callers should use print.h.
 */

#include "print.h"

void ScreenIO_Print_Log_Thermalization(const long nGen);

void ScreenIO_Print_Log_FullRun(const long nGen, const int run_cycle);

void FileIO_WriteTo_MCMoveFile(const char* filename, const long nGen, const float fMCTemp);

void FileIO_AppendEnergyTo_EnergyFile(const char* fileNameStr, const long nGen);

void FileIO_Trajectory_AppendFrame(const char* fileNameStr, const int run_it, const long nGen);

#endif // _PRINT_INTERNAL_H_
