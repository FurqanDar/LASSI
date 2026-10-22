#include "print.h"
#include "initialize.h"
#include "cluster.h"
#include "energy.h"
#include "structure.h"


/// ForPrinting_GetReportState. Calculates if the current step is a multiple of the reporting frequency
/// for this particular report.
/// Convenient wrapper around some modulo arithmetic.
/// \param nGen
/// \param thisReport
/// \return
char ForPrinting_GetReportState(const long nGen, const long thisReport)
{
    char dum_log = (char) ((nGen % thisReport) == 0);
    dum_log      = dum_log ? 1 : 0;
    return dum_log;
}


/// DataPrinting_Thermalization - helper function for printing out data during the thermalization cycle.
/// No data analysis is performed during the thermalization procedure.
/// This function decides if, given the MCStep, it is time to print the following:
/// 1. The Log - to the screen.
/// 2. Trajectory - to the file (or saved if in total format)
/// 3. Energy - to the file.
/// 4. MCMove - to the file.
/// \param nGen
void DataPrinting_Thermalization(const long nGen)
{

    char cFlagForEnCal = 0;
    char cLogFlag      = 0;
    char cEnergyFlag   = 0;
    char cAccFlag      = 0;
    char cConfigFlag   = 0;

    if (g_reports.naReportFreqs[REPORT_LOG])
        {
            cLogFlag = ForPrinting_GetReportState(nGen, g_reports.naReportFreqs[REPORT_LOG]);
            if (cLogFlag)
                {
                    PerformRuntimeSanityChecks(nGen, -1);
                    Energy_Total_System();
                    cFlagForEnCal = 1;
                    ScreenIO_Print_Log_Thermalization(nGen);
                }
        }

    if (nGen)
        {
            if (g_reports.naReportFreqs[REPORT_CONFIG])
                {
                    cConfigFlag = ForPrinting_GetReportState(nGen, g_reports.naReportFreqs[REPORT_CONFIG]);
                    if (cConfigFlag)
                        {
                            FileIO_Trajectory_AppendFrame(g_reports.strFileTraj, -1, nGen);
                        }
                }

            if (g_reports.naReportFreqs[REPORT_ENERGY])
                {
                    // DO ENERGY FLAGS
                    cEnergyFlag = ForPrinting_GetReportState(nGen, g_reports.naReportFreqs[REPORT_ENERGY]);
                    if (cEnergyFlag)
                        {
                            if (! cFlagForEnCal)
                                {
                                    Energy_Total_System();
                                    cFlagForEnCal = 1;
                                }
                            FileIO_AppendEnergyTo_EnergyFile(g_reports.strFileEnergy, nGen);
                        }
                }

            if (g_reports.naReportFreqs[REPORT_MCMOVE])
                {
                    // DO MC_ACC PRINTING
                    cAccFlag = ForPrinting_GetReportState(nGen, g_reports.naReportFreqs[REPORT_MCMOVE]);
                    if (cAccFlag)
                        {
                            FileIO_WriteTo_MCMoveFile(g_reports.strFileMCMove, nGen, g_mc.fCuTemp);
                        }
                }
        }
}


/// DataPrinting_DuringRunCycles
/// \param nGen
/// \param run_it
void DataPrinting_DuringRunCycles(const long nGen, const int run_it)
{
    char cFlagForEnCal = 0;
    char cLogFlag      = 0;
    char cEnergyFlag   = 0;
    char cAccFlag      = 0;
    char cConfigFlag   = 0;

    if (g_reports.naReportFreqs[REPORT_LOG])
        {
            cLogFlag = ForPrinting_GetReportState(nGen, g_reports.naReportFreqs[REPORT_LOG]);
            if (cLogFlag)
                {
                    PerformRuntimeSanityChecks(nGen, run_it);
                    Energy_Total_System();
                    cFlagForEnCal = 1;
                    ScreenIO_Print_Log_FullRun(nGen, run_it);
                }
        }

    if (nGen)
        {
            if (g_reports.naReportFreqs[REPORT_CONFIG])
                {
                    cConfigFlag = ForPrinting_GetReportState(nGen, g_reports.naReportFreqs[REPORT_CONFIG]);
                    if (cConfigFlag)
                        {
                            FileIO_Trajectory_AppendFrame(g_reports.strFileTraj, run_it, nGen);
                        }
                }

            if (g_reports.naReportFreqs[REPORT_ENERGY])
                {
                    // DO ENERGY PRINTING
                    cEnergyFlag = ForPrinting_GetReportState(nGen, g_reports.naReportFreqs[REPORT_ENERGY]);
                    if (cEnergyFlag)
                        {
                            if (! cFlagForEnCal)
                                {
                                    Energy_Total_System();
                                    cFlagForEnCal = 1;
                                }
                            FileIO_AppendEnergyTo_EnergyFile(g_reports.strFileEnergy, nGen);
                        }
                }

            if (g_reports.naReportFreqs[REPORT_MCMOVE])
                {
                    // DO MC_ACC PRINTING
                    cAccFlag = ForPrinting_GetReportState(nGen, g_reports.naReportFreqs[REPORT_MCMOVE]);
                    if (cAccFlag)
                        {
                            FileIO_WriteTo_MCMoveFile(g_reports.strFileMCMove, nGen, g_mc.fCuTemp);
                        }
                }
        }
}


/// DataAnalysis_DuringRunCycles
/// \param nGen
/// \param run_it
void DataAnalysis_DuringRunCycles(const long nGen, const int run_it)
{

    //    char cRDF_flag  = 0;
    //    char cCOM_flag  = 0;
    //    char cCLUS_flag = 0;

    if (g_reports.naReportFreqs[REPORT_RDFTOT])
        { // SysProp is printed outside of this function in main.c, lol
            if (nGen % g_reports.naReportFreqs[REPORT_RDFTOT] == 0)
                {
                    // TODO: Update RDF Calculations
                    RDF_ComponentWise_Avg();
                }
        }
    if (g_reports.naReportFreqs[REPORT_COMDEN])
        { // SysProp is printed outside of this function in main.c, lol
            if (nGen % g_reports.naReportFreqs[REPORT_COMDEN] == 0)
                {
                    RadialDensityAnalysis_Perform_Analysis();
                }
        }
    if (g_reports.naReportFreqs[REPORT_NETWORK])
        { // SysProp is printed outside of this function in main.c, lol
            if (nGen % g_reports.naReportFreqs[REPORT_NETWORK] == 0)
                {
                    ClusAnalysis_Perform_Analysis(nClusteringMode_glb);
                    // TODO: Update Gyration calculations
                    GyrTensor_GyrRad_Avg();
                }
        }
}

