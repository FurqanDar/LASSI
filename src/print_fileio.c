#include "print_internal.h"
#include "structure.h"

/* File-local helpers (not in print.h / print_internal.h). */
static long TrajArr_Index(const int beadID, const int nFrameNumber, const int beadProp);
static void Write_ClusterDist(char* filename, long nGen);
static void Write_GyrTen(char* filename, long nGen);
static void FileIOUtil_WriteHeader_ForMCMove(const char* fileName);
static void FileIOUtil_WriteHeader_ForRDF(const char* strFileName);
static void FileIOUtil_WriteHeader_ForCOMDen(const char* strFileName);
static void FileIOUtil_WriteHeader_ForGyrRad(const char* strFileName);
static void FileIOUtil_WriteHeader_ForCLUS(const char* strFileName);
static void FileIOUtil_WriteHeader_ForMolClus(const char* strFileName);
static void FileIOUtil_WriteHeader_ForEnergy(const char* fileName);
static void FileIOUtil_Traj_Txt_AppendFrame_ToFile(const char* filename, const long nGen);
static void TrajUtil_SubselectDataFromBeadInfo(int** restrict sub_beads, int** restrict input_beads);
static void FileIOUtil_Traj_Bin_AppendFrame_ToFile(const char* filename, const long nGen);
static void Save_Trajectory(const long nGen, const long curFrame);
static void Write_Saved_Trajectory(char* filename, const int run_it);
static void Write_RDF_ComponentWise(char* filename, long nGen);
static void FileIO_WriteTo_TopFile(const char* filename);
static void Write_SysProp(char* filename);
static void FileIO_WriteTo_RDFFile_TotFromGLB(const int run_it);
static void FileIO_WriteTo_COMDenFile_TotFromGLB(const int run_it);
static void FileIO_WriteTo_ClusFile_TotFromGLB(const int run_it);
static void FileIO_WriteTo_MolClusFile_TotFromGLB(const int run_it);
static void FileIO_WriteTo_GyrRadFile_TotFromGLB(const int run_it);
static void FileIO_Write_TotalSysProp_TotFromGLB(const int run_it);
static void FileIOUtil_CreateFile_Overwrite(const char* fileName);
static void FileIOUtil_CreateFile_Binary_Overwrite(const char* fileName);
static void CopyData_All(const int run_it);
static void CopyData_RDF(const int run_it);
static void CopyData_COMDen(const int run_it);
static void CopyData_Clus(const int run_it);
static void FileIOUtil_AppendTo_GyrRadFile_ForRun(const int run_it);
static void FileIOUtil_AppendTo_ClusFile_ForRun(const int run_it);
static void FileIOUtil_AppendTo_MolClusFile_ForRun(const int run_it);
static void FileIOUtil_AppendTo_RDFFile_ForRun(const int run_it);
static void FileIOUtil_AppendTo_COMDenFile_ForRun(const int run_it);

/* initialize.h is main/parse-layer; forward-declare the one allocator we use
 * so I/O does not take a hard include edge upward (Phase 3.1). */
int** Create2DInt(const size_t xDim, const size_t yDim, const char* ArrName);

/// TrajArr_Index
/// \param beadID
/// \param nFrameNumber
/// \param beadProp
/// \return
static long TrajArr_Index(const int beadID, const int nFrameNumber, const int beadProp)
{
    return beadProp + BEADINFO_MAX * (beadID + g_system.tot_beads * nFrameNumber);
}

/// Write_ClusterDist - write the cluster histogram to a separate file. NOT USED
/// ANYMORE, but might find a use later Always appends to the file for this run.
/// Stopped using it because the IO load was slowing things down at our cluster
/// Would be a good way to gather proper statistics on clusters as the runs went
/// on \param filename \param nGen
static void Write_ClusterDist(char* filename, long nGen)
{
    FILE* fp;
    int i;
    if (nGen == -1)
    {
        fp = fopen(filename, "w"); // overwrite
    }
    else
    {
        fp = fopen(filename, "a");
    }

    if (nGen == -1)
    { // title
        fprintf(fp, "#Step followed by histogram\n");
    }
    else
    {
        fprintf(fp, "#%ld\n", nGen);
        for (i = 0; i <= g_system.tot_chains; i++)
        {
            fprintf(fp, "%ld\t", g_analysis.cycle.laClusHistList[i]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

/// Write_GyrTen - write the total gyration tensor to a file. Remember that we
/// only have 6 things in a symmetric 3x3 tensor. Always appends to the file for
/// this run. Stopped using it because the IO load was slowing things down at
/// our cluster Would be a good way to gather proper statistics on the shapes
/// and such as the runs went on \param filename \param nGen
static void Write_GyrTen(char* filename, long nGen)
{
    FILE* fp;
    int i;
    if (nGen == -1)
    {
        fp = fopen(filename, "w"); // overwrite
    }
    else
    {
        fp = fopen(filename, "a");
    }

    if (nGen == -1)
    { // title
        fprintf(fp, "#Step followed by Gyration Tensor Vals: xx yy zz xy yz xz\n");
    }
    else
    {
        fprintf(fp, "#%ld\n", nGen);
        for (i = 0; i < 7; i++)
        {
            if (i == 5)
            {
                continue;
            } // Not smart enough to make a better indexing function
            fprintf(fp, "%f\t", g_analysis.cycle.faGyrTensor[i]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

/// FileIOUtil_WriteHeader_ForMCMove - write the header for the MC move acceptance file.
/// \param fileName
static void FileIOUtil_WriteHeader_ForMCMove(const char* fileName)
{
    FILE* fp = fopen(fileName, "a");
    fprintf(fp, "#Steps, and Temp, are followed by (rej, acc) for each MC Move.\n");
    fprintf(
        fp,
        "#(nStep, T)             "
        "STROT                   "
        "LOCAL                   "
        "COLOCAL                 "
        "MTLOCAL                 "
        "SNAKE                   "
        "TRANS                   "
        "SM_CLST                 "
        "CLST                    "
        "PIVOT                   "
        "BRROT                   "
        "DBPVT                   "
        "PR_CLST                 "
        "\n");
    fclose(fp);
}

/// FileIO_WriteTo_MCMoveFile - writes the acceptance and rejection ratios of the various moves.
/// Just keeps appending to that file. Can be used to track the 'dynamics' during a simulation.
/// ALSO: Zeros out all the acceptance arrays. \param filename \param nGen \param fMCTemp
void FileIO_WriteTo_MCMoveFile(const char* filename, const long nGen, float const fMCTemp)
{
    FILE* fp = fopen(filename, "a+");
    int i;                                           // Iterator
    fprintf(fp, "%-10ld  %-10.2f  ", nGen, fMCTemp); // Step and Temp
    for (i = 1; i < MAX_MV; i++)
    {
        fprintf(fp, "%-10ld  %-10ld  ", g_mc.cycle.naMCAccepMat[0][i], g_mc.cycle.naMCAccepMat[1][i]);
    }
    fprintf(fp, "\n");

    for (i = 1; i < MAX_MV; i++)
    {
        g_mc.cycle.naMCAccepMat[0][i] = 0;
        g_mc.cycle.naMCAccepMat[1][i] = 0;
        // This way the print function will zero out the matrix every time
        // we print to a file!
    }

    fclose(fp);
}

static void FileIOUtil_WriteHeader_ForRDF(const char* strFileName)
{
    FILE* fp = fopen(strFileName, "w");

    fprintf(
        fp,
        "#Split RDFs. ALL-ALL; DIAGONALS (0-0, 1-1, ..) and then off-diagonals (0-1, 0-2, "
        "..., 1-2, ...)\n");

    fclose(fp);
}

static void FileIOUtil_WriteHeader_ForCOMDen(const char* strFileName)
{
    FILE* fp = fopen(strFileName, "w");

    fprintf(fp, "#Density Distribution from the COM outwards. Order is ChainType.\n");

    fclose(fp);
}

static void FileIOUtil_WriteHeader_ForGyrRad(const char* strFileName)
{
    FILE* fp = fopen(strFileName, "w");

    fprintf(fp, "#Gyration Radius Of ENTIRE System & BoxSize.\n");

    fclose(fp);
}

static void FileIOUtil_WriteHeader_ForCLUS(const char* strFileName)
{
    FILE* fp = fopen(strFileName, "w");

    fprintf(fp, "#Cluster size histograms. First is largest cluster, then n=1, n=2, and so on.\n");

    fclose(fp);
}

static void FileIOUtil_WriteHeader_ForMolClus(const char* strFileName)
{
    FILE* fp = fopen(strFileName, "w");

    fprintf(
        fp,
        "#Cluster size histograms, split by molecule type. Each row is a ChainType, and "
        "each column is the "
        "cluster size starting from 1.\n");

    fclose(fp);
}

/// FileIOUtil_WriteHeader_ForEnergy - write the header for the energy file.
/// \param fileName
static void FileIOUtil_WriteHeader_ForEnergy(const char* fileName)
{
    FILE* fp = fopen(fileName, "a+");
    fprintf(
        fp,
        "#"
        "STEP         "
        "TOT         "
        "OVLP        "
        "CONT        "
        "SC-SC       "
        "FSOL        "
        "T_IND       "
        "STIFF       "
        "BIAS        "
        "\n");
    fclose(fp);
}

/// FileIO_AppendEnergyTo_EnergyFile - write the decomposed energy of the system to a file.
/// Just appends to the file for this run.
/// \param fileNameStr
/// \param nGen
void FileIO_AppendEnergyTo_EnergyFile(const char* fileNameStr, const long nGen)
{
    FILE* fp = fopen(fileNameStr, "a");

    int i;
    fprintf(fp, "%-12ld", nGen);
    for (i = 0; i < (MAX_E); i++)
    {
        fprintf(fp, "  %-10.2E", g_energy.faCurrEn[i]);
    }
    fprintf(fp, "\n");

    fclose(fp);
}

/// FileIO_Trajectory_AppendFrame -
/// \param fileNameStr
/// \param run_it
/// \param nGen
void FileIO_Trajectory_AppendFrame(const char* fileNameStr, const int run_it, const long nGen)
{
    if (g_reports.nTrajMode == 1)
    {
        FileIOUtil_Traj_Bin_AppendFrame_ToFile(fileNameStr, 0);
    }
    else
    {
        if (run_it == -1) // Thermalization cycle
        {
            FileIOUtil_Traj_Txt_AppendFrame_ToFile(fileNameStr, nGen);
        }
        else
        {
            FileIOUtil_Traj_Txt_AppendFrame_ToFile(
                fileNameStr,
                nGen + g_mc.nMCStepsForTherm + run_it * g_mc.nMCStepsPerCycle);
        }
    }
}

/// FileIOUtil_Traj_Txt_AppendFrame_ToFile -  write a LAMMPS style trajectory file for the system.
/// Appends to the file for this run.
/// \param filename
/// \param nGen
static void FileIOUtil_Traj_Txt_AppendFrame_ToFile(const char* filename, const long nGen)
{
    // Writes the trajectory in LAMMPS format. To be viewed with VMD (hopefully). Read the LAMMPS
    // dump documentation for the actual format of the file
    FILE* fp = fopen(filename, "a"); // We will append to the file made for this run;

    int i; // Looping index
    fprintf(fp, "ITEM: TIMESTEP\n");
    fprintf(fp, "%ld\n", nGen); // Timestep

    fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fp, "%ld\n", (size_t)g_system.tot_beads); // Total atom number

    fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n"); // BCs are always periodic for now
    fprintf(
        fp,
        "0 %d\n0 %d\n0 %d\n",
        g_system.naBoxSize[0],
        g_system.naBoxSize[1],
        g_system.naBoxSize[2]); // Box dimensions

    fprintf(fp, "ITEM: ATOMS id type mol x y z bP\n"); // What we are printing

    for (i = 0; i < g_system.tot_beads; i++)
    {
        fprintf(
            fp,
            "%d %d %d %d %d %d %d\n",
            i,
            g_system.bead_info[i][BEAD_TYPE],
            g_system.bead_info[i][BEAD_CHAINID],
            g_system.bead_info[i][POS_X],
            g_system.bead_info[i][POS_Y],
            g_system.bead_info[i][POS_Z],
            g_system.bead_info[i][BEAD_FACE]);
    }

    fclose(fp);
}

/// TrajUtil_SubselectDataFromBeadInfo - Assuming that sub_beads is a 2D array with enough space
/// (nAtoms * 4), we extract POS_X, POS_Y, POS_Z, and BEAD_FACE from each bead. \param sub_beads
/// \param input_beads
static void TrajUtil_SubselectDataFromBeadInfo(int** restrict sub_beads, int** restrict input_beads)
{
    int i;

    for (i = 0; i < g_system.tot_beads; i++)
    {
        sub_beads[i][0] = input_beads[i][POS_X];
        sub_beads[i][1] = input_beads[i][POS_Y];
        sub_beads[i][2] = input_beads[i][POS_Z];
        sub_beads[i][3] = input_beads[i][BEAD_FACE];
    }
}

/// FileIOUtil_Traj_Bin_AppendFrame_ToFile - Appends the current frame in binary format. The format
/// is flat, and for every bead, in order of beadID, we write X, Y, Z, bP \param filename \param
/// nGen
static void FileIOUtil_Traj_Bin_AppendFrame_ToFile(const char* filename, const long nGen)
{
    const int nCrds = 4;

    int** subBeadInfo = Create2DInt(nCrds, g_system.tot_beads, "TrajBinSubArr");

    TrajUtil_SubselectDataFromBeadInfo(subBeadInfo, g_system.bead_info);

    FILE* fp = fopen(filename, "ab");

    const int nBeads = (int)g_system.tot_beads;

    fwrite(&nBeads, sizeof(int), 1, fp);
    fwrite(subBeadInfo[0], sizeof(int), nCrds * g_system.tot_beads, fp);
    fclose(fp);

    free(subBeadInfo[0]);
    free(subBeadInfo);
}

/// Save_Trajectory -  saves the current position onto the total array
/// \param nGen
/// \param curFrame
static void Save_Trajectory(const long nGen, const long curFrame)
{
    int i, j;
    long ar_idx;
    for (i = 0; i < g_system.tot_beads; i++)
    {
        for (j = 0; j < BEADINFO_MAX; j++)
        {
            ar_idx = TrajArr_Index(i, curFrame, j);
            g_reports.cycle.naTOTTRAJ_Arr[ar_idx] = g_system.bead_info[i][j];
        }
    }
}

/// Write_Saved_Trajectory -  write a LAMMPS style trajectory file for the
/// system. This funcion writes all the stored frames for this temperature cycle
/// \param filename
/// \param nGen
static void Write_Saved_Trajectory(char* filename, const int run_it)
{
    // Writes the trajectory in LAMMPS format. To be viewed with VMD
    // (Hopefully). Read the LAMMPS dump documentation for the actual formate of
    // the file
    FILE* fp;
    fp = fopen(filename, "a"); // We will append to the file made for this run

    int i, j, k; // Looping index
    int ar_idx;

    for (i = 0; i < g_reports.cycle.nTrajCurFrame; i++)
    {
        fprintf(fp, "ITEM: TIMESTEP\n");
        if (run_it >= 0)
        {
            fprintf(
                fp,
                "%ld\n",
                g_mc.nMCStepsForTherm + run_it * g_mc.nMCStepsPerCycle +
                    i * g_reports.naReportFreqs[REPORT_CONFIG]); // Timestep
        }
        else
        {
            fprintf(fp, "%ld\n", i * g_reports.naReportFreqs[REPORT_CONFIG]); // Timestep
        }

        fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
        fprintf(fp, "%ld\n", (size_t)g_system.tot_beads); // Total atom number

        fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n"); // BCs are always periodic for now
        fprintf(
            fp,
            "0 %d\n0 %d\n0 %d\n",
            g_system.naBoxSize[0],
            g_system.naBoxSize[1],
            g_system.naBoxSize[2]); // Box dimensions

        fprintf(fp, "ITEM: ATOMS id type mol x y z bP\n"); // What we are printing

        for (j = 0; j < g_system.tot_beads; j++)
        {
            fprintf(fp, "%d", j);
            ar_idx = TrajArr_Index(j, i, BEAD_TYPE);
            fprintf(fp, " %d", g_reports.cycle.naTOTTRAJ_Arr[ar_idx]);
            ar_idx = TrajArr_Index(j, i, BEAD_CHAINID);
            fprintf(fp, " %d", g_reports.cycle.naTOTTRAJ_Arr[ar_idx]);
            ar_idx = TrajArr_Index(j, i, POS_X);
            fprintf(fp, " %d", g_reports.cycle.naTOTTRAJ_Arr[ar_idx]);
            ar_idx = TrajArr_Index(j, i, POS_Y);
            fprintf(fp, " %d", g_reports.cycle.naTOTTRAJ_Arr[ar_idx]);
            ar_idx = TrajArr_Index(j, i, POS_Z);
            fprintf(fp, " %d", g_reports.cycle.naTOTTRAJ_Arr[ar_idx]);
            ar_idx = TrajArr_Index(j, i, BEAD_FACE);
            fprintf(fp, " %d", g_reports.cycle.naTOTTRAJ_Arr[ar_idx]);
            fprintf(fp, "\n");
        }
    }

    fclose(fp);
}

/// Write_RDF_ComponentWise - old implementation of printing the RDF, component
/// by component. Always appends to the file for this run. Stopped using it
/// because the IO load was slowing things down at our cluster Would be a good
/// way to gather proper statistics on clusters as the runs went on.
/// TODO: update this to work with the new indexing
/// \param filename
/// \param nGen
static void Write_RDF_ComponentWise(char* filename, long nGen)
{
    FILE* fp;
    if (nGen == -1)
    {
        fp = fopen(filename, "w"); // overwrite
    }
    else
    {
        fp = fopen(filename, "a");
    }
    int i, k;
    if (nGen == -1)
    { // title
        fprintf(fp, "#Row-by-row split RDF. dr = 1/4\n");
    }
    else
    {
        for (k = 0; k < g_analysis.nRDF_TotComps; k++)
        {
            for (i = 0; i < g_analysis.nRDF_TotBins; i++)
            {
                fprintf(fp, "%.5Lf\t", g_analysis.cycle.ldaRDF_Arr[RDFArr_Index(0, k, i)]);
            }
            fprintf(fp, "\n");
        }
    }

    fclose(fp);
}

/// FileIO_WriteTo_TopFile - write a LAMMPS DATA file, which can be read in VMD to store
/// topological information for the system Only writes the topology once at the
/// beginning of the runs.
/// \param filename
static void FileIO_WriteTo_TopFile(const char* filename)
{
    /*
    Writes a topology file for VMD. Since the trajectory is saved in the LAMMPS
    format, the topology file is also in the LAMMPS format. The format for this
    file is one which is used as input data for LAMMPS. This function copies the
    format given in read_data for LAMMPS:
    https://lammps.sandia.gov/doc/read_data.html
    */
    FILE* fp;
    fp = fopen(filename, "w"); // Just overwrite over if last file exists.
    int i, j, k;               // Loop iterators.
    int numBonds;              // Used to count total number of bonds!
    printf("Writing the topology file!\n");
    fprintf(fp, "LAMMPS Description\n");                       // The file must start with this.
    fprintf(fp, "\n");                                         // Empty line.
    fprintf(fp, "\t%ld\tatoms\n", (size_t)g_system.tot_beads); // Listing total number of atoms
    numBonds = 0;
    for (i = 0; i < g_system.tot_beads; i++)
    {
        for (j = 0; j < MAX_BONDS; j++)
        {
            if (g_system.topo_info[i][j] != -1)
            { // There is a bond between i and g_system.topo_info[i][j]
                numBonds++;
            }
        }
    }

    fprintf(fp, "\t%d\tbonds\n", numBonds);               // Listing number of bonds
    fprintf(fp, "\t0\tangles\n");                         // Systems don't have angle depenece yet
    fprintf(fp, "\t0\tdihedrals\n");                      // Systems don't have dihedrals  yet
    fprintf(fp, "\t0\timpropers\n");                      // Systems don't have imporopers yet
    fprintf(fp, "\n");                                    // Empty line.
    fprintf(fp, "%d\tatom types\n", g_energy.nBeadTypes); // This many bead-types
    fprintf(fp, "1\tbond types\n");                       // System can have multiple bond-lengths
    fprintf(fp, "0\tangle types\n");                      // Systems don't have any angular forces yet
    fprintf(fp, "\n");                                    // Empty line.
    fprintf(fp, " 0 %d xlo xhi\n", g_system.naBoxSize[0]);
    fprintf(fp, " 0 %d ylo yhi\n", g_system.naBoxSize[1]);
    fprintf(fp, " 0 %d zlo zhi\n", g_system.naBoxSize[2]);
    fprintf(fp, "\n");       // Empty line.
    fprintf(fp, "Masses\n"); // These don't really mean anything
    fprintf(fp, "\n");       // Empty line.
    for (i = 0; i < g_energy.nBeadTypes; i++)
    {
        fprintf(fp, "%d\t1.0\n", i);
    }

    fprintf(fp, "\n");             // Empty line.
    fprintf(fp, "Atoms # bond\n"); // Signifying the beginning of atom coordinates.
    fprintf(fp, "\n");             // Empty line.
    for (i = 0; i < g_system.tot_beads; i++)
    {
        fprintf(
            fp,
            "%d %d %d %d %d %d\n",
            i,
            g_system.bead_info[i][BEAD_CHAINID],
            g_system.bead_info[i][BEAD_TYPE],
            g_system.bead_info[i][POS_X],
            g_system.bead_info[i][POS_Y],
            g_system.bead_info[i][POS_Z]);
    } // Done with the coordinates
    fprintf(fp, "\n");      // Empty line.
    fprintf(fp, "Bonds\n"); // Signifying the beginning of atom coordinates.
    fprintf(fp, "\n");      // Empty line.
    k = 0;                  // This guy counts bondIDs
    for (i = 0; i < g_system.tot_beads; i++)
    {
        for (j = 0; j < MAX_BONDS; j++)
        {
            if (g_system.topo_info[i][j] != -1)
            { // There is a bond between i and g_system.topo_info[i][j]
                fprintf(
                    fp,
                    "%d %d %d %d\n",
                    k,
                    g_system.linker_len[i][j] / (int)g_system.fLinkerLength,
                    i,
                    g_system.topo_info[i][j]);
                k++;
            }
        }
    } // Done with the coordinates
    fprintf(fp, "\n"); // Empty line.

    fclose(fp);
}

/// Write_SysProp - writes the RDF, CLUS and GyrTen avg rrays to a file. OLD implementation that
/// should not be used yet Also, TODO: update this to work with the new indexing \param filename
static void Write_SysProp(char* filename)
{
    FILE* fp;
    int i, j;
    fp = fopen(filename, "w");
    fprintf(fp, "#Contains various averaged quatities.\n");
    // Gyration Radius
    fprintf(
        fp,
        "#Total Rg\n%f\t%f\n#Cluster Hist\n",
        g_analysis.cycle.faSysGyrRad / (float)g_analysis.cycle.nTotGyrRadCounter,
        (float)g_system.naBoxSize[0] / 2.);
    // Cluster Histogram
    fprintf(fp, "%f\t", (float)g_analysis.cycle.nLargestClusterRightNow / (float)g_analysis.cycle.nTotClusCounter);
    for (i = 1; i <= g_system.tot_chains; i++)
    {
        fprintf(fp, "%f\t", (float)g_analysis.cycle.laClusHistList[i] / (float)g_analysis.cycle.nTotClusCounter);
    }
    // Split RDFs
    fprintf(fp, "\n#Split RDFs. ALL-ALL; DIAGONALS and then from 0 onwards \n");
    for (j = 0; j < g_analysis.nRDF_TotComps; j++)
    {
        for (i = 0; i < g_analysis.nRDF_TotBins; i++)
        {
            // fprintf(fp, "%LE\t", ldRDF_ARR[j][i] / (float)g_analysis.cycle.nTotRDFCounter);
            fprintf(
                fp,
                "%LE\t",
                g_analysis.cycle.ldaRDF_Arr[RDFArr_Index(0, j, i)] / (float)g_analysis.cycle.nTotRDFCounter);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n#Done");
}

/// FileIO_WriteTo_RDFFile_TotFromGLB - writes the RDFTot file.
/// \param run_it Number of temperature cycles.
static void FileIO_WriteTo_RDFFile_TotFromGLB(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_RDF.dat",
            g_reports.strReportPrefix); // Name Of the RDF Files
    fp = fopen(dumFile, "w");
    fprintf(fp, "#Split RDFs. ALL-ALL; DIAGONALS and then from 0-0, 0-1, and onwards \n");
    for (i = 0; i < run_it; i++)
    {
        fprintf(fp, "#Run_Cycle = %d\n", i);
        for (j = 0; j < g_analysis.nRDF_TotComps; j++)
        {
            for (k = 0; k < g_analysis.nRDF_TotBins; k++)
            {
                fprintf(fp, "%LE\t", g_analysis.totals.ldaTOTRDF_Arr[RDFArr_Index(i, j, k)]);
            }
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_WriteTo_COMDenFile_TotFromGLB - writes the COMDenTot file.
/// \param run_it Number of temperature cycles.
static void FileIO_WriteTo_COMDenFile_TotFromGLB(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_COMDen.dat",
            g_reports.strReportPrefix); // Name Of the COM Density Distribution
    fp = fopen(dumFile, "w");
    fprintf(fp, "#Density distribution from the COM outwards. Order is ChainType.\n");
    for (i = 0; i < run_it; i++)
    {
        fprintf(fp, "#Run_Cycle = %d\n", i);
        for (j = 0; j < g_analysis.nRadDen_TotComps; j++)
        {
            for (k = 0; k < g_analysis.nRDF_TotBins; k++)
            {
                fprintf(fp, "%LE\t", g_analysis.totals.ldaTOTRadDen_Arr[RadDenArr_Index(i, j, k)]);
            }
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_WriteTo_ClusFile_TotFromGLB - writes the CLUS file.
/// \param run_it Number of temperature cycles.
static void FileIO_WriteTo_ClusFile_TotFromGLB(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_CLUS.dat",
            g_reports.strReportPrefix); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "w");
    fprintf(
        fp,
        "#Cluster Histograms: 1st column is largest cluster, and then clusters of size 1, "
        "2, and so on\n");
    for (i = 0; i < run_it; i++)
    {
        fprintf(fp, "#Run_Cycle = %d\n", i);
        for (k = 0; k <= g_system.tot_chains; k++)
        {
            fprintf(fp, "%LE\t", g_analysis.totals.ldaTOTCLUS_Arr[i][k]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_WriteTo_MolClus - writes the MolClus file
/// \param run_it Number of temperature cycles.
static void FileIO_WriteTo_MolClusFile_TotFromGLB(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_MolClus.dat",
            g_reports.strReportPrefix); // Name Of the COM Density Distribution
    fp = fopen(dumFile, "w");
    fprintf(
        fp,
        "#Cluster Histograms: 1st column is largest cluster, and then clusters of size 1, "
        "2, and so on. Each "
        "row is a different moltype\n");
    for (i = 0; i < run_it; i++)
    {
        fprintf(fp, "#Run_Cycle = %d\n", i);
        for (j = 0; j < g_system.tot_chain_types; j++)
        {
            for (k = 0; k < g_system.tot_chains; k++)
            {
                fprintf(fp, "%LE\t", g_analysis.totals.ldaTOTMOLCLUS_Arr[MolClusArr_Index(i, j, k)]);
            }
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_WriteTo_GyrRad - writes the GR file.
/// \param run_it Number of temperature cycles.
static void FileIO_WriteTo_GyrRadFile_TotFromGLB(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_GR.dat",
            g_reports.strReportPrefix); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "w");
    fprintf(
        fp,
        "#Cluster Histograms: 1st column is largest cluster, and then "
        "clusters of size 1, 2, and so on\n");
    for (i = 0; i < run_it; i++)
    {
        fprintf(fp, "#Run_Cycle = %d\n", i);
        fprintf(fp, "%LE\t", g_analysis.totals.ldaTOTRg_Arr[i][0]);
        fprintf(fp, "%LE\n", g_analysis.totals.ldaTOTRg_Arr[i][1]);
    }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_Write_TotalSysProp_TotFromGLB - writes the RDF, CLUS and GyrTen arrays to their
/// respective files. These arrays store the data over the course of an ENTIRE
/// run, or over all the temp cycles.
/// \param run_it
static void FileIO_Write_TotalSysProp_TotFromGLB(const int run_it)
{
    /* This function writes one large file with all the averaged values from
     * each run_cycle. run_it let's the function know how many cycles to write.
     * As opposed to Write_SysProp, each of the relevant averaged parameters
     * will be in their appropriately named files. The naming convention shall
     * be: $REPORTNAME_*.dat where * will be GR, CLUS and RDF for the different
     * measured quantities. Within each file, each run_cycle shall be written
     * one-by-one in the style of Write_SysProp
     */
    if (g_reports.naReportFreqs[REPORT_RDFTOT] != 0)
    {
        FileIO_WriteTo_RDFFile_TotFromGLB(run_it);
    }

    if (g_reports.naReportFreqs[REPORT_COMDEN] != 0)
    {
        FileIO_WriteTo_COMDenFile_TotFromGLB(run_it);
    }

    if (g_reports.naReportFreqs[REPORT_NETWORK] != 0)
    {
        FileIO_WriteTo_ClusFile_TotFromGLB(run_it);
        FileIO_WriteTo_MolClusFile_TotFromGLB(run_it);
        FileIO_WriteTo_GyrRadFile_TotFromGLB(run_it);
    }
}

/// FileIOUtil_CreateFile_Overwrite - Creates a new overwritten file with the given name.
/// \param fileName
static void FileIOUtil_CreateFile_Overwrite(const char* fileName)
{
    FILE* fp = fopen(fileName, "w+");
    fclose(fp);
}

/// FileIOUtil_CreateFile_Overwrite - Creates a new overwritten file with the given name.
/// \param fileName
static void FileIOUtil_CreateFile_Binary_Overwrite(const char* fileName)
{
    FILE* fp = fopen(fileName, "wb");
    fclose(fp);
}

/// FileIO_CreateRunningDataFiles - Creates the necessary files that are continuously written to
/// over the course of a simulation.
void FileIO_CreateRunningDataFiles(void)
{
    // Trajectory
    if (g_reports.naReportFreqs[REPORT_CONFIG])
    {
        sprintf(g_reports.strFileTraj, "%s_topo.lammpstrj",
                g_reports.strReportPrefix);            // Name of the topology file
        FileIO_WriteTo_TopFile(g_reports.strFileTraj); // Write the topology file. Only need to write once

        if (g_reports.nTrajMode == 1)
        {
            sprintf(g_reports.strFileTraj, "%s_EQ_trj.lassi", g_reports.strReportPrefix);
            FileIOUtil_CreateFile_Binary_Overwrite(g_reports.strFileTraj); // Opens a new binary file.
        }
        else
        {
            sprintf(
                g_reports.strFileTraj,
                "%s_EQ_trj.lammpstrj",
                g_reports.strReportPrefix);                         // Naming convention for trajectory files.
            FileIOUtil_CreateFile_Overwrite(g_reports.strFileTraj); // This opens a new trajectory file; each run_it
                                                                    // will have its own
        }
    }

    // Energy
    if (g_reports.naReportFreqs[REPORT_ENERGY])
    {
        sprintf(g_reports.strFileEnergy, "%s_energy.dat", g_reports.strReportPrefix);
        FileIOUtil_CreateFile_Overwrite(
            g_reports.strFileEnergy); // Open a new energy file; each run_it will have its own
        FileIOUtil_WriteHeader_ForEnergy(g_reports.strFileEnergy);
    }

    // MC Move Acceptance
    if (g_reports.naReportFreqs[REPORT_MCMOVE])
    {
        sprintf(g_reports.strFileMCMove, "%s_mcmove.dat", g_reports.strReportPrefix);
        FileIOUtil_CreateFile_Overwrite(
            g_reports.strFileMCMove); // Open a new MCInfo file; each run_it will have its own
        FileIOUtil_WriteHeader_ForMCMove(g_reports.strFileMCMove);
    }

    char strFileName[512];

    if (g_reports.naReportFreqs[REPORT_COMDEN])
    {
        sprintf(strFileName, "%s_COMDen.dat", g_reports.strReportPrefix);
        FileIOUtil_CreateFile_Overwrite(strFileName);
        FileIOUtil_WriteHeader_ForCOMDen(strFileName);
    }

    if (g_reports.naReportFreqs[REPORT_NETWORK])
    {
        sprintf(strFileName, "%s_CLUS.dat", g_reports.strReportPrefix);
        FileIOUtil_CreateFile_Overwrite(strFileName);
        FileIOUtil_WriteHeader_ForCLUS(strFileName);

        sprintf(strFileName, "%s_MolClus.dat", g_reports.strReportPrefix);
        FileIOUtil_CreateFile_Overwrite(strFileName);
        FileIOUtil_WriteHeader_ForMolClus(strFileName);

        sprintf(strFileName, "%s_GR.dat", g_reports.strReportPrefix);
        FileIOUtil_CreateFile_Overwrite(strFileName);
        FileIOUtil_WriteHeader_ForGyrRad(strFileName);
    }

    if (g_reports.naReportFreqs[REPORT_RDFTOT])
    {
        sprintf(strFileName, "%s_RDF.dat", g_reports.strReportPrefix);
        FileIOUtil_CreateFile_Overwrite(strFileName);
        FileIOUtil_WriteHeader_ForRDF(strFileName);
    }
}

void FileIO_PrintCrashSnapshot(void)
{
    FILE* fp = fopen("crash_snapshot.txt", "w+");

    int i, x, y, z, bP;
    fputs("#CRASH SNAPSHOT\n#X, Y, Z, bondPartner\n", fp);
    for (i = 0; i < g_system.tot_beads; i++)
    {
        x = g_system.bead_info[i][POS_X];
        y = g_system.bead_info[i][POS_Y];
        z = g_system.bead_info[i][POS_Z];
        bP = g_system.bead_info[i][BEAD_FACE];
        fprintf(fp, "%d %d %d %d\n", x, y, z, bP);
    }

    fclose(fp);
}

/// FileIOUtil_PreCycle_Init
/// \param run_it
void FileIOUtil_PreCycle_Init(const int run_it)
{
    if (g_reports.naReportFreqs[REPORT_CONFIG])
    {
        if (g_reports.nTrajMode == 1)
        {
            sprintf(g_reports.strFileTraj, "%s_trj.lassi", g_reports.strReportPrefix);
            if (run_it == 0)
            {
                FileIOUtil_CreateFile_Binary_Overwrite(g_reports.strFileTraj); // Opens a new binary file.
            }
        }
        else
        {
            sprintf(g_reports.strFileTraj, "%s_trj.lammpstrj", g_reports.strReportPrefix);
            if (run_it == 0)
            {
                FileIOUtil_CreateFile_Overwrite(g_reports.strFileTraj);
            }
        }
    }

    if (g_reports.naReportFreqs[REPORT_ENERGY])
    {
        sprintf(g_reports.strFileEnergy, "%s_%d_energy.dat", g_reports.strReportPrefix, run_it);
        FileIOUtil_CreateFile_Overwrite(
            g_reports.strFileEnergy); // Open a new Energy file; each run_it will have its own
        FileIOUtil_WriteHeader_ForEnergy(g_reports.strFileEnergy);
    }
    if (g_reports.naReportFreqs[REPORT_MCMOVE])
    {
        sprintf(g_reports.strFileMCMove, "%s_%d_mcmove.dat", g_reports.strReportPrefix, run_it);
        FileIOUtil_CreateFile_Overwrite(
            g_reports.strFileMCMove); // Open a new MCInfo file; each run_it will have its own
        FileIOUtil_WriteHeader_ForMCMove(g_reports.strFileMCMove);
    }
}

/// FileIO_PostCycle_WriteSystemRestart
/// \param run_it
void FileIO_PostCycle_WriteSystemRestart(const int run_it)
{
    sprintf(
        g_reports.strFileTraj,
        "%s_%d_restart.lammpstrj",
        g_reports.strReportPrefix,
        run_it); // Naming convention for trajectory files.
    FileIOUtil_CreateFile_Overwrite(g_reports.strFileTraj);
    FileIOUtil_Traj_Txt_AppendFrame_ToFile(
        g_reports.strFileTraj,
        g_mc.nMCStepsForTherm + (run_it + 1) * g_mc.nMCStepsPerCycle);
}

/// FileIO_WriteRestart_ForThermalization
void FileIO_WriteRestart_ForThermalization(void)
{
    if (g_mc.nMCStepsForTherm)
    {
        sprintf(
            g_reports.strFileTraj,
            "%s_EQ_restart.lammpstrj",
            g_reports.strReportPrefix); // Naming convention for trajectory files.
        FileIOUtil_CreateFile_Overwrite(g_reports.strFileTraj);
        FileIOUtil_Traj_Txt_AppendFrame_ToFile(g_reports.strFileTraj, g_mc.nMCStepsForTherm);
    }
}

/// CopyData_All - copies data from run_it specific data arrays to the overall
/// global data arrays that are printed later. Also note that the averaging, or
/// dividing by the frequency of acquisitions, occurs here.
/// \param run_it
static void CopyData_All(const int run_it)
{
    if (g_reports.naReportFreqs[REPORT_RDFTOT])
    {
        CopyData_RDF(run_it);
    }
    if (g_reports.naReportFreqs[REPORT_COMDEN])
    {
        CopyData_COMDen(run_it);
    }
    if (g_reports.naReportFreqs[REPORT_NETWORK])
    {
        CopyData_Clus(run_it);
    }
}

/// CopyData_RDF: Copies g_analysis.cycle.ldaRDF_Arr into g_analysis.totals.ldaTOTRDF_Arr
/// \param run_it: which run cycle we are on. Should be >= 0.
static void CopyData_RDF(const int run_it)
{
    int i, j;
    const long double normFactor = 1.0 / (long double)g_analysis.cycle.nTotRDFCounter;
    for (i = 0; i < g_analysis.nRDF_TotComps; i++)
    {
        for (j = 0; j < g_analysis.nRDF_TotBins; j++)
        {
            g_analysis.totals.ldaTOTRDF_Arr[RDFArr_Index(run_it, i, j)] =
                g_analysis.cycle.ldaRDF_Arr[RDFArr_Index(0, i, j)] * normFactor;
        }
    }
}

/// CopyData_RDF: Copies g_analysis.cycle.ldaRadDen_Arr into g_analysis.totals.ldaTOTRadDen_Arr
/// \param run_it: which run cycle we are on. Should be >= 0.
static void CopyData_COMDen(const int run_it)
{
    int i, j;
    const long double normFactor = 1.0 / (long double)g_analysis.cycle.nTotRadDenCounter;
    for (i = 0; i < g_analysis.nRadDen_TotComps; i++)
    {
        for (j = 0; j < g_analysis.nRDF_TotBins; j++)
        {
            g_analysis.totals.ldaTOTRadDen_Arr[RadDenArr_Index(run_it, i, j)] =
                g_analysis.cycle.ldaRadDen_Arr[RadDenArr_Index(0, i, j)] * normFactor;
        }
    }
}

/// CopyData_Clus - Copies ldMOLClUS_Arr into g_analysis.totals.ldaTOTMOLCLUS_Arr. Also copies
/// g_analysis.cycle.laClusHistList into g_analysis.totals.ldaTOTCLUS_Arr. And stores the Gyration radius. \param
/// run_it: which run cycle we are on. Should be >= 0.
static void CopyData_Clus(const int run_it)
{
    int i, j;
    const long double normFactor = 1.0 / (long double)g_analysis.cycle.nTotClusCounter;
    g_analysis.totals.ldaTOTCLUS_Arr[run_it][0] = (long double)g_analysis.cycle.nLargestClusterRightNow * normFactor;
    for (i = 1; i <= g_system.tot_chains; i++)
    {
        g_analysis.totals.ldaTOTCLUS_Arr[run_it][i] = (long double)g_analysis.cycle.laClusHistList[i] * normFactor;
    }
    for (i = 0; i < g_system.tot_chains; i++)
    {
        for (j = 0; j < g_system.tot_chain_types; j++)
        {
            g_analysis.totals.ldaTOTMOLCLUS_Arr[MolClusArr_Index(run_it, j, i)] =
                g_analysis.cycle.ldaMOLCLUS_Arr[MolClusArr_Index(0, j, i)] * normFactor;
        }
    }
    g_analysis.totals.ldaTOTRg_Arr[run_it][0] =
        (long double)g_analysis.cycle.faSysGyrRad / (long double)g_analysis.cycle.nTotGyrRadCounter;
    g_analysis.totals.ldaTOTRg_Arr[run_it][1] = (long double)g_system.naBoxSize[0] / 2.0;
}

///
/// \param run_it
static void FileIOUtil_AppendTo_GyrRadFile_ForRun(const int run_it)
{
    FILE* fp;
    char dumFile[512];
    sprintf(dumFile, "%s_GR.dat", g_reports.strReportPrefix); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "a");

    fprintf(fp, "#Run_Cycle = %d; Avg Over %d Samples\n", run_it, g_analysis.cycle.nTotClusCounter);
    fprintf(fp, "%LE\t", (long double)g_analysis.cycle.faSysGyrRad / (long double)g_analysis.cycle.nTotGyrRadCounter);
    fprintf(fp, "%LE\n", (long double)g_system.naBoxSize[0] / 2.0);

    fclose(fp);
}

///
/// \param run_it
static void FileIOUtil_AppendTo_ClusFile_ForRun(const int run_it)
{
    FILE* fp;
    int i;
    char dumFile[512];
    sprintf(dumFile, "%s_CLUS.dat",
            g_reports.strReportPrefix); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "a");
    fprintf(fp, "#Run_Cycle = %d; Avg Over %d Samples\n", run_it, g_analysis.cycle.nTotClusCounter);
    const double normFactor = 1.0 / (double)g_analysis.cycle.nTotClusCounter;

    g_analysis.cycle.laClusHistList[0] = g_analysis.cycle.nLargestClusterRightNow;

    for (i = 0; i <= g_system.tot_chains; i++)
    {
        fprintf(fp, "%E\t", (double)g_analysis.cycle.laClusHistList[i] * normFactor);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

///
/// \param run_it
static void FileIOUtil_AppendTo_MolClusFile_ForRun(const int run_it)
{
    FILE* fp;
    int i, j;
    char dumFile[512];
    sprintf(dumFile, "%s_MolClus.dat",
            g_reports.strReportPrefix); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "a");
    fprintf(fp, "#Run_Cycle = %d; Avg Over %d Samples\n", run_it, g_analysis.cycle.nTotClusCounter);

    const long double normFactor = 1.0 / (long double)g_analysis.cycle.nTotClusCounter;

    for (j = 0; j < g_system.tot_chain_types; j++)
    {
        for (i = 0; i < g_system.tot_chains; i++)
        {
            fprintf(fp, "%LE\t", g_analysis.cycle.ldaMOLCLUS_Arr[MolClusArr_Index(0, j, i)] * normFactor);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

///
/// \param run_it
static void FileIOUtil_AppendTo_RDFFile_ForRun(const int run_it)
{
    FILE* fp;
    int i, j;
    char dumFile[512];
    sprintf(dumFile, "%s_RDF.dat", g_reports.strReportPrefix); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "a");
    fprintf(fp, "#Run_Cycle = %d; Avg Over %d Samples\n", run_it, g_analysis.cycle.nTotRDFCounter);

    const long double normFactor = 1.0 / (long double)g_analysis.cycle.nTotRDFCounter;

    for (j = 0; j < g_analysis.nRDF_TotComps; j++)
    {
        for (i = 0; i < g_analysis.nRDF_TotBins; i++)
        {
            fprintf(fp, "%LE\t", g_analysis.cycle.ldaRDF_Arr[RDFArr_Index(0, j, i)] * normFactor);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

///
/// \param run_it
static void FileIOUtil_AppendTo_COMDenFile_ForRun(const int run_it)
{
    FILE* fp;
    int i, j;
    char dumFile[512];
    sprintf(dumFile, "%s_COMDen.dat", g_reports.strReportPrefix);
    fp = fopen(dumFile, "a");
    fprintf(fp, "#Run_Cycle = %d; Avg Over %d Samples\n", run_it, g_analysis.cycle.nTotRadDenCounter);

    const long double normFactor = 1.0 / (long double)g_analysis.cycle.nTotRadDenCounter;

    for (j = 0; j < g_analysis.nRadDen_TotComps; j++)
    {
        for (i = 0; i < g_analysis.nRDF_TotBins; i++)
        {
            fprintf(fp, "%LE\t", g_analysis.cycle.ldaRadDen_Arr[RadDenArr_Index(0, j, i)] * normFactor);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

///
/// \param run_it
void FileIO_PostCycle_WriteCycleAvgData(const int run_it)
{
    if (g_reports.naReportFreqs[REPORT_RDFTOT])
    {
        FileIOUtil_AppendTo_RDFFile_ForRun(run_it);
    }
    if (g_reports.naReportFreqs[REPORT_COMDEN])
    {
        FileIOUtil_AppendTo_COMDenFile_ForRun(run_it);
    }
    if (g_reports.naReportFreqs[REPORT_NETWORK])
    {
        FileIOUtil_AppendTo_ClusFile_ForRun(run_it);
        FileIOUtil_AppendTo_MolClusFile_ForRun(run_it);
        FileIOUtil_AppendTo_GyrRadFile_ForRun(run_it);
    }
}
