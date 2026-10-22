#include "mcmove.h"
#include "cluster.h"
#include "energy.h"
#include "structure.h"

/// MC_Step - given the MC move frequencies, picks which move to try.
/// \param fMCTemp
/// \return bAccept - the acceptance state where MCMove=bAccept/12 and
/// MCAccept=bAccept%2.
int MC_Step(float fMCTemp)
{
    int mode = 0; // Which move to do
    int nAccept;  // Used in MC steps
    float prob = (float) rand() / (float) RAND_MAX;
    int i; // Internal iterator, and int.
    for (i = 0; i < MAX_MV; i++)
        {
            if (prob < g_mc.faMCFreq[i])
                { // Pick this move!
                    mode = i;
                    break;
                }
        }
    nAccept = 0;

    // Performing the appropriate move.
    switch (mode)
        {
            // Translation
            case MV_TRANS:
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_Trans(i, fMCTemp);
                break;

                // Cluster translation
            case MV_CLSTR:
                nAccept = Move_Clus_Network(fMCTemp);
                break;

                // Cluster translation iff ClusSize <= 5
            case MV_SMCLSTR:
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_SmallClus_Network(i, fMCTemp);
                break;

                // Change Rotational State
            case MV_STROT:
                i       = rand() % tot_beads_glb;
                nAccept = Move_Rot(i, fMCTemp);
                break;

                // Local Move
            case MV_LOCAL:
                i       = rand() % tot_beads_glb;
                nAccept = Move_Local(i, fMCTemp);
                break;

                // Slithering Snake
            case MV_SNAKE:
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_Snake(i, fMCTemp);
                break;

                // Double Pivot
            case MV_DBPVT:
                i       = rand() % tot_beads_glb; // Pick a random bead
                nAccept = Move_DbPvt(i, fMCTemp);
                break;

                // Co-Local
            case MV_COLOCAL:
                i       = rand() % tot_beads_glb; // Pick a random bead
                nAccept = Move_CoLocal(i, fMCTemp);
                break;

                // Shake Move
            case MV_MTLOCAL:
                i       = rand() % tot_beads_glb;
                nAccept = Move_MultiLocal(i, fMCTemp);
                break;

                // Pivot
            case MV_PIVOT:
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_Pivot(i, fMCTemp);
                break;

                // Branched Rotation
            case MV_BRROT:
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_BranchedRot(i, fMCTemp);
                break;
            case MV_PR_SMCLSTR:
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_SmallClus_Proximity(i, fMCTemp);
                break;

            case MV_PR_CLSTR:
                nAccept = Move_Clus_Proximity(fMCTemp);
            default:
                nAccept = 0;
                break;
        }

    g_mc.naMCAccepMat[nAccept][mode]++; // Just adding which move got accepted/rejected.
    return mode * 12 + nAccept;
}

/// MC_Step_Therm - given the MC move frequencies, picks non-anisotropic
/// variants of the moves. \param fMCTemp \return bAccept - the acceptance state
/// where MCMove=bAccept/12 and MCAccept=bAccept%2.
int MC_Step_Therm(float fMCTemp)
{
    int mode = 0; // Which move to do
    int nAccept;  // Used in MC steps
    float prob = (float) rand() / (float) RAND_MAX;
    int i; // Internal iterator, and int.
    for (i = 0; i < MAX_MV; i++)
        {
            if (prob < g_mc.faMCFreq[i])
                { // Pick this move!
                    mode = i;
                    break;
                }
        }
    // Start by assuming failure.
    nAccept = 0;

    // Performing the appropriate move.
    switch (mode)
        {
            // translation of a random chain
            case MV_TRANS:
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_Trans_Equil(i, fMCTemp);
                break;

                // cluster translation: moves largest cluster to another spot.
            case MV_CLSTR:                         // In equil, just a translation of the chain
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_Trans_Equil(i, fMCTemp);
                break;

            case MV_SMCLSTR:                       // In equil, just a translation of the chain
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_Trans_Equil(i, fMCTemp);
                break;

                // face change
            case MV_STROT: // In equil, just a local move
                i       = rand() % tot_beads_glb;
                nAccept = Move_Local_Equil(i, fMCTemp);
                break;

                // local moves
            case MV_LOCAL:
                i       = rand() % tot_beads_glb;
                nAccept = Move_Local_Equil(i, fMCTemp);
                break;

                // slithering snake
            case MV_SNAKE:
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_Snake_Equil(i, fMCTemp);
                break;

                // double pivot
            case MV_DBPVT:
                i       = rand() % tot_beads_glb; // Pick a random bead
                nAccept = Move_DbPvt(i, fMCTemp);
                break;

                // co-local
            case MV_COLOCAL: // In equil, just a local move
                i       = rand() % tot_beads_glb;
                nAccept = Move_Local_Equil(i, fMCTemp);
                break;

                // shake
            case MV_MTLOCAL:
                i       = rand() % tot_beads_glb;
                nAccept = Move_MultiLocal_Equil(i, fMCTemp);
                break;

                // pivot
            case MV_PIVOT:
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_Pivot_Equil(i, fMCTemp);
                break;

                // branched rotate
            case MV_BRROT:
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_BranchedRot_Equil(i, fMCTemp);
                break;

            case MV_PR_SMCLSTR:
                // In equil, just a translation of the chain
                i       = rand() % tot_beads_glb;  // Pick random bead
                i       = bead_info_glb[i][BEAD_CHAINID]; // Get the chainID
                nAccept = Move_Trans_Equil(i, fMCTemp);
                break;

            default:
                nAccept = 0;
                break;
        }

    g_mc.naMCAccepMat[nAccept][mode]++; // Just adding which move got accepted/rejected.
    return mode * 12 + nAccept;
}
