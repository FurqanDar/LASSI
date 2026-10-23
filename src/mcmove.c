#include "mcmove_internal.h"

/*
 * MV_* strategy tables. Each entry is a thin runner that preserves the old
 * switch's bead/chain selection and Move_* call. MC_Step / MC_Step_Therm only
 * pick a mode from faMCFreq and invoke the table. Moves that already match
 * McMoveRunner (e.g. Move_Clus_Network) are registered directly.
 */

typedef int (*McMoveRunner)(float fMCTemp);

static int McRun_Reject(float fMCTemp)
{
    (void)fMCTemp;
    return 0;
}

static int McRun_PickBead_Rot(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    return Move_Rot(i, fMCTemp);
}

static int McRun_PickBead_Local(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    return Move_Local(i, fMCTemp);
}

static int McRun_PickBead_Local_Equil(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    return Move_Local_Equil(i, fMCTemp);
}

static int McRun_PickBead_CoLocal(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    return Move_CoLocal(i, fMCTemp);
}

static int McRun_PickBead_MultiLocal(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    return Move_MultiLocal(i, fMCTemp);
}

static int McRun_PickBead_MultiLocal_Equil(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    return Move_MultiLocal_Equil(i, fMCTemp);
}

static int McRun_PickBead_DbPvt(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    return Move_DbPvt(i, fMCTemp);
}

static int McRun_PickChain_Snake(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    i = g_system.bead_info[i][BEAD_CHAINID];
    return Move_Snake(i, fMCTemp);
}

static int McRun_PickChain_Snake_Equil(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    i = g_system.bead_info[i][BEAD_CHAINID];
    return Move_Snake_Equil(i, fMCTemp);
}

static int McRun_PickChain_Trans(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    i = g_system.bead_info[i][BEAD_CHAINID];
    return Move_Trans(i, fMCTemp);
}

static int McRun_PickChain_Trans_Equil(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    i = g_system.bead_info[i][BEAD_CHAINID];
    return Move_Trans_Equil(i, fMCTemp);
}

static int McRun_PickChain_SmallClus_Network(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    i = g_system.bead_info[i][BEAD_CHAINID];
    return Move_SmallClus_Network(i, fMCTemp);
}

static int McRun_PickChain_Pivot(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    i = g_system.bead_info[i][BEAD_CHAINID];
    return Move_Pivot(i, fMCTemp);
}

static int McRun_PickChain_Pivot_Equil(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    i = g_system.bead_info[i][BEAD_CHAINID];
    return Move_Pivot_Equil(i, fMCTemp);
}

static int McRun_PickChain_BranchedRot(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    i = g_system.bead_info[i][BEAD_CHAINID];
    return Move_BranchedRot(i, fMCTemp);
}

static int McRun_PickChain_BranchedRot_Equil(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    i = g_system.bead_info[i][BEAD_CHAINID];
    return Move_BranchedRot_Equil(i, fMCTemp);
}

static int McRun_PickChain_SmallClus_Proximity(float fMCTemp)
{
    int i = rand() % g_system.tot_beads;
    i = g_system.bead_info[i][BEAD_CHAINID];
    return Move_SmallClus_Proximity(i, fMCTemp);
}

/* Old switch: case MV_PR_CLSTR called Move_Clus_Proximity then fell through to
 * default which set nAccept = 0. Preserve that accept accounting. */
static int McRun_Clus_Proximity_ThenReject(float fMCTemp)
{
    (void)Move_Clus_Proximity(fMCTemp);
    return 0;
}

/* Indexed by MV_* / mode from faMCFreq selection. */
static const McMoveRunner s_mcMoves_Prod[MAX_MV] = {
    [MV_NULL] = McRun_Reject,
    [MV_STROT] = McRun_PickBead_Rot,
    [MV_LOCAL] = McRun_PickBead_Local,
    [MV_COLOCAL] = McRun_PickBead_CoLocal,
    [MV_MTLOCAL] = McRun_PickBead_MultiLocal,
    [MV_SNAKE] = McRun_PickChain_Snake,
    [MV_TRANS] = McRun_PickChain_Trans,
    [MV_SMCLSTR] = McRun_PickChain_SmallClus_Network,
    [MV_CLSTR] = Move_Clus_Network,
    [MV_PIVOT] = McRun_PickChain_Pivot,
    [MV_BRROT] = McRun_PickChain_BranchedRot,
    [MV_DBPVT] = McRun_PickBead_DbPvt,
    [MV_PR_SMCLSTR] = McRun_PickChain_SmallClus_Proximity,
    [MV_PR_CLSTR] = McRun_Clus_Proximity_ThenReject,
};

static const McMoveRunner s_mcMoves_Therm[MAX_MV] = {
    [MV_NULL] = McRun_Reject,
    [MV_STROT] = McRun_PickBead_Local_Equil,
    [MV_LOCAL] = McRun_PickBead_Local_Equil,
    [MV_COLOCAL] = McRun_PickBead_Local_Equil,
    [MV_MTLOCAL] = McRun_PickBead_MultiLocal_Equil,
    [MV_SNAKE] = McRun_PickChain_Snake_Equil,
    [MV_TRANS] = McRun_PickChain_Trans_Equil,
    [MV_SMCLSTR] = McRun_PickChain_Trans_Equil,
    [MV_CLSTR] = McRun_PickChain_Trans_Equil,
    [MV_PIVOT] = McRun_PickChain_Pivot_Equil,
    [MV_BRROT] = McRun_PickChain_BranchedRot_Equil,
    [MV_DBPVT] = McRun_PickBead_DbPvt,
    [MV_PR_SMCLSTR] = McRun_PickChain_Trans_Equil,
    [MV_PR_CLSTR] = McRun_Reject,
};

static int MC_PickMode(void)
{
    int mode = 0;
    float prob = (float)rand() / (float)RAND_MAX;
    int i;
    for (i = 0; i < MAX_MV; i++)
    {
        if (prob < g_mc.faMCFreq[i])
        {
            mode = i;
            break;
        }
    }
    return mode;
}

/// MC_Step - given the MC move frequencies, picks which move to try.
/// \param fMCTemp
/// \return bAccept - the acceptance state where MCMove=bAccept/12 and
/// MCAccept=bAccept%2.
int MC_Step(float fMCTemp)
{
    const int mode = MC_PickMode();
    const int nAccept = s_mcMoves_Prod[mode](fMCTemp);
    g_mc.cycle.naMCAccepMat[nAccept][mode]++;
    return mode * 12 + nAccept;
}

/// MC_Step_Therm - given the MC move frequencies, picks non-anisotropic
/// variants of the moves. \param fMCTemp \return bAccept - the acceptance state
/// where MCMove=bAccept/12 and MCAccept=bAccept%2.
int MC_Step_Therm(float fMCTemp)
{
    const int mode = MC_PickMode();
    const int nAccept = s_mcMoves_Therm[mode](fMCTemp);
    g_mc.cycle.naMCAccepMat[nAccept][mode]++;
    return mode * 12 + nAccept;
}
