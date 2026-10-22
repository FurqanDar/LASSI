#include "mcmove.h"
#include "energy.h"
#include "structure.h"

/// OP_Rotation - given the rotation operation PivotM, rotate beadID using
/// tmpR[POS_MAX] as the origin. Note the global array g_system.naTempR[POS_MAX] stores
/// the locations of the post-rotated beads. Furthermore, if PivotM is 10, the
/// null operation is still performed. Mathematically, we have that
/// f$vec{r}^prime = hat{R}_i(theta)vec{r}f$ where f$vec{r}f$ is the new
/// position, f$hat{R}_i(theta)f$ is the standard rotation matrix where i=X,Y,Z
/// and $theta=45^circ,180^circ,270^circ$. For now, only 90 degree rotations are
/// implemented where the Rot_{#1}_{#2} functions/sub-routines below explicitly
/// calculate what the new vector is given the axis ({#1}) and angle ({#2}).
/// Aribitrary rotations are a little harder on lattices because some rotations
/// do not fully overlap with points on the lattice being used. In our case, it
/// is a simple cubic lattice.
/// \param PivotM - this is the operation to perform
/// \param beadID - the ID of the bead to move
/// \param tmpR   - the location of the anchor point, or, origin around which to
/// rotate
void OP_Rotation(int PivotM, int beadID, int* tmpR)
{
    int j;
    switch (PivotM)
        {
            case 1:
                Rot_X_90(beadID, tmpR);
                break;

            case 2:
                Rot_X_180(beadID, tmpR);
                break;

            case 3:
                Rot_X_270(beadID, tmpR);
                break;

            case 4:
                Rot_Y_90(beadID, tmpR);
                break;

            case 5:
                Rot_Y_180(beadID, tmpR);
                break;

            case 6:
                Rot_Y_270(beadID, tmpR);
                break;

            case 7:
                Rot_Z_90(beadID, tmpR);
                break;

            case 8:
                Rot_Z_180(beadID, tmpR);
                break;

            case 9:
                Rot_Z_270(beadID, tmpR);
                break;

            default:
                // printf("How\n");
                for (j = 0; j < POS_MAX; j++)
                    {
                        g_system.naTempR[j] = g_system.bead_info[beadID][j]; // Does nothing; should make the move fail.
                    }
                break;
        }
}

void Rot_X_90(int beadID, const int tmpR[POS_MAX])
{

    int j;
    // Performs a rotation by 90 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            g_system.naTempR[j] = g_system.bead_info[beadID][j];
        }
    g_system.naTempR[POS_Y] = tmpR[POS_Y] + tmpR[POS_Z] - g_system.bead_info[beadID][POS_Z];
    g_system.naTempR[POS_Z] = g_system.bead_info[beadID][POS_Y] - tmpR[POS_Y] + tmpR[POS_Z];
    // Adjusting for periodic boundaries
    g_system.naTempR[POS_Y] = (g_system.naTempR[POS_Y] + g_system.naBoxSize[POS_Y]) % g_system.naBoxSize[POS_Y];
    g_system.naTempR[POS_Z] = (g_system.naTempR[POS_Z] + g_system.naBoxSize[POS_Z]) % g_system.naBoxSize[POS_Z];
}

void Rot_X_180(int beadID, const int tmpR[POS_MAX])
{

    int j;
    // Performs a rotation by 180 degrees around x-axis. tmpR is the new origin
    // for the rotation.
    for (j = 0; j < POS_MAX; j++)
        {
            g_system.naTempR[j] = g_system.bead_info[beadID][j];
        }
    g_system.naTempR[POS_Y] = 2 * tmpR[POS_Y] - g_system.bead_info[beadID][POS_Y];
    g_system.naTempR[POS_Z] = 2 * tmpR[POS_Z] - g_system.bead_info[beadID][POS_Z];
    // Adjusting for periodic boundaries
    g_system.naTempR[POS_Y] = (g_system.naTempR[POS_Y] + g_system.naBoxSize[POS_Y]) % g_system.naBoxSize[POS_Y];
    g_system.naTempR[POS_Z] = (g_system.naTempR[POS_Z] + g_system.naBoxSize[POS_Z]) % g_system.naBoxSize[POS_Z];
}

void Rot_X_270(int beadID, const int tmpR[POS_MAX])
{

    int j;
    // Performs a rotation by 270 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            g_system.naTempR[j] = g_system.bead_info[beadID][j];
        }
    g_system.naTempR[POS_Y] = g_system.bead_info[beadID][POS_Z] + tmpR[POS_Y] - tmpR[POS_Z];
    g_system.naTempR[POS_Z] = tmpR[POS_Y] + tmpR[POS_Z] - g_system.bead_info[beadID][POS_Y];
    // Adjusting for periodic boundaries
    g_system.naTempR[POS_Y] = (g_system.naTempR[POS_Y] + g_system.naBoxSize[POS_Y]) % g_system.naBoxSize[POS_Y];
    g_system.naTempR[POS_Z] = (g_system.naTempR[POS_Z] + g_system.naBoxSize[POS_Z]) % g_system.naBoxSize[POS_Z];
}

void Rot_Y_90(int beadID, const int tmpR[POS_MAX])
{
    int j;
    // Performs a rotation by 90 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            g_system.naTempR[j] = g_system.bead_info[beadID][j];
        }
    g_system.naTempR[POS_X] = (g_system.bead_info[beadID][POS_Z] - tmpR[POS_Z] + tmpR[POS_X]);
    g_system.naTempR[POS_Z] = (tmpR[POS_X] + tmpR[POS_Z] - g_system.bead_info[beadID][POS_X]);
    // Adjusting for periodic boundaries
    g_system.naTempR[POS_X] = (g_system.naTempR[POS_X] + g_system.naBoxSize[POS_X]) % g_system.naBoxSize[POS_X];
    g_system.naTempR[POS_Z] = (g_system.naTempR[POS_Z] + g_system.naBoxSize[POS_Z]) % g_system.naBoxSize[POS_Z];
}

void Rot_Y_180(int beadID, const int tmpR[POS_MAX])
{
    int j;
    // Performs a rotation by 270 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            g_system.naTempR[j] = g_system.bead_info[beadID][j];
        }
    g_system.naTempR[POS_X] = (2 * tmpR[POS_X] - g_system.bead_info[beadID][POS_X]);
    g_system.naTempR[POS_Z] = (2 * tmpR[POS_Z] - g_system.bead_info[beadID][POS_Z]);
    // Adjusting for periodic boundaries
    g_system.naTempR[POS_X] = (g_system.naTempR[POS_X] + g_system.naBoxSize[POS_X]) % g_system.naBoxSize[POS_X];
    g_system.naTempR[POS_Z] = (g_system.naTempR[POS_Z] + g_system.naBoxSize[POS_Z]) % g_system.naBoxSize[POS_Z];
}

void Rot_Y_270(int beadID, const int tmpR[POS_MAX])
{
    int j;
    // Performs a rotation by 270 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            g_system.naTempR[j] = g_system.bead_info[beadID][j];
        }
    g_system.naTempR[POS_X] = (tmpR[POS_X] + tmpR[POS_Z] - g_system.bead_info[beadID][POS_Z]);
    g_system.naTempR[POS_Z] = (g_system.bead_info[beadID][POS_X] - tmpR[POS_X] + tmpR[POS_Z]);
    // Adjusting for periodic boundaries
    g_system.naTempR[POS_X] = (g_system.naTempR[POS_X] + g_system.naBoxSize[POS_X]) % g_system.naBoxSize[POS_X];
    g_system.naTempR[POS_Z] = (g_system.naTempR[POS_Z] + g_system.naBoxSize[POS_Z]) % g_system.naBoxSize[POS_Z];
}

void Rot_Z_90(int beadID, const int tmpR[POS_MAX])
{
    int j;
    // Performs a rotation by 270 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            g_system.naTempR[j] = g_system.bead_info[beadID][j];
        }
    g_system.naTempR[POS_X] = (tmpR[POS_X] + tmpR[POS_Y] - g_system.bead_info[beadID][POS_Y]);
    g_system.naTempR[POS_Y] = (g_system.bead_info[beadID][POS_X] - tmpR[POS_X] + tmpR[POS_Y]);
    // Adjusting for periodic boundaries
    g_system.naTempR[POS_X] = (g_system.bead_info[beadID][POS_X] + g_system.naBoxSize[POS_X]) % g_system.naBoxSize[POS_X];
    g_system.naTempR[POS_Y] = (g_system.bead_info[beadID][POS_Y] + g_system.naBoxSize[POS_Y]) % g_system.naBoxSize[POS_Y];
}

void Rot_Z_180(int beadID, const int tmpR[POS_MAX])
{
    int j;
    // Performs a rotation by 270 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            g_system.naTempR[j] = g_system.bead_info[beadID][j];
        }
    g_system.naTempR[POS_X] = (2 * tmpR[POS_X] - g_system.bead_info[beadID][POS_X]);
    g_system.naTempR[POS_Y] = (2 * tmpR[POS_Y] - g_system.bead_info[beadID][POS_Y]);
    // Adjusting for periodic boundaries
    g_system.naTempR[POS_X] = (g_system.bead_info[beadID][POS_X] + g_system.naBoxSize[POS_X]) % g_system.naBoxSize[POS_X];
    g_system.naTempR[POS_Y] = (g_system.bead_info[beadID][POS_Y] + g_system.naBoxSize[POS_Y]) % g_system.naBoxSize[POS_Y];
}

void Rot_Z_270(int beadID, const int tmpR[POS_MAX])
{
    int j;
    // Performs a rotation by 270 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            g_system.naTempR[j] = g_system.bead_info[beadID][j];
        }
    g_system.naTempR[POS_X] = (g_system.bead_info[beadID][POS_Y] + tmpR[POS_X] - tmpR[POS_Y]);
    g_system.naTempR[POS_Y] = (tmpR[POS_X] + tmpR[POS_Y] - g_system.bead_info[beadID][POS_X]);
    // Adjusting for periodic boundaries
    g_system.naTempR[POS_X] = (g_system.bead_info[beadID][POS_X] + g_system.naBoxSize[POS_X]) % g_system.naBoxSize[POS_X];
    g_system.naTempR[POS_Y] = (g_system.bead_info[beadID][POS_Y] + g_system.naBoxSize[POS_Y]) % g_system.naBoxSize[POS_Y];
}

/// OP_ShuffleRotIndecies - randomly shuffles the incdecies to sample the
/// rotational states around a bead. Ensures that we fully, but randomly, sample
/// the states. The implementation was taken from the provided URL.
void OP_ShuffleRotIndecies(void)
{
    // Algorithm taken from
    // https://www.w3resource.com/c-programming-exercises/array/c-array-exercise-77.php
    // Takes the array and shuffles the numbers in it. Used to randomly, but
    // fully, sample the rotational states.

    int i, j;
    int i_val, j_val;

    for (i = MAX_ROTSTATES - 2; i > 0; i--)
        {
            j                   = rand() % (i + 1);
            i_val               = g_system.naRot_IndArr[i];
            j_val               = g_system.naRot_IndArr[j];
            g_system.naRot_IndArr[i] = j_val;
            g_system.naRot_IndArr[j] = i_val;
        }
}

/// OP_ShuffleArray: Given an array of at-most size arr_size, we shuffle the
/// elements in dum_arr. Note that if arr_size is smaller than the actual size
/// of the array, only the elements less than arr_size shall be shuffled. So
/// this _could_ be used to shuffle a starting subset of the array.
/// \param arr_size
/// \param dum_arr
void OP_ShuffleArray(const int arr_size, int* dum_arr)
{
    // Algorithm taken from
    // https://www.w3resource.com/c-programming-exercises/array/c-array-exercise-77.php
    // Takes the array and shuffles the numbers in it. Used to randomly, but
    // fully, sample arrays.
    int i;
    int jVal, j;
    for (i = arr_size - 1; i > 0; i--)
        {
            j          = rand() % (i + 1);
            jVal       = dum_arr[i];
            dum_arr[i] = dum_arr[j];
            dum_arr[j] = jVal;
        }
}

/// Check_RotStates_wNeighList: Given this bead and the supplied neighbor-list
/// and numbers, we see which of the neighbors are viable anisotropic partners.
/// For each viable bead, we also record the Rosenbluth factors, and store which
/// bead it is.
/// \param beadID
/// \param resi
/// \param neighList
/// \param neighNum
/// \return The number of beads that beadID could bond to.
int Check_RotStates_wNeighList(int const beadID, int const resi, const int* neighList, int const neighNum)
{
    int CandNums = 0;
    int tmpBead;
    int resj, i;

    for (i = 0; i < neighNum; i++)
        {
            tmpBead = neighList[i];
            resj    = g_system.bead_info[tmpBead][BEAD_TYPE];

            if (g_energy.faEnergy[resi][resj][E_SC_SC] != 0.f)
                {
                    if (g_system.bead_info[tmpBead][BEAD_FACE] == -1 || g_system.bead_info[tmpBead][BEAD_FACE] == beadID)
                        {
                            g_energy.ldaBoltzFac[CandNums]   = g_energy.ldaBoltzFacNorm[resi][resj];
                            g_energy.naRotTrial[0][CandNums] = tmpBead;
                            CandNums++;
                        }
                }
        }

    return CandNums;
}

int Check_RotStates_wNeighList_ForChains(int const beadID, int const resi, const int* neighList, int const neighNum)
{
    int CandNums = 0;
    int tmpBead, tmpfaceID, firstbeadID, lastbeadID, chainID;
    int resj, i;

    chainID     = g_system.bead_info[beadID][BEAD_CHAINID];
    firstbeadID = g_system.chain_info[chainID][CHAIN_START]; // First bead of the chain where beadID is.
    lastbeadID  = firstbeadID + g_system.chain_info[chainID][CHAIN_LENGTH];

    for (i = 0; i < neighNum; i++)
        {
            tmpBead = neighList[i];
            resj    = g_system.bead_info[tmpBead][BEAD_TYPE];
            // Check if the tmpBead is a preceding part of the chain.
            // In this case, the tmpBead already decided that it did not bond to any bead.
            // So, it can not be a candidate.
            if (tmpBead >= firstbeadID && tmpBead < beadID) 
                {
                    continue;
                }
            if (g_energy.faEnergy[resi][resj][E_SC_SC] == 0.f)
                {
                    continue;
                }
            
            tmpfaceID = g_system.bead_info[tmpBead][BEAD_FACE];
            // Check if the tmpfaceBead is a remained part of the chain.
            // if the tmpfaceBead is a later checked bead, then tmpBead can be a candidate
            // because the bond between tmpBead and tmpfaceBead was not formed 
            // before deciding the beadID's bond.
            if (tmpfaceID == -1 || ( (tmpfaceID >= beadID) && (tmpfaceID < lastbeadID) ) )
                {
                    g_energy.ldaBoltzFac[CandNums]   = g_energy.ldaBoltzFacNorm[resi][resj];
                    g_energy.naRotTrial[0][CandNums] = tmpBead;
                    CandNums++;
                }
        }

    return CandNums;
}

int Check_RotStates_wNeighList_ForLists(int const beadIdx, const int listSize, const int beadList[MAX_BONDS + 1], int const resi, const int* neighList, int const neighNum)
{
    int CandNums = 0;
    int tmpBead, tmpfaceID;
    int resj, i, j;
    int beadID = beadList[beadIdx];

    for (i = 0; i < neighNum; i++)
        {
            tmpBead = neighList[i];
            resj    = g_system.bead_info[tmpBead][BEAD_TYPE];
            // Check if the bead is already in the list.
            // If it is, then it cannot be a candidate because it
            // already decided that it did not bond to any bead.
            for (j = 0; j < beadIdx; j++) 
                {
                    if (tmpBead == beadList[j])
                        {
                            break;
                        }
                }
            if (j < beadIdx)
                {
                    continue;
                }

            if (g_energy.faEnergy[resi][resj][E_SC_SC] == 0.f)
                {
                    continue; // No interaction between these two types.
                }
           
            tmpfaceID = g_system.bead_info[tmpBead][BEAD_FACE];

            // check if the tmpfaceBead is a remained part of the list.
            for(j= beadIdx; j < listSize; j++) 
                {
                    if (tmpfaceID == beadList[j])
                        {
                            break;
                        }
                }
            // j < listSize means that the tmpfaceID is a remained part of the list.
            // if the tmpfaceID is a later checked bead, then tmpBead can be a candidate
            if (tmpfaceID == -1 || j < listSize )  
                {
                    g_energy.ldaBoltzFac[CandNums]   = g_energy.ldaBoltzFacNorm[resi][resj];
                    g_energy.naRotTrial[0][CandNums] = tmpBead;
                    CandNums++;
                }
        }

    return CandNums;
}

/// @brief this function is used to scan the roational states of the beadID 
///        in the sequencial order of the given range.
/// @param beadID 
/// @param resi 
/// @param smallestBead means the bead which has the smallest beadID in the Range
/// @param largestBead means the bead which has the largest beadID in the Range
/// @param neighList 
/// @param neighNum 
/// @return 
int Check_RotState_wNeighList_ForRange(int const beadID, int const resi, int const smallestBead, int const largestBead, const int* neighList, int const neighNum)
{
    int CandNums = 0;
    int tmpBead, tmpfaceID;
    int resj, i;

    for (i = 0; i < neighNum; i++)
        {
            tmpBead = neighList[i];
            resj    = g_system.bead_info[tmpBead][BEAD_TYPE];
            // Check if the tmpBead is a preceding part of the range.
            // In this case, the tmpBead already decided that it did not bond to any bead.
            // So, it can not be a candidate.
            if (tmpBead >= smallestBead && tmpBead < beadID) 
                {
                    continue;
                }
            if (g_energy.faEnergy[resi][resj][E_SC_SC] == 0.f)
                {
                    continue;
                }
            
            tmpfaceID = g_system.bead_info[tmpBead][BEAD_FACE];
            // Check if the tmpfaceBead is a remained part of the range.
            // if the tmpfaceBead is a later checked bead, then tmpBead can be a candidate
            // because the bond between tmpBead and tmpfaceBead was not formed
            // before deciding the beadID's bond.
            if (tmpfaceID == -1 || ( (tmpfaceID >= beadID) && (tmpfaceID <= largestBead) ) )
                {
                    g_energy.ldaBoltzFac[CandNums]   = g_energy.ldaBoltzFacNorm[resi][resj];
                    g_energy.naRotTrial[0][CandNums] = tmpBead;
                    CandNums++;
                }
        }

    return CandNums;
}


/// Check_RotStatesOld -- goes around beadID, of type resi to built the
/// boltzmann distribution if the rotational states. Used in the old location.
/// Note that this implementation is one where the solvent interaction does not
/// exist for the rotational state. The commented version below is one where the
/// solvent also has a contribution to the boltzmann distribution of the
/// rotational states.
/// \param beadID
/// \param resi
/// \param MyTemp
/// \return The total number of possible candidates
int Check_RotStatesOld(int const beadID, int const resi, float const MyTemp)
{

    int i, j, k, tmpBead;
    int CandNums = 0;
    int r_pos0[POS_MAX], r_search[POS_MAX];
    for (j = 0; j < POS_MAX; j++)
        {
            r_pos0[j] = g_system.bead_info[beadID][j];
        }
    for (k = 0; k < MAX_ROTSTATES - 1; k++)
        {
            i = g_system.naRot_IndArr[k];

            LatPos_add_wPBC(r_search, r_pos0, g_system.naLocalArr[i]);

            tmpBead = Lat_Ind_FromVec(r_search);
            tmpBead = g_system.naTotLattice[tmpBead];
            if (tmpBead != -1)
                {
                    j = g_system.bead_info[tmpBead][BEAD_TYPE];
                    if (g_energy.faEnergy[resi][j][E_SC_SC] != 0 &&
                        (g_system.bead_info[tmpBead][BEAD_FACE] == -1 || g_system.bead_info[tmpBead][BEAD_FACE] == beadID))
                        {
                            g_energy.ldaBoltzFac[CandNums]   = g_energy.ldaBoltzFacNorm[resi][j];
                            g_energy.naRotTrial[0][CandNums] = tmpBead;
                            CandNums++;
                        }
                }
        }

    return CandNums;
}

/// Check_RotStatesNew -- goes around beadID, of type resi to built the
/// boltzmann distribution if the rotational states. Used in the new proposed
/// location. Note that this implementation is one where the solvent interaction
/// does not exist for the rotaitonal state. The commented version below is one
/// where the solvent also has a contribution to the boltzmann distribution of
/// the rotational states.
/// \param beadID
/// \param resi
/// \param MyTemp
/// \return The total number of possible candidates
int Check_RotStatesNew(int const beadID, int const resi, float const MyTemp)
{

    int i, j, k, tmpBead;
    int CandNums = 0;
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    for (j = 0; j < POS_MAX; j++)
        {
            tmpR[j] = g_system.bead_info[beadID][j];
        }
    for (k = 0; k < MAX_ROTSTATES - 1; k++)
        {
            i = g_system.naRot_IndArr[k];
            for (j = 0; j < POS_MAX; j++)
                {
                    tmpR2[j] = (tmpR[j] + g_system.naLocalArr[i][j] + g_system.naBoxSize[j]) % g_system.naBoxSize[j];
                }
            tmpBead = Lat_Ind_FromVec(tmpR2);
            tmpBead = g_system.naTotLattice[tmpBead];
            if (tmpBead != -1)
                {
                    j = g_system.bead_info[tmpBead][BEAD_TYPE];
                    if (g_energy.faEnergy[resi][j][E_SC_SC] != 0 &&
                        (g_system.bead_info[tmpBead][BEAD_FACE] == -1 || g_system.bead_info[tmpBead][BEAD_FACE] == beadID))
                        {
                            g_energy.ldaBoltzFac[CandNums]   = g_energy.ldaBoltzFacNorm[resi][j];
                            g_energy.naRotTrial[0][CandNums] = tmpBead;
                            CandNums++;
                        }
                }
        }

    return CandNums;
}

/// OP_NormalizeRotState -- normalizes the rotational boltzmann distribution,
/// and generates the cumulative distribution so that we can sample from it.
/// \param beadVal - g_energy.ldaBoltzNorm[MAX_VALENCY] stores the total Rosenbluth weights
/// for moves that sample multiple stickers. Note that g_energy.ldaBoltzFac[] is now the
/// cumulative boltzmann distribution, which is sampled to propose bonds.
/// \param CandNums - total possible bonding partners for this particular sticker case.
/// A dissociated state integrated into the cumulative distribution implicitly.
/// The normalized Boltzmann factor for the dissociated state is considered to be
/// appended to the end of the g_energy.ldaBoltzFac[] array. This is because the
/// dissociated state is always a possible state.
void OP_NormalizeRotState(const int beadVal, const int CandNum)
{
    int i;

    g_energy.ldaBoltzNorm[beadVal] = 1.; // Default value is a dissociated state

    for (i = 0; i < CandNum; i++)
        {
            g_energy.ldaBoltzNorm[beadVal] += g_energy.ldaBoltzFac[i]; // Rosenbluth weights for beadVal
        }
    for (i = 0; i < CandNum; i++)
        {
            g_energy.ldaBoltzFac[i] /= g_energy.ldaBoltzNorm[beadVal];
        }
    for (i = 1; i < CandNum; i++)
        {
            g_energy.ldaBoltzFac[i] += g_energy.ldaBoltzFac[i - 1];
        }
}


/// \param CandNum - the number of possible bonding candidates.
/// \return The new proposed bonding partner.
int OP_PickRotState(int CandNum)
{
    int newRot = -1;
    int i;

    lLDub fProb;

    fProb = (lLDub) rand() / (lLDub) RAND_MAX;
    for (i = 0; i < CandNum; i++)
        {
            if (fProb < g_energy.ldaBoltzFac[i])
                {
                    break;
                }
        }
    // if i == CandNum, it means that all the states are unchosen and the
    // last state, which is the dissociated state, is selected.
    if (i < CandNum) 
    {
        newRot = g_energy.naRotTrial[0][i];
    }

    return newRot;

}
