/****************************************************************

****************************************************************/

#include <iostream> // cout
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <iostream> // cout
#include <cstring>  // memset
#include <vector>

#include <StdCout.hpp>
#include <Assert.hpp>

#include "LibPotentials.hpp"
#include "Constants.hpp"
#include "Code_Functions_Declarations.hpp"

#include "Potentials_Coulomb.hpp"


using namespace libpotentials;


// **************************************************************
// ********** Global variable (to this cpp file) ****************
// ********** These are used as global parameters for  **********
// ********** the different potentials shapes.  They are ********
// ********** set/calculated once at initialization, and ********
// ********** treated as constants afterwards. ******************
// **************************************************************


// **************************************************************
fdouble Coulomb_Potential(const fdouble kQ, const fdouble r)
/**
 * Calculate the Coulomb potential due to a point charge.
 * @param   kQ      Charge times Coulomb constant [V.m]
 * @param   r       Distance where the potential is to be evaluated [m]
 * @return  Coulombic potential [V]
 */
{
    if (r > std::numeric_limits<fdouble>::min())
        return kQ / r;
    else
        return std::numeric_limits<fdouble>::max();
}

// **************************************************************
void Set_Coulomb_Field(const fdouble phi12, fdouble E[3],
                       const fdouble dr[3], const fdouble r2)
/**
 * Calculate the Coulomb field due to a point charge.
 * @param   phi12   Input:  Potential at position 1 due to charge 2 [V]
 * @param   E       Output: Electrostatic field at position 1 [V.m^-1]
 * @param   dr      Input:  Vector from position 2 to position 1 [m]
 * @param   r2      Input:  Length squared of vector from position 2 to position 1 [m^2]
 */
{
    if (r2 > std::numeric_limits<fdouble>::min())
    {
        // E = r[:] . (V / |r[:]|^2)
        const fdouble phi12_over_dr2 = phi12 / r2;
        for (int d = 0 ; d < 3 ; d++)
            E[d] += dr[d] * phi12_over_dr2;
    }
    else
    {
        for (int d = 0 ; d < 3 ; d++)
            E[d] = std::numeric_limits<fdouble>::max();
    }
}

// **************************************************************
void Potentials_Set_Parameters_Coulomb(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);

    // We only add field if other particle has a charge not equal 0
    int charge_state2 = Get_Charge_State(p2);

    if (charge_state2 != 0)
    {
        potparams.kQ2 = one_over_4Pieps0 * fdouble(charge_state2) * e0;
    }
    else
    {
        potparams.kQ2 = 0.0;
    }
}

// **************************************************************
fdouble Calculate_Potential_Cutoff_Coulomb(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    return Coulomb_Potential(potparams.kQ2, potparams.r);
}

// **************************************************************
void Set_Field_Cutoff_Coulomb(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3])
{
    Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);
}

// ********** End of file ***************************************
