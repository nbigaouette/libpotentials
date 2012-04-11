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

#include "Potentials_Simple.hpp"
#include "Potentials_Coulomb.hpp"


using namespace libpotentials;


// **************************************************************
void Initialize_Simple(const fdouble cutoff_base_potential, const fdouble cutoff_radius)
{
    if (cutoff_base_potential > 0.0)
    {
        libpotentials_private::cutoff_radius         = one / (cutoff_base_potential * libpotentials::eV_to_Eh)  * libpotentials::bohr_to_m;
    }
    else if (cutoff_radius > 0.0)
    {
        libpotentials_private::cutoff_base_potential = one / (cutoff_radius * libpotentials::m_to_bohr) * libpotentials::Eh_to_eV;
    }
}

// **************************************************************
void Potentials_Set_Parameters_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams)
/**
 * p1 is the particle of interest, p2 is "the other" particle
 * acting on p1.
 *
 */
{
    Check_if_LibPotentials_is_initialized();

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);

    // we only add field if other particle has a charge not equal 0
    fdouble Q2 = Get_Charge(p2);
    int charge_state2 = Get_Charge_State(p2);
    if (charge_state2 != 0)
    {
        potparams.kQ2 = one_over_4Pieps0 * Q2;
    }
    else
    {
        potparams.kQ2 = 0.0;
    }
}

// **************************************************************
fdouble Calculate_Potential_Cutoff_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    fdouble potential;

    if (potparams.r > libpotentials_private::cutoff_radius)
        potential = Coulomb_Potential(potparams.kQ2, potparams.r);
    else
        potential = Coulomb_Potential(potparams.kQ2, libpotentials_private::cutoff_radius);

    return potential;
}

// **************************************************************
void Set_Field_Cutoff_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3])
{
    Check_if_LibPotentials_is_initialized();

    if (potparams.r > libpotentials_private::cutoff_radius)
        Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);
    // Simple cutoff means the potential is constant for distance
    // less then the cutoff radius. If the potential is constant, then
    // their is not field.
}

// ********** End of file ***************************************
