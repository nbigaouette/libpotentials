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

#include "Potentials_Simple_CJ.hpp"
#include "Potentials_Coulomb.hpp"


using namespace libpotentials;


// **************************************************************
void Initialize_Simple_CJ(const fdouble cutoff_base_potential, const fdouble cutoff_radius)
{
    if (cutoff_base_potential > 0.0)
    {
        std_cout
            << "ERROR: Cannot initialize potential library using a potential depth for\n"
            <<  "      the \"SimpleCJ\" potential shape!\n"
            << "       Use a cutoff radius instead."<< std::endl;
        std_cout.Flush();
        abort();
    }
    else if (cutoff_radius > 0.0)
    {
        libpotentials_private::cutoff_base_potential = one / (cutoff_radius * libpotentials::m_to_bohr) * libpotentials::Eh_to_eV;
    }
}

// **************************************************************
void Potentials_Set_Parameters_Simple_CJ(
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

    potparams.scj_cs = std::max(std::abs(Get_Charge_State(p1)), std::abs(charge_state2));
    potparams.scj_cr = bohr_to_m * ((m_to_bohr*libpotentials_private::cutoff_radius) + half * (fdouble(potparams.scj_cs)-one) * one_over_seven);
}

// **************************************************************
fdouble Calculate_Potential_Cutoff_Simple_CJ(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    fdouble potential;

    if (potparams.r > potparams.scj_cr)
        potential = Coulomb_Potential(potparams.kQ2, potparams.r);
    else
        potential = Coulomb_Potential(potparams.kQ2, potparams.scj_cr);

    return potential;
}

// **************************************************************
void Set_Field_Cutoff_Simple_CJ(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3])
{
    Check_if_LibPotentials_is_initialized();

    if (potparams.r > potparams.scj_cr)
        Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);
    // Simple cutoff means the potential is constant for distance
    // less then the cutoff radius. If the potential is constant, then
    // their is not field.
}

// ********** End of file ***************************************
