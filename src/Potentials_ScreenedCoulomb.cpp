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

#include "Potentials_ScreenedCoulomb.hpp"
#include "Potentials_Coulomb.hpp"


using namespace libpotentials;

// **************************************************************
void Potentials_Set_Parameters_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    std_cout << "FIXME: sc_alpha is not set in Potentials_Set_Parameters_ScreenedCoulomb()\n";
    std_cout.Flush();
    abort();

    Check_if_LibPotentials_is_initialized();

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);
    potparams.kQ2 = one_over_4Pieps0 * Get_Charge(p2);
}

// **************************************************************
fdouble Calculate_Potential_Cutoff_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    const fdouble one_over_r_plus_alpha = libpotentials::one / (potparams.r + sc_alpha);
    const fdouble phi = potparams.kQ2 * (libpotentials::two * potparams.r + sc_alpha) * one_over_r_plus_alpha * one_over_r_plus_alpha * libpotentials::half;
    return phi;
}

// **************************************************************
void Set_Field_Cutoff_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3])
{
    Check_if_LibPotentials_is_initialized();

    const fdouble Eabs = potparams.kQ2 / std::pow(potparams.r + sc_alpha,3);
    for (int d = 0 ; d < 3 ; d++)
        E[d] += potparams.dr[d] * Eabs;
}


// ********** End of file ***************************************
