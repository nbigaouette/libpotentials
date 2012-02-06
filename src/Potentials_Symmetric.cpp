// ***************************************************************
//
// Functions for the potential and field of a two gaussian-shaped
// charge distribution particles.
//
// See "Notes on gaussian particles.pdf" equations (38) and (56).
//
// The "symmetric" potential was calculated for two gaussian-shaped
// charge distribution interacting. The potential in that case is
// the potential _energy_ divided by the total charge of a
// particle. If these two particles have the same shape (same
// sigma), then the "GaussianDistribution" potential gives
// the exact same thing if sigma_symmetric = sqrt(2)/2 * sigma_gaussiandistribution
// See "Notes on potential shapes.pdf" equation (12) and (15)
//
// To make sure both "Symmetric" and "GaussianDistribution" potentials
// gives the same thing, sigma_symmetric is multiplied by "sqrt(2)/2"
// See line ~63 here.
//
// ***************************************************************


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

#include "Potentials_Symmetric.hpp"
#include "Potentials_Coulomb.hpp"


namespace Potentials_Symmetric
{
    fdouble sigma;
}

using namespace libpotentials;
using namespace Potentials_Symmetric;


void Initialize_Symmetric(const fdouble cutoff_base_potential, const fdouble cutoff_radius)
{
    if (cutoff_base_potential > 0.0)
    {
        sigma = one_over_4Pieps0 * e0 / (cutoff_base_potential * sqrt_Pi);
        libpotentials_private::cutoff_radius = eight * sigma;
    }
    else if (cutoff_radius > 0.0)
    {
        // NOTE: "sqrt_2/two" because here, the calculation is done for a single gaussian distribution particle.
        sigma = cutoff_radius / eight * (sqrt_2/two);
        libpotentials_private::cutoff_base_potential = one_over_4Pieps0 * e0 / (sigma * sqrt_Pi);
    }

    assert(sigma > 1.0e-25 * bohr_to_m);
}

// **************************************************************
void Potentials_Set_Parameters_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    potparams.sym_cs1 = Get_Charge_State(p1);
    potparams.sym_cs2 = Get_Charge_State(p2);

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);

    // It does not make sense to calculate effect of particle 2 if
    // its charge is 0
    if (potparams.sym_cs2 != 0)
    {
        potparams.kQ2 = one_over_4Pieps0 * fdouble(potparams.sym_cs2) * e0;
    }
    else
    {
        potparams.kQ2 = 0.0;
    }
}

// **************************************************************
fdouble Calculate_Potential_Cutoff_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams)
/**
 * The electrostatic potential for two charge distribution interaction
 * does not make sense. The electrostatic potential is only calculated
 * to later calculate the electrostatic potential _energy_.
 * Because the codes (treecode, MD) are expecting point charges interactions,
 * they will multiply the electrostatic potential by the charge of
 * the current particle to obtain its electrostatic potential energy.
 * To be exact, the codes would need to integrate the charge density
 * times the potential over the whole space, which is impractical.
 *
 * So the "potential" calculated here is really the potential _energy_,
 * but divided by the charge of the current particle. The potential
 * energy, or the work, is given by the equation (38) of the document
 * "Notes on gaussian particles".
 *
 * The "real" electrostatic potential of equation (14) cannot be used
 * here.
 */
{
    Check_if_LibPotentials_is_initialized();

    fdouble potential = 0.0;

    // It does not make sense to calculate effect of particle 2 if
    // its charge is 0
    if (potparams.sym_cs2 != 0)
    {
        if (potparams.r > libpotentials_private::cutoff_radius)
        {
            // If the distance is not less then the shielding radius, get the
            // normal Coulomb potential.

            potential = Coulomb_Potential(potparams.kQ2, potparams.r);
        }
        else
        {
            // Else, use the lookup table for erf(x)/x
            // Note that y = x / sqrt(2) where x is the unitless distance from GaussianDistribution
            const fdouble y = potparams.r / ( two * sigma );
            Assert_isinf_isnan(y);
            potential = potparams.kQ2 / (two * sigma) * libpotentials_private::lut_potential.read(y);

            //printf("kQ2 = %20.15g  sigma = %20.15g  potential = %20.15g  E = %20.15g   x=%20.15g\n",
            //       potparams.kQ2, sigma*si_to_au_length, potential*si_to_au_pot, 0*si_to_au_field, y);
        }
    }

    return potential;
}

// **************************************************************
void Set_Field_Cutoff_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3])
{
    Check_if_LibPotentials_is_initialized();

    // It does not make sense to calculate effect of particle 2 if
    // its charge is 0
    if (potparams.sym_cs2 != 0)
    {
        if (potparams.r > libpotentials_private::cutoff_radius)
        {
            // If the distance is not less then the shielding radius, get the
            // normal Coulomb potential.

            Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);
        }
        else
        {
            // Else, use the lookup table
            // Get E/r
            // Note that y = x / sqrt(2) where x is the unitless distance from GaussianDistribution
            const fdouble y = potparams.r / (two * sigma);
            const fdouble E_over_r = potparams.kQ2 / (four * sigma * sigma * sigma) * libpotentials_private::lut_field.read(y);

            Assert_isinf_isnan(E_over_r);

            for (int d = 0 ; d < 3 ; d++)
            {
                // We multiply E/r by r to get E
                E[d]  += E_over_r * potparams.dr[d];
            }

            //printf("kQ2 = %20.15g  sigma = %20.15g  potential = %20.15g  E = %20.15g   x=%20.15g\n",
            //       potparams.kQ2, sigma*si_to_au_length, 0, E[0]*si_to_au_field, y);
        }
    }
}


// ********** End of file ***************************************
