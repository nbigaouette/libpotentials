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

#include "Potentials_SuperGaussian.hpp"
#include "Potentials_Coulomb.hpp"


using namespace libpotentials;

fdouble sg_m;
fdouble sg_two_m;
fdouble sg_one_over_two_m;
fdouble sg_m_pow_one_over_two_m;
fdouble sg_exp_one_over_two_m;

// **************************************************************
void Initialize_SuperGaussian(const int &m,
                              const fdouble cutoff_base_potential,
                              const fdouble cutoff_radius)
/**
 * Initialize super-gaussian like potential
 * @param   m       Order of the Super-Gaussian (m=1 is a simple gaussian) [-]
 */
{
    sg_m = fdouble(m);
    sg_two_m = fdouble(2 * m);

    if (sg_m >= 1)
    {
        sg_one_over_two_m         = libpotentials::one / fdouble(sg_two_m);
        sg_m_pow_one_over_two_m   = std::pow(sg_m, sg_one_over_two_m);
        sg_exp_one_over_two_m     = std::exp(sg_one_over_two_m);
    }
    else
    {
        sg_one_over_two_m         = std::numeric_limits<fdouble>::max();
        sg_m_pow_one_over_two_m   = 0.0;
        sg_exp_one_over_two_m     = 0.0;
    }

    if (cutoff_base_potential > 0.0)
    {
        libpotentials_private::cutoff_radius = one_over_4Pieps0 * e0 / cutoff_base_potential * sg_exp_one_over_two_m;
    }
    else if (cutoff_radius > 0.0)
    {
        libpotentials_private::cutoff_base_potential = one_over_4Pieps0 * e0 / cutoff_radius * sg_exp_one_over_two_m;
    }
}

// **************************************************************
void Potentials_Set_Parameters_SuperGaussian(
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

    if( charge_state2 != 0 )
    {
        // B is the electrostatic potential right on top of particle 2
        // Its units are Volts (since its a potential, not to be
        // confused with potential _energy_). The Ip is the
        // potential energy we are choosing for a point charge of "e0"
        // if located at r=0 (on to of particle 2). So we first transform
        // the Ip, being in eV, from eV to J, and divide
        // by "e0" (the charge of the point charge) to get the potential
        // from the potential energy. This result in multipling
        // by 1.0, which is not done.
        if (Q2 > std::numeric_limits<fdouble>::min())
        {   // Positive charge (ion), positive potential
            potparams.B = libpotentials_private::cutoff_base_potential*fdouble(charge_state2);
            //potparams.B *= eV_to_J * (1.0 / e0); // *= 1.0
        }
        else
        {   // Negative charge (electron), negative potential
            // B is the (negative of the)
            // ionization potential of the ion.
            potparams.B = -libpotentials_private::cutoff_base_potential;
            //potparams.B *= eV_to_J * (1.0 / e0); // *= 1.0
        }

        potparams.kQ2 = one_over_4Pieps0 * Q2;
        potparams.kQ2_over_B = potparams.kQ2 / potparams.B;

        // Width of Super-Gaussian
        if (potparams.r <= libpotentials_private::cutoff_radius)
        {
            potparams.sg_sigma = potparams.kQ2_over_B
                                    * sg_m_pow_one_over_two_m
                                    * sg_exp_one_over_two_m;
            potparams.sg_r_over_sigma_two_m = std::pow(potparams.r/potparams.sg_sigma, sg_two_m);
            potparams.sg_exp_half_r_over_sigma_two_m =
                    std::exp( -libpotentials::half * potparams.sg_r_over_sigma_two_m );
        }
    }
    else
    {
        potparams.kQ2           = 0.0;
        potparams.kQ2_over_B    = 0.0;
        potparams.B             = 0.0;
    }
}

// **************************************************************
fdouble Calculate_Potential_Cutoff_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    fdouble phi12;   // Electrostatic potential

    if (potparams.r <= libpotentials_private::cutoff_radius)
    {
        // If the distance between two bodys is less than the shielding
        // radius, we use the special Super-Gaussian potential instead
        // of the Coulomb potential.
        phi12 = potparams.B * potparams.sg_exp_half_r_over_sigma_two_m;
    }
    else
    {
        // If the distance is not less then the shielding radius, get the
        // normal Coulomb potential.
        phi12 = Coulomb_Potential(potparams.kQ2, potparams.r);
    }
    return phi12;
}

// **************************************************************
void Set_Field_Cutoff_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3])
{
    Check_if_LibPotentials_is_initialized();

    if (potparams.r > libpotentials_private::cutoff_radius)
        Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);
    else
    {
        // The derivative of a super gaussian w.r. to r is:
        fdouble diff_sg = -(potparams.B*sg_m*potparams.one_over_r)
                            * potparams.sg_r_over_sigma_two_m
                            * potparams.sg_exp_half_r_over_sigma_two_m;

        // We have the norm of the gradient of the potential (diff_sg),
        // we need to multiply this by the unit vector, to get
        // the electric field.
        fdouble unit_dr[3];
//         MULVS(unit_dr, dr, one_over_distance);  // Calculate unit vector
//         ADDMULVS(E, unit_dr, -diff_sg);         // Add to the electric field
//                                                 // the gradient of the potential.
        for (int d = 0 ; d < 3 ; d++)
        {
            unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
            E[d]  += unit_dr[d] * (-diff_sg);
        }
    }
}

// ********** End of file ***************************************
