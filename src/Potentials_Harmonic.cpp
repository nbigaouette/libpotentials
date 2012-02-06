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

#include "Potentials_Harmonic.hpp"
#include "Potentials_Coulomb.hpp"


using namespace libpotentials;


// **************************************************************
void Initialize_Harmonic(const fdouble cutoff_base_potential, const fdouble cutoff_radius)
{
    if (cutoff_base_potential > 0.0)
    {
        libpotentials_private::cutoff_radius = three_over_two * libpotentials::one_over_4Pieps0 * libpotentials::e0 / cutoff_base_potential;
    }
    else if (cutoff_radius > 0.0)
    {
        libpotentials_private::cutoff_base_potential = three_over_two * libpotentials::one_over_4Pieps0 * libpotentials::e0 / libpotentials_private::cutoff_radius;
    }
}

// **************************************************************
void Potentials_Set_Parameters_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams)
/**
 * p1 is the particle of interest, p2 is "the other" particle
 * acting on p1.
 *
 */
{
    std_cout << "FIXME: Potentials_Set_Parameters_Harmonic() is probably broken.\n";
    std_cout.Flush();
    abort();

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
            potparams.B = -libpotentials_private::cutoff_base_potential;
            //potparams.B *= (1.0 / e0) * eV_to_J; // *= 1.0
        }

        potparams.kQ2 = one_over_4Pieps0 * Q2;
        potparams.kQ2_over_B = potparams.kQ2 / potparams.B;
    }
    else
    {
        potparams.kQ2           = 0.0;
        potparams.kQ2_over_B    = 0.0;
        potparams.B             = 0.0;
    }
}

// **************************************************************
fdouble Calculate_Potential_Cutoff_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams)
/**
 * Calculate the electrostatic potential energy "phi1" of particle "p1"
 * due to another particle/cell "p2" at distance sqrt(dr2)
 */
{
    Check_if_LibPotentials_is_initialized();

    fdouble phi12;   // Electrostatic potential

    if (potparams.r <= libpotentials_private::cutoff_radius)
    {
        // const fdouble A  = ( 4.0 * B*B*B) / (27.0 * kQ*kQ);
        potparams.h_A  = (four_over_twenty_seven * potparams.B)
                            / (potparams.kQ2_over_B*potparams.kQ2_over_B);
        potparams.h_A_r = potparams.h_A*potparams.r;

        // Get partial potential
        phi12 = -potparams.h_A_r*potparams.r + potparams.B;
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
void Set_Field_Cutoff_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3])
{
    Check_if_LibPotentials_is_initialized();

    if (potparams.r > libpotentials_private::cutoff_radius)
        Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);
    else
    {
        // The derivative of a -Ar^2+B w.r. to r is:
        const fdouble diff_x2 = -libpotentials::two*potparams.h_A*potparams.r;

        // We have the norm of the gradient of the potential (diff_x2),
        // we need to multiply this by the unit vector, to get
        // the electric field.
        fdouble unit_dr[3];
        //MULVS(unit_dr, dr, one_over_distance); // Calculate unit vector
        //ADDMULVS(E, unit_dr, -diff_x2);     // Add to the electric field
        //                                    // the gradient of the potential.
        for (int d = 0 ; d < 3 ; d++)
        {
            unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
            E[d]  += unit_dr[d] * -diff_x2;
        }
    }
}

// ********** End of file ***************************************
