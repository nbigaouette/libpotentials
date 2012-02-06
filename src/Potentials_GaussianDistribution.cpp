// ***************************************************************
//
// Functions for the potential and field of a SINGLE gaussian-shaped
// charge distribution particle.
//
// See "Notes on potential shapes.pdf" equation (12) and (15)
//
// Note that the implementations are overcomplicated.
// Use instead the "Symmetric" potential which gives the EXACT
// same thing but is a lot more simple.
//
// Actually, the two are not exactly the same. The "symmetric"
// potential was calculated for two gaussian-shaped charge distribution
// interacting. The potential in that case is the potential _energy_
// divided by the total charge of a particle. If these two particles
// have the same shape (same sigma), then the "Symmetric" case gives
// the exact same thing if sigma_symmetric = sqrt(2)/2 * sigma
// See "Notes on gaussian particles.pdf" equations (38) and (56).
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


#include "Potentials_GaussianDistribution.hpp"
#include "Potentials_Coulomb.hpp"


using namespace libpotentials;


void Initialize_GaussianDistribution(const fdouble cutoff_base_potential, const fdouble cutoff_radius)
{
    if (cutoff_base_potential > 0.0)
    {
        Potentials_Symmetric::sigma                 = one_over_4Pieps0 * e0 / cutoff_base_potential * sqrt_2_over_pi;
        libpotentials_private::cutoff_radius        = eight * Potentials_Symmetric::sigma;
    }
    else if (cutoff_radius > 0.0)
    {
        Potentials_Symmetric::sigma                 = cutoff_radius / eight;
        libpotentials_private::cutoff_base_potential= one_over_4Pieps0 * e0 / Potentials_Symmetric::sigma * sqrt_2_over_pi;
    }

}

// **************************************************************
void Potentials_Set_Parameters_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams)
/**
 * p1 is the particle of interest, p2 is "the other" particle
 * acting on p1.
 */
{
    Check_if_LibPotentials_is_initialized();

// std::cout << "chg1 = " << Get_Charge(p1) <<" chg2 = " << Get_Charge(p2) << std::endl;
// std::cout << "chgs1 = " << Get_Charge_State(p1) <<" chgs2 = " << Get_Charge_State(p2) << std::endl;

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);
// #ifdef YDEBUG
//     if (potparams.r <= 1.0e-200 || isnan(potparams.r))
//     {
//         Print_Particles(p1, 1);
//         Print_Particles(p2, 1);
//         assert(potparams.r > 1.0e-200);
//         abort();
//     }
// #endif // #ifdef YDEBUG
//     assert(potparams.r > 1.0e-200);

    // we only add field if other particle has a charge not equal 0
    fdouble Q = Get_Charge(p2);
    int distribution1_charge_state = 1,distribution2_charge_state = 1;
    int charge_state2 = Get_Charge_State(p2);
    int charge_state1 = Get_Charge_State(p1);
    int factor = charge_state2;
    //if the first particle is neutral then skip the rest otherwise when the second is an electron
    //one would get negative width for the electron!
    if (charge_state1 == 0 and charge_state2 < 0) {
        potparams.kQ2           = 0.0;
        potparams.kQ2_over_B    = 0.0;
        potparams.B             = 0.0;
        potparams.gd_sigma      = 0.0;
    }
    else if( charge_state2 != 0 )
    {
      //fdouble Ip =element.IpsLowest[std::abs(charge_state2)];
      fdouble Ip = libpotentials_private::cutoff_base_potential;
      //charges are not equal
      //take the higher charge as the distribution
      //Thus set the parameters for the guassian
      //using the values of the larger charge
      //NOTE: it will always do this for p2=electron
      // unless p1 is also an electron
        if ((charge_state2 < charge_state1) && (charge_state1 != 0)){
          Ip = libpotentials_private::cutoff_base_potential;
          if (charge_state2 < 0)
            Q = -Get_Charge(p1);
          else
            Q = Get_Charge(p1);
          //save the charge states so kQ2 = KQ1 * charge_state2/charge_state1
          //which gives the correct scale. Needed for ion-ion
          distribution1_charge_state = charge_state1;
          distribution2_charge_state = std::abs(charge_state2);
          factor = charge_state1;
        }
        //Make the Ip linearly deeper for ions.
        //not used for e- e-
        Ip*=fdouble(factor);

        // B is the electrostatic potential right on top of particle 2
        // Its units are Volts (since its a potential, not to be
        // confused with potential _energy_). The Ip is the
        // potential energy we are choosing for a point charge of "e0"
        // if located at r=0 (on to of particle 2). So we first transform
        // the Ip, being in eV, from eV to J, and divide
        // by "e0" (the charge of the point charge) to get the potential
        // from the potential energy. This result in multipling
        // by 1.0, which is not done.
        if (charge_state2 > 0)
        {   // Positive charge (ion), positive potential
            // This is Max( Ip(p1),Ip(p2) )
            //potparams.B = element.IpsLowest[charge_state2];
            potparams.B = Ip;
            //potparams.B *= eV_to_J * (1.0 / e0); // *= 1.0
        }
        else  // This means the p2 is an electron
        {     // and the ion is 1+ or p1 is an electron

          // If p1 is an electron use Ip[0]
          if (charge_state1 < 0 )
            potparams.B = -libpotentials_private::cutoff_base_potential;
          else //else p1 is an ion
            potparams.B = -Ip;
        }

        potparams.kQ2 = one_over_4Pieps0 * Q;
        potparams.kQ2_over_B = potparams.kQ2 / potparams.B;

        // Make sure the other particle's charge is normalized by the current
        // particle's charge state so forces are symmetric between the two.
        potparams.kQ2 /= fdouble(distribution1_charge_state);
        potparams.kQ2 *= fdouble(distribution2_charge_state);
        Assert_isinf_isnan(potparams.kQ2);

        potparams.gd_sigma = potparams.kQ2_over_B * sqrt_2_over_pi;
/*        std_cout << "Id(p1)="<<Get_Id(p1)<<"  B="<<potparams.B<<"  Cs(p1)="<<Get_Charge_State(p1)<<" Cs(p2)="<<Get_Charge_State(p2)<<" kQ2_over_B="<<potparams.kQ2_over_B<<" well="<<libpotentials_private::cutoff_base_potential<<"\n";*/
//         potparams.cutoff_radius = 4.0 * potparams.gd_sigma;
//         potparams.cutoff_radius = 2.5 * potparams.gd_sigma;
        //potparams.cutoff_radius = two * potparams.gd_sigma;

    }
    else
    {
        potparams.kQ2           = 0.0;
        potparams.kQ2_over_B    = 0.0;
        potparams.B             = 0.0;
        potparams.gd_sigma      = 0.0;
    }
}

// **************************************************************
fdouble Calculate_Potential_Cutoff_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    fdouble phi12=0.0;   // Electrostatic potential

    if      (potparams.r > libpotentials_private::cutoff_radius)
    {
        // If the distance is not less then the shielding radius, get the
        // normal Coulomb potential.

        phi12 = Coulomb_Potential(potparams.kQ2, potparams.r);
    }
    else
    {
        // Else, use the lookup table for erf(x)/x
        phi12 = potparams.kQ2 / (potparams.gd_sigma * sqrt_2) *
                    libpotentials_private::lut_potential.read(potparams.r / (potparams.gd_sigma * sqrt_2));
    }

    return phi12;
}

// **************************************************************
void Set_Field_Cutoff_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3])
{
    Check_if_LibPotentials_is_initialized();

    if      (potparams.r > libpotentials_private::cutoff_radius)
    {
        // If the distance is not less then the shielding radius, get the
        // normal Coulomb potential.

        Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);

        //std_cout << "High range field:    Expansion: dr = (" << m_to_bohr*potparams.dr[0] << ", " << m_to_bohr*potparams.dr[1] << ", " << m_to_bohr*potparams.dr[2] << ")   E = (" << E[0] <<", "<< E[1] <<", "<< E[2] << ")\n";
    }
    else
    {
        // Else, use the lookup table
        // Get E/r
        const fdouble E_over_r = potparams.kQ2 / (sqrt_2 * potparams.gd_sigma * potparams.gd_sigma * potparams.gd_sigma) *
                                        libpotentials_private::lut_field.read(potparams.r / (potparams.gd_sigma * sqrt_2));

        Assert_isinf_isnan(E_over_r);

        for (int d = 0 ; d < 3 ; d++)
        {
            // We multiply E/r by r to get E
            E[d]  += E_over_r * potparams.dr[d];
        }
        //std_cout << "Low range field:    Expansion: dr = (" << m_to_bohr*potparams.dr[0] << ", " << m_to_bohr*potparams.dr[1] << ", " << m_to_bohr*potparams.dr[2] << ")   E = (" << E[0] <<", "<< E[1] <<", "<< E[2] << ")    E_over_r = " << E_over_r << "     E_over_r.dr = (" << E_over_r*potparams.dr[0] <<", "<< E_over_r*potparams.dr[1] <<", "<< E_over_r*potparams.dr[2] << ")\n";

    }
}


// ********** End of file ***************************************
