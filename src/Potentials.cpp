/****************************************************************

****************************************************************/

#include <iostream> // cout
#ifdef __SUNPRO_CC
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>  // memset
#include <assert.h>
#else // #ifdef __SUNPRO_CC
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <iostream> // cout
#include <cstring>  // memset
#include <cassert>
#endif // #ifdef __SUNPRO_CC

#include "Constants.hpp"
#include "Potentials.hpp"
#include "General.hpp"
#include "Code_Functions_Declarations.hpp"
#include "Global.hpp"

using namespace libpotentials;

// **************************************************************
// ********** Global variable (to this cpp file) ****************
// ********** These are used as global parameters for  **********
// ********** the different potentials shapes.  They are ********
// ********** set/calculated once at initialization, and ********
// ********** treated as constants afterwards. ******************
// **************************************************************
double s_cutoff_radius_from_input_file;
double sg_m;
double sg_two_m;

double sg_one_over_two_m;
double sg_m_pow_one_over_two_m;
double sg_exp_one_over_two_m;

double hs_min_rad;

const double ps_A = 0.1   * angstrom_to_m;
const double ps_B = 0.45  * angstrom_to_m;
const double ps_C = 0.358 * angstrom_to_m;
const double ps_D = 0.01  * angstrom_to_m*angstrom_to_m;

const double ps_A2          = ps_A * ps_A;
const double ps_B2          = ps_B * ps_B;
const double ps_A_A_minus2B = ps_A * (ps_A - 2.0*ps_B);
const double ps_A_minus_B   = ps_A - ps_B;
const double ps_A_minus_B2  = ps_A_minus_B*ps_A_minus_B;
const double ps_A_minus_B3  = ps_A_minus_B*ps_A_minus_B*ps_A_minus_B;


// **************************************************************
// ********** Local functions prototypes ************************
// **************************************************************
void Set_Coulomb_Field(const double phi12, double E[3], const double dr[3],
                   const double dr2);

// **************************************************************
// ********** Accessible functions implementations **************
// **************************************************************

void Get_r21(double r1[3], double r2[3], double r21[3])
/**
 * r21 is the vector groing from position r2 to position r1
 * r21 = r1 - r2
 * r2 + r21 = r1
 *
 *     2
 *    / \
 *   /   \ dr
 *  /    _\/
 * /------>1
 *
 * @param  r1   Input:  Position r1 [any units]
 * @param  r2   Input:  Position r2 [any units]
 * @param  r21  Output: Vector from to 2 to 1 [same units]
 */
{
    for (int d = 0 ; d < 3 ; d++) r21[d] = r1[d] - r2[d];
}

// **************************************************************
double Get_Distance_Squared(double r1[3], double r2[3])
/**
 * r21 is the vector groing from position r2 to position r1
 * r21 = r1 - r2
 * r2 + r21 = r1
 *
 *     2
 *    / \
 *   /   \ dr
 *  /    _\/
 * /------>1
 *
 * @param  r1    Position r1 [any units]
 * @param  r2    Position r2 [any units]
 * @return r212  Length of vector from to 2 to 1 squared [same units^2]
 */
{
    double r212 = 0.0;
    double r21[3];
    Get_r21(r1, r2, r21);
    for (int d = 0 ; d < 3 ; d++) r212 += r21[d]*r21[d];
    return r212;
}

// **************************************************************
double Get_Distance(double r1[3], double r2[3])
/**
 * r21 is the vector groing from position r2 to position r1
 * r21 = r1 - r2
 * r2 + r21 = r1
 *
 *     2
 *    / \
 *   /   \ dr
 *  /    _\/
 * /------>1
 *
 * @param  r1    Position r1 [any units]
 * @param  r2    Position r2 [any units]
 * @return r212  Length of vector from to 2 to 1 [same units]
 */
{
    return sqrt(Get_Distance_Squared(r1, r2));
}

// **************************************************************
void set_vector_between_particles(
        double r1[3], double r2[3],
        double r21[3], double &r212, double &r, double &one_over_r)
/**
 * r21 is the vector groing from position r2 to position r1
 * r21 = r1 - r2
 * r2 + r21 = r1
 *
 *     2
 *    / \
 *   /   \ dr
 *  /    _\/
 * /------>1
 *
 * @param  r1   Input:  Position 1 [any units]
 * @param  r2   Input:  Position 2 [any units]
 * @param  r21  Output: Vector from to 2 to 1 [same units]
 * @param  r212 Output: Length squared of vector r12 [same units^2]
 * @param  r    Output: Length of vector r12 [same units]
 * @param  one_over_r Output: Inverse of the length of vector r12 [same units^-1]
 */
{
    Get_r21(r1, r2, r21);
    r212       = Get_Distance_Squared(r1, r2);
    r          = sqrt(r212);
    one_over_r = 1.0 / r;
}

// **************************************************************
void Initialize_Simple(const double &minr)
/**
 * Initialize simple cutoff potential
 * @param   minr    Minimum radius where Coulombic potential is
 *                  calculated. For smaller distances, the simple
 *                  cutoff will be used. [m]
 */
{
    s_cutoff_radius_from_input_file = minr;
}

// **************************************************************
void Initialize_SuperGaussian(const int &m)
/**
 * Initialize super-gaussian like potential
 * @param   m       Order of the Super-Gaussian (m=1 is a simple gaussian) [-]
 */
{
    sg_m = m;
    sg_two_m = 2 * m;
    if (sg_m >= 1)
    {
        sg_one_over_two_m         = 1.0 / double(sg_two_m);
        sg_m_pow_one_over_two_m   = pow(sg_m, sg_one_over_two_m);
        sg_exp_one_over_two_m     = exp(sg_one_over_two_m);
    }
    else
    {
        sg_one_over_two_m         = DBL_MAX;
        sg_m_pow_one_over_two_m   = 0.0;
        sg_exp_one_over_two_m     = 0.0;
    }
}

// **************************************************************
void Initialize_HS(const int &m, const double &min_rad)
/**
 * Initialize super-gaussian like potential
 * @param   m       Order of the Super-Gaussian (m=1 is a simple gaussian) [-]
 * @param   min_rad Minimum radius where Herman-Skillman (HS) potential is
 *                  calculated. For smaller distances, a hard cutoffis used. [m]
 */
{
    hs_min_rad = min_rad;
    if (hs_min_rad <= (0.0729 * au_to_si_length))
    {
        std_cout << "##############################################\n";
        DEBUGP("Initialize_HS() called with a minimum too small radius\n");
        std_cout << "The value "<<hs_min_rad * au_to_si_length<<" au ("<<hs_min_rad<<" m) should not be lower than 0.073 a.u. ("<<0.073 * au_to_si_length<<" m)\n";
        std_cout << "Note that this distance should be set in SI units and NOT in atomic units.\n";
        std_cout << "Exiting\n";
        abort();
    }
    Initialize_SuperGaussian(m);
}

// **************************************************************
double Coulomb_Potential(const double kQ, const double r)
/**
 * Calculate the Coulomb potential due to a point charge.
 * @param   kQ      Charge times Coulomb constant [V.m]
 * @param   r       Distance where the potential is to be evaluated [m]
 * @return  Coulombic potential [V]
 */
{
    if (r > DBL_MIN)
        return kQ / r;
    else
        return DBL_MAX;
}

// **************************************************************
void Set_Coulomb_Field(const double phi12, double E[3], const double dr[3],
                   const double r2)
/**
 * Calculate the Coulomb field due to a point charge.
 * @param   phi12   Input:  Potential at position 1 due to charge 2 [V]
 * @param   E       Output: Electrostatic field at position 1 [V.m^-1]
 * @param   dr      Input:  Vector from position 2 to position 1 [m]
 * @param   r2      Input:  Length squared of vector from position 2 to position 1 [m^2]
 */
{
    if (r2 > DBL_MIN)
    {
        // E = r[:] . (V / |r[:]|^2)
        const double phi12_over_dr2 = phi12 / r2;
        for (int d = 0 ; d < 3 ; d++)
            E[d] += dr[d] * phi12_over_dr2;
    }
}

// **************************************************************
double deriv_genericHSfit(const double *par, double x){
    // The 1/2 factor is becuase HS outputs the
    // potential as 2V(x) and that's how they were fit
    return 0.5*(
              par[0]*par[5]*pow(x,par[5]-1)/(pow(x,par[5])-par[1])/(pow(x,par[5])-par[1])
            + par[2]*par[4]/pow(x,par[4]+1)
            + par[6]*par[3]*pow(x,par[6]-1)
        );
}

// **************************************************************
// f(x) = -a/(x^n-b)-b/x^m +d*x^o where the parameters are passed alphabetically
double genericHSfit(const double *par, double x){
    //the 1/2 factor is becuase HS outputs the potential as 2V and that's how they were fit
    return -0.5*(
            - par[0]/(pow(x,par[5])-par[1])
            - par[2]/pow(x,par[4])
            + par[3]*pow(x,par[6])
        );
}

// **************************************************************
// Function pointers for...
// ...setting the parameters of the potential/field calculation
void   (*Potentials_Set_Parameters)(void *p1, void *p2, potential_paramaters &potparams) = NULL;
// ...calculating the potential
double (*Calculate_Potential)(      void *p1, void *p2, potential_paramaters &potparams) = NULL;
// ...setting the electric field
void   (*Set_Field)(                void *p1, void *p2, potential_paramaters &potparams, double &phi, double E[3]) = NULL;

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
    Check_if_initialized();

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);

    // we only add field if other particle has a charge not equal 0
    double Q2 = Get_Charge(p2);
    int charge_state2 = Get_Charge_State(p2);
    if (charge_state2 != 0)
    {
        potparams.kQ2 = one_over_4Pieps0 * Q2;

            potparams.cutoff_radius = tmp_get_shieldr_2(charge_state2, Get_Charge_State(p1));
//         else
//             potparams.cutoff_radius = 0.0;
    }
    else
    {
        potparams.cutoff_radius = 0.0;
        potparams.kQ2           = 0.0;
    }
}

// **************************************************************
double Calculate_Potential_Cutoff_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_initialized();

    double potential;

    if (potparams.r > potparams.cutoff_radius)
        potential = Coulomb_Potential(potparams.kQ2, potparams.r);
    else
        potential = Coulomb_Potential(potparams.kQ2, potparams.cutoff_radius);

    return potential;
}

// **************************************************************
void Set_Field_Cutoff_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3])
{
    Check_if_initialized();

    if (potparams.r > potparams.cutoff_radius)
        Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);
    // Simple cutoff means the potential is constant for distance
    // less then the cutoff radius. If the potential is constant, then
    // their is not field.
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
    Check_if_initialized();

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);

    // we only add field if other particle has a charge not equal 0
    double Q2 = Get_Charge(p2);
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
        if (Q2 > DBL_MIN)
        {   // Positive charge (ion), positive potential
            potparams.B = base_pot_well_depth*charge_state2;
            //potparams.B *= eV_to_J * (1.0 / e0); // *= 1.0
        }
        else
        {   // Negative charge (electron), negative potential
            potparams.B = -base_pot_well_depth;
            //potparams.B *= (1.0 / e0) * eV_to_J; // *= 1.0
        }

        potparams.kQ2 = one_over_4Pieps0 * Q2;
        potparams.kQ2_over_B = potparams.kQ2 / potparams.B;

        // Radius where the Coulomb potential and its first derivative are
        // equal to a potential in a -A r^2 + B shape and its first derivative.
        potparams.cutoff_radius = three_over_two * potparams.kQ2_over_B;
    }
    else
    {
        potparams.cutoff_radius = 0.0;
        potparams.kQ2           = 0.0;
        potparams.kQ2_over_B    = 0.0;
        potparams.B             = 0.0;
    }
}

// **************************************************************
double Calculate_Potential_Cutoff_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams)
/**
 * Calculate the electrostatic potential energy "phi1" of particle "p1"
 * due to another particle/cell "p2" at distance sqrt(dr2)
 */
{
    Check_if_initialized();

    double phi12;   // Electrostatic potential

    if (potparams.r <= potparams.cutoff_radius)
    {
        // const double A  = ( 4.0 * B*B*B) / (27.0 * kQ*kQ);
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
    double &phi, double E[3])
{
    Check_if_initialized();

    if (potparams.r > potparams.cutoff_radius)
        Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);
    else
    {
        // The derivative of a -Ar^2+B w.r. to r is:
        const double diff_x2 = -2.0*potparams.h_A*potparams.r;

        // We have the norm of the gradient of the potential (diff_x2),
        // we need to multiply this by the unit vector, to get
        // the electric field.
        double unit_dr[3];
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
    Check_if_initialized();

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);

    // we only add field if other particle has a charge not equal 0
    double Q2 = Get_Charge(p2);
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
        if (Q2 > DBL_MIN)
        {   // Positive charge (ion), positive potential
            potparams.B = base_pot_well_depth*charge_state2;
            //potparams.B *= eV_to_J * (1.0 / e0); // *= 1.0
        }
        else
        {   // Negative charge (electron), negative potential
            // B is the (negative of the)
            // ionization potential of the ion.
            potparams.B = -base_pot_well_depth;
            //potparams.B *= eV_to_J * (1.0 / e0); // *= 1.0
        }

        potparams.kQ2 = one_over_4Pieps0 * Q2;
        potparams.kQ2_over_B = potparams.kQ2 / potparams.B;

        // Radius where the Coulomb potential and its first derivative are
        // equal to a the super-gaussian potential.
        potparams.cutoff_radius = potparams.kQ2_over_B * sg_exp_one_over_two_m;

        // Width of Super-Gaussian
        if (potparams.r <= potparams.cutoff_radius)
        {
            potparams.sg_sigma = potparams.kQ2_over_B
                                    * sg_m_pow_one_over_two_m
                                    * sg_exp_one_over_two_m;
            potparams.sg_r_over_sigma_two_m = pow(potparams.r/potparams.sg_sigma, sg_two_m);
            potparams.sg_exp_half_r_over_sigma_two_m =
                    exp( -0.5 * potparams.sg_r_over_sigma_two_m );
        }
    }
    else
    {
        potparams.cutoff_radius = 0.0;
        potparams.kQ2           = 0.0;
        potparams.kQ2_over_B    = 0.0;
        potparams.B             = 0.0;
    }
}

// **************************************************************
double Calculate_Potential_Cutoff_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_initialized();

    double phi12;   // Electrostatic potential

    if (potparams.r <= potparams.cutoff_radius)
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
    double &phi, double E[3])
{
    Check_if_initialized();

    if (potparams.r > potparams.cutoff_radius)
        Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);
    else
    {
        // The derivative of a super gaussian w.r. to r is:
        double diff_sg = -(potparams.B*sg_m*potparams.one_over_r)
                            * potparams.sg_r_over_sigma_two_m
                            * potparams.sg_exp_half_r_over_sigma_two_m;

        // We have the norm of the gradient of the potential (diff_sg),
        // we need to multiply this by the unit vector, to get
        // the electric field.
        double unit_dr[3];
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

// **************************************************************

void Potentials_Set_Parameters_HS_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams)
/**
 * p1 is the particle of interest, p2 is "the other" particle
 * acting on p1.
 *
 */
{
    Check_if_initialized();

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);

    potparams.kQ2           = one_over_4Pieps0 * Get_Charge(p2);
    potparams.hs_type2      = Get_Type(p2);
    potparams.hs_cs2        = Get_Charge_State(p2);

    // Default values
    potparams.cutoff_radius = 0.073 * au_to_si_length;
    potparams.kQ2_over_B    = 0.0;

    // If p2 is not a body (thus it is a cell), no approximation is needed
    // because it should be far away, so it should already be approximated.
    // FIXME: HS+BODY
    //if (potparams.hs_type2 == BODY)
    {
        // If we have an electron, set the correct parameters for the
        // potential cutoff. All other cases (neutral atom or ion) are
        // treated as HS potential (neutral to 7+) and coulomb (8+ and up)
        if ( (potparams.hs_cs2 == -1) or (potparams.hs_cs2 >= 8) )
        {

            // Negative charge (electron), negative potential
            // B is the (negative of the)
            // ionization potential of the ion.
            potparams.B = -base_pot_well_depth;

            potparams.kQ2_over_B = potparams.kQ2 / potparams.B;

            // Radius where the Coulomb potential and its first derivative are
            // equal to a the super-gaussian potential.
            potparams.cutoff_radius = potparams.kQ2_over_B * sg_exp_one_over_two_m;

            // Width of Super-Gaussian
            potparams.sg_sigma = potparams.kQ2_over_B
                                    * sg_m_pow_one_over_two_m
                                    * sg_exp_one_over_two_m;

            potparams.sg_r_over_sigma_two_m = pow(potparams.r / potparams.sg_sigma, sg_two_m);
            potparams.sg_exp_half_r_over_sigma_two_m =
                                exp( -0.5 * potparams.sg_r_over_sigma_two_m );
        }
    }
}

// **************************************************************
double Calculate_Potential_Cutoff_HS_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams)
/**
 * Output in Volts units. The potential outputs the potential energy
 * (the potential in volts needs to be multiplied
 * by charge_state*e, charge of the electron)
 * Volts = J/C and the fits are in atomic units (Hartree) thus:
 * Hatree * Eh_to_eV(eV/Hartree) * e0(J/eV) * 1/e0(1/C) = Hatree * Eh_to_eV (J/C=Volts)
 * hs_min_rad SHOULDN'T be < 0.073 since the functions aren't fitted there
 */
{
    Check_if_initialized();

    double phi12 = 0.0;   // Electrostatic potential

    // Fits are in atomic units
    const double distance_au = potparams.r * si_to_au_length;

    // Ions are given a potential inside the electron cloud.
    // The last two doubles of the fit_lessthan_R array is the range
    // the fit is valid for below the smallest range we have a simple
    // hard cutoff V(r<r_0) = V(r_0)

    // ___________________________________________________________________ ...
    // |        |             |                |       |
    // |  CP    |     R1      |     R2         |  R3   |  Coulomb
    // |________|_____________|________________|_______|__________________ ...

    const int cs = potparams.hs_cs2;

    // FIXME: HS+BODY
    //if ((cs == 0) && ( potparams.hs_type2 == BODY))
    if (cs == 0)
    {
        if (distance_au >= fit_lt_R3[0][8])    /* In Coulomb */
        {
            // Potential outside the electron cloud goes to 0
            // exponentially using f(x)=h*exp(-v*x+k)
            phi12 = 1.93775072943628 * exp(
                -0.533297816151*distance_au
                -0.7486357665822807
            ); /*std_cout << "0. HS in C\n";*/
        }
        else if (   (distance_au <  fit_lt_R3[0][8]) &&
                    (distance_au >= fit_lt_R3[0][7]))               /* In R3 */
        {
            phi12 = genericHSfit(&(fit_lt_R3[0][0]),distance_au);
            /*std_cout << "1. HS in R3\n";*/
        }
        else if (   (distance_au <  fit_lt_R2[0][8]) &&
                    (distance_au >= fit_lt_R2[0][7]))               /* In R2 */
        {
            phi12 = genericHSfit(&(fit_lt_R2[0][0]),distance_au);
            /*std_cout << "2. HS in R2\n";*/
        }
        else if (   (distance_au <  fit_lt_R1[0][8]) &&
                    (distance_au >= fit_lt_R1[0][7]))               /* In R1 */
        {
            phi12 = genericHSfit(&(fit_lt_R1[0][0]),distance_au);
            /*std_cout << "3. HS in R1\n";*/
        }
        else    // Hard cutoff
        {
            phi12 = genericHSfit(&(fit_lt_R1[0][0]), hs_min_rad);
            /*std_cout << "4. HS inside\n";*/
        }
        phi12 *= Eh_to_eV;
    }
    // FIXME: HS+BODY
    //else if ( (cs < 10) && (cs > 0) && ( potparams.hs_type2 == BODY) )
    else if ( (cs < 10) && (cs > 0) )
    { // If charge state is between 0 and 9 (inclusive)...

        if          (distance_au >= fit_lt_R3[cs][8])    /* In Coulomb */
        {
            phi12 = (potparams.kQ2 / potparams.r);   // Outside electron cloud: Coulombic pot.
        }
        else if (   (distance_au <  fit_lt_R3[cs][8]) &&
                    (distance_au >= fit_lt_R3[cs][7]))   /* In R3 */
        {
            phi12 = genericHSfit(&(fit_lt_R3[cs][0]),distance_au) * Eh_to_eV;
        }
        else if (   (distance_au < fit_lt_R2[cs][8]) &&
                    (distance_au >= fit_lt_R2[cs][7]))   /* In R2 */
        {
            phi12 = genericHSfit(&(fit_lt_R2[cs][0]),distance_au) * Eh_to_eV;
        }
        else if (   (distance_au <  fit_lt_R1[cs][8]) &&
                    (distance_au >= fit_lt_R1[cs][7]))   /* In R1 */
        {
            phi12 = genericHSfit(&(fit_lt_R1[cs][0]),distance_au) * Eh_to_eV;
        }
        else                                   /* In CP (constant potential) */
        {
            phi12 = genericHSfit(&(fit_lt_R1[cs][0]),hs_min_rad) * Eh_to_eV;
        }
    }
    else if( abs(cs) >= 1 )
    {
        if (potparams.r <= potparams.cutoff_radius)
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
    }
//     std_cout << "HS: phi12 = " << phi12 << "\n";

    return phi12;
}

// **************************************************************
void Set_Field_Cutoff_HS_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3])
{
    Check_if_initialized();

    // Fits are in atomic units
    const double distance_au = potparams.r * si_to_au_length;

    const int cs = potparams.hs_cs2;
    double Ef;

    // Ions are given a potential inside the electron cloud.
    // The last two doubles of the fit_lessthan_R array is the range
    // the fit is valid for below the smallest range we have a simple
    // hard cutoff V(r<r_0) = V(r_0)
    // FIXME: HS+BODY
    //if ((cs == 0) && ( potparams.hs_type2 == BODY))
    if (cs == 0)
    {
        if (distance_au >= fit_lt_R3[0][8])
        {
            // Potential outside the electron cloud goes to 0
            // exponentially using f(x)=h*exp(-v*x+k)
            Ef = 0.272*(    -1.93775072943628
                    * (-0.533297816151*distance_au - 0.7486357665822807) *(-0.533297816151)
                    * exp(-0.533297816151*distance_au - 0.7486357665822807)
                    )* Eh_to_eV;
            double unit_dr[3];
//             MULVS(unit_dr, dr, one_over_distance);  // Calculate unit vector
//             ADDMULVS(E, unit_dr, Ef/a0);
            for (int d = 0 ; d < 3 ; d++)
            {
                unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
                E[d]  += unit_dr[d] * (Ef/a0);
            }

        }
        else if (   (distance_au < fit_lt_R3[0][8]) &&
                    (distance_au >= fit_lt_R3[0][7])){
            Ef = deriv_genericHSfit(&(fit_lt_R3[0][0]),distance_au) * Eh_to_eV;
            double unit_dr[3];
//             MULVS(unit_dr, dr, one_over_distance);  // Calculate unit vector
//             ADDMULVS(E, unit_dr, Ef/a0);
            for (int d = 0 ; d < 3 ; d++)
            {
                unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
                E[d]  += unit_dr[d] * (Ef/a0);
            }
        }
        else if (   (distance_au < fit_lt_R2[0][8]) &&
                    (distance_au >= fit_lt_R2[0][7])){
            Ef = deriv_genericHSfit(&(fit_lt_R2[0][0]),distance_au) * Eh_to_eV;
            double unit_dr[3];
//             MULVS(unit_dr, dr, one_over_distance);  // Calculate unit vector
//             ADDMULVS(E, unit_dr, Ef/a0);
            for (int d = 0 ; d < 3 ; d++)
            {
                unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
                E[d]  += unit_dr[d] * (Ef/a0);
            }
        }
        else if (   (distance_au < fit_lt_R1[0][8]) &&
                    (distance_au >= fit_lt_R1[0][7])){
            Ef = deriv_genericHSfit(&(fit_lt_R1[0][0]),distance_au) * Eh_to_eV;
            double unit_dr[3];
//             MULVS(unit_dr, dr, one_over_distance);  // Calculate unit vector
//             ADDMULVS(E, unit_dr, Ef/a0);
            for (int d = 0 ; d < 3 ; d++)
            {
                unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
                E[d]  += unit_dr[d] * (Ef/a0);
            }
         }
        //else
            // No field because of hard cutoff (-grad(constant) = 0 )
    }
    // FIXME: HS+BODY
    //else if ( (cs < 8) && (cs > 0) && ( potparams.hs_type2 == BODY) )
    else if ( (cs < 8) && (cs > 0) )
    { //for 8+ use this instead of Coulomb potential

        // Find what range the radial distance_au is in and use the appropriate potential
        // ___________________________________________________________________ ...
        // |        |             |                |       |
        // |  CP    |     R1      |     R2         |  R3   |  Coulomb
        // |________|_____________|________________|_______|__________________ ...
        if (distance_au >= fit_lt_R3[cs][8])
        {
            //potential outside the electron cloud is Coulombic
            double unit_dr[3];
//             MULVS(unit_dr, dr, one_over_distance);  // Calculate unit vector
//             ADDMULVS(E, unit_dr, kQ2 / distance / distance);
            for (int d = 0 ; d < 3 ; d++)
            {
                unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
                E[d]  += unit_dr[d] * (potparams.kQ2 / potparams.r / potparams.r);
            }

        }
        else if (   (distance_au <  fit_lt_R3[cs][8]) &&
                    (distance_au >= fit_lt_R3[cs][7])
        ) {
            double unit_dr[3];
//             MULVS(unit_dr, dr, one_over_distance);  // Calculate unit vector
//             ADDMULVS(E, unit_dr,
//                 deriv_genericHSfit(&(fit_lt_R3[cs][0]),distance_au) * Eh_to_eV/a0);
            double tmp = deriv_genericHSfit(&(fit_lt_R3[cs][0]),distance_au) * Eh_to_eV/a0;
            for (int d = 0 ; d < 3 ; d++)
            {
                unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
                E[d]  += unit_dr[d] * tmp;
            }
        }
        else if (   (distance_au <  fit_lt_R2[cs][8]) &&
                    (distance_au >= fit_lt_R2[cs][7])
        ) {
            double unit_dr[3];
//             MULVS(unit_dr, dr, one_over_distance);  // Calculate unit vector
//             ADDMULVS(E, unit_dr,
//                 deriv_genericHSfit(&(fit_lt_R2[cs][0]),distance_au) * Eh_to_eV/a0);
            double tmp = deriv_genericHSfit(&(fit_lt_R2[cs][0]),distance_au) * Eh_to_eV/a0;
            for (int d = 0 ; d < 3 ; d++)
            {
                unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
                E[d]  += unit_dr[d] * tmp;
            }
        }
        else if (   (distance_au <  fit_lt_R1[cs][8]) &&
                    (distance_au >= fit_lt_R1[cs][7])
        ) {
            double unit_dr[3];
//             MULVS(unit_dr, dr, one_over_distance);  // Calculate unit vector
//             ADDMULVS(E, unit_dr,
//                 deriv_genericHSfit(&(fit_lt_R1[cs][0]),distance_au) * Eh_to_eV/a0);
            double tmp = deriv_genericHSfit(&(fit_lt_R1[cs][0]),distance_au) * Eh_to_eV/a0;
            for (int d = 0 ; d < 3 ; d++)
            {
                unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
                E[d]  += unit_dr[d] * tmp;
            }
        }
//         else
            // No field because of hard cutoff (-grad(constant) = 0 )
    }
    else if( abs(cs) >= 1 )
    {
        if (potparams.r > potparams.cutoff_radius)
        {
            Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);
        }
        else
        {
            potparams.sg_r_over_sigma_two_m = pow(potparams.r / potparams.sg_sigma, sg_two_m);
            potparams.sg_exp_half_r_over_sigma_two_m =
                            exp( -0.5 * potparams.sg_r_over_sigma_two_m );

            // The derivative of a super gaussian w.r. to r is:
            double diff_sg =
                -(potparams.B*sg_m*potparams.one_over_r) * potparams.sg_r_over_sigma_two_m
                * potparams.sg_exp_half_r_over_sigma_two_m;

            // We have the norm of the gradient of the potential (diff_sg),
            // we need to multiply this by the unit vector, to get
            // the electric field.
            double unit_dr[3];
//             MULVS(unit_dr, dr, one_over_distance);  // Calculate unit vector
//             ADDMULVS(E, unit_dr, -diff_sg);         // Add to the electric field
//                                                     // the gradient of the  potential.
            for (int d = 0 ; d < 3 ; d++)
            {
                unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
                E[d]  += unit_dr[d] * -diff_sg;
            }
        }
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
    Check_if_initialized();

//   if (Get_Id(p1) == Get_Id(p2)){
/*  std::cout << "Id1=" << Get_Id(p1) << " Id2=" << Get_Id(p2)
    << "  potparams.r=" << potparams.r << std::endl;
   std::cout << " Type1=" << Get_Type(p1) << " Type1=" << Get_Type(p2) <<  std::endl;
   std::cout << " Get_Position(p1)=" << Get_Position(p1)[0] << "," << Get_Position(p1)[1] << "," << Get_Position(p1)[2]
<< "\nGet_Position(p1)=" << Get_Position(p2)[0] << "," << Get_Position(p2)[1] << "," << Get_Position(p2)[2]
<<  std::endl;*/
//    }
//    assert(Get_Id(p1) != Get_Id(p2));

    // Well depth
    //const int depth = 0;

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);
#ifdef YDEBUG
    if (potparams.r <= 1.0e-200 || isnan(potparams.r))
    {
        Print_Partiles(p1, NULL, 1, 1, 0, 0);
        Print_Partiles(p2, NULL, 1, 1, 0, 0);
        assert(potparams.r > 1.0e-200);
        abort();
    }
#endif // #ifdef YDEBUG
    assert(potparams.r > 1.0e-200);

    // we only add field if other particle has a charge not equal 0
    double Q = Get_Charge(p2);
    int distribution1_charge_state = 1,distribution2_charge_state = 1;
    int charge_state2 = Get_Charge_State(p2);
    int charge_state1 = Get_Charge_State(p1);
    int factor = charge_state2;
    if( charge_state2 != 0 )
    {
      //double Ip =element.IpsLowest[abs(charge_state2)];
      double Ip = base_pot_well_depth;
      //charges are not equal
      //take the higher charge as the distribution
      //Thus set the parameters for the guassian
      //using the values of the larger charge
      //NOTE: it will always do this for p2=electron
      // unless p1 is also an electron
        if ((charge_state2 < charge_state1) && (charge_state1 != 0)){
          Ip = base_pot_well_depth;
          if (charge_state2 < 0)
            Q = -Get_Charge(p1);
          else
            Q = Get_Charge(p1);
          //save the charge states so kQ2 = KQ1 * charge_state2/charge_state1
          //which gives the correct scale. Needed for ion-ion
          distribution1_charge_state = charge_state1;
          distribution2_charge_state = abs(charge_state2);
          factor = charge_state1;
        }
        //Make the Ip linearly deeper for ions.
        //not used for e- e-
        Ip*=(factor);

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
            potparams.B = -base_pot_well_depth;
          else //else p1 is an ion
            potparams.B = -Ip;
        }

        potparams.kQ2 = one_over_4Pieps0 * Q;
        potparams.kQ2_over_B = potparams.kQ2 / potparams.B;

        // Make sure the other particle's charge is normalized by the current
        // particle's charge state so forces are symmetric between the two.
        potparams.kQ2 /= double(distribution1_charge_state);
        potparams.kQ2 *= double(distribution2_charge_state);
        assert(!isnan(potparams.kQ2));
#ifndef __SUNPRO_CC
        assert(!isinf(potparams.kQ2));
#endif // #ifndef __SUNPRO_CC

        potparams.gd_sigma = potparams.kQ2_over_B * sqrt_2_over_pi;
/*        std_cout << "Id(p1)="<<Get_Id(p1)<<"  B="<<potparams.B<<"  Cs(p1)="<<Get_Charge_State(p1)<<" Cs(p2)="<<Get_Charge_State(p2)<<" kQ2_over_B="<<potparams.kQ2_over_B<<" well="<<base_pot_well_depth<<"\n";*/
        // Radius where the Coulomb potential and its first derivative are
        // equal to a the gaussian charge distribution potential.
        potparams.cutoff_radius = 4.0 * potparams.gd_sigma;

    }
    else
    {
        potparams.cutoff_radius = 0.0;
        potparams.kQ2           = 0.0;
        potparams.kQ2_over_B    = 0.0;
        potparams.B             = 0.0;
        potparams.gd_sigma      = 0.0;
    }
}

// **************************************************************
double Calculate_Potential_Cutoff_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_initialized();

    double phi12=0.0;   // Electrostatic potential

    if (potparams.r <= potparams.cutoff_radius)
    {
        // If the distance between two bodys is less than the shielding
        // radius, we use the special potential calculated from a
        // gaussian charge distribution instead of the Coulomb potential.

        const double r_over_sigma_sqrt_2 = potparams.r/ (sqrt_2 * potparams.gd_sigma);

        // Get partial potential
        phi12 = potparams.kQ2 * (potparams.one_over_r * my_erf(r_over_sigma_sqrt_2));
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
void Set_Field_Cutoff_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3])
{
    Check_if_initialized();

    if (potparams.r > potparams.cutoff_radius)
        Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);
    // A debug statement to catch very small distances where the true formular should be inf - inf
    //This rarely cancels well and thus divergers giving a huge kick to the particles so
    //we test for small r and set the field to 0 (a true expansion would be costly for little return)
    else if (potparams.r < 1.48e-12 ) //r-vale determined by code trial and error for base_potential of 0.5
    {
      return; //don''t add any field making it effectively 0
    }
    else
    {
        const double r_over_sigma_sqrt_2 = potparams.r / (sqrt_2 * potparams.gd_sigma);
        // The radial component of the gradient of the potential w.r. to r is:
        double grad_cd = potparams.kQ2 * (
            (my_erf(r_over_sigma_sqrt_2) / potparams.r2)
            - sqrt_2_over_pi * exp(- r_over_sigma_sqrt_2 * r_over_sigma_sqrt_2)
                / (potparams.gd_sigma * potparams.r)
        );
        // We have the radial component of the gradient of the potential,
        // we need to multiply this by the unit vector, to get
        // the electric field.
        double unit_dr[3];
//         MULVS(unit_dr, dr, one_over_distance);  // Calculate unit vector
//         ADDMULVS(E, unit_dr, grad_cd);          // Add to the electric field
//                                                 // the gradient of the potential.
        for (int d = 0 ; d < 3 ; d++)
        {
            unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
            E[d]  += unit_dr[d] * grad_cd;
        }
    }
}

// **************************************************************
// **************************************************************

// **************************************************************
void Potentials_Set_Parameters_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_initialized();

    potparams.sym_cs1 = Get_Charge_State(p1);
    potparams.sym_cs2 = Get_Charge_State(p2);

    // It does not make sense to calculate effect of particle 2 if
    // its charge is 0

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);
    if (potparams.sym_cs2 == 0) return;
#ifdef YDEBUG
    // Failsafe
    if (potparams.r <= 1.0e-200 || isnan(potparams.r))
    {
        std_cout << "Error in Potentials_Set_Parameters_ChargeDistribution_Symmetric()\n";
        std_cout << "Particles p1 ("<<p1<<") and p2 ("<<p2<<") are too close.\n";
        std_cout << "Distance = " << potparams.r * si_to_au_length << "\n";
        Print_Partiles(p1, NULL, 1, 1, 0, 0);
        Print_Partiles(p2, NULL, 1, 1, 0, 0);
        assert(potparams.r > 1.0e-200);
        abort();
    }
#endif // #ifdef YDEBUG
    assert(potparams.r > 1.0e-200);

    if (potparams.sym_cs1 != 0 &&  potparams.sym_cs2 != 0)
    {
        // FIXME: Should it be "sym_cs - 1" or plain "sym_cs" for the IP?
        // FIXME: The jump in energy before and after ionization occurs
        // seems to be affected by the depth of the well.
        //double Ip1 = element.IpsLowest[max(potparams.sym_cs1-1, 0)];
        //double Ip2 = element.IpsLowest[max(potparams.sym_cs2-1, 0)];
        double Ip1 = base_pot_well_depth*2.0;
        double Ip2 = base_pot_well_depth*2.0;

        potparams.kQ2 = one_over_4Pieps0 * double(potparams.sym_cs2) * e0;
//         potparams.sym_sigma1 = one_over_4Pieps0 * double(abs(potparams.sym_cs1)) * e0 * sqrt_2_over_pi / Ip1;
//         potparams.sym_sigma2 = one_over_4Pieps0 * double(abs(potparams.sym_cs2)) * e0 * sqrt_2_over_pi / Ip2;
        potparams.sym_sigma1 = one_over_4Pieps0 * e0 * sqrt_2_over_pi / Ip1;
        potparams.sym_sigma2 = one_over_4Pieps0 * e0 * sqrt_2_over_pi / Ip2;
    }
    else
    {
        potparams.kQ2         = 0.0;
        potparams.sym_sigma1  = 0.0;
        potparams.sym_sigma2  = 0.0;
    }
}

// **************************************************************
double Calculate_Potential_Cutoff_ChargeDistribution_Symmetric(
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
    Check_if_initialized();

    double potential = 0.0;
    // It does not make sense to calculate effect of particle 2 if
    // its charge is 0
    if (potparams.sym_cs2 != 0)
    {
        if (potparams.sym_cs1 == 0)
        {
            potential = Coulomb_Potential(potparams.kQ2, potparams.one_over_r);
        }
        else
        {
            potential = (potparams.kQ2 * potparams.one_over_r)
                * my_erf(
                    potparams.r / sqrt(
                        2.0*potparams.sym_sigma1*potparams.sym_sigma1
                        + 2.0*potparams.sym_sigma2*potparams.sym_sigma2
                    )
                );
        }
    }
    return potential;
}

// **************************************************************
void Set_Field_Cutoff_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3])
{
    Check_if_initialized();

    // It does not make sense to calculate effect of particle 2 if
    // its charge is 0
    if (potparams.sym_cs2 != 0)
    {
        if (potparams.sym_cs1 == 0)
        {
            Set_Coulomb_Field(phi, E, potparams.dr, potparams.r2);
        }
        else
        {
            const double r_over_sqrt = potparams.r / sqrt(
                      (2.0*potparams.sym_sigma1*potparams.sym_sigma1)
                    + (2.0*potparams.sym_sigma2*potparams.sym_sigma2)
                );
            const double absE = potparams.kQ2 * potparams.one_over_r * potparams.one_over_r * (
                my_erf(r_over_sqrt) - two_over_sqrt_Pi * r_over_sqrt * exp(-r_over_sqrt*r_over_sqrt)
            );

            double f;
            double unit_vector;
            for (int d = 0 ; d < 3 ; d++)
            {
                unit_vector = potparams.dr[d] * potparams.one_over_r;
                f = unit_vector * absE;
                assert(!isnan(f));
#ifndef __SUNPRO_CC
                assert(!isinf(f));
#endif // #ifndef __SUNPRO_CC
                assert(f < 1.0e200);
                E[d] += f;
            }

        }
    }
}

// **************************************************************
void Potentials_Set_Parameters_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_initialized();

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);
    potparams.kQ2 = one_over_4Pieps0 * Get_Charge(p2);
}

// **************************************************************
double Calculate_Potential_Cutoff_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_initialized();

    const double one_over_r_plus_alpha = 1.0 / (potparams.r + sc_alpha);
    const double phi = potparams.kQ2 * (2.0 * potparams.r + sc_alpha) * one_over_r_plus_alpha * one_over_r_plus_alpha * 0.5;
    return phi;
}

// **************************************************************
void Set_Field_Cutoff_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3])
{
    Check_if_initialized();

    const double Eabs = potparams.kQ2 / pow(potparams.r + sc_alpha,3);
    for (int d = 0 ; d < 3 ; d++)
        E[d] += potparams.dr[d] * Eabs;
}

// **************************************************************
void Potentials_Set_Parameters_PseudoParticles(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_initialized();

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);
    potparams.kQ2 = one_over_4Pieps0 * Get_Charge(p2);
}

// **************************************************************
double Calculate_Potential_PseudoParticles(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_initialized();

    double phi = - potparams.kQ2 * (
          ( ps_A2 / ( 2.0 * ps_A_minus_B * pow(ps_A + potparams.r, 2) ) )
        - ( ps_A_A_minus2B / ( ps_A_minus_B2 * (ps_A + potparams.r) ) )
        - ( ps_B2 * log( (ps_A+potparams.r)/(ps_B+potparams.r) ) / ps_A_minus_B3 )
        - ( ps_C / ( 2.0 * ( pow(potparams.r, 2) + ps_D) ) )
    );
    return phi;
}

// **************************************************************
void Set_Field_PseudoParticles(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3])
{
    Check_if_initialized();

    const double Eabs = potparams.kQ2 * (
          ( potparams.r / ( pow(potparams.r + ps_A, 3) * (potparams.r + ps_B) ) )
        + ( ps_C / pow( pow(potparams.r, 2) + ps_D, 2) )
    );
    for (int d = 0 ; d < 3 ; d++)
        E[d] += potparams.dr[d] * Eabs;
}

// **************************************************************
double tmp_get_shieldr(const int chg_st, const char *message)
{
    if( chg_st == 0 )
    {
        std_cout << "error in gecutoff_radius(): chg_st == 0 from="<<message<<"\n";
        abort();
    }

    return (s_cutoff_radius_from_input_file + 0.5 * (abs(chg_st)-1)*one_over_seven * a0);
//    return(shieldr);
}

// **************************************************************
double tmp_get_shieldr_2(const int chg_st_1, const int chg_st_2)
{
    if( chg_st_1 == 0 )
    {
        std_cout << "error in gecutoff_radius(): chg_st == 0 from=gecutoff_radius_2\n";
        abort();
    }

    int abs_chg_st_1 = abs(chg_st_1);
    int abs_chg_st_2 = abs(chg_st_2);
    int max_chg_st = (abs_chg_st_1 > abs_chg_st_2 ) ? abs_chg_st_1 : abs_chg_st_2;

    return tmp_get_shieldr(max_chg_st,"gecutoff_radius_2");
}

// ********** End of file ***************************************
