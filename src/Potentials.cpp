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
#include "Global.hpp"
#include "Vectors.hpp"

using namespace libpotentials;

bool is_libpotentials_initialized = false;
bool USING_HS;

// **************************************************************
// ********** Global variable (to this cpp file) ****************
// ********** These are used as global parameters for  **********
// ********** the different potentials shapes.  They are ********
// ********** set/calculated once at initialization, and ********
// ********** treated as constants afterwards. ******************
// **************************************************************
fdouble s_cutoff_radius_from_input_file;
fdouble sg_m;
fdouble sg_two_m;

fdouble sg_one_over_two_m;
fdouble sg_m_pow_one_over_two_m;
fdouble sg_exp_one_over_two_m;

std::vector<fdouble> hs_min_rad;
std::vector<LookUpTable<fdouble> > hs_lut_potential;
const int max_hs_cs = 7;

const fdouble ps_A = fdouble(0.1)   * angstrom_to_m;
const fdouble ps_B = fdouble(0.45)  * angstrom_to_m;
const fdouble ps_C = fdouble(0.358) * angstrom_to_m;
const fdouble ps_D = fdouble(0.01)  * angstrom_to_m*angstrom_to_m;

const fdouble ps_A2          = ps_A * ps_A;
const fdouble ps_B2          = ps_B * ps_B;
const fdouble ps_A_A_minus2B = ps_A * (ps_A - fdouble(2.0)*ps_B);
const fdouble ps_A_minus_B   = ps_A - ps_B;
const fdouble ps_A_minus_B2  = ps_A_minus_B*ps_A_minus_B;
const fdouble ps_A_minus_B3  = ps_A_minus_B*ps_A_minus_B*ps_A_minus_B;

// ********** Herman-Skillman (HS) potential fit parameters *****
// The fit function is f(x)=-a/(x^n-b)-b/x^m +d*x^o and the fit
// parameters are in alphabetical order there is a different
// array for each radial distance, where the cutoff radial
// distance is the last element of the array.
const fdouble fit_lt_R1[max_hs_cs+1][9] = {
/*                                                                                                                                                  r_min               r_max */
    {fdouble(-39.3117),  fdouble(-0.23822),       fdouble(1137.15),    fdouble(1093.87),    fdouble(0.926033),   fdouble(1.35102),    fdouble(-0.902534),  fdouble(0.073),  fdouble(1.0)},
    {fdouble(-50.9699),  fdouble(-0.249349),      fdouble(1190.93),    fdouble(1137.92),    fdouble(0.934615),   fdouble(1.26191),    fdouble(-0.915538),  fdouble(0.073),  fdouble(1.0)},
    {fdouble(-103318.0), fdouble(-0.0025085),     fdouble(109884.0),   fdouble(6808.38),    fdouble(0.452195),   fdouble(0.453217),   fdouble(-0.452051),  fdouble(0.073),  fdouble(0.6)},
    {fdouble(-103309.0), fdouble(-0.00222854),    fdouble(109893.0),   fdouble(6799.44),    fdouble(0.462727),   fdouble(0.46361),    fdouble(-0.462679),  fdouble(0.073),  fdouble(0.75)},
    {fdouble(-106539.0), fdouble(-0.00253429),    fdouble(106375.0),   fdouble(90.6738),    fdouble(0.576876),   fdouble(0.577581),   fdouble(-1.31088),   fdouble(0.073),  fdouble(0.75)},
    {fdouble(-106552.0), fdouble(-0.00285323),    fdouble(106363.0),   fdouble(97.8488),    fdouble(0.555262),   fdouble(0.556083),   fdouble(-1.26473),   fdouble(0.073),  fdouble(0.75)},
    {fdouble(-106572.0), fdouble(-0.00333483),    fdouble(106342.0),   fdouble(106.134),    fdouble(0.523154),   fdouble(0.524146),   fdouble(-1.19394),   fdouble(0.073),  fdouble(0.75)},
    {fdouble(-156137.0), fdouble(-9.52875),       fdouble(16279.0),    fdouble(-30.664),    fdouble(0.00362412), fdouble(4.27026),    fdouble(-1.32769),   fdouble(0.073),  fdouble(0.35)}
};
const fdouble fit_lt_R2[max_hs_cs+1][9] = {
/*                                                                                                                                                                                                              r_min           r_max */
    {fdouble(-106456.523613218), fdouble(-0.00434541093317553),   fdouble(106457.47661029),    fdouble(449.688681389621),   fdouble(1.05645523648719),   fdouble(1.05644674944298),   fdouble(-2.10055725950707),  fdouble(1.0),    fdouble(3.0)},
    {fdouble(-103240.467920728), fdouble(-0.000208924735834372),  fdouble(109961.532079643),   fdouble(6730.47792027321),   fdouble(0.935905881633714),  fdouble(0.935947358356231),  fdouble(-0.93589486453368),  fdouble(1.0),    fdouble(5.0)},
    {fdouble(-7.43911046752643), fdouble(-7.49680170563087),      fdouble(83544.7086195816),   fdouble(83531.3679985203),   fdouble(2.4477467508823),    fdouble(6.7666148229704),    fdouble(-2.44780121116816),  fdouble(0.6),    fdouble(2.0)},
    {fdouble(-106458.718251124), fdouble(-0.000545743677370998),  fdouble(106455.682016356),   fdouble(42.3237633727421),   fdouble(1.00559843304636),   fdouble(1.00563829185779),   fdouble(-1.95708048412661),  fdouble(0.75),   fdouble(4.6)},
    {fdouble(-106453.495071328), fdouble(-0.00399495917548577),   fdouble(106460.925622402),   fdouble(418.039392846222),   fdouble(2.49073456323941),   fdouble(2.4909137590075),    fdouble(-4.9738722919108),   fdouble(0.75),   fdouble(1.49)},
    {fdouble(-106455.157115451), fdouble(-0.00456229507833856),   fdouble(106460.145184005),   fdouble(475.327392665337),   fdouble(2.27725233310332),   fdouble(2.27744553473598),   fdouble(-4.5488537847976),   fdouble(0.75),   fdouble(1.5)},
    {fdouble(-106452.180354907), fdouble(-0.00278122155186969),   fdouble(106461.721043604),   fdouble(291.588888724572),   fdouble(3.59580873362151),   fdouble(3.59615037646864),   fdouble(-7.17813960396325),  fdouble(0.75),   fdouble(1.1)},
    {fdouble(-156436.219173519), fdouble(-13.360177523064),       fdouble(10907.4536590735),   fdouble(-0.0178811573295934),fdouble(0.0295757482829108), fdouble(0.398808602998421),  fdouble(-5.49402342863045),  fdouble(0.35),   fdouble(0.96)}
};

const fdouble fit_lt_R3[max_hs_cs+1][9] = {
/*                                                                                                                                                                                                          r_min           r_max */
    {fdouble(-106344.499357271), fdouble(-0.0870404156519758),    fdouble(106379.969770542),   fdouble(8916.02780769541),   fdouble(2.34571347967461),   fdouble(2.34558512875328),   fdouble(-4.64724093315347),  fdouble(3.0),    fdouble(6.0)},
    {fdouble(-103237.178865962), fdouble(-6.19966863330973e-05),  fdouble(109964.821133342),   fdouble(6727.38883676891),   fdouble(0.990416309150688),  fdouble(0.990415990770504),  fdouble(-0.990490798848876), fdouble(5.0),    fdouble(12.0)},
    {fdouble(-106453.321357016), fdouble(-0.0233720244005975),    fdouble(106447.424341854),   fdouble(2423.61663166259),   fdouble(1.69020647850117),   fdouble(1.69030805063035),   fdouble(-3.36829845029172),  fdouble(2.0),    fdouble(6.0)},
    {fdouble(-106457.189833221), fdouble(-0.000453936408454839),  fdouble(106457.21043453),    fdouble(42.3245989284602),   fdouble(0.499881464943715),  fdouble(0.499881437435555),  fdouble(-0.999349099198404), fdouble(4.6),    fdouble(12.0)},
    {fdouble(-106478.807529316), fdouble(-0.00470475292274558),   fdouble(106435.613140363),   fdouble(443.194839747241),   fdouble(0.192492936878364),  fdouble(0.1932366392085),    fdouble(-0.192465568481317), fdouble(1.49),   fdouble(2.0)},
    {fdouble(-106455.157115451), fdouble(-0.00456229507833856),   fdouble(106460.145184005),   fdouble(475.327392665337),   fdouble(2.27725233310332),   fdouble(2.27744553473598),   fdouble(-4.5488537847976),   fdouble(1.49),   fdouble(1.49)},
    {fdouble(-106452.180354907), fdouble(-0.00278122155186969),   fdouble(106461.721043604),   fdouble(291.588888724572),   fdouble(3.59580873362151),   fdouble(3.59615037646864),   fdouble(-7.17813960396325),  fdouble(1.1),    fdouble(1.1)},
    {fdouble(-156436.219173519), fdouble(-13.360177523064),       fdouble(10907.4536590735),   fdouble(-0.0178811573295934),fdouble(0.0295757482829108), fdouble(0.398808602998421),  fdouble(-5.49402342863045),  fdouble(0.96),   fdouble(0.96)}
};

// **************************************************************
// ********** Local functions prototypes ************************
// **************************************************************
void Set_Coulomb_Field(const fdouble phi12, fdouble E[3], const fdouble dr[3],
                   const fdouble dr2);

// Herman-Skillman (HS) potential fit functions
fdouble deriv_genericHSfit(const fdouble *par, fdouble x);
fdouble genericHSfit(const fdouble *par, fdouble x);

fdouble tmp_get_shieldr_2(const int chg_st_1, const int chg_st_2);
fdouble tmp_get_shieldr(const int chg_st, const char *message);

// **************************************************************
// ********** Accessible functions implementations **************
// **************************************************************


// **************************************************************
void Check_if_LibPotentials_is_initialized(void)
{
#ifdef YDEBUG
    if (!is_libpotentials_initialized)
    {
        std_cout << "ERROR!!!\n";
        std_cout << "is_libpotentials_initialized = " << (is_libpotentials_initialized ? "yes" : "no") << "\n";
        std_cout << "Potentials library is not initialized, please call Potentials_Initialize()\n";
        std_cout << "Exiting\n";
        abort();
    }
#endif
}

// **************************************************************
void Initialize_Simple(const fdouble &minr)
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
}

// **************************************************************
void Initialize_HS(const int &m, const fdouble &base_potential)
/**
 * Initialize super-gaussian like potential
 * @param   m               Order of the Super-Gaussian (m=1 is a simple gaussian) [-]
 * @param   base_potential  Potential depth wanted [eV]
 */
{
    hs_min_rad.resize(max_hs_cs+1); // +1 for the neutral.
    potential_paramaters potparams;

    // Find the radius where the HS potential is equal to "base_potential"
    // by doing a bisection, for all supported charge states.
    fdouble r_left, r_right, found_r, pot; // [Bohr]
    for (int cs = 0 ; cs <= max_hs_cs ; cs++)
    {
        // Initial conditions
        hs_min_rad[cs] = fit_lt_R1[cs][7];
        r_left  = fit_lt_R1[cs][7]; // Minimum radius of fit
        r_right = fdouble(10.0);             // At 10 bohr, it should be coulombic

        // Start bisection!
        // See http://en.wikipedia.org/wiki/Bisection_method#Practical_considerations
        found_r = r_right + (r_left - r_right) / libpotentials::two;

        potparams.hs_cs2 = cs;
        potparams.kQ2 = fdouble(cs) * e0;

        while (std::abs(found_r - r_left) > 1.0e-100 && std::abs(found_r - r_right) > 1.0e-100)
        {
            potparams.r = found_r * bohr_to_m;
            pot = Calculate_Potential_Cutoff_HS_SuperGaussian(NULL, NULL, potparams);
            //printf("base_potential = %10.5g   r_left = %10.5g   r = %10.5g   r_right = %10.5g   HS(r) = %10.5g\n", base_potential, r_left, found_r, r_right, pot);
            if (pot <= fdouble(std::max(1,cs))*base_potential)
            {
                r_right = found_r;
            }
            else
            {
                r_left = found_r;
            }
            found_r = r_right + (r_left - r_right) /libpotentials::two;
        }
        std_cout << "Bisection end: cs = " << cs << "  HS(r="<<found_r<<") = " << pot << "\n";

        Assert_isinf_isnan(found_r);
        Assert_isinf_isnan(pot);
        assert(found_r > 0.0);
        //assert(std::abs(pot - std::max(1,cs)*base_potential) < 1.0e-3);

        if (found_r * si_to_au_length <= 0.99999*hs_min_rad[cs])
        {
            std_cout << "##############################################\n";
            DEBUGP("Initialize_HS() called with a potential depth too deep.\n");
            std_cout << "The value found " << found_r * si_to_au_length << " Bohr (" << found_r << " m)\n";
            std_cout << "for charge state " << cs << " should not be lower than " << hs_min_rad[cs] << " Bohr (" << hs_min_rad[cs] * au_to_si_length<<" m)\n";
            std_cout << "Potential depth wanted: " << base_potential << " eV (" << base_potential*eV_to_Eh << " Eh)\n";
            std_cout << "Exiting\n";
            abort();
        }

        hs_min_rad[cs] = found_r;
    }

    Initialize_SuperGaussian(m);
}

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
void Set_Coulomb_Field(const fdouble phi12, fdouble E[3], const fdouble dr[3],
                   const fdouble r2)
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
}

// **************************************************************
fdouble deriv_genericHSfit(const fdouble *par, fdouble x){
    // The 1/2 factor is becuase HS outputs the
    // potential as 2V(x) and that's how they were fit
    return libpotentials::half*(
              par[0]*par[5]*std::pow(x,par[5]-1)/(std::pow(x,par[5])-par[1])/(std::pow(x,par[5])-par[1])
            + par[2]*par[4]/std::pow(x,par[4]+1)
            + par[6]*par[3]*std::pow(x,par[6]-1)
        );
}

// **************************************************************
// f(x) = -a/(x^n-b)-b/x^m +d*x^o where the parameters are passed alphabetically
fdouble genericHSfit(const fdouble *par, fdouble x){
    //the 1/2 factor is becuase HS outputs the potential as 2V and that's how they were fit
    return -libpotentials::half*(
            - par[0]/(std::pow(x,par[5])-par[1])
            - par[2]/std::pow(x,par[4])
            + par[3]*std::pow(x,par[6])
        );
}

// **************************************************************
fdouble LibPotentialErf(fdouble x)
{
    fdouble erf_value = 0.0;
    //erf_value = nr::int_erf(x);
    erf_value = fdouble(erf(x));
    //erf_value = nr::erff(x);
    //erf_value = nr::python_erf(x);
    //erf_value = gsl_sf_erf (x);
//     if (x < libpotentials_private::tl_Rmax)
//     {
//         const int base = int(x * libpotentials_private::tl_one_over_dR);
//         const fdouble gain = fdouble(x * libpotentials_private::tl_one_over_dR) - fdouble(base);
//         erf_value = libpotentials_private::tl_erf[base] * (1.0 - gain) + gain*libpotentials_private::tl_erf[base+1];
//     } else {
//         erf_value = 1.0;
//     }

    return erf_value;
}

// **************************************************************
// Function pointers for...
// ...setting the parameters of the potential/field calculation
void   (*Potentials_Set_Parameters)(void *p1, void *p2, potential_paramaters &potparams) = NULL;
// ...calculating the potential
fdouble (*Calculate_Potential)(      void *p1, void *p2, potential_paramaters &potparams) = NULL;
// ...setting the electric field
void   (*Set_Field)(                void *p1, void *p2, potential_paramaters &potparams, fdouble &phi, fdouble E[3]) = NULL;

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
fdouble Calculate_Potential_Cutoff_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    fdouble potential;

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
    fdouble &phi, fdouble E[3])
{
    Check_if_LibPotentials_is_initialized();

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
            potparams.B = libpotentials_private::base_pot_well_depth*fdouble(charge_state2);
            //potparams.B *= eV_to_J * (1.0 / e0); // *= 1.0
        }
        else
        {   // Negative charge (electron), negative potential
            potparams.B = -libpotentials_private::base_pot_well_depth;
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

    if (potparams.r <= potparams.cutoff_radius)
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

    if (potparams.r > potparams.cutoff_radius)
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
            potparams.B = libpotentials_private::base_pot_well_depth*fdouble(charge_state2);
            //potparams.B *= eV_to_J * (1.0 / e0); // *= 1.0
        }
        else
        {   // Negative charge (electron), negative potential
            // B is the (negative of the)
            // ionization potential of the ion.
            potparams.B = -libpotentials_private::base_pot_well_depth;
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
            potparams.sg_r_over_sigma_two_m = std::pow(potparams.r/potparams.sg_sigma, sg_two_m);
            potparams.sg_exp_half_r_over_sigma_two_m =
                    std::exp( -libpotentials::half * potparams.sg_r_over_sigma_two_m );
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
fdouble Calculate_Potential_Cutoff_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    fdouble phi12;   // Electrostatic potential

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
    fdouble &phi, fdouble E[3])
{
    Check_if_LibPotentials_is_initialized();

    if (potparams.r > potparams.cutoff_radius)
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
    Check_if_LibPotentials_is_initialized();

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);

    potparams.kQ2           = one_over_4Pieps0 * Get_Charge(p2);
    potparams.hs_cs2        = Get_Charge_State(p2);

    // If we have an electron, set the correct parameters for the
    // potential cutoff. All other cases (neutral atom or ion) are
    // treated as HS potential (neutral to 7+) and coulomb (8+ and up)
    if ( (potparams.hs_cs2 == -1) or (potparams.hs_cs2 > max_hs_cs) )
    {
        Potentials_Set_Parameters_ChargeDistribution_Symmetric(p1, p2, potparams);
    }
}

// **************************************************************
fdouble Calculate_Potential_Cutoff_HS_SuperGaussian(
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
    Check_if_LibPotentials_is_initialized();

    fdouble phi12 = 0.0;   // Electrostatic potential

    // Fits are in atomic units
    fdouble distance_au = potparams.r * si_to_au_length;

    const int cs = potparams.hs_cs2;

    // Ions are given a potential inside the electron cloud.
    // The last two fdoubles of the fit_lessthan_R array is the range
    // the fit is valid for below the smallest range we have a simple
    // hard cutoff V(r<r_0) = V(r_0)

    // ___________________________________________________________________ ...
    // |        |             |                |       |
    // |  CP    |     R1      |     R2         |  R3   |  Coulomb
    // |________|_____________|________________|_______|__________________ ...

    if (cs == -1 or cs > max_hs_cs)
    {
        // Electron or high charge state ion
        assert(p1 != NULL);
        assert(p2 != NULL);
        phi12 = Calculate_Potential_Cutoff_ChargeDistribution_Symmetric(p1, p2, potparams);
    }
    else
    {
        // Make sure that for distances less then the hard cutoff, we
        // use that cutoff;
        assert(cs < int(hs_min_rad.size()));
        if (distance_au < hs_min_rad[cs])
        {
            distance_au = hs_min_rad[cs];
        }

        if          (distance_au >= fit_lt_R3[cs][8])    // In Coulomb
        {
            if (cs == 0)
            {
                // Potential outside the electron cloud goes to 0
                // exponentially using f(x)=h*exp(-v*x+k)
                phi12 = fdouble(1.93775072943628) * std::exp(fdouble(-0.533297816151)*distance_au -fdouble(0.7486357665822807));
            }
            else
            {
                phi12 = (fdouble(cs) / distance_au);   // Outside electron cloud: Coulombic pot.
            }
        }
        else if (   (distance_au <  fit_lt_R3[cs][8]) &&
                    (distance_au >= fit_lt_R3[cs][7]))   // In R3
        {
            phi12 = genericHSfit(&(fit_lt_R3[cs][0]),distance_au);
        }
        else if (   (distance_au <  fit_lt_R2[cs][8]) &&
                    (distance_au >= fit_lt_R2[cs][7]))   // In R2
        {
            phi12 = genericHSfit(&(fit_lt_R2[cs][0]),distance_au);
        }
        else if (   (distance_au <  fit_lt_R1[cs][8]) &&
                    (distance_au >= fit_lt_R1[cs][7]))   // In R1
        {
            phi12 = genericHSfit(&(fit_lt_R1[cs][0]),distance_au);
        }
        else                                   // In CP (constant potential)
        {
            phi12 = genericHSfit(&(fit_lt_R1[cs][0]),hs_min_rad[cs]);
        }
        phi12 *= Eh_to_eV;
    }

    if (p2 != NULL)
        phi12 /= fdouble(std::max(1,Get_Charge_State(p2)));

    return phi12;
}

// **************************************************************
void Set_Field_Cutoff_HS_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3])
{
    Check_if_LibPotentials_is_initialized();

    // Fits are in atomic units
    fdouble distance_au = potparams.r * si_to_au_length;

    const int cs = potparams.hs_cs2;
    fdouble Ef, unit_dr[3];

    // Ions are given a potential inside the electron cloud.
    // The last two fdoubles of the fit_lessthan_R array is the range
    // the fit is valid for below the smallest range we have a simple
    // hard cutoff V(r<r_0) = V(r_0)

    // ___________________________________________________________________ ...
    // |        |             |                |       |
    // |  CP    |     R1      |     R2         |  R3   |  Coulomb
    // |________|_____________|________________|_______|__________________ ...

    bool field_vector_multiplication = true;

    if (cs == -1 or cs > max_hs_cs)
    {
        // Electron or high charge state ion
        assert(p1 != NULL);
        assert(p2 != NULL);
        Set_Field_Cutoff_ChargeDistribution_Symmetric(p1, p2, potparams, phi, E);
    }
    else
    {
        if (distance_au < hs_min_rad[cs]) /* In CP (constant potential) */
        {
            // No field in hard cutoff region
            Ef = 0.0;
            field_vector_multiplication = false;
        }
        else if (distance_au >= fit_lt_R3[cs][8])    /* In Coulomb */
        {
            if (cs == 0)
            {
                // Potential outside the electron cloud goes to 0
                // exponentially using f(x)=h*exp(-v*x+k)
                Ef = fdouble(-0.272) * (
                             (fdouble(-1.93775072943628))
                        *    (fdouble(-0.533297816151))
                        *    (fdouble(-0.533297816151)*distance_au - fdouble(0.7486357665822807))
                        * std::exp(fdouble(-0.533297816151)*distance_au - fdouble(0.7486357665822807))
                    );
            }
            else
            {
                Ef = fdouble(cs) / (distance_au*distance_au);
            }
        }
        else if (   (distance_au <  fit_lt_R3[cs][8]) &&
                    (distance_au >= fit_lt_R3[cs][7]))   /* In R3 */
        {
            Ef = deriv_genericHSfit(&(fit_lt_R3[cs][0]),distance_au);
        }
        else if (   (distance_au <  fit_lt_R2[cs][8]) &&
                    (distance_au >= fit_lt_R2[cs][7]))   /* In R2 */
        {
            Ef = deriv_genericHSfit(&(fit_lt_R2[cs][0]),distance_au);
        }
        else if (   (distance_au <  fit_lt_R1[cs][8]) &&
                    (distance_au >= fit_lt_R1[cs][7]))   /* In R1 */
        {
            Ef = deriv_genericHSfit(&(fit_lt_R1[cs][0]),distance_au);
        }
        else                                   /* In CP (constant potential) */
        {
            printf("CAN'T BE HERE!!!\n");
            abort();
        }

        if (p2 != NULL)
            Ef /= fdouble(std::max(1,Get_Charge_State(p2)));

        Ef *= au_to_si_field;
        if (field_vector_multiplication)
        {
            for (int d = 0 ; d < 3 ; d++)
            {
                unit_dr[d] = potparams.dr[d] * potparams.one_over_r;
                E[d]  += unit_dr[d] * Ef;
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
        potparams.cutoff_radius = 0.0;
        potparams.kQ2           = 0.0;
        potparams.kQ2_over_B    = 0.0;
        potparams.B             = 0.0;
        potparams.gd_sigma      = 0.0;
    }
    else if( charge_state2 != 0 )
    {
      //fdouble Ip =element.IpsLowest[std::abs(charge_state2)];
      fdouble Ip = libpotentials_private::base_pot_well_depth;
      //charges are not equal
      //take the higher charge as the distribution
      //Thus set the parameters for the guassian
      //using the values of the larger charge
      //NOTE: it will always do this for p2=electron
      // unless p1 is also an electron
        if ((charge_state2 < charge_state1) && (charge_state1 != 0)){
          Ip = libpotentials_private::base_pot_well_depth;
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
            potparams.B = -libpotentials_private::base_pot_well_depth;
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
/*        std_cout << "Id(p1)="<<Get_Id(p1)<<"  B="<<potparams.B<<"  Cs(p1)="<<Get_Charge_State(p1)<<" Cs(p2)="<<Get_Charge_State(p2)<<" kQ2_over_B="<<potparams.kQ2_over_B<<" well="<<libpotentials_private::base_pot_well_depth<<"\n";*/
        // Radius where the Coulomb potential and its first derivative are
        // equal to a the gaussian charge distribution potential.
        potparams.cutoff_radius = libpotentials::eight* potparams.gd_sigma;
//         potparams.cutoff_radius = 4.0 * potparams.gd_sigma;
//         potparams.cutoff_radius = 2.5 * potparams.gd_sigma;
        //potparams.cutoff_radius = libpotentials::two * potparams.gd_sigma;

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
fdouble Calculate_Potential_Cutoff_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    fdouble phi12=0.0;   // Electrostatic potential

    if      (potparams.r > potparams.cutoff_radius)
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

    if      (potparams.r > potparams.cutoff_radius)
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

// **************************************************************
// **************************************************************

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
        potparams.gd_sigma = one_over_4Pieps0 * e0 / (libpotentials_private::base_pot_well_depth * sqrt_Pi);
        potparams.cutoff_radius = libpotentials::eight* potparams.gd_sigma;
    }
    else
    {
        potparams.kQ2 = 0.0;
        potparams.gd_sigma = 0.0;
        potparams.cutoff_radius = 0.0;
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
        if (potparams.r > potparams.cutoff_radius)
        {
            // If the distance is not less then the shielding radius, get the
            // normal Coulomb potential.

            potential = Coulomb_Potential(potparams.kQ2, potparams.r);
        }
        else
        {
            // Else, use the lookup table for erf(x)/x
            // Note that y = x / sqrt(2) where x is the unitless distance from GaussianDistribution
            const fdouble y = potparams.r / ( two * potparams.gd_sigma );
            Assert_isinf_isnan(y);
            potential = potparams.kQ2 / (two * potparams.gd_sigma) * libpotentials_private::lut_potential.read(y);

            //printf("kQ2 = %20.15g  sigma = %20.15g  potential = %20.15g  E = %20.15g   x=%20.15g\n",
            //       potparams.kQ2, potparams.gd_sigma*si_to_au_length, potential*si_to_au_pot, 0*si_to_au_field, y);
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
        if (potparams.r > potparams.cutoff_radius)
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
            const fdouble y = potparams.r / (two * potparams.gd_sigma);
            const fdouble E_over_r = potparams.kQ2 / (four * potparams.gd_sigma * potparams.gd_sigma * potparams.gd_sigma) * libpotentials_private::lut_field.read(y);

            Assert_isinf_isnan(E_over_r);

            for (int d = 0 ; d < 3 ; d++)
            {
                // We multiply E/r by r to get E
                E[d]  += E_over_r * potparams.dr[d];
            }

            //printf("kQ2 = %20.15g  sigma = %20.15g  potential = %20.15g  E = %20.15g   x=%20.15g\n",
            //       potparams.kQ2, potparams.gd_sigma*si_to_au_length, 0, E[0]*si_to_au_field, y);
        }
    }
}

// **************************************************************
void Potentials_Set_Parameters_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
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

// **************************************************************
void Potentials_Set_Parameters_PseudoParticles(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);
    potparams.kQ2 = one_over_4Pieps0 * Get_Charge(p2);
}

// **************************************************************
fdouble Calculate_Potential_PseudoParticles(
    void *p1, void *p2,
    potential_paramaters &potparams)
{
    Check_if_LibPotentials_is_initialized();

    fdouble phi = - potparams.kQ2 * (
          ( ps_A2 / ( libpotentials::two * ps_A_minus_B * std::pow(ps_A + potparams.r, 2) ) )
        - ( ps_A_A_minus2B / ( ps_A_minus_B2 * (ps_A + potparams.r) ) )
        - ( ps_B2 * std::log( (ps_A+potparams.r)/(ps_B+potparams.r) ) / ps_A_minus_B3 )
        - ( ps_C / ( libpotentials::two * ( std::pow(potparams.r, 2) + ps_D) ) )
    );
    return phi;
}

// **************************************************************
void Set_Field_PseudoParticles(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3])
{
    Check_if_LibPotentials_is_initialized();

    const fdouble Eabs = potparams.kQ2 * (
          ( potparams.r / ( std::pow(potparams.r + ps_A, 3) * (potparams.r + ps_B) ) )
        + ( ps_C / std::pow( std::pow(potparams.r, 2) + ps_D, 2) )
    );
    for (int d = 0 ; d < 3 ; d++)
        E[d] += potparams.dr[d] * Eabs;
}

// **************************************************************
fdouble tmp_get_shieldr(const int chg_st, const char *message)
{
    if( chg_st == 0 )
    {
        std_cout << "error in gecutoff_radius(): chg_st == 0 from="<<message<<"\n";
        abort();
    }

    return (s_cutoff_radius_from_input_file + libpotentials::half * fdouble(std::abs(chg_st)-1)*one_over_seven * a0);
//    return(shieldr);
}

// **************************************************************
fdouble tmp_get_shieldr_2(const int chg_st_1, const int chg_st_2)
{
    if( chg_st_1 == 0 )
    {
        std_cout << "error in gecutoff_radius(): chg_st == 0 from=gecutoff_radius_2\n";
        abort();
    }

    int abs_chg_st_1 = std::abs(chg_st_1);
    int abs_chg_st_2 = std::abs(chg_st_2);
    int max_chg_st = (abs_chg_st_1 > abs_chg_st_2 ) ? abs_chg_st_1 : abs_chg_st_2;

    return tmp_get_shieldr(max_chg_st,"gecutoff_radius_2");
}

// ********** End of file ***************************************
