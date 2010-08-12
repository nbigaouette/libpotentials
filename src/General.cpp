/****************************************************************

    General functions

****************************************************************/

#ifdef __SUNPRO_CC
#include <string.h>
#else
#include <cstring>
#endif

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <cfloat>
#include <climits> // CHAR_BIT
#include <stdint.h> // (u)int64_t

#include "Constants.hpp"
#include "General.hpp"
#include "Potentials.hpp"
#include "NR.hpp"
#include "Global.hpp"
#include "Code_Functions_Declarations.hpp"
#include "Std_Cout.hpp"
#include "Version.hpp"


using namespace libpotentials;

// Is the library initialized?
int is_initialized = -1;

// Error function lookup table
// erf(R)
double *tl_erf;             // Table
const int    tl_n = 200;    // Number of points of the lookup table
const double tl_Rmax = 4.0; // Maximum value of R: erf(4) = 0.999999984582742
const double tl_dR = tl_Rmax / (tl_n-1); // Step
const double tl_one_over_dR = 1.0 / tl_dR;
//
void initialize_erf_lookup_table();

// **************************************************************
// ********** Accessible functions implementations **************
// **************************************************************

void Check_if_initialized(void)
{
#ifdef YDEBUG
    if (!is_initialized)
    {
        std_cout << "ERROR!!!\n";
        std_cout << "is_initialized = " << (is_initialized ? "yes" : "no") << "\n";
        std_cout << "Potentials library is not initialized, please call Potentials_Initialize()\n";
        std_cout << "Exiting\n";
        abort();
    }
#endif
}

// Compare two numbers. If their relative difference is
// bigger than a certain precision (1e-14) then assert
// (exit with error)
void assert_diff(double a, double b)
{
    if ((a+b) > DBL_MIN)
        assert((fabs( (a-b) / (a+b) )) <= 1.0e-14);
}

// **************************************************************
void Potentials_Initialize(const std::string potential_shape,
                           const double base_potential_depth,
                           const double input_s_rmin,
                           const double input_sg_m
                          )
{
    is_initialized = true;

    base_pot_well_depth = base_potential_depth; //in eV

    Potentials_Set_Parameters = NULL;
    Calculate_Potential       = NULL;
    Set_Field                 = NULL;

    std_cout << "###############################################\n";
    std_cout << "### Potentials library initialization...    ###\n";
    std_cout << "###############################################\n";

    if (potential_shape == "PureCoulomb")
    {
        double r_temp = 0.0;
        std_cout << "### Using a pure Coulomb interaction        ###\n";
        Initialize_Simple(r_temp);

        Potentials_Set_Parameters = &Potentials_Set_Parameters_Simple;
        Calculate_Potential       = &Calculate_Potential_Cutoff_Simple;
        Set_Field                 = &Set_Field_Cutoff_Simple;
    }
    else if (potential_shape == "Simple")
    {
        std_cout << "### Using a simple cutoff                   ###\n";
        std_cout << "### for close range interaction             ###\n";
        Initialize_Simple(input_s_rmin);

        Potentials_Set_Parameters = &Potentials_Set_Parameters_Simple;
        Calculate_Potential       = &Calculate_Potential_Cutoff_Simple;
        Set_Field                 = &Set_Field_Cutoff_Simple;
    }
    else if (potential_shape == "Harmonic")
    {
        std_cout << "### Using an harmonic potential             ###\n";
        std_cout << "### for close range interaction             ###\n";

        Potentials_Set_Parameters = &Potentials_Set_Parameters_Harmonic;
        Calculate_Potential       = &Calculate_Potential_Cutoff_Harmonic;
        Set_Field                 = &Set_Field_Cutoff_Harmonic;
    }
    else if (potential_shape == "SuperGaussian")
    {
        std_cout << "### Using a super-gaussian potential        ###\n";
        std_cout << "### for close range interaction             ###\n";
        std_cout << "### with m = " << input_sg_m << "\n";
        Initialize_SuperGaussian(input_sg_m);

        Potentials_Set_Parameters = &Potentials_Set_Parameters_SuperGaussian;
        Calculate_Potential       = &Calculate_Potential_Cutoff_SuperGaussian;
        Set_Field                 = &Set_Field_Cutoff_SuperGaussian;
    }
    else if (potential_shape == "GaussianDistribution")
    {
        std_cout << "### Using a gaussian charge distribution    ###\n";
        std_cout << "### potential for close range interaction   ###\n";
        initialize_erf_lookup_table();

        Potentials_Set_Parameters = &Potentials_Set_Parameters_GaussianDistribution;
        Calculate_Potential       = &Calculate_Potential_Cutoff_GaussianDistribution;
        Set_Field                 = &Set_Field_Cutoff_GaussianDistribution;
    }
    else if (potential_shape == "HermanSkillman")
    {
        std_cout << "### Using the Herman-Skillman (HS) potential ##\n";
        std_cout << "### for close range interaction             ###\n";
        std_cout << "### and Super-Gaussian for electrons and 8+ ###\n";
        std_cout << "### and up ions (m = " << input_sg_m << ")                    ###\n";
        Initialize_HS(input_sg_m, input_s_rmin);

        Potentials_Set_Parameters = &Potentials_Set_Parameters_HS_SuperGaussian;
        Calculate_Potential       = &Calculate_Potential_Cutoff_Harmonic;
        Set_Field                 = &Set_Field_Cutoff_HS_SuperGaussian;
    }
    else if (potential_shape == "Symmetric")
    {
        std_cout << "### Using the symmetric two charge          ###\n";
        std_cout << "### distribution interaction                ###\n";
        initialize_erf_lookup_table();

        Potentials_Set_Parameters = &Potentials_Set_Parameters_ChargeDistribution_Symmetric;
        Calculate_Potential       = &Calculate_Potential_Cutoff_ChargeDistribution_Symmetric;
        Set_Field                 = &Set_Field_Cutoff_ChargeDistribution_Symmetric;
    }
    else if (potential_shape == "ScreenedCoulomb")
    {
        std_cout << "### Using the Coulomb potential screened with##\n";
        std_cout << "### parameter alpha = " << sc_alpha * libpotentials::m_to_angstrom << "                 ###\n";

        Potentials_Set_Parameters = &Potentials_Set_Parameters_ScreenedCoulomb;
        Calculate_Potential       = &Calculate_Potential_Cutoff_ScreenedCoulomb;
        Set_Field                 = &Set_Field_Cutoff_ScreenedCoulomb;
    }
    else
    {
        std_cout << "### ERROR\n";
        std_cout << "    ### Potentials_Initialize() called without a:\n";
        std_cout << "    ### valid potential shape (" << potential_shape << ")\n";
        std_cout << "    ### Exiting\n";
        abort();
    }

    std_cout << "### Git versioning:"   << std::endl
             << "###     build_time:   "<< build_time << std::endl
             << "###     build_sha:    "<< build_sha << std::endl
             << "###     build_branch: "<< build_branch << std::endl;
    std_cout << "### Ionization library initialization done. ###\n";
    std_cout << "###############################################\n";
}

// **************************************************************
void Potentials_Finalize()
{
    free(tl_erf);
}

// **************************************************************
const char* Potentials_Get_Git_Commit()
{
    return build_sha;
}

// **************************************************************
void initialize_erf_lookup_table()
{
    tl_erf = (double *) calloc_and_check(tl_n, sizeof(double));
    if (tl_erf == NULL)
    {
        std_cout << "Allocationg for tl_erf failed!\nAborting.\n";
        abort();
    }
    // Populating the lookup table
    for (int i = 0 ; i < tl_n ; i++)
    {
        tl_erf[i] = nr::int_erf(i*tl_dR);
    }
}

// **************************************************************
int factorial(int x)
/**
 * Calculate the factorial of a number.
 */
{
    if ( x == 1 || x == 0 ) return 1;
    else if (x < 0)
    {
        std_cout << "ERROR:\nCall of function \"int factorial(int x)\" with negative value.\n";
    }

    int res = 1;
    for (int i = 2 ; i <= x ; i++) res *= i;

    return res;
}

// **************************************************************
double factoriald(int x)
{
    return double(factorial(x));
}

// **************************************************************
double my_erf(double x)
{
    double erf_value = 0.0;
    //erf_value = nr::int_erf(x);
    //erf_value = erf(x);
    //erf_value = nr::erff(x);
    //erf_value = nr::python_erf(x);
    //erf_value = gsl_sf_erf (x);
    if (x < tl_Rmax)
    {
        const int base = int(x * tl_one_over_dR);
        const double gain = double(x * tl_one_over_dR) - double(base);
        erf_value = tl_erf[base] * (1.0 - gain) + gain*tl_erf[base+1];
    } else {
        erf_value = 1.0;
    }

    return erf_value;
}

// **************************************************************
void Print_Partiles(void *al, void *el,
                    const int &Nb_atoms, const int &Nb_atoms_max,
                    const int &Nb_electrons, const int &Nb_electrons_max)
{
    void *pv;

    std_cout << "     Id ";
    std_cout << " Address";
    std_cout << "    (pos bohr)";
    std_cout << "                         (vel a.u.)";
    std_cout << "                         CS";
    std_cout << " Type  K [eV]    U [eV] ";
    std_cout << " Mass [au]";
    std_cout << " E [V/m]";
    std_cout << "                  Potential [V]";
    std_cout << "\n";
    if (Nb_atoms_max > 0)
        std_cout << "Atoms: (" << Nb_atoms << "/" << Nb_atoms_max << ")\n";
    for (int j = 0 ; j < Nb_atoms ; j++)
    {
        pv = get_voidp(al, j);
        std_cout << "a"<<j<<" Id." << Get_Id(pv);
        std_cout << " ("<<pv<<")";
        std_cout << " ("<<Get_Position(pv)[0] * m_to_bohr<<","<<Get_Position(pv)[1] * m_to_bohr<<","<<Get_Position(pv)[2] * m_to_bohr<<")";
        std_cout << " ("<<Get_Velocity(pv)[0] * si_to_au_vel<<","<<Get_Velocity(pv)[1] * si_to_au_vel<<","<<Get_Velocity(pv)[2] * si_to_au_vel<<")";
        std_cout << " "<<Get_Charge_State(pv)<<" ";
        std_cout << " "<<Get_Type(pv)<<" ";
        std_cout << " " << 0.5 * Get_Mass(pv) * (
            Get_Velocity(pv)[0]*Get_Velocity(pv)[0] + Get_Velocity(pv)[1]*Get_Velocity(pv)[1] + Get_Velocity(pv)[2]*Get_Velocity(pv)[2]
        ) * J_to_eV;
        std_cout << " " << Get_Charge_State(pv) * Get_Potential(pv);
        //std_cout << " " << Get_Mass(pv) * (1000.0 * Na);
        std_cout << " " << Get_Mass(pv) * si_to_au_mass;
        std_cout << " ("<<Get_E(pv)[0]<<","<<Get_E(pv)[1]<<","<<Get_E(pv)[2]<<")";
        std_cout << " " << Get_Potential(pv);
        std_cout << "\n";
    }
    if (Nb_electrons > 0)
    {
    std_cout << "Electrons: (" << Nb_electrons << "/" << Nb_electrons_max << ")\n";
//     for (int j = 0 ; j < Nb_electrons_max ; j++)
    for (int j = 0 ; j < Nb_electrons ; j++)
    {
        pv = get_voidp(el, j);
        std_cout << "e"<<j<<" Id." << Get_Id(pv);
        std_cout << " ("<<pv<<")";
        std_cout << " ("<<Get_Position(pv)[0] * m_to_bohr<<","<<Get_Position(pv)[1] * m_to_bohr<<","<<Get_Position(pv)[2] * m_to_bohr<<")";
        std_cout << " ("<<Get_Velocity(pv)[0] * si_to_au_vel<<","<<Get_Velocity(pv)[1] * si_to_au_vel<<","<<Get_Velocity(pv)[2] * si_to_au_vel<<")";
        std_cout << " "<<Get_Charge_State(pv)<<" ";
        std_cout << " "<<Get_Type(pv)<<" ";
        std_cout << " " << 0.5 * Get_Mass(pv) * (
            Get_Velocity(pv)[0]*Get_Velocity(pv)[0] + Get_Velocity(pv)[1]*Get_Velocity(pv)[1] + Get_Velocity(pv)[2]*Get_Velocity(pv)[2]
        ) * J_to_eV;
        std_cout << " " << Get_Charge_State(pv) * Get_Potential(pv);
        //std_cout << " " << Get_Mass(pv) * (1000.0 * Na);
        std_cout << " " << Get_Mass(pv) * si_to_au_mass;
        std_cout << " ("<<Get_E(pv)[0]<<","<<Get_E(pv)[1]<<","<<Get_E(pv)[2]<<")";
        std_cout << " " << Get_Potential(pv);
        std_cout << "\n";
    }
    }
    std_cout << "--------\n\n";
}

// ********** End of file ***************************************

