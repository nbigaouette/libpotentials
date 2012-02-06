
#include <StdCout.hpp>

#include "LibPotentials.hpp"
#include "Constants.hpp"
#include "Code_Functions_Declarations.hpp"
#include "Git_Diff.hpp"
#include "Version.hpp"

extern void Initialize_SuperGaussian(const int &m);
extern void Initialize_HS(const fdouble &base_potential);

std::string io_basename;

bool is_libpotentials_initialized = false;

namespace libpotentials_private
{
    fdouble cutoff_base_potential;
    fdouble cutoff_radius;
    bool cutoff_radius_linear_in_CS;

//     // Error function lookup table erf(R)
//     fdouble *tl_erf;
//     const int    tl_n = 1000;    // Number of points of the lookup table
//     const fdouble tl_Rmax = 6.0; // Maximum value of R: erf(4) = 0.999999984582742
//     const fdouble tl_dR = tl_Rmax / (tl_n-1); // Step
//     const fdouble tl_one_over_dR = 1.0 / tl_dR;
//
//     // **********************************************************
//     void initialize_erf_lookup_table()
//     {
//         tl_erf = (fdouble *) calloc_and_check(tl_n, sizeof(fdouble));
//         // Populating the lookup table
//         for (int i = 0 ; i < tl_n ; i++)
//         {
//             tl_erf[i] = nr::int_erf(i*tl_dR);
//         }
//     }

    LookUpTable<fdouble> lut_potential;
    LookUpTable<fdouble> lut_field;
}

// **************************************************************
// ********** Function pointers for... **************************
// ...setting the parameters of the potential/field calculation
void   (*Potentials_Set_Parameters)(void *p1, void *p2, potential_paramaters &potparams);
// ...calculating the potential
fdouble (*Calculate_Potential)(      void *p1, void *p2, potential_paramaters &potparams);
// ...setting the electric field
void   (*Set_Field)(                void *p1, void *p2, potential_paramaters &potparams, fdouble &phi, fdouble E[3]);

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
void Potentials_Initialize(const std::string _io_basename,
                           const std::string potential_shape,
                           const fdouble cutoff_base_potential,
                           const fdouble cutoff_radius,
                           const int input_sg_m)
/**
 * Initialize potentials library.
 * @param   potential_shape         String containing the type of close
 *                                  range potential to use
 * @param   cutoff_base_potential   Potential depth of a 1+ ion [eV] (== -1 if cutoff radius is wanted)
 * @param   cutoff_radius           Radius at which Coulomb/HS potential becomes the smoothed value (== -1 if a base potential is wanted) [m]
 * @param   input_sg_m              Super-Gaussian "m" parameter
 */
{
    is_libpotentials_initialized = true;

    io_basename = _io_basename;

    if (cutoff_base_potential <= 0.0 and cutoff_radius <= 0.0)
    {
        std_cout
            << "Error in processing cutoff_base_potential=" << cutoff_base_potential << " or cutoff_radius=" << cutoff_radius << "\n"
            << "Only one of these two should be negative (negative value enables the other one).\n"
            << "Aborting\n";
        std_cout.Flush();
        abort();
    }
    else if (cutoff_base_potential > 0.0 and cutoff_radius > 0.0)
    {
        std_cout
            << "Error in processing cutoff_base_potential=" << cutoff_base_potential << " or cutoff_radius=" << cutoff_radius << "\n"
            << "Only one of these two should be negative (negative value enables the other one).\n"
            << "Aborting\n";
        std_cout.Flush();
        abort();
    }

    libpotentials_private::cutoff_base_potential = cutoff_base_potential;
    libpotentials_private::cutoff_radius         = cutoff_radius;

    Potentials_Set_Parameters = NULL;
    Calculate_Potential       = NULL;
    Set_Field                 = NULL;

    std_cout
        << "######################################################################\n"
        << "###             Potentials library initialization                  ###\n"
        << "###----------------------------------------------------------------###\n";

    if (potential_shape == "PureCoulomb")
    {
        std_cout << "### Using a pure Coulomb interaction                               ###\n";

        libpotentials_private::cutoff_base_potential = 0.0;
        libpotentials_private::cutoff_radius         = 0.0;

        Potentials_Set_Parameters = &Potentials_Set_Parameters_Simple;
        Calculate_Potential       = &Calculate_Potential_Cutoff_Simple;
        Set_Field                 = &Set_Field_Cutoff_Simple;
    }
    else if (potential_shape == "Simple")
    {
        std_cout << "### Using a simple cutoff                                          ###\n";
        std_cout << "### for close range interaction                                    ###\n";

        Initialize_Simple(cutoff_base_potential, cutoff_radius);

        std_cout << "### Cutoff radius  = " << libpotentials_private::cutoff_radius*libpotentials::m_to_bohr << " bohr\n";
        std_cout << "### Base potential = " << libpotentials_private::cutoff_base_potential << " eV";
        std_cout <<                   " = " << libpotentials_private::cutoff_base_potential*libpotentials::eV_to_Eh << " Eh\n";

        Potentials_Set_Parameters = &Potentials_Set_Parameters_Simple;
        Calculate_Potential       = &Calculate_Potential_Cutoff_Simple;
        Set_Field                 = &Set_Field_Cutoff_Simple;
    }
    else if (potential_shape == "Harmonic")
    {
        std_cout << "### Using an harmonic potential                                    ###\n";
        std_cout << "### for close range interaction                                    ###\n";

        Initialize_Harmonic(cutoff_base_potential, cutoff_radius);

        Potentials_Set_Parameters = &Potentials_Set_Parameters_Harmonic;
        Calculate_Potential       = &Calculate_Potential_Cutoff_Harmonic;
        Set_Field                 = &Set_Field_Cutoff_Harmonic;
    }
    else if (potential_shape == "SuperGaussian")
    {
        std_cout << "### Using a super-gaussian potential                               ###\n";
        std_cout << "### for close range interaction                                    ###\n";
        std_cout << "### with m = " << input_sg_m << "\n";

        Initialize_SuperGaussian(input_sg_m, cutoff_base_potential, cutoff_radius);

        Potentials_Set_Parameters = &Potentials_Set_Parameters_SuperGaussian;
        Calculate_Potential       = &Calculate_Potential_Cutoff_SuperGaussian;
        Set_Field                 = &Set_Field_Cutoff_SuperGaussian;
    }
    else if (potential_shape == "GaussianDistribution")
    {
        std_cout << "### Using a gaussian charge distribution                           ###\n";
        std_cout << "### potential for close range interaction                          ###\n";
        std_cout << "### Initializing the lookup tables...                              ###\n" << std::flush;
        libpotentials_private::lut_potential.Initialize(erf_over_x,                 0.0, fdouble(4.5*std::sqrt(2.0)), 10000, "Potential LookUpTable");
        libpotentials_private::lut_field.Initialize(erf_over_x3_minus_exp_over_x2,  0.0, fdouble(4.5*std::sqrt(2.0)), 10000, "Field LookUpTable");
        std_cout << "### Initializing the lookup tables done.                           ###\n" << std::flush;

        Initialize_GaussianDistribution(cutoff_base_potential, cutoff_radius);

        Potentials_Set_Parameters = &Potentials_Set_Parameters_GaussianDistribution;
        Calculate_Potential       = &Calculate_Potential_Cutoff_GaussianDistribution;
        Set_Field                 = &Set_Field_Cutoff_GaussianDistribution;
    }
    else if (potential_shape == "HermanSkillman")
    {
        std_cout << "### Using the Herman-Skillman (HS) potential                        ##\n";
        std_cout << "### for close range interaction                                    ###\n";

        Initialize_HermanSkillman(cutoff_base_potential, cutoff_radius);

        // Same as symmetric. Necessary for ion-ion interactions
        libpotentials_private::lut_potential.Initialize(erf_over_x,                 0.0, fdouble(4.5*std::sqrt(2.0)), 10000, "Potential LookUpTable");
        libpotentials_private::lut_field.Initialize(erf_over_x3_minus_exp_over_x2,  0.0, fdouble(4.5*std::sqrt(2.0)), 10000, "Field LookUpTable");

        Potentials_Set_Parameters = &Potentials_Set_Parameters_HS;
        Calculate_Potential       = &Calculate_Potential_Cutoff_HS;
        Set_Field                 = &Set_Field_Cutoff_HS;
    }
    else if (potential_shape == "Symmetric")
    {
        std_cout << "### Using the symmetric two charge                                 ###\n";
        std_cout << "### distribution interaction                                       ###\n";
        std_cout << "### Initializing the lookup tables...                              ###\n" << std::flush;
        libpotentials_private::lut_potential.Initialize(erf_over_x,                 0.0, fdouble(4.5*std::sqrt(2.0)), 10000, "Potential LookUpTable");
        libpotentials_private::lut_field.Initialize(erf_over_x3_minus_exp_over_x2,  0.0, fdouble(4.5*std::sqrt(2.0)), 10000, "Field LookUpTable");
        std_cout << "### Initializing the lookup tables done.                           ###\n" << std::flush;
        std_cout <<                                            "Done!...                ###\n";

        Initialize_Symmetric(cutoff_base_potential, cutoff_radius);

        Potentials_Set_Parameters = &Potentials_Set_Parameters_ChargeDistribution_Symmetric;
        Calculate_Potential       = &Calculate_Potential_Cutoff_ChargeDistribution_Symmetric;
        Set_Field                 = &Set_Field_Cutoff_ChargeDistribution_Symmetric;
    }
    else if (potential_shape == "ScreenedCoulomb")
    {
        std_cout << "### Using the Coulomb potential screened with                       ##\n";
        std_cout << "### parameter alpha = " << sc_alpha * libpotentials::m_to_angstrom << "                                        ###\n";

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

    Log_Git_Info(io_basename);

    std_cout << "###----------------------------------------------------------------###\n"
             << "###            Potentials library initialization done.             ###\n"
             << "######################################################################\n";
}

// **************************************************************
void Potentials_Finalize()
{
}

// **************************************************************
void Print_Particles(void *list, const int &N)
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
    for (int j = 0 ; j < N ; j++)
    {
        pv = get_voidp(list, j);
        std_cout << "a"<<j<<" Id." << Get_Id(pv);
        std_cout << " ("<<pv<<")";
        std_cout << " ("<<Get_Position(pv)[0] * libpotentials::m_to_bohr<<","<<Get_Position(pv)[1] * libpotentials::m_to_bohr<<","<<Get_Position(pv)[2] * libpotentials::m_to_bohr<<")";
        std_cout << " ("<<Get_Velocity(pv)[0] * libpotentials::si_to_au_vel<<","<<Get_Velocity(pv)[1] * libpotentials::si_to_au_vel<<","<<Get_Velocity(pv)[2] * libpotentials::si_to_au_vel<<")";
        std_cout << " "<<Get_Charge_State(pv)<<" ";
        std_cout << " " << 0.5 * Get_Mass(pv) * (
            Get_Velocity(pv)[0]*Get_Velocity(pv)[0] + Get_Velocity(pv)[1]*Get_Velocity(pv)[1] + Get_Velocity(pv)[2]*Get_Velocity(pv)[2]
        ) * libpotentials::J_to_eV;
        std_cout << " " << fdouble(Get_Charge_State(pv)) * Get_Potential(pv);
        //std_cout << " " << Get_Mass(pv) * (1000.0 * Na);
        std_cout << " " << Get_Mass(pv) * libpotentials::si_to_au_mass;
        std_cout << " ("<<Get_E(pv)[0]<<","<<Get_E(pv)[1]<<","<<Get_E(pv)[2]<<")";
        std_cout << " " << Get_Potential(pv);
        std_cout << "\n";
    }
    std_cout << "--------\n\n";
}
fdouble erf_over_x(fdouble x)
/**
 * Calculate electrostatic potential:
 *      V(x) / (k*Q / (sqrt(2)*sigma))
 * where x = r/(sqrt(2)*sigma)
 *
 * Used to construct the lookup table.
 *
 * See doc/expansions/expansions.pdf and scripts/expansions.py
 *
 * http://www.wolframalpha.com/input/?i=erf%28x%29%2Fx
 */
{
    fdouble value;

    if (x < libpotentials::one)
    {
        value = (
               libpotentials::two
            - (libpotentials::two * x*x)                                                                / fdouble(3.0)
            +        x*x*x*x                                                                            / fdouble(5.0)
            -        x*x*x*x*x*x                                                                        / fdouble(21.0)
            +        x*x*x*x*x*x*x*x                                                                    / fdouble(108.0)
            -        x*x*x*x*x*x*x*x*x*x                                                                / fdouble(660.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x                                                            / fdouble(4680.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                        / fdouble(37800.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                    / fdouble(342720.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                / fdouble(3447360.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                            / fdouble(38102400.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                        / fdouble(459043200.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                    / fdouble(5987520000.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                / fdouble(84064780800.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                            / fdouble(1264085222400.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                        / fdouble(20268952704000.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                    / fdouble(345226033152000.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                / fdouble(6224529991680000.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x            / fdouble(118443913555968000.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x        / fdouble(2372079457972224000.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x    / fdouble(49874491167621120000.0)
            //   +O(x^41)
        );
        value /= libpotentials::sqrt_Pi;
    }
    else
    {
        value = fdouble(erf(x)) / x;
    }

    return value;
}

// **************************************************************
fdouble erf_over_x3_minus_exp_over_x2(fdouble x)
/**
 * Calculate electrostatic field:
 *      (E(x) / x) / (k*Q / (sqrt(2)*sigma^3))
 * where x = r/(sqrt(2)*sigma)
 *
 * Used to construct the lookup table.
 *
 * See doc/expansions/expansions.pdf and scripts/expansions.py
 *
 * http://www.wolframalpha.com/input/?i=expansion+erf%28x%29%2F%282*x**3%29-1%2Fsqrt%28pi%29*exp%28-x**2%29%2F%28x**2%29
 */
{
    fdouble value;

    if (x < 1.0)
    {
        // http://www.wolframalpha.com/input/?i=expansion+erf%28x%29%2F%282*x**3%29-1%2Fsqrt%28pi%29*exp%28-x**2%29%2F%28x**2%29
        value = (
               libpotentials::two                                                                               / fdouble(3.0)
            - (libpotentials::two * x*x)                                                                        / fdouble(5.0)
            +        x*x*x*x                                                                                    / fdouble(7.0)
            -        x*x*x*x*x*x                                                                                / fdouble(27.0)
            +        x*x*x*x*x*x*x*x                                                                            / fdouble(132.0)
            -        x*x*x*x*x*x*x*x*x*x                                                                        / fdouble(780.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x                                                                    / fdouble(5400.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                                / fdouble(42840.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                            / fdouble(383040.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                        / fdouble(3810240.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                    / fdouble(41731200.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                / fdouble(498960000.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                            / fdouble(6466521600.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                        / fdouble(90291801600.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                    / fdouble(1351263513600.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                / fdouble(21576627072000.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                            / fdouble(366148823040000.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                        / fdouble(6580217419776000.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                    / fdouble(124846287261696000.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                / fdouble(2493724558381056000.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x            / fdouble(52307393175797760000.0)
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x        / fdouble(1149546198863462400000.0)
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x    / fdouble(26414017102773780480000.0)
            // +O(x^45)
        ) / libpotentials::sqrt_Pi;
    }
    else
    {
        value = fdouble(erf(x))/(libpotentials::two*x*x*x) - libpotentials::one/libpotentials::sqrt_Pi*std::exp(-x*x)/(x*x);
    }

    return value;
}


// ********** End of file ***************************************
