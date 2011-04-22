
#include <StdCout.hpp>

#include "LibPotentials.hpp"
#include "Constants.hpp"
#include "Potentials.hpp"
#include "Global.hpp"
#include "Code_Functions_Declarations.hpp"

extern void Initialize_Simple(const fdouble &minr);
extern void Initialize_SuperGaussian(const int &m);
extern void Initialize_HS(const int &input_sg_m, const fdouble &min_rad);

// **************************************************************
void Potentials_Initialize(const std::string potential_shape,
                           const fdouble base_potential_depth,
                           const fdouble input_s_rmin,
                           const int input_sg_m)
/**
 * Initialize potentials library.
 * @param   potential_shape         String containing the type of close
 *                                  range potential to use
 * @param   base_potential_depth    Potential depth of a 1+ ion [eV]
 * @param   input_s_rmin            Simple cutoff distance [m]
 * @param   input_sg_m              Super-Gaussian "m" parameter
 */
{
    is_libpotentials_initialized = true;

    libpotentials_private::base_pot_well_depth = base_potential_depth;

    Potentials_Set_Parameters = NULL;
    Calculate_Potential       = NULL;
    Set_Field                 = NULL;

    std_cout
        << "######################################################################\n"
        << "###             Potentials library initialization                  ###\n"
        << "###----------------------------------------------------------------###\n";

    if (potential_shape == "PureCoulomb")
    {
        fdouble r_temp = 0.0;
        std_cout << "### Using a pure Coulomb interaction                               ###\n";
        Initialize_Simple(r_temp);

        Potentials_Set_Parameters = &Potentials_Set_Parameters_Simple;
        Calculate_Potential       = &Calculate_Potential_Cutoff_Simple;
        Set_Field                 = &Set_Field_Cutoff_Simple;
    }
    else if (potential_shape == "Simple")
    {
        std_cout << "### Using a simple cutoff                                          ###\n";
        std_cout << "### for close range interaction                                    ###\n";
        Initialize_Simple(input_s_rmin);

        Potentials_Set_Parameters = &Potentials_Set_Parameters_Simple;
        Calculate_Potential       = &Calculate_Potential_Cutoff_Simple;
        Set_Field                 = &Set_Field_Cutoff_Simple;
    }
    else if (potential_shape == "Harmonic")
    {
        std_cout << "### Using an harmonic potential                                    ###\n";
        std_cout << "### for close range interaction                                    ###\n";

        Potentials_Set_Parameters = &Potentials_Set_Parameters_Harmonic;
        Calculate_Potential       = &Calculate_Potential_Cutoff_Harmonic;
        Set_Field                 = &Set_Field_Cutoff_Harmonic;
    }
    else if (potential_shape == "SuperGaussian")
    {
        std_cout << "### Using a super-gaussian potential                               ###\n";
        std_cout << "### for close range interaction                                    ###\n";
        std_cout << "### with m = " << input_sg_m << "\n";
        Initialize_SuperGaussian(input_sg_m);

        Potentials_Set_Parameters = &Potentials_Set_Parameters_SuperGaussian;
        Calculate_Potential       = &Calculate_Potential_Cutoff_SuperGaussian;
        Set_Field                 = &Set_Field_Cutoff_SuperGaussian;
    }
    else if (potential_shape == "GaussianDistribution")
    {
        std_cout << "### Using a gaussian charge distribution                           ###\n";
        std_cout << "### potential for close range interaction                          ###\n";
        std_cout << "### Initializing the lookup tables...                              ###\n" << std::flush;
        libpotentials_private::lut_potential.Initialize(erf_over_x,                 0.0, 4.5*std::sqrt(2.0), 10000, "Potential LookUpTable");
        libpotentials_private::lut_field.Initialize(erf_over_x3_minus_exp_over_x2,  0.0, 4.5*std::sqrt(2.0), 10000, "Field LookUpTable");
        std_cout << "### Initializing the lookup tables done.                           ###\n" << std::flush;

        Potentials_Set_Parameters = &Potentials_Set_Parameters_GaussianDistribution;
        Calculate_Potential       = &Calculate_Potential_Cutoff_GaussianDistribution;
        Set_Field                 = &Set_Field_Cutoff_GaussianDistribution;
    }
    else if (potential_shape == "HermanSkillman")
    {
        std_cout << "### Using the Herman-Skillman (HS) potential                        ##\n";
        std_cout << "### for close range interaction                                    ###\n";
        std_cout << "### and Super-Gaussian for electrons and 8+                        ###\n";
        std_cout << "### and up ions (m = " << input_sg_m << ")                                           ###\n";
        Initialize_HS(input_sg_m, input_s_rmin);

        Potentials_Set_Parameters = &Potentials_Set_Parameters_HS_SuperGaussian;
        Calculate_Potential       = &Calculate_Potential_Cutoff_HS_SuperGaussian;
        Set_Field                 = &Set_Field_Cutoff_HS_SuperGaussian;
    }
    else if (potential_shape == "Symmetric")
    {
        std_cout << "### Using the symmetric two charge                                 ###\n";
        std_cout << "### distribution interaction                                       ###\n";
        std_cout << "### Initializing the lookup tables...                              ###\n" << std::flush;
        libpotentials_private::lut_potential.Initialize(erf_over_x,                 0.0, 4.5*std::sqrt(2.0), 10000, "Potential LookUpTable");
        libpotentials_private::lut_field.Initialize(erf_over_x3_minus_exp_over_x2,  0.0, 4.5*std::sqrt(2.0), 10000, "Field LookUpTable");
        std_cout << "### Initializing the lookup tables done.                           ###\n" << std::flush;
        std_cout <<                                            "Done!...                ###\n";

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

    std_cout << "### Git versioning:                                                ###\n"
             << "###     build_time:   "<< libpotentials::build_time << std::endl
             << "###     build_sha:    "<< libpotentials::build_sha << std::endl
             << "###     build_branch: "<< libpotentials::build_branch << std::endl
             << "###----------------------------------------------------------------###\n"
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
        std_cout << " " << Get_Charge_State(pv) * Get_Potential(pv);
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

    if (x < 1.0)
    {
        value = (
               2.0
            - (2.0 * x*x)                                                                               / 3.0
            +        x*x*x*x                                                                            / 5.0
            -        x*x*x*x*x*x                                                                        / 21.0
            +        x*x*x*x*x*x*x*x                                                                    / 108.0
            -        x*x*x*x*x*x*x*x*x*x                                                                / 660.0
            +        x*x*x*x*x*x*x*x*x*x*x*x                                                            / 4680.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                        / 37800.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                    / 342720.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                / 3447360.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                            / 38102400.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                        / 459043200.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                    / 5987520000.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                / 84064780800.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                            / 1264085222400.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                        / 20268952704000.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                    / 345226033152000.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                / 6224529991680000.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x            / 118443913555968000.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x        / 2372079457972224000.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x    / 49874491167621120000.0
            //   +O(x^41)
        );
        value /= libpotentials::sqrt_Pi;
    }
    else
    {
        value = erf(x) / x;
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
               2.0                                                                                              / 3.0
            - (2.0 * x*x)                                                                                       / 5.0
            +        x*x*x*x                                                                                    / 7.0
            -        x*x*x*x*x*x                                                                                / 27.0
            +        x*x*x*x*x*x*x*x                                                                            / 132.0
            -        x*x*x*x*x*x*x*x*x*x                                                                        / 780.0
            +        x*x*x*x*x*x*x*x*x*x*x*x                                                                    / 5400.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                                / 42840.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                            / 383040.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                        / 3810240.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                    / 41731200.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                / 498960000.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                            / 6466521600.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                        / 90291801600.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                    / 1351263513600.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                / 21576627072000.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                            / 366148823040000.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                        / 6580217419776000.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                    / 124846287261696000.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                / 2493724558381056000.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x            / 52307393175797760000.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x        / 1149546198863462400000.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x    / 26414017102773780480000.0
            // +O(x^45)
        ) / libpotentials::sqrt_Pi;
    }
    else
    {
        value = erf(x)/(2.0*x*x*x) - 1.0/libpotentials::sqrt_Pi*std::exp(-x*x)/(x*x);
    }

    return value;
}


// ********** End of file ***************************************
