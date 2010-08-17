
#include "LibPotentials.hpp"
#include "Constants.hpp"
#include "Potentials.hpp"
#include "Global.hpp"
#include "Std_Cout.hpp"
#include "Code_Functions_Declarations.hpp"

extern void Initialize_Simple(const double &minr);
extern void Initialize_SuperGaussian(const int &m);
extern void Initialize_HS(const int &input_sg_m, const double &min_rad);

// **************************************************************
void Potentials_Initialize(const std::string potential_shape,
                           const double base_potential_depth,
                           const double input_s_rmin,
                           const double input_sg_m)
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
        double r_temp = 0.0;
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
        libpotentials_private::initialize_erf_lookup_table();

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
        Calculate_Potential       = &Calculate_Potential_Cutoff_Harmonic;
        Set_Field                 = &Set_Field_Cutoff_HS_SuperGaussian;
    }
    else if (potential_shape == "Symmetric")
    {
        std_cout << "### Using the symmetric two charge                                 ###\n";
        std_cout << "### distribution interaction                                       ###\n";
        libpotentials_private::initialize_erf_lookup_table();

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
    free(libpotentials_private::tl_erf);
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

// ********** End of file ***************************************
