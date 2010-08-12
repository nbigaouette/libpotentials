
#include <cstdio>
#include <iostream>

#include "Constants.hpp"
#include "Std_Cout.hpp"

#define F std_cout.Format(20,14,'e')

using namespace libpotentials;

void Print_Constants()
{
    std_cout << "###########################################\n";
    std_cout << "### Constants: ############################\n";
    std_cout << "###########################################\n";
    std_cout << "Pi                  = "; F; std_cout << Pi << "\n";
    std_cout << "twoPi               = "; F; std_cout << twoPi << "\n";
    std_cout << "sqrt_Pi             = "; F; std_cout << sqrt_Pi << "\n";
    std_cout << "sqrt_2              = "; F; std_cout << sqrt_2 << "\n";
    std_cout << "sqrt_2_over_pi      = "; F; std_cout << sqrt_2_over_pi << "\n";
    std_cout << "two_over_sqrt_Pi    = "; F; std_cout << two_over_sqrt_Pi << "\n";
    std_cout << "one_over_Pi         = "; F; std_cout << one_over_Pi << "\n";
    std_cout << "deg_to_rad          = "; F; std_cout << deg_to_rad << "\n";
    std_cout << "rad_to_deg          = "; F; std_cout << rad_to_deg << "\n";

    std_cout << "two_over_three      = "; F; std_cout << two_over_three << "\n";
    std_cout << "three_over_two      = "; F; std_cout << three_over_two << "\n";

    std_cout << "one_over_three      = "; F; std_cout << one_over_three << "\n";
    std_cout << "one_over_four       = "; F; std_cout << one_over_four << "\n";
    std_cout << "one_over_seven      = "; F; std_cout << one_over_seven << "\n";

    std_cout << "one_over_four_Pi    = "; F; std_cout << one_over_four_Pi << "\n";
    std_cout << "sqrt_one_over_four_Pi= "; F; std_cout << sqrt_one_over_four_Pi << "\n";

    std_cout << "four_over_twenty_seven = "; F; std_cout << four_over_twenty_seven << "\n";
    std_cout << "eight_over_three    = "; F; std_cout << eight_over_three << "\n";
    std_cout << "three_over_sixteen  = "; F; std_cout << three_over_sixteen << "\n";

    std_cout << "SI:\n";
    std_cout << "co                  = "; F; std_cout << co << "\n";
    std_cout << "co2                 = "; F; std_cout << co2 << "\n";
    std_cout << "inv_co              = "; F; std_cout << inv_co << "\n";
    std_cout << "inv_co2             = "; F; std_cout << inv_co2 << "\n";
    std_cout << "mu0                 = "; F; std_cout << mu0 << "\n";
    std_cout << "eps0                = "; F; std_cout << eps0 << "\n";
    std_cout << "eta0                = "; F; std_cout << eta0 << "\n";
    std_cout << "one_over_4Pieps0    = "; F; std_cout << one_over_4Pieps0 << "\n";

    std_cout << "e0                  = "; F; std_cout << e0 << "\n";
    std_cout << "qe                  = "; F; std_cout << qe << "\n";
    std_cout << "qp                  = "; F; std_cout << qp << "\n";
    std_cout << "mp                  = "; F; std_cout << mp << "\n";
    std_cout << "mp_eV               = "; F; std_cout << mp_eV << "\n";
    std_cout << "mp_MeV              = "; F; std_cout << mp_MeV << "\n";
    std_cout << "me                  = "; F; std_cout << me << "\n";
    std_cout << "me_eV               = "; F; std_cout << me_eV << "\n";
    std_cout << "me_MeV              = "; F; std_cout << me_MeV << "\n";
    std_cout << "one_over_me         = "; F; std_cout << one_over_me << "\n";
    std_cout << "one_over_mp         = "; F; std_cout << one_over_mp << "\n";

    std_cout << "planck              = "; F; std_cout << planck << "\n";
    std_cout << "hbar                = "; F; std_cout << hbar << "\n";
    std_cout << "a0                  = "; F; std_cout << a0 << "\n";

    std_cout << "Na                  = "; F; std_cout << Na << "\n";

    std_cout << "Atomic units:\n";
    std_cout << "Eh                  = "; F; std_cout << Eh << "\n";
    std_cout << "au_energy           = "; F; std_cout << au_energy << "\n";
    std_cout << "au_time             = "; F; std_cout << au_time << "\n";
    std_cout << "au_electric_field   = "; F; std_cout << au_electric_field << "\n";
    std_cout << "w_au                = "; F; std_cout << w_au << "\n";
    std_cout << "alpha               = "; F; std_cout << alpha << "\n";

    std_cout << "Direct Units Conversions:\n";
    std_cout << "Eh_to_eV            = "; F; std_cout << Eh_to_eV << "\n";
    std_cout << "eV_to_Eh            = "; F; std_cout << eV_to_Eh << "\n";
    std_cout << "eV_to_J             = "; F; std_cout << eV_to_J << "\n";
    std_cout << "J_to_eV             = "; F; std_cout << J_to_eV << "\n";
    std_cout << "J_to_Eh             = "; F; std_cout << J_to_Eh << "\n";
    std_cout << "Eh_to_J             = "; F; std_cout << Eh_to_J << "\n";

    std_cout << "m_to_angstrom       = "; F; std_cout << m_to_angstrom << "\n";
    std_cout << "angstrom_to_m       = "; F; std_cout << angstrom_to_m << "\n";
    std_cout << "m_to_bohr           = "; F; std_cout << m_to_bohr << "\n";
    std_cout << "bohr_to_m           = "; F; std_cout << bohr_to_m << "\n";
    std_cout << "m_to_cm             = "; F; std_cout << m_to_cm << "\n";
    std_cout << "cm_to_m             = "; F; std_cout << cm_to_m << "\n";
    std_cout << "m_to_nm             = "; F; std_cout << m_to_nm << "\n";
    std_cout << "nm_to_m             = "; F; std_cout << nm_to_m << "\n";

    std_cout << "fs_to_s             = "; F; std_cout << fs_to_s << "\n";
    std_cout << "s_to_fs             = "; F; std_cout << s_to_fs << "\n";
    std_cout << "as_to_s             = "; F; std_cout << as_to_s << "\n";
    std_cout << "s_to_as             = "; F; std_cout << s_to_as << "\n";
    std_cout << "deg_to_rad          = "; F; std_cout << deg_to_rad << "\n";
    std_cout << "rad_to_deg          = "; F; std_cout << rad_to_deg << "\n";

    std_cout << "SI  <--> Atomic Units Conversions:\n";
    std_cout << "au_to_si_mass       = "; F; std_cout << au_to_si_mass << "\n";
    std_cout << "si_to_au_mass       = "; F; std_cout << si_to_au_mass << "\n";
    std_cout << "au_to_si_length     = "; F; std_cout << au_to_si_length << "\n";
    std_cout << "si_to_au_length     = "; F; std_cout << si_to_au_length << "\n";
    std_cout << "au_to_si_charge     = "; F; std_cout << au_to_si_charge << "\n";
    std_cout << "si_to_au_charge     = "; F; std_cout << si_to_au_charge << "\n";
    std_cout << "au_to_si_energy     = "; F; std_cout << au_to_si_energy << "\n";
    std_cout << "si_to_au_energy     = "; F; std_cout << si_to_au_energy << "\n";
    std_cout << "au_to_si_k          = "; F; std_cout << au_to_si_k << "\n";
    std_cout << "si_to_au_k          = "; F; std_cout << si_to_au_k << "\n";
    std_cout << "au_to_si_time       = "; F; std_cout << au_to_si_time << "\n";
    std_cout << "si_to_au_time       = "; F; std_cout << si_to_au_time << "\n";
    std_cout << "au_to_si_field      = "; F; std_cout << au_to_si_field << "\n";
    std_cout << "si_to_au_field      = "; F; std_cout << si_to_au_field << "\n";
    std_cout << "au_to_si_force      = "; F; std_cout << au_to_si_force << "\n";
    std_cout << "si_to_au_force      = "; F; std_cout << si_to_au_force << "\n";
    std_cout << "au_to_si_pot        = "; F; std_cout << au_to_si_pot << "\n";
    std_cout << "si_to_au_pot        = "; F; std_cout << si_to_au_pot << "\n";

    std_cout << "###########################################\n";
    std_cout << "###########################################\n";
}

// ********** End of File ***************************************
