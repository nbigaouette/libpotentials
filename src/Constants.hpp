#ifndef INC_LIBPOTENTIALS_CONSTANTS_hpp
#define INC_LIBPOTENTIALS_CONSTANTS_hpp

#ifdef __SUNPRO_CC
#include <math.h>
#else
#include <cmath>
#endif

#include "FloatType.hpp"

/**
 * The constants are defined as double precision.
 * This is NECESSARY to prevent single precision (floats) overflowing!
 * The calculations are done in double in the libpotentials_double
 * namespace, then casted to single in the libpotentials namespace.
 */


void Print_Constants();

namespace libpotentials_double
{
    // **************************************************************

    const double zero                   = 0.0;
    const double half                   = 0.5;
    const double one                    = 1.0;
    const double two                    = 2.0;
    const double three                  = 3.0;
    const double four                   = 4.0;
    const double six                    = 6.0;
    const double seven                  = 7.0;
    const double eight                  = 8.0;

    const double Pi                     = std::acos(-one);
    const double twoPi                  = two * Pi;
    const double sqrt_Pi                = std::sqrt(Pi);
    const double sqrt_2                 = std::sqrt(two);
    const double sqrt_2_over_pi         = std::sqrt(two / Pi);
    const double one_over_Pi            = one / Pi;
    const double two_over_sqrt_Pi       = two / sqrt_Pi;

    const double two_over_three         = two / three;
    const double three_over_two         = three / two;

    const double one_over_three         = one / three;
    const double one_over_four          = one / four;
    const double one_over_seven         = one / seven;

    const double one_over_four_Pi       = one_over_four * one_over_Pi;
    const double sqrt_one_over_four_Pi  = std::sqrt(one_over_four_Pi);

    const double four_over_twenty_seven = four / 27.0;
    const double eight_over_three       = eight / three;
    const double three_over_eight       = three / eight;
    const double three_over_sixteen     = three / 16.0;

    // **************************************************************
    // ********** http://physics.nist.gov/cuu/Constants *************
    // **************************************************************


    // **************************************************************
    // ********** SI

    const double co                     = 299792458.0;          // Speed of ligth in vacuum [m/s]
    const double co2                    = co*co;                // [m^2/s^2]
    const double inv_co                 = one / co;             // 1/co [s/m]
    const double inv_co2                = one / (co * co);      // 1/co^2 [s^2/m^2]
    const double mu0                    = 4.0e-7 * Pi;          // Magnetic constant (vacuum) [N/A^2]
    const double eps0                   = inv_co2 / mu0;        // Electric constant (vacuum) [F/m]
    const double eta0                   = mu0 * co;             // Impedance of vacuum [Ohm]
    const double one_over_4Pieps0       = one_over_four_Pi / eps0;// Coulomb's law constant [V.m.C^-1]

    const double e0                     = 1.602176462e-19;      // Unitary charge [C]
    const double qe                     = -e0;                  // Electron charge [C]
    const double qp                     =  e0;                  // Proton charge [C]
    const double mp                     = 1.67262158e-27;       // Proton mass [kg]
    const double mp_eV                  = mp * co*co / e0;      // Proton mass energy (mc^2) [eV]
    const double mp_MeV                 = mp_eV / 1.0e6;        // Proton mass energy (mc^2) [MeV]
    const double me                     = 9.10938188e-31;       // Electron mass [kg]
    const double me_eV                  = me * co*co / e0;      // Electron mass energy (mc^2) [eV]
    const double me_MeV                 = me_eV / 1.0e6;        // Electron mass energy (mc^2) [MeV]
    const double one_over_me            = one / me;             // [1/kg]
    const double one_over_mp            = one / mp;             // [1/kg]

    const double planck                 = 6.62606896e-34;       // Planck constant [J . s]
    const double hbar                   = planck / twoPi;       // Planck cst / 2 pi [J . s]
    const double a0                     = hbar*hbar / (one_over_4Pieps0 * me * e0*e0);
                                                                // Bohr radius [m]
    const double Na                     = 6.02214179e23;        // Avogadro constant [mol^-1]

    const double kB                     = 1.3806504e-23;        // Boltzmann constant  J . K^-1

    const double alpha                  = 7.29735257e-3;        // Fine-structure constant [-]

    // **************************************************************
    // ********** Atomic Units

    const double Eh                     = hbar*hbar/(me*a0*a0); // ...of energy (Hartree) [J]
    const double au_energy              = Eh;                   // ...of energy (Hartree) [J]
    const double au_time                = hbar / Eh;            // ...of time [s]
    const double au_electric_field      = Eh / (e0 * a0);       // ...of electric field [V . m^-1]
    const double w_au                   = one / au_time;        // ...of frequency [???]

    // **************************************************************
    // ********** Conversions

    const double Eh_to_eV               = Eh / e0;              // eV . Eh^-1
    const double eV_to_Eh               = one / Eh_to_eV;       // Eh . eV^-1
    const double eV_to_J                = e0;                   // J . eV^-1
    const double J_to_eV                = one / eV_to_J;        // eV . J^-1
    const double J_to_Eh                = J_to_eV * eV_to_Eh;   // Eh . J^-1
    const double Eh_to_J                = one / J_to_Eh;        // J . Eh^-1

    const double m_to_angstrom          = 1.0e10;               // Å . m^-1
    const double angstrom_to_m          = one / m_to_angstrom;  // m . Å^-1
    const double bohr_to_m              = a0;                   // m . bohr^-1
    const double m_to_bohr              = one / bohr_to_m;      // bohr . m^-1
    const double m_to_cm                = 100.0;                // cm . m^-1
    const double cm_to_m                = one / m_to_cm;        // cm . m^-1
    const double m_to_nm                = 1.0e9;                // nm . m^-1
    const double nm_to_m                = one / m_to_nm;        // m . nm^-1
    const double bohr_to_angstrom       = bohr_to_m * m_to_angstrom;  // Å . bohr^-1
    const double angstrom_to_bohr       = one / bohr_to_angstrom;     // bohr . Å^-1

    const double as_to_s                = 1.0e-18;              // s . as^-1
    const double s_to_as                = one / as_to_s;        // as . s^-1
    const double fs_to_s                = 1.0e-15;              // s . fs^-1
    const double s_to_fs                = one / fs_to_s;        // fs . s^-1
    const double ps_to_s                = 1.0e-12;              // s . ps^-1
    const double s_to_ps                = one / ps_to_s;        // ps . s^-1
    const double ns_to_s                = 1.0e-9;               // s . ns^-1
    const double s_to_ns                = one / ns_to_s;        // ns . s^-1
    const double mus_to_s               = 1.0e-6;               // s . mus^-1
    const double s_to_mus               = one / mus_to_s;       // mus . s^-1
    const double ms_to_s                = 1.0e-3;               // s . ms^-1
    const double s_to_ms                = one / ms_to_s;        // ms . s^-1

    const double deg_to_rad             = two * Pi / 360.0;     // radians . degrees^-1
    const double rad_to_deg             = one / deg_to_rad;     // degrees . radians^-1

    const double byte_to_KiB            = one / 1024.0;
    const double byte_to_MiB            = byte_to_KiB / 1024.0;
    const double byte_to_GiB            = byte_to_MiB / 1024.0;

    // ********** SI  <--> Atomic Units Conversions
    const double au_to_si_mass          = me;                   // kg . au^-1
    const double si_to_au_mass          = one / au_to_si_mass;  // au . kg^-1

    const double au_to_si_length        = a0;                   // m . bohr^-1
    const double si_to_au_length        = one / au_to_si_length;// bohr . m^-1

    const double au_to_si_charge        = e0;                   // C . au^-1
    const double si_to_au_charge        = one / au_to_si_charge;// au . C^-1

    const double au_to_si_energy        = Eh;                   // J . au^-1
    const double si_to_au_energy        = one / au_to_si_energy;// au . J^-1

    const double au_to_si_k             = one_over_4Pieps0;     // C^-2 . N . m^2 . au^-1
    const double si_to_au_k             = one / au_to_si_energy;// au . C^2 . N^-1 . m^-2

    const double au_to_si_time          = hbar / Eh;            // s . au^-1
    const double si_to_au_time          = one / au_to_si_time;  // au . s^-1

    const double au_to_si_force         = Eh / a0;              // N . au^-1
    const double si_to_au_force         = one / au_to_si_force; // au . N^-1

    const double au_to_si_field         = au_electric_field;    // V . m^-1 . au^-1
    const double si_to_au_field         = one / au_to_si_field; // au . m . V^-1

    const double au_to_si_pot           = au_electric_field * a0; // V . au^-1
    const double si_to_au_pot           = one / au_to_si_pot;   // au . V^-1

    const double au_to_si_vel           = Eh * a0 / hbar;       // m.s^-1 . bohr . aut^-1
    const double si_to_au_vel           = one / au_to_si_vel;   // au . m.s^-1

    const double au_to_si_acc           = au_to_si_vel / au_to_si_time; // m.s^-2 . bohr^-1.aut^2
    const double si_to_au_acc           = one / au_to_si_acc;   // au . m.s^-1

    const double au_to_si_vpot          = au_to_si_pot * au_to_si_time / au_to_si_length; // (V.s.m^-1) . au^-1
    const double si_to_au_vpot          = 1.0 / au_to_si_vpot;  // au . (V.s.m^-1)^-1

    // ********** Other Usefull Conversions
    // Cross sections (area): Barn (b), Megabarnes (Mb) http://en.wikipedia.org/wiki/Barn_(unit)
    const double b_to_m2                = 1.0e-28;              // m^2 . b^-1
    const double m2_to_b                = one / b_to_m2;        // b . m^-2
    const double Mb_to_m2               = 1.0e6 * b_to_m2;      // m^2 . Mb^-1
    const double m2_to_Mb               = one / Mb_to_m2;       // Mb . m^-2
    const double m2_to_atomic_area      = si_to_au_length * si_to_au_length; // bohr^2 . m^-2
    const double atomic_area_to_m2      = one / m2_to_atomic_area;           // m^2 . bohr^-2
    const double Mb_to_atomic_area      = Mb_to_m2 * m2_to_atomic_area;      // bohr^2 . Mb^-1
    const double atomic_area_to_Mb      = one / Mb_to_atomic_area;           // Mb . bohr^-2

    // sin^4 intensity profile. NOTE: This period contains two pulses!
    // See md.git/scripts/sin4_laser.py and notes.git/xournal/20111104_14h22_Sine4_Intensity_Profile.xoj
    const double acos_pow_half_fourth   = std::acos(std::pow(0.5, 0.25));
    const double sin4_fwhm_to_period    = Pi / acos_pow_half_fourth;
}

namespace libpotentials
{
    // To compare with the return value of strcmp()
    // 0 = same strings, 1 = different strings
    // http://www.cplusplus.com/reference/clibrary/cstring/strcmp/
    const int strcmp_success = 0;

    const fdouble zero                  = fdouble(libpotentials_double::zero                   );
    const fdouble half                  = fdouble(libpotentials_double::half                   );
    const fdouble one                   = fdouble(libpotentials_double::one                    );
    const fdouble two                   = fdouble(libpotentials_double::two                    );
    const fdouble three                 = fdouble(libpotentials_double::three                  );
    const fdouble four                  = fdouble(libpotentials_double::four                   );
    const fdouble six                   = fdouble(libpotentials_double::six                    );
    const fdouble seven                 = fdouble(libpotentials_double::seven                  );
    const fdouble eight                 = fdouble(libpotentials_double::eight                  );
    const fdouble Pi                    = fdouble(libpotentials_double::Pi                     );
    const fdouble twoPi                 = fdouble(libpotentials_double::twoPi                  );
    const fdouble sqrt_Pi               = fdouble(libpotentials_double::sqrt_Pi                );
    const fdouble sqrt_2                = fdouble(libpotentials_double::sqrt_2                 );
    const fdouble sqrt_2_over_pi        = fdouble(libpotentials_double::sqrt_2_over_pi         );
    const fdouble one_over_Pi           = fdouble(libpotentials_double::one_over_Pi            );
    const fdouble two_over_sqrt_Pi      = fdouble(libpotentials_double::two_over_sqrt_Pi       );
    const fdouble two_over_three        = fdouble(libpotentials_double::two_over_three         );
    const fdouble three_over_two        = fdouble(libpotentials_double::three_over_two         );
    const fdouble one_over_three        = fdouble(libpotentials_double::one_over_three         );
    const fdouble one_over_four         = fdouble(libpotentials_double::one_over_four          );
    const fdouble one_over_seven        = fdouble(libpotentials_double::one_over_seven         );
    const fdouble one_over_four_Pi      = fdouble(libpotentials_double::one_over_four_Pi       );
    const fdouble sqrt_one_over_four_Pi = fdouble(libpotentials_double::sqrt_one_over_four_Pi  );
    const fdouble four_over_twenty_seven= fdouble(libpotentials_double::four_over_twenty_seven );
    const fdouble eight_over_three      = fdouble(libpotentials_double::eight_over_three       );
    const fdouble three_over_eight      = fdouble(libpotentials_double::three_over_eight       );
    const fdouble three_over_sixteen    = fdouble(libpotentials_double::three_over_sixteen     );
    const fdouble co                    = fdouble(libpotentials_double::co                     );
    const fdouble co2                   = fdouble(libpotentials_double::co2                    );
    const fdouble inv_co                = fdouble(libpotentials_double::inv_co                 );
    const fdouble inv_co2               = fdouble(libpotentials_double::inv_co2                );
    const fdouble mu0                   = fdouble(libpotentials_double::mu0                    );
    const fdouble eps0                  = fdouble(libpotentials_double::eps0                   );
    const fdouble eta0                  = fdouble(libpotentials_double::eta0                   );
    const fdouble one_over_4Pieps0      = fdouble(libpotentials_double::one_over_4Pieps0       );
    const fdouble e0                    = fdouble(libpotentials_double::e0                     );
    const fdouble qe                    = fdouble(libpotentials_double::qe                     );
    const fdouble qp                    = fdouble(libpotentials_double::qp                     );
    const fdouble mp                    = fdouble(libpotentials_double::mp                     );
    const fdouble mp_eV                 = fdouble(libpotentials_double::mp_eV                  );
    const fdouble mp_MeV                = fdouble(libpotentials_double::mp_MeV                 );
    const fdouble me                    = fdouble(libpotentials_double::me                     );
    const fdouble me_eV                 = fdouble(libpotentials_double::me_eV                  );
    const fdouble me_MeV                = fdouble(libpotentials_double::me_MeV                 );
    const fdouble one_over_me           = fdouble(libpotentials_double::one_over_me            );
    const fdouble one_over_mp           = fdouble(libpotentials_double::one_over_mp            );
    const fdouble planck                = fdouble(libpotentials_double::planck                 );
    const fdouble hbar                  = fdouble(libpotentials_double::hbar                   );
    const fdouble a0                    = fdouble(libpotentials_double::a0                     );
    const fdouble Na                    = fdouble(libpotentials_double::Na                     );
    const fdouble kB                    = fdouble(libpotentials_double::kB                     );
    const fdouble Eh                    = fdouble(libpotentials_double::Eh                     );
    const fdouble au_energy             = fdouble(libpotentials_double::au_energy              );
    const fdouble au_time               = fdouble(libpotentials_double::au_time                );
    const fdouble au_electric_field     = fdouble(libpotentials_double::au_electric_field      );
    const fdouble w_au                  = fdouble(libpotentials_double::w_au                   );
    const fdouble alpha                 = fdouble(libpotentials_double::alpha                  );
    const fdouble Eh_to_eV              = fdouble(libpotentials_double::Eh_to_eV               );
    const fdouble eV_to_Eh              = fdouble(libpotentials_double::eV_to_Eh               );
    const fdouble eV_to_J               = fdouble(libpotentials_double::eV_to_J                );
    const fdouble J_to_eV               = fdouble(libpotentials_double::J_to_eV                );
    const fdouble J_to_Eh               = fdouble(libpotentials_double::J_to_Eh                );
    const fdouble Eh_to_J               = fdouble(libpotentials_double::Eh_to_J                );
    const fdouble m_to_angstrom         = fdouble(libpotentials_double::m_to_angstrom          );
    const fdouble angstrom_to_m         = fdouble(libpotentials_double::angstrom_to_m          );
    const fdouble bohr_to_m             = fdouble(libpotentials_double::bohr_to_m              );
    const fdouble m_to_bohr             = fdouble(libpotentials_double::m_to_bohr              );
    const fdouble m_to_cm               = fdouble(libpotentials_double::m_to_cm                );
    const fdouble cm_to_m               = fdouble(libpotentials_double::cm_to_m                );
    const fdouble m_to_nm               = fdouble(libpotentials_double::m_to_nm                );
    const fdouble nm_to_m               = fdouble(libpotentials_double::nm_to_m                );
    const fdouble bohr_to_angstrom      = fdouble(libpotentials_double::bohr_to_angstrom       );
    const fdouble angstrom_to_bohr      = fdouble(libpotentials_double::angstrom_to_bohr       );
    const fdouble as_to_s               = fdouble(libpotentials_double::as_to_s                );
    const fdouble s_to_as               = fdouble(libpotentials_double::s_to_as                );
    const fdouble fs_to_s               = fdouble(libpotentials_double::fs_to_s                );
    const fdouble s_to_fs               = fdouble(libpotentials_double::s_to_fs                );
    const fdouble ps_to_s               = fdouble(libpotentials_double::ps_to_s                );
    const fdouble s_to_ps               = fdouble(libpotentials_double::s_to_ps                );
    const fdouble ns_to_s               = fdouble(libpotentials_double::ns_to_s                );
    const fdouble s_to_ns               = fdouble(libpotentials_double::s_to_ns                );
    const fdouble mus_to_s              = fdouble(libpotentials_double::mus_to_s               );
    const fdouble s_to_mus              = fdouble(libpotentials_double::s_to_mus               );
    const fdouble ms_to_s               = fdouble(libpotentials_double::ms_to_s                );
    const fdouble s_to_ms               = fdouble(libpotentials_double::s_to_ms                );
    const fdouble deg_to_rad            = fdouble(libpotentials_double::deg_to_rad             );
    const fdouble rad_to_deg            = fdouble(libpotentials_double::rad_to_deg             );
    const fdouble byte_to_KiB           = fdouble(libpotentials_double::byte_to_KiB            );
    const fdouble byte_to_MiB           = fdouble(libpotentials_double::byte_to_MiB            );
    const fdouble byte_to_GiB           = fdouble(libpotentials_double::byte_to_GiB            );
    const fdouble au_to_si_mass         = fdouble(libpotentials_double::au_to_si_mass          );
    const fdouble si_to_au_mass         = fdouble(libpotentials_double::si_to_au_mass          );
    const fdouble au_to_si_length       = fdouble(libpotentials_double::au_to_si_length        );
    const fdouble si_to_au_length       = fdouble(libpotentials_double::si_to_au_length        );
    const fdouble au_to_si_charge       = fdouble(libpotentials_double::au_to_si_charge        );
    const fdouble si_to_au_charge       = fdouble(libpotentials_double::si_to_au_charge        );
    const fdouble au_to_si_energy       = fdouble(libpotentials_double::au_to_si_energy        );
    const fdouble si_to_au_energy       = fdouble(libpotentials_double::si_to_au_energy        );
    const fdouble au_to_si_k            = fdouble(libpotentials_double::au_to_si_k             );
    const fdouble si_to_au_k            = fdouble(libpotentials_double::si_to_au_k             );
    const fdouble au_to_si_time         = fdouble(libpotentials_double::au_to_si_time          );
    const fdouble si_to_au_time         = fdouble(libpotentials_double::si_to_au_time          );
    const fdouble au_to_si_force        = fdouble(libpotentials_double::au_to_si_force         );
    const fdouble si_to_au_force        = fdouble(libpotentials_double::si_to_au_force         );
    const fdouble au_to_si_field        = fdouble(libpotentials_double::au_to_si_field         );
    const fdouble si_to_au_field        = fdouble(libpotentials_double::si_to_au_field         );
    const fdouble au_to_si_pot          = fdouble(libpotentials_double::au_to_si_pot           );
    const fdouble si_to_au_pot          = fdouble(libpotentials_double::si_to_au_pot           );
    const fdouble au_to_si_vel          = fdouble(libpotentials_double::au_to_si_vel           );
    const fdouble si_to_au_vel          = fdouble(libpotentials_double::si_to_au_vel           );
    const fdouble au_to_si_acc          = fdouble(libpotentials_double::au_to_si_acc           );
    const fdouble si_to_au_acc          = fdouble(libpotentials_double::si_to_au_acc           );
    const fdouble au_to_si_vpot         = fdouble(libpotentials_double::au_to_si_vpot          );
    const fdouble si_to_au_vpot         = fdouble(libpotentials_double::si_to_au_vpot          );
    const fdouble b_to_m2               = fdouble(libpotentials_double::b_to_m2                );
    const fdouble m2_to_b               = fdouble(libpotentials_double::m2_to_b                );
    const fdouble Mb_to_m2              = fdouble(libpotentials_double::Mb_to_m2               );
    const fdouble m2_to_Mb              = fdouble(libpotentials_double::m2_to_Mb               );
    const fdouble m2_to_atomic_area     = fdouble(libpotentials_double::m2_to_atomic_area      );
    const fdouble atomic_area_to_m2     = fdouble(libpotentials_double::atomic_area_to_m2      );
    const fdouble Mb_to_atomic_area     = fdouble(libpotentials_double::Mb_to_atomic_area      );
    const fdouble atomic_area_to_Mb     = fdouble(libpotentials_double::atomic_area_to_Mb      );
    const fdouble sin4_fwhm_to_period   = fdouble(libpotentials_double::sin4_fwhm_to_period    );
    const fdouble acos_pow_half_fourth  = fdouble(libpotentials_double::acos_pow_half_fourth   );
}


#endif // INC_LIBPOTENTIALS_CONSTANTS_hpp

// *******f*** End of File ***************************************
