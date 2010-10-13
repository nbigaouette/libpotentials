#ifndef INC_LIBPOTENTIALS_CONSTANTS_hpp
#define INC_LIBPOTENTIALS_CONSTANTS_hpp

#ifdef __SUNPRO_CC
#include <math.h>
#else
#include <cmath>
#endif

#include "FloatType.hpp"

namespace libpotentials
{

    // The base potential well giving the bottom of the potential well.
    //const fdouble base_pot_well_depth = 26.0; //[eV]

    // **************************************************************

    const fdouble zero              = fdouble(0.0);
    const fdouble one               = fdouble(1.0);
    const fdouble two               = fdouble(2.0);
    const fdouble three             = fdouble(3.0);
    const fdouble four              = fdouble(4.0);
    const fdouble seven             = fdouble(7.0);
    const fdouble eight             = fdouble(8.0);

    const fdouble Pi             = std::acos(-one);
    const fdouble twoPi          = two * Pi;
    const fdouble sqrt_Pi        = std::sqrt(Pi);
    const fdouble sqrt_2         = std::sqrt(two);
    const fdouble sqrt_2_over_pi = std::sqrt(two / Pi);
    const fdouble one_over_Pi    = one / Pi;
    const fdouble two_over_sqrt_Pi = two / sqrt_Pi;

    const fdouble two_over_three = two / three;
    const fdouble three_over_two = three / two;

    const fdouble one_over_three = one / three;
    const fdouble one_over_four  = one / four;
    const fdouble one_over_seven = one / seven;

    const fdouble one_over_four_Pi = one_over_four * one_over_Pi;
    const fdouble sqrt_one_over_four_Pi = std::sqrt(one_over_four_Pi);

    const fdouble four_over_twenty_seven = four / fdouble(27.0);
    const fdouble eight_over_three       = eight / three;
    const fdouble three_over_eight       = three / eight;
    const fdouble three_over_sixteen     = three / fdouble(16.0);

    const int ZERO              = 0;
    const int ONE               = 1;
    const int TWO               = 2;

    // **************************************************************
    // ********** http://physics.nist.gov/cuu/Constants *************
    // **************************************************************


    // **************************************************************
    // ********** SI

    const fdouble co             = fdouble(299792458.0);          // Speed of ligth in vacuum [m/s]
    const fdouble co2            = co*co;                // [m^2/s^2]
    const fdouble inv_co         = one / co;             // 1/co [s/m]
    const fdouble inv_co2        = one / (co * co);      // 1/co^2 [s^2/m^2]
    const fdouble mu0            = fdouble(4.0e-7) * Pi;          // Magnetic constant (vacuum) [N/A^2]
    const fdouble eps0           = inv_co2 / mu0;        // Electric constant (vacuum) [F/m]
    const fdouble eta0           = mu0 * co;             // Impedance of vacuum [Ohm]
    const fdouble one_over_4Pieps0                       // Coulomb's law
                                = one_over_four_Pi / eps0;// constant
                                                        // [V.m.C^-1]

    const fdouble e0             = fdouble(1.602176462e-19);      // Unitary charge [C]
    const fdouble qe             = -e0;                  // Electron charge [C]
    const fdouble qp             =  e0;                  // Proton charge [C]
    const fdouble mp             = fdouble(1.67262158e-27);       // Proton mass [kg]
    const fdouble mp_eV          = mp * co*co / e0;      // Proton mass energy (mc^2) [eV]
    const fdouble mp_MeV         = mp_eV / fdouble(1.0e6);        // Proton mass energy (mc^2) [MeV]
    const fdouble me             = fdouble(9.10938188e-31);       // Electron mass [kg]
    const fdouble me_eV          = me * co*co / e0;      // Electron mass energy (mc^2) [eV]
    const fdouble me_MeV         = me_eV / fdouble(1.0e6);        // Electron mass energy (mc^2) [MeV]
    const fdouble one_over_me    = one / me;             // [1/kg]
    const fdouble one_over_mp    = one / mp;             // [1/kg]

    const fdouble planck         = fdouble(6.62606896e-34);       // Planck constant [J . s]
    const fdouble hbar           = planck / twoPi;       // Planck cst / 2 pi [J . s]
    const fdouble a0             = hbar*hbar / (one_over_4Pieps0 * me * e0*e0);
                                                        // Bohr radius [m]
    const fdouble Na             = fdouble(6.02214179e23);        // Avogadro constant [mol^-1]

    const fdouble kB             = fdouble(1.3806504e-23);        // Boltzmann constant  J . K^-1

    // **************************************************************
    // ********** Atomic Units

    const fdouble Eh             = hbar*hbar/(me*a0*a0); // ...of energy (Hartree) [J]
    const fdouble au_energy      = Eh;                   // ...of energy (Hartree) [J]
    const fdouble au_time        = hbar / Eh;            // ...of time [s]
    const fdouble au_electric_field = Eh / (e0 * a0);    // ...of electric field [V . m^-1]
    const fdouble w_au           = one / au_time;        // ...of frequency [???]
    const fdouble alpha          = fdouble(7.29735257e-3);        //fine-structure constant

    // **************************************************************
    // ********** Conversions

    const fdouble Eh_to_eV       = Eh / e0;              // eV . Eh^-1
    const fdouble eV_to_Eh       = one / Eh_to_eV;       // Eh . eV^-1
    const fdouble eV_to_J        = e0;                   // J . eV^-1
    const fdouble J_to_eV        = one / eV_to_J;        // eV . J^-1
    const fdouble J_to_Eh        = J_to_eV * eV_to_Eh;   // Eh . J^-1
    const fdouble Eh_to_J        = one / J_to_Eh;        // J . Eh^-1

    const fdouble m_to_angstrom  = fdouble(1.0e10);               // Å . m^-1
    const fdouble angstrom_to_m  = one / m_to_angstrom;  // m . Å^-1
    const fdouble bohr_to_m      = a0;                   // m . bohr^-1
    const fdouble m_to_bohr      = one / bohr_to_m;      // bohr . m^-1
    const fdouble m_to_cm        = fdouble(100.0);                // cm . m^-1
    const fdouble cm_to_m        = one / m_to_cm;        // cm . m^-1
    const fdouble m_to_nm        = fdouble(1.0e9);                // nm . m^-1
    const fdouble nm_to_m        = one / m_to_nm;        // m . nm^-1
    const fdouble bohr_to_angstrom = bohr_to_m * m_to_angstrom;  // Å . bohr^-1
    const fdouble angstrom_to_bohr = one / bohr_to_angstrom;     // bohr . Å^-1

    const fdouble as_to_s        = fdouble(1.0e-18);              // s . as^-1
    const fdouble s_to_as        = one / as_to_s;        // as . s^-1
    const fdouble fs_to_s        = fdouble(1.0e-15);              // s . fs^-1
    const fdouble s_to_fs        = one / fs_to_s;        // fs . s^-1
    const fdouble ps_to_s        = fdouble(1.0e-12);              // s . ps^-1
    const fdouble s_to_ps        = one / ps_to_s;        // ps . s^-1
    const fdouble ns_to_s        = fdouble(1.0e-9);               // s . ns^-1
    const fdouble s_to_ns        = one / ns_to_s;        // ns . s^-1
    const fdouble mus_to_s       = fdouble(1.0e-6);               // s . mus^-1
    const fdouble s_to_mus       = one / mus_to_s;       // mus . s^-1
    const fdouble ms_to_s        = fdouble(1.0e-3);               // s . ms^-1
    const fdouble s_to_ms        = one / ms_to_s;        // ms . s^-1

    const fdouble deg_to_rad     = two * Pi / fdouble(360.0);     // radians . degrees^-1
    const fdouble rad_to_deg     = one / deg_to_rad;     // degrees . radians^-1

    const fdouble byte_to_KiB    = one / fdouble(1024.0);
    const fdouble byte_to_MiB    = byte_to_KiB / fdouble(1024.0);
    const fdouble byte_to_GiB    = byte_to_MiB / fdouble(1024.0);

    // ********** SI  <--> Atomic Units Conversions
    const fdouble au_to_si_mass  = me;                   // kg . au^-1
    const fdouble si_to_au_mass  = one / au_to_si_mass;  // au . kg^-1

    const fdouble au_to_si_length= a0;                   // m . bohr^-1
    const fdouble si_to_au_length= one / au_to_si_length;// bohr . m^-1

    const fdouble au_to_si_charge= e0;                   // C . au^-1
    const fdouble si_to_au_charge= one / au_to_si_charge;// au . C^-1

    const fdouble au_to_si_energy= Eh;                   // J . au^-1
    const fdouble si_to_au_energy= one / au_to_si_energy;// au . J^-1

    const fdouble au_to_si_k     = one_over_4Pieps0;     // C^-2 . N . m^2 . au^-1
    const fdouble si_to_au_k     = one / au_to_si_energy;// au . C^2 . N^-1 . m^-2

    const fdouble au_to_si_time  = hbar / Eh;            // s . au^-1
    const fdouble si_to_au_time  = one / au_to_si_time;  // au . s^-1

    const fdouble au_to_si_force = Eh / a0;              // N . au^-1
    const fdouble si_to_au_force = one / au_to_si_force; // au . N^-1

    const fdouble au_to_si_field = au_electric_field;    // V . m^-1 . au^-1
    const fdouble si_to_au_field = one / au_to_si_field; // au . m . V^-1

    const fdouble au_to_si_pot   = au_electric_field * a0; // V . au^-1
    const fdouble si_to_au_pot   = one / au_to_si_pot;   // au . V^-1

    const fdouble au_to_si_vel   = Eh * a0 / hbar;       // m.s^-1 . au^-1
    const fdouble si_to_au_vel   = one / au_to_si_vel;   // au . m.s^-1


    // ********** Other Usefull Conversions
    // Cross sections (area): Barn (b), Megabarnes (Mb) http://en.wikipedia.org/wiki/Barn_(unit)
    const fdouble b_to_m2        = fdouble(1.0e-28);              // m^2 . b^-1
    const fdouble m2_to_b        = one / b_to_m2;        // b . m^-2
    const fdouble Mb_to_m2       = fdouble(1.0e6) * b_to_m2;      // m^2 . Mb^-1
    const fdouble m2_to_Mb       = one / Mb_to_m2;       // Mb . m^-2
    const fdouble m2_to_atomic_area = si_to_au_length * si_to_au_length; // bohr^2 . m^-2
    const fdouble atomic_area_to_m2 = one / m2_to_atomic_area;           // m^2 . bohr^-2
    const fdouble Mb_to_atomic_area = Mb_to_m2 * m2_to_atomic_area;      // bohr^2 . Mb^-1
    const fdouble atomic_area_to_Mb = one / Mb_to_atomic_area;           // Mb . bohr^-2

    // **************************************************************
    // **********  Print all constants
    void Print_Constants();
}

#endif // INC_LIBPOTENTIALS_CONSTANTS_hpp

// ********** End of File ***************************************
