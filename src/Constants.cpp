
#include <cstdio>
#include <iostream>

#include "Constants.hpp"
#include <StdCout.hpp>

// Quote something, usefull to quote a macro's value
#define _QUOTEME(x) #x
#define QUOTEME(x) _QUOTEME(x)

template <class T>
inline void Assert_isinf_isnan(T value)
{
#ifndef __SUNPRO_CC_
     if (std::isinf(value))
     {
         std_cout << "value is inf!!! value = " << value << "\n";
         std_cout << "Aborting\n" << std::flush;
     }
    assert(!std::isinf(value));
#endif // #ifndef __SUNPRO_CC
     if (std::isnan(value))
     {
         std_cout << "value is NaN!!! value = " << value << "\n";
         std_cout << "Aborting\n" << std::flush;
     }
    assert(!std::isnan(value));
}


#define Print(x)                                                                    \
    std_cout.Format(24,0,'d','l');  std_cout << QUOTEME(x) << " = ";                \
    std_cout.Format(23,16,'e');     std_cout << libpotentials_double::x << "  (";   \
    std_cout.Format(15,8,'e');      std_cout << libpotentials::x << ")" << "\n";    \
    Assert_isinf_isnan(libpotentials_double::x);                                    \
    Assert_isinf_isnan(libpotentials::x);

void Print_Constants()
{
    std_cout << "######################################################################\n";
    std_cout << "########################### Constants ################################\n";
    std_cout << "######################################################################\n";
    Print(Pi);
    Print(zero);
    Print(one);
    Print(half);
    Print(two);
    Print(three);
    Print(four);
    Print(seven);
    Print(eight);
    Print(Pi);
    Print(twoPi);
    Print(sqrt_Pi);
    Print(sqrt_2);
    Print(sqrt_2_over_pi);
    Print(one_over_Pi);
    Print(two_over_sqrt_Pi);
    Print(two_over_three);
    Print(three_over_two);
    Print(one_over_three);
    Print(one_over_four);
    Print(one_over_seven);
    Print(one_over_four_Pi);
    Print(sqrt_one_over_four_Pi);
    Print(four_over_twenty_seven);
    Print(eight_over_three);
    Print(three_over_eight);
    Print(three_over_sixteen);
    std_cout << "----------------------------------------------------------------------\n";
    std_cout << "                            SI                                        \n";
    Print(co);
    Print(co2);
    Print(inv_co);
    Print(inv_co2);
    Print(mu0);
    Print(eps0);
    Print(eta0);
    Print(one_over_4Pieps0);
    Print(e0);
    Print(qe);
    Print(qp);
    Print(mp);
    Print(mp_eV);
    Print(mp_MeV);
    Print(me);
    Print(me_eV);
    Print(me_MeV);
    Print(one_over_me);
    Print(one_over_mp);
    Print(planck);
    Print(hbar);
    Print(a0);
    Print(Na);
    Print(kB);
    std_cout << "----------------------------------------------------------------------\n";
    std_cout << "                        Atomic units                                  \n";
    Print(Eh);
    Print(au_energy);
    Print(au_time);
    Print(au_electric_field);
    Print(w_au);
    Print(alpha);
    std_cout << "----------------------------------------------------------------------\n";
    std_cout << "                    Direct Units Conversions                          \n";
    Print(Eh_to_eV);
    Print(eV_to_Eh);
    Print(eV_to_J);
    Print(J_to_eV);
    Print(J_to_Eh);
    Print(Eh_to_J);
    Print(m_to_angstrom);
    Print(angstrom_to_m);
    Print(bohr_to_m);
    Print(m_to_bohr);
    Print(m_to_cm);
    Print(cm_to_m);
    Print(m_to_nm);
    Print(nm_to_m);
    Print(bohr_to_angstrom);
    Print(angstrom_to_bohr);
    Print(as_to_s);
    Print(s_to_as);
    Print(fs_to_s);
    Print(s_to_fs);
    Print(ps_to_s);
    Print(s_to_ps);
    Print(ns_to_s);
    Print(s_to_ns);
    Print(mus_to_s);
    Print(s_to_mus);
    Print(ms_to_s);
    Print(s_to_ms);
    Print(deg_to_rad);
    Print(rad_to_deg);
    Print(byte_to_KiB);
    Print(byte_to_MiB);
    Print(byte_to_GiB);
    std_cout << "----------------------------------------------------------------------\n";
    std_cout << "                SI  <--> Atomic Units Conversions                     \n";
    Print(au_to_si_mass);
    Print(si_to_au_mass);
    Print(au_to_si_length);
    Print(si_to_au_length);
    Print(au_to_si_charge);
    Print(si_to_au_charge);
    Print(au_to_si_energy);
    Print(si_to_au_energy);
    Print(au_to_si_k);
    Print(si_to_au_k);
    Print(au_to_si_time);
    Print(si_to_au_time);
    Print(au_to_si_force);
    Print(si_to_au_force);
    Print(au_to_si_field);
    Print(si_to_au_field);
    Print(au_to_si_pot);
    Print(si_to_au_pot);
    Print(au_to_si_vel);
    Print(si_to_au_vel);
    Print(au_to_si_acc);
    Print(au_to_si_acc);
    std_cout << "----------------------------------------------------------------------\n";
    std_cout << " Laser Pulse: E(t) = A*sin(pi/T*t)^2, I(t) = E(t)^2 = A*sin(pi/T*t)^4 \n";
    Print(sin4_fwhm_to_period);
    std_cout << "######################################################################\n";
    std_cout << "######################################################################\n";
}

// ********** End of File ***************************************
