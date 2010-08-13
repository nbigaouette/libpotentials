#ifndef INC_POTENTIAL_GENERAL_hpp
#define INC_POTENTIAL_GENERAL_hpp

#include <iostream>
#include <cstring>

#include <cassert>
#include <cmath>
#include <cstdlib> // calloc()
#include <cstdio>

#include "Std_Cout.hpp"

#define DEBUGP(x)           std_cout << __FILE__ << ":" << __LINE__ << ":\n    " << (x)

void Check_if_LibPotentials_is_initialized(void);
void Potentials_Initialize(const std::string potential_shape,
                           const double base_potential_depth,
                           const double input_s_rmin,
                           const double input_sg_m);
void Potentials_Finalize();

void Print_Particles(void *al, void *el,
                    const int &Nb_atoms, const int &Nb_atoms_max,
                    const int &Nb_electrons, const int &Nb_electrons_max);

double my_erf(double x);

#endif // INC_POTENTIAL_GENERAL_hpp

// ********** End of file ***************************************
