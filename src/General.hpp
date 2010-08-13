#ifndef INC_GENERAL_hpp
#define INC_GENERAL_hpp

#include <iostream>
#include <cstring>

#include <cassert>
#include <cmath>
#include <cstdlib> // calloc()
#include <cstdio>

#include "Std_Cout.hpp"

#define DEBUGP(x)           std_cout << __FILE__ << ":" << __LINE__ << ":\n    " << (x)

const int Potential_Shapes_Simple               = 1;
const int Potential_Shapes_Harmonic             = 2;
const int Potential_Shapes_SuperGaussian        = 3;
const int Potential_Shapes_GaussianDistribution = 4;
const int Potential_Shapes_HS_SuperGaussian     = 5;
const int Potential_Shapes_ChargeDistribution_Symmetric = 6;
const int Potential_Shapes_ScreenedCoulomb      = 7;


void Check_if_initialized(void);
void Potentials_Initialize(const int arg_called_from,
                           const std::string potential_shape,
                           const double base_potential_depth);
void Potentials_Finalize();

void Print_Partiles(void *al, void *el,
                    const int &Nb_atoms, const int &Nb_atoms_max,
                    const int &Nb_electrons, const int &Nb_electrons_max);

double my_erf(double x);


#endif // INC_GENERAL_hpp

// ********** End of file ***************************************
