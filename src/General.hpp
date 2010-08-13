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

// **************************************************************
template <class Integer>
static inline void * calloc_and_check(Integer nb, size_t s, const std::string msg = "")
{
    void *p = NULL;
    p = calloc(nb, s);
    if (p == NULL)
    {
        double nb_s = double(nb) * double(s);
        std_cout << "Allocation of " << nb << " x " << s << " bytes = " << nb_s << " bytes failed" << std::endl;
        std_cout << "(" << double(nb_s) / (1024.0) << " KiB, "
                         << double(nb_s) / (1024.0*1024.0) << " MiB, "
                         << double(nb_s) / (1024.0*1024.0*1024.0) << " GiB)" << std::endl;
        std_cout << "p = " << p << "\n";
        if (msg != "")
        {
            std_cout << "Comment: " << msg << std::endl;
        }
        std_cout << "Aborting\n" << std::flush;
        abort();
    }
    return p;
}

// **************************************************************
template <class Integer>
static inline void * malloc_and_check(Integer nb, size_t s, const std::string msg = "")
{
    void *p = NULL;
    p = malloc(nb * s);
    if (p == NULL)
    {
        const long unsigned int slu  = (long unsigned int) s;
        const long unsigned int slun = (long unsigned int) s * nb;
        std_cout << "Allocation of " << nb << " x " << slu << " bytes = " << slun << " bytes failed" << std::endl;
        std_cout << "(" << double(slun) / (1024.0) << " KiB, "
                         << double(slun) / (1024.0*1024.0) << " MiB, "
                         << double(slun) / (1024.0*1024.0*1024.0) << " GiB)" << std::endl;
        std_cout << "p = " << p << "\n";
        if (msg != "")
        {
            std_cout << "Comment: " << msg << std::endl;
        }
        std_cout << "Aborting\n" << std::flush;
        abort();
    }
    return p;
}

// **************************************************************
template <class Pointer>
void free_me(Pointer &p)
{
    if (p != NULL)
    {
        free(p);
    }
    p = NULL;
}

#endif // INC_GENERAL_hpp

// ********** End of file ***************************************
