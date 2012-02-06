#ifndef INC_LIBPOTENTIALS_hpp
#define INC_LIBPOTENTIALS_hpp

#include <string>

#include "FloatType.hpp"

#include "Vectors.hpp"
#include "Structure_Potentials.hpp"
#include "Constants.hpp"
#include "Global.hpp"

#include "Potentials_Coulomb.hpp"
#include "Potentials_Symmetric.hpp"
#include "Potentials_GaussianDistribution.hpp"
#include "Potentials_Harmonic.hpp"
#include "Potentials_HermanSkillman.hpp"
#include "Potentials_ScreenedCoulomb.hpp"
#include "Potentials_Simple.hpp"
#include "Potentials_SuperGaussian.hpp"


#ifndef DEBUGP
#define DEBUGP(x)           std_cout << __FILE__ << ":" << __LINE__ << ":\n    " << (x)
#endif // #ifndef DEBUGP

extern bool is_libpotentials_initialized;

void Check_if_LibPotentials_is_initialized(void);
void Potentials_Initialize(const std::string _io_basename,
                           const std::string potential_shape,
                           const fdouble cutoff_base_potential,
                           const fdouble cutoff_radius,
                           const int input_sg_m);
void Potentials_Finalize();

void Set_LibPotentials_as_initialized(void);

void Print_Particles(void *list, const int &N);

bool Is_HS_used();

fdouble erf_over_x(fdouble x);
fdouble erf_over_x3_minus_exp_over_x2(fdouble x);


// **************************************************************
// ********** Function pointers for... **************************
// ...setting the parameters of the potential/field calculation
extern void   (*Potentials_Set_Parameters)(void *p1, void *p2, potential_paramaters &potparams);
// ...calculating the potential
extern fdouble (*Calculate_Potential)(      void *p1, void *p2, potential_paramaters &potparams);
// ...setting the electric field
extern void   (*Set_Field)(                void *p1, void *p2, potential_paramaters &potparams, fdouble &phi, fdouble E[3]);

#endif // INC_LIBPOTENTIALS_hpp

// ********** End of file ***************************************
