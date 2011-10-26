#ifndef INC_LIBPOTENTIALS_hpp
#define INC_LIBPOTENTIALS_hpp

#include <string>

#include "FloatType.hpp"

#include "Potentials.hpp"
#include "Vectors.hpp"
#include "Structure_Potentials.hpp"
#include "Version.hpp"
#include "Constants.hpp"

// Public interface of potentials library

#ifndef DEBUGP
#define DEBUGP(x)           std_cout << __FILE__ << ":" << __LINE__ << ":\n    " << (x)
#endif // #ifndef DEBUGP

void Potentials_Initialize(const std::string _io_basename,
                           const std::string potential_shape,
                           const fdouble base_potential_depth,
                           const fdouble input_s_rmin,
                           const int input_sg_m);
void Potentials_Finalize();

void Set_LibPotentials_as_initialized(void);

void Print_Particles(void *list, const int &N);

bool Is_HS_used();

fdouble erf_over_x(fdouble x);
fdouble erf_over_x3_minus_exp_over_x2(fdouble x);

#endif // INC_LIBPOTENTIALS_hpp

// ********** End of file ***************************************
