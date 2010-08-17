#ifndef INC_LIBPOTENTIAL_hpp
#define INC_LIBPOTENTIAL_hpp

#include <string>

#include "Potentials.hpp"
#include "Memory.hpp"
#include "Vectors.hpp"
#include "Structure_Potentials.hpp"
#include "Std_Cout.hpp"
#include "Version.hpp"
#include "Constants.hpp"
// Public interface of potentials library

#ifndef DEBUGP
#define DEBUGP(x)           std_cout << __FILE__ << ":" << __LINE__ << ":\n    " << (x)
#endif // #ifndef DEBUGP

void Potentials_Initialize(const std::string potential_shape,
                           const double base_potential_depth,
                           const double input_s_rmin,
                           const double input_sg_m);
void Potentials_Finalize();

void Set_LibPotentials_as_initialized(void);

void Print_Particles(void *list, const int &N);

#endif // INC_LIBPOTENTIAL_hpp

// ********** End of file ***************************************
