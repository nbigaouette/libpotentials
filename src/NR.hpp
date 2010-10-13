#ifndef INC_NR_hpp
#define INC_NR_hpp

#include "FloatType.hpp"

/*
    Functions taken from
    Numerical Recipes in C
*/

namespace nr
{
fdouble python_erf(fdouble x);
fdouble Gamma(fdouble x);
fdouble Gamma_NaturalLogarithm(fdouble x);
fdouble gammln(fdouble xx);
fdouble gammp(fdouble a, fdouble x);
void gcf(fdouble *gammcf, fdouble a, fdouble x, fdouble *gln);
void gser(fdouble *gamser, fdouble a, fdouble x, fdouble *gln);
void gcf(fdouble *gammcf, fdouble a, fdouble x, fdouble *gln);
fdouble erff(fdouble x);
fdouble int_erf(fdouble x);
}

#endif // #ifndef INC_NR_hpp

// ********** End of file ***************************************
