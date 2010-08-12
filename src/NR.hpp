#ifndef INC_NR_hpp
#define INC_NR_hpp

/*
    Functions taken from
    Numerical Recipes in C
*/

namespace nr
{
double python_erf(double x);
double Gamma(double x);
double Gamma_NaturalLogarithm(double x);
double gammln(double xx);
double gammp(double a, double x);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);
double erff(double x);
double int_erf(double x);
}

#endif // #ifndef INC_NR_hpp

// ********** End of file ***************************************
