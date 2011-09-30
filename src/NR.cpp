#ifdef __SUNPRO_CC
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>  // memset
#include <assert.h>
#else // #ifdef __SUNPRO_CC
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <iostream> // cout
#include <cstring>  // memset
#include <cassert>
#endif // #ifdef __SUNPRO_CC


#include <StdCout.hpp>
#include "Constants.hpp"
#include "NR.hpp"

namespace nr
{

// **************************************************************
fdouble Gamma_NaturalLogarithm(fdouble xx)
/**
 * Natural logarithm of the (approximate) gamma function. See:
 * See "Numerical Recipes in C", 2nd edition, page 214
 * Gives the same result as gammln() from the treecode (2008 09 25)
 */
{
    int i;
    fdouble cof[6],stp;
    fdouble x,y,tmp,ser;

    cof[0] = 76.18009172947146;
    cof[1] = -86.50532032941677;
    cof[2] = 24.01409824083091;
    cof[3] = -1.231739572450155;
    cof[4] = 0.1208650973866179E-2;
    cof[5] = -0.5395239384953E-5;

    stp = 2.5066282746310005;

    x = xx;
    y = x;
    tmp = x + 5.5;
    tmp = (x+fhalf)*std::log(tmp) - tmp;

    ser = 1.000000000190015;

    for( i = 0; i < 6; i++ )
    {
        y += 1;
        ser += cof[i]/y;
    }

    return (tmp + std::log(stp*ser/x) );
}

// **************************************************************
fdouble Gamma(fdouble xx)
/**
 * Approximate gamma function.
 * See "Numerical Recipes in C", 2nd edition, page 214
 */
{
    return std::exp( Gamma_NaturalLogarithm(xx) );
}

// **************************************************************
fdouble gammln(fdouble xx)
/**
 * Returns the value ln[Γ(xx)] for xx > 0.
 * See "Numerical Recipes in C", 2nd edition, page 214
 */
{
    // Internal arithmetic will be done in fdouble precision,
    // a nicety that you can omit if ﬁve-ﬁgure accuracy is good enough.
    fdouble x, y, tmp, ser;
    static fdouble cof[6] = {
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5
    };

    int j;
    y = x = xx;
    tmp = x + 5.5;
    ser = 1.000000000190015;
    tmp -= (x + fhalf) * std::log(tmp);
    for (j = 0 ; j <= 5 ; j++) ser += cof[j] / ++y;

    return -tmp + std::log(fdouble(2.5066282746310005) * ser / x);
}

// **************************************************************
fdouble gammp(fdouble a, fdouble x)
/**
 * Returns the incomplete gamma function P (a, x).
 * See "Numerical Recipes in C", 2nd edition, page 218
 */
{
    fdouble gamser, gammcf, gln;

    if (x < 0.0 || a <= 0.0)
    {
        std_cout << "Invalid arguments in routine gammp\n";
        abort();
    }

    if (x < (a+1.0))
    {
        //Use the series representation.
        gser(&gamser, a, x, &gln);
        return gamser;
    } else {
        // Use the continued fraction representation and take its complement.
        gcf(&gammcf, a, x, &gln);
        return 1.0 - gammcf;
    }
}

// **************************************************************
void gser(fdouble *gamser, fdouble a, fdouble x, fdouble *gln)
/**
 * Returns the incomplete gamma function P (a, x) evaluated by
 * its series representation as gamser.
 * Also returns ln Γ(a) as gln.
 * See "Numerical Recipes in C", 2nd edition, page 218
 */
{
    const int    ITMAX = 100;
    const fdouble EPS   = 3.0e-7;

    int n;
    fdouble sum, del, ap;
    *gln = gammln(a);

    if (x <= 0.0)
    {
        if (x < 0.0)
        {
            std_cout << "x less than 0 in routine gser\n";
            abort();
        }

        *gamser=0.0;
        return;
    } else {
        ap = a;
        del = sum = 1.0 / a;
        for (n = 1 ; n <= ITMAX ; n++)
        {
            ++ap;
            del *= x / ap;
            sum += del;
            if (std::abs(del) < std::abs(sum)*EPS)
            {
                *gamser = sum * std::exp(-x + a * std::log(x) - (*gln));
                return;
            }
        }
        std_cout << "a too large, ITMAX too small in routine gser\n";
        abort();
    }
}

// **************************************************************
void gcf(fdouble *gammcf, fdouble a, fdouble x, fdouble *gln)
/**
 * Returns the incomplete gamma function Q(a, x) evaluated by
 * its continued fraction representation as gammcf. Also
 * returns ln Γ(a) as gln.
 * See "Numerical Recipes in C", 2nd edition, page 219
 */
{
    int i;
    fdouble an, b, c, d, del, h;

    const int    ITMAX = 100;       // Maximum allowed number of iterations.
    const fdouble EPS   = 3.0e-7;    // Relative accuracy.
    const fdouble FPMIN = 1.0e-30;   // Number near the smallest representable ﬂoating-point number.

    *gln = gammln(a);

    // Set up for evaluating continued fraction by modiﬁed
    // Lentz's method (§5.2) with b0 = 0.
    b = x + 1.0 - a;
    c=1.0 / FPMIN;
    d=1.0 / b;

    h = d;
    for (i = 1 ; i <= ITMAX ; i++)
    {
        // Iterate to convergence.
        an = -i*(i-a);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (std::abs(del-fone) < EPS) break;
    }
    if (i > ITMAX)
    {
        std_cout << "a too large, ITMAX too small in gcf\n";
        abort();
    }

    // Put factors in front.
    *gammcf = std::exp(-x + a * std::log(x) - (*gln)) * h;
}

// **************************************************************
fdouble erff(fdouble x)
/**
 * Returns the error function erf(x).
 * See "Numerical Recipes in C", 2nd edition, page 220
 */
{
    return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}

// **************************************************************
fdouble int_erf(fdouble x)
/**
 * "Pretty accurate simpson integration for the error function"
 * Taken from MDCluster's mdcluster_dyn.c line 27
 */
{
    int loop,n=1000; /* number of points */
    fdouble I,h;
    h=(fdouble)x/(n-1);
    I=0;
    I+=fdouble(3.0/8.0)  *h;
    I+=fdouble(7.0/6.0)  *h*std::exp(-h*h);
    I+=fdouble(23.0/24.0)*h*std::exp(-fdouble(4.0)*h*h);
    for (loop=3;loop<n-3;loop++) I+=h*std::exp(-fdouble(loop*loop)*h*h);
    I+=fdouble(23.0/24.0)*h*std::exp(-fdouble(loop*loop)*h*h);   loop++;
    I+=fdouble(7.0/6.0)*h*std::exp(-fdouble(loop*loop)*h*h);   loop++;
    I+=fdouble(3.0/8.0)*h*std::exp(-fdouble(loop*loop)*h*h);
    I*=fdouble(2.0)/std::sqrt(libpotentials::Pi);
    return(I);
}

// **************************************************************
fdouble python_erf(fdouble x)
/**
 * Taken from http://www.johndcook.com/blog/2009/01/19/stand-alone-error-function-erf/
 */
{
    // constants
    const fdouble a1 =  0.254829592;
    const fdouble a2 = -0.284496736;
    const fdouble a3 =  1.421413741;
    const fdouble a4 = -1.453152027;
    const fdouble a5 =  1.061405429;
    const fdouble p  =  0.3275911;

    // Save the sign of x
    fdouble sign = 1.0;
    if (x < 0.0) sign = -1.0;
    const fdouble absx = std::abs(x);

    // A & S 7.1.26
    const fdouble t = 1.0 / (1.0 + p * absx);
    const fdouble y = fone - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t * std::exp(-absx*absx);

    return sign*y;
}

}

// ********** End of file ***************************************
