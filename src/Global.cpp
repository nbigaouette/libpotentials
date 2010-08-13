/****************************************************************

    File containing global variable declarations. This needs
    to be done to prevent multiple declaration.
    All other files who need these variables will have to
    include "Global.hpp" which declares as extern the variables.

****************************************************************/

#include "Global.hpp"
#include "Memory.hpp"
#include "NR.hpp"

namespace libpotentials_private
{
    double base_pot_well_depth;

    // Error function lookup table erf(R)
    double *tl_erf;
    const int    tl_n = 200;    // Number of points of the lookup table
    const double tl_Rmax = 4.0; // Maximum value of R: erf(4) = 0.999999984582742
    const double tl_dR = tl_Rmax / (tl_n-1); // Step
    const double tl_one_over_dR = 1.0 / tl_dR;

    // **********************************************************
    void initialize_erf_lookup_table()
    {
        tl_erf = (double *) calloc_and_check(tl_n, sizeof(double));
        // Populating the lookup table
        for (int i = 0 ; i < tl_n ; i++)
        {
            tl_erf[i] = nr::int_erf(i*tl_dR);
        }
    }
}

// ********** End of File ***************************************
