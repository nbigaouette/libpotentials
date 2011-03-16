/****************************************************************

    File containing global variable declarations. This needs
    to be done to prevent multiple declaration.
    All other files who need these variables will have to
    include "Global.hpp" which declares as extern the variables.

****************************************************************/

#include <Memory.hpp>
#include "Global.hpp"
#include "NR.hpp"

namespace libpotentials_private
{
    fdouble base_pot_well_depth;

//     // Error function lookup table erf(R)
//     fdouble *tl_erf;
//     const int    tl_n = 1000;    // Number of points of the lookup table
//     const fdouble tl_Rmax = 6.0; // Maximum value of R: erf(4) = 0.999999984582742
//     const fdouble tl_dR = tl_Rmax / (tl_n-1); // Step
//     const fdouble tl_one_over_dR = 1.0 / tl_dR;
//
//     // **********************************************************
//     void initialize_erf_lookup_table()
//     {
//         tl_erf = (fdouble *) calloc_and_check(tl_n, sizeof(fdouble));
//         // Populating the lookup table
//         for (int i = 0 ; i < tl_n ; i++)
//         {
//             tl_erf[i] = nr::int_erf(i*tl_dR);
//         }
//     }

    LookUpTable<fdouble> lut_potential;
    LookUpTable<fdouble> lut_field;
}

// ********** End of File ***************************************
