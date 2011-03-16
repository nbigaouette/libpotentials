#ifndef INC_GLOBAL_hpp
#define INC_GLOBAL_hpp

#include "FloatType.hpp"
#include "LookUpTable.hpp"

namespace libpotentials_private
{
    // The base potential well giving the bottom of the potential well.
    extern fdouble base_pot_well_depth; // [eV]

//     extern fdouble *tl_erf;             // Error function lookup table
//     extern const int    tl_n;           // Number of points of the lookup table
//     extern const fdouble tl_Rmax;       // Maximum value of R: erf(4) = 0.999999984582742
//     extern const fdouble tl_dR;         // Step
//     extern const fdouble tl_one_over_dR;
//
//     void initialize_erf_lookup_table();

    extern LookUpTable<fdouble> lut_potential;
    extern LookUpTable<fdouble> lut_field;
}

#endif // INC_GLOBAL_hpp

// ********** End of File ***************************************
