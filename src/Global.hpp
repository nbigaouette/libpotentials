#ifndef INC_GLOBAL_hpp
#define INC_GLOBAL_hpp

namespace libpotentials_private
{
    // The base potential well giving the bottom of the potential well.
    extern double base_pot_well_depth;  // [eV]

    extern double *tl_erf;              // Error function lookup table
    extern const int    tl_n;           // Number of points of the lookup table
    extern const double tl_Rmax;        // Maximum value of R: erf(4) = 0.999999984582742
    extern const double tl_dR;          // Step
    extern const double tl_one_over_dR;

    void initialize_erf_lookup_table();
}

#endif // INC_GLOBAL_hpp

// ********** End of File ***************************************
