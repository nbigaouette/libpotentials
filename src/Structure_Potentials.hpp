#ifndef INC_STRUCTURES_POTENTIALS_HPP
#define INC_STRUCTURES_POTENTIALS_HPP

#include "FloatType.hpp"

// Parameter preventing the Coulomb singularity at r = 0
// r / (r+alpha)^3                          [m]
const fdouble sc_alpha = fdouble(0.1 * 1.0e-10); // See Mijoule2006, PRA 73, 033203

struct potential_paramaters {

    // ********** Pontential parameters *****************************
    // See "Notes on Potential Shapes" (NPS) and "Notes on Gaussian
    // Particles" (NGP) documents.
    //     Name             Comment                                     [Units]

    // Used by all
    fdouble  dr[3];          // Vector from particle 2 to particle 1     [m]
    fdouble  r;              // Distance between 2 particles             [m]
    fdouble  r2;             // Distance between 2 particles (squared)   [m]
    fdouble  one_over_r;     // Inverse of distance                      [m^-1]

    fdouble  kQ2;            // Charge of particle 2 times k             [V.m]
    fdouble  B;              // Potential right on top of particle 2     [V]
    fdouble  kQ2_over_B;     // kQ2 divided B                            [m]


    // ********** Simple pontential parameters **********************
    // None

    // ********** Harmonic pontential parameters ********************
    fdouble  h_A;            // Amplitude of potential (NPS eq. 6)       [V.m^-2]
    fdouble  h_A_r;          // Amplitude times distance                 [V.m^-1]

    // ********** Super Gaussian parameters *************************
    fdouble  sg_sigma;       // Width of SG (NPS eq. 11)                 [m]
    fdouble  sg_r_over_sigma_two_m;          // (r/sigma)^(2m)           [-]
    fdouble  sg_exp_half_r_over_sigma_two_m; // exp(-0.5(r/sigma)^(2m))  [-]

    // ********** Herman-Skillman (HS) parameters *******************
    int     hs_cs1;         // Charge state of particle 1               [-]
    int     hs_cs2;         // Charge state of particle 2               [-]
    unsigned int hs_lut_i;  // Lookup table index

    // ********** Gaussian Distribution parameters ******************
    fdouble  gd_sigma;       // Width of gaussian distribution (NPS eq. 21) [m]

    // ********** Symmetric Charge distribution parameters **********
    int     sym_cs1;        // Charge state of particle 1               [-]
    int     sym_cs2;        // Charge state of particle 2               [-]

    // ********** Screened Coulomb **********************************
    // None

    // Constructor: Set everything to 0
    potential_paramaters()
    {
        dr[0] = 0.0; dr[1] = 0.0; dr[2] = 0.0;
        r = 0.0;
        r2 = 0.0;
        one_over_r = 0.0;
        kQ2 = 0.0;
        B = 0.0;
        kQ2_over_B = 0.0;
        h_A = 0.0;
        h_A_r = 0.0;
        sg_sigma = 0.0;
        sg_r_over_sigma_two_m = 0.0;
        sg_exp_half_r_over_sigma_two_m = 0.0;
        hs_cs1 = 0;
        hs_cs2 = 0;
        hs_lut_i = 0;
        gd_sigma = 0.0;
        sym_cs1 = 0;
        sym_cs2 = 0;
    }
};

#endif // #ifndef INC_STRUCTURES_POTENTIALS_HPP

// ********** End of file ***************************************
