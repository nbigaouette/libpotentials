#ifndef INC_STRUCTURES_POTENTIALS_HPP
#define INC_STRUCTURES_POTENTIALS_HPP

// Parameter preventing the Coulomb singularity at r = 0
// r / (r+alpha)^3                          [m]
const double sc_alpha = 0.1 * 1.0e-10; // See Mijoule2006, PRA 73, 033203

struct potential_paramaters {

    // ********** Pontential parameters *****************************
    // See "Notes on Potential Shapes" (NPS) and "Notes on Gaussian
    // Particles" (NGP) documents.
    //     Name             Comment                                     [Units]

    // Used by all
    double  dr[3];          // Vector from particle 2 to particle 1     [m]
    double  r;              // Distance between 2 particles             [m]
    double  r2;             // Distance between 2 particles (squared)   [m]
    double  one_over_r;     // Inverse of distance                      [m^-1]
    double  cutoff_radius;  // Minimum radius where Coulombic           [m]

    double  kQ2;            // Charge of particle 2 times k             [V.m]
    double  B;              // Potential right on top of particle 2     [V]
    double  kQ2_over_B;     // kQ2 divided B                            [m]


    // ********** Simple pontential parameters **********************
    // None

    // ********** Harmonic pontential parameters ********************
    double  h_A;            // Amplitude of potential (NPS eq. 6)       [V.m^-2]
    double  h_A_r;          // Amplitude times distance                 [V.m^-1]

    // ********** Super Gaussian parameters *************************
    double  sg_sigma;       // Width of SG (NPS eq. 11)                 [m]
    double  sg_r_over_sigma_two_m;          // (r/sigma)^(2m)           [-]
    double  sg_exp_half_r_over_sigma_two_m; // exp(-0.5(r/sigma)^(2m))  [-]

    // ********** Herman-Skillman (HS) parameters *******************
    int     hs_cs2;         // Charge state of particle 2               [-]

    // ********** Gaussian Distribution parameters ******************
    double  gd_sigma;       // Width of gaussian distribution (NPS eq. 21) [m]

    // ********** Symmetric Charge distribution parameters **********
    double  sym_sigma1;     // Width of gaussian distribution 1         [m]
    double  sym_sigma2;     // Width of gaussian distribution 2         [m]
    int     sym_cs1;        // Charge state of particle 1               [-]
    int     sym_cs2;        // Charge state of particle 2               [-]

    // ********** Screened Coulomb **********************************
    // None
};

#endif // #ifndef INC_STRUCTURES_POTENTIALS_HPP

// ********** End of file ***************************************
