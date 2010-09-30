#ifndef INC_POTENTIALS_HPP
#define INC_POTENTIALS_HPP

#include "Structure_Potentials.hpp"


extern bool is_libpotentials_initialized;

void Check_if_LibPotentials_is_initialized(void);

double Coulomb_Potential(const double kQ, const double r);
void   Coulomb_Field(const double phi12, double E[3], const double dr[3], const double dr2);

double LibPotentialErf(double x);

// **************************************************************
// ********** Function pointers for... **************************

// ...setting the parameters of the potential/field calculation
extern void   (*Potentials_Set_Parameters)(void *p1, void *p2, potential_paramaters &potparams);
// ...calculating the potential
extern double (*Calculate_Potential)(      void *p1, void *p2, potential_paramaters &potparams);
// ...setting the electric field
extern void   (*Set_Field)(                void *p1, void *p2, potential_paramaters &potparams, double &phi, double E[3]);


// **************************************************************
// ********** Real functions ************************************

// ********** Simple pontential *********************************
void Potentials_Set_Parameters_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Harmonic pontential *******************************
void Potentials_Set_Parameters_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Super Gaussian ************************************
void Potentials_Set_Parameters_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Herman-Skillman (HS) ******************************
void Potentials_Set_Parameters_HS_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_HS_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_HS_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Gaussian Distribution parameters ******************
void Potentials_Set_Parameters_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Symmetric Charge distribution *********************
void Potentials_Set_Parameters_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Screened Coulomb **********************************
void Potentials_Set_Parameters_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Pseudo-particles (QFD) ****************************
void Potentials_Set_Parameters_PseudoParticles(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_PseudoParticles(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_PseudoParticles(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);


#endif // #ifndef INC_POTENTIALS_HPP

// ********** End of file ***************************************
