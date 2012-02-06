#ifndef INC_POTENTIALS_HARMONIC_HPP
#define INC_POTENTIALS_HARMONIC_HPP

#include "FloatType.hpp"

#include "Structure_Potentials.hpp"

// ********** Harmonic pontential *******************************
void Initialize_Harmonic(const fdouble cutoff_base_potential, const fdouble cutoff_radius);

void Potentials_Set_Parameters_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams);
fdouble Calculate_Potential_Cutoff_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3]);

#endif // #ifndef INC_POTENTIALS_HARMONIC_HPP

// ********** End of file ***************************************
