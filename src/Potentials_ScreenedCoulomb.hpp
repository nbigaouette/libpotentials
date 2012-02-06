#ifndef INC_POTENTIALS_SCREENEDCOULOMB_HPP
#define INC_POTENTIALS_SCREENEDCOULOMB_HPP

#include "FloatType.hpp"

#include "Structure_Potentials.hpp"

// ********** Screened Coulomb **********************************
void Potentials_Set_Parameters_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams);
fdouble Calculate_Potential_Cutoff_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3]);

#endif // #ifndef INC_POTENTIALS_SCREENEDCOULOMB_HPP

// ********** End of file ***************************************
