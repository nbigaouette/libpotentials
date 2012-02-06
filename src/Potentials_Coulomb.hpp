#ifndef INC_POTENTIALS_COULOMB_HPP
#define INC_POTENTIALS_COULOMB_HPP

#include "FloatType.hpp"

#include "Structure_Potentials.hpp"

// **************************************************************

fdouble Coulomb_Potential(const fdouble kQ, const fdouble r);
void Set_Coulomb_Field(const fdouble phi12, fdouble E[3],
                       const fdouble dr[3], const fdouble r2);

void Potentials_Set_Parameters_Coulomb(
    void *p1, void *p2,
    potential_paramaters &potparams);
fdouble Calculate_Potential_Cutoff_Coulomb(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_Coulomb(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3]);

#endif // #ifndef INC_POTENTIALS_COULOMB_HPP

// ********** End of file ***************************************
