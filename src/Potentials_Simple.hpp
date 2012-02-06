#ifndef INC_POTENTIALS_SIMPLE_HPP
#define INC_POTENTIALS_SIMPLE_HPP

#include "FloatType.hpp"

#include "Structure_Potentials.hpp"

void Initialize_Simple(const fdouble cutoff_base_potential, const fdouble cutoff_radius);

void Potentials_Set_Parameters_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams);
fdouble Calculate_Potential_Cutoff_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3]);


#endif // #ifndef INC_POTENTIALS_SIMPLE_HPP

// ********** End of file ***************************************
