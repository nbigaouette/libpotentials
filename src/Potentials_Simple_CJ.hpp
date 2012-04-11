#ifndef INC_POTENTIALS_SIMPLE_CJ_HPP
#define INC_POTENTIALS_SIMPLE_CJ_HPP

#include "FloatType.hpp"

#include "Structure_Potentials.hpp"

void Initialize_Simple_CJ(const fdouble cutoff_base_potential, const fdouble cutoff_radius);

void Potentials_Set_Parameters_Simple_CJ(
    void *p1, void *p2,
    potential_paramaters &potparams);
fdouble Calculate_Potential_Cutoff_Simple_CJ(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_Simple_CJ(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3]);


#endif // #ifndef INC_POTENTIALS_SIMPLE_CJ_HPP

// ********** End of file ***************************************
