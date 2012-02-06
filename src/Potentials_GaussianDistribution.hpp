#ifndef INC_POTENTIALS_GAUSSIANDISTRIBUTION_HPP
#define INC_POTENTIALS_GAUSSIANDISTRIBUTION_HPP

#include "FloatType.hpp"

#include "Potentials_GaussianDistribution.hpp"


void Initialize_GaussianDistribution(const fdouble cutoff_base_potential, const fdouble cutoff_radius);

void Potentials_Set_Parameters_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams);
fdouble Calculate_Potential_Cutoff_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3]);


#endif // #ifndef INC_POTENTIALS_GAUSSIANDISTRIBUTION_HPP

// ********** End of file ***************************************
