#ifndef INC_POTENTIALS_CHARGEDISTRIBUTION_SYMMETRIC_HPP
#define INC_POTENTIALS_CHARGEDISTRIBUTION_SYMMETRIC_HPP

#include "FloatType.hpp"

#include "Structure_Potentials.hpp"


namespace Potentials_Symmetric
{
    extern fdouble sigma;
}

void Initialize_Symmetric(const fdouble cutoff_base_potential, const fdouble cutoff_radius);

void Potentials_Set_Parameters_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams);
fdouble Calculate_Potential_Cutoff_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3]);


#endif // #ifndef INC_POTENTIALS_CHARGEDISTRIBUTION_SYMMETRIC_HPP

// ********** End of file ***************************************
