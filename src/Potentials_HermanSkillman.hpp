#ifndef INC_POTENTIALS_HS_HPP
#define INC_POTENTIALS_HS_HPP

#include <vector>

#include "FloatType.hpp"

#include "Structure_Potentials.hpp"

// Includes neutrals
const int HS_Xe_MaxNbCS = 7;

extern std::vector<fdouble> hs_min_rad;
extern std::vector<LookUpTable<fdouble> > hs_lut_potential;
extern std::vector<LookUpTable<fdouble> > hs_lut_field;

bool Is_HS_used();
void Initialize_HermanSkillman(const fdouble cutoff_base_potential, const fdouble cutoff_radius);
void Initialize_HS_Base_Potential(const fdouble &base_potential_eV);
void Initialize_HS_Cutoff_Radius(const fdouble &cutoff_radius_m, const bool scale_radius_to_symmetric = true);
void Set_HermanSkillman_Lookup_Tables_Xe(std::vector<LookUpTable<fdouble> > &lut_pot,
                                         std::vector<LookUpTable<fdouble> > &lut_field);

void Potentials_Set_Parameters_HS(
    void *p1, void *p2,
    potential_paramaters &potparams);
fdouble Calculate_Potential_Cutoff_HS(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_HS(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3]);

#endif // #ifndef INC_POTENTIALS_HS_HPP

// ********** End of file ***************************************
