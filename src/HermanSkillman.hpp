#ifndef INC_LIBPOTENTIALS_HS_hpp
#define INC_LIBPOTENTIALS_HS_hpp

#include <vector>

#include <LookUpTable.hpp>

#include "FloatType.hpp"

// Includes neutrals
const int HS_Xe_MaxNbCS = 7;

extern std::vector<fdouble> hs_min_rad;
extern std::vector<LookUpTable<fdouble> > hs_lut_potential;
extern std::vector<LookUpTable<fdouble> > hs_lut_field;

void Initialize_HS_Base_Potential(const fdouble &base_potential_eV);
void Set_HermanSkillman_Lookup_Tables_Xe(std::vector<LookUpTable<fdouble> > &lut_pot,
                                         std::vector<LookUpTable<fdouble> > &lut_field);


#endif // INC_LIBPOTENTIALS_HS_hpp

// ********** End of file ***************************************
