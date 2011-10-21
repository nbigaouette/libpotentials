#ifndef INC_LIBPOTENTIALS_HS_hpp
#define INC_LIBPOTENTIALS_HS_hpp

#include <vector>

#include <LookUpTable.hpp>

#include "FloatType.hpp"

// Includes electrons (0) and neutral (1).
const int HS_Xe_MaxNbCS = 8;

extern std::vector<fdouble> hs_min_rad;
extern std::vector<LookUpTable<fdouble> > hs_lut_potential;
extern std::vector<LookUpTable<fdouble> > hs_lut_field;


void Set_HermanSkillman_Lookup_Tables_Xe(std::vector<LookUpTable<fdouble> > &lut_pot,
                                         std::vector<LookUpTable<fdouble> > &lut_field);


#endif // INC_LIBPOTENTIALS_HS_hpp

// ********** End of file ***************************************
