#ifndef INC_LIBPOTENTIALS_HS_hpp
#define INC_LIBPOTENTIALS_HS_hpp

#include <vector>

#include <LookUpTable.hpp>

#include "FloatType.hpp"

void Set_HermanSkillman_Lookup_Tables_Xe(std::vector<LookUpTable<fdouble> > &lut_pot,
                                         std::vector<LookUpTable<fdouble> > &lut_field);


#endif // INC_LIBPOTENTIALS_HS_hpp

// ********** End of file ***************************************
