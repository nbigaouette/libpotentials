#ifndef INC_LIBPOTENTIALS_HS_hpp
#define INC_LIBPOTENTIALS_HS_hpp

#include "FloatType.hpp"

std::vector<std::pair<fdouble,fdouble> > Load_HermanSkillman_Xe_1();
std::vector<std::pair<fdouble,fdouble> > Load_HermanSkillman_Xe_2();
std::vector<std::pair<fdouble,fdouble> > Load_HermanSkillman_Xe_3();
std::vector<std::pair<fdouble,fdouble> > Load_HermanSkillman_Xe_4();
std::vector<std::pair<fdouble,fdouble> > Load_HermanSkillman_Xe_5();

inline std::vector<std::pair<fdouble,fdouble> >Load_HermanSkillman_Xe(const int cs)
{
    if      (cs == -1)
        return Load_HermanSkillman_Xe_1();
    else if (cs == 0)
        return Load_HermanSkillman_Xe_1();
    else if (cs == 1)
        return Load_HermanSkillman_Xe_1();
    else if (cs == 2)
        return Load_HermanSkillman_Xe_2();
    else if (cs == 3)
        return Load_HermanSkillman_Xe_3();
    else if (cs == 4)
        return Load_HermanSkillman_Xe_4();
    else if (cs == 5)
        return Load_HermanSkillman_Xe_5();
    else
        abort();
}


#endif // INC_LIBPOTENTIALS_HS_hpp

// ********** End of file ***************************************
