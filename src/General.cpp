/****************************************************************

    General functions

****************************************************************/

#include <cstring>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <cfloat>
#include <climits> // CHAR_BIT
#include <stdint.h> // (u)int64_t

#include "Constants.hpp"
#include "General.hpp"
#include "Potentials.hpp"
#include "NR.hpp"
#include "Global.hpp"
#include "Code_Functions_Declarations.hpp"
#include "Std_Cout.hpp"
#include "Version.hpp"
#include "Memory.hpp"


using namespace libpotentials;
//


// **************************************************************
// ********** Accessible functions implementations **************
// **************************************************************


// **************************************************************
const char* Potentials_Get_Git_Commit()
{
    return build_sha;
}

// ********** End of file ***************************************

