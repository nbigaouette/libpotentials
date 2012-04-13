#ifndef _CODE_FUNCTIONS_DECLARATIONS_HPP
#define _CODE_FUNCTIONS_DECLARATIONS_HPP

/*
    This file contain declaration of functions that should
    be defined in the actual code, not the library, since they
    are code dependant.
*/

#include "FloatType.hpp"

extern size_t  Get_Sizeof_particle();

static inline void * get_voidp(void *p, int i)
{
    return (void *) (((unsigned char*) p)+(i*Get_Sizeof_particle()));
}

// **************************************************************
// Obtained running:
//      sed -n '/.*my_pattern1.*/,/.*my_pattern2.*/p' treestructures.h | sed -e "s|inline ||g" -e "s| {.*}|;|g"
// in same directory as "treestructures.h"

extern char *  Get_Name(void *b);
extern void    Print_name(void *b);

extern int     Get_Id(void *b);
extern void    Set_Id(void *b, int id);

extern int     Get_Is_on(void *b);
extern void    Set_Is_on(void *b, int i);

extern int     Get_Charge_State(void *b);
extern void    Set_Charge_State(void *b, int cs);
extern fdouble Get_Charge(void *b);
extern void    Set_Charge(void *b, fdouble c);
extern fdouble Get_Chargea(void *b);
extern void    Set_Chargea(void *b, fdouble c);

extern fdouble Get_Mass(void *b);
extern void    Set_Mass(void *b, fdouble m);

extern fdouble*Get_Position(void *b);
extern void    Set_Position(void *b, fdouble p[3]);
extern fdouble*Get_Velocity(void *b);
extern void    Set_Velocity(void *b, fdouble p[3]);

extern fdouble Get_Potential(void *b);
extern void    Set_Potential(void *b, fdouble m);

extern fdouble*Get_E(void *b);
extern void    Set_E(void *b, fdouble E[3]);

extern int     Get_NextIon(void *b);
extern void    Set_NextIon(void *b, fdouble n);

extern fdouble Get_ClosestIon_Distance(void *b);
extern void    Set_ClosestIon_Distance(void *b, fdouble d);

extern int     Get_Nghb_id(void *b);
extern void    Set_Nghb_id(void *b, int i);

extern int     Get_Rcb_mark(void *b);
extern void    Set_Rcb_mark(void *b, int i);

extern int     Get_Rcb_indx(void *b);
extern void    Set_Rcb_indx(void *b, int i);

extern int     Get_Imp_mark(void *b);
extern void    Set_Imp_mark(void *b, int i);

extern char *  Get_Name(void *b);
extern void    Print_name(void *b);

extern fdouble*Get_Efld(void *b);

// **************************************************************

#endif // #ifndef _CODE_FUNCTIONS_DECLARATIONS_HPP

// ********** End of file ***************************************
