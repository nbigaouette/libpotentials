#ifndef INC_MDSTRUCTURES_ACCESS_HPP
#define INC_MDSTRUCTURES_ACCESS_HPP

#include "limits.h"

#include <FloatType.hpp>
#include <Assert.hpp>

// **************************************************************
// Function definitions to access structure members
// **************************************************************

void    Print(void *b);

void    Copy_Particle_Data(void *src, void *dest);
void    Clear_Particle_Data(void *p);

fdouble* Get_E(void *b);
void    Set_E(void *b, fdouble E[3]);

int     Get_Id(void *b);
void    Set_Id(void *b, int id);

int     Get_Indx(void *b);
void    Set_Indx(void *b, int i);

bool    Get_Is_on(void *b);
void    Set_Is_on(void *b, int i);

int     Get_Charge_State(void *b);
void    Set_Charge_State(void *b, int cs);
fdouble  Get_Charge(void *b);
void    Set_Charge(void *b, fdouble c);
fdouble  Get_Chargea(void *b);
void    Set_Chargea(void *b, fdouble c);

fdouble  Get_Mass(void *b);
void    Set_Mass(void *b, fdouble m);

fdouble* Get_Position(void *b);
void    Set_Position(void *b, fdouble pos[3]);

fdouble* Get_Velocity(void *b);
void    Set_Velocity(void *b, fdouble vel[3]);

fdouble  Get_Potential(void *b);
void    Set_Potential(void *b, fdouble m);

fdouble  Get_NextIon(void *b);
void    Set_NextIon(void *b, fdouble n);

int     Get_ClosestIon(void *b);
void    Set_ClosestIon(void *b, int i);

fdouble  Get_ClosestIon_Distance(void *b);
void    Set_ClosestIon_Distance(void *b, fdouble d);

int     Get_Nghb_id(void *b);
void    Set_Nghb_id(void *b, int i);

int     Get_Nghb_indx(void *b);
void    Set_Nghb_indx(void *b, int i);

int     Get_Rcb_mark(void *b);
void    Set_Rcb_mark(void *b, int i);

int     Get_Rcb_indx(void *b);
void    Set_Rcb_indx(void *b, int i);

int     Get_Imp_mark(void *b);
void    Set_Imp_mark(void *b, int i);

int     Get_Coll(void *b);
void    Set_Coll(void *b, int i);

int   * Get_orbiting_electrons_Pointer(void *b);
void    Set_orbiting_electrons_Pointer(void *b, int *p);

int     Get_orbiting_electrons_Size();

short int Get_Excited_State_index(void *b);
void    Set_Excited_State_index(void *b, short int i);
void    Set_Excited_State_index(void *b, int i);

#endif // INC_MDSTRUCTURES_ACCESS_HPP

// ********** End of file ***************************************
