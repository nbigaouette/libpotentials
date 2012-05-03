/*
    This file defines functions used by the library but which
    requires information from the structure.
*/

#include <cstring> // memset()
#include <cstdlib> // abort()
#include <limits> // http://www.cplusplus.com/reference/std/limits/numeric_limits/

#include "MDStructures_Access.hpp"


size_t  Get_Sizeof_particle() { return 0; }

// **************************************************************
// Function definitions to access structure members
// **************************************************************

void    Print(void *b) { }

void    Clear_Particle_Data(void *p)              { }
void    Copy_Particle_Data(void *src, void *dest) { }

void    Print_Particle(void *p) { }

fdouble*Get_E(void *b) { return NULL; }
void    Set_E(void *b, fdouble E[3]) { }

int     Get_Id(void *b) { return -1; }
void    Set_Id(void *b, int id) { }

int     Get_Indx(void *b) { return -1; }
void    Set_Indx(void *b, int i) { }

bool    Get_Is_on(void *b) { return false; }
void    Set_Is_on(void *b, int i) {  }

int     Get_Charge_State(void *b) { return 999999999; }
void    Set_Charge_State(void *b, int cs) { }
fdouble Get_Charge(void *b) { return fdouble(-9999999.9); }
void    Set_Charge(void *b, fdouble c) {  }
fdouble Get_Chargea(void *b) { return fdouble(-9999999.9); }
void    Set_Chargea(void *b, fdouble c) {  }

fdouble Get_Mass(void *b) { return fdouble(-9999999.9); }
void    Set_Mass(void *b, fdouble m) { }

fdouble*Get_Position(void *b) { return NULL; }
void    Set_Position(void *b, fdouble pos[3]) { }

fdouble*Get_Velocity(void *b) { return NULL; }
void    Set_Velocity(void *b, fdouble vel[3]) { }

fdouble Get_Potential(void *b) { return 0.0; }
void    Set_Potential(void *b, fdouble m) { }

fdouble Get_NextIon(void *b) { return 0.0; }
void    Set_NextIon(void *b, fdouble n) { }

int     Get_ClosestIon(void *b) { return 0; }
void    Set_ClosestIon(void *b, int i) { }

fdouble Get_ClosestIon_Distance(void *b) { return 0.0; }
void    Set_ClosestIon_Distance(void *b, fdouble d) { }

int     Get_Nghb_id(void *b) { return 0; }
void    Set_Nghb_id(void *b, int i) { }

int     Get_Nghb_indx(void *b) { return 0; }
void    Set_Nghb_indx(void *b, int i) { }

int     Get_Rcb_mark(void *b) { return 0; }
void    Set_Rcb_mark(void *b, int i) { }

int     Get_Rcb_indx(void *b) { return 0; }
void    Set_Rcb_indx(void *b, int i) { }

int     Get_Imp_mark(void *b) { return 0; }
void    Set_Imp_mark(void *b, int i) { }

int     Get_Coll(void *b) { return 0; }
void    Set_Coll(void *b, int i) { }

int   * Get_orbiting_electrons_Pointer(void *b) { return NULL; }
void    Set_orbiting_electrons_Pointer(void *b, int *p) { }

int     Get_orbiting_electrons_Size() { return 0; }

short int Get_Excited_State_index(void *b)            { return 0; }
void    Set_Excited_State_index(void *b, short int i) { }
void    Set_Excited_State_index(void *b, int i)       { }

// ********** End of File ***************************************
