/***************************************************************
 *
 *
 ***************************************************************/

#include <cstdlib>
#include <iostream>
#include <vector>
#include <stdint.h> // uint64_t

#include <LibPotentials.hpp>
#include "FloatType.hpp"

class Particle
{
    public:
    int id;
    int charge_state;
    fdouble mass;
    fdouble vel[3];
    fdouble pos[3];
    fdouble E[3];
    fdouble potential;
    fdouble charge() { return fdouble(charge_state)*libpotentials::e0; }
};

int         Get_Id(void *b)         { return ((Particle *)b)->id; }
fdouble     Get_Charge(void *b)     { return ((Particle *)b)->charge(); }
fdouble     Get_Mass(void *b)       { return ((Particle *)b)->mass; }
fdouble*    Get_Velocity(void *b)   { return ((Particle *)b)->vel; }
fdouble     Get_Potential(void *b)  { return ((Particle *)b)->potential; }
fdouble*    Get_Position(void *b)   { return ((Particle *)b)->pos; }
fdouble*    Get_E(void *b)          { return ((Particle *)b)->E; }
int         Get_Charge_State(void *b) { return ((Particle *)b)->charge_state; }
size_t      Get_Sizeof_particle()   { return sizeof(Particle); }

// **************************************************************
int main(int argc, char *argv[])
{
    //Potentials_Initialize("HermanSkillman",
    Potentials_Initialize("Symmetric",
                            fdouble(1.0 * libpotentials::Eh_to_eV),     // base potential
                            fdouble(0.5 * libpotentials::bohr_to_m),    // Simple cutoff radius
                            1);                                 // Super Gaussian order (m=1 for gaussian)

    Particle p0;
    Particle p1;
    potential_paramaters potparams;

    fdouble E_at_p0_from_p1[3]  = {libpotentials::zero, libpotentials::zero, libpotentials::zero};
    fdouble potential_at_p0_from_p1;

    for (int d = 0 ; d < 3 ; d++)
    {
        p0.pos[d] = libpotentials::zero;
        p1.pos[d] = libpotentials::zero;
    }
    p0.charge_state = libpotentials::one;
    p1.charge_state = libpotentials::one;


    const int N = 1000;
    const fdouble xmin =  0.001 * libpotentials::bohr_to_m;
    const fdouble xmax = 20.000 * libpotentials::bohr_to_m;
    const fdouble dx = (xmax - xmin) / fdouble(N);

    fdouble r;
    for (int i = 0 ; i < N ; i++)
    {
        r = fdouble(i) * dx + xmin;
        p1.pos[0] = r;
        for (int d = 0 ; d < 3 ; d++)
        {
            E_at_p0_from_p1[d]  = libpotentials::zero;
        }
        potential_at_p0_from_p1 = libpotentials::zero;

        // Set parameters, calculate potential and field
        Potentials_Set_Parameters((void *) &p0, (void *) &p1, potparams);
        potential_at_p0_from_p1 = Calculate_Potential((void *) &p0, (void *) &p1, potparams);
        Set_Field((void *) &p0, (void *) &p1, potparams, potential_at_p0_from_p1, E_at_p0_from_p1);

        //std::cerr << r * libpotentials::m_to_bohr << "  " << potential_at_p0_from_p1 * libpotentials::si_to_au_pot << "\n";
        std::cerr << r * libpotentials::m_to_bohr << "  " << E_at_p0_from_p1[0] * libpotentials::si_to_au_field << "\n";
    }


    return EXIT_SUCCESS;
}
