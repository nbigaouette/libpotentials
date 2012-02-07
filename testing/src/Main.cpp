/***************************************************************
 *
 *
 ***************************************************************/

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <vector>
#include <stdint.h> // uint64_t
#include <stdio.h>
#include <fstream>


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

    std::vector<std::string> potential_shapes;
    potential_shapes.push_back("HermanSkillman");
    potential_shapes.push_back("Symmetric");
    potential_shapes.push_back("Simple");
    potential_shapes.push_back("GaussianDistribution");
    potential_shapes.push_back("SuperGaussian");
    //potential_shapes.push_back("PureCoulomb");

    for (unsigned int pi = 0 ; pi < potential_shapes.size() ; pi++)
    {
        const std::string potential_shape = potential_shapes[pi];

        const std::string cmd = std::string("mkdir -p output/") + potential_shape;
        system(cmd.c_str());

        fdouble cutoff_base_potential = -1.5 * libpotentials::Eh_to_eV;
        fdouble cutoff_radius         = +1.0 * libpotentials::bohr_to_m;
        Potentials_Initialize("output",
                                potential_shape,
                                cutoff_base_potential, // base potential (negative to ignore) [eV]
                                cutoff_radius, // Cutoff radius (negative to ignore) [m]
                                1);  // Super Gaussian order (m=1 for gaussian)

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

        const int N = 10000;
        const fdouble xmin = fdouble(0.001 * libpotentials::bohr_to_m);
        const fdouble xmax = fdouble(5.00 * libpotentials::bohr_to_m);
        const fdouble dx = (xmax - xmin) / fdouble(N);

        fdouble r;

        std::vector<int> charge_states;
        charge_states.push_back(0);
        charge_states.push_back(1);
        charge_states.push_back(2);
        charge_states.push_back(3);
        //charge_states.push_back(4);
        //charge_states.push_back(5);
        //charge_states.push_back(6);
        //charge_states.push_back(7);
        //charge_states.push_back(8);
        //charge_states.push_back(9);
        //charge_states.push_back(10);
        //charge_states.push_back(11);
        //charge_states.push_back(12);

        for (unsigned int csi = 0 ; csi < charge_states.size() ; csi++)
        {
            const int cs = charge_states[csi];

            //p0.charge_state = 1;    // p0 is a 1+ feeling the ion p1
            p0.charge_state = -1;   // p0 is an electron feeling the ion p1
            p1.charge_state = cs;

            char filename[1024];
            if (cs == -1)
                sprintf(filename, "output/%s/field_-%1d.csv", potential_shape.c_str(), std::abs(cs));
            else
                sprintf(filename, "output/%s/field_%02d.csv", potential_shape.c_str(), cs);
            std::ofstream f_field(filename);
            if (cs == -1)
                sprintf(filename, "output/%s/poten_-%1d.csv", potential_shape.c_str(), std::abs(cs));
            else
                sprintf(filename, "output/%s/poten_%02d.csv", potential_shape.c_str(), cs);
            std::ofstream f_poten(filename);

            assert(f_poten.is_open());
            assert(f_field.is_open());

            for (int i = 0 ; i < N ; i++)
            {
                r = fdouble(i) * dx + xmin;

                // Put first particle at right of second. The particle creating the field is at the left,
                // so a positively charged p1 will create a positive field.
                p0.pos[0] = r;

                // Reset values
                potential_at_p0_from_p1 = libpotentials::zero;
                for (int d = 0 ; d < 3 ; d++)
                    E_at_p0_from_p1[d]  = libpotentials::zero;

                // Set parameters, calculate potential and field
                Potentials_Set_Parameters((void *) &p0, (void *) &p1, potparams);
                potential_at_p0_from_p1 = Calculate_Potential((void *) &p0, (void *) &p1, potparams);
                Set_Field((void *) &p0, (void *) &p1, potparams, potential_at_p0_from_p1, E_at_p0_from_p1);

                // Save potential ENERGY and field
                f_poten << r * libpotentials::m_to_bohr << ", " << fdouble(p0.charge_state) * potential_at_p0_from_p1 * libpotentials::si_to_au_pot << "\n";
                f_field << r * libpotentials::m_to_bohr << ", " << E_at_p0_from_p1[0] * libpotentials::si_to_au_field << "\n";
            }

            f_poten.close();
            f_field.close();
        }
    }

    return EXIT_SUCCESS;
}

