#ifndef _ELEMENTSTRUCTURE_HPP
#define _ELEMENTSTRUCTURE_HPP

#define MAX_IPS               55

#include <cstdlib> // free()

#include "FloatType.hpp"

struct elem_props
{
    char name[20];
    char symbol[20];
    int atomicNumber;
    fdouble mass;            // kg
    fdouble next_nghb;       // m
    fdouble next_nghb_05;    // m
    int nb_avail_ips;       // Number of available IPs
    int max_Z;              // Max charge state
    int max_n;              // Max shell (n) considered
    int max_l;              // Max subshell (l) considered
    int *l;
    int *m;
    // Orbitals
    // For each Z, there is one "slice":
    //        | l=0 |  1  |  2  |
    //   n=0  | --- | --- | --- |
    //   1    |     | --- | --- |
    //   2    |     |     | --- |
    //   3    |     |     |     |
    //   4    |     |     |     |
    //   5    |     |     |     |
    //          (s)   (p)   (d)
    // "---" are un-allowed states
    // array[l][n][Z] = array[l + n*max_l + Z*max_l*max_n]
    // n between 1 and max_n       (total: max_n)
    // l between 0 and max_n-1     (total: max_l)
    // Because C/C++ arrays start at 0, a column for n=0 will be
    // allocated too, but not used. This will keep code consistent
    int *occ;               // Occupied orbitals
    fdouble *Ips;            // Ionization potentials (eV)
    fdouble *IpsCummul;      // Cumulative ionization potentials (eV)
    fdouble *IpsLowest;      // Lowest ionization potentials for each charge states (eV)
    // aIp removed (use IpsCummul instead)
    fdouble *Lotz_a;        // cm^2 . (eV)^2
    fdouble *Lotz_aq;       // cm^2 . (eV)^2
    fdouble *Lotz_b;        // none
    fdouble *Lotz_c;        // none

    void set_pointers_to_null()
    {
        l         = NULL;
        m         = NULL;
        occ       = NULL;
        Ips       = NULL;
        IpsCummul = NULL;
        IpsLowest = NULL;
        Lotz_a    = NULL;
        Lotz_aq   = NULL;
        Lotz_b    = NULL;
        Lotz_c    = NULL;
    }

    void free_arrays()
    {
        free(l);
        free(m);
        free(occ);
        free(Ips);
        free(IpsCummul);
        free(IpsLowest);
        free(Lotz_a);
        free(Lotz_aq);
        free(Lotz_b);
        free(Lotz_c);

        set_pointers_to_null();
    }
};

#endif // #ifndef _ELEMENTSTRUCTURE_HPP

// ********** End of file ***************************************
