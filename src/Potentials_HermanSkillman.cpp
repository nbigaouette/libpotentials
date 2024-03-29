/****************************************************************

****************************************************************/

#include <iostream> // cout
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <iostream> // cout
#include <cstring>  // memset
#include <vector>

#include <StdCout.hpp>
#include <Assert.hpp>

#include "LibPotentials.hpp"
#include "Constants.hpp"
#include "Code_Functions_Declarations.hpp"
#include "LookUpTable.hpp"

#include "Potentials_Symmetric.hpp"
#include "Potentials_HermanSkillman.hpp"
#include "Potentials_Coulomb.hpp"



using namespace libpotentials;

bool USING_HS = false;


std::vector<fdouble> hs_min_rad;
std::vector<LookUpTable<fdouble> > hs_lut_potential;
std::vector<LookUpTable<fdouble> > hs_lut_field;

// Distance (in Bohr) where HS potential reaches the Coulombic values
const double HS_Xe_rmax[HS_Xe_MaxNbCS] = {
//  Neutral   1+           2+                   3+                  4+                   5+                6+
    12.0, 3.20671, 3.4419260284993465, 2.5786383389900931, 2.3069998849855851, 2.1424931697942897, 1.9288708307182563
};

// Function parameters: -a*exp(-r*b+c) - d*exp(-r*e+f) - g*exp(-r*h+i)
const double HS_Xe_parameters[HS_Xe_MaxNbCS][10] = {
//                   a                  b               c               d               e               f               g               h               i               j
/* Neutral  */  {3.14083382,        2.23690529,     2.45999159,     1.04922253,     0.758964055,    0.916259659,    6.43949225,     5.57500308,     3.46700894,     -2.22471387e-03},
/* 1+ */        {6.68943368,        5.46613511,     3.50223818,     9.84854609,     16.33928308,    4.55862051,     2.05998631,     1.79149357,     2.67105113,     -0.21651473},
/* 2+ */        {1.08862692,        1.12845509,     2.40634711,     7.5231977,      11.80857511,    4.40029841,     3.5171341,      4.02105327,     3.70863489,     -0.33244088},
/* 3+ */        {1.26082951,        1.24346292,     2.48202614,     7.60391482,     14.20436211,    4.63976132,     5.75320941,     4.57482337,     3.45935956,     -0.55091234},
/* 4+ */        {8.33659368,        15.53383795,    4.69224278,     5.66740119,     4.93161199,     3.58851214,     1.33122023,     1.36069086,     2.65699251,     -0.90941801},
/* 5+ */        {8.13621709,        15.39455048,    4.6973397,      1.33881001,     1.40783802,     2.72036815,     5.60695758,     4.96351559,     3.59035494,     -1.33283627},
/* 6+ */        {7.52331956,        15.56584267,    4.77821787,     2.17218048,     1.51817071,     2.38100923,     5.09462365,     5.11830058,     3.70739486,     -1.84326541},
};

double Adapt_HS_Cutoff_Radius_to_Match_Symmetric(const fdouble &cutoff_radius, const int cs)
/**
 * The cutoff radius from input file will be used to find where a symmetric field
 * would be maximum (r = Symmetric_r1). Then, a bisection is used to find the HS
 * cutoff radius that will place the maximum of the HS field at the same radius r1.
 */
{
    // Find the radius at which the derivative of the symmetric field is 0.
    // Setting y = r/(sqrt(2)*sigma)
    // http://www.wolframalpha.com/input/?i=diff%28erf%28y%29%2F%282*y^2%29+-+exp%28-y^2%29%2F%28sqrt%28pi%29+*+y%29%2Cy%29
    // But note that the symmetric potential was modified to have sigma *= sqrt_2/2 so
    // it would be the same as the Gaussian Distribution. Thus divided r1 here by sqrt2/2.
    const double Symmetric_y1 = 0.9678571637866073;
    const double Symmetric_r1 = (libpotentials_double::m_to_bohr * Potentials_Symmetric::sigma)
                                * (libpotentials_double::sqrt_2 * Symmetric_y1)
                                / (sqrt_2/two);

    // Do a bisection to find the right cutoff radius that will match "r1", the distance
    // where the HS field is maximum, to "Symmetric_r1", the distance where the
    // symmetric field is maximum.
    double cutoff_r_left, cutoff_r_right, cutoff_found_r; // [Bohr]
    // Initial conditions
    cutoff_r_left  = 0.0; // [Bohr]
    cutoff_r_right = 10.0;// [Bohr]

    // Start bisection!
    // See http://en.wikipedia.org/wiki/Bisection_method#Practical_considerations
    cutoff_found_r = cutoff_r_right + (cutoff_r_left - cutoff_r_right) / 2.0;

    while (std::abs(cutoff_found_r - cutoff_r_left) > 1.0e-100 && std::abs(cutoff_found_r - cutoff_r_right) > 1.0e-100)
    {
        // Find index of cutoff radius
        const int index = hs_lut_potential[cs].Get_i_from_x(fdouble(cutoff_found_r));
        assert(index > 0);

        const double dr = hs_lut_potential[cs].Get_dx();
        const double r0 = cutoff_found_r;

        // Get the last two points of the potential curve before cutoff
        //const double H0 = hs_lut_potential[cs].Table(index);
        //const double H1 = hs_lut_potential[cs].Table(index+1);
        const double F0 = hs_lut_field[cs].Table(index)   *  r0;
        const double F1 = hs_lut_field[cs].Table(index+1) * (r0 + dr);

        const double diffF0 = (F1 - F0) / dr;

        // r1 is the distance where the field is maximum
        const double r1 = (
                2.0*F0*r0 - r0*r0*diffF0
            ) / (
                2.0*F0 - 2.0*r0*diffF0
            );
        assert(r1 > 0.0);
        assert(r1 < r0);

        //log("cs=%2d   [r_left,r,r_right] = [%10.5g,%10.5g,%10.5g]     (r1,r0)=(%10.5g,%10.5g)", cs, cutoff_r_left, cutoff_found_r, cutoff_r_right, r1, r0);

        if (r1 >= Symmetric_r1)
        {
            cutoff_r_right = cutoff_found_r;
        }
        else
        {
            cutoff_r_left = cutoff_found_r;
        }
        cutoff_found_r = cutoff_r_right + (cutoff_r_left - cutoff_r_right) / 2.0;
    }
    return cutoff_found_r;
    //std_cout << "Bisection end: cs = " << cs << "  HS(r="<<cutoff_found_r<<")\n";
}


// **************************************************************
bool Is_HS_used()
{
    return USING_HS;
}

// **************************************************************
double HS_Fitting_Function_Xe_Potential(const double r, const int cs)
{
    const double a = HS_Xe_parameters[cs][0];
    const double b = HS_Xe_parameters[cs][1];
    const double c = HS_Xe_parameters[cs][2];
    const double d = HS_Xe_parameters[cs][3];
    const double e = HS_Xe_parameters[cs][4];
    const double f = HS_Xe_parameters[cs][5];
    const double g = HS_Xe_parameters[cs][6];
    const double h = HS_Xe_parameters[cs][7];
    const double i = HS_Xe_parameters[cs][8];
    const double j = HS_Xe_parameters[cs][9];
    const double potential =
        - a * std::exp(-r*b + c)
        - d * std::exp(-r*e + f)
        - g * std::exp(-r*h + i)
        + j;

    return potential;
}

// **************************************************************
double HS_Fitting_Function_Xe_Field(const double r, const int cs)
{
    const double a = HS_Xe_parameters[cs][0];
    const double b = HS_Xe_parameters[cs][1];
    const double c = HS_Xe_parameters[cs][2];
    const double d = HS_Xe_parameters[cs][3];
    const double e = HS_Xe_parameters[cs][4];
    const double f = HS_Xe_parameters[cs][5];
    const double g = HS_Xe_parameters[cs][6];
    const double h = HS_Xe_parameters[cs][7];
    const double i = HS_Xe_parameters[cs][8];
    //const double j = HS_Xe_parameters[cs][9];
    const double field = -(
          a * b * std::exp(-r*b + c)
        + d * e * std::exp(-r*e + f)
        + g * h * std::exp(-r*h + i));

    return field;
}

// **************************************************************
template <class Integer>
inline std::string IntToStr(const Integer integer, const int width = 0, const char fill = ' ')
{
    std::ostringstream MyStream;
    if (width != 0)
    {
        MyStream << std::setw(width);
        MyStream << std::setfill(fill);
    }
    MyStream << integer << std::flush;
    return (MyStream.str());
}

// **************************************************************
void Initialize_HermanSkillman(const fdouble cutoff_base_potential, const fdouble cutoff_radius)
{
    USING_HS = true;

    // We'll need one lookup table per charge state
    std_cout << "FIXME: Dynamically choose between atom types for HS (" << __FILE__ << ", line " << __LINE__ << ")\n";
    // LUTs stored in atomic units
    Set_HermanSkillman_Lookup_Tables_Xe(hs_lut_potential, hs_lut_field);

    hs_min_rad.resize(hs_lut_potential.size());

    // HS will fallback to symmetric if charge state is not HS fitted or if distance is more then
    // the HS cutoff.
    Initialize_Symmetric(cutoff_base_potential, cutoff_radius);

    // If both radius and base potential are positive, we reloaded from a snapshot. Call the
    // radius cutoff function with the scaling disabled.
    if      (cutoff_base_potential > 0.0 and cutoff_radius > 0.0)
        Initialize_HS_Cutoff_Radius(cutoff_radius, false);
    else if (cutoff_base_potential > 0.0)
        Initialize_HS_Base_Potential(cutoff_base_potential);
    else if (cutoff_radius > 0.0)
        Initialize_HS_Cutoff_Radius(cutoff_radius);

}

// **************************************************************
void Set_HermanSkillman_Lookup_Tables_Xe(std::vector<LookUpTable<fdouble> > &lut_pot,
                                         std::vector<LookUpTable<fdouble> > &lut_field)
{
    const int lut_n = 10000;

    lut_pot.resize(HS_Xe_MaxNbCS);
    lut_field.resize(HS_Xe_MaxNbCS);

    // Set potential lookup table
    // cs: 0 == neutral, 1 == 1+, etc.
    for (int cs = 0 ; cs < HS_Xe_MaxNbCS ; cs++)
    {
        // Distance range for the HS.
        const double xmin = 0.0;    // [Bohr]
        const double xmax = 15.0;   // [Bohr]
        lut_pot[cs].Initialize(  NULL, fdouble(xmin), fdouble(xmax), lut_n, "Initialize_HS() LookUpTable (lut_pot, cs=" + IntToStr(cs) + ")");
        lut_field[cs].Initialize(NULL, fdouble(xmin), fdouble(xmax), lut_n, "Initialize_HS() LookUpTable (lut_field, cs=" + IntToStr(cs) + ")");

        // FIXME: Dynamically choose between atom types for HS
        // Start: use function pointers
        double (*HS_Fitting_Function_Potential)(const double r, const int cs);
        double (*HS_Fitting_Function_Field)(const double r, const int cs);

        HS_Fitting_Function_Potential = &HS_Fitting_Function_Xe_Potential;
        HS_Fitting_Function_Field     = &HS_Fitting_Function_Xe_Field;

        // Populate the lookup tables.
        // The HS potential is fitted. But the resulting field has a discontinuity with Coulomb.
        // To fix, scale the HS field at rmax so it fits the Coulomb field. This add
        // a multiplicative constant to the field. But it needs to be added to the potential too!
        // Then, to make sure the potential fits Coulomb at rmax, add an additive constant, which
        // is not included in the field (due to E=-grad(V)).
        double r;
        const double HS_E_xmax      = -HS_Fitting_Function_Field(HS_Xe_rmax[cs], cs);
        const double Coulomb_E_xmax = double(cs) / (HS_Xe_rmax[cs]*HS_Xe_rmax[cs]);
        const double HS_U_xmax      = HS_Fitting_Function_Potential(HS_Xe_rmax[cs], cs);
        const double Coulomb_U_xmax = -double(cs) / (HS_Xe_rmax[cs]);
              double HS_E_scaling_factor = Coulomb_E_xmax / HS_E_xmax;
              double HS_U_add_factor     = Coulomb_U_xmax - HS_E_scaling_factor*HS_U_xmax;
        // Prevent neutral being scaled by 0
        if (cs == 0)
        {
            HS_U_add_factor = -HS_U_xmax;
            HS_E_scaling_factor = 1.0;
        }

        //std_cout << "cs="<<cs<<"  HS_E_scaling_factor = " <<HS_E_scaling_factor << "\n";
        for (int i = 0 ; i < lut_n ; i++)
        {
            r = lut_pot[cs].Get_x_from_i(i);  // [Bohr]

            double HS_U         = 0.0;  // [Hartree]
            double HS_E_over_r  = 0.0;  // [au/Bohr]

            // Electrostatic potential
            if (r >= HS_Xe_rmax[cs])
            {
                // Use Coulomb
                HS_U = -double(cs) / r;
            }
            else
            {
                // Use HS
                HS_U = HS_Fitting_Function_Potential(r, cs);
                // Scale to prevent discontinuities
                HS_U *= HS_E_scaling_factor;
                HS_U += HS_U_add_factor;
            }

            // Electrostatic field
            if (r < HS_Xe_rmax[cs])
            {
                HS_E_over_r = -HS_Fitting_Function_Field(r, cs);
                // Scale to prevent discontinuities
                HS_E_over_r *= HS_E_scaling_factor;
            }
            else
                HS_E_over_r = double(cs) / (r*r);


            // We store E/r, not E
            if (r > 1.0e-10) // Bohr
                HS_E_over_r /= r;

            Assert_isinf_isnan(HS_U);

            // Save it (in atomic units)
            lut_pot[cs].Set(  i, fdouble(HS_U));
            lut_field[cs].Set(i, fdouble(HS_E_over_r));
        }
    }
}

// **************************************************************
void Initialize_HS_Cutoff_Radius(const fdouble &cutoff_radius_m, const bool scale_radius_to_symmetric)
/**
 * Initialize Herman-Skillman potential with a given cutoff radius
 * @param   cutoff_radius   Cutoff radius under which a smoothing is wanted [m]
 */
{
    const fdouble cutoff_radius = cutoff_radius_m*libpotentials::m_to_bohr;

    std_cout
        << "Calling Initialize_HS_Cutoff_Radius():\n"
        << "    The cutoff radius from input file (" << cutoff_radius << " bohr) will be used\n"
        << "    to find where a symmetric field would be maximum (r = Symmetric_r1). Then, a\n"
        << "    bisection is performed to find the HS cutoff radius that will place the maximum\n"
        << "    of the HS field at the same radius r1." << std::endl;

    // Change lookup tables values using a smoothed potential/field inside cutoff radius
    for (int cs = 0 ; cs < int(hs_lut_potential.size()) ; cs++)
    {
        // Set every charge states cutoff radius to the wanted value
        hs_min_rad[cs] = cutoff_radius;

        // Flag to enable/disable the cutoff setting.
        // If false:
        //      The cutoff radius from input file (called "r0") will separate the regions where the
        //      potential/field are either HS (r>r0) or a smoothed harmonic value (r<r0).
        // If true:
        //      The cutoff radius from input file will be used to find where a symmetric field
        //      would be maximum (r = Symmetric_r1). Then, a bisection is used to find the HS
        //      cutoff radius that will place the maximum of the HS field at the same radius r1.
        //const bool hs_do_bisection_r0_matching_symmetric = false;
        const bool hs_do_bisection_r0_matching_symmetric = true;

        // Set the HS's cutoff radius to a value which will put the maximum of the electric field
        // at the same radius as the Symmetric case.
        if (hs_do_bisection_r0_matching_symmetric and scale_radius_to_symmetric)
        {
            hs_min_rad[cs] = fdouble( Adapt_HS_Cutoff_Radius_to_Match_Symmetric(cutoff_radius, cs) );
        }

        // Find index of cutoff radius
        const int index = hs_lut_potential[cs].Get_i_from_x(hs_min_rad[cs]);
        assert(index > 0);

        const double dr = hs_lut_potential[cs].Get_dx();
        const double r0 = hs_min_rad[cs];

        // Get the last two points of the potential curve before cutoff
        const double H0 = hs_lut_potential[cs].Table(index);
        //const double H1 = hs_lut_potential[cs].Table(index+1);
        const double F0 = hs_lut_field[cs].Table(index)   *  r0;
        const double F1 = hs_lut_field[cs].Table(index+1) * (r0 + dr);

        const double diffF0 = (F1 - F0) / dr;

        // r1 is the distance where the field is maximum
        const double r1 = (
                2.0*F0*r0 - r0*r0*diffF0
            ) / (
                2.0*F0 - 2.0*r0*diffF0
            );
        assert(r1 > 0.0);
        assert(r1 < r0);
        const double B = F0 - (r0-r1)*diffF0/2.0;
        const double A = B / (r1*r1);
        const double C = H0 + A*r0*r0*r0/3.0 - A*r0*r0*r1;

        /*
        if (cs == 1)
        {
            std_cout << "index = " << index << "\n";
            std_cout << "HS_Fitting_Function_Xe_Potential(1.0 bohr, 1+) = " << HS_Fitting_Function_Xe_Potential(1.0, 1) << "\n";
            std_cout << "hs_lut_potential[cs].Table(index) = " << hs_lut_potential[cs].Table(index) << "\n";
            std_cout << "hs_lut_potential[cs].Table(index+1) = " << hs_lut_potential[cs].Table(index+1) << "\n";
            std_cout << "r1 = " << r1 << "\n";

            std_cout << "A = " << A << "\n";
            std_cout << "B = " << B << "\n";
            std_cout << "C = " << C << "\n";

            std_cout << "Hs(0.5 bohr) = " << hs_lut_potential[cs].read(0.5) << "\n";
            std_cout << "Hs(1.0 bohr) = " << hs_lut_potential[cs].read(1.0) << "\n";
            std_cout << "HsF(0.5 bohr) = " << hs_lut_field[cs].read(0.5) << "\n";
            std_cout << "HsF(1.0 bohr) = " << hs_lut_field[cs].read(1.0) << "\n";
            std_cout << "H0 = " << H0 << "\n";
            std_cout << "H1 = " << H1 << "\n";
            std_cout << "F0 = " << F0 << "\n";
            std_cout << "F1 = " << F1 << "\n";
            assert(A > 0.0);
            assert(B > 0.0);
            assert(C < 0.0);
            assert(hs_lut_potential[cs].read(r0) < 0.0);
            assert(hs_lut_field[cs].read(r0)     > 0.0);
        }
        */

        for (int i = 0 ; i < index ; i++)
        {
            const double r = hs_lut_potential[cs].Get_x_from_i(i);
            //hs_lut_potential[cs].Set(i, fdouble(-A*std::pow(r,3)/3.0 - A*std::pow(r,2)*r1 + C));
            hs_lut_potential[cs].Set(i, fdouble(-A*std::pow(r,3)/3.0 - A*std::pow(r1,2)*r + A*r*r*r1 + B*r + C));
            if (i == 0)
                hs_lut_field[cs].Set(i, 0.0 ); // We store E/r but diverge at r == 0
            else
                hs_lut_field[cs].Set(i, fdouble((-A*std::pow(r-r1,2) + B) / r) ); // We store E/r
        }
    }

    // Store the potential of a 1+ at r == 0 as being the cutoff_base_potential
    libpotentials_private::cutoff_base_potential = Eh_to_eV * std::abs(hs_lut_potential[1].Table(0));
    libpotentials_private::cutoff_radius = hs_min_rad[1] * bohr_to_m;

    log("Initialize_HS_Cutoff_Radius()  libpotentials_private::cutoff_base_potential = %15.13g eV\n", libpotentials_private::cutoff_base_potential);
    log("Initialize_HS_Cutoff_Radius()  libpotentials_private::cutoff_radius = %15.13g bohr\n", libpotentials_private::cutoff_radius * libpotentials::m_to_bohr);

}

// **************************************************************
void Initialize_HS_Base_Potential(const fdouble &base_potential_eV)
/**
 * Initialize Herman-Skillman potential with a given potential depth
 * @param   base_potential_eV  Potential depth wanted [eV]
 */
{
    const fdouble base_potential = -std::abs(base_potential_eV)*libpotentials::eV_to_Eh;

    // Find the radius where the HS potential is equal to "base_potential"
    // by doing a bisection, for all supported charge states.
    fdouble r_left, r_right, found_r, pot; // [Bohr]
    for (int cs = 0 ; cs < int(hs_lut_potential.size()) ; cs++)
    {
        // Initial conditions
        hs_min_rad[cs] = hs_lut_potential[cs].Get_x_from_i(0);
        r_left  = hs_min_rad[cs];     // Minimum radius of fit
        r_right = hs_lut_potential[cs].Get_XMax();

        // Start bisection!
        // See http://en.wikipedia.org/wiki/Bisection_method#Practical_considerations
        found_r = r_right + (r_left - r_right) / fdouble(2.0);

        while (std::abs(found_r - r_left) > 1.0e-100 && std::abs(found_r - r_right) > 1.0e-100)
        {
            pot = hs_lut_potential[cs].read(found_r); // Potential energy [au]
            assert(pot < 0.0);

            // We want a 2+'s base potential to be twice as deep as a 1+. We will
            // thus compare the base potential with the potential divided by the
            // charge state, but only for 1+ and up.

            //printf("cs=%2d  base_potential = %10.5g Hartree   r_left = %10.5g   r = %10.5g   r_right = %10.5g   HS(r) = %10.5g Hartree\n", cs, base_potential, r_left, found_r, r_right, pot);
            Assert_isinf_isnan(pot);
            if (pot >= base_potential*fdouble(std::max(1,cs)))
            {
                r_right = found_r;
            }
            else
            {
                r_left = found_r;
            }
            found_r = r_right + (r_left - r_right) / fdouble(2.0);
        }
        std_cout << "Bisection end: cs = " << cs << "  HS(r="<<found_r<<") = " << double(std::max(1,cs))*pot << "\n";

        Assert_isinf_isnan(found_r);
        Assert_isinf_isnan(pot);
        assert(found_r > 0.0);
        //assert(std::abs(pot - std::max(1,cs)*base_potential) < 1.0e-3);

        if (found_r <= 0.99999*hs_min_rad[cs])
        {
            std_cout << "##############################################\n";
            DEBUGP("Initialize_HS() called with a potential depth too deep.\n");
            std_cout << "The value found " << found_r << " Bohr\n";
            std_cout << "for charge state " << cs << " should not be lower than " << hs_min_rad[cs] << " Bohr\n";
            std_cout << "Potential depth wanted: " << base_potential << " Hartree (" << base_potential_eV << " eV)\n";
            std_cout << "Exiting\n";
            abort();
        }

        hs_min_rad[cs] = found_r;
    }

    // Now that the cutting radius is found for each charge states, change lookup tables values
    // This is a "hard" cutoff: the field inside hs_min_rad will be 0, and the
    // potential will be the base_potential.
    for (int cs = 0 ; cs < int(hs_lut_potential.size()) ; cs++)
    {
        // Set neutral's charge state to 1, so it does not clear the lookup tables.
        const fdouble cs_factor = fdouble(std::max(1, cs));
        const int lut_n = hs_lut_potential[cs].Get_n();
        for (int i = 0 ; i < lut_n ; i++)
        {
            assert(base_potential < 0.0);
            if (hs_lut_potential[cs].Table(i) < base_potential*cs_factor)
            {
                hs_lut_potential[cs].Set(i, base_potential*cs_factor);
                hs_lut_field[cs].Set(i, 0.0);
            }
        }
    }

#ifdef HSSPLINE
    // Do a cubic spline interpolation on the field to prevent the drop from the maximum
    // of field to 0 at cutoff radius. This hard cutoff introduce a lot of numerical heating.
    // The spline should make it smooth and thus no more heating.
    for (int cs = 0 ; cs < int(hs_lut_potential.size()) ; cs++)
    {
        const int lut_n = hs_lut_potential[cs].Get_n();
        // Find index of maximum field
        int iEmax = -1;
        fdouble max_field = 0.0;
        for (int i = 0 ; i < lut_n ; i++)
        {
            // The field's lut stores E/r, not E.
            const fdouble r_i = hs_lut_field[cs].Get_x_from_i(i);
            if (r_i*hs_lut_field[cs].Table(i) > max_field)
            {
                max_field = r_i*hs_lut_field[cs].Table(i);
                iEmax = i;
            }
        }
        assert(iEmax != -1);

        // Do a bisection to find the proper factor to use for a given potential depth and charge state
        fdouble f_up, f_down, found_f;
        // Initial conditions
        f_down = 2.0f;
        f_up   = 10.0f;

        // Start bisection!
        // See http://en.wikipedia.org/wiki/Bisection_method#Practical_considerations
        found_f = f_down + (f_up - f_down) / fdouble(2.0);
        while (std::abs(found_f - f_down) > 1.0e-100 && std::abs(found_f - f_up) > 1.0e-100)
        {
            // Find the index where the electric field is found_f times less then the max
            int iE = -1;
            for (int i = iEmax ; i < lut_n ; i++)
            {
                // The field's lut stores E/r, not E.
                const fdouble r_i   = hs_lut_field[cs].Get_x_from_i(i);
                const fdouble r_ip1 = hs_lut_field[cs].Get_x_from_i(i+1);
                const fdouble E_i   = r_i  *hs_lut_field[cs].Table(i);
                const fdouble E_ip1 = r_ip1*hs_lut_field[cs].Table(i+1);
                if ( (E_i >= max_field/found_f) and (max_field/found_f > E_ip1) )
                {
                    iE = i;
                    break;
                }
            }
            assert(iE != -1);

            // Store the 3 points used for the cubic spline
            const int n = 2;
            const int intervals[n+1] = {0, iE, iE+10};
            //std_cout
            //    << "intervals = ("<<intervals[0] << "," << intervals[1] << "," << intervals[2] << ")    "
            //    << "x = (" << hs_lut_field[cs].Get_x_from_i(intervals[0]) << "," << hs_lut_field[cs].Get_x_from_i(intervals[1]) << "," << hs_lut_field[cs].Get_x_from_i(intervals[1]) << ")\n";
            const fdouble p0[2] = {hs_lut_field[cs].Get_x_from_i(intervals[0]), hs_lut_field[cs].Get_x_from_i(intervals[0])*hs_lut_field[cs].Table(intervals[0])};
            const fdouble p1[2] = {hs_lut_field[cs].Get_x_from_i(intervals[1]), hs_lut_field[cs].Get_x_from_i(intervals[1])*hs_lut_field[cs].Table(intervals[1])};
            const fdouble p2[2] = {hs_lut_field[cs].Get_x_from_i(intervals[2]), hs_lut_field[cs].Get_x_from_i(intervals[2])*hs_lut_field[cs].Table(intervals[2])};

            const fdouble r[n+1]= {p0[0], p1[0], p2[0]};
            const fdouble y[n+1]= {p0[1], p1[1], p2[1]};

            const fdouble h[n]  = {r[1] - r[0], r[2] - r[1]};

            // Construct cubic spline coefficients. See doc/cubic_splines/cubic_spline.pdf
            const fdouble two   = libpotentials::two;
            const fdouble three = libpotentials::three;
            const fdouble four  = libpotentials::four;
            const fdouble six   = libpotentials::six;
            // For a 3 points "natural" spline, the solution "m" to "A.m = b" is easy:
            const fdouble m[n+1] = {0.0, (three / (h[0]-h[1])) * ( (y[2]-y[1])/h[1] - (y[1]-y[0])/h[0] ), 0.0};
            // Get the coefficients from "m"
            fdouble a[n];
            fdouble b[n];
            fdouble c[n];
            fdouble d[n];
            for (int i = 0 ; i < n ; i++)
            {
                a[i] = y[i];
                b[i] = (y[i+1] - y[i])/h[i] - h[i]*m[i]/two - h[i]/six*(m[i+1] - m[i]);
                c[i] = m[i] / two;
                d[i] = (m[i+1] - m[i]) / (six * h[i]);
            }
            //for (int i = 0 ; i < n ; i++)
            //{
            //    std_cout << "    i="<<i<<" (a,b,c,d,m) = ("<<a[i]<<","<<b[i]<<","<<c[i]<<","<<d[i]<<","<<m[i]<<")\n";
            //}

            // For the two regions defined by the three points, calculate the cubic spline
            fdouble pot_at_end_section = 0.0;
            fdouble x = 0.0;
            // Set neutral's charge state to 1, so it does not clear the lookup tables.
            const fdouble cs_factor = fdouble(std::max(1, cs));
            fdouble last_IntV = base_potential*cs_factor;
            for (int interval = 0 ; interval < n ; interval++)
            {
                x = hs_lut_field[cs].Get_x_from_i(intervals[interval+1]);
                pot_at_end_section = a[interval]*(x - r[interval]) + b[interval]*std::pow(x - r[interval],2)/two   + c[interval]*std::pow(x - r[interval], 3)/three + d[interval]*std::pow(x - r[interval], 4)/four + last_IntV;
                last_IntV = pot_at_end_section;
                //std_cout << "intervals[interval="<<interval<<"] = " << intervals[interval] << "  ->  (x,r) = (" << x << "," << r[interval] << ")   pot_at_end_section = " << pot_at_end_section
                //    << "    (a,b,c,d,m) = ("<<a[interval]<<","<<b[interval]<<","<<c[interval]<<","<<d[interval]<<","<<m[interval]<<")"
                //    << "\n";
            }

            //printf("cs=%2d  f_up = %5.3g   f = %5.3g   f_down = %5.3g int[n]=%d  Table(int[n])=%6.3g  r(int[n])=%5.3g  last_IntV=%5.3g\n",
            //       cs, f_up, found_f, f_down, intervals[n], hs_lut_potential[cs].Table(intervals[n]), hs_lut_potential[cs].Get_x_from_i(intervals[n]), last_IntV);
            if (last_IntV >= hs_lut_potential[cs].Table(intervals[n]))
            {
                f_down = found_f;
            }
            else
            {
                f_up = found_f;
            }
            found_f = f_down + (f_up - f_down) / fdouble(2.0);
        }

        std_cout << "cs=" << cs << "  found_f = " << found_f << "\n";

        // Find the index where the electric field is 6.67 times less then the max
        const fdouble factor = found_f;
        int iE = -1;
        for (int i = iEmax ; i < lut_n ; i++)
        {
            // The field's lut stores E/r, not E.
            const fdouble r_i   = hs_lut_field[cs].Get_x_from_i(i);
            const fdouble r_ip1 = hs_lut_field[cs].Get_x_from_i(i+1);
            const fdouble E_i   = r_i  *hs_lut_field[cs].Table(i);
            const fdouble E_ip1 = r_ip1*hs_lut_field[cs].Table(i+1);
            if ( (E_i >= max_field/factor) and (max_field/factor > E_ip1) )
            {
                iE = i;
                break;
            }
        }
        assert(iE != -1);

        // Store the 3 points used for the cubic spline
        const int n = 2;
        const int intervals[2*n] = {0, iE, iE+10};
        const fdouble p0[2] = {hs_lut_field[cs].Get_x_from_i(intervals[0]), hs_lut_field[cs].Get_x_from_i(intervals[0])*hs_lut_field[cs].Table(intervals[0])};
        const fdouble p1[2] = {hs_lut_field[cs].Get_x_from_i(intervals[1]), hs_lut_field[cs].Get_x_from_i(intervals[1])*hs_lut_field[cs].Table(intervals[1])};
        const fdouble p2[2] = {hs_lut_field[cs].Get_x_from_i(intervals[2]), hs_lut_field[cs].Get_x_from_i(intervals[2])*hs_lut_field[cs].Table(intervals[2])};

        const fdouble r[n+1]= {p0[0], p1[0], p2[0]};
        const fdouble y[n+1]= {p0[1], p1[1], p2[1]};

        const fdouble h[n]  = {r[1] - r[0], r[2] - r[1]};

        // Construct cubic spline coefficients. See doc/cubic_splines/cubic_spline.pdf
        const fdouble two   = libpotentials::two;
        const fdouble three = libpotentials::three;
        const fdouble four  = libpotentials::four;
        const fdouble six   = libpotentials::six;
        // For a 3 points "natural" spline, the solution "m" to "A.m = b" is easy:
        const fdouble m[n+1] = {0.0, (three / (h[0]-h[1])) * ( (y[2]-y[1])/h[1] - (y[1]-y[0])/h[0] ), 0.0};
        // Get the coefficients from "m"
        fdouble a[n];
        fdouble b[n];
        fdouble c[n];
        fdouble d[n];
        for (int i = 0 ; i < n ; i++)
        {
            a[i] = y[i];
            b[i] = (y[i+1] - y[i])/h[i] - h[i]*m[i]/two - h[i]/six*(m[i+1] - m[i]);
            c[i] = m[i] / two;
            d[i] = (m[i+1] - m[i]) / (six * h[i]);
        }

        // For the two regions defined by the three points, calculate the cubic spline
        fdouble new_field = 0.0;
        fdouble new_pot   = 0.0;
        fdouble old_pot   = 0.0;
        fdouble x = 0.0;
        // Set neutral's charge state to 1, so it does not clear the lookup tables.
        const fdouble cs_factor = fdouble(std::max(1, cs));
        fdouble last_IntV = base_potential*cs_factor;
        for (int interval = 0 ; interval < n ; interval++)
        {
            for (int i = intervals[interval] ; i <= intervals[interval+1] ; i++)
            {
                x = hs_lut_field[cs].Get_x_from_i(i);
                new_field = a[interval]                   + b[interval]*        (x - r[interval])         + c[interval]*std::pow(x - r[interval], 2)       + d[interval]*std::pow(x - r[interval], 3);
                new_pot   = a[interval]*(x - r[interval]) + b[interval]*std::pow(x - r[interval],2)/two   + c[interval]*std::pow(x - r[interval], 3)/three + d[interval]*std::pow(x - r[interval], 4)/four + last_IntV;
                old_pot   = hs_lut_potential[cs].Table(i);

                // As soon as the spline crosses the old potential, stop using it.
                new_pot = std::max(new_pot, old_pot);
                // We store E/r, not E
                if (x > 0.0)
                    new_field /= x;

                hs_lut_field[cs].Set(i, new_field);
                hs_lut_potential[cs].Set(i, new_pot);
            }
            last_IntV = new_pot;
        }
    }
#endif // #ifdef HSSPLINE

    //hs_lut_potential[1].Print_Table();
    //hs_lut_field[1].Print_Table();
    //exit(0);

    /*
    // Print lookup table for verification
    // To run: ./mdgit 2> lut_hs.dat
    const int max_lut = HS_Xe_MaxNbCS;
    std::string filename("lut_hs.dat");
    std::string gnuplot_command("");
    gnuplot_command += "#set term wxt 3; plot ";
    int row = 1;
    for (int l = 0 ; l < max_lut ; l++)
    {
        gnuplot_command += "\"" + filename + "\" using " + IntToStr(row++) + ":";
        gnuplot_command += IntToStr(row++) + "  title \"LUT(V(" + IntToStr(l) + ")) (HS)\" with lines lw 3";
        if (l == max_lut-1)
            gnuplot_command += "\n";
        else
            gnuplot_command += ", ";
    }
    gnuplot_command += "#set term wxt 4; plot ";
    for (int l = 0 ; l < max_lut ; l++)
    {
        gnuplot_command += "\"" + filename + "\" using " + IntToStr(row++) + ":";
        gnuplot_command += IntToStr(row++) + "  title \"cs*LUT(V(" + IntToStr(l) + ")) (HS)\" with lines lw 3, ";
    }
    for (int l = 1 ; l <= max_lut-2 ; l++)
    {
        gnuplot_command += "-" + IntToStr(l) + "/x with points";
        if (l == max_lut-2)
            gnuplot_command += "\n";
        else
            gnuplot_command += ", ";
    }
    gnuplot_command += "#set term wxt 5; plot ";
    for (int l = 0 ; l < max_lut ; l++)
    {
        gnuplot_command += "\"" + filename + "\" using " + IntToStr(row++) + ":";
        gnuplot_command += IntToStr(row++) + "  title \"LUT(E/r(" + IntToStr(l) + ")) (HS)\" with lines lw 3";
        if (l == max_lut-1)
            gnuplot_command += "\n";
        else
            gnuplot_command += ", ";
    }
    gnuplot_command += "#set term wxt 6; plot ";
    for (int l = 0 ; l < max_lut ; l++)
    {
        gnuplot_command += "\"" + filename + "\" using " + IntToStr(row++) + ":";
        gnuplot_command += IntToStr(row++) + "  title \"r*cs*LUT(E/r(" + IntToStr(l) + ")) (HS)\" with lines lw 3, ";
    }
    for (int l = 1 ; l <= max_lut-2 ; l++)
    {
        gnuplot_command += "-" + IntToStr(l) + "/x**2 with points";
        if (l == max_lut-2)
            gnuplot_command += "\n";
        else
            gnuplot_command += ", ";
    }
    fprintf(stderr, "%s", gnuplot_command.c_str());
    std::string::size_type pos = 0;
    std::string searchString("#");
    while ((pos = gnuplot_command.find(searchString, pos)) != std::string::npos)
    {
        gnuplot_command.replace(pos, searchString.size(), "");
        pos++;
    }
    std_cout << "\n" << gnuplot_command << "\n";
    fprintf(stderr, "#");
    for (int l = 0 ; l < max_lut ; l++)
    {
        fprintf(stderr, "%20s ", "r (bohr)");
        fprintf(stderr, "%20s ", std::string("V(" + IntToStr(l) + ")").c_str());
    }
    for (int l = 0 ; l < max_lut ; l++)
    {
        fprintf(stderr, "%20s ", "r (bohr)");
        fprintf(stderr, "%20s ", std::string("cs*V(" + IntToStr(l) + ")").c_str());
    }
    for (int l = 0 ; l < max_lut ; l++)
    {
        fprintf(stderr, "%20s ", "r (bohr)");
        fprintf(stderr, "%20s ", std::string("E/r(" + IntToStr(l) + ")").c_str());
    }
    for (int l = 0 ; l < max_lut ; l++)
    {
        fprintf(stderr, "%20s ", "r (bohr)");
        fprintf(stderr, "%20s ", std::string("cs*E(" + IntToStr(l) + ")").c_str());
    }
    fprintf(stderr, "\n");

    const int lut_n = hs_lut_potential[0].Get_n();
    for (int i = 0 ; i < lut_n ; i++)
    {
        for (int lut_index = 0 ; lut_index < max_lut ; lut_index++)
        {
            const fdouble xmin   = hs_lut_potential[lut_index].Get_XMin();
            const float distance = float(i)/hs_lut_potential[lut_index].Get_inv_dx() + xmin;
            fprintf(stderr, "%20.15g ", distance);
            fprintf(stderr, "%20.15g ", hs_lut_potential[lut_index].read(distance));
        }
        for (int lut_index = 0 ; lut_index < max_lut ; lut_index++)
        {
            const fdouble xmin   = hs_lut_potential[lut_index].Get_XMin();
            const float distance = float(i)/hs_lut_potential[lut_index].Get_inv_dx() + xmin;
            fprintf(stderr, "%20.15g ", distance);
            int cs = lut_index-1;
            if (cs == 0)
                cs = 1;
            fprintf(stderr, "%20.15g ", fdouble(cs)*hs_lut_potential[lut_index].read(distance));
        }
        for (int lut_index = 0 ; lut_index < max_lut ; lut_index++)
        {
            const fdouble xmin   = hs_lut_potential[lut_index].Get_XMin();
            const float distance = float(i)/hs_lut_potential[lut_index].Get_inv_dx() + xmin;
            fprintf(stderr, "%20.15g ", distance);
            fprintf(stderr, "%20.15g ", hs_lut_field[lut_index].read(distance));
        }
        for (int lut_index = 0 ; lut_index < max_lut ; lut_index++)
        {
            const fdouble xmin   = hs_lut_potential[lut_index].Get_XMin();
            const float distance = float(i)/hs_lut_potential[lut_index].Get_inv_dx() + xmin;
            int cs = lut_index-1;
            if (cs == 0)
                cs = 1;
            fprintf(stderr, "%20.15g ", distance);
            fprintf(stderr, "%20.15g ", distance*fdouble(cs)*hs_lut_field[lut_index].read(distance));
        }
        fprintf(stderr, "\n");
    }
    std_cout << std::flush;
//     hs_lut_field[0].Print_Table();
//     hs_lut_potential[0].Print_Table();
    exit(0);
    */
}

// **************************************************************
void Potentials_Set_Parameters_HS(
    void *p1, void *p2,
    potential_paramaters &potparams)
/**
 * p1 is the particle of interest, p2 is "the other" particle
 * acting on p1.
 *
 */
{
    Check_if_LibPotentials_is_initialized();

    set_vector_between_particles(Get_Position(p1), Get_Position(p2),
                                  potparams.dr, potparams.r2,
                                  potparams.r, potparams.one_over_r);

    potparams.kQ2           = one_over_4Pieps0 * Get_Charge(p2);
    potparams.hs_cs1        = Get_Charge_State(p1);
    potparams.hs_cs2        = Get_Charge_State(p2);

    // Fits are in atomic units
    const fdouble distance_au = potparams.r * si_to_au_length;

    // Find the right lookup table to use.
    // LUT indices: 0 == neutral, 1 == 1+, etc.
    // The max() will pick up an ion if present. If two electrons, abs() will make it a 1+
    potparams.hs_lut_i = std::abs(std::max(potparams.hs_cs1, potparams.hs_cs2));

    if (
        (potparams.hs_cs1 >= 0 and potparams.hs_cs2 >= 0) or                // Two ions
        (potparams.hs_lut_i >= hs_lut_potential.size())   or                // Too high charge state
        (distance_au >= hs_lut_potential[potparams.hs_lut_i].Get_XMax())    // Outside lookup table range
    )
    {
        // Coulombic
        Potentials_Set_Parameters_ChargeDistribution_Symmetric(p1, p2, potparams);
    }
}

// **************************************************************
fdouble Calculate_Potential_Cutoff_HS(
    void *p1, void *p2,
    potential_paramaters &potparams)
/**
 */
{
    Check_if_LibPotentials_is_initialized();

    fdouble phi12 = 0.0;   // Electrostatic potential [Volt]

    // Fits are in atomic units
    const fdouble distance_au = potparams.r * si_to_au_length;

    if (
        (potparams.hs_cs1 >= 0 and potparams.hs_cs2 >= 0) or                // Two ions
        (potparams.hs_lut_i >= hs_lut_potential.size())   or                // Too high charge state
        (distance_au >= hs_lut_potential[potparams.hs_lut_i].Get_XMax())    // Outside lookup table range
    )
    {
        // Coulombic
        phi12 = Calculate_Potential_Cutoff_ChargeDistribution_Symmetric(p1, p2, potparams);
    }
    else
    {
        // Get the lookup table's value. NOTE: The potential ENERGY of an electron in an atomic potential is stored [Hartree].
        phi12 = std::abs(hs_lut_potential[potparams.hs_lut_i].read(distance_au)); // [Hartree]

        // We get |LUT|/|cs1*cs2|
        if (potparams.hs_cs1 > 0)
            phi12 /= fdouble(potparams.hs_cs1);
        if (potparams.hs_cs2 > 0)
            phi12 /= fdouble(potparams.hs_cs2);

        phi12 *= fdouble(potparams.hs_cs2);
        phi12 *= au_to_si_pot;
    }

    //std_cout << "(cs1,cs2)=(" << Get_Charge_State(p1) << ", " << Get_Charge_State(p2) << ")  r = " << distance_au << " Bohr   phi12 = " << phi12 << " Volt   (Coulomb: cs1*cs2/r = " << (fdouble(Get_Charge_State(p1)*Get_Charge_State(p2))/distance_au)*au_to_si_pot << ")     hs_lut_potential[potparams.hs_lut_i].Get_XMax() = " << hs_lut_potential[potparams.hs_lut_i].Get_XMax() << "\n";

    return phi12;
}

// **************************************************************
void Set_Field_Cutoff_HS(
    void *p1, void *p2,
    potential_paramaters &potparams,
    fdouble &phi, fdouble E[3])
{
    Check_if_LibPotentials_is_initialized();

    fdouble E_over_r = 0.0;   // Electrostatic field over distance r

    // Fits are in atomic units
    fdouble distance_au = potparams.r * si_to_au_length;

    if (
        (potparams.hs_cs1 >= 0 and potparams.hs_cs2 >= 0) or                // Two ions
        (potparams.hs_lut_i >= hs_lut_potential.size())   or                // Too high charge state
        (distance_au >= hs_lut_potential[potparams.hs_lut_i].Get_XMax())    // Outside lookup table range
    )
    {
        // Electron or high charge state ion
        assert(p1 != NULL);
        assert(p2 != NULL);
        Set_Field_Cutoff_ChargeDistribution_Symmetric(p1, p2, potparams, phi, E);
    }
    else
    {
        E_over_r = hs_lut_field[potparams.hs_lut_i].read(distance_au);

        // We get LUT/|cs1*cs2|
        if (potparams.hs_cs1 > 0)
            E_over_r /= fdouble(potparams.hs_cs1);
        if (potparams.hs_cs2 > 0)
            E_over_r /= fdouble(potparams.hs_cs2);

        E_over_r *= fdouble(potparams.hs_cs2);
        E_over_r *= au_to_si_field;
        E_over_r /= au_to_si_length;

        for (int d = 0 ; d < 3 ; d++)
        {
            E[d]  += potparams.dr[d] * E_over_r;
        }
    }
}

// ********** End of file ***************************************
