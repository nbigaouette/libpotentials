
#include <cmath>

#include <Assert.hpp>

#include "HermanSkillman.hpp"

const int HS_Xe_MaxNbCS = 7;

// Distance (in Bohr) where HS potential reaches the Coulombic values
const double HS_Xe_rmax[HS_Xe_MaxNbCS] = {
    12.0, 3.20671, 3.4419260284993465, 2.5786383389900931, 2.3069998849855851, 2.1424931697942897, 1.9288708307182563
};

const double HS_Xe_parameters[HS_Xe_MaxNbCS][10] = {
    {3.14083382e+00,    2.23690529e+00, 2.45999159e+00, 1.04922253e+00, 7.58964055e-01, 9.16259659e-01, 6.43949225e+00, 5.57500308e+00, 3.46700894e+00, -2.22471387e-03},
    {6.68943368,        5.46613511,     3.50223818,     9.84854609,     16.33928308,    4.55862051,     2.05998631,     1.79149357,     2.67105113,     -0.21651473},
    {1.08862692,        1.12845509,     2.40634711,     7.5231977,      11.80857511,    4.40029841,     3.5171341,      4.02105327,     3.70863489,     -0.33244088},
    {1.26082951,        1.24346292,     2.48202614,     7.60391482,     14.20436211,    4.63976132,     5.75320941,     4.57482337,     3.45935956,     -0.55091234},
    {8.33659368,        15.53383795,    4.69224278,     5.66740119,     4.93161199,     3.58851214,     1.33122023,     1.36069086,     2.65699251,     -0.90941801},
    {8.13621709,        15.39455048,    4.6973397,      1.33881001,     1.40783802,     2.72036815,     5.60695758,     4.96351559,     3.59035494,     -1.33283627},
    {7.52331956,        15.56584267,    4.77821787,     2.17218048,     1.51817071,     2.38100923,     5.09462365,     5.11830058,     3.70739486,     -1.84326541},
};

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

void Set_HermanSkillman_Lookup_Tables_Xe(std::vector<LookUpTable<fdouble> > &lut_pot,
                                         std::vector<LookUpTable<fdouble> > &lut_field)
{
    const int lut_n = 10000;

    lut_pot.resize(HS_Xe_MaxNbCS);
    lut_field.resize(HS_Xe_MaxNbCS);

    // Set potential lookup table
    for (int cs_i = 0 ; cs_i < HS_Xe_MaxNbCS ; cs_i++)
    {
        const int cs = int(cs_i) - 1;
        const double xmin = 0.0;
        const double xmax = 1.1 * HS_Xe_rmax[cs_i];
        lut_pot[cs_i].Initialize(NULL,   fdouble(xmin), fdouble(xmax), lut_n, "Initialize_HS() LookUpTable (lut_pot, cs=" + IntToStr(cs) + ")");
        lut_field[cs_i].Initialize(NULL, fdouble(xmin), fdouble(xmax), lut_n, "Initialize_HS() LookUpTable (lut_field, cs=" + IntToStr(cs) + ")");

        // FIXME: Dynamically choose between atom types for HS

        // Now populate the lookup tables.
        double r;
        for (int i = 0 ; i <= lut_n ; i++)
        {
            r = lut_pot[cs_i].Get_x_from_i(i);

            double HS_U_r = 0.0;
            double HS_E_r = 0.0;

            if (r >= HS_Xe_rmax[cs_i])
            {
                // Use Coulomb
                HS_U_r = -double(cs) / r;
                HS_E_r = double(cs) / (r*r);
            }
            else
            {
//                 if (cs == 0)
//                 {
//                     // FIXME: We set the potential of a neutral to be 0.
//                     //        In real HS, neutrals do have a potential.
//                     HS_U_r = 0.0;
//                     HS_E_r = 0.0;
//                 }
//                 else
                {
                    // Use the 1+ for the electron by taking the absolute value of the charge state.
                    // Electron has same potential as 1+. This is necessary to have a 1+ and an electron
                    // on top of each other be seen by other particles as a neutral.
                    HS_U_r = HS_Fitting_Function_Xe_Potential(r, std::abs(cs));
                    HS_E_r = HS_Fitting_Function_Xe_Field(    r, std::abs(cs));
                }
            }

            // Scale with the charge state
            if (cs != 0)
            {
                HS_U_r /= double(cs);
                HS_E_r /= double(cs);
            }

            Assert_isinf_isnan(HS_U_r);

            // Save it
            lut_pot[cs_i].Set(  i, fdouble(HS_U_r));
            lut_field[cs_i].Set(i, fdouble(HS_E_r));
        }
    }
}

