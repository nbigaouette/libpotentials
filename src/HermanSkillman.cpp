
#include <cmath>

#include <Assert.hpp>

#include "HermanSkillman.hpp"
#include "Constants.hpp"

std::vector<fdouble> hs_min_rad;
std::vector<LookUpTable<fdouble> > hs_lut_potential;
std::vector<LookUpTable<fdouble> > hs_lut_field;

// Distance (in Bohr) where HS potential reaches the Coulombic values
const double HS_Xe_rmax[HS_Xe_MaxNbCS] = {
//  El.     Neutral   1+           2+                   3+                  4+                   5+                6+
    3.20671, 12.0, 3.20671, 3.4419260284993465, 2.5786383389900931, 2.3069998849855851, 2.1424931697942897, 1.9288708307182563
};

// Function parameters: -a*exp(-r*b+c) - d*exp(-r*e+f) - g*exp(-r*h+i)
const double HS_Xe_parameters[HS_Xe_MaxNbCS][10] = {
//                   a                  b               c               d               e               f               g               h               i               j
/* Electron */  {6.68943368,        5.46613511,     3.50223818,     9.84854609,     16.33928308,    4.55862051,     2.05998631,     1.79149357,     2.67105113,     -0.21651473},
/* Neutral  */  {3.14083382,        2.23690529,     2.45999159,     1.04922253,     0.758964055,    0.916259659,    6.43949225,     5.57500308,     3.46700894,     -2.22471387e-03},
/* 1+ */        {6.68943368,        5.46613511,     3.50223818,     9.84854609,     16.33928308,    4.55862051,     2.05998631,     1.79149357,     2.67105113,     -0.21651473},
/* 2+ */        {1.08862692,        1.12845509,     2.40634711,     7.5231977,      11.80857511,    4.40029841,     3.5171341,      4.02105327,     3.70863489,     -0.33244088},
/* 3+ */        {1.26082951,        1.24346292,     2.48202614,     7.60391482,     14.20436211,    4.63976132,     5.75320941,     4.57482337,     3.45935956,     -0.55091234},
/* 4+ */        {8.33659368,        15.53383795,    4.69224278,     5.66740119,     4.93161199,     3.58851214,     1.33122023,     1.36069086,     2.65699251,     -0.90941801},
/* 5+ */        {8.13621709,        15.39455048,    4.6973397,      1.33881001,     1.40783802,     2.72036815,     5.60695758,     4.96351559,     3.59035494,     -1.33283627},
/* 6+ */        {7.52331956,        15.56584267,    4.77821787,     2.17218048,     1.51817071,     2.38100923,     5.09462365,     5.11830058,     3.70739486,     -1.84326541},
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
    // cs_i: 0 == electron, 1 == neutral, 2 == 1+, etc.
    for (int cs_i = 0 ; cs_i < HS_Xe_MaxNbCS ; cs_i++)
    {
        const int cs = cs_i - 1;

        // Distance range for the HS.
        const double xmin = 0.0;        // [Bohr]
        double xmax = HS_Xe_rmax[cs_i]; // [Bohr]
        if (cs != 0)
            xmax *= 2.0;
        lut_pot[cs_i].Initialize(  NULL, fdouble(xmin), fdouble(xmax), lut_n, "Initialize_HS() LookUpTable (lut_pot, cs=" + IntToStr(cs) + ")");
        lut_field[cs_i].Initialize(NULL, fdouble(xmin), fdouble(xmax), lut_n, "Initialize_HS() LookUpTable (lut_field, cs=" + IntToStr(cs) + ")");

        // FIXME: Dynamically choose between atom types for HS
        // Start: use function pointers
        double (*HS_Fitting_Function_Potential)(const double r, const int cs);
        double (*HS_Fitting_Function_Field)(const double r, const int cs);

        HS_Fitting_Function_Potential = &HS_Fitting_Function_Xe_Potential;
        HS_Fitting_Function_Field     = &HS_Fitting_Function_Xe_Field;

        // Populate the lookup tables.
        double r;
        for (int i = 0 ; i <= lut_n ; i++)
        {
            r = lut_pot[cs_i].Get_x_from_i(i);  // [Bohr]

            double HS_U_r = 0.0;
            double HS_E_r = 0.0;

            if (r >= HS_Xe_rmax[cs_i])
            {
                // Use Coulomb
                HS_U_r = -double(std::abs(cs)) / r;
            }
            else
            {
                // Use the 1+ for the electron by taking the absolute value of the charge state.
                // Electron has same potential as 1+. This is necessary to have a 1+ and an electron
                // on top of each other be seen by other particles as a neutral.
                HS_U_r = HS_Fitting_Function_Potential(r, cs_i);
            }

            // The fitting is performed on the potential, not the field. By doing the same test
            // over "r" for the field as for the potential, a discontinuity always show up in the field
            // close to where the potential becomes Coulombic. What is important is the field: it must
            // be smooth. The fit will cross the Coulomb field at some point, a bit farther then
            // the HS_Xe_rmax distance. Detect that crossing by taking the minumum value between
            // the Coulomb field and the fit. This enforce a continuous field.
            if (cs == 0)
            {
                if (r < HS_Xe_rmax[cs_i])
                    HS_E_r = HS_Fitting_Function_Field(r, cs_i);
            }
            else
            {
                HS_E_r = std::min(
                    -double(std::abs(cs)) / (r*r),
                    HS_Fitting_Function_Field(r, cs_i)
                );
            }

            // Scale with the charge state
            if (cs != 0)
            {
                HS_U_r /= double(cs);
                HS_E_r /= double(cs);
            }
//             else
//             {
//                 // FIXME: We set the potential of a neutral to be 0.
//                 //        In real HS, neutrals do have a potential.
//                 HS_U_r = 0.0;
//                 HS_E_r = 0.0;
//             }

            // We store E/r, not E
            HS_E_r /= r;

            Assert_isinf_isnan(HS_U_r);

            // Save it
            lut_pot[cs_i].Set(  i, fdouble(HS_U_r));
            lut_field[cs_i].Set(i, fdouble(HS_E_r));
        }
    }
}

// **************************************************************
void Initialize_HS(const fdouble &base_potential)
/**
 * Initialize Herman-Skillman potential
 * @param   base_potential  Potential depth wanted [eV]
 */
{
    // We'll need one lookup table per charge state
    std_cout << "FIXME: Dynamically choose between atom types for HS (" << __FILE__ << ", line " << __LINE__ << ")\n";
    // LUTs stored in atomic units
    Set_HermanSkillman_Lookup_Tables_Xe(hs_lut_potential, hs_lut_field);

    hs_min_rad.resize(hs_lut_potential.size());

    // Find the radius where the HS potential is equal to "base_potential"
    // by doing a bisection, for all supported charge states.
    fdouble r_left, r_right, found_r, pot; // [Bohr]
    for (int cs_i = 0 ; cs_i < int(hs_lut_potential.size()) ; cs_i++)
    {
        const int cs = cs_i - 1;

        // Initial conditions
        hs_min_rad[cs_i] = hs_lut_potential[cs_i].Get_x_from_i(0);
        r_left  = hs_min_rad[cs_i];     // Minimum radius of fit
        r_right = hs_lut_potential[cs_i].Get_XMax();

        // Start bisection!
        // See http://en.wikipedia.org/wiki/Bisection_method#Practical_considerations
        found_r = r_right + (r_left - r_right) / fdouble(2.0);

        while (std::abs(found_r - r_left) > 1.0e-100 && std::abs(found_r - r_right) > 1.0e-100)
        {
            pot = std::abs(hs_lut_potential[cs_i].read(found_r));
            //printf("base_potential = %10.5g   r_left = %10.5g   r = %10.5g   r_right = %10.5g   HS(r) = %10.5g\n", base_potential, r_left, found_r, r_right, pot);
            Assert_isinf_isnan(pot);
            if (pot <= base_potential)
            {
                r_right = found_r;
            }
            else
            {
                r_left = found_r;
            }
            found_r = r_right + (r_left - r_right) / fdouble(2.0);
        }
        std_cout << "Bisection end: cs = " << cs << "  HS(r="<<found_r<<") = " << (cs == 0 ? 1 : fdouble(std::abs(cs)))*pot << "\n";

        Assert_isinf_isnan(found_r);
        Assert_isinf_isnan(pot);
        assert(found_r > 0.0);
        //assert(std::abs(pot - std::max(1,cs)*base_potential) < 1.0e-3);

        if (found_r <= 0.99999*hs_min_rad[cs_i])
        {
            std_cout << "##############################################\n";
            DEBUGP("Initialize_HS() called with a potential depth too deep.\n");
            std_cout << "The value found " << found_r << " Bohr\n";
            std_cout << "for charge state " << cs << " should not be lower than " << hs_min_rad[cs_i] << " Bohr\n";
            std_cout << "Potential depth wanted: " << base_potential << " eV (" << base_potential*libpotentials::eV_to_Eh << " Eh)\n";
            std_cout << "Exiting\n";
            abort();
        }

        hs_min_rad[cs_i] = found_r;
    }

    // Now that the cutting radius is found for each charge states, change lookup tables values
    for (int cs_i = 0 ; cs_i < int(hs_lut_potential.size()) ; cs_i++)
    {
        const int lut_n = hs_lut_potential[cs_i].Get_n();
        for (int i = 0 ; i <= lut_n ; i++)
        {
            if (hs_lut_potential[cs_i].Table(i) < -base_potential*libpotentials::eV_to_Eh)
            {
                hs_lut_potential[cs_i].Set(i, -base_potential*libpotentials::eV_to_Eh);
                hs_lut_field[cs_i].Set(i, 0.0);
            }
            if (hs_lut_potential[cs_i].Table(i) > base_potential*libpotentials::eV_to_Eh)
            {
                hs_lut_potential[cs_i].Set(i, base_potential*libpotentials::eV_to_Eh);
                hs_lut_field[cs_i].Set(i, 0.0);
            }
        }
    }

    /*
    // Print lookup table for verification
    const int max_lut = 8;
    std::string filename("lut_hs.dat");
    std::string gnuplot_command("");
    gnuplot_command += "#set term wxt 3; plot ";
    int row = 1;
    for (int l = 0 ; l < max_lut ; l++)
    {
        gnuplot_command += "\"" + filename + "\" using " + IntToStr(row++) + ":";
        gnuplot_command += IntToStr(row++) + "  title \"LUT(V(" + IntToStr(l-1) + ")) (HS)\" with lines lw 3";
        if (l == max_lut-1)
            gnuplot_command += "\n";
        else
            gnuplot_command += ", ";
    }
    gnuplot_command += "#set term wxt 4; plot ";
    for (int l = 0 ; l < max_lut ; l++)
    {
        gnuplot_command += "\"" + filename + "\" using " + IntToStr(row++) + ":";
        gnuplot_command += IntToStr(row++) + "  title \"cs*LUT(V(" + IntToStr(l-1) + ")) (HS)\" with lines lw 3, ";
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
        gnuplot_command += IntToStr(row++) + "  title \"LUT(E/r(" + IntToStr(l-1) + ")) (HS)\" with lines lw 3";
        if (l == max_lut-1)
            gnuplot_command += "\n";
        else
            gnuplot_command += ", ";
    }
    gnuplot_command += "#set term wxt 6; plot ";
    for (int l = 0 ; l < max_lut ; l++)
    {
        gnuplot_command += "\"" + filename + "\" using " + IntToStr(row++) + ":";
        gnuplot_command += IntToStr(row++) + "  title \"r*cs*LUT(E/r(" + IntToStr(l-1) + ")) (HS)\" with lines lw 3, ";
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
        fprintf(stderr, "%20s ", std::string("V(" + IntToStr(l-1) + ")").c_str());
    }
    for (int l = 0 ; l < max_lut ; l++)
    {
        fprintf(stderr, "%20s ", "r (bohr)");
        fprintf(stderr, "%20s ", std::string("cs*V(" + IntToStr(l-1) + ")").c_str());
    }
    for (int l = 0 ; l < max_lut ; l++)
    {
        fprintf(stderr, "%20s ", "r (bohr)");
        fprintf(stderr, "%20s ", std::string("E/r(" + IntToStr(l-1) + ")").c_str());
    }
    for (int l = 0 ; l < max_lut ; l++)
    {
        fprintf(stderr, "%20s ", "r (bohr)");
        fprintf(stderr, "%20s ", std::string("cs*E(" + IntToStr(l-1) + ")").c_str());
    }
    fprintf(stderr, "\n");

    const int lut_n      = hs_lut_potential[0].Get_n();
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

