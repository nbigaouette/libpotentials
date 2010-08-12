#ifndef INC_POTENTIALS_HPP
#define INC_POTENTIALS_HPP

#include "General.hpp"
#include "Structure_Potentials.hpp"

void Initialize_Simple(const double &minr);
void Initialize_SuperGaussian(const int &m);
void Initialize_HS(const int &input_sg_m, const double &min_rad);

void Get_r21(double r1[3], double r2[3], double r21[3]);
double Get_Distance_Squared(double r1[3], double r2[3]);
double Get_Distance(double r1[3], double r2[3]);
void set_vector_between_particles(
        double r1[3], double r2[3],
        double r21[3], double &r212, double &r, double &one_over_r);

double Coulomb_Potential(const double kQ, const double r);
void   Coulomb_Field(const double phi12, double E[3], const double dr[3],
    const double dr2);

double tmp_get_shieldr_2(const int chg_st_1, const int chg_st_2);
double tmp_get_shieldr(const int chg_st, const char *message);

// **************************************************************
// ********** Function pointers for... **************************

// ...setting the parameters of the potential/field calculation
extern void   (*Potentials_Set_Parameters)(void *p1, void *p2, potential_paramaters &potparams);
// ...calculating the potential
extern double (*Calculate_Potential)(      void *p1, void *p2, potential_paramaters &potparams);
// ...setting the electric field
extern void   (*Set_Field)(                void *p1, void *p2, potential_paramaters &potparams, double &phi, double E[3]);


// **************************************************************
// ********** Real functions ************************************

// ********** Simple pontential *********************************
void Potentials_Set_Parameters_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_Simple(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Harmonic pontential *******************************
void Potentials_Set_Parameters_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_Harmonic(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Super Gaussian ************************************
void Potentials_Set_Parameters_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Herman-Skillman (HS) ******************************
void Potentials_Set_Parameters_HS_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_HS_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_HS_SuperGaussian(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Gaussian Distribution parameters ******************
void Potentials_Set_Parameters_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_GaussianDistribution(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Symmetric Charge distribution *********************
void Potentials_Set_Parameters_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_ChargeDistribution_Symmetric(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Screened Coulomb **********************************
void Potentials_Set_Parameters_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_Cutoff_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_Cutoff_ScreenedCoulomb(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// ********** Pseudo-particles (QFD) ****************************
void Potentials_Set_Parameters_PseudoParticles(
    void *p1, void *p2,
    potential_paramaters &potparams);
double Calculate_Potential_PseudoParticles(
    void *p1, void *p2,
    potential_paramaters &potparams);
void Set_Field_PseudoParticles(
    void *p1, void *p2,
    potential_paramaters &potparams,
    double &phi, double E[3]);

// **************************************************************
// **************************************************************


// ********** Herman-Skillman (HS) potential fit functions ******
double deriv_genericHSfit(const double *par, double x);
double genericHSfit(const double *par, double x);

// ********** Herman-Skillman (HS) potential fit parameters *****
// The fit function is f(x)=-a/(x^n-b)-b/x^m +d*x^o and the fit
// parameters are in alphabetical order there is a different
// array for each radial distance, where the cutoff radial
// distance is the last element of the array.
const double fit_lt_R1[10][9]={
    {-39.3117,-0.23822 ,1137.15,1093.87,0.926033 ,1.35102,-0.902534,0.073,1.0},
    {-50.9699,-0.249349,1190.93,1137.92,0.934615,1.26191,-0.915538,0.073,1.0},
    {-103318.0,-0.0025085,109884.0,6808.38,0.452195,0.453217,-0.452051,0.073,0.6},
    {-103309.0,-0.00222854,109893.0,6799.44,0.462727,0.46361,-0.462679,0.073,0.75},
    {-106539.0,-0.00253429,106375.0,90.6738,0.576876 ,0.577581,-1.31088,0.073,0.75},
    {-106552.0,-0.00285323,106363.0,97.8488,0.555262,0.556083,-1.26473,0.073,0.75},
    {-106572.0,-0.00333483,106342.0,106.134,0.523154,0.524146,-1.19394,0.073,0.75},
    {-156137.0, -9.52875,16279.0,-30.664,0.00362412,4.27026,-1.32769,0.073,0.35},
    {506.586882065287, -5.70042070187032, -44.6577028237441,-1.00674477808791, -58.7532076356085, -1.00659228186075,0.0,0.02},
    {525.638029177344, -4.17758373083591, -44.8158100544942, -1.00624596615361, -58.9113171412815,-1.00622039392236,0.0,0.02}
};
const double fit_lt_R2[10][9]={
    {-106456.523613218,-0.00434541093317553,106457.47661029,
        449.688681389621,1.05645523648719,1.05644674944298,-2.10055725950707,
        1.0,3.0},
    {-103240.467920728,-0.000208924735834372,109961.532079643,
        6730.47792027321,0.935905881633714,0.935947358356231,-0.93589486453368,
        1.0,5.0},
    {-7.43911046752643,-7.49680170563087,83544.7086195816,
        83531.3679985203,2.4477467508823,6.7666148229704,-2.44780121116816,
        0.6,2.0},
    {-106458.718251124,-0.000545743677370998,106455.682016356,
        42.3237633727421,1.00559843304636,1.00563829185779,-1.95708048412661,
        0.75,4.6},
    {-106453.495071328,-0.00399495917548577,106460.925622402,
        418.039392846222,2.49073456323941,2.4909137590075,-4.9738722919108,
        0.75,1.49},
    {-106455.157115451,-0.00456229507833856,106460.145184005,
        475.327392665337,2.27725233310332,2.27744553473598,-4.5488537847976,
        0.75,1.5},
    {-106452.180354907,-0.00278122155186969,106461.721043604,
        291.588888724572,3.59580873362151,3.59615037646864,-7.17813960396325,
        0.75,1.1},
    {-156436.219173519,-13.360177523064,10907.4536590735,
        -0.0178811573295934,0.0295757482829108,0.398808602998421,-5.49402342863045,
        0.35,0.96},
    {290.185723442269, 2.47652109557943, -88.9477236879255, -1.03429895973465,  111.942600145551,  0.285566497558354,0.02,0.2},
    {286.879451804095, 8.00573092007186,  -33060.3534705734, -0.760161832369098,  33083.8629659416, -0.757992817776906,0.02,0.2}
};

const double fit_lt_R3[10][9]={
    {-106344.499357271,-0.0870404156519758,106379.969770542,
        8916.02780769541,2.34571347967461,2.34558512875328,-4.64724093315347,
        3.0,6.0},
    {-103237.178865962,-6.19966863330973e-05,109964.821133342,
        6727.38883676891,0.990416309150688,0.990415990770504,-0.990490798848876,
        5.0,12.0},
    {-106453.321357016,-0.0233720244005975,106447.424341854,
        2423.61663166259,1.69020647850117,1.69030805063035,-3.36829845029172,
        2.0,6.0},
    {-106457.189833221,-0.000453936408454839,106457.21043453,
        42.3245989284602,0.499881464943715,0.499881437435555,-0.999349099198404,
        4.6,12.0},
    {-106478.807529316,-0.00470475292274558,106435.613140363,
        443.194839747241,0.192492936878364,0.1932366392085,-0.192465568481317,
        1.49,2.0},
    {-106455.157115451,-0.00456229507833856,106460.145184005,
        475.327392665337,2.27725233310332,2.27744553473598,-4.5488537847976,
        1.49,1.49},
    {-106452.180354907,-0.00278122155186969,106461.721043604,
        291.588888724572,3.59580873362151,3.59615037646864,-7.17813960396325,
        1.1,1.1},
    {-156436.219173519,-13.360177523064,10907.4536590735,
        -0.0178811573295934,0.0295757482829108,0.398808602998421,-5.49402342863045,
        0.96,0.96},
    {301.086976771718, 3.09640164288097, -33085.4767378537, -0.704929036919971, 33058.7587426746,  -0.702937583476305,0.2,0.9},
    {273.866594715728, 5.00998369874365, -33079.2658494451, -0.969158166151665,  33064.9547364488, -0.968029709305527,0.2,0.86}
};

#endif // #ifndef INC_POTENTIALS_HPP

// ********** End of file ***************************************
