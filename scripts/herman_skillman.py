#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np

max_hs_cs = 7

# Taken from libpotentials.git's Potential.cpp, variables fit_lt_R1, fit_lt_R2 and fit_lt_R3 (247a4be84f2b7e84504a723eeeb07b8c6c2c3537)
cutoffs = np.array([
        [0.000, 0.00, 0.00, 0.000], # Electron
        [0.073, 1.00, 3.00,  6.00], # Neutral
        [0.073, 1.00, 5.00, 12.00], # 1+
        [0.073, 0.60, 2.00,  6.00], # 2+
        [0.073, 0.75, 4.60, 12.00], # 3+
        [0.073, 0.75, 1.49,  2.00], # 4+
        [0.073, 0.75, 1.50,  1.49], # 5+
        [0.073, 0.75, 1.10,  1.10], # 6+
        [0.073, 0.35, 0.96,  0.96], # 7+
    ])

assert(max_hs_cs+2 == cutoff.shape[0])

def genericHSfit(par, r):
    # The 1/2 factor is becuase HS outputs the potential as 2V and that's how they were fit
    return -0.5*(
            - par[0] / (pow(r,par[5]) - par[1])
            - par[2] /  pow(r,par[4])
            + par[3] *  pow(r,par[6])
        )
#

def deriv_genericHSfit(par, r):
    # The 1/2 factor is becuase HS outputs the
    # potential as 2V(x) and that's how they were fit
    return 0.5*(
              par[0]*par[5]*pow(x,par[5]-1)/(pow(x,par[5])-par[1])/(pow(x,par[5])-par[1])
            + par[2]*par[4]/pow(x,par[4]+1)
            + par[6]*par[3]*pow(x,par[6]-1)
        )
#

hs_min_rad = 0.073

# ********** Herman-Skillman (HS) potential fit parameters *****
# The fit function is f(x)=-a/(x^n-b)-b/x^m +d*x^o and the fit
# parameters are in alphabetical order there is a different
# array for each radial distance, where the cutoff radial
# distance is the last element of the array.
fit_lt_R1 = [
    [-39.3117,          -0.23822,           1137.15,            1093.87,            0.926033,           1.35102,        -0.902534,  0.073, 1.0],
    [-50.9699,          -0.249349,          1190.93,            1137.92,            0.934615,           1.26191,        -0.915538,  0.073, 1.0],
    [-103318.0,         -0.0025085,         109884.0,           6808.38,            0.452195,           0.453217,       -0.452051,  0.073, 0.6],
    [-103309.0,         -0.00222854,        109893.0,           6799.44,            0.462727,           0.46361,        -0.462679,  0.073, 0.75],
    [-106539.0,         -0.00253429,        106375.0,           90.6738,            0.576876,           0.577581,       -1.31088,   0.073, 0.75],
    [-106552.0,         -0.00285323,        106363.0,           97.8488,            0.555262,           0.556083,       -1.26473,   0.073, 0.75],
    [-106572.0,         -0.00333483,        106342.0,           106.134,            0.523154,           0.524146,       -1.19394,   0.073, 0.75],
    [-156137.0,         -9.52875,           16279.0,            -30.664,            0.00362412,         4.27026,        -1.32769,   0.073, 0.35]
]

fit_lt_R2 = [
    [-106456.523613218,-0.00434541093317553, 106457.47661029,   449.688681389621,    1.05645523648719,    1.05644674944298, -2.10055725950707,  1.0,  3.0],
    [-103240.467920728,-0.000208924735834372,109961.532079643,  6730.47792027321,    0.935905881633714,   0.935947358356231,-0.93589486453368,  1.0,  5.0],
    [-7.43911046752643,-7.49680170563087,    83544.7086195816,  83531.3679985203,    2.4477467508823,     6.7666148229704,  -2.44780121116816,  0.6,  2.0],
    [-106458.718251124,-0.000545743677370998,106455.682016356,  42.3237633727421,    1.00559843304636,    1.00563829185779, -1.95708048412661,  0.75, 4.6],
    [-106453.495071328,-0.00399495917548577, 106460.925622402,  418.039392846222,    2.49073456323941,    2.4909137590075,  -4.9738722919108,   0.75, 1.49],
    [-106455.157115451,-0.00456229507833856, 106460.145184005,  475.327392665337,    2.27725233310332,    2.27744553473598, -4.5488537847976,   0.75, 1.5],
    [-106452.180354907,-0.00278122155186969, 106461.721043604,  291.588888724572,    3.59580873362151,    3.59615037646864, -7.17813960396325,  0.75, 1.1],
    [-156436.219173519,-13.360177523064,     10907.4536590735,  -0.0178811573295934, 0.0295757482829108,  0.398808602998421,-5.49402342863045,  0.35, 0.96]
]

fit_lt_R3 = [
    [-106344.499357271,-0.0870404156519758,  106379.969770542, 8916.02780769541,    2.34571347967461,   2.34558512875328,   -4.64724093315347,  3.0,    6.0],
    [-103237.178865962,-6.19966863330973e-05,109964.821133342, 6727.38883676891,    0.990416309150688,  0.990415990770504,  -0.990490798848876, 5.0,    12.0],
    [-106453.321357016,-0.0233720244005975,  106447.424341854, 2423.61663166259,    1.69020647850117,   1.69030805063035,   -3.36829845029172,  2.0,    6.0],
    [-106457.189833221,-0.000453936408454839,106457.21043453,  42.3245989284602,    0.499881464943715,  0.499881437435555,  -0.999349099198404, 4.6,    12.0],
    [-106478.807529316,-0.00470475292274558, 106435.613140363, 443.194839747241,    0.192492936878364,  0.1932366392085,    -0.192465568481317, 1.49,   2.0],
    [-106455.157115451,-0.00456229507833856, 106460.145184005, 475.327392665337,    2.27725233310332,   2.27744553473598,   -4.5488537847976,   1.49,   1.49],
    [-106452.180354907,-0.00278122155186969, 106461.721043604, 291.588888724572,    3.59580873362151,   3.59615037646864,  -7.17813960396325,   1.1,    1.1],
    [-156436.219173519,-13.360177523064,     10907.4536590735, -0.0178811573295934, 0.0295757482829108, 0.398808602998421, -5.49402342863045,   0.96,   0.96]
]

assert(max_hs_cs+1 == shape(fit_lt_R1)[0])
assert(max_hs_cs+1 == shape(fit_lt_R2)[0])
assert(max_hs_cs+1 == shape(fit_lt_R3)[0])

def HS(r, cs):

    if (cs == 0):

        if   (r >= fit_lt_R3[0][8]):
            # Potential outside the electron cloud goes to 0
            # exponentially using f(x)=h*exp(-v*x+k)
            phi12 = 1.93775072943628 * np.exp(-0.533297816151*r -0.7486357665822807)
        elif (r <  fit_lt_R3[0][8] and r >= fit_lt_R3[0][7]):
            # R3
            phi12 = genericHSfit(fit_lt_R3[0], r)
        elif (r <  fit_lt_R2[0][8] and r >= fit_lt_R2[0][7]):
            # R2
            phi12 = genericHSfit(fit_lt_R2[0], r)
        elif (r <  fit_lt_R1[0][8] and r >= fit_lt_R1[0][7]):
            # R1
            phi12 = genericHSfit(fit_lt_R1[0],r)
        else:
            # Hard cutoff
            phi12 = genericHSfit(fit_lt_R1[0], hs_min_rad)
    elif ( (cs < 10) and (cs > 0) ):
        # If charge state is between 0 and 9 (inclusive)...
        # Constant potential
        phi12 = np.ones(len(r)) * genericHSfit(fit_lt_R1[cs], hs_min_rad)
        # Coulomb
        # Outside electron cloud: Coulombic pot.
        i = np.nonzero(r >= fit_lt_R3[cs][8])
        phi12[i] = float(cs) / r[i]
        # R3
        i = np.nonzero((r <  fit_lt_R3[cs][8]) & (r >= fit_lt_R3[cs][7]))
        phi12[i] = genericHSfit(fit_lt_R3[cs], r[i])
        # R2
        i = np.nonzero((r <  fit_lt_R2[cs][8]) & (r >= fit_lt_R2[cs][7]))
        phi12[i] = genericHSfit(fit_lt_R2[cs], r[i])
        # R1
        i = np.nonzero((r <  fit_lt_R1[cs][8]) & (r >= fit_lt_R1[cs][7]))
        phi12[i] = genericHSfit(fit_lt_R1[cs], r[i])

    else:
        phi12 = 0.0

    return phi12



