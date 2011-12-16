#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np

# Distance (in Bohr) where HS potential reaches the Coulombic values
HS_Xe_rmax = np.array([
#  Neutral   1+           2+                   3+                  4+                   5+                6+
    12.0, 3.20671, 3.4419260284993465, 2.5786383389900931, 2.3069998849855851, 2.1424931697942897, 1.9288708307182563
])

# Function parameters: -a*exp(-r*b+c) - d*exp(-r*e+f) - g*exp(-r*h+i)
HS_Xe_parameters = np.array([
#                   a                  b               c               d               e               f               g               h               i               j
                [3.14083382,        2.23690529,     2.45999159,     1.04922253,     0.758964055,    0.916259659,    6.43949225,     5.57500308,     3.46700894,     -2.22471387e-03],   # Neutral
                [6.68943368,        5.46613511,     3.50223818,     9.84854609,     16.33928308,    4.55862051,     2.05998631,     1.79149357,     2.67105113,     -0.21651473],       # 1+
                [1.08862692,        1.12845509,     2.40634711,     7.5231977,      11.80857511,    4.40029841,     3.5171341,      4.02105327,     3.70863489,     -0.33244088],       # 2+
                [1.26082951,        1.24346292,     2.48202614,     7.60391482,     14.20436211,    4.63976132,     5.75320941,     4.57482337,     3.45935956,     -0.55091234],       # 3+
                [8.33659368,        15.53383795,    4.69224278,     5.66740119,     4.93161199,     3.58851214,     1.33122023,     1.36069086,     2.65699251,     -0.90941801],       # 4+
                [8.13621709,        15.39455048,    4.6973397,      1.33881001,     1.40783802,     2.72036815,     5.60695758,     4.96351559,     3.59035494,     -1.33283627],       # 5+
                [7.52331956,        15.56584267,    4.77821787,     2.17218048,     1.51817071,     2.38100923,     5.09462365,     5.11830058,     3.70739486,     -1.84326541]])      # 6+

#print "HS_Xe_rmax.shape =", HS_Xe_rmax.shape
#print "HS_Xe_parameters.shape =", HS_Xe_parameters.shape

def HS_Fitting_Function_Xe_Potential(r, cs):
    a = HS_Xe_parameters[cs,0]
    b = HS_Xe_parameters[cs,1]
    c = HS_Xe_parameters[cs,2]
    d = HS_Xe_parameters[cs,3]
    e = HS_Xe_parameters[cs,4]
    f = HS_Xe_parameters[cs,5]
    g = HS_Xe_parameters[cs,6]
    h = HS_Xe_parameters[cs,7]
    i = HS_Xe_parameters[cs,8]
    j = HS_Xe_parameters[cs,9]
    potential = - a * np.exp(-r*b + c) - d * np.exp(-r*e + f) - g * np.exp(-r*h + i) + j
    return potential

def HS_Fitting_Function_Xe_Field(r, cs):
    a = HS_Xe_parameters[cs,0]
    b = HS_Xe_parameters[cs,1]
    c = HS_Xe_parameters[cs,2]
    d = HS_Xe_parameters[cs,3]
    e = HS_Xe_parameters[cs,4]
    f = HS_Xe_parameters[cs,5]
    g = HS_Xe_parameters[cs,6]
    h = HS_Xe_parameters[cs,7]
    i = HS_Xe_parameters[cs,8]
    #j = HS_Xe_parameters[cs,9]
    field = -(
          a * b * np.exp(-r*b + c)
        + d * e * np.exp(-r*e + f)
        + g * h * np.exp(-r*h + i))

    return field


def main():
    import matplotlib.pyplot as plt
    import on_key

    r = np.linspace(0.0, 10.0, 1000)

    css = []
    css.append(0)
    css.append(1)
    css.append(2)
    css.append(3)
    css.append(4)
    css.append(5)
    css.append(6)

    fig = on_key.figure()
    axprops = dict()
    ax1 = fig.add_subplot(211, **axprops)
    axprops['sharex'] = ax1
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = fig.add_subplot(212, **axprops)
    plt.subplots_adjust(hspace=0.0)
    ax1.grid(True)
    ax2.grid(True)

    for cs in css:
        HS_U = HS_Fitting_Function_Xe_Potential(r, cs)
        HS_E = HS_Fitting_Function_Xe_Field(r, cs)

        indices = np.where(r > HS_Xe_rmax[cs])
        HS_U[indices] = -float(cs) / r[indices]
        HS_E[indices] = -float(cs) / (r[indices]*r[indices])
        #HS_U[indices] = 0.0
        #HS_E[indices] = 0.0

        ax1.plot(r, HS_U, label = str(cs) + "+")
        ax2.plot(r, HS_E, label = str(cs) + "+")


    ax2.legend(loc='best')
    ax1.set_ylabel("Potentiel energy (Hartree)")
    ax2.set_ylabel("Electric field (atomic units)")
    ax2.set_xlabel("Distance (Bohr)")
    ax1.set_ylim((-5.0, 0.0))
    ax2.set_ylim((-0.6, 0.0))
    plt.show()


if __name__ == "__main__":
    main()

