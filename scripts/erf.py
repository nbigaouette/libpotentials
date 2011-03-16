#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import math
import numpy as np
import scipy.special as special
import matplotlib.pyplot as plt

import on_key

sqrt2 = math.sqrt(2.0)
sqrtPi = math.sqrt(math.pi)

def lut_potential(x):
    return special.erf(x) / x

def lut_potential_expanded(x):
    return (
               2.0
            - (2.0 * x*x)                                                                               / 3.0
            +        x*x*x*x                                                                            / 5.0
            -        x*x*x*x*x*x                                                                        / 21.0
            +        x*x*x*x*x*x*x*x                                                                    / 108.0
            -        x*x*x*x*x*x*x*x*x*x                                                                / 660.0
            +        x*x*x*x*x*x*x*x*x*x*x*x                                                            / 4680.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                        / 37800.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                    / 342720.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                / 3447360.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                            / 38102400.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                        / 459043200.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                    / 5987520000.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                / 84064780800.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                            / 1264085222400.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                        / 20268952704000.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                    / 345226033152000.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                / 6224529991680000.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x            / 118443913555968000.0
            -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x        / 2372079457972224000.0
            +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x    / 49874491167621120000.0
    ) / sqrtPi

def lut_field(x):
    return (
        special.erf(x) / (2.0 * x**3)
        - (1.0/sqrtPi) * np.exp(-x**2)/x**2
    )

def lut_field_expanded(x):
    # http://www.wolframalpha.com/input/?i=expansion+erf%28x%29%2F%282*x**3%29-1%2Fsqrt%28pi%29*exp%28-x**2%29%2F%28x**2%29
    return (
           2.0                                                                                              / 3.0
        - (2.0 * x*x)                                                                                       / 5.0
        +        x*x*x*x                                                                                    / 7.0
        -        x*x*x*x*x*x                                                                                / 27.0
        +        x*x*x*x*x*x*x*x                                                                            / 132.0
        -        x*x*x*x*x*x*x*x*x*x                                                                        / 780.0
        +        x*x*x*x*x*x*x*x*x*x*x*x                                                                    / 5400.0
        -        x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                                / 42840.0
        +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                            / 383040.0
        -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                        / 3810240.0
        +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                    / 41731200.0
        -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                                / 498960000.0
        +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                            / 6466521600.0
        -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                        / 90291801600.0
        +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                    / 1351263513600.0
        -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                                / 21576627072000.0
        +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                            / 366148823040000.0
        -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                        / 6580217419776000.0
        +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                    / 124846287261696000.0
        -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x                / 2493724558381056000.0
        +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x            / 52307393175797760000.0
        -        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x        / 1149546198863462400000.0
        +        x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x*x    / 26414017102773780480000.0
        # +O(x^45)
    ) / sqrtPi


# ************************************************************************************************************************************

x = np.linspace(0.0, 6.0, 1000)

# ************************************************************************************************************************************
fig = on_key.figure()
axprops = dict()
axes_xmin  = 0.1
axes_sizex = 0.8
axes_sizey = 0.45
axes_ymin  = [0.5, 0.5 - axes_sizey]

ax1 = fig.add_axes([axes_xmin, axes_ymin[0], axes_sizex, axes_sizey], **axprops)
plt.plot(x, lut_field(x),               'b-',  label=r'$f(x) = \frac{\rm{erf}(x)}{x^{3}} - \frac{1}{\sqrt{\pi}} \frac{\rm{exp}(-x^2)}{x^{2}}$')
plt.plot(x, lut_field_expanded(x),      'b--', label='expansion')
plt.plot(x, 1.0/x**3,                   'b:',  label=r'$x^{-2}/x$')
plt.plot(x, lut_potential(x),           'r-',  label=r'$g(x) = \frac{\rm{erf}(x)}{x}$')
plt.plot(x, lut_potential_expanded(x),  'r--', label='expansion')
plt.plot(x, 1.0/x,                      'r:',  label=r'$x^{-1}$')
plt.grid()
plt.legend(loc='lower right')
axprops['sharex'] = ax1
ax1.set_xlim((x.min(), x.max()))
ax1.set_ylim((0.0, 1.2))
plt.setp(ax1.get_xticklabels(), visible=False)


ax2 = fig.add_axes([axes_xmin, axes_ymin[1], axes_sizex, axes_sizey], **axprops)
plt.semilogy(x, abs(lut_field(x)     - lut_field_expanded(x)),      'b-', label=r'$|f(x) - \rm{expansion}|$')
plt.semilogy(x, abs(lut_potential(x) - lut_potential_expanded(x)),  'r-',  label=r'$|g(x) - \rm{expansion}|$')
plt.grid()
plt.legend(loc='center left')
plt.xlabel("x")


plt.show()
