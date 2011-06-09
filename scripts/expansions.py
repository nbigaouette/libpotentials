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


def lut_potential_log_expanded(x):
    # http://www.wolframalpha.com/input/?i=erf%28exp%28a%29%29+%2F+%28exp%28a%29%29
    a = np.log(x)
    e = math.exp(1.0)
    erf1 = special.erf(1.0)
    return (
                                erf1
        +               a     * (2.0 / (e * sqrtPi)-erf1)
        +               a**2  * (erf1/2.0 - 3.0/(e * sqrtPi))
        +               a**3  * (1.0/(e * sqrtPi)-(erf1)/6.0)
        + 1.0/24.0  *   a**4  * (erf1+10.0/(e * sqrtPi))
        + 1.0/120.0 *   a**5  * (22.0/(e * sqrtPi)-erf1)
        + 1.0/720.0 *   a**6  * (erf1-150.0/(e * sqrtPi))
        -               a**7  * (e*erf1+1002.0/sqrtPi)/(5040.0 * e)
        +               a**8  * (erf1-1302.0/(e * sqrtPi))/40320.0
        +               a**9  * (26902.0/(e * sqrtPi)-erf1)/362880.0
        +               a**10 * (e*erf1+246506.0/sqrtPi)/(3628800.0 * e)
        +               a**11 * (599318.0/sqrtPi - e*erf1)/(39916800.0 * e)
        -               a**12 * (9528598.0/sqrtPi - e*erf1)/(479001600.0 * e)
        -               a**13 * (e*erf1+135723754.0/sqrtPi)/(6227020800.0 * e)
        -               a**14 * (692208918.0/sqrtPi - e*erf1)/(87178291200.0 * e)
        +               a**15 * (4302456086.0/sqrtPi - e*erf1)/(1307674368000.0 * e)
        +               a**16 * (e*erf1+124593552106.0/sqrtPi)/(20922789888000.0 * e)
        +               a**17 * (1178087671062.0/sqrtPi - e*erf1)/(355687428096000.0 * e)
        +               a**18 * (e*erf1+1086500420330.0/sqrtPi)/(6402373705728000.0 * e)
        -               a**19 * (e*erf1+147087633201898.0/sqrtPi)/(121645100408832000.0 * e)
        -               a**20 * (2517298379179286.0/sqrtPi - e*erf1)/(2432902008176640000.0 * e)
        -               a**21 * (e*erf1+18070149954321130.0/sqrtPi)/(51090942171709440000.0 * e)
        +               a**22 * (e*erf1+132038987875400426.0/sqrtPi)/(1124000727777607680000.0 * e)
        +               a**23 * (5953921562136495382.0/sqrtPi - e*erf1)/(25852016738884976640000.0 * e)
        +               a**24 * (e*erf1+87873383278392163050.0/sqrtPi)/(620448401733239439360000.0 * e)
        +               a**25 * (443078107750850315542.0/sqrtPi - e*erf1)/(15511210043330985984000000.0 * e)
        -               a**26 * (11743714906826509075734.0/sqrtPi - e*erf1)/(403291461126605635584000000.0 * e)
        -               a**27 * (e*erf1+370523876199881176562410.0/sqrtPi)/(10888869450418352160768000000.0 * e)
        -               a**28 * (5218119697123425490257174.0/sqrtPi - e*erf1)/(304888344611713860501504000000.0 * e)
        -               a**29 * (e*erf1+15610098398647240048556778.0/sqrtPi)/(8841761993739701954543616000000.0 * e)
        +               a**30 * (e*erf1+1220358267742920850482969322.0/sqrtPi)/(265252859812191058636308480000000.0 * e)
        +               a**31 * (35791740664812835109132518678.0/sqrtPi - e*erf1)/(8222838654177922817725562880000000.0 * e)
        +               a**32 * (e * erf1+512227721799073043023911762666/sqrtPi)/(263130836933693530167218012160000000.0 * e)
        +               a**33 * (805820224276882734389243364630/sqrtPi - e*erf1)/(8683317618811886495518194401280000000.0 * e)
        -               a**34 * (172352783672515496324431963178262/sqrtPi - e*erf1)/(295232799039604140847618609643520000000.0 * e)
        -               a**35 * (e * erf1+5218497864994500452628855878892266/sqrtPi)/(10333147966386144929666651337523200000000.0 * e)
        -               a**36 * (80273838077660983402587126114501910/sqrtPi - e*erf1)/(371993326789901217467999448150835200000000.0 * e)
        -               a**37 * (e * erf1+119597113017073805379403882993597162/sqrtPi)/(13763753091226345046315979581580902400000000.0 * e)
        +               a**38 * (e * erf1+33167031420062725477719844729643907818/sqrtPi)/(523022617466601111760007224100074291200000000.0 * e)
        +               a**39 * (1105452898018999810213082930577589355798/sqrtPi - e*erf1)/(20397882081197443358640281739902897356800000000.0 * e)
        +               a**40 * (e * erf1+19194819359843402980378759970280384281322/sqrtPi)/(815915283247897734345611269596115894272000000000.0 * e)
        +               a**41 * (59528695268406316842612281655736656676118/sqrtPi - e*erf1)/(33452526613163807108170062053440751665152000000000.0 * e)
        -               a**42 * (8420623623139728132573770040775432651883798/sqrtPi - e*erf1)/(1405006117752879898543142606244511569936384000000000.0 * e)
        -               a**43 * (e * erf1+326300612751223971102684361133439100841409258/sqrtPi)/(60415263063373835637355132068513997507264512000000000.0 * e)
        -               a**44 * (6636648699735983413339758349915275122587944214/sqrtPi - e*erf1)/(2658271574788448768043625811014615890319638528000000000.0 * e)
        -               a**45 * (e * erf1+41454073463196484471219345470169757674258873066/sqrtPi)/(119622220865480194561963161495657715064383733760000000000.0 * e)
        +               a**46 * (e * erf1+2674222695660923324574048214427267778419633470186/sqrtPi)/(5502622159812088949850305428800254892961651752960000000000.0 * e)
        # +O(a^47)
    )

# ************************************************************************************************************************************

x = np.linspace(0.0, 6.0, 1000)

#fig = on_key.figure()
#ax1 = plt.subplot(111)
##plt.plot(x, lut_potential_log_expanded(np.log(x)))
#a = np.log(x)
#plt.plot(x, special.erf(np.exp(a))/np.exp(a))
#plt.plot(x, lut_potential_log_expanded(x))
#ax1.set_ylim((0.0, 1.2))
#plt.show()
#sys.exit(0)

# ************************************************************************************************************************************
fig = on_key.figure()
axprops = dict()
axes_xmin  = 0.1
axes_sizex = 0.8
axes_sizey = 0.42
axes_ymin  = [0.5, 0.5 - axes_sizey]

ax1 = fig.add_axes([axes_xmin, axes_ymin[0], axes_sizex, axes_sizey], **axprops)
plt.plot(x, lut_field(x),               'b-',  label=r'$F(x) = \frac{\rm{erf}(x)}{x^{3}} - \frac{1}{\sqrt{\pi}} \frac{\rm{exp}(-x^2)}{x^{2}}$')
plt.plot(x, lut_field_expanded(x),      'b--', label='expansion')
plt.plot(x, 1.0/x**3,                   'b:',  label=r'$x^{-2}/x$')
plt.plot(x, lut_potential(x),           'r-',  label=r'$G(x) = \frac{\rm{erf}(x)}{x}$')
plt.plot(x, lut_potential_expanded(x),  'r--', label='Expansion of $G(x)$')
plt.plot(x, lut_potential_log_expanded(x), 'm--', label='Expansion of $G(log(x))$')
plt.plot(x, 1.0/x,                      'r:',  label=r'$x^{-1}$')
plt.grid()
plt.legend(loc='lower right')
axprops['sharex'] = ax1
ax1.set_xlim((x.min(), x.max()))
ax1.set_ylim((0.0, 1.2))
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = fig.add_axes([axes_xmin, axes_ymin[1], axes_sizex, axes_sizey], **axprops)
plt.semilogy(x, abs(lut_field(x)     - lut_field_expanded(x)),      'b-', label=r'$|F(x) - \rm{expansion}|$')
plt.semilogy(x, abs(lut_potential(x) - lut_potential_expanded(x)),  'r-',  label=r'$|G(x) - \rm{expansion}|$')
plt.semilogy(x, abs(lut_potential(x) - lut_potential_log_expanded(x)), 'm-',  label=r'$|G(x) - \rm{expansion(log)}|$')
plt.grid()
plt.legend(loc='center left')
plt.xlabel("x")


plt.show()
