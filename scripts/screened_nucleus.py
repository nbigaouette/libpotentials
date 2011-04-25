#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy.special import erf

import on_key


r = np.linspace(1.0e-1, 10.0, 1000) # Bohr

#Z = 54
#sigmas = [1.0, 0.5, 0.25]
sigmas = [1.0, 0.5, 0.25]

V1 = 1.0 / r
Vn = len(sigmas) * V1
Ve = [
        -erf(r/(m.sqrt(2) * sigmas[0]))/r,
        -erf(r/(m.sqrt(2) * sigmas[1]))/r,
        -erf(r/(m.sqrt(2) * sigmas[2]))/r
    ]


fig = on_key.figure()

plt.plot(r, V1,     '--', label='$V_{+1}$')
plt.plot(r, Vn,     label='$V_{+' + str(len(sigmas)) + '}$')
plt.plot(r, Ve[0],  label='$V_{e,0}$')
plt.plot(r, Ve[1],  label='$V_{e,1}$')
plt.plot(r, Ve[2],  label='$V_{e,2}$')
plt.plot(r, Vn+Ve[0]+Ve[1]+Ve[2],  label='$V_{+' + str(len(sigmas)) + '} + \sum V_e$')

plt.legend()
plt.show()
