#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy.special import erf

import on_key
import herman_skillman

r = np.linspace(1.0e-1, 20.0, 100) # Bohr

s = 1.75

Ve = -erf(r/(m.sqrt(2) * s))/r
VHS = herman_skillman.HS(r, 1)

fig = on_key.figure()
ax1 = plt.subplot(111)

plt.plot(r, VHS,        label='$V_{HS}$')
for i in xrange(10):
    plt.plot(r, float(i)/r,  '--k')
    #plt.plot(r, float(i)/r,  '--k', label='$V_{+'+str(i)+'}$')

plt.plot(r, Ve,         label='$V_{e}$')
#plt.plot(r, V2+Ve,      label='$V_{+2} + V_e$')
plt.plot(r, 10.0/r+9.0*Ve,      label='$V_{+10} + 9*V_e$')

# Plot HS's cuttofs
cs = 1
for hsi in xrange(herman_skillman.cutoffs.shape[1]):
    plt.plot([herman_skillman.cutoffs[cs+1,hsi], herman_skillman.cutoffs[cs+1,hsi], herman_skillman.cutoffs[cs+1,hsi]], [-1.0, 1.0e-5, 20.0], '--k')

ax1.set_yscale('log')
ax1.set_ylim((1.0e-2, 10.0))
#ax1.set_ylim((-2.0, 10.0))
plt.legend()
plt.show()


