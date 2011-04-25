#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.optimize import leastsq

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
z = 8.0
plt.plot(r, z/r+(z-1)*Ve,      label='$V_{+'+str('%d' % z) + '} + '+str('%d' % (z-1))+'*V_e$')

# Plot HS's cuttofs
cs = 1
for hsi in xrange(herman_skillman.cutoffs.shape[1]):
    plt.plot([herman_skillman.cutoffs[cs+1,hsi], herman_skillman.cutoffs[cs+1,hsi], herman_skillman.cutoffs[cs+1,hsi]], [-1.0, 1.0e-5, 20.0], '--k')

# Least-square fitting
# http://www.scipy.org/Cookbook/FittingData
#fitfunc = lambda p, r: -erf(r / (m.sqrt(2)*p[0])) / r # Target function
#errfunc = lambda p, r, x: fitfunc(p, r) - x # Distance to the target function
#p0 = [-15., 0.8, 0., -1.] # Initial guess for the parameters
#p1, success = optimize.leastsq(errfunc, p0[:], args=(Tx, tX))


ax1.set_yscale('log')
ax1.set_ylim((1.0e-2, 10.0))
#ax1.set_ylim((-2.0, 10.0))
plt.legend()
plt.show()


