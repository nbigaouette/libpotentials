#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.optimize import leastsq

import on_key
import herman_skillman

def Velectron(r, sigma):
    return -erf(r / (m.sqrt(2) * sigma)) / r

r = np.linspace(1.0, 50.0, 1000) # Bohr
#r = np.linspace(herman_skillman.cutoffs[2,0], herman_skillman.cutoffs[2,1], 100) # Bohr
#r = np.linspace(herman_skillman.cutoffs[2,1], herman_skillman.cutoffs[2,2], 100) # Bohr
#r = np.linspace(herman_skillman.cutoffs[2,2], herman_skillman.cutoffs[2,3], 100) # Bohr

#r = np.linspace(herman_skillman.cutoffs[2,1], herman_skillman.cutoffs[2,3], 100) # Bohr
#r = np.linspace(0.5, herman_skillman.cutoffs[2,3], 100) # Bohr


VHS = herman_skillman.HS(r, 1)

fig = on_key.figure()
ax1 = plt.subplot(111)

# Plot HS's cuttofs
#cs = 1
#for hsi in xrange(herman_skillman.cutoffs.shape[1]):
    #plt.plot([herman_skillman.cutoffs[cs+1,hsi], herman_skillman.cutoffs[cs+1,hsi], herman_skillman.cutoffs[cs+1,hsi]], [-1.0, 1.0e-5, 20.0], '--k')

def Vtot(sigmas, r):
    Vtot = float(Ne+1)/r
    for e in xrange(Ne):
        Vtot += Velectron(r, sigmas[e])
    return Vtot
def Vtot_err(p, r, x):
    return (Vtot(p, r) - x)

#for Ne in xrange(1, 40, 5):
#for Ne in [14, 16, 18]:
for Ne in [14]:
    r = np.linspace(1.0, 50.0, 1000) # Bohr

    sigmas0 = [1.75]*Ne
    sigmas1, success = leastsq(Vtot_err, sigmas0[:], args=(r, VHS))
    sigmas1.sort()

    r = np.linspace(1.0e-5, 50.0, 1000) # Bohr
    plt.plot(r, Vtot(sigmas1, r), label='Fit ($N_e =$' + str(Ne) + ')')
    print "sigmas1 =", sigmas1


plt.plot(r, herman_skillman.HS(r, 1), '--k', label='$V_{HS}$')


for e in xrange(len(sigmas1)):
    plt.plot(r, Velectron(r,sigmas1[e]), label='$V_e(\sigma = '+str(sigmas1[e])+')$')

#ax1.set_yscale('log')
#ax1.set_ylim((1.0e-2, 10.0))
ax1.set_ylim((-2.0, 10.0))
plt.legend()
plt.show()


