#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys, os, glob
import numpy
import matplotlib.pyplot as plt

import on_key

globber = glob.glob(os.path.join("output", "*"))


# Taken from libpotentials.git's Potential.cpp, variables fit_lt_R1, fit_lt_R2 and fit_lt_R3 (247a4be84f2b7e84504a723eeeb07b8c6c2c3537)
HS_cutoffs = numpy.array([
        [0.000, 0.00, 0.00, 0.000], # Electron
        [0.073, 1.00, 3.00,  6.00], # Neutral
        [0.073, 1.00, 5.00, 12.00], # 1+
        [0.073, 0.60, 2.00,  6.00], # 2+...
        [0.073, 0.75, 4.60, 12.00],
        [0.073, 0.75, 1.49,  2.00],
        [0.073, 0.75, 1.50,  1.49],
        [0.073, 0.75, 1.10,  1.10],
        [0.073, 0.35, 0.96,  0.96],
        [0.020, 0.00, 0.00,  0.00]
    ])




fig = on_key.figure()
axprops = dict()
ax1 = fig.add_axes([0.15, 0.51, 0.6, 0.4], **axprops)
axprops['sharex'] = ax1
plt.setp(ax1.get_xticklabels(), visible=False)
ax2 = fig.add_axes([0.15, 0.1, 0.6, 0.4], **axprops)

colors  = ['b', 'r', 'm', 'c', 'g', 'y']
symbols = ['-', '--', ':', '-.']
line_width = 2

fi = 0
for folder in globber:

    # Skip files
    if (os.path.isfile(folder)):
        continue

    potential_shape = folder.replace("output/", "")

    field_files = glob.glob(os.path.join(folder, "field_*"))
    field_files.sort()

    pot_files = glob.glob(os.path.join(folder, "pot*"))
    pot_files.sort()

    assert(len(field_files) == len(pot_files))
    nb_cs = len(field_files)

    for cs in xrange(nb_cs):
        data = numpy.loadtxt(pot_files[cs], delimiter=',', skiprows=0, dtype=float)
        r     = data[:,0]
        pot   = data[:,1]

        charge_state = int(pot_files[cs].replace(folder,"").replace("/poten_", "").replace(".csv", ""))
        ax1.plot(r, pot, symbols[fi]+colors[cs], label = str(charge_state) + " " + potential_shape, lw=line_width)

        # Plot HS's cuttofs
        if (potential_shape == "HermanSkillman"):
            for hsi in xrange(1,4):
                ax2.plot([HS_cutoffs[cs,hsi], HS_cutoffs[cs,hsi], HS_cutoffs[cs,hsi]], [-1.0, 1.0e-5, 20.0], '--' + colors[cs])

    for cs in xrange(nb_cs):
        data = numpy.loadtxt(field_files[cs], delimiter=',', skiprows=0, dtype=float)
        r     = data[:,0]
        field = data[:,1]

        charge_state = int(field_files[cs].replace(folder,"").replace("/field_", "").replace(".csv", ""))
        ax2.plot(r, field, symbols[fi]+colors[cs], label = str(charge_state) + " " + potential_shape, lw=line_width)

        # Plot HS's cuttofs
        if (potential_shape == "HermanSkillman"):
            for hsi in xrange(1,4):
                ax2.plot([HS_cutoffs[cs,hsi], HS_cutoffs[cs,hsi], HS_cutoffs[cs,hsi]], [-1.0, 1.0e-5, 20.0], '--' + colors[cs])

    fi += 1


ax1.set_ylabel("Potential (au)")
ax2.set_ylabel("Field (au)")
ax2.set_xlabel("r (Bohr)")

ax1.set_ylim((-2.0, 11.0))
#ax2.set_ylim((-1.0, 5.0))
#ax2.set_yscale('log')
ax1.grid(True)
ax2.grid(True)

plt.legend()

plt.show()

