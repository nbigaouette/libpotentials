#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys, os, glob
import numpy
import matplotlib.pyplot as plt

import on_key

globber = glob.glob(os.path.join("output", "*"))

fig = on_key.figure()
axprops = dict()
#ax1 = plt.subplot(211)
#ax2 = plt.subplot(212)
ax1 = fig.add_axes([0.15, 0.51, 0.6, 0.4], **axprops)
axprops['sharex'] = ax1
plt.setp(ax1.get_xticklabels(), visible=False)
ax2 = fig.add_axes([0.15, 0.1, 0.6, 0.4], **axprops)

for folder in globber:

    # Skip files
    if (os.path.isfile(folder)):
        continue

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
        ax1.plot(r, pot, label = str(charge_state) + " " + folder.replace("output/", ""))

    for cs in xrange(nb_cs):
        data = numpy.loadtxt(field_files[cs], delimiter=',', skiprows=0, dtype=float)
        r     = data[:,0]
        field = data[:,1]

        charge_state = int(field_files[cs].replace(folder,"").replace("/field_", "").replace(".csv", ""))
        ax2.plot(r, field, label = str(charge_state) + " " + folder.replace("output/", ""))


ax1.set_ylabel("Potential (au)")
ax2.set_ylabel("Field (au)")
ax2.set_xlabel("r (Bohr)")

ax1.set_ylim((-2.0, 11.0))
#ax2.set_ylim((-1.0, 5.0))

plt.legend()

plt.show()

