#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys, os, glob
import numpy
import matplotlib.pyplot as plt

import on_key

globber = glob.glob(os.path.join("output", "*"))

# Electron's curve?
p0_cs = -1
p0_name = "electron"
# Ion's curve?
#p0_cs = 1
#p0_name = "ion (" + str(p0_cs) + "+)"


fig = on_key.figure()
axprops = dict()
ax1 = fig.add_subplot(211, **axprops)
axprops['sharex'] = ax1
plt.setp(ax1.get_xticklabels(), visible=False)
ax2 = fig.add_subplot(212, **axprops)
plt.subplots_adjust(hspace=0.0)

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
        p1_cs = int(pot_files[cs].replace(folder,"").replace("/poten_", "").replace(".csv", ""))
        r     = data[:,0]
        pot   = data[:,1]
        Coulomb_U = p0_cs*p1_cs/r
        Umax = +10.0
        Umin = -10.0
        Coulomb_U[numpy.where(Coulomb_U > Umax)] = Umax
        Coulomb_U[numpy.where(Coulomb_U < Umin)] = Umin

        ax1.plot(r, pot, symbols[fi]+colors[cs%len(colors)], label = str(p1_cs) + "+ " + potential_shape, lw=line_width)
        ax1.plot(r, Coulomb_U, ':'+colors[cs%len(colors)], lw=line_width)

    for cs in xrange(nb_cs):
        data = numpy.loadtxt(field_files[cs], delimiter=',', skiprows=0, dtype=float)
        p1_cs = int(field_files[cs].replace(folder,"").replace("/field_", "").replace(".csv", ""))
        r     = data[:,0]
        field = data[:,1]
        Coulomb_E = p1_cs/(r*r)
        Emax = +10.0
        Emin = -10.0
        Coulomb_E[numpy.where(Coulomb_E > Emax)] = Emax
        Coulomb_E[numpy.where(Coulomb_E < Emin)] = Emin

        ax2.plot(r, field, symbols[fi]+colors[cs%len(colors)], label = str(p1_cs) + "+ " + potential_shape, lw=line_width)
        ax2.plot(r, Coulomb_E, ':'+colors[cs%len(colors)], lw=line_width)

    fi += 1


ax1.set_ylabel("Potential Energy (Hartree)")
ax2.set_ylabel("Field (au)")
ax2.set_xlabel("r (Bohr)")

#ax1.set_yscale('log')
#ax2.set_yscale('log')
#ax1.set_xscale('log')
#ax2.set_xscale('log')
ax1.grid(True)
ax2.grid(True)


#ax1.set_xlim((0.0, 4.0))
#ax1.set_ylim((-7.0, 0.5))
#ax2.set_ylim((-2.0, 10.0))

plt.suptitle("What an " + p0_name + " feels")
plt.legend(loc="best")

plt.show()

