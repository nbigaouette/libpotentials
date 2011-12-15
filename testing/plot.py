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

def cubic_spline(cs, rV, V, rE, E):
    # http://people.math.sfu.ca/~stockie/teaching/macm316/notes/splines.pdf

    # Find where field is maximum
    index = numpy.where(E == E.max())[0][0]
    Emax = E[index]
    rEmax = rE[index]
    #print "index =", index
    #print "Emax =", Emax
    #print "rEmax =", rEmax
    f = 6.67
    #index = numpy.where((Emax/3.9 <= E) & (E <= Emax/3.8))[0][0]
    index = numpy.where((0.99*Emax/f <= E) & (E <= 1.01*Emax/f))[0][0]
    p1 = (rE[index], E[index])
    p2 = (rE[index+10], E[index+10])
    #print "index =", index
    #print "rE[index] =", rE[index]
    #print "E[index] =", E[index]
    #sys.exit(0)

    #V0 = -1.5
    V0 = V[0]
    ## Cubic spline with 3 points
    n = 2
    p0 = (0.0, 0.0)
    # This gives the right potential
    #p1 = (2.1876497612049692, 0.60793404213156088)
    #p2 = (2.1981567182814885, 0.59642039133345082)
    r = numpy.array([p0[0], p1[0], p2[0]])
    y = numpy.array([p0[1], p1[1], p2[1]])

    A = numpy.zeros((3,3))
    b = numpy.zeros((3))
    h = numpy.array([r[1]-r[0], r[2]-r[1]])
    b[0] = 0.0
    b[1] = 6.0 * ( (y[2] - y[1])/h[1] - (y[1] - y[0])/h[0] )
    b[2] = 0.0
    ## First line
    A[0,0] = 1.0
    ## Second line
    A[1,0] = h[0]
    A[1,1] = 2.0 * (h[0] + h[1])
    A[1,2] = h[1]
    ## Third line
    A[2,2] = 1.0


    ## Cubic spline with 4 points
    #n = 3
    #p0 = (0.0, 0.0)
    #p1 = (0.90464476213767697, 1.2121277423263328)
    #p2 = (1.4765745007680491, 2.376101321585903)
    #p3 = (1.5341781874039935, 2.0952643171806171)
    #r = numpy.array([p0[0], p1[0], p2[0], p3[0]])
    #y = numpy.array([p0[1], p1[1], p2[1], p3[1]])
    ##r = numpy.array([0.0, 0.75049374588545104, 1.5, 1.5487339576829691])
    ##y = numpy.array([-1.5, -1.5, -1.4468473596608349, -1.3381362921849833])
    #A = numpy.zeros((n+1,n+1))
    #b = numpy.zeros((n+1))
    #h = numpy.array([r[1]-r[0], r[2]-r[1], r[3]-r[2]])
    #b[0] = 0.0
    #b[1] = 6.0 * ( (y[2] - y[1])/h[1] - (y[1] - y[0])/h[0] )
    #b[2] = 6.0 * ( (y[3] - y[2])/h[2] - (y[2] - y[1])/h[1] )
    #b[3] = 0.0
    ## First line
    #A[0,0] = 1.0
    ## Second line
    #A[1,0] = h[0]                   # First column
    #A[1,1] = 2.0 * (h[0] + h[1])    # Second column
    #A[1,2] = h[1]                   # Third column
    ## Third line
    #A[2,1] = h[1]                   # Second column
    #A[2,2] = 2.0 * (h[1] + h[2])    # Third column
    #A[2,2] = h[2]                   # Fourth column
    ## Fourth line
    #A[3,3] = 1.0                    # Fifth column

    m = numpy.linalg.solve(A, b)
    #print "A =", A
    #print "b =", b
    #print "h =", h
    #print "m =", m
    #print "A.dot(m) - b =", A.dot(m) - b
    # Spline coefficients:
    ai = numpy.zeros((n))
    bi = numpy.zeros((n))
    ci = numpy.zeros((n))
    di = numpy.zeros((n))
    for i in xrange(n):
        ai[i] = y[i]
        bi[i] = (y[i+1] - y[i]) / h[i] - h[i]/2.0*m[i] - h[i]/6.0*(m[i+1]-m[i])
        ci[i] = m[i] / 2.0
        di[i] = (m[i+1] - m[i]) / (6.0 * h[i])
    new_r = numpy.linspace(r[0], r[-1], 100)
    new_V = numpy.zeros((len(new_r)))
    #new_dV = numpy.zeros((len(new_r)))
    new_IntV = numpy.zeros((len(new_r)))
    last_IntV = V0
    #last_IntV = 0.0
    for i in xrange(n):
        indices = numpy.nonzero((r[i] <= new_r) & (new_r <= r[i+1]))[0]
        new_V[indices]  = ai[i] + bi[i]*(new_r[indices] - r[i]) +     ci[i]*(new_r[indices]-r[i])**2  +     di[i]*(new_r[indices]-r[i])**3
        #new_dV[indices] =         bi[i]                         + 2.0*ci[i]*(new_r[indices]-r[i])     + 3.0*di[i]*(new_r[indices]-r[i])**2
        #new_IntV[indices] = V0 + ai[i]*(new_r[indices]-r[0]) + bi[i]/3.0*(new_r[indices]-r[0])**2 + ci[i]/4.0*(new_r[indices]-r[0])**3 + di[i]/5.0*(new_r[indices]-r[0])**4
        new_IntV[indices] = ai[i]*(new_r[indices]-r[i]) + bi[i]/2.0*(new_r[indices]-r[i])**2 + ci[i]/3.0*(new_r[indices]-r[i])**3 + di[i]/4.0*(new_r[indices]-r[i])**4 + last_IntV
        last_IntV = new_IntV[indices][-1]

    return new_r, new_V, new_IntV, r, y

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
        data        = numpy.loadtxt(pot_files[cs], delimiter=',', skiprows=0, dtype=float)
        p1_cs       = int(pot_files[cs].replace(folder,"").replace("/poten_", "").replace(".csv", ""))
        rU          = data[:,0]
        U           = data[:,1]
        Coulomb_U   = p0_cs*p1_cs/rU
        Umax = +2.0
        Umin = -2.0
        Coulomb_U[numpy.where(Coulomb_U > Umax)] = Umax
        Coulomb_U[numpy.where(Coulomb_U < Umin)] = Umin

        data        = numpy.loadtxt(field_files[cs], delimiter=',', skiprows=0, dtype=float)
        rE          = data[:,0]
        E           = data[:,1]
        Coulomb_E   = p1_cs/(rE*rE)
        Emax = +2.0
        Emin = -2.0
        Coulomb_E[numpy.where(Coulomb_E > Emax)] = Emax
        Coulomb_E[numpy.where(Coulomb_E < Emin)] = Emin


        try:
            if (potential_shape == "HermanSkillman"):
                new_r, new_E, new_IntE, pts_r, pts_E = cubic_spline(p1_cs, rU, U, rE, E)
                ax1.plot(new_r, new_IntE, colors[cs%len(colors)])
                ax2.plot(new_r, new_E, colors[cs%len(colors)])
                ax2.plot(pts_r, pts_E, colors[cs%len(colors)] + 'x', ms=10, markeredgewidth=3)
        except:
            pass

        ax1.plot(rU, U, symbols[fi]+colors[cs%len(colors)], label = str(p1_cs) + "+ " + potential_shape, lw=line_width)
        ax1.plot(rU, Coulomb_U, ':'+colors[cs%len(colors)], lw=line_width)
        ax2.plot(rE, E, symbols[fi]+colors[cs%len(colors)], label = str(p1_cs) + "+ " + potential_shape, lw=line_width)
        ax2.plot(rE, Coulomb_E, ':'+colors[cs%len(colors)], lw=line_width)

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

