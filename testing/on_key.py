#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pylab, sys

###################################################################
def on_key(event):
    #print 'you pressed', event.key, event.xdata, event.ydata
    if   (event.key == 'q'):
        sys.exit(0)
    #
    elif (event.key == 'w'):
        pylab.close(event.canvas.figure)
    #
    elif (event.key == 'd'):
        print "##############################################################"
        print "Please click two points to get the distance and slope."
        from matplotlib.widgets import Cursor
        cursor = Cursor(event.inaxes, useblit=False, color='red', linewidth=1)
        points = pylab.ginput(n=2, show_clicks=True, timeout=0)
        xs  = pylab.array([points[0][0], points[1][0]])
        ys  = pylab.array([points[0][1], points[1][1]])
        dx = xs[1] - xs[0]
        dy = ys[1] - ys[0]
        dy_log = pylab.log10(ys[0]) - pylab.log10(ys[1])
        dy_ln = pylab.log(ys[0]) - pylab.log(ys[1])
        angle = pylab.arctan2(dy, dx)
        print "points = ", points
        print "distance (x) =", dx
        print "distance (y) =", dy
        print "distance =", pylab.sqrt( dx**2 + dy**2 )
        print "Ratio: x0/x1 =", xs[0] / xs[1], "   x1/x0 =", xs[1] / xs[0], "  y0/y1 =", ys[0] / ys[1], "  y1/y0 =", ys[1] / ys[0]
        print "dy/y0 = ", dy/ys[0]
        print "dx/x0 = ", dx/xs[0]
        print "Angle: theta = atan2(y/x) =", angle, "rad =", angle*180.0/pylab.pi, "deg"
        print "Slope: ", dy/dx, "  1/slope =", dx/dy
        print "Slope (log10 scale):", dy_log/dx
        print "Slope (ln scale):",    dy_ln/dx
        event.inaxes.plot([points[0][0], points[1][0]], [points[0][1], points[1][1]], '--r', lw=1.0)
        event.inaxes.plot([points[0][0], points[1][0]], [points[0][1], points[1][1]],  '+r', lw=1.0)
        pylab.draw()
    #
    elif (event.key == 'a'):
        print "##############################################################"
        print "Please click a point to get the position."
        from matplotlib.widgets import Cursor
        cursor = Cursor(event.inaxes, useblit=False, color='red', linewidth=1)
        point = pylab.ginput(n=1, show_clicks=False, timeout=0)
        print "point = ", point
        event.inaxes.plot(point[0][0], point[0][1],  '+r', lw=1.0)
        pylab.draw()

###################################################################

def figure():
    fig = pylab.figure()
    fig.canvas.mpl_connect('key_press_event', on_key)
    return fig

