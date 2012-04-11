#!/usr/bin/env python2
# -*- coding: utf-8 -*-


if __name__ == "__main__":
    import sys
    print "Don't run this file directly. It is used by other scripts."
    sys.exit(0)



colors = ['b', 'r', 'g', 'm', 'c']
symbols = ['-', '--', '-.', ':']


def color(i):
    return colors[i % len(colors)%len(symbols)]

def symbol(i):
    return symbols[int(i/len(colors))%len(symbols)]

def symb_col(i):
    return color(i) + symbol(i)

