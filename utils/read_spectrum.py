#!/usr/bin/env python3

import numpy as np


def read_ttmult_spectrum(fn):
    with open(fn) as fp:
        x = list()
        y = list()
        for line in fp:
            if 'Sticks' in line:
                break
            tokens = line.split()
            if tokens:
                x.append(float(tokens[0]))
                y.append(float(tokens[1]))

        x = np.array(x)
        y = np.array(y)

        next(fp)

        xsticks = list()
        ysticks = list()
        for line in fp:
            tokens = line.split()
            if tokens:
                xsticks.append(float(tokens[0]))
                ysticks.append(float(tokens[1]))

        xsticks = np.array(xsticks)
        ysticks = np.array(ysticks)

    return x, y, xsticks, ysticks


def read_quanty_spectrum(fn, column=2):
    data = np.loadtxt(fn, skiprows=5)
    x = data[:, 0]
    y = data[:, column]
    return x, y
