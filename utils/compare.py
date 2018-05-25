#!/usr/bin/env python3
# conding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import sys


def read_ttmult_spectrum(fn):
    with open(fn) as fp:
        values = list()
        for line in fp:
            if 'Sticks' in line:
                break
            tokens = line.split()
            if tokens:
                values.append([float(token) for token in tokens])

        next(fp)

        sticks = list()
        for line in fp:
            tokens = line.split()
            sticks.append([float(token) for token in tokens])

    return np.array(values), np.array(sticks)


def read_quanty_spectrum(fn, column=2):
    data = np.loadtxt(fn, skiprows=5)
    x = data[:, 0]
    y = -data[:, column]
    return np.array([x, y]).transpose()


def main():
    fig, ax = plt.subplots()

    if len(sys.argv) > 1:
        experiment = sys.argv[1]

    npoints = 8192

    data, _ = read_ttmult_spectrum('TTMult_without_Bander/input_iso.xy')
    x = data[:, 0]
    y = data[:, 1] / 3
    xn = np.linspace(min(x), max(x), npoints)
    y1 = np.interp(xn, x, y)
    ax.plot(xn, y1, label='TTMult_without_Bander')

    data, _ = read_ttmult_spectrum('TTMult/input_iso.xy')
    x = data[:, 0]
    y = data[:, 1] / 3
    y2 = np.interp(xn, x, y)
    ax.plot(xn, y2, label='TTMult')

    print('Largest difference between Racah and Bander: {0:.2E}'.format(np.max(np.abs(y1 - y2))))

    data = read_quanty_spectrum('Quanty/input.spec')
    x = data[:, 0]
    y = data[:, 1]
    y3 = np.interp(xn, x, y)
    ax.plot(xn, y3, label='Quanty')

    ax.legend()
    plt.show()


if __name__ == '__main__':
    main()
