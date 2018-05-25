#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np


def read_spectrum(fn):
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

    values = np.array(values)
    sticks = np.array(sticks)

    if 'iso' in fn:
        values[:, 1] = values[:, 1] / 3.
    return values, sticks

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--files', nargs='+')
    parser.add_argument('-s', '--sticks', action='store_true')

    args = parser.parse_args()

    fig, ax = plt.subplots()

    for fn in args.files:
        values, sticks = read_spectrum(fn)

        x, y = np.split(values, 2, 1)
        line, = ax.plot(x, y, label=fn)

        if args.sticks:
            color = line.get_color()
            x, y = np.split(sticks, 2, 1)
            ax.vlines(x, np.zeros_like(y), y, color=color)

    ax.grid()
    ax.legend()
    plt.show()


if __name__ == '__main__':
    main()

