"""
Given the output of getgbt4mmcal.py this will 
generate the zenith opacity at the desired frequency.
"""

import argparse
import numpy as np

import weather_forecast as wf


def parse_args():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str,
                        help='Output of getgbt4mmcal.py'
                             'Example: AGMVA21B_164_02_gbt4mmcal.dat')
    parser.add_argument('-f', '--freq', type=float,
                        help='Frequency in GHz')
    parser.add_argument('-o', '--output', type=str,
                        help='Output file name.')
    parser.add_argument('-c', '--column', type=int,
                        default=6,
                        help='Column with the MJDs in the input file.'
                             'Default: 6')

    args = parser.parse_args()

    return args


def get_tau0(mjds, frequency):
    """
    """

    tau0 = np.zeros(len(mjds), dtype=float)

    for i,mjd in enumerate(mjds):

        tau0[i] = wf.get_forecast_tau0(mjd, frequency)

    return tau0


if __name__ == '__main__':

    args = parse_args()

    mjds = np.loadtxt(args.input_file, usecols=(args.column))

    tau0 = get_tau0(mjds, args.freq)

    np.savetxt(args.output, np.c_[mjds, tau0], header='MJD tau0')
