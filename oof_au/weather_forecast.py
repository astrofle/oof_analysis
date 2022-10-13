"""
Wrapper around getForecastValues.tcl
"""

import subprocess

tcl_path = '/home/cleoversions/Cleo/mainscreens/getForecastValues.tcl'


def get_forecast_tau0(mjd, frequency):
    """
    Returns the zenith opacity for a given MJD and frequency.
    """

    args = [tcl_path, '-freqList', f'{frequency}',
            '-timeList', f'{mjd}', '-type', 'Opacity']

    result = subprocess.run(args, stdout=subprocess.PIPE)

    tau0 = float(result.stdout.decode('ascii').strip().split('=')[1])

    return tau0


def get_forecast_tamt(mjd, frequency):
    """
    Returns the atmospheric temperature for a given MJD and frequency.
    """

    args = [tcl_path, '-freqList', f'{frequency}',
            '-timeList', f'{mjd}', '-type', 'Tatm']

    result = subprocess.run(args, stdout=subprocess.PIPE)

    tatm = float(result.stdout.decode('ascii').strip().split('=')[1])

    return tatm


def get_forecast_atmtsys(mjd, frequency, elevation):
    """
    Returns the atmospheres contribution to the system temperature
    for a given MJD, frequency and elevation.
    """

    args = [tcl_path, '-freqList', f'{frequency}',
            '-timeList', f'{mjd}', '-elev', f'{elevation}', 
            '-type', 'AtmTsys']

    result = subprocess.run(args, stdout=subprocess.PIPE)

    atmtsys = float(result.stdout.decode('ascii').strip().split('=')[1])

    return atmtsys
