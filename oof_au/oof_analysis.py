

import numpy as np

from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy import constants as ac


import alma_flux
import get_tau0
import scan_log
import ruze


time_fmt = "%d-%b-%Y"


class OOFAnalysis():
    """
    """

    def __init__(self, path, project, session):

        self.path = path
        self.project = project
        self.session = session
        self.projid = f"{project}_{session}"


    def compute_surface_rms(self, scan, ta, verbose=False, source=False):
        """
        """

        #if np.isnan(ta):
        #    return {"rms": np.nan, "tau0": np.nan}

        # Open the scan log and retrieve the date and paths to the GO and Antenna fits files.
        sl = scan_log.ScanLog(self.path, self.projid)
        scan_date = sl.get_scan_dateobs(scan)
        go_file = sl.get_scan_files(scan, "GO")
        an_file = sl.get_scan_files(scan, "Antenna")

        # Get the antenna elevation.
        an_hdu = fits.open(an_file)
        el = an_hdu[0].header["SOBSC_EL"]

        # Get the observing frequency, source name and receiver.
        go_hdu = fits.open(go_file)
        freq = go_hdu[0].header["SKYFREQ"]
        if not source:
            source = go_hdu[0].header["OBJECT"]
        rx = go_hdu[0].header["RECEIVER"]

        if verbose:
            print(f"Observed {source} at an elevation of {el} and frequency {freq:g} Hz with {rx}")

        # Query the ALMA catalogue for the flux density of the source at the time of the peak scan.
        tab = alma_flux.get_alma_flux(Time(scan_date).strftime(time_fmt), freq=freq, source=f"J{source}")
        if verbose:
            print(f"ALMA flux for the source:")
            print(tab['FluxDensity'])

        # Correct the antenna temperature for the atmospheric opacity.
        airm = 1/np.sin(np.deg2rad(el))
        mjd = Time(scan_date).mjd
        tau0 = get_tau0.get_tau0([mjd], freq*1e-9)
        # Only do this if the Rx is not Argus.
        if rx != "RcvrArray75_115":
            ta_p = ta*np.exp(tau0*airm)
        else:
            ta_p = ta

        if verbose:
            print(f"Using a zenith opacity of {tau0} to conver from Ta to Ta'={ta_p} K")

        # Compute the aperture efficiency.
        eta_a = ruze.aperture_efficiency(ta_p, tab["FluxDensity"])

        # Turn the aperture efficiency into a surface rms.
        lmbd = (ac.c/(freq*u.Hz)).to("m")
        surf_rms = ruze.eta_to_surface_rms(eta_a[0], lmbd, eta_0=0.71)

        try:
            surf_rms = surf_rms.to("um")
        except AttributeError:
            pass

        output = {"rms": surf_rms, "tau0": tau0, "mjd": mjd, "El": el, "object": source}

        return output

