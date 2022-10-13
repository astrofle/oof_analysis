"""
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits


def find_vane_peak_scans(sl):
    """
    """

    scans = np.unique(sl.scans)

    vane_scans = []
    sky_scans = []
    peak_scans = []
    focus_scans = []
    peak_rx = []
    peak_src = []
    focus_rx = []
    focus_src = []

    for scan in scans:
        
        try:
            go_file = sl.get_scan_files(scan, "GO")
        except TypeError:
            continue
        
        # Check if file exists before opening it.
        # Sometimes a file is present in the ScanLog,
        # but not in the archive.
        pfile = Path(go_file) 
        if pfile.is_file():
            hdu = fits.open(go_file)
        else:
            continue

        head = hdu[0].header

        procname = head["PROCNAME"]
        proctype = head["PROCTYPE"]
        procscan = head["PROCSCAN"]
        
        if procname == "VANECAL":
            if procscan == "SKY":
                sky_scans.append(scan)
            elif procscan == "CALIBRATION":
                vane_scans.append(scan)
        
        elif procname == "Peak":
            peak_scans.append(scan)
            peak_rx.append(head["RECEIVER"])
            peak_src.append(head["OBJECT"])

        elif procname == "FocusSubreflector":
            focus_scans.append(scan)
            focus_rx.append(head["RECEIVER"])
            focus_src.append(head["OBJECT"])

    return {"vane": vane_scans, "sky": sky_scans, 
            "peak": peak_scans, "focus": focus_scans,
            "peak_rx": peak_rx, "peak_src": peak_src,
            "focus_rx": focus_rx, "focus_src": focus_src}


class AutoPeak():
    """
    """

    def __init__(self, vane_scan, sky_scan, peak_scans, sl):
        """
        """
        self.vane_scan = vane_scan
        self.sky_scan = sky_scan
        self.peak_scans = peak_scans
        self.sl = sl


    def compute_gain(self, rx):
        """
        """
        
        rx_file_vane = self.sl.get_scan_files(self.vane_scan, rx)
        rx_file_sky = self.sl.get_scan_files(self.sky_scan, rx)

        # Get calibrator load temperature, in C.
        hdu = fits.open(rx_file_vane)
        twarm_vane = hdu[0].header["TWARM"] + 273.15
        hdu = fits.open(rx_file_sky)
        twarm_sky = hdu[0].header["TWARM"] + 273.15

        # Get power from DCR.
        # Data for two beams with Argus.
        dcr_file_vane = self.sl.get_scan_files(self.vane_scan, "DCR")
        hdu = fits.open(dcr_file_vane)
        dcr_pow_vane = np.nanmedian(hdu[3].data["DATA"], axis=(0,2))
        dcr_file_sky = self.sl.get_scan_files(self.sky_scan, "DCR")
        hdu = fits.open(dcr_file_sky)
        dcr_pow_sky = np.nanmedian(hdu[3].data["DATA"], axis=(0,2))

        # Compute the gain.
        gains = twarm_vane / (dcr_pow_vane[0] - dcr_pow_sky[0]), twarm_vane / (dcr_pow_vane[1] - dcr_pow_sky[1])
        tsys = gains * dcr_pow_sky

        self.gain = gains
        self.tsys = tsys

    
    def calibrate(self, peak_scan, rx):
        """
        """

        # Check that the Rx is the correct one.
        go_file_peak = self.sl.get_scan_files(peak_scan, "GO")
        hdu = fits.open(go_file_peak)
        peak_rx = hdu[0].header["RECEIVER"]
        if peak_rx != rx:
            return np.array([[np.nan,np.nan]]) 

        # Get the power from the first peak scan.
        dcr_file_peak = self.sl.get_scan_files(peak_scan, "DCR")
        hdu = fits.open(dcr_file_peak)
        dcr_pow_peak = hdu[3].data["DATA"]

        # Calibrate to antenna temperature.
        ta_sig = self.gain[0]*(dcr_pow_peak[:,0] - np.nanmedian(dcr_pow_peak[:,0]))
        ta_ref = self.gain[1]*(dcr_pow_peak[:,1] - np.nanmedian(dcr_pow_peak[:,1]))

        return ta_sig - ta_ref

    
    def calibrate_all(self, rx):
        """
        """
        
        self.ta_peak = []

        for i,ps in enumerate(self.peak_scans):
            self.ta_peak.append(self.calibrate(ps, rx))
       
    
    def get_peak_ta(self, rx):
        """
        """

        self.compute_gain(rx)
        self.calibrate_all(rx)

        self.peak_ta = {}

        for i,ps in enumerate(self.peak_scans):
            if len(self.ta_peak[i][:,0]) == 0 \
               or np.all(np.isnan(self.ta_peak[i])):
                self.peak_ta[ps] = np.nan
            else:
                self.peak_ta[ps] = np.nanmax(self.ta_peak[i][:,0]) - \
                                   np.nanmedian(self.ta_peak[i][:,0])

        return self.peak_ta


    def plot_peak_ta(self):
        """
        """
    
        fig = plt.figure(dpi=150, frameon=False)

        nx = 2
        ny = 2

        for i,ps in enumerate(self.peak_scans[0:nx*ny]):
            ax = fig.add_subplot(2,2,i+1)
            ax.plot(self.ta_peak[i][:,0])
            ax.text(0.1, 0.85, 
                    f"{self.peak_ta[ps]:.2f}", 
                    transform=ax.transAxes)
            ax.set_xlabel("Sample")
            #ax.set_ylabel("")
            

        if len(self.peak_scans) > nx*ny:
            fig = plt.figure(dpi=150, frameon=False)
            ax = fig.add_subplot(111)
            ax.plot(self.ta_peak[nx*ny][:,0])
            ax.text(0.1, 0.85,
                    f"{self.peak_ta[self.peak_scans[-1]]:.2f}",
                    transform=ax.transAxes)
