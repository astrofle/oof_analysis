
import os

from astropy.time import Time
from astropy.io import fits


class ScanLog(object):
    """
    """

    def __init__(self, path, projid):

        self.filename = f"{path}/{projid}/ScanLog.fits"
        self.hdu = fits.open(self.filename)
        self.scans = self.hdu[1].data["SCAN"]
        self.path = path

    
    def scan_mask(self, scan):
        """
        """

        if self.scans.max() < scan:
            raise ValueError(f"Scan requested {scan} is larger than largest scan in ScanLog {self.scans.max()}")

        mask = ( self.scans == scan )

        return mask


    def get_end_time(self, scan):
        """
        """

        mask = self.scan_mask(scan)
        date = Time(self.hdu[1].data[mask][-1][-1].split(' ')[3], format='mjd')
        time = self.hdu[1].data[mask][-1][-1].split(' ')[-1]

        return f"{date.to_value('iso').split(' ')[0]}T{time}"


    def get_start_time(self, scan):
        """
        """

        mask = self.scan_mask(scan)
        date = Time(self.hdu[1].data[mask][-2][-1].split(' ')[3], format='mjd')
        time = self.hdu[1].data[mask][-2][-1].split(' ')[-1]

        return f"{date.to_value('iso').split(' ')[0]}T{time}"


    def get_scan_files(self, scan, manager):
        """
        manager : string
            One of the M&C managers, e.g., Antenna, GO, DCR
        returns : string
            The first occurrence of the manager
        """

        mask = self.scan_mask(scan)
        files = self.hdu[1].data["FILEPATH"][mask]
    #    try:
        mngr_file = next((s for s in files if manager in s and "2016_01_26_00:00:00" not in s), None)[1:]
    #    except TypeError:
    #        mngr_file = ""


        return f"{self.path}/{mngr_file}" 


    def get_scan_dateobs(self, scan):
        """
        """

        mask = self.scan_mask(scan)

        return self.hdu[1].data[mask]["DATE-OBS"][0]
