
import numpy as np


def eta_to_surface_rms(eta_a, lmbd, eta_0=0.71):
    """
    Returns the surface rms given an aperture efficiency.

    eta_a : float
        Aperture efficiency.
    lmbd : float, Quantity
        Wavelength of the observation.
    eta_0 : float
        Aperture efficiency at long wavelengths.
        For the GBT 0.71
    returns : float, Quantity
        Surface rms at `lmbd` in units of `lmbd`.
    """
    
    return np.sqrt(-np.log(eta_a/eta_0))/(4*np.pi)*lmbd 


def aperture_efficiency(ta, snu):
    """
    Eq. (7) in Frayer et al. (2019)

    https://library.nrao.edu/public/memos/gbt/GBT_302.pdf
    """

    return 0.352*ta/snu
