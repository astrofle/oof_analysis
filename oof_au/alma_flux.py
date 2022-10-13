
import warnings

from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning


def get_alma_flux(date, freq, source):
    """
    Queries the ALMA calibrator source catalogue: https://almascience.nrao.edu/sc/

    date : 
        Date in the format DD-MONTH-YEAR,
        e.g., "10-DEC-2021"
    freq : float
        Frequency in Hz, e.g., 90e9
    source : string
        Source name, e.g., "J0006-0623"
    
    Details of the flux service: 
    https://almascience.nrao.edu/documents-and-tools/cycle8/flux-service-of-the-alma-source-catalogue
    """
    
    url = "https://almascience.nrao.edu/sc/flux"
    query = f"{url}?DATE={date}&FREQUENCY={freq}&NAME={source}"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=AstropyWarning)
        table = Table.read(query)
    
    return table
