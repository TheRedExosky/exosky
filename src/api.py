import ssl
ssl._create_default_https_context = ssl._create_unverified_context

from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from dataclasses import dataclass
import numpy as np

from calc import spherical_to_cartesian

@dataclass
class StarObject():
    """ Star object to draw in matplotlib """
    x: float
    y: float
    z: float
    luminosity: float
    radius: float
    temperature: float


def fetch_api(ra=280, dec=-60):
    """
    Fetch AstroQuery to fetch star data given a `ra` and `dec` position.
    """
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    width = u.Quantity(0.1, u.deg)
    height = u.Quantity(0.1, u.deg)

    # fetch star table
    stars = Gaia.query_object_async(coordinate=coord, width=width, height=height)
    stars.pprint(max_lines=12, max_width=130)

    # collect relevant data for futher plotting from star table
    objs = []
    for idx in range(len(stars)):
        row = stars[idx]

        ra = row['ra']
        dec = row['dec']
        parallax = row['parallax']
        luminosity = row['phot_g_mean_mag']
        bp_rp = row['bp_rp']
        temperature = row['teff_gspphot']

        if not temperature:
             temperature = 5601 * (0.4 * bp_rp + 1) ** (-1.6)

        distance = 1000 / parallax
        stefan_boltzmann_constant = 5.67e-8  # W/m^2/K^4
        luminosity_in_watt = luminosity * 3.828e26
        radius = np.sqrt(luminosity_in_watt / (4 * np.pi * stefan_boltzmann_constant * temperature**4))

        x, y, z = spherical_to_cartesian(ra, dec, distance)
        drawobject = StarObject(x, y, z, luminosity, radius, temperature)

        objs.append(drawobject)
    return objs
