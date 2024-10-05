import ssl
ssl._create_default_https_context = ssl._create_unverified_context

from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from dataclasses import dataclass

from calc import spherical_to_cartesian

import numpy as np

@dataclass
class DrawObject():
    """ Star object to draw """
    x: float
    y: float
    z: float
    luminosity: float
    radius: float
    temp: float

def run_api(ra=280, dec=-60):
    """ do things! """

    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    width = u.Quantity(0.1, u.deg)
    height = u.Quantity(0.1, u.deg)

    stars = Gaia.query_object_async(coordinate=coord, width=width, height=height)

    stars.pprint(max_lines=12, max_width=130)

    objs = []
    for idx in range(len(stars)):
        row = stars[idx]

        ra = row['ra']
        dec = row['dec']
        parallax = row['parallax']
        luminosity = row['phot_g_mean_mag']
        bp_rp = row['bp_rp']
        temp = row['teff_gspphot']

        if not(temp):
            temp = 5601 * (0.4 * bp_rp + 1) ** (-1.6)

        # Calculate distance
        distance = 1000 / parallax

        # Stefan-Boltzmann-Konstante
        sigma = 5.67e-8  # W/m^2/K^4

        # Calculation of luminosity in watts (solar luminosity = 3.828e26 W)
        luminosity_in_watt = luminosity * 3.828e26

        # Calculate radius
        radius = np.sqrt(luminosity_in_watt / (4 * np.pi * sigma * temp**4))

        x, y, z = spherical_to_cartesian(ra, dec, distance)

        drawobject = DrawObject(x, y, z, luminosity, radius, temp)
        objs.append(drawobject)
    return objs


# obs = run_api()
# print(obs)