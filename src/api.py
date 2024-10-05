import ssl
ssl._create_default_https_context = ssl._create_unverified_context

from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from dataclasses import dataclass

from calc import spherical_to_cartesian

@dataclass
class DrawObject():
    """ Star object to draw """
    x: float
    y: float
    z: float
    brightness: float

def run_api(ra=280, dec=-60):
    """ do things! """

    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    width = u.Quantity(0.1, u.deg)
    height = u.Quantity(0.1, u.deg)

    stars = Gaia.query_object_async(coordinate=coord, width=width, height=height)

    objs = []
    for idx in range(len(stars)):
        row = stars[idx]

        ra = row['ra']
        dec = row['dec']
        parallax = row['parallax']
        # XXX: 
        brightness = 1.0

        distance = 1000 / parallax
        x, y, z = spherical_to_cartesian(ra, dec, distance)

        drawobject = DrawObject(x, y, z, brightness)
        objs.append(drawobject)
    return objs


# obs = run_api()
# print(obs)