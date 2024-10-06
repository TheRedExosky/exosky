"""
Everything related to fetching star data from AstroQuery.
"""

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


def fetch_api(ra=280, dec=-60, limit=100, min_brightness=21):
    """
    Fetch a random subset of stars using the `random_index` column from Gaia.
    Args:
        ra (float): Right Ascension of the center of the search region.
        dec (float): Declination of the center of the search region.
        limit (int): Number of random stars to retrieve.
    Returns:
        List[StarObject]: A list of star objects for plotting.
    """
    # Split into 12 patches to not run into timeouts
    num_patches = 8
    radius_deg = u.Quantity(30, u.deg)
    ra_centers = np.linspace(0, 360, num_patches, endpoint=False)

    # Iterate over all starting points & query with smaller angles
    for ra_center in ra_centers:
        adql_query = f"""
        SELECT TOP {limit} *
        FROM gaiadr3.gaia_source
        WHERE CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra_center}, {dec}, {radius_deg.value})
        ) = 1
        AND parallax IS NOT NULL AND parallax > 0
        AND phot_g_mean_mag <= {min_brightness}
        ORDER BY random_index
        """

        # Führe die Abfrage durch und gib die Ergebnisse zurück
        job = Gaia.launch_job_async(adql_query)
        stars = job.get_results()

        # Sterne für die Darstellung sammeln
        objs = []
        for idx in range(len(stars)):
            row = stars[idx]

            ra = row['ra']
            dec = row['dec']
            parallax = row['parallax']

            if not parallax or parallax <= 0:
                continue

            luminosity = row['phot_g_mean_mag']
            bp_rp = row['bp_rp']
            temperature = row['teff_gspphot']

            if not temperature:
                temperature = 5601 * (0.4 * bp_rp + 1) ** (-1.6)

            # Berechnung der Distanz und des Radius
            distance = 1000 / parallax
            stefan_boltzmann_constant = 5.67e-8  # W/m^2/K^4
            luminosity_in_watt = luminosity * 3.828e26
            radius = np.sqrt(luminosity_in_watt / (4 * np.pi * stefan_boltzmann_constant * temperature**4))

            # Konvertiere in kartesische Koordinaten
            x, y, z = spherical_to_cartesian(ra, dec, distance)
            drawobject = StarObject(x, y, z, luminosity, radius, temperature)
            objs.append(drawobject)

    return objs
