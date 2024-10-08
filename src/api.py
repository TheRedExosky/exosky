"""
Everything related to fetching star data from AstroQuery.
"""

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import astropy.units as u
from astroquery.gaia import Gaia
from dataclasses import dataclass
import numpy as np

from calc import spherical_to_cartesian
from debug import debug

@dataclass
class StarObject:
    """ Star object to draw in matplotlib """
    x: float
    y: float
    z: float
    luminosity: float
    radius: float
    temperature: float
    designation: str


def fetch_api(ra: float, dec: float, limit: int = 150, min_brightness: int = 21):
    """
    Fetch a random subset of stars using the `random_index` column from Gaia.
    Args:
        ra (float): Right Ascension of the center of the search region.
        dec (float): Declination of the center of the search region.
        limit (int): Number of random stars to retrieve.
        min_brightness (int): Minimum brightness of the stars. Lower number means brighter.
    Returns:
        List[StarObject]: A list of star objects for plotting.
    """
    # Split into 12 patches to not run into timeouts
    num_patches = 4
    patch_size = round(limit / num_patches)
    radius_deg = u.Quantity(90, u.deg)

    jobs = []

    # Iterate over all starting points & query with smaller angles
    for index in range(num_patches):
        offset = index * patch_size
        adql_query = f"""
        SELECT TOP {patch_size} ra, dec, parallax, phot_g_mean_mag, bp_rp, teff_gspphot, designation
        FROM gaiadr3.gaia_source
        WHERE CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra}, {dec}, {radius_deg.value})
        ) = 1
        AND parallax IS NOT NULL AND parallax > 0
        AND phot_g_mean_mag <= {min_brightness}
        ORDER BY random_index
        OFFSET {offset}
        """
        print(offset, patch_size, ra, dec, radius_deg.value, min_brightness)

        debug("API", f"Starting query {index + 1} for stars...")

        job = Gaia.launch_job_async(adql_query)
        jobs.append(job)
    
    debug("API", "Fetched all stars, starting to process...")

    objs = []
    for job in jobs:
        stars = job.get_results()

        for idx in range(len(stars)):
            row = stars[idx]

            ra = row['ra']
            dec = row['dec']
            parallax = row['parallax']

            # skip stars without a parallax value
            if not parallax or parallax <= 0:
                continue

            luminosity = row['phot_g_mean_mag']
            bp_rp = row['bp_rp']
            temperature = row['teff_gspphot']

            # calculate temperature, if not provided by request
            if not temperature:
                temperature = 5601 * (0.4 * bp_rp + 1) ** (-1.6)

            distance = 1000 / parallax
            stefan_boltzmann_constant = 5.67e-8  # W/m^2/K^4
            luminosity_in_watt = luminosity * 3.828e26
            nom = luminosity_in_watt
            denom = 4 * np.pi * stefan_boltzmann_constant * temperature**4
            # skip stars with invalid radius
            if not denom:
                continue
            div = nom / denom
            radius = np.sqrt(div)

            x, y, z = spherical_to_cartesian(ra, dec, distance)
            designation = row['DESIGNATION']
            drawobject = StarObject(x, y, z, luminosity, radius, temperature, designation)
            objs.append(drawobject)

    debug("API", f"Done processing, returning {len(objs)} stars")
    return objs