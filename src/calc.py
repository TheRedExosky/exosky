"""
Some general calculation work for the simulation.
"""

import numpy as np

from astropy.coordinates import Angle
import astropy.units as u

def spherical_to_cartesian(ra: float, dec: float, distance: float):
    """
    Convert `ra` and `dec` coordinates to cartesian.
    """
    ra_rad = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)

    # Convert to cartesian coordinates
    x = distance * np.cos(dec_rad) * np.cos(ra_rad)
    y = distance * np.cos(dec_rad) * np.sin(ra_rad)
    z = distance * np.sin(dec_rad)
    
    return x, y, z


def hms_to_degrees(ra_hms: float):
    """
    Convert RA in HMS format to decimal degrees.
    """
    ra_angle = Angle(ra_hms, unit=u.hour)
    return ra_angle.deg


def dms_to_degrees(dec_dms: float):
    """
    Convert DEC in DMS format to decimal degrees.
    """
    dec_angle = Angle(dec_dms, unit=u.deg)
    return dec_angle.deg


def temperature_to_color(temperature: float):
    """
    Temperature [K] to light wave length [nm].
    """
    return round(0.002898 / temperature * 1000000000)


def wavelength_to_rgb(wavelength: int):
    """
    Light wave length [nm] to tuple of rgb values.
    """
    gamma = 0.8
    intensity_max = 255
    
    if 380 <= wavelength <= 440:
        r = -(wavelength - 440) / (440 - 380)
        g = 0.0
        b = 1.0
    elif 440 < wavelength <= 490:
        r = 0.0
        g = (wavelength - 440) / (490 - 440)
        b = 1.0
    elif 490 < wavelength <= 510:
        r = 0.0
        g = 1.0
        b = -(wavelength - 510) / (510 - 490)
    elif 510 < wavelength <= 580:
        r = (wavelength - 510) / (580 - 510)
        g = 1.0
        b = 0.0
    elif 580 < wavelength <= 645:
        r = 1.0
        g = -(wavelength - 645) / (645 - 580)
        b = 0.0
    elif 645 < wavelength <= 780:
        r = 1.0
        g = 0.0
        b = 0.0
    else:
        # Wavelength outside the visible spectrum
        r = g = b = 0.0

    # Let the intensity fall off near the vision limits
    if 380 <= wavelength <= 420:
        factor = 0.3 + 0.7 * (wavelength - 380) / (420 - 380)
    elif 645 < wavelength <= 780:
        factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 645)
    else:
        factor = 1.0

    r = int(intensity_max * (r * factor) ** gamma)
    g = int(intensity_max * (g * factor) ** gamma)
    b = int(intensity_max * (b * factor) ** gamma)

    return r, g, b