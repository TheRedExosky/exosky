import numpy as np

def spherical_to_cartesian(ra, dec, distance):
    """
    ra und dec Koordinaten in Kartesische Koordinaten umrechnen
    """
    ra_rad = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)

    # Umrechnung in kartesische Koordinaten
    x = distance * np.cos(dec_rad) * np.cos(ra_rad)
    y = distance * np.cos(dec_rad) * np.sin(ra_rad)
    z = distance * np.sin(dec_rad)
    
    return x, y, z
