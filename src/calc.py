import numpy as np

import pandas as pd
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

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





# Beispiel-Datenrahmen erstellen
data = {
    'dist': [1.0, 2.0, 3.0],
    'solution_id': [1, 2, 3],
    'DESIGNATION': ['A', 'B', 'C'],
    'SOURCE_ID': [1001, 1002, 1003],
    'random_index': [10, 20, 30],
    'ref_epoch': [2000, 2001, 2002],
    'ra': [10.0, 20.0, 30.0],
    'ra_error': [0.1, 0.2, 0.3],
    'dec': [10.0, 20.0, 30.0],
    'dec_error': [0.1, 0.2, 0.3],
    'parallax': [0.1, 0.2, 0.3],
    'parallax_error': [0.01, 0.02, 0.03],
    'parallax_over_error': [10, 10, 10],
    'pm': [1.0, 2.0, 3.0],
    'pmra': [0.1, 0.2, 0.3],
    'pmra_error': [0.01, 0.02, 0.03],
    'pmdec': [0.1, 0.2, 0.3],
    'pmdec_error': [0.01, 0.02, 0.03],
    'ra_dec_corr': [0.1, 0.2, 0.3],
    'ra_parallax_corr': [0.1, 0.2, 0.3],
    'ra_pmra_corr': [0.1, 0.2, 0.3],
    'ra_pmdec_corr': [0.1, 0.2, 0.3],
    'dec_parallax_corr': [0.1, 0.2, 0.3],
    'dec_pmra_corr': [0.1, 0.2, 0.3],
    'dec_pmdec_corr': [0.1, 0.2, 0.3],
    'parallax_pmra_corr': [0.1, 0.2, 0.3],
    'parallax_pmdec_corr': [0.1, 0.2, 0.3],
    'pmra_pmdec_corr': [0.1, 0.2, 0.3],
    'astrometric_n_obs_al': [10, 20, 30],
    'astrometric_n_obs_ac': [10, 20, 30],
    'astrometric_n_good_obs_al': [10, 20, 30],
    'astrometric_n_bad_obs_al': [0, 0, 0],
    'astrometric_gof_al': [1.0, 2.0, 3.0],
    'astrometric_chi2_al': [1.0, 2.0, 3.0],
    'astrometric_excess_noise': [0.1, 0.2, 0.3],
    'astrometric_excess_noise_sig': [0.1, 0.2, 0.3],
    'astrometric_params_solved': [1, 1, 1],
    'astrometric_primary_flag': [True, True, True],
    'nu_eff_used_in_astrometry': [1.0, 2.0, 3.0],
    'pseudocolour': [1.0, 2.0, 3.0],
    'pseudocolour_error': [0.1, 0.2, 0.3],
    'ra_pseudocolour_corr': [0.1, 0.2, 0.3],
    'dec_pseudocolour_corr': [0.1, 0.2, 0.3],
    'parallax_pseudocolour_corr': [0.1, 0.2, 0.3],
    'pmra_pseudocolour_corr': [0.1, 0.2, 0.3],
    'pmdec_pseudocolour_corr': [0.1, 0.2, 0.3],
    'astrometric_matched_transits': [10, 20, 30],
    'visibility_periods_used': [10, 20, 30],
    'astrometric_sigma5d_max': [1.0, 2.0, 3.0],
    'matched_transits': [10, 20, 30],
    'new_matched_transits': [10, 20, 30],
    'matched_transits_removed': [0, 0, 0],
    'ipd_gof_harmonic_amplitude': [1.0, 2.0, 3.0],
    'ipd_gof_harmonic_phase': [1.0, 2.0, 3.0],
    'ipd_frac_multi_peak': [0.1, 0.2, 0.3],
    'ipd_frac_odd_win': [0.1, 0.2, 0.3],
    'ruwe': [1.0, 2.0, 3.0],
    'scan_direction_strength_k1': [1.0, 2.0, 3.0],
    'scan_direction_strength_k2': [1.0, 2.0, 3.0],
    'scan_direction_strength_k3': [1.0, 2.0, 3.0],
    'scan_direction_strength_k4': [1.0, 2.0, 3.0],
    'scan_direction_mean_k1': [1.0, 2.0, 3.0],
    'scan_direction_mean_k2': [1.0, 2.0, 3.0],
    'scan_direction_mean_k3': [1.0, 2.0, 3.0],
    'scan_direction_mean_k4': [1.0, 2.0, 3.0],
    'duplicated_source': [False, False, False],
    'phot_g_n_obs': [10, 20, 30],
    'phot_g_mean_flux': [1.0, 2.0, 3.0],
    'phot_g_mean_flux_error': [0.1, 0.2, 0.3],
    'phot_g_mean_flux_over_error': [10, 10, 10],
    'phot_g_mean_mag': [1.0, 2.0, 3.0],
    'phot_bp_n_obs': [10, 20, 30],
    'phot_bp_mean_flux': [1.0, 2.0, 3.0],
    'phot_bp_mean_flux_error': [0.1, 0.2, 0.3],
    'phot_bp_mean_flux_over_error': [10, 10, 10],
    'phot_bp_mean_mag': [1.0, 2.0, 3.0],
    'phot_rp_n_obs': [10, 20, 30],
    'phot_rp_mean_flux': [1.0, 2.0, 3.0],
    'phot_rp_mean_flux_error': [0.1, 0.2, 0.3],
    'phot_rp_mean_flux_over_error': [10, 10, 10],
    'phot_rp_mean_mag': [1.0, 2.0, 3.0],
    'phot_bp_rp_excess_factor': [1.0, 2.0, 3.0],
    'phot_bp_n_contaminated_transits': [0, 0, 0]
}

df = pd.DataFrame(data)

def berechne_luminositaet(temperatur, radius):
    sigma = 5.67e-8  # Stefan-Boltzmann-Konstante
    luminositaet = 4 * 3.14159 * (radius ** 2) * sigma * (temperatur ** 4)
    return luminositaet

def farbe_aus_luminositaet(luminositaet, min_lum, max_lum):
    norm_lum = (luminositaet - min_lum) / (max_lum - min_lum)
    farbe = sns.color_palette("viridis", as_cmap=True)(norm_lum)
    return mcolors.to_hex(farbe)

# Beispielwerte für Temperatur (in Kelvin) und Radius (in Metern)
temperatur = 5778  # Temperatur der Sonne
radius = 6.96e8  # Radius der Sonne

# Berechnung der Luminosität für jeden Eintrag im DataFrame
df['luminositaet'] = df.apply(lambda row: berechne_luminositaet(temperatur, radius), axis=1)
