import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astroquery.gaia import Gaia
import pandas as pd

# Funktion zur Umrechnung von Kugelkoordinaten in kartesische Koordinaten
def spherical_to_cartesian(ra, dec, distance):
    ra_rad = np.deg2rad(ra)  # Rektaszension in Radiant umwandeln
    dec_rad = np.deg2rad(dec)  # Deklination in Radiant umwandeln

    # Umrechnung in kartesische Koordinaten
    x = distance * np.cos(dec_rad) * np.cos(ra_rad)
    y = distance * np.cos(dec_rad) * np.sin(ra_rad)
    z = distance * np.sin(dec_rad)
    
    return x, y, z

# Funktion zur Erstellung einer 3D-Sternkarte mit einem Planeten in der Mitte
def plot_3d_star_map_with_planet(ra, dec, parallax, exoplanet_name):
    # Parallaxe in Distanz umrechnen (in Parsecs)
    distance = 1000 / parallax  # Umrechnung der Parallaxwerte

    # Filtere sehr große Distanzen heraus (realistische Grenzen setzen)
    valid_indices = np.isfinite(distance) & (distance < 5000) & (distance > 0)
    ra = ra[valid_indices]
    dec = dec[valid_indices]
    distance = distance[valid_indices]

    # Umrechnung der Kugelkoordinaten in kartesische Koordinaten
    x, y, z = spherical_to_cartesian(ra, dec, distance)

    # 3D-Abbildung erstellen
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Sterne als Punktwolke darstellen
    ax.scatter(x, y, z, s=10, c='white', alpha=0.8)

    # Hintergrund und Titel setzen
    ax.set_facecolor('black')
    ax.set_title(f'3D Star Map from {exoplanet_name}', color='white')

    # Achsen entfernen für eine klare Ansicht
    ax.set_axis_off()

    # Einen größeren Planeten in der Mitte platzieren
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
    r_planet = 200  # Größe des Planeten anpassen
    x_planet = r_planet * np.cos(u) * np.sin(v)
    y_planet = r_planet * np.sin(u) * np.sin(v)
    z_planet = r_planet * np.cos(v)

    # Planeten in der Mitte plotten
    ax.plot_surface(x_planet, y_planet, z_planet, color='blue', alpha=0.7, rstride=5, cstride=5)

    # Achsenlimits für bessere Darstellung setzen
    ax.set_xlim([-2000, 2000])
    ax.set_ylim([-2000, 2000])
    ax.set_zlim([-2000, 2000])

    # Überflüssige Ränder entfernen
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # Fenster maximieren
    mng = plt.get_current_fig_manager()
    try:
        mng.window.state('zoomed')  # Maximiert das Fenster auf macOS
    except:
        mng.full_screen_toggle()  # Alternativ auf Linux/Windows

    # Interaktive Rotation ermöglichen
    plt.show()

# Funktion zur Abfrage von Gaia-Sterndaten
def fetch_gaia_stars_3d(ra, dec, radius=5):
    query = f"""
    SELECT TOP 1500 ra, dec, parallax
    FROM gaiadr3.gaia_source
    WHERE CONTAINS(POINT('ICRS', ra, dec),
                   CIRCLE('ICRS', {ra}, {dec}, {radius}))=1
    AND parallax > 0
    """
    job = Gaia.launch_job(query)
    result = job.get_results()
    return result['ra'], result['dec'], result['parallax']

# Funktion zum Auswählen eines Exoplaneten und Erstellen der Sternkarte
def on_select_exoplanet_3d():
    selected_index = exoplanet_listbox.curselection()
    if not selected_index:
        messagebox.showwarning("Selection Error", "Please select an exoplanet.")
        return

    # Details des ausgewählten Exoplaneten holen
    selected_exoplanet = exoplanet_data.iloc[selected_index[0]]
    exoplanet_name = selected_exoplanet['75 Cet b']  # Exoplanetenname
    exoplanet_ra = selected_exoplanet['38.0391598']  # Dezimal-RA-Spalte
    exoplanet_dec = selected_exoplanet['-1.0350269']  # Dezimal-Dec-Spalte

    # Gaia-Sterne um die Koordinaten des Exoplaneten abfragen
    ra, dec, parallax = fetch_gaia_stars_3d(exoplanet_ra, exoplanet_dec, radius=10)

    # 3D-Sternkarte erstellen
    plot_3d_star_map_with_planet(ra, dec, parallax, exoplanet_name)

# Exoplanetendaten aus der CSV-Datei laden
file_path = 'PSCompPars_2024.10.04_08.31.39.csv'
exoplanet_data = pd.read_csv(file_path, skiprows=45)

# Spaltennamen bereinigen
exoplanet_data.columns = exoplanet_data.columns.str.strip()  # Entfernt überflüssige Leerzeichen

# Tkinter GUI für die Exoplanetenliste
root = tk.Tk()
root.title("3D Star Map Viewer")

# Listbox für die Exoplanetenanzeige
exoplanet_listbox = tk.Listbox(root, width=50, height=20)
exoplanet_listbox.pack(padx=10, pady=10)

# Exoplanetennamen zur Liste hinzufügen
for name in exoplanet_data['75 Cet b']:  # Exoplanetennamen
    exoplanet_listbox.insert(tk.END, name)

# Button zur Generierung der Sternkarte
select_button = tk.Button(root, text="Show 3D Star Map", command=on_select_exoplanet_3d)
select_button.pack(pady=10)

# Tkinter-Hauptschleife starten
root.mainloop()
