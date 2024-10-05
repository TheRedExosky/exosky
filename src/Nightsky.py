import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astroquery.gaia import Gaia
from matplotlib.colors import to_rgba
import pandas as pd
import sys

from api import DrawObject, run_api

def plot_3d_star_map_with_planet(drawobjects: [DrawObject]):
    # 3D-Abbildung erstellen
    fig = plt.figure(figsize=(10, 8))

    ax = fig.add_subplot(111, projection='3d')

    print(len(drawobjects))

    # clamp lumonisoty from 0 to 1
    min_lum = sys.float_info.max
    max_lum = sys.float_info.min
    for obj in drawobjects:
        if obj.luminosity > max_lum:
            max_lum = obj.luminosity
            print("new max lum", max_lum)
        if obj.luminosity < min_lum:
            min_lum = obj.luminosity
            print("new min lum", min_lum)

    def translate_clamp(value, min_lum, max_lum):
        nom = (value - min_lum)
        denom = (max_lum - min_lum)
        if denom == 0.0:
            return 1.0
        else:
            return nom / denom

    print("Min lum", min_lum)
    print("Max lum", max_lum)
    
    # Durch die Objekte in drawobject iterieren und für jedes die Helligkeit berücksichtigen
    for obj in drawobjects:
        print("clamped", translate_clamp(obj.luminosity,min_lum, max_lum))
        # Helligkeit des Objekts beeinflusst seine Darstellung
        color_with_luminosity = to_rgba('white', alpha=translate_clamp(obj.luminosity,min_lum, max_lum))  # Helligkeit für Transparenz
        size_with_luminosity = 10 * obj.luminosity # Größe proportional zur Helligkeit

        # Objekt als Punktwolke darstellen
        ax.scatter(obj.x, obj.y, obj.z, s=size_with_luminosity, c=[color_with_luminosity])

    # Hintergrund und Titel setzen
    ax.set_facecolor('black')

    # Optional: Einen speziellen Planeten in der Mitte der Szene
    ax.set_title('3D Star Map', color='yellow')

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

    # Achsen entfernen für eine klare Ansicht
    ax.set_axis_off()

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

    # 3D-Sternkarte erstellen
    plot_3d_star_map_with_planet(run_api(exoplanet_ra, exoplanet_dec))

# Exoplanetendaten aus der CSV-Datei laden
file_path = './../PSCompPars_2024.10.04_08.31.39.csv'
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
