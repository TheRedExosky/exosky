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
from astropy.coordinates import Angle
import astropy.units as u
from calc import temperature_to_color, wavelength_to_rgb

def plot_3d_star_map_with_planet(drawobjects: (DrawObject)):

    # Create 3D illustration
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
    
    # Iterate through the objects in drawobject and take the brightness into account for each one
    for obj in drawobjects:
        print("clamped", translate_clamp(obj.luminosity,min_lum, max_lum))
        # Brightness of the object affects its appearance
        size = obj.radius * 0.00000001
        wavelength = temperature_to_color(obj.temp)
        rgb_value = wavelength_to_rgb(wavelength)
        normalized_rgb = tuple([x / 255.0 for x in rgb_value])

        color_with_luminosity = to_rgba(normalized_rgb, alpha=translate_clamp(obj.luminosity, min_lum, max_lum))

        # Display object as a point cloud
        ax.scatter(obj.x, obj.y, obj.z, s=size, c=[color_with_luminosity])

    # Set background and title
    ax.set_facecolor('black')

    # Optional: A special planet in the middle of the scene
    ax.set_title('3D Star Map', color='yellow')

    # Place a larger planet in the center
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
    r_planet = 200  # Adjust the size of the planet
    x_planet = r_planet * np.cos(u) * np.sin(v)
    y_planet = r_planet * np.sin(u) * np.sin(v)
    z_planet = r_planet * np.cos(v)

    # Plot planets in the middle
    ax.plot_surface(x_planet, y_planet, z_planet, color='blue', alpha=0.7, rstride=5, cstride=5)

    # Set axis limits for better display
    ax.set_xlim([-2000, 2000])
    ax.set_ylim([-2000, 2000])
    ax.set_zlim([-2000, 2000])

    # Remove axes for a clear view
    ax.set_axis_off()

    # Remove unnecessary edges
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # Maximize windows
    mng = plt.get_current_fig_manager()
    try:
        mng.window.state('zoomed')  # Maximizes the window on macOS
    except:
        mng.full_screen_toggle()  # Enable interactive rotation

    # Enable interactive rotation
    selected_index = exoplanet_listbox.curselection()
    selected_exoplanet = exoplanet_data.iloc[selected_index[0]]
    exoplanet_name = selected_exoplanet[0]

    fig.suptitle("Exoplanet: " + exoplanet_name, color="white")
    plt.show()


# Function to convert RA in HMS format to decimal degrees
def hms_to_degrees(ra_hms):
    ra_angle = Angle(ra_hms, unit=u.hour)
    return ra_angle.deg

# Function to convert DEC in DMS format to decimal degrees
def dms_to_degrees(dec_dms):
    dec_angle = Angle(dec_dms, unit=u.deg)
    return dec_angle.deg

# Function to select an exoplanet and create the star map
def on_select_exoplanet_3d():
    index_name = 0
    index_ra = 28
    index_dec = 30

    selected_index = exoplanet_listbox.curselection()
    if not selected_index:
        messagebox.showwarning("Selection Error", "Please select an exoplanet.")
        return

    # Get details of the selected exoplanet
    selected_exoplanet = exoplanet_data.iloc[selected_index[0]]
    exoplanet_name = selected_exoplanet[index_name]  # Get exoplanet name
    print(exoplanet_name)

    exoplanet_ra_hms = selected_exoplanet[index_ra]  # RA column in HMS format
    exoplanet_dec_dms = selected_exoplanet[index_dec]  # DEC column in DMS format
    
    # Conversion to decimal values
    exoplanet_ra = hms_to_degrees(exoplanet_ra_hms)  # RA to decimal degrees
    exoplanet_dec = dms_to_degrees(exoplanet_dec_dms)  # DEC to decimal degrees

    # Create 3D star map
    plot_3d_star_map_with_planet(run_api(exoplanet_ra, exoplanet_dec))

# Load exoplanet data from CSV file
file_path = './../PSCompPars_2024.10.04_08.31.39.csv'
exoplanet_data = pd.read_csv(file_path, skiprows=45)

# Clean up column names
exoplanet_data.columns = exoplanet_data.columns.str.strip()  # Removes unnecessary spaces

# Tkinter GUI for the exoplanet list
root = tk.Tk()
root.title("3D Star Map Viewer")

# Listbox for the exoplanet display
exoplanet_listbox = tk.Listbox(root, width=50, height=20)
exoplanet_listbox.pack(padx=10, pady=10)

# Add exoplanet names to the list
for name in exoplanet_data['75 Cet b']:  # Exoplanet names
    exoplanet_listbox.insert(tk.END, name)

# Button for generating the star map
select_button = tk.Button(root, text="Show 3D Star Map", command=on_select_exoplanet_3d)
select_button.pack(pady=10)

# Start Tkinter main loop
root.mainloop()
