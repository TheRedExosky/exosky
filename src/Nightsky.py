"""
Main Entry Point to the whole simulation.
Just run this as `python3 NightSky.py`.
"""

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import tkinter as tk
from tkinter import messagebox
from tkinter import HORIZONTAL
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import pandas as pd
import sys
from typing import List

from api import StarObject, fetch_api
from calc import *

def plot_stars(stars: List[StarObject]):
    """
    Plot all stars, retrieved by the API.
    """
    print("Amount start", len(stars))

    # setup plot and subplot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')
    ax.set_title('3D Star Map', color='yellow')
    ax.set_xlim([-2000, 2000])
    ax.set_ylim([-2000, 2000])
    ax.set_zlim([-2000, 2000])
    ax.set_axis_off()
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # exoplanet name
    selected_index = exoplanet_listbox.curselection()
    selected_exoplanet = exoplanet_data.iloc[selected_index[0]]
    exoplanet_name = selected_exoplanet[0]
    fig.suptitle("Exoplanet: " + exoplanet_name, color="white")

    # clamp luminosity from 0 to 1 for alpha value in plot
    min_lum = sys.float_info.max
    max_lum = sys.float_info.min
    for star in stars:
        if star.luminosity > max_lum:
            max_lum = star.luminosity
        if star.luminosity < min_lum:
            min_lum = star.luminosity

    def alpha_clamp(value, min_lum, max_lum):
        """
        Clamp values from `min_lum` to `max_lum` to 0 - 1.
        """
        nom = (value - min_lum)
        denom = (max_lum - min_lum)
        if denom == 0.0:
            return 1.0
        else:
            return nom / denom

    print("Min lum", min_lum)
    print("Max lum", max_lum)
    
    # draw each star
    for star in stars:
        size = star.radius * 0.00000001
        wavelength = temperature_to_color(star.temperature)
        rgb_value = wavelength_to_rgb(wavelength)
        normalized_rgb = tuple([x / 255.0 for x in rgb_value])

        alpha_clamp_value = alpha_clamp(star.luminosity, min_lum, max_lum) 
        if not alpha_clamp_value:
            alpha_clamp_value = 0

        color_with_luminosity = to_rgba('white', alpha=alpha_clamp_value)

        # display object as a point cloud
        ax.scatter(star.x, star.y, star.z, s=size, c=[color_with_luminosity])

    # draw planet in center
    if True:
        u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
        r_planet = 200  # Adjust the size of the planet
        x_planet = r_planet * np.cos(u) * np.sin(v)
        y_planet = r_planet * np.sin(u) * np.sin(v)
        z_planet = r_planet * np.cos(v)
        ax.plot_surface(x_planet, y_planet, z_planet, color='blue', alpha=0.7, rstride=5, cstride=5)

    # interactive rotation
    mng = plt.get_current_fig_manager()
    try:
        mng.window.state('zoomed')
    except:
        mng.full_screen_toggle()

    # show all to screen
    plt.show()


def select_exoplanet():
    """
    Select and exoplanet and parse `ra` and `dec` data.
    """
    index_name = 0
    index_ra = 28
    index_dec = 30

    selected_index = exoplanet_listbox.curselection()
    if not selected_index:
        messagebox.showwarning("Selection Error", "Please select an exoplanet.")
        return

    selected_exoplanet = exoplanet_data.iloc[selected_index[0]]
    exoplanet_name = selected_exoplanet[index_name]
    print("selected planet", exoplanet_name)

    # convert to hms
    exoplanet_ra_hms = selected_exoplanet[index_ra]
    exoplanet_dec_dms = selected_exoplanet[index_dec]
    
    # convert to decimal values
    exoplanet_ra = hms_to_degrees(exoplanet_ra_hms)
    exoplanet_dec = dms_to_degrees(exoplanet_dec_dms)

    return exoplanet_ra,exoplanet_dec


def run_simulation():
    """
    Run star simulation
    """
    ra, dec = select_exoplanet()
    plot_stars(fetch_api(ra, dec))


if __name__ == "__main__":
    file_path = './../PSCompPars_2024.10.04_08.31.39.csv'
    exoplanet_data = pd.read_csv(file_path, skiprows=45)
    exoplanet_data.columns = exoplanet_data.columns.str.strip()

    # selection interface
    root = tk.Tk()
    root.title("Star Viewer")
    exoplanet_listbox = tk.Listbox(root, width=50, height=20)
    exoplanet_listbox.pack(padx=10, pady=10)

    # add exoplantes to selection
    for name in exoplanet_data['75 Cet b']:
        exoplanet_listbox.insert(tk.END, name)

    # number of stars slider
    nr_stars_label = tk.Label(root, text="Nr. of Stars")
    nr_stars_label.pack()
    nr_stars_slider = tk.Scale(root, orient=HORIZONTAL, length=400, from_=150, to=2000)
    nr_stars_slider.pack()

    # add start button
    select_button = tk.Button(root, text="Show stars ⭐", command=run_simulation)
    select_button.pack(pady=10)

    root.bind('<Escape>', lambda e: root.destroy())
    root.mainloop()
