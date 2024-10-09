"""
Main Entry Point to the whole simulation.
Just run this as `python3 NightSky.py`.
"""

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import tkinter as tk
from tkinter import messagebox
from tkinter import HORIZONTAL
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import pandas as pd
import sys
from typing import List
import mplcursors

from api import StarObject, fetch_api
from calc import *
from debug import debug

if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "agg":
        import matplotlib
        matplotlib.use('Qt5Agg')


def plot_stars(stars: List[StarObject]):
    """
    Plot all stars, retrieved by the API.
    """
    # setup plot and subplot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')
    ax.set_title('3D Star Map', color='black')
    ax.set_xlim([-2000, 2000])
    ax.set_ylim([-2000, 2000])
    ax.set_zlim([-2000, 2000])
    ax.set_axis_off()
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # exoplanet name
    selected_index = exoplanet_listbox.curselection()
    selected_exoplanet = exoplanet_data.iloc[selected_index[0]]
    exoplanet_name = selected_exoplanet.iloc[0]
    fig.suptitle("Exoplanet: " + exoplanet_name, color="white")

    # clamp luminosity from 0 to 1 for alpha value in plot
    min_lum = sys.float_info.max
    max_lum = sys.float_info.min
    for star in stars:
        if star.luminosity > max_lum:
            max_lum = star.luminosity
        if star.luminosity < min_lum:
            min_lum = star.luminosity

    def alpha_clamp(value, min_lum_x, max_lum_x):
        """
        Clamp values from `min_lum` to `max_lum` to 0 - 1.
        """
        nom = (value - min_lum_x)
        denom = (max_lum_x - min_lum_x)
        if denom == 0.0:
            return 1.0
        else:
            return nom / denom

    debug("Minimum found luminosity", str(min_lum))
    debug("Maximum found luminosity", str(max_lum))
    
    # draw each star
    xs = []
    ys = []
    zs = []
    colors = []
    radia = []
    sizes = []
    temperatures = []
    for idx,star in enumerate(stars):
        size = star.radius * 1e-9
        temperatures.append(star.temperature)
        wavelength = temperature_to_color(star.temperature)
        rgb_value = wavelength_to_rgb(wavelength)
        normalized_rgb = tuple([x / 255.0 for x in rgb_value])

        alpha_clamp_value = alpha_clamp(star.luminosity, min_lum, max_lum) 
        if not alpha_clamp_value:
            alpha_clamp_value = 0
        # higher luminosity implies lower brightness, resulting in lower alpha
        alpha_clamp_value = 1 - alpha_clamp_value

        color_with_luminosity = to_rgba('white', alpha=alpha_clamp_value)

        # display object as a point cloud
        stretch = 10
        xs.append(star.x * stretch)
        ys.append(star.y * stretch)
        zs.append(star.z * stretch)
        radia.append(star.radius)
        sizes.append(size)
        colors.append(color_with_luminosity)

    # show stars to screen
    scatter = ax.scatter(xs,ys,zs, s=sizes, c=colors)
    names = [s.designation for s in stars]

    # add hover information
    mplcursors.cursor(scatter, hover=2).connect("add", lambda sel: sel.annotation.set_text(f"Name: {names[sel.index]}\nRadius: {radia[sel.index]}m\nTemperatur: {temperatures[sel.index]}K"))

    # draw planet in center
    if True:
        c, v = np.mgrid[0:2 * np.pi:50j, 0:np.pi:25j]
        r_planet = 50
        x_planet = r_planet * np.cos(c) * np.sin(v)
        y_planet = r_planet * np.sin(c) * np.sin(v)
        z_planet = r_planet * np.cos(v) * 1.3
        ax.plot_surface(x_planet, y_planet, z_planet, color='#006994', alpha=1, rstride=5, cstride=5)

    # interactive rotation
    mng = plt.get_current_fig_manager()
    try:
        mng.window.state('zoomed')
    except:
        mng.full_screen_toggle()
    
    mng.set_window_title("TheRedExoplanet")

    # show all to screen
    plt.show()


def select_exoplanet():
    """
    Select and exoplanet and parse `ra` and `dec` data.
    """
    index_name = 0
    index_ra = 1
    index_dec = 2

    selected_index = exoplanet_listbox.curselection()
    if not selected_index:
        messagebox.showwarning("Selection Error", "Please select an exoplanet.")
        return

    selected_exoplanet = exoplanet_data.iloc[selected_index[0]]
    exoplanet_name = selected_exoplanet.iloc[index_name]
    debug("Selected planet", exoplanet_name)

    # extract right ascension & declination
    exoplanet_ra = selected_exoplanet.iloc[index_ra]
    exoplanet_dec = selected_exoplanet.iloc[index_dec]

    return exoplanet_ra,exoplanet_dec


def run_simulation():
    """
    Run star simulation
    """
    ra, dec = select_exoplanet()
    plot_stars(fetch_api(ra, dec, round(nr_stars_slider.get()), round(min_brightness_slider.get())))


def update_approximation(ignored):
    """
    Callback method for the nr stars slider,
    provides the current value as string, ignored because useless,
    sets the approximation text
    """
    seconds = round((27 / 7122500) * nr_stars_slider.get()**2 + (-391 / 142450) * nr_stars_slider.get() + (114890 /2849))
    approximation_label.config(text=f"Approximately {seconds} seconds")


if __name__ == "__main__":
    file_path = './../exoplanets.csv'
    exoplanet_data = pd.read_csv(file_path, skiprows=0)
    exoplanet_data.columns = exoplanet_data.columns.str.strip()

    # selection interface
    root = tk.Tk()
    root.title("Star Viewer")
    exoplanet_listbox = tk.Listbox(root, width=50, height=20)
    exoplanet_listbox.pack(padx=10, pady=10)

    # add exoplanets to selection
    for name in exoplanet_data['pl_name']:
        exoplanet_listbox.insert(tk.END, name)

    # number of stars slider
    nr_stars_slider = tk.Scale(root, orient=HORIZONTAL, length=400, from_=150, to=4000, command=update_approximation)
    nr_stars_slider.pack()
    nr_stars_label = tk.Label(root, text="Nr. of Stars")
    nr_stars_label.pack()

    # approximation text
    approximation_label = tk.Label(root, text="Approximately 40 seconds")
    approximation_label.pack()

    # min brightness slider
    min_brightness_slider = tk.Scale(root, orient=HORIZONTAL, length=400, from_=3, to=21)
    min_brightness_slider.set(21)
    min_brightness_slider.pack()
    min_brightness_label = tk.Label(root, text="Minimum Star Luminosity (Small means brighter)")
    min_brightness_label.pack()

    # add start button
    select_button = tk.Button(root, text="Show stars ⭐", command=run_simulation)
    select_button.pack(pady=10)

    root.bind('<Escape>', lambda e: root.destroy())
    root.mainloop()
