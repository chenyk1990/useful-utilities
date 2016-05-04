#!/bin/sh

gnuplot plot_thermocline.gnu

gnuplot plot_points_per_wavelength_histogram.gnu

python generate_shallow_water_profile_script_by_Paul_Cristini.py


