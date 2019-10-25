import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from leptonic_version_results import (electron_distribution,
                                      co_moving_lumionsity,
                                      electron_Lorentz_factor,
                                      photon_frequency)


# The function below creates combined subplots of the
# temporal evolution of the co-moving electron and photon
# (i.e. spectra) distributions within the plasmoid. The user
# provides the number of equally-spaced times to plot. Note
# that color coding (going from purple to red) denotes increasing
# time. Here, the spectra is in units of erg/s while for the
# electron distribution, we plot the electron distribution times
# the electron Lorentz factor to the injected power "p" (i.e.
# g^p * N^e).
def single_plasmoid_results(num):
    # Determine the spacing in which to plot the time-evolving
    # distributions.
    viz_index = int(np.shape(electron_distribution)[0]/num)

    # Array of temporal indices.
    viz_array = viz_index * np.arange(1, num)

    # Color map of both plots.
    colors = cm.rainbow(np.linspace(0, 1, num))

    # Define the figure and subplots.
    fig, axs = plt.subplots(1, 2, figsize=(17, 8))

    # For each of the number of curves with a given subplot,
    # we plots the electron and photon distributions.
    for i in range(num - 1):
        # Plot the electron distribution
        axs[0].loglog(10**electron_Lorentz_factor,
                      10**electron_distribution[viz_array[i]], color=colors[i])

        # Plot the spectra
        axs[1].loglog(10**photon_frequency,
                      10**co_moving_lumionsity[viz_array[i]], color=colors[i])

        # x and y-axis ranges
        axs[1].set_ylim([10**35., 10**41.])
        axs[1].set_xlim([10**8., 10**32.])

        # Labels of the first subplot (electron distribution)
        axs[0].set_xlabel(r'$\gamma$', fontsize=15)
        axs[0].set_ylabel(r'$\gamma^p$ N$^e$($\gamma , t$)', fontsize=15)

        # Labels of the second subplot (photon distribution)
        axs[1].set_xlabel(r'$\nu$ (Hz)', fontsize=15)
        axs[1].set_ylabel(r'$\nu \, L_\nu$ (erg/s)', fontsize=15)

    plt.savefig('leptonic_version_results/co_moving_electron_photon_distribution.png')


# We define the number of curves, equally spaced in time,
# to plot the co-moving electron and photon distributions.
num = 10

# Call the function to create plots of the elctron and
# photon distributions.
single_plasmoid_results(num)
