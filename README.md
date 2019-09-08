This repository contains a single-zone radiative code used for computing the emission from a time-evolving relativistic, non-thermal particle distribution. Currently, this code contains all relevent radiative processes for a leptonic (i.e. electron and positron) emission model (for a more apt description, see []). Here, the emitting region is assumed to be a plasmoid, i.e. a quasi-spherical blob of plasma containing magnetic fields and relativistic particles (for detailed description of plasmoids and their astrophysical application, see []). As such, the code is coupled to results of particle-in-cell simulations of relativistic reconnection in electron-positron (i.e. pair) plasma which govern the evolution and dynamics of plasmoids (for description of PIC and relativistic reconnection, see []). 

If using this code (and any of its contents), please make proper references of the following:

i) Christie I. M., Petropoulou M., Sironi L., Giannios D., 2019, MNRAS, 482, 65 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.482...65C/abstract)

ii) Sironi L., Giannios D., Petropoulou M., 2016, MNRAS, 462, 48 (https://ui.adsabs.harvard.edu/abs/2016MNRAS.462...48S/abstract)

iii) Petropoulou M., Giannios D., Sironi L., 2016, MNRAS, 462, 3325 (https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.3325P/abstract)

To run the contents of the code, in the directoy containing `leptonic.c` type `gcc leptonic.c -o leptonic -O3` to compile and `./leptonic` to run.



# Introduction
An over-arching description of the code is as follows:

i) Set up the initial conditions required to compute the emission from the evolving particle distribution. For this scenario, it includes: `magnetic field of plasma far upstream from the reconnection layer` (see []), `half-length of the reconnection layer` (see []), `slope of the injected particle distribution` (see eqn. and table in []), `maximum Lorentz factor of the injected particle distribution` (see Appendix in []). 

ii) Compute all source and loss rate terms for the evolving particle and photon distributions. 

iii) Numerically solves the continuity equation for the evolving the particle and photon distributions. Solving the former requires the use of the tri-diagonal matrix algorithm [].

iv) Writes to a file the log10 of the following: `particle distribution`, `particle Lorentz factors`, `photon frequencis` (in Hz), and `photon spectrum` (i.e. `nu L_nu` in erg/s).

# More to Come

# Citations

i) 

ii)

iii)
