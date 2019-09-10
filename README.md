This repository contains a single-zone radiative code used for computing the emission from a time-evolving relativistic, non-thermal particle distribution. Currently, this code contains all relevent radiative processes for a leptonic (i.e. electron and positron) emission model (for a more apt description, see [[3][3]]). Here, the emitting region is assumed to be a plasmoid, i.e. a quasi-spherical blob of plasma containing magnetic fields and relativistic particles (for detailed description of plasmoids and their astrophysical application, see [4]). As such, the code is coupled to results of particle-in-cell simulations of relativistic reconnection in electron-positron (i.e. pair) plasma which govern the evolution and dynamics of plasmoids (for description of PIC and relativistic reconnection, see [6]). 

If using this code (and any of its contents), please make proper references of the following:

i) Christie I. M., Petropoulou M., Sironi L., Giannios D., 2019, MNRAS, 482, 65 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.482...65C/abstract)

ii) Sironi L., Giannios D., Petropoulou M., 2016, MNRAS, 462, 48 (https://ui.adsabs.harvard.edu/abs/2016MNRAS.462...48S/abstract)

To run the contents of the code, in the directory containing `leptonic.c` type `gcc leptonic.c -o leptonic -O3` to compile and `./leptonic` to run. Possible issues might arise and may require linking to the math library. If this errors occurs, compile via `gcc leptonic.c -o leptonic -O3 -lm`.



# Introduction
An over-arching description of the code is as follows:

i) Set up the initial conditions required to compute the emission from the evolving particle distribution. For this scenario, it includes: `magnetic field of plasma far upstream from the reconnection layer` (see Fig. 5 in [6]), `half-length of the reconnection layer` (see Sec 2 and table 1 [3]), `slope of the injected particle distribution` (see eqn. 5 and table 1 in [3]), `maximum Lorentz factor of the injected particle distribution` (see Appendix A in [3]). 

ii) Computes all source and loss rate terms for the evolving particle and photon distributions. 

iii) Numerically solves the continuity equation for the evolving the particle and photon distributions. Solving the former requires the use of the tri-diagonal matrix algorithm [1, 2, 5].

iv) Writes to a file the log10 of the following: `particle distribution` (i.e. `gamma^p * N(gamma, t)`, where `p` is the slope of the injected particle distribution), `particle Lorentz factors`, `photon frequencies` (in Hz), and `photon spectrum` (i.e. `nu L_nu` in erg/s).

# More to Come

Development of a lepto-hadronic code which allows for the computation of the potential neutrino flux from blazars!

# Citations

(1) Chang, J.S. & Cooper, G., 1970, Journal of Computational Physics 6, 1 (https://ui.adsabs.harvard.edu/abs/1970JCoPh...6....1C/abstract)

(2) Chiaberge M., Ghisellini G., 1999, MNRAS, 306, 551 (https://ui.adsabs.harvard.edu/abs/1999MNRAS.306..551C/abstract)

(3) Christie I. M., Petropoulou M., Sironi L., Giannios D., 2019, MNRAS, 482, 65 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.482...65C/abstract)

(4) Petropoulou M., Giannios D., Sironi L., 2016, MNRAS, 462, 3325 (https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.3325P/abstract)

(5) Press, W.H. et al., 1989, Numerical Recipes in Fortran, Cambridge University Press (https://ui.adsabs.harvard.edu/abs/1989nrpa.book.....P/abstract)

(6) Sironi L., Giannios D., Petropoulou M., 2016, MNRAS, 462, 48 (https://ui.adsabs.harvard.edu/abs/2016MNRAS.462...48S/abstract)

[1]: https://ui.adsabs.harvard.edu/abs/1970JCoPh...6....1C/abstract
[2]: https://ui.adsabs.harvard.edu/abs/1999MNRAS.306..551C/abstract
[3]: https://ui.adsabs.harvard.edu/abs/2019MNRAS.482...65C/abstract
[4]: https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.3325P/abstract
[5]: https://ui.adsabs.harvard.edu/abs/1989nrpa.book.....P/abstract
[6]: https://ui.adsabs.harvard.edu/abs/2016MNRAS.462...48S/abstract

