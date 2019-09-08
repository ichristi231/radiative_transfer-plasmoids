# Results
This directory contains text files of the results from the radiative transfer code. In it, we have the log-10 values of:

i) Photon frequency array (in Hz)

ii) Particle (i.e. electron/positron) Lorentz factor array

iii) Temporal evolutino of photon spectrum (i.e. `nu*L_nu` in erg/s)

iv) Temporal evolution of particle distribution (i.e. `gamma^p * N_e(gamma, t)`, where `gamma` is the particle's Lorentz factor and `p` is the slope of the injected power-law distribution).

v) The cumulative number of injected particles, determined from the particle injection rate, along with the total number of particles within the plasmoid, determined by integrating the particle distribution, at every timestep. This is a check to make sure that the code conserves particle number.
