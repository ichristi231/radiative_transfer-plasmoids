# Plasmoid Properties
This directory contains text files of the plasmoid properties required to calculate the particle evolution and resulting photon spectrum. These properties are the result of particle-in-cell (PIC) simulations of electron-positron pair plasmas, for a plasma magnetization of 10, as reported in `Sironi L., Giannios D., Petropoulou M., 2016, MNRAS, 462, 48` [[1][1]]. As such, any use of these results should require this citation as well as the ones listed below.

Contained within this directory are the plamoid's:

i) Transverse diameter normalized to the half-length of the reconnection layer [[1][1]]. This is used to compute the photon escape timescale and the timestep `delta_t` [[2][2]].

ii) Co-moving (i.e. in the rest frame of the plasmoid) volume (for definition, see Sec. 2 [[2][2]]) normalized to the half-length of the reconnection layer cubed [[1][1]]. This quantitiy is used in several source and loss rate terms for photons and particles. 

iii) The co-moving temporal derivative of the plasmoid's volume. This quantity is used to set the injection rate of particles (see eqn. 5 in [[2][2]]).

iv) Lorentz factor as measured in the reconnection layer's rest frame (for additional description, see [[3][3]]). During the radiative transfer caculation, this quantity is used in the particle injection rate (see eqn. 5 in [[2][2]]) and for determining the co-moving energy of any photon fields external to the plasmoid and/or relativistic blazar jet (e.g. Broad-line region, Infrared Torus).



# Citations

[1] Sironi L., Giannios D., Petropoulou M., 2016, MNRAS, 462, 48 (https://ui.adsabs.harvard.edu/abs/2016MNRAS.462...48S/abstract)

[2] Christie I. M., Petropoulou M., Sironi L., Giannios D., 2019, MNRAS, 482, 65 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.482...65C/abstract)

[3] Petropoulou M., Giannios D., Sironi L., 2016, MNRAS, 462, 3325 (https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.3325P/abstract)

[1]: https://ui.adsabs.harvard.edu/abs/2016MNRAS.462...48S/abstract
[2]: https://ui.adsabs.harvard.edu/abs/2019MNRAS.482...65C/abstract
[3]: https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.3325P/abstract
