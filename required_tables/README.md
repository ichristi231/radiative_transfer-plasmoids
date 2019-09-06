# Required Tables
This directory contains text files needed to compute the different production or loss rates for a variety of radiative processes. Include in here is the following (with appropriate reference and description):

i) `single_part_emissivity_x_store.txt` & `single_part_emissivity_f_store.txt`:

The first file stores the log-10 values of the ratio of the photon frequency to the critical synchrotron frequency (see eqn. 6.17c in [1]) used in the definition of `F(x)` (eqn. 6.31c in [1]) for the single-particle synchrotron emissivity (eqn. 6.33 in [1]). The second file stores the log-10 values of `F(x)` (eqn. 6.31c in [1]) for the values store in `x`. 

The values in `x` are chosen such that the code will not need to manually calculate `F(x)`. Instead, we interpolate `F(x)` such that we can maually select the values from the array.

ii) `..._production_numerical_values_photo_hadronic_interactions.txt`:

All files ending with the text listed above contain the production rate of differen particle species as produced from photo-hadronic interactions (see [2] for more detail). The `...` can contain any of the following particle species: `electron_neutrino`, `anti-electron_neutrino`, `muon_neutrino`, `anti_muon_neutrino`, `electron`, `positron`, or `photon`. Each of these text files contain four columns, with their variable representation described in detail in [2]. We note, that the fourth and final column are the log-10 of the production rate of that particular species. Moreover, these will only be used in the lepto_hadronic.c version of the code.


# Citations

[1] Rybicki G. B., Lightman A. P., 1986, Radiative Processes in Astrophysics. p. 400

[2] Kelner, S. R., & Aharonian, F. A. 2008, PhRvD, 78, 034013
(https://arxiv.org/pdf/0803.0688.pdf)
