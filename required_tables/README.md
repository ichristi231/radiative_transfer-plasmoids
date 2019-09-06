# Required Tables
This directory contains text files needed to compute the differen production or loss rates for a variety of processes. Include in here is the following (with appropriate reference and description):

i) `single_part_emissivity_x_store.txt` & `single_part_emissivity_f_store.txt`:

The first file stores the log-10 values of the ratio of the photon frequency to the critical synchrotron frequency (see eqn. 6.17c in [1]) used in the definition of `F(x)` (eqn. 6.31c in [1]) for the single-particle synchrotron emissivity (eqn. 6.33 in [1]). The second file stores the log-10 values of `F(x)` (eqn. 6.31c in [1]) for the values store in `x`. 

The values in `x` are chosen such that the code will not need to manually calculate `F(x)`. Instead, we interpolate `F(x)` such that we can maually select the values from the array.

ii) 
