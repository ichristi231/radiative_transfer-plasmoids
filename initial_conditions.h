/************************ Insert initial parameters below ************************/
/*********************************************************************************/

/* Number of grid points in the electron energy range (in log space). */
#define kmax 350

/* Number of grid points in the photon frequency range (in log space). */
#define lmax 350

/* Magnetization of un-reconnected (upstream) plasma. Defined as the ratio of the
 magnetic energy density to particle eneregy density far upstream from the 
 reconnection layer. */
#define magnetization_sigma 10.0

/* Half-length of the reconnection layer in cm. */
#define half_length_L 5.e16

/* Magnetic field strength (in G) and magnetic energy desntiy within the radiating 
 blob. The magnetic field far upstream from the reconnecting plasma is set first.
 The magnetic field in the plasmoid is then sqrt(2) times larger (see Fig. 5 in 
 Sironi et al. 2016). */
#define B_up 1.5/3.3

/* Number of pairs to ions within the jet. */
#define pair_multiplicity_N_pm 1.

/* Bulk Lorentz factor and dimensionless velocity (i.e. normalized
 to the speed of light) of the jet. */
#define Gamma_j 12.

/* Power-law index of injected electron distribution. */
#define electron_injec_slope_p_e 2.1

/* Characteristic temperature in Kelvin of blackbody external 
 radiation field (if applicable). */
#define BLR_T 1.E4

/* Maximum Lorentz factor of injected electron distribution. */
#define gamma_e_max 5.E4

/* The angle (in degrees) between the blazar jet's access and the plasmoid's 
 direction of motion (i.e. with respect to the jet's co-moving frame). */
#define theta_prime 0.0

/* Fraction of the disk's luminosity which is reprocessed by the BLR. 
 Typicall value is 0.1 (i.e. 10%). */
#define frac_BLR 0.1

/*********************************************************************************/
/*********************************************************************************/