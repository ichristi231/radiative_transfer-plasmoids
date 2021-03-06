/* 
This is the leptonic version of a radiative transfer code developed by I. Christie 
(in collaboration with M. Petropoulou). If using this code, make reference to the
following:

i) Christie I. M., Petropoulou M., Sironi L., Giannios D., 2019, MNRAS, 482, 65 
(https://ui.adsabs.harvard.edu/abs/2019MNRAS.482...65C/abstract)

ii) Sironi L., Giannios D., Petropoulou M., 2016, MNRAS, 462, 48 
(https://ui.adsabs.harvard.edu/abs/2016MNRAS.462...48S/abstract)

This version uses the area-averaged quantities of plasmoids as obtained from 
particle-in-cell (PIC) simulations of relativistic magnetic reconnection. Given several 
initial free conditions (e.g. half length of the reconnection layer, magentic field upstream 
from the layer), it computes the co-moving spectra of a single plasmoid. The code will save the 
following to the 'results' directory: i) log10 of the photon frequency values (in Hz),
ii) log10 of the particle Lorentz factors, iii) log10 of the temporal evolution of the particle 
distribution, iv) log10 of the temporal evolution of the photon distribution, and v) a check 
that the particle number within the plasmoid is conserved.

Different processes can be switched on and off by use of the flags presented at the 
beginning of the main() function. 

Throughout this work, we have adopted the work of several studies which are referenced
below and which are referenced with appropriate eqn. # throughout the code.

i) Christie I. M., Petropoulou M., Sironi L., Giannios D., 2019, MNRAS, 482, 65 
(https://ui.adsabs.harvard.edu/abs/2019MNRAS.482...65C/abstract)

ii) Coppi & Blandford '90 (MNRAS, 245) 
(https://ui.adsabs.harvard.edu/abs/1990MNRAS.245..453C/abstract)

iii) Chiaberge & Ghisellini '99 (MNRAS, 306)
(http://adsabs.harvard.edu/abs/1999MNRAS.306..551C
iii) Ghisellini & Madau '96 (MNRAS, 280): http://adsabs.harvard.edu/abs/1996MNRAS.280...67G)

iv) Mastichiadis & Kirk '95 (A&A, 295) 
(http://adsabs.harvard.edu/abs/1995A%26A...295..613M)

v) Rybicki & Lightman '86
(http://adsabs.harvard.edu/abs/1986rpa..book.....R)

vi) Sironi et al. '16
(https://ui.adsabs.harvard.edu/abs/2016MNRAS.462...48S/abstract)

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "physical_constants.h"
#include "initial_conditions_leptonic.h"


/* Declaration of photon-photon pair production reaction 
   rate (see eqn. in 4.7 in Coppi & Blandford '90). */

double photon_annihilation_rate (double o);	

/* Linear interpolation scheme. */ 
double interpolation(double xx,double x1, double x2, double y1, double y2);

/* Radiative Transfer function for a single plasmoid. */
void single_plasmoid_calculation(double magnetization, double half_length, double B_upstream, 
  double plasmoid_avg_num_den, double pair_multiplicity, double Gamma_jet,
  double electron_injec_slope, double BLR_temp, double ge_max, double theta_p,
  double f_BLR, int plasmoid_number, int plasmoid_iteration);



int main()
{

  /* Set the number of plasmoids to perform the radiative transfer calculation. */
  int number_of_plasmoids = 1;

  /* Initialize arrays to store: i) the number of time-steps for each plasmoid
     in which the radiative section is called and ii) the area-averaged particle
     number density per plasmoids. The latter is determined directly from PIC (see
     Fig. 5 in Sironi at al. '16). */
  int plasmoid_iteration_number[number_of_plasmoids];
  double plasmoid_avg_num_den[number_of_plasmoids];


  /* Open and read in files of the number of time-steps and the area-averaged 
     particle number density per plasmoid. */
  FILE *pls_avg_num_den = fopen("plasmoid_properties/plasmoid_avg_particle_num_density.txt", "r");
  FILE *pls_iteration = fopen("plasmoid_properties/number_of_iterations_per_plasmoid.txt", "r");

  for (int i = 0; i < number_of_plasmoids; ++i)
  {
    fscanf(pls_iteration, "%d", &plasmoid_iteration_number[i]);
    fscanf(pls_avg_num_den, "%lf", &plasmoid_avg_num_den[i]);
  }

  fclose(pls_avg_num_den);
  fclose(pls_iteration);

  /* For each plasmoid, we call the radiative transfer function to calculate
     the co-moving photon and particle distributions for each plasmoid using
     the free parameters chosen in the initial_conditions.h file. */
  for (int i = 0; i < number_of_plasmoids; i++)
  {
    single_plasmoid_calculation(magnetization_sigma, half_length_L, B_up, plasmoid_avg_num_den[i],
    pair_multiplicity_N_pm, Gamma_j, electron_injec_slope_p_e, BLR_T, gamma_e_max, theta_prime,
    frac_BLR, i + 1, plasmoid_iteration_number[i]);
  }
  
  return 0;
}






/****** Radiative Transfer Calculation For Single Plasmoid ******/
/****************************************************************/
void single_plasmoid_calculation(double magnetization, double half_length, double B_upstream, 
  double plasmoid_avg_num_den, double pair_multiplicity, double Gamma_jet,
  double electron_injec_slope, double BLR_temp, double ge_max, double theta_p,
  double f_BLR, int plasmoid_number, int plasmoid_iteration) 
{
  printf("Beginning Plasmoid number %d. \n", plasmoid_number);

  double B, UB, Beta_jet;

  /* From the inputed magnetic field strength of the upstream plasma,
     we determine the magnetic field strength and energy density within
     the plasmoid. */
  B = sqrt(2.) * B_upstream; 
  UB = B * B / (8. * PI);

  /* From the inputed Lorentz factor of the jet, we compute the 
     dimensionless velocity, normalized to the speed of light. */
  Beta_jet = sqrt(1. - pow(Gamma_jet, -2.));


  /*********** Radiative Processes Switches **********/
  /***************************************************/
  /* Below, we create switches which permits to turn on/off a certain 
     radiative process. This can be used to test the code or explore 
     the emission from one or several different processes (named below). 
     If the switch is equal to unity, then the process is on, if zero, 
     it's off. */
  int sw_syn_e, sw_ict, sw_ickn, sw_ssa, sw_ggee, sw_ext, sw_rtc;

  /* Synchrotron emission */
  sw_syn_e = 1; 
  /* Photon annihilation to pair production */
  sw_ggee = 1; 
  /* Inverse Compton (Klein-Nishina regime) emission */
  sw_ickn = 1; 
  /* Inverse Compton (Thomson regime) emission */
  sw_ict = 1; 
  /* Synchrotron self-absorption */
  sw_ssa = 1; 
  /* The addition of radiation fields external to the jet. 
     Here, the radiation field is taken to be the BLR of 
     the jet, assumed to be a blackbody source with 
     temperature set in the initial_conditions.h . */
  sw_ext = 0; 
  /* The radiative transfer caculation. Used for checking
     the initalization of all quantities prior to the start
     of the radiative section. */
  sw_rtc = 1;
  /***************************************************/
  /***************************************************/



  /********** Importing Plasmoid Properties from PIC Simulations **********/
  /************************************************************************/
  /* Here, we import the plasmoid properties, derived from PIC results (Sironi
     et al. '16'), required to evolve the particle distribution within a 
     plasmoid and compute its photon spectrum. The quantities required are: 
     i) plasmoid's tranverse size (i.e. perpendicular to the plasmoid's 
     direction of motion) normalized to the half-length of the reconnection
     layer
     ii) plasmoid's Lorentz factor as measured in the reconnection frame.
     iii) plasmoid's co-moving volume 
     iv) the co-moving temporal derivative of the plasmoid's co-moving
     volume (normalized the  half-length of the reconnection layer cubed)

     Note that all quantities read in are, except the plasmoid's Lorentz
     factor, the log-10 values. 
  */
  int i = 0;


  char Lorentz_file_name[60], Volume_file_name[60];
  char Volume_derivative_file_name[60], size_file_name[60];

  snprintf(Lorentz_file_name, sizeof(Lorentz_file_name), 
    "plasmoid_properties/plasmoid_%d_Lorentz_factor.txt", plasmoid_number);
  snprintf(Volume_file_name, sizeof(Volume_file_name),
   "plasmoid_properties/plasmoid_%d_volume.txt", plasmoid_number);
  snprintf(Volume_derivative_file_name, sizeof(Volume_derivative_file_name),
   "plasmoid_properties/plasmoid_%d_volume_derivative.txt", plasmoid_number);
  snprintf(size_file_name, sizeof(size_file_name),
   "plasmoid_properties/plasmoid_%d_transverse_size.txt", plasmoid_number);

  FILE *plasmoid_Lorentz_file = fopen(Lorentz_file_name,"r");
  FILE *plasmoid_volume_file = fopen(Volume_file_name,"r");
  FILE *plasmoid_volume_derivative_file = fopen(Volume_derivative_file_name,"r");
  FILE *plamoid_size_file = fopen(size_file_name,"r");

  double plasmoid_Lorentz_factor[plasmoid_iteration];
  double plasmoid_volume[plasmoid_iteration];
  double plasmoid_volume_derivative[plasmoid_iteration];
  double plasmoid_size[plasmoid_iteration];



  for (i = 0; i < plasmoid_iteration; i++)
  {
    fscanf(plasmoid_Lorentz_file, "%lf", &plasmoid_Lorentz_factor[i]);
    fscanf(plasmoid_volume_file, "%lf", &plasmoid_volume[i]);
    fscanf(plasmoid_volume_derivative_file, "%lf", &plasmoid_volume_derivative[i]);
    fscanf(plamoid_size_file, "%lf", &plasmoid_size[i]);
  }

  fclose(plasmoid_Lorentz_file);
  fclose(plasmoid_volume_file);
  fclose(plasmoid_volume_derivative_file);
  fclose(plamoid_size_file);
  /************************************************************************/
  /************************************************************************/


  /************************ Import x & F(x) tables below ************************/
  /******************************************************************************/
  /* Below, we declare the tables of the x and F(x) values as governed by Rybicki & Lightman '86 
  eqn. 6.31c used for the single particle synchrotron emissivity. Note that the x values
  range (in log-10 space) from -20 to 2. Beyond these limits, F(x) evaluates roughly to zero.
  Note that both files are in log-10 space. */
  double x_store_table[400], f_store_table[400];

  FILE *xt; FILE *Ft;
  xt = fopen("required_tables/single_part_emissivity_x_store.txt", "r"); 
  Ft = fopen("required_tables/single_part_emissivity_f_store.txt", "r"); 
  for (i = 0; i < 400; i++) 
  {
    fscanf(xt, "%lf", &x_store_table[i]); 
    fscanf(Ft, "%lf", &f_store_table[i]);
  }
  fclose(xt); fclose(Ft);
  /******************************************************************************/
  /******************************************************************************/


  /************ Photon Frequency Array in Hz (& Dimensionless Frequency) ************/
  /**********************************************************************************/
  /* Below, we introduce the photon frequency array and the dimensionless photon 
     frequency, defined as the ratio of the photon energy to rest mass electron 
     energy. */
  FILE *nu_save = fopen("leptonic_version_results/log_10_photon_frequency_array.txt", "w");
  int l, ll;
  double nu[lmax], x[lmax], delta_x;

  for (l = 0; l < lmax; l++) 
  {
    nu[l] = pow(10, 7. + (l / (lmax - 1.)) * (40. - 7.)); 
    x[l] = planck_const * nu[l] / (electron_mass * c * c); 
    fprintf(nu_save, "%lf\n", log10(nu[l]));
  }
  fclose(nu_save);

  /* The difference of the log-space spacing between photon energies.
  This will be used when integrating over the photon energy 
  (e.g. inverse Compton scattering). */
  delta_x = log(x[1]) - log(x[0]);
  /**********************************************************************************/
  /**********************************************************************************/


  /*********** Electron Injection & Energy Distribution Parameters ***********/
  /***************************************************************************/
  FILE *gamma_e_save = fopen("leptonic_version_results/log_10_electron_Lorentz_factor_array.txt", "w");

  double ge_min, ge[kmax], ge_minus[kmax], ge_plus[kmax], delta_ge[kmax], delta_ge_int, Qe_inj[kmax];
  int k;

  /* Minimum Lorentz factor of injected electron distribution. Determine from 
     eqn. A4 in Christie et al. 2019. */ 
  ge_min = (electron_injec_slope - 2.) * (1. + 4. * magnetization * proton_mass /
      (plasmoid_avg_num_den * electron_mass * pair_multiplicity)) / (electron_injec_slope - 1.);
  printf("Minimum Lorentz factor: %f. Maximum Lorentz factor: %f. \n", ge_min, ge_max);

  /* Lorentz factor of the electron distribution, 
  evaluated at the half grid pints, and the spacing between 
  each of them as governed by Chiaberge & Ghisellini '99. */
  double max_ge_array = 1.E7;

  for (k = 1; k < kmax + 1; k++)
  {
  	ge[k - 1] = 1.01 * pow(max_ge_array / 1.01, (k - 1.) / (kmax - 1.));
    ge_minus[k - 1] = 1.01 * pow(max_ge_array / 1.01, (k - 0.5 - 1.) / (kmax - 1.));
    ge_plus[k - 1] = 1.01 * pow(max_ge_array / 1.01, (k + 0.5 - 1.) / (kmax - 1.));
    delta_ge[k - 1] = ge_plus[k - 1] - ge_minus[k - 1];
    fprintf(gamma_e_save, "%lf\n", log10(ge[k - 1]));
  }
  fclose(gamma_e_save);

  /* The difference of the log-space spacing between particle energies.
  This will be used when integrating over the particle Lorentz factor 
  (e.g. synchrotron emission). */
  delta_ge_int = log(ge[1]) - log(ge[0]);
  /***************************************************************************/
  /***************************************************************************/


  /************************ Inverse Compton (Thomson regime) Parameters ************************/
  /*********************************************************************************************/
  int ict_bound_index[lmax], ict_bound_loss_plus[kmax], ict_bound_loss_minus[kmax];
  int ict_Ne_gamma[lmax][lmax], ict_Ne_gamma_index[lmax][lmax];
  double ict_bound[lmax];
  double sum;

  /* Determine the upper bound of the integral in eqn. 44  of Mastichiadis & Kirk '95 for the inverse
  Compton scattering (Thomson regime) of photons off an electron distribution. For those upper bounds,
  we find the closest dimensionless photon frequency x-value within the array such that we can easily 
  perform the summation of the inverse Compton (Thomson regime) source term. */
  for (l = 0; l < lmax; l++)
  {
    if (3. * x[l] / 4. > 3. / (4. * x[l])) 
    {
      ict_bound[l] = 3. / (4. * x[l]);
    }
    else 
    {
      ict_bound[l] = 3. * x[l] / 4.;
    }

    sum = 0;
    for (ll = 0; ll < lmax; ll++) 
    {
      sum += (ict_bound[l] >= x[ll]);
    }
    ict_bound_index[l] = sum;
  }

  /* Below, we determine the index of the dimensionless photon frequency array which is closest
     to the value of 3/4*gamma, i.e. the upper bound of the photon energy density in eqn. 43 of Mastichiadis
     & Kirk '95. We find the values, evaluated at the half-mesh points in the electron energy density. These
     indices will then be used as the limit in which the summation (integration) is carried to. */
  for (k = 0; k < kmax; k++)
  {
  	ict_bound_loss_plus[k] = 0; 
    ict_bound_loss_minus[k] = 0;

  	for (l = 0; l < lmax; l++) 
    {
      ict_bound_loss_plus[k] += (3. / (4. * ge_plus[k])) >= x[l]; 
      ict_bound_loss_minus[k] += (3. / (4. * ge_minus[k])) >= x[l];
    }
    
  }

  /* Below, we determine for which values of the dimensionless photon frequency fall within the range 
  of the minimum (i.e. 1) & maximum Lorentz factor of the electron distribution. This is used
  for evaluating the electron distribution as goverened by eqn. 44 in Mastichiadis & Kirk '95.*/
  for (l = 0; l < lmax; l++)
  {
  	for (ll = 0; ll < lmax; ll++)
  	{
  		sum = 0;
  		ict_Ne_gamma[l][ll] = (sqrt(3. * x[l] / (4. * x[ll])) >= 1.) && (sqrt(3. * x[l] / (4. * x[ll])) <= ge_max);
  		if (ict_Ne_gamma[l][ll] == 1)
  		{
  			for (k = 0; k < kmax; k++) 
        {
          sum += (sqrt(3. * x[l] / (4. * x[ll])) >= ge[k]);
        }
  			ict_Ne_gamma_index[l][ll] = sum;
  		}
  	}
  }
  /*********************************************************************************************/
  /*********************************************************************************************/


  /************************ Inverse Compton (Klein Nishina regime) Parameters ************************/
  /***************************************************************************************************/
  int ickn_bound_loss[kmax], ickn_bound_source[lmax], ickn_gamma_source[lmax];

  /* Below, we determine the appropriate index in the dimensionless photon frequency array in order to 
  determine the lower bound on the Klein-Nishina loss term integral defined in eqn. 45 of Mastichiadis 
  & Kirk '95. */
  for (k = 0; k < kmax; k++)
  {
  	ickn_bound_loss[k] = 0;
  	for (l = 0; l < lmax; l++) 
    {
      ickn_bound_loss[k] += (3. / (4. * ge[k])) >= x[l];
    }
  }

  /* Below, we determine the index in the dimensionless photon frequency array in order to determine
  the lower bound on the Klein-Nishina photon source term defined in eqn. 46 of Mastichiadis
  & Kirk '95. */
  for (l = 0; l < lmax; l++)
  {
  	ickn_bound_source[l] = 0;
  	for (ll = 0; ll < lmax; ll++) 
    {
      ickn_bound_source[l] += (3. / (4. * x[l])) >= x[ll];
    }
  }

  /* Below, we determine the index within the electron Lorentz factor array which is closest to
  the dimensionless frequency value. This will then be evaluated within the electron distribution
  function as governed by eqn. 46 in Mastichiadis & Kirk '95. */
  for (l = 0; l < lmax; l++)
  {
  	ickn_gamma_source[l] = 0;
  	if (x[l] >= 1. && x[l] <= ge_max) 
    {
      for (k = 0; k < kmax; k++) 
      {
        ickn_gamma_source[l] += x[l] >= ge[k];
      }
    }
  }
  /***************************************************************************************************/
  /***************************************************************************************************/


  /****************** Single particle electron Synchrotron emissivity in cgs. ******************/
  /*********************************************************************************************/
  double const_syn_e, x_store, f_store, pe_single[lmax][kmax];
  int x_index;

  /* Unitless constant which belongs in front of the electron synchrotron cooling rate. */
  const_syn_e = 4. * sigmaT * UB / (3. * electron_mass *c);

  /* Below, we determine the electron single particle synchrotron emissivity. In doing so,
     we first compute the value of x (i.e. ratio of photon frequency to the critical 
     synchrotron cooling frequency, see Rybicki & Lightman '86 just following eqn. 6.31c) using 
     our values of the photon frequency nu and the electron Lorentz factor gamma. We then, for 
     every value of x which falls in the range of the imported x-table (see above), we compute 
     the synchroton emissivity using eqn. 6.33 in Rybicki & Lightman '86. */
  for (l = 0; l < lmax; l++)
  {
    f_store = 0;
    for (k = 0; k < kmax; k++)
    {
    	x_store = nu[l] * 4. * PI * electron_mass * c /(3. * electric_charge * B * ge[k] * ge[k]);

    	if (x_store >= 1E-20 && x_store <= 100)
    	{
    		x_index = 0;
    		for (i = 0; i < 400; i++) 
        {
          x_index += log10(x_store) >= x_store_table[i];
        }
    		f_store = pow(10., f_store_table[x_index]);
    	}
    	else 
      {
        f_store = 0;
      }

    	pe_single[l][k] = sqrt(3.) * pow(electric_charge, 3.) * B * f_store / (electron_mass * c * c) ;
    }
  }
  /*********************************************************************************************/
  /*********************************************************************************************/

  /************************ Photon-Photon Pair Production Parameters ***************************/
  /*********************************************************************************************/
  int gg_ee_gamma_bound_source[kmax], gg_ee_bound_loss[lmax], gg_ee_bound_source[kmax];

  /* Below, we determine the index in the dimensionless frequency array in which two times the 
  electron Lorentz factor is closest to the dimensionless frequency. This is used to evaluate the
  photon distribution at two times the electron Lorentz factor as stated in eqn. 57 in 
  Mastichiadis & Kirk '95. */
  for (k = 0; k < kmax; k++)
  {
    gg_ee_gamma_bound_source[k] = 0;
    for (l = 0; l < lmax; l++) 
    {
      gg_ee_gamma_bound_source[k] += 2. * ge[k] >= x[l];
    }
  }

  /* Below, we determine the index in the dimensionless photon frequency array in which 1/x = x, where
  x is the dimensionless photon frequency. This index is used for the lower bound of the loss term as 
  governed by eqn. 54 in Mastichiadis & Kirk '95. */
  for (l = 0; l < lmax; l++)
  {
    gg_ee_bound_loss[l] = 0;
    for (ll = 0; ll < lmax; ll++) 
    {
      gg_ee_bound_loss[l] += (1. / x[l]) >= x[ll];
    }
  }

  /* Below, we determine the index in the dimensionless photon frequency array in which the inverse of 
  two times the electron Lorentz factor is closest to the dimensionless photon frequency. This is used 
  in the lower bound of the source term integral as governed by eqn. 57 in Mastichiadis & Kirk '95. */
  for (k = 0; k < kmax; k++)
  {
    gg_ee_bound_source[k] = 0;
    for (l = 0; l < lmax; l++) 
    {
      gg_ee_bound_source[k] += (1. / (2. * ge[k])) >= x[l];
    }
  }
  /*********************************************************************************************/
  /*********************************************************************************************/


  /******************* Radiative Transfer Section *****************/
  /****************************************************************/
  /* Below, we declare all variables, particle species distributions, and all 
     source/loss terms used in the radiative transfer section. Additionally, we
     create txt files used in saving the temporal evolution of all particle
     species distributions. */

  char nuLnu_file_name[150], Ne_file_name[150], particle_conservation_file_name[150];

  snprintf(nuLnu_file_name, sizeof(nuLnu_file_name), 
    "leptonic_version_results/co_moving_nu_l_nu_plasmoid_%d.txt", plasmoid_number);
  snprintf(Ne_file_name, sizeof(Ne_file_name), 
    "leptonic_version_results/electron_distribution_plasmoid_%d.txt", plasmoid_number);
  snprintf(particle_conservation_file_name, sizeof(particle_conservation_file_name), 
    "leptonic_version_results/particle_conservation_check_plasmoid_%d.txt", plasmoid_number);

  /* Opens files in which to save the log-10 values of the co-moving
     particle and photon distributions. */
  FILE *nu_L_nu_save = fopen(nuLnu_file_name, "w");
  FILE *N_e_save = fopen(Ne_file_name, "w");

  /* Opens file to save the log-10 values of: i) the cumulative number 
     of injected particles and ii) the total number of particles 
     determined from the updated electron distribution. This allows
     to check particle number conservation. */
  FILE *particle_conservation = fopen(particle_conservation_file_name, "w");

  /* Declares the particle and photon distributions. */
  double N_e[kmax], N_x[lmax];

  /* Declares all source terms Q. */
  double Q_syn_e[lmax], Q_ICT_e[lmax], Q_ICKN_e[lmax],  Q_gg_ee[kmax];

  /* Declares all loss terms L. */
  double L_ICKN_e[kmax], L_ICT_e_plus[kmax], L_ICT_e_minus[kmax], L_gg_ee[lmax], L_ssa[lmax];

  /* Declares the v2 and v3 coefficients needed to solve the tri-
     diagonal matrix algorithm for updating the particle distribution
     (see eqn. 10 in Chiaberge & Ghisellini '99). */
  double v2_e[kmax], v3_e[kmax];

  /* Declartion of the photon escape timescale and plasmoid volume. 
     These quantities are updated each timestep as the plasmoid size
     changes. */
  double t_esc, vol;

  /* The value of delta_t is determined from the smallest value of the 
     plasmoid's size (i.e. its initial size when born in the layer). */
  double delta_t = pow(10., plasmoid_size[0]) * half_length / (2. * c);

  /* Declaration of quantities needed to compute the number of particles
     after each time-step. At the end of each timestep, we compare the 
     cumulative injected particles to the number of particles within the 
     electron distribution. */
  double N_tot_inj = 0;
  double N_tot = 0;
  double sum_inj = 0;

  /* Declaration of the plasmoid's Lorentz factor as measured in the
     rest frame of the SMBH. This quantity is used to boost the energy
     densities of photon fields external to the jet and is updated each 
     time step. */
  double plasmoid_Lorentz_factor_SMBH = 0;

  /* Declaring the photon distribution external radiation fields. */
  double N_external[lmax];



  /******* Initialization of all species distributions at initial timestep. *******/
  /********************************************************************************/
  for (k = 0; k < kmax; k++)
  {
    N_e[k] = 0;
  } 
  for (l = 0; l < lmax; l++) 
  {
    N_x[l] = 0;
    N_external[l] = 0;
  }
  /********************************************************************************/
  /********************************************************************************/


  /* Begin time-loop for radiative calulations of plasmoid. */
  for (i = 0; i < (plasmoid_iteration - 1) * sw_rtc; i++)
  {
    /************* Update plasmoid quantities *************/
    /******************************************************/
    /*  For each time-step, we re-compute the escape time, the particle
        injection rate, and the plasmoid volume as the plasmoid's size changes. 
    */
    for (k = 0; k < kmax; k++) 
    {
      /* Particle injection rate (see eqn. 5 in Christie et al. '19). */
      Qe_inj[k] = (ge[k] >= ge_min) * (ge[k] <= ge_max) * (1. - electron_injec_slope) * pow(ge[k], -electron_injec_slope) / 
      (pow(ge_max, 1. - electron_injec_slope) - pow(ge_min, 1. - electron_injec_slope)) * (10. * 
        plasmoid_avg_num_den / 4.) * (pow(B_upstream, 2.) / (4. * PI * proton_mass * c * c * 
          magnetization)) * c * plasmoid_Lorentz_factor[i] * pow(10, plasmoid_volume_derivative[i]) * 
            half_length * half_length;
    }
    /* Photon escape timescale, in seconds. */
    t_esc = pow(10., plasmoid_size[i]) * half_length / (2. * c);

    /* Plasmoid volume in cm^3. */
    vol = pow(10., plasmoid_volume[i]) * pow(half_length, 3.);

    /* Plasmoid Lorentz factor with respect to the BH frame. */
    plasmoid_Lorentz_factor_SMBH = Gamma_jet * plasmoid_Lorentz_factor[i] * (1. + 
      Beta_jet * sqrt(1. - pow(plasmoid_Lorentz_factor[i], -2.)) * cos(theta_p * PI / 180.));

    for (l = 0; l < lmax; l++)
    {
      /* Updating the external photon field from the BLR. */
      N_external[l] = 0.265 * f_BLR * (15. / pow(boltzman_const * BLR_temp * plasmoid_Lorentz_factor_SMBH * PI, 4.)) *
      x[l] * x[l] * vol * pow(electron_mass * c * c, 3.) / (exp(x[l] * electron_mass * c * c / 
        (boltzman_const * BLR_temp * plasmoid_Lorentz_factor_SMBH)) - 1.) * pow(plasmoid_Lorentz_factor_SMBH, 2.);

      if (N_external[l] < 1.E-200) 
      {
        N_external[l] = 1.E-200;
      }
    }
    /******************************************************/
    /******************************************************/



  	/************* Synchrotron Source Terms (Electron) *************/
  	/***************************************************************/
    for (l = 0; l < lmax; l++)
  	{
  		Q_syn_e[l] = 0;
  	  for (k = 0; k < kmax; k++) 
      {
        Q_syn_e[l] += pe_single[l][k] * ge[k] * N_e[k] * delta_ge_int / (planck_const * x[l]);
      }
  	}
  	/***************************************************************/
  	/***************************************************************/



  	/************* Electron Inverse Compton (Thomson Regime) Source Term *************/
  	/*********************************************************************************/
  	/* As determined from eqn. 44 in Mastichiadis & Kirk '95. */
  	for (l = 0; l < lmax; l++)
  	{
  		Q_ICT_e[l] = 0;
  		if (ict_bound[l] >= x[0])
  		{
  			for (ll = 0; ll < lmax; ll++)
  			{
  				if (ict_Ne_gamma[l][ll] == 1) 
          {
            Q_ICT_e[l] += c * sigmaT * sqrt(3.0 / x[l]) / (4. * vol) * delta_x * sqrt(x[ll]) * 
            (N_x[ll] + N_external[ll] * sw_ext) * N_e[ict_Ne_gamma_index[l][ll]];
          }
  			}
      }
    }
  	/*********************************************************************************/
  	/*********************************************************************************/



  	/************* Electron Inverse Compton (Thomson Regime) Loss Terms *************/
  	/********************************************************************************/
  	/* As determined from eqns. 42 & 43 in Mastichiadis & Kirk '95. However, we note that the
  	 photon energy density, determined from eqn. 43, is dependent upon the electron Lorentz factor
  	 as is brought within the derivative in eqn. 42. As such, we are required to evaluate it at
  	 the half grid points in the electron Lorentz factor. */
  	for (k = 0; k < kmax; k++)
  	{
  		L_ICT_e_plus[k] = 0; L_ICT_e_minus[k] = 0;

  		//for (l = 0; l < ict_bound_loss_minus[k]; l++)
      for (l = 0; l < lmax; l++) 
      {
        L_ICT_e_plus[k] += delta_x * x[l] * x[l] * (N_x[l] + N_external[l] * sw_ext) * 4. * sigmaT * c / (3. * vol);
      } 
  		//for (l = 0; l < ict_bound_loss_minus[k]; l++) 
      for (l = 0; l < lmax; l++) 
      {
        L_ICT_e_minus[k] += delta_x * x[l] * x[l] * (N_x[l] + N_external[l] * sw_ext) * 4. * sigmaT * c / (3. * vol);
      } 
  	}
  	/********************************************************************************/
  	/********************************************************************************/



  	/*********** Electron Inverse Comptron (Klein-Nishina Regime) Loss Term ***********/
  	/**********************************************************************************/
  	/* Determined from eqn. 45 in Mastichiadis & Kirk '95. */
  	for (k = 0; k < kmax; k++)
  	{
  		L_ICKN_e[k] = 0;
  		for (l = ickn_bound_loss[k]; l < lmax; l++) 
      {
        L_ICKN_e[k] += (sigmaT * c / vol) * (N_e[k] / ge[k]) * delta_x * (N_x[l] + N_external[l] * sw_ext);
      }
  	}
  	/**********************************************************************************/
  	/**********************************************************************************/



  	/*********** Electron Inverse Comptron (Klein-Nishina Regime) Source Term ***********/
  	/************************************************************************************/
  	/* Determined from eqn. 46 in Mastichiadis & Kirk '95. */
  	for (l = 0; l < lmax; l++)
  	{
  		Q_ICKN_e[l] = 0;
  		if (ickn_gamma_source[l] != 0 && ickn_bound_source[l] != 0) 
      {
        for (ll = ickn_bound_source[l]; ll < lmax; ll++) 
        {
          Q_ICKN_e[l] += (sigmaT * c / vol) * (N_e[ickn_gamma_source[l]] / x[l]) * delta_x * (N_x[ll] 
            + N_external[ll] * sw_ext);
        }
      }
  	}
  	/************************************************************************************/
  	/************************************************************************************/



  	/*********** Electron Distribution Update Coefficients ***********/
  	/*****************************************************************/
  	/* Coefficients used in updating the electron particle distributions
  	as determined from Chiaberge & Gisellini '99, eqn. 10. */
  	for (k = 0; k < kmax; k++)
  	{
      v2_e[k] = 1. + (delta_t / delta_ge[k]) * (ge_minus[k] * ge_minus[k]) * (const_syn_e * sw_syn_e + 
        L_ICT_e_minus[k] * sw_ict);

      v3_e[k] = -(delta_t / delta_ge[k]) * (ge_plus[k] * ge_plus[k]) * (const_syn_e * sw_syn_e + 
        L_ICT_e_plus[k] * sw_ict);
  	}
  	/*****************************************************************/
  	/*****************************************************************/



  	/************* Photon Distribution Update *************/
  	/******************************************************/
  	for (l = 0; l < lmax; l++)
  	{
      N_x[l] = N_x[l] * (1. - delta_t / t_esc) + delta_t * (Q_syn_e[l] * sw_syn_e + Q_ICT_e[l] * 
        sw_ict + Q_ICKN_e[l] * sw_ickn);
  		
      if (N_x[l] < 0) 
      {
        N_x[l] = 1.E-200;
      }
  	}
  	/******************************************************/
  	/******************************************************/

  	

  	/*********** Photon-Photon Pair Production Loss Term ***********/
  	/***************************************************************/
  	/* Determined from eqn. 54 in Mastichiadis & Kirk '95. However, we
  	integrate over the dimensionless frequency range of the soft photons 
  	from 1/x to Max[x], where x is the dimensionless frequency. */
  	for (l = 0; l < lmax; l++)
  	{
  		L_gg_ee[l] = 0;
  		if (gg_ee_bound_loss[l] != lmax)
  		{
  			for (ll = gg_ee_bound_loss[l]; ll < lmax; ll++) 
        {
          L_gg_ee[l] += (1. / vol) * (N_x[l] + N_external[l] * sw_ext) * delta_x
          * x[ll] * (N_x[ll] + N_external[ll] * sw_ext) * photon_annihilation_rate(x[l]
            * x[ll]);
        }
  		}
  	}
  	/***************************************************************/
  	/***************************************************************/



  	/*********** Photon-Photon Pair Production Source Term ***********/
  	/*****************************************************************/
  	/* Determined from eqn. 57 in Mastichiadis & Kirk '95. However, we
  	integrate over the dimensionless frequency range of the soft photons 
  	from 1/(2*gamma) to Max[x], where x is the dimensionless frequency and 
  	gamma is the Lorentz factor of the produced electron/positron. */
  	for (k = 0; k < kmax; k++)
  	{
  		Q_gg_ee[k] = 0;
  		for (l = gg_ee_bound_source[k]; l < lmax; l++)
  		{
        Q_gg_ee[k] += (4. / vol) * (N_x[gg_ee_gamma_bound_source[k]] + N_external[gg_ee_gamma_bound_source[k]] * sw_ext)
        * delta_x * x[l] * photon_annihilation_rate(2. * ge[k] * x[l]) * (N_x[l] + N_external[l] * sw_ext);
  		}
  	}
  	/*****************************************************************/
  	/*****************************************************************/



  	/************* Photon Distribution Update *************/
  	/******************************************************/
  	for (l = 0; l < lmax; l++)
  	{
  		N_x[l] -= delta_t * L_gg_ee[l] * (L_gg_ee[l] >= 0.) * sw_ggee;
  		if (N_x[l] < 0) 
      {
        N_x[l] = 1.E-200;
      }
  	}
  	/******************************************************/
  	/******************************************************/



  	/************* Electron Distribution Update *************/
  	/********************************************************/
  	for (k = kmax - 1; k >= 0; k--)
  	{
  		if (k == kmax - 1) 
      {
        N_e[k] = (N_e[k] + delta_t * (Qe_inj[k] - L_ICKN_e[k] * sw_ickn + Q_gg_ee[k] * sw_ggee)) / v2_e[k];
      }
  		if (k != kmax - 1) 
      {
        N_e[k] = (N_e[k] + delta_t * (Qe_inj[k] - L_ICKN_e[k] * sw_ickn + Q_gg_ee[k] * sw_ggee) - 
          v3_e[k] * N_e[k + 1]) / v2_e[k];
      }
  		if (N_e[k] <= 0) 
      {
        N_e[k] = 1.E-20;
      }
  	}

  	for (k = 0; k < kmax; k++)
  	{
  		if (k != kmax - 1) 
      {
        fprintf(N_e_save, "%lf ", log10(pow(ge[k],electron_injec_slope) * N_e[k]));
      }
  		else 
      {
        fprintf(N_e_save, "%lf\n", log10(pow(ge[k],electron_injec_slope) * N_e[k]));
      }
  	}
  	/********************************************************/
  	/********************************************************/



  	/************* Synchrotron Self Absorption *************/
  	/*******************************************************/
  	/* Here, we follow the prescription for the loss term due to synchrotron
  	self absorption as specified in eqn. 40 in Mastichiadis & Kirk '95 however,
  	we adopt the full expression of the absoprtion coefficient as governed by 
  	eqn. 6.50 in Rybicki & Lightman. */
  	for (l = 0; l < lmax; l++)
  	{
  		L_ssa[l] = 0.;
  		for (k = 0; k < kmax - 1; k++) 
      {
        L_ssa[l] += -(planck_const * planck_const / (8. * PI * electron_mass * pow(electron_mass * 
          c * c, 2.))) * (pow(x[l], -2.)) * (c * (N_x[l] + N_external[l] * sw_ext) / vol) 
        * delta_ge_int * pe_single[l][k] * ((N_e[k + 1] - N_e[k]) / (delta_ge_int) - 2. * N_e[k]);
      }
  		if (L_ssa[l] <= 0.) 
      {
        L_ssa[l] = 1.E-20;
      }
  	}
  	/*******************************************************/
  	/*******************************************************/



  	/************* Photon Distribution Update *************/
  	/******************************************************/
  	/* We update the photon distribution to remove those photons 
  	lossed due to synchrotron self absorption. */
  	for (l = 0; l < lmax; l++)
  	{
  		N_x[l] -= delta_t * L_ssa[l] * (L_ssa[l] >= 0.) * sw_ssa;

  		if (N_x[l] < 0) 
      {
        N_x[l] = 1.E-200;
      }

  		if (l != lmax - 1) 
      {
        fprintf(nu_L_nu_save, "%lf ", log10(N_x[l] * pow(planck_const * nu[l], 2.) / 
          (electron_mass * c *c * t_esc)));
      }
  		else 
      {
        fprintf(nu_L_nu_save, "%lf\n", log10(N_x[l]*pow(planck_const * nu[l], 2.) / 
          (electron_mass * c * c * t_esc)));
      }
  	}
  	/******************************************************/
  	/******************************************************/



    /******************* Checks Particle Number Conservation **************/
    /**********************************************************************/
    sum_inj = 0; N_tot = 0;

    for (k = 0; k < kmax; k++)
    {
      N_tot += ge[k] * N_e[k] * delta_ge_int;
      sum_inj += ge[k] * Qe_inj[k] * delta_ge_int;
    }
    N_tot_inj += delta_t * sum_inj;

    fprintf(particle_conservation, "%lf %lf\n", log10(N_tot_inj), log10(N_tot));
    /**********************************************************************/
    /**********************************************************************/



    /* Prints the interation value every 100 interations. */
    if(i % 100 == 0) 
    {
      printf("Finished %d out of %d.\n", i, plasmoid_iteration);
    }
  }

  fclose(nu_L_nu_save); 
  fclose(N_e_save); 
  fclose(particle_conservation);
  /************************************************************************/
  /******************** End Radiative Transfer Section ********************/

  printf("Ending Plasmoid number %d. \n \n", plasmoid_number);
}
/***************************************************************/
/***************************************************************/



/*********** Reaction Rate of Photon-Photon Pair Production ***********/
/**********************************************************************/
/* See eqn. 4.7 in Coppi & Blandford '90. */
double photon_annihilation_rate(double o)
{
  if (o >= 1.)
  {
    return 0.652 * sigmaT * c * log(o) * (o * o - 1.) / pow(o, 3.);
  }
  else 
  {
    return 0.;
  }
}
/**********************************************************************/
/**********************************************************************/



/********************* Linear Interpolation Scheme ********************/
/**********************************************************************/
double interpolation(double xx,double x1, double x2, double y1, double y2)
{
  return y1 + (xx - x1) * (y2 - y1) / (x2 - x1); 
}
/**********************************************************************/
/**********************************************************************/
