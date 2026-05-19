/** @file background.c Documented background module
 *
 * * Julien Lesgourgues, 17.04.2011
 * * routines related to ncdm written by T. Tram in 2011
 *
 * Deals with the cosmological background evolution.
 * This module has two purposes:
 *
 * - at the beginning, to initialize the background, i.e. to integrate
 *    the background equations, and store all background quantities
 *    as a function of conformal time inside an interpolation table.
 *
 * - to provide routines which allow other modules to evaluate any
 *    background quantity for a given value of the conformal time (by
 *    interpolating within the interpolation table), or to find the
 *    correspondence between redshift and conformal time.
 *
 *
 * The overall logic in this module is the following:
 *
 * 1. most background parameters that we will call {A}
 * (e.g. rho_gamma, ..) can be expressed as simple analytical
 * functions of a few variables that we will call {B} (in simplest
 * models, of the scale factor 'a'; in extended cosmologies, of 'a'
 * plus e.g. (phi, phidot) for quintessence, or some temperature for
 * exotic particles, etc...).
 *
 * 2. in turn, quantities {B} can be found as a function of conformal
 * time by integrating the background equations.
 *
 * 3. some other quantities that we will call {C} (like e.g. the
 * sound horizon or proper time) also require an integration with
 * respect to time, that cannot be inferred analytically from
 * parameters {B}.
 *
 * So, we define the following routines:
 *
 * - background_functions() returns all background
 *    quantities {A} as a function of quantities {B}.
 *
 * - background_solve() integrates the quantities {B} and {C} with
 *    respect to conformal time; this integration requires many calls
 *    to background_functions().
 *
 * - the result is stored in the form of a big table in the background
 *    structure. There is one column for conformal time 'tau'; one or
 *    more for quantities {B}; then several columns for quantities {A}
 *    and {C}.
 *
 * Later in the code, if we know the variables {B} and need some
 * quantity {A}, the quickest and most precise way is to call directly
 * background_functions() (for instance, in simple models, if we want
 * H at a given value of the scale factor). If we know 'tau' and want
 * any other quantity, we can call background_at_tau(), which
 * interpolates in the table and returns all values. Finally it can be
 * useful to get 'tau' for a given redshift 'z': this can be done with
 * background_tau_of_z(). So if we are somewhere in the code, knowing
 * z and willing to get background quantities, we should call first
 * background_tau_of_z() and then background_at_tau().
 *
 *
 * In order to save time, background_at_tau() can be called in three
 * modes: short_info, normal_info, long_info (returning only essential
 * quantities, or useful quantities, or rarely useful
 * quantities). Each line in the interpolation table is a vector whose
 * first few elements correspond to the short_info format; a larger
 * fraction contribute to the normal format; and the full vector
 * corresponds to the long format. The guideline is that short_info
 * returns only geometric quantities like a, H, H'; normal format
 * returns quantities strictly needed at each step in the integration
 * of perturbations; long_info returns quantities needed only
 * occasionally.
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# background_init() at the beginning
 * -# background_at_tau(), background_tau_of_z() at any later time
 * -# background_free() at the end, when no more calls to the previous functions are needed
 */

#include "background.h"

/**
 * Background quantities at given conformal time tau.
 *
 * Evaluates all background quantities at a given value of
 * conformal time by reading the pre-computed table and interpolating.
 *
 * @param pba           Input: pointer to background structure (containing pre-computed table)
 * @param tau           Input: value of conformal time
 * @param return_format Input: format of output vector (short, normal, long)
 * @param intermode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

int background_at_tau(
                      struct background *pba,
                      double tau,
                      short return_format,
                      short intermode,
                      int * last_index,
                      double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                      ) {

  /** Summary: */

  /** - define local variables */

  /* size of output vector, controlled by input parameter return_format */
  int pvecback_size;

  /** - check that tau is in the pre-computed range */

  class_test(tau < pba->tau_table[0],
             pba->error_message,
             "out of range: tau=%e < tau_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->tau_table[0]);

  class_test(tau > pba->tau_table[pba->bt_size-1],
             pba->error_message,
             "out of range: tau=%e > tau_max=%e\n",tau,pba->tau_table[pba->bt_size-1]);

  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }

  /** - interpolate from pre-computed table with array_interpolate()
      or array_interpolate_growing_closeby() (depending on
      interpolation mode) */

  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->tau_table,
                                        pba->bt_size,
                                        pba->background_table,
                                        pba->d2background_dtau2_table,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->tau_table,
                                                        pba->bt_size,
                                                        pba->background_table,
                                                        pba->d2background_dtau2_table,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }

  return _SUCCESS_;
}

/**
 * Conformal time at given redshift.
 *
 * Returns tau(z) by interpolation from pre-computed table.
 *
 * @param pba Input: pointer to background structure
 * @param z   Input: redshift
 * @param tau Output: conformal time
 * @return the error status
 */

int background_tau_of_z(
                        struct background *pba,
                        double z,
                        double * tau
                        ) {

  /** Summary: */

  /** - define local variables */

  /* necessary for calling array_interpolate(), but never used */
  int last_index;

  /** - check that \f$ z \f$ is in the pre-computed range */
  class_test(z < pba->z_table[pba->bt_size-1],
             pba->error_message,
             "out of range: z=%e < z_min=%e\n",z,pba->z_table[pba->bt_size-1]);

  class_test(z > pba->z_table[0],
             pba->error_message,
             "out of range: a=%e > a_max=%e\n",z,pba->z_table[0]);

  /** - interpolate from pre-computed table with array_interpolate() */
  class_call(array_interpolate_spline(
                                      pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      pba->d2tau_dz2_table,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;
}

/**
 * Background quantities at given \f$ a \f$.
 *
 * Function evaluating all background quantities which can be computed
 * analytically as a function of {B} parameters such as the scale factor 'a'
 * (see discussion at the beginning of this file). In extended
 * cosmological models, the pvecback_B vector contains other input parameters than
 * just 'a', e.g. (phi, phidot) for quintessence, some temperature of
 * exotic relics, etc...
 *
 * @param pba           Input: pointer to background structure
 * @param pvecback_B    Input: vector containing all {B} type quantities (scale factor, ...)
 * @param return_format Input: format of output vector
 * @param pvecback      Output: vector of background quantities (assumed to be already allocated)
 * @return the error status
 */

int background_functions(
                         struct background *pba,
                         double * pvecback_B, /* Vector containing all {B} quantities. */
                         short return_format,
                         double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                         ) {

  /** Summary: */

  /** - define local variables */

  /* total density */
  double rho_tot;
  /* total pressure */
  double p_tot;
  /* total relativistic density */
  double rho_r;
  /* total non-relativistic density */
  double rho_m;
 /* total dark sector density */
  double rho_ds;
  /* scale factor relative to scale factor today */
  double a_rel;
  /* background ncdm quantities */
  double rho_ncdm,p_ncdm,pseudo_p_ncdm;
  /* index for n_ncdm species */
  int n_ncdm;
  /* scale factor */
  double a;
  /* scalar field quantities */
  double phi, phi_prime;
  

  /** - initialize local variables */
  a = pvecback_B[pba->index_bi_a];
  rho_tot = 0.;
  p_tot = 0.;
  rho_r=0.;
  rho_m=0.;
  rho_ds = 0.;

  a_rel = a / pba->a_today;

  class_test(a_rel <= 0.,
             pba->error_message,
             "a = %e instead of strictly positive",a_rel);

  /** - pass value of \f$ a\f$ to output */
  pvecback[pba->index_bg_a] = a;

  /** - compute each component's density and pressure */

  /* photons */
  pvecback[pba->index_bg_rho_g] = pba->Omega0_g * pow(pba->H0,2) / pow(a_rel,4);
  rho_tot += pvecback[pba->index_bg_rho_g];
  p_tot += (1./3.) * pvecback[pba->index_bg_rho_g];
  rho_r += pvecback[pba->index_bg_rho_g];

  /* baryons */
  pvecback[pba->index_bg_rho_b] = pba->Omega0_b * pow(pba->H0,2) / pow(a_rel,3);
  rho_tot += pvecback[pba->index_bg_rho_b];
  p_tot += 0;
  rho_m += pvecback[pba->index_bg_rho_b];

  /* cdm */
  if (pba->has_cdm == _TRUE_) {
    pvecback[pba->index_bg_rho_cdm] = pba->Omega0_cdm * pow(pba->H0,2) / pow(a_rel,3);
    rho_tot += pvecback[pba->index_bg_rho_cdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_cdm];
  }


/* ds */
  if (pba->has_ds == _TRUE_) {
	double H0 = pba->H0;
	double HSb,fR0,om0,Om0,Ol,Or0,w_fR,rho_fR,b,l,N,PI,G;
	//double a = a_rel;  //PUÒ ESSERE QUI IL PROBLEMA??
	Om0 = pba->Omega0_b + pba->Omega0_cdm;
	om0 = pba->Omega0_b + pba->Omega0_cdm;
	Or0 = pba->Omega0_ur + pba->Omega0_g;
	Ol  = pba->Omega0_lambda;
	l   = pba->l;
	fR0 = 1 - pba->fR0; 
  	if (l == 0){		//lcdm
  		static int printed = 0;
  		if (!printed) {
      		printf("DEBUG: BACKGROUND MODEL = LCDM (Case 0) SELECTED\n");
     		printed = 1;          // non stampare più
   		}
   		
		pvecback[pba->index_bg_rho_ds] = rho_fR * pba->Omega0_ds * pow(pba->H0,2);
       		pvecback[pba->index_bg_p_ds]   = w_fR *pvecback[pba->index_bg_rho_ds];
  		
    		}
  	if (l > 0){
    		if (l == 1){
    			static int printed = 0;
    			if (!printed) {
    			printf("-------------------------------------------\n");
      			printf("DEBUG: f(R) model = HS1 (Case 1) SELECTED\n");
      			printf("-------------------------------------------\n");
      			fflush(stdout);	
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = (2 - 2*fR0)*pow(-3.0/2.0*Om0 - 2*Or0 + 2, 2)/pow(Ol, 2);

   			w_fR = (2.0/3.0)*pow(Ol, 2)*pow(a, 4)*pow(b, 2)*(122880*pow(Ol, 6)*Om0*pow(a, 20) + 49152*pow(Ol, 6)*Or0*pow(a, 19) + 239616*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17) + 53248*pow(Ol, 5)*Om0*Or0*pow(a, 16) - 8192*pow(Ol, 5)*pow(Or0, 2)*pow(a, 15) + 225792*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14) + 363008*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13) - 243712*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12) - 55344*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11) + 417664*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10) + 503296*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9) - 195084*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8) - 447296*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7) - 244352*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6) + 19944*Ol*pow(Om0, 6)*pow(a, 5) + 32960*Ol*pow(Om0, 5)*Or0*pow(a, 4) + 13328*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3) + 327*pow(Om0, 7)*pow(a, 2) + 520*pow(Om0, 6)*Or0*a + 196*pow(Om0, 5)*pow(Or0, 2))/pow(4*Ol*pow(a, 3) + Om0, 9) + (4.0/3.0)*Ol*pow(a, 2)*b*(72*pow(Ol, 2)*Om0*pow(a, 7) + 32*pow(Ol, 2)*Or0*pow(a, 6) + 63*Ol*pow(Om0, 2)*pow(a, 4) + 88*Ol*Om0*Or0*pow(a, 3) - 9*pow(Om0, 3)*a - 7*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 4) - 1
;
   			
   			// dark energy density
   			rho_fR =  exp(-3*pow(b, 2)*((-2621440*pow(Ol, 7)*Om0*pow(a, 21) - 786432*pow(Ol, 7)*Or0*pow(a, 20) - 4849664*pow(Ol, 6)*pow(Om0, 2)*pow(a, 18) - 1048576*pow(Ol, 6)*Om0*Or0*pow(a, 17) + 65536*pow(Ol, 6)*pow(Or0, 2)*pow(a, 16) - 4030464*pow(Ol, 5)*pow(Om0, 3)*pow(a, 15) - 2768896*pow(Ol, 5)*pow(Om0, 2)*Or0*pow(a, 14) + 1441792*pow(Ol, 5)*Om0*pow(Or0, 2)*pow(a, 13) - 964352*pow(Ol, 4)*pow(Om0, 4)*pow(a, 12) - 2801664*pow(Ol, 4)*pow(Om0, 3)*Or0*pow(a, 11) - 1966080*pow(Ol, 4)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 10) + 639488*pow(Ol, 3)*pow(Om0, 5)*pow(a, 9) + 1307648*pow(Ol, 3)*pow(Om0, 4)*Or0*pow(a, 8) + 630784*pow(Ol, 3)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 7) + 9024*pow(Ol, 2)*pow(Om0, 6)*pow(a, 6) + 26624*pow(Ol, 2)*pow(Om0, 5)*Or0*pow(a, 5) + 12544*pow(Ol, 2)*pow(Om0, 4)*pow(Or0, 2)*pow(a, 4) - 352*Ol*pow(Om0, 7)*pow(a, 3) - 11*pow(Om0, 8))/(25165824*pow(Ol, 8)*pow(a, 24) + 50331648*pow(Ol, 7)*Om0*pow(a, 21) + 44040192*pow(Ol, 6)*pow(Om0, 2)*pow(a, 18) + 22020096*pow(Ol, 5)*pow(Om0, 3)*pow(a, 15) + 6881280*pow(Ol, 4)*pow(Om0, 4)*pow(a, 12) + 1376256*pow(Ol, 3)*pow(Om0, 5)*pow(a, 9) + 172032*pow(Ol, 2)*pow(Om0, 6)*pow(a, 6) + 12288*Ol*pow(Om0, 7)*pow(a, 3) + 384*pow(Om0, 8)) - (-2621440*pow(Ol, 7)*Om0 - 786432*pow(Ol, 7)*Or0 - 4849664*pow(Ol, 6)*pow(Om0, 2) - 1048576*pow(Ol, 6)*Om0*Or0 + 65536*pow(Ol, 6)*pow(Or0, 2) - 4030464*pow(Ol, 5)*pow(Om0, 3) - 2768896*pow(Ol, 5)*pow(Om0, 2)*Or0 + 1441792*pow(Ol, 5)*Om0*pow(Or0, 2) - 964352*pow(Ol, 4)*pow(Om0, 4) - 2801664*pow(Ol, 4)*pow(Om0, 3)*Or0 - 1966080*pow(Ol, 4)*pow(Om0, 2)*pow(Or0, 2) + 639488*pow(Ol, 3)*pow(Om0, 5) + 1307648*pow(Ol, 3)*pow(Om0, 4)*Or0 + 630784*pow(Ol, 3)*pow(Om0, 3)*pow(Or0, 2) + 9024*pow(Ol, 2)*pow(Om0, 6) + 26624*pow(Ol, 2)*pow(Om0, 5)*Or0 + 12544*pow(Ol, 2)*pow(Om0, 4)*pow(Or0, 2) - 352*Ol*pow(Om0, 7) - 11*pow(Om0, 8))/(25165824*pow(Ol, 8) + 50331648*pow(Ol, 7)*Om0 + 44040192*pow(Ol, 6)*pow(Om0, 2) + 22020096*pow(Ol, 5)*pow(Om0, 3) + 6881280*pow(Ol, 4)*pow(Om0, 4) + 1376256*pow(Ol, 3)*pow(Om0, 5) + 172032*pow(Ol, 2)*pow(Om0, 6) + 12288*Ol*pow(Om0, 7) + 384*pow(Om0, 8))) - 3*b*((-192*pow(Ol, 2)*Om0*pow(a, 6) - 64*pow(Ol, 2)*Or0*pow(a, 5) - 132*Ol*pow(Om0, 2)*pow(a, 3) - 112*Ol*Om0*Or0*pow(a, 2) - 3*pow(Om0, 3))/(1536*pow(Ol, 3)*pow(a, 9) + 1152*pow(Ol, 2)*Om0*pow(a, 6) + 288*Ol*pow(Om0, 2)*pow(a, 3) + 24*pow(Om0, 3)) - (-192*pow(Ol, 2)*Om0 - 64*pow(Ol, 2)*Or0 - 132*Ol*pow(Om0, 2) - 112*Ol*Om0*Or0 - 3*pow(Om0, 3))/(1536*pow(Ol, 3) + 1152*pow(Ol, 2)*Om0 + 288*Ol*pow(Om0, 2) + 24*pow(Om0, 3))))
;	

	pvecback[pba->index_bg_rho_ds] = rho_fR * pba->Omega0_ds * pow(pba->H0,2);
        pvecback[pba->index_bg_p_ds]   = w_fR *pvecback[pba->index_bg_rho_ds];
    		
  
	}if (l == 2){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = HS4 (Case 2) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = (1.0/2.0)*pow(2, 1.0/4.0)*pow(pow(H0, 10)*(fR0 - 1)*pow(3*Om0 + 4*Or0 - 4, 5), 1.0/4.0)/pow(pow(H0, 2)*Ol, 5.0/4.0);

   			w_fR = (4.0/3.0)*pow(Ol, 8)*pow(a, 22)*pow(b, 8)*(3932160*pow(Ol, 6)*Om0*pow(a, 20) + 589824*pow(Ol, 6)*Or0*pow(a, 19) + 12201984*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17) - 2289664*pow(Ol, 5)*Om0*Or0*pow(a, 16) - 65536*pow(Ol, 5)*pow(Or0, 2)*pow(a, 15) + 1301760*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14) + 54906880*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13) - 10833920*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12) - 46567680*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11) - 11008000*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10) + 61829120*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9) - 25864560*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8) - 95334400*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7) - 67425280*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6) + 16514640*Ol*pow(Om0, 6)*pow(a, 5) + 31165456*Ol*pow(Om0, 5)*Or0*pow(a, 4) + 14671360*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3) + 48021*pow(Om0, 7)*pow(a, 2) + 93104*pow(Om0, 6)*Or0*a + 45056*pow(Om0, 5)*pow(Or0, 2))/pow(4*Ol*pow(a, 3) + Om0, 15) + (4.0/3.0)*pow(Ol, 4)*pow(a, 11)*pow(b, 4)*(576*pow(Ol, 2)*Om0*pow(a, 7) + 128*pow(Ol, 2)*Or0*pow(a, 6) + 288*Ol*pow(Om0, 2)*pow(a, 4) + 784*Ol*Om0*Or0*pow(a, 3) - 369*pow(Om0, 3)*a - 352*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 7) - 1;
   			
   			// dark energy density
   			rho_fR = exp(-3*pow(b, 8)*((-171798691840*pow(Ol, 13)*Om0*pow(a, 39) - 19327352832*pow(Ol, 13)*Or0*pow(a, 38) - 545729282048*pow(Ol, 12)*pow(Om0, 2)*pow(a, 36) + 16642998272*pow(Ol, 12)*Om0*Or0*pow(a, 35) + 1073741824*pow(Ol, 12)*pow(Or0, 2)*pow(a, 34) - 564687536128*pow(Ol, 11)*pow(Om0, 3)*pow(a, 33) - 705112834048*pow(Ol, 11)*pow(Om0, 2)*Or0*pow(a, 32) + 129922760704*pow(Ol, 11)*Om0*pow(Or0, 2)*pow(a, 31) + 120420564992*pow(Ol, 10)*pow(Om0, 4)*pow(a, 30) - 322927853568*pow(Ol, 10)*pow(Om0, 3)*Or0*pow(a, 29) - 506940358656*pow(Ol, 10)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 28) + 286218256384*pow(Ol, 9)*pow(Om0, 5)*pow(a, 27) + 634652721152*pow(Ol, 9)*pow(Om0, 4)*Or0*pow(a, 26) + 311116693504*pow(Ol, 9)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 25) - 12924092416*pow(Ol, 8)*pow(Om0, 6)*pow(a, 24) + 2122317824*pow(Ol, 8)*pow(Om0, 5)*Or0*pow(a, 23) + 1073741824*pow(Ol, 8)*pow(Om0, 4)*pow(Or0, 2)*pow(a, 22) - 3992322048*pow(Ol, 7)*pow(Om0, 7)*pow(a, 21) - 873320448*pow(Ol, 6)*pow(Om0, 8)*pow(a, 18) - 145553408*pow(Ol, 5)*pow(Om0, 9)*pow(a, 15) - 18194176*pow(Ol, 4)*pow(Om0, 10)*pow(a, 12) - 1654016*pow(Ol, 3)*pow(Om0, 11)*pow(a, 9) - 103376*pow(Ol, 2)*pow(Om0, 12)*pow(a, 6) - 3976*Ol*pow(Om0, 13)*pow(a, 3) - 71*pow(Om0, 14))/(105553116266496*pow(Ol, 14)*pow(a, 42) + 369435906932736*pow(Ol, 13)*Om0*pow(a, 39) + 600333348765696*pow(Ol, 12)*pow(Om0, 2)*pow(a, 36) + 600333348765696*pow(Ol, 11)*pow(Om0, 3)*pow(a, 33) + 412729177276416*pow(Ol, 10)*pow(Om0, 4)*pow(a, 30) + 206364588638208*pow(Ol, 9)*pow(Om0, 5)*pow(a, 27) + 77386720739328*pow(Ol, 8)*pow(Om0, 6)*pow(a, 24) + 22110491639808*pow(Ol, 7)*pow(Om0, 7)*pow(a, 21) + 4836670046208*pow(Ol, 6)*pow(Om0, 8)*pow(a, 18) + 806111674368*pow(Ol, 5)*pow(Om0, 9)*pow(a, 15) + 100763959296*pow(Ol, 4)*pow(Om0, 10)*pow(a, 12) + 9160359936*pow(Ol, 3)*pow(Om0, 11)*pow(a, 9) + 572522496*pow(Ol, 2)*pow(Om0, 12)*pow(a, 6) + 22020096*Ol*pow(Om0, 13)*pow(a, 3) + 393216*pow(Om0, 14)) - (-171798691840*pow(Ol, 13)*Om0 - 19327352832*pow(Ol, 13)*Or0 - 545729282048*pow(Ol, 12)*pow(Om0, 2) + 16642998272*pow(Ol, 12)*Om0*Or0 + 1073741824*pow(Ol, 12)*pow(Or0, 2) - 564687536128*pow(Ol, 11)*pow(Om0, 3) - 705112834048*pow(Ol, 11)*pow(Om0, 2)*Or0 + 129922760704*pow(Ol, 11)*Om0*pow(Or0, 2) + 120420564992*pow(Ol, 10)*pow(Om0, 4) - 322927853568*pow(Ol, 10)*pow(Om0, 3)*Or0 - 506940358656*pow(Ol, 10)*pow(Om0, 2)*pow(Or0, 2) + 286218256384*pow(Ol, 9)*pow(Om0, 5) + 634652721152*pow(Ol, 9)*pow(Om0, 4)*Or0 + 311116693504*pow(Ol, 9)*pow(Om0, 3)*pow(Or0, 2) - 12924092416*pow(Ol, 8)*pow(Om0, 6) + 2122317824*pow(Ol, 8)*pow(Om0, 5)*Or0 + 1073741824*pow(Ol, 8)*pow(Om0, 4)*pow(Or0, 2) - 3992322048*pow(Ol, 7)*pow(Om0, 7) - 873320448*pow(Ol, 6)*pow(Om0, 8) - 145553408*pow(Ol, 5)*pow(Om0, 9) - 18194176*pow(Ol, 4)*pow(Om0, 10) - 1654016*pow(Ol, 3)*pow(Om0, 11) - 103376*pow(Ol, 2)*pow(Om0, 12) - 3976*Ol*pow(Om0, 13) - 71*pow(Om0, 14))/(105553116266496*pow(Ol, 14) + 369435906932736*pow(Ol, 13)*Om0 + 600333348765696*pow(Ol, 12)*pow(Om0, 2) + 600333348765696*pow(Ol, 11)*pow(Om0, 3) + 412729177276416*pow(Ol, 10)*pow(Om0, 4) + 206364588638208*pow(Ol, 9)*pow(Om0, 5) + 77386720739328*pow(Ol, 8)*pow(Om0, 6) + 22110491639808*pow(Ol, 7)*pow(Om0, 7) + 4836670046208*pow(Ol, 6)*pow(Om0, 8) + 806111674368*pow(Ol, 5)*pow(Om0, 9) + 100763959296*pow(Ol, 4)*pow(Om0, 10) + 9160359936*pow(Ol, 3)*pow(Om0, 11) + 572522496*pow(Ol, 2)*pow(Om0, 12) + 22020096*Ol*pow(Om0, 13) + 393216*pow(Om0, 14))) - 3*pow(b, 4)*((-49152*pow(Ol, 5)*Om0*pow(a, 15) - 8192*pow(Ol, 5)*Or0*pow(a, 14) - 43008*pow(Ol, 4)*pow(Om0, 2)*pow(a, 12) - 32768*pow(Ol, 4)*Om0*Or0*pow(a, 11) - 3840*pow(Ol, 3)*pow(Om0, 3)*pow(a, 9) - 720*pow(Ol, 2)*pow(Om0, 4)*pow(a, 6) - 72*Ol*pow(Om0, 5)*pow(a, 3) - 3*pow(Om0, 6))/(3145728*pow(Ol, 6)*pow(a, 18) + 4718592*pow(Ol, 5)*Om0*pow(a, 15) + 2949120*pow(Ol, 4)*pow(Om0, 2)*pow(a, 12) + 983040*pow(Ol, 3)*pow(Om0, 3)*pow(a, 9) + 184320*pow(Ol, 2)*pow(Om0, 4)*pow(a, 6) + 18432*Ol*pow(Om0, 5)*pow(a, 3) + 768*pow(Om0, 6)) - (-49152*pow(Ol, 5)*Om0 - 8192*pow(Ol, 5)*Or0 - 43008*pow(Ol, 4)*pow(Om0, 2) - 32768*pow(Ol, 4)*Om0*Or0 - 3840*pow(Ol, 3)*pow(Om0, 3) - 720*pow(Ol, 2)*pow(Om0, 4) - 72*Ol*pow(Om0, 5) - 3*pow(Om0, 6))/(3145728*pow(Ol, 6) + 4718592*pow(Ol, 5)*Om0 + 2949120*pow(Ol, 4)*pow(Om0, 2) + 983040*pow(Ol, 3)*pow(Om0, 3) + 184320*pow(Ol, 2)*pow(Om0, 4) + 18432*Ol*pow(Om0, 5) + 768*pow(Om0, 6))));	
    		
    		pvecback[pba->index_bg_rho_ds] = rho_fR * pba->Omega0_ds * pow(pba->H0,2);
        	pvecback[pba->index_bg_p_ds]   = w_fR *pvecback[pba->index_bg_rho_ds];
    		}
    		
    		if (l == 3){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = Star1 (Case 3) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = (1.0/2.0)*sqrt(pow(H0, 6)*(fR0 - 1)*pow(3*Om0 + 4*Or0 - 4, 3))/pow(pow(H0, 2)*Ol, 3.0/2.0);

   			w_fR = (2.0/3.0)*pow(Ol, 4)*pow(a, 10)*pow(b, 4)*(884736*pow(Ol, 6)*Om0*pow(a, 20) + 229376*pow(Ol, 6)*Or0*pow(a, 19) + 2016768*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17) - 57344*pow(Ol, 5)*Om0*Or0*pow(a, 16) - 32768*pow(Ol, 5)*pow(Or0, 2)*pow(a, 15) + 1936896*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14) + 5477888*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13) - 2050048*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12) - 2307840*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11) + 3363200*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10) + 6776320*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9) - 2905440*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8) - 7610144*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7) - 4563584*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6) + 641490*Ol*pow(Om0, 6)*pow(a, 5) + 1157960*Ol*pow(Om0, 5)*Or0*pow(a, 4) + 519200*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3) + 4845*pow(Om0, 7)*pow(a, 2) + 8844*pow(Om0, 6)*Or0*a + 4000*pow(Om0, 5)*pow(Or0, 2))/pow(4*Ol*pow(a, 3) + Om0, 11) + (2.0/3.0)*pow(Ol, 2)*pow(a, 5)*pow(b, 2)*(384*pow(Ol, 2)*Om0*pow(a, 7) + 128*pow(Ol, 2)*Or0*pow(a, 6) + 300*Ol*pow(Om0, 2)*pow(a, 4) + 496*Ol*Om0*Or0*pow(a, 3) - 111*pow(Om0, 3)*a - 100*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 5) - 1;
   			
   			// dark energy density
   			rho_fR = exp(-3*pow(b, 4)*((-37748736*pow(Ol, 9)*Om0*pow(a, 27) - 7340032*pow(Ol, 9)*Or0*pow(a, 26) - 85491712*pow(Ol, 8)*pow(Om0, 2)*pow(a, 24) - 5767168*pow(Ol, 8)*Om0*Or0*pow(a, 23) + 524288*pow(Ol, 8)*pow(Or0, 2)*pow(a, 22) - 84541440*pow(Ol, 7)*pow(Om0, 3)*pow(a, 21) - 73433088*pow(Ol, 7)*pow(Om0, 2)*Or0*pow(a, 20) + 24117248*pow(Ol, 7)*Om0*pow(Or0, 2)*pow(a, 19) - 12369920*pow(Ol, 6)*pow(Om0, 4)*pow(a, 18) - 61358080*pow(Ol, 6)*pow(Om0, 3)*Or0*pow(a, 17) - 53772288*pow(Ol, 6)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 16) + 21082112*pow(Ol, 5)*pow(Om0, 5)*pow(a, 15) + 44582912*pow(Ol, 5)*pow(Om0, 4)*Or0*pow(a, 14) + 21708800*pow(Ol, 5)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 13) - 169600*pow(Ol, 4)*pow(Om0, 6)*pow(a, 12) + 411648*pow(Ol, 4)*pow(Om0, 5)*Or0*pow(a, 11) + 204800*pow(Ol, 4)*pow(Om0, 4)*pow(Or0, 2)*pow(a, 10) - 53760*pow(Ol, 3)*pow(Om0, 7)*pow(a, 9) - 5040*pow(Ol, 2)*pow(Om0, 8)*pow(a, 6) - 280*Ol*pow(Om0, 9)*pow(a, 3) - 7*pow(Om0, 10))/(805306368*pow(Ol, 10)*pow(a, 30) + 2013265920*pow(Ol, 9)*Om0*pow(a, 27) + 2264924160*pow(Ol, 8)*pow(Om0, 2)*pow(a, 24) + 1509949440*pow(Ol, 7)*pow(Om0, 3)*pow(a, 21) + 660602880*pow(Ol, 6)*pow(Om0, 4)*pow(a, 18) + 198180864*pow(Ol, 5)*pow(Om0, 5)*pow(a, 15) + 41287680*pow(Ol, 4)*pow(Om0, 6)*pow(a, 12) + 5898240*pow(Ol, 3)*pow(Om0, 7)*pow(a, 9) + 552960*pow(Ol, 2)*pow(Om0, 8)*pow(a, 6) + 30720*Ol*pow(Om0, 9)*pow(a, 3) + 768*pow(Om0, 10)) - (-37748736*pow(Ol, 9)*Om0 - 7340032*pow(Ol, 9)*Or0 - 85491712*pow(Ol, 8)*pow(Om0, 2) - 5767168*pow(Ol, 8)*Om0*Or0 + 524288*pow(Ol, 8)*pow(Or0, 2) - 84541440*pow(Ol, 7)*pow(Om0, 3) - 73433088*pow(Ol, 7)*pow(Om0, 2)*Or0 + 24117248*pow(Ol, 7)*Om0*pow(Or0, 2) - 12369920*pow(Ol, 6)*pow(Om0, 4) - 61358080*pow(Ol, 6)*pow(Om0, 3)*Or0 - 53772288*pow(Ol, 6)*pow(Om0, 2)*pow(Or0, 2) + 21082112*pow(Ol, 5)*pow(Om0, 5) + 44582912*pow(Ol, 5)*pow(Om0, 4)*Or0 + 21708800*pow(Ol, 5)*pow(Om0, 3)*pow(Or0, 2) - 169600*pow(Ol, 4)*pow(Om0, 6) + 411648*pow(Ol, 4)*pow(Om0, 5)*Or0 + 204800*pow(Ol, 4)*pow(Om0, 4)*pow(Or0, 2) - 53760*pow(Ol, 3)*pow(Om0, 7) - 5040*pow(Ol, 2)*pow(Om0, 8) - 280*Ol*pow(Om0, 9) - 7*pow(Om0, 10))/(805306368*pow(Ol, 10) + 2013265920*pow(Ol, 9)*Om0 + 2264924160*pow(Ol, 8)*pow(Om0, 2) + 1509949440*pow(Ol, 7)*pow(Om0, 3) + 660602880*pow(Ol, 6)*pow(Om0, 4) + 198180864*pow(Ol, 5)*pow(Om0, 5) + 41287680*pow(Ol, 4)*pow(Om0, 6) + 5898240*pow(Ol, 3)*pow(Om0, 7) + 552960*pow(Ol, 2)*pow(Om0, 8) + 30720*Ol*pow(Om0, 9) + 768*pow(Om0, 10))) - 3*pow(b, 2)*((-512*pow(Ol, 3)*Om0*pow(a, 9) - 128*pow(Ol, 3)*Or0*pow(a, 8) - 392*pow(Ol, 2)*pow(Om0, 2)*pow(a, 6) - 320*pow(Ol, 2)*Om0*Or0*pow(a, 5) - 16*Ol*pow(Om0, 3)*pow(a, 3) - pow(Om0, 4))/(6144*pow(Ol, 4)*pow(a, 12) + 6144*pow(Ol, 3)*Om0*pow(a, 9) + 2304*pow(Ol, 2)*pow(Om0, 2)*pow(a, 6) + 384*Ol*pow(Om0, 3)*pow(a, 3) + 24*pow(Om0, 4)) - (-512*pow(Ol, 3)*Om0 - 128*pow(Ol, 3)*Or0 - 392*pow(Ol, 2)*pow(Om0, 2) - 320*pow(Ol, 2)*Om0*Or0 - 16*Ol*pow(Om0, 3) - pow(Om0, 4))/(6144*pow(Ol, 4) + 6144*pow(Ol, 3)*Om0 + 2304*pow(Ol, 2)*pow(Om0, 2) + 384*Ol*pow(Om0, 3) + 24*pow(Om0, 4))))
;	
    		
    		pvecback[pba->index_bg_rho_ds] = rho_fR * pba->Omega0_ds * pow(pba->H0,2);
        pvecback[pba->index_bg_p_ds]   = w_fR *pvecback[pba->index_bg_rho_ds];
        
         		
    		}if (l == 4){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = Star2 (Case 4) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0) DA MATHEMATICA
   			b = (1.0/2.0)*pow(2, 1.0/4.0)*pow(pow(H0, 10)*(fR0 - 1)*pow(3*Om0 + 4*Or0 - 4, 5), 1.0/4.0)/pow(pow(H0, 2)*Ol, 5.0/4.0);

   			w_fR = -4*pow(Ol, 6)*pow(a, 17)*pow(b, 6)*(768*pow(Ol, 2)*Om0*pow(a, 7) + 128*pow(Ol, 2)*Or0*pow(a, 6) + 132*Ol*pow(Om0, 2)*pow(a, 4) + 1072*Ol*Om0*Or0*pow(a, 3) - 771*pow(Om0, 3)*a - 748*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 9) + (4.0/3.0)*pow(Ol, 4)*pow(a, 11)*pow(b, 4)*(576*pow(Ol, 2)*Om0*pow(a, 7) + 128*pow(Ol, 2)*Or0*pow(a, 6) + 288*Ol*pow(Om0, 2)*pow(a, 4) + 784*Ol*Om0*Or0*pow(a, 3) - 369*pow(Om0, 3)*a - 352*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 7) - 1;
   			
   			// dark energy density
   			rho_fR =exp(-3*pow(b, 6)*((393216*pow(Ol, 7)*Om0*pow(a, 21) + 49152*pow(Ol, 7)*Or0*pow(a, 20) + 377856*pow(Ol, 6)*pow(Om0, 2)*pow(a, 18) + 270336*pow(Ol, 6)*Om0*Or0*pow(a, 17) + 57344*pow(Ol, 5)*pow(Om0, 3)*pow(a, 15) + 17920*pow(Ol, 4)*pow(Om0, 4)*pow(a, 12) + 3584*pow(Ol, 3)*pow(Om0, 5)*pow(a, 9) + 448*pow(Ol, 2)*pow(Om0, 6)*pow(a, 6) + 32*Ol*pow(Om0, 7)*pow(a, 3) + pow(Om0, 8))/(100663296*pow(Ol, 8)*pow(a, 24) + 201326592*pow(Ol, 7)*Om0*pow(a, 21) + 176160768*pow(Ol, 6)*pow(Om0, 2)*pow(a, 18) + 88080384*pow(Ol, 5)*pow(Om0, 3)*pow(a, 15) + 27525120*pow(Ol, 4)*pow(Om0, 4)*pow(a, 12) + 5505024*pow(Ol, 3)*pow(Om0, 5)*pow(a, 9) + 688128*pow(Ol, 2)*pow(Om0, 6)*pow(a, 6) + 49152*Ol*pow(Om0, 7)*pow(a, 3) + 1536*pow(Om0, 8)) - (393216*pow(Ol, 7)*Om0 + 49152*pow(Ol, 7)*Or0 + 377856*pow(Ol, 6)*pow(Om0, 2) + 270336*pow(Ol, 6)*Om0*Or0 + 57344*pow(Ol, 5)*pow(Om0, 3) + 17920*pow(Ol, 4)*pow(Om0, 4) + 3584*pow(Ol, 3)*pow(Om0, 5) + 448*pow(Ol, 2)*pow(Om0, 6) + 32*Ol*pow(Om0, 7) + pow(Om0, 8))/(100663296*pow(Ol, 8) + 201326592*pow(Ol, 7)*Om0 + 176160768*pow(Ol, 6)*pow(Om0, 2) + 88080384*pow(Ol, 5)*pow(Om0, 3) + 27525120*pow(Ol, 4)*pow(Om0, 4) + 5505024*pow(Ol, 3)*pow(Om0, 5) + 688128*pow(Ol, 2)*pow(Om0, 6) + 49152*Ol*pow(Om0, 7) + 1536*pow(Om0, 8))) - 3*pow(b, 4)*((-49152*pow(Ol, 5)*Om0*pow(a, 15) - 8192*pow(Ol, 5)*Or0*pow(a, 14) - 43008*pow(Ol, 4)*pow(Om0, 2)*pow(a, 12) - 32768*pow(Ol, 4)*Om0*Or0*pow(a, 11) - 3840*pow(Ol, 3)*pow(Om0, 3)*pow(a, 9) - 720*pow(Ol, 2)*pow(Om0, 4)*pow(a, 6) - 72*Ol*pow(Om0, 5)*pow(a, 3) - 3*pow(Om0, 6))/(3145728*pow(Ol, 6)*pow(a, 18) + 4718592*pow(Ol, 5)*Om0*pow(a, 15) + 2949120*pow(Ol, 4)*pow(Om0, 2)*pow(a, 12) + 983040*pow(Ol, 3)*pow(Om0, 3)*pow(a, 9) + 184320*pow(Ol, 2)*pow(Om0, 4)*pow(a, 6) + 18432*Ol*pow(Om0, 5)*pow(a, 3) + 768*pow(Om0, 6)) - (-49152*pow(Ol, 5)*Om0 - 8192*pow(Ol, 5)*Or0 - 43008*pow(Ol, 4)*pow(Om0, 2) - 32768*pow(Ol, 4)*Om0*Or0 - 3840*pow(Ol, 3)*pow(Om0, 3) - 720*pow(Ol, 2)*pow(Om0, 4) - 72*Ol*pow(Om0, 5) - 3*pow(Om0, 6))/(3145728*pow(Ol, 6) + 4718592*pow(Ol, 5)*Om0 + 2949120*pow(Ol, 4)*pow(Om0, 2) + 983040*pow(Ol, 3)*pow(Om0, 3) + 184320*pow(Ol, 2)*pow(Om0, 4) + 18432*Ol*pow(Om0, 5) + 768*pow(Om0, 6))));	
    		
    		
    	pvecback[pba->index_bg_rho_ds] = rho_fR * pba->Omega0_ds * pow(pba->H0,2);
        pvecback[pba->index_bg_p_ds]   = w_fR *pvecback[pba->index_bg_rho_ds];
    		
    		}if (l == 5){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = Exp (Case 5) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = -1.0/2.0*(fR0 - 1)*pow(3*Om0 + 4*Or0 - 4, 2)/pow(Ol, 2); // GIÀ CAMBIATO

   			w_fR = (1.0/3.0)*pow(Ol, 2)*pow(a, 4)*pow(b, 2)*(344064*pow(Ol, 6)*Om0*pow(a, 20) + 131072*pow(Ol, 6)*Or0*pow(a, 19) + 654336*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17) + 266240*pow(Ol, 5)*Om0*Or0*pow(a, 16) - 16384*pow(Ol, 5)*pow(Or0, 2)*pow(a, 15) + 536832*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14) + 839680*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13) - 487424*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12) - 104160*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11) + 859392*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10) + 1006592*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9) - 395640*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8) - 896128*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7) - 488704*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6) + 38412*Ol*pow(Om0, 6)*pow(a, 5) + 64816*Ol*pow(Om0, 5)*Or0*pow(a, 4) + 26656*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3) + 543*pow(Om0, 7)*pow(a, 2) + 940*pow(Om0, 6)*Or0*a + 392*pow(Om0, 5)*pow(Or0, 2))/pow(4*Ol*pow(a, 3) + Om0, 9) + (4.0/3.0)*Ol*pow(a, 2)*b*(72*pow(Ol, 2)*Om0*pow(a, 7) + 32*pow(Ol, 2)*Or0*pow(a, 6) + 63*Ol*pow(Om0, 2)*pow(a, 4) + 88*Ol*Om0*Or0*pow(a, 3) - 9*pow(Om0, 3)*a - 7*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 4) - 1;
   			
   			// dark energy density
   			rho_fR = exp(-3*pow(b, 2)*((-3670016*pow(Ol, 7)*Om0*pow(a, 21) - 1048576*pow(Ol, 7)*Or0*pow(a, 20) - 6701056*pow(Ol, 6)*pow(Om0, 2)*pow(a, 18) - 1966080*pow(Ol, 6)*Om0*Or0*pow(a, 17) + 65536*pow(Ol, 6)*pow(Or0, 2)*pow(a, 16) - 5259264*pow(Ol, 5)*pow(Om0, 3)*pow(a, 15) - 3522560*pow(Ol, 5)*pow(Om0, 2)*Or0*pow(a, 14) + 1441792*pow(Ol, 5)*Om0*pow(Or0, 2)*pow(a, 13) - 1365760*pow(Ol, 4)*pow(Om0, 4)*pow(a, 12) - 3063808*pow(Ol, 4)*pow(Om0, 3)*Or0*pow(a, 11) - 1966080*pow(Ol, 4)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 10) + 570880*pow(Ol, 3)*pow(Om0, 5)*pow(a, 9) + 1265664*pow(Ol, 3)*pow(Om0, 4)*Or0*pow(a, 8) + 630784*pow(Ol, 3)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 7) + 3072*pow(Ol, 2)*pow(Om0, 6)*pow(a, 6) + 24064*pow(Ol, 2)*pow(Om0, 5)*Or0*pow(a, 5) + 12544*pow(Ol, 2)*pow(Om0, 4)*pow(Or0, 2)*pow(a, 4) - 608*Ol*pow(Om0, 7)*pow(a, 3) - 19*pow(Om0, 8))/(25165824*pow(Ol, 8)*pow(a, 24) + 50331648*pow(Ol, 7)*Om0*pow(a, 21) + 44040192*pow(Ol, 6)*pow(Om0, 2)*pow(a, 18) + 22020096*pow(Ol, 5)*pow(Om0, 3)*pow(a, 15) + 6881280*pow(Ol, 4)*pow(Om0, 4)*pow(a, 12) + 1376256*pow(Ol, 3)*pow(Om0, 5)*pow(a, 9) + 172032*pow(Ol, 2)*pow(Om0, 6)*pow(a, 6) + 12288*Ol*pow(Om0, 7)*pow(a, 3) + 384*pow(Om0, 8)) - (-3670016*pow(Ol, 7)*Om0 - 1048576*pow(Ol, 7)*Or0 - 6701056*pow(Ol, 6)*pow(Om0, 2) - 1966080*pow(Ol, 6)*Om0*Or0 + 65536*pow(Ol, 6)*pow(Or0, 2) - 5259264*pow(Ol, 5)*pow(Om0, 3) - 3522560*pow(Ol, 5)*pow(Om0, 2)*Or0 + 1441792*pow(Ol, 5)*Om0*pow(Or0, 2) - 1365760*pow(Ol, 4)*pow(Om0, 4) - 3063808*pow(Ol, 4)*pow(Om0, 3)*Or0 - 1966080*pow(Ol, 4)*pow(Om0, 2)*pow(Or0, 2) + 570880*pow(Ol, 3)*pow(Om0, 5) + 1265664*pow(Ol, 3)*pow(Om0, 4)*Or0 + 630784*pow(Ol, 3)*pow(Om0, 3)*pow(Or0, 2) + 3072*pow(Ol, 2)*pow(Om0, 6) + 24064*pow(Ol, 2)*pow(Om0, 5)*Or0 + 12544*pow(Ol, 2)*pow(Om0, 4)*pow(Or0, 2) - 608*Ol*pow(Om0, 7) - 19*pow(Om0, 8))/(25165824*pow(Ol, 8) + 50331648*pow(Ol, 7)*Om0 + 44040192*pow(Ol, 6)*pow(Om0, 2) + 22020096*pow(Ol, 5)*pow(Om0, 3) + 6881280*pow(Ol, 4)*pow(Om0, 4) + 1376256*pow(Ol, 3)*pow(Om0, 5) + 172032*pow(Ol, 2)*pow(Om0, 6) + 12288*Ol*pow(Om0, 7) + 384*pow(Om0, 8))) - 3*b*((-192*pow(Ol, 2)*Om0*pow(a, 6) - 64*pow(Ol, 2)*Or0*pow(a, 5) - 132*Ol*pow(Om0, 2)*pow(a, 3) - 112*Ol*Om0*Or0*pow(a, 2) - 3*pow(Om0, 3))/(1536*pow(Ol, 3)*pow(a, 9) + 1152*pow(Ol, 2)*Om0*pow(a, 6) + 288*Ol*pow(Om0, 2)*pow(a, 3) + 24*pow(Om0, 3)) - (-192*pow(Ol, 2)*Om0 - 64*pow(Ol, 2)*Or0 - 132*Ol*pow(Om0, 2) - 112*Ol*Om0*Or0 - 3*pow(Om0, 3))/(1536*pow(Ol, 3) + 1152*pow(Ol, 2)*Om0 + 288*Ol*pow(Om0, 2) + 24*pow(Om0, 3))));	
    		
    	pvecback[pba->index_bg_rho_ds] = rho_fR * pba->Omega0_ds * pow(pba->H0,2);
        pvecback[pba->index_bg_p_ds]   = w_fR *pvecback[pba->index_bg_rho_ds];
    		
    		}if (l == 6){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = Tsjk (Case 6) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = -1.0/2.0*sqrt(pow(H0, 6)*(fR0 - 1)*pow(3*Om0 + 4*Or0 - 4, 3))/pow(pow(H0, 2)*Ol, 3.0/2.0);

   			w_fR = (2.0/3.0)*pow(Ol, 4)*pow(a, 10)*pow(b, 3)*(32*pow(Ol, 2)*pow(a, 7) - 32*Ol*Om0*pow(a, 4) - 16*Ol*Or0*pow(a, 3) - 37*pow(Om0, 2)*a - 40*Om0*Or0)*(384*pow(Ol, 2)*Om0*pow(a, 7) + 128*pow(Ol, 2)*Or0*pow(a, 6) + 300*Ol*pow(Om0, 2)*pow(a, 4) + 496*Ol*Om0*Or0*pow(a, 3) - 111*pow(Om0, 3)*a - 100*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 9) + (2.0/3.0)*pow(Ol, 2)*pow(a, 5)*pow(b, 2)*(384*pow(Ol, 2)*Om0*pow(a, 7) + 128*pow(Ol, 2)*Or0*pow(a, 6) + 300*Ol*pow(Om0, 2)*pow(a, 4) + 496*Ol*Om0*Or0*pow(a, 3) - 111*pow(Om0, 3)*a - 100*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 5) - 1
;
   			
   			// dark energy density
   			rho_fR = 1;	
    		
    		
    	pvecback[pba->index_bg_rho_ds] = rho_fR * pba->Omega0_ds * pow(pba->H0,2);
        pvecback[pba->index_bg_p_ds]   = w_fR *pvecback[pba->index_bg_rho_ds];
    		
    		}if (l == 7){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = Log (Case 7) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = (1.0/2.0)*(fR0 - 1)*(3*Om0 + 4*Or0 - 4)*pow(log((-3*Om0 - 4*Or0 + 4)/Ol), 2)/Ol;

   			w_fR = pow(b, 2)*(147456*pow(Ol, 6)*Om0*pow(a, 20)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 344064*pow(Ol, 6)*Om0*pow(a, 20)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 122880*pow(Ol, 6)*Om0*pow(a, 20)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 114688*pow(Ol, 6)*Or0*pow(a, 19)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 32768*pow(Ol, 6)*Or0*pow(a, 19)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 6144*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + 235008*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 664320*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 192768*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 152064*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 28672*pow(Ol, 5)*Om0*Or0*pow(a, 16)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + 215040*pow(Ol, 5)*Om0*Or0*pow(a, 16)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 200704*pow(Ol, 5)*Om0*Or0*pow(a, 16)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) - 159744*pow(Ol, 5)*Om0*Or0*pow(a, 16)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 16384*pow(Ol, 5)*pow(Or0, 2)*pow(a, 15)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) - 15360*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + 101952*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 575808*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 357312*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 17280*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 580608*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*log(4 + Om0/(Ol*pow(a, 3))) - 1024*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + 88320*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 601088*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 1310208*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 1465344*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 152576*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) - 466944*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) - 405504*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 14016*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) - 144*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 105264*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) - 398472*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 1181952*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 1143072*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*log(4 + Om0/(Ol*pow(a, 3))) + 217728*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11) - 12032*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + 123968*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 876416*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 1301184*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 155520*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 1669248*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10)*log(4 + Om0/(Ol*pow(a, 3))) + 111104*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 779520*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 1790208*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 1423872*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 3792*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) - 50112*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) - 245472*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) - 707844*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 890568*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 258552*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*log(4 + Om0/(Ol*pow(a, 3))) + 435456*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8) - 4160*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) - 99024*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) - 526208*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) - 1553520*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 2379456*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 1378944*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7)*log(4 + Om0/(Ol*pow(a, 3))) + 435456*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7) - 47488*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) - 271296*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) - 822528*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 1444608*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 1088640*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6)*log(4 + Om0/(Ol*pow(a, 3))) - 348*Ol*pow(Om0, 6)*pow(a, 5)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) - 513*Ol*pow(Om0, 6)*pow(a, 5)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 11832*Ol*pow(Om0, 6)*pow(a, 5)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 63876*Ol*pow(Om0, 6)*pow(a, 5)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 192888*Ol*pow(Om0, 6)*pow(a, 5)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 303912*Ol*pow(Om0, 6)*pow(a, 5)*log(4 + Om0/(Ol*pow(a, 3))) + 217728*Ol*pow(Om0, 6)*pow(a, 5) - 448*Ol*pow(Om0, 5)*Or0*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) - 2724*Ol*pow(Om0, 5)*Or0*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 14000*Ol*pow(Om0, 5)*Or0*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 110796*Ol*pow(Om0, 5)*Or0*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 346248*Ol*pow(Om0, 5)*Or0*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 576072*Ol*pow(Om0, 5)*Or0*pow(a, 4)*log(4 + Om0/(Ol*pow(a, 3))) + 435456*Ol*pow(Om0, 5)*Or0*pow(a, 4) - 1952*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 2640*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 47448*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 153792*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 272160*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3)*log(4 + Om0/(Ol*pow(a, 3))) + 217728*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3) - 3*pow(Om0, 7)*pow(a, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + 63*pow(Om0, 7)*pow(a, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 582*pow(Om0, 7)*pow(a, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 1260*pow(Om0, 7)*pow(a, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 1296*pow(Om0, 7)*pow(a, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 8*pow(Om0, 6)*Or0*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) - 8*pow(Om0, 6)*Or0*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 780*pow(Om0, 6)*Or0*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 2556*pow(Om0, 6)*Or0*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 2592*pow(Om0, 6)*Or0*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 64*pow(Om0, 5)*pow(Or0, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 192*pow(Om0, 5)*pow(Or0, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 1296*pow(Om0, 5)*pow(Or0, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 1296*pow(Om0, 5)*pow(Or0, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2))/(49152*pow(Ol, 7)*pow(a, 23)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 86016*pow(Ol, 6)*Om0*pow(a, 20)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 64512*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 26880*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 6720*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 1008*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 84*Ol*pow(Om0, 6)*pow(a, 5)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 3*pow(Om0, 7)*pow(a, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9)) + b*(192*pow(Ol, 2)*Om0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 192*pow(Ol, 2)*Om0*pow(a, 7)*log(4 + Om0/(Ol*pow(a, 3))) + 128*pow(Ol, 2)*Or0*pow(a, 6)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 168*Ol*pow(Om0, 2)*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 204*Ol*pow(Om0, 2)*pow(a, 4)*log(4 + Om0/(Ol*pow(a, 3))) - 108*Ol*pow(Om0, 2)*pow(a, 4) + 208*Ol*Om0*Or0*pow(a, 3)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 288*Ol*Om0*Or0*pow(a, 3)*log(4 + Om0/(Ol*pow(a, 3))) + 3*pow(Om0, 3)*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 42*pow(Om0, 3)*a*log(4 + Om0/(Ol*pow(a, 3))) - 108*pow(Om0, 3)*a + 8*pow(Om0, 2)*Or0*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 36*pow(Om0, 2)*Or0*log(4 + Om0/(Ol*pow(a, 3))) - 108*pow(Om0, 2)*Or0)/(192*pow(Ol, 3)*pow(a, 10)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 144*pow(Ol, 2)*Om0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 36*Ol*pow(Om0, 2)*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 3*pow(Om0, 3)*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 4)) - 1;
   			
   			// dark energy density
   			rho_fR = exp(-3*pow(b, 2)*((18144*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) + 36288*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) + 36288*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 7) + 18144*Ol*pow(Om0, 5)*pow(a, 5) + 36288*Ol*pow(Om0, 4)*Or0*pow(a, 4) + 18144*Ol*pow(Om0, 3)*pow(Or0, 2)*pow(a, 3) + (-4096*pow(Ol, 6)*pow(a, 20) - 6144*pow(Ol, 5)*Om0*pow(a, 17) - 3840*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) - 1280*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) - 240*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) - 24*Ol*pow(Om0, 5)*pow(a, 5) - pow(Om0, 6)*pow(a, 2))*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + (-24192*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) - 33696*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) - 69120*pow(Ol, 3)*pow(Om0, 2)*Or0*pow(a, 10) + 16848*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) - 29376*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 7) - 44928*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 6) + 26352*Ol*pow(Om0, 5)*pow(a, 5) + 51408*Ol*pow(Om0, 4)*Or0*pow(a, 4) + 25056*Ol*pow(Om0, 3)*pow(Or0, 2)*pow(a, 3))*log(4 + Om0/(Ol*pow(a, 3))) + (4096*pow(Ol, 6)*pow(a, 20) + 4096*pow(Ol, 5)*Om0*pow(a, 17) + 512*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) - 2048*pow(Ol, 4)*Om0*Or0*pow(a, 13) + 384*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) - 512*pow(Ol, 3)*pow(Om0, 2)*Or0*pow(a, 10) + 496*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) + 384*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 7) + 152*Ol*pow(Om0, 5)*pow(a, 5) + 160*Ol*pow(Om0, 4)*Or0*pow(a, 4) + 14*pow(Om0, 6)*pow(a, 2) + 16*pow(Om0, 5)*Or0*a)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + (-15360*pow(Ol, 5)*Om0*pow(a, 17) - 49344*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) + 9216*pow(Ol, 4)*Om0*Or0*pow(a, 13) - 69552*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) - 71424*pow(Ol, 3)*pow(Om0, 2)*Or0*pow(a, 10) + 18432*pow(Ol, 3)*Om0*pow(Or0, 2)*pow(a, 9) - 16728*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) - 73440*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 7) - 55872*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 6) + 17040*Ol*pow(Om0, 5)*pow(a, 5) + 32760*Ol*pow(Om0, 4)*Or0*pow(a, 4) + 15696*Ol*pow(Om0, 3)*pow(Or0, 2)*pow(a, 3) + 144*pow(Om0, 6)*pow(a, 2) + 288*pow(Om0, 5)*Or0*a + 144*pow(Om0, 4)*pow(Or0, 2))*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + (4096*pow(Ol, 6)*pow(a, 20) - 37888*pow(Ol, 5)*Om0*pow(a, 17) - 4096*pow(Ol, 5)*Or0*pow(a, 16) - 79808*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) - 22528*pow(Ol, 4)*Om0*Or0*pow(a, 13) - 64928*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) - 43008*pow(Ol, 3)*pow(Om0, 2)*Or0*pow(a, 10) + 21504*pow(Ol, 3)*Om0*pow(Or0, 2)*pow(a, 9) - 19528*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) - 42880*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 7) - 23616*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 6) + 5800*Ol*pow(Om0, 5)*pow(a, 5) + 11864*Ol*pow(Om0, 4)*Or0*pow(a, 4) + 6048*Ol*pow(Om0, 3)*pow(Or0, 2)*pow(a, 3) + 168*pow(Om0, 6)*pow(a, 2) + 360*pow(Om0, 5)*Or0*a + 192*pow(Om0, 4)*pow(Or0, 2))*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + (13312*pow(Ol, 6)*pow(a, 20) - 1536*pow(Ol, 5)*Om0*pow(a, 17) - 14336*pow(Ol, 5)*Or0*pow(a, 16) - 19008*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) - 22528*pow(Ol, 4)*Om0*Or0*pow(a, 13) + 1024*pow(Ol, 4)*pow(Or0, 2)*pow(a, 12) - 11872*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) - 11072*pow(Ol, 3)*pow(Om0, 2)*Or0*pow(a, 10) + 7168*pow(Ol, 3)*Om0*pow(Or0, 2)*pow(a, 9) - 2148*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) - 6560*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 7) - 3072*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 6) + 1530*Ol*pow(Om0, 5)*pow(a, 5) + 2636*Ol*pow(Om0, 4)*Or0*pow(a, 4) + 1216*Ol*pow(Om0, 3)*pow(Or0, 2)*pow(a, 3) + 97*pow(Om0, 6)*pow(a, 2) + 160*pow(Om0, 5)*Or0*a + 64*pow(Om0, 4)*pow(Or0, 2))*pow(log(4 + Om0/(Ol*pow(a, 3))), 4))/((24576*pow(Ol, 6)*pow(a, 20) + 36864*pow(Ol, 5)*Om0*pow(a, 17) + 23040*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) + 7680*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) + 1440*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) + 144*Ol*pow(Om0, 5)*pow(a, 5) + 6*pow(Om0, 6)*pow(a, 2))*pow(log(4 + Om0/(Ol*pow(a, 3))), 8)) - (18144*pow(Ol, 3)*pow(Om0, 3) + 36288*pow(Ol, 2)*pow(Om0, 4) + 36288*pow(Ol, 2)*pow(Om0, 3)*Or0 + 18144*Ol*pow(Om0, 5) + 36288*Ol*pow(Om0, 4)*Or0 + 18144*Ol*pow(Om0, 3)*pow(Or0, 2) + (-4096*pow(Ol, 6) - 6144*pow(Ol, 5)*Om0 - 3840*pow(Ol, 4)*pow(Om0, 2) - 1280*pow(Ol, 3)*pow(Om0, 3) - 240*pow(Ol, 2)*pow(Om0, 4) - 24*Ol*pow(Om0, 5) - pow(Om0, 6))*pow(log(4 + Om0/Ol), 6) + (-24192*pow(Ol, 4)*pow(Om0, 2) - 33696*pow(Ol, 3)*pow(Om0, 3) - 69120*pow(Ol, 3)*pow(Om0, 2)*Or0 + 16848*pow(Ol, 2)*pow(Om0, 4) - 29376*pow(Ol, 2)*pow(Om0, 3)*Or0 - 44928*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2) + 26352*Ol*pow(Om0, 5) + 51408*Ol*pow(Om0, 4)*Or0 + 25056*Ol*pow(Om0, 3)*pow(Or0, 2))*log(4 + Om0/Ol) + (4096*pow(Ol, 6) + 4096*pow(Ol, 5)*Om0 + 512*pow(Ol, 4)*pow(Om0, 2) - 2048*pow(Ol, 4)*Om0*Or0 + 384*pow(Ol, 3)*pow(Om0, 3) - 512*pow(Ol, 3)*pow(Om0, 2)*Or0 + 496*pow(Ol, 2)*pow(Om0, 4) + 384*pow(Ol, 2)*pow(Om0, 3)*Or0 + 152*Ol*pow(Om0, 5) + 160*Ol*pow(Om0, 4)*Or0 + 14*pow(Om0, 6) + 16*pow(Om0, 5)*Or0)*pow(log(4 + Om0/Ol), 5) + (-15360*pow(Ol, 5)*Om0 - 49344*pow(Ol, 4)*pow(Om0, 2) + 9216*pow(Ol, 4)*Om0*Or0 - 69552*pow(Ol, 3)*pow(Om0, 3) - 71424*pow(Ol, 3)*pow(Om0, 2)*Or0 + 18432*pow(Ol, 3)*Om0*pow(Or0, 2) - 16728*pow(Ol, 2)*pow(Om0, 4) - 73440*pow(Ol, 2)*pow(Om0, 3)*Or0 - 55872*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2) + 17040*Ol*pow(Om0, 5) + 32760*Ol*pow(Om0, 4)*Or0 + 15696*Ol*pow(Om0, 3)*pow(Or0, 2) + 144*pow(Om0, 6) + 288*pow(Om0, 5)*Or0 + 144*pow(Om0, 4)*pow(Or0, 2))*pow(log(4 + Om0/Ol), 2) + (4096*pow(Ol, 6) - 37888*pow(Ol, 5)*Om0 - 4096*pow(Ol, 5)*Or0 - 79808*pow(Ol, 4)*pow(Om0, 2) - 22528*pow(Ol, 4)*Om0*Or0 - 64928*pow(Ol, 3)*pow(Om0, 3) - 43008*pow(Ol, 3)*pow(Om0, 2)*Or0 + 21504*pow(Ol, 3)*Om0*pow(Or0, 2) - 19528*pow(Ol, 2)*pow(Om0, 4) - 42880*pow(Ol, 2)*pow(Om0, 3)*Or0 - 23616*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2) + 5800*Ol*pow(Om0, 5) + 11864*Ol*pow(Om0, 4)*Or0 + 6048*Ol*pow(Om0, 3)*pow(Or0, 2) + 168*pow(Om0, 6) + 360*pow(Om0, 5)*Or0 + 192*pow(Om0, 4)*pow(Or0, 2))*pow(log(4 + Om0/Ol), 3) + (13312*pow(Ol, 6) - 1536*pow(Ol, 5)*Om0 - 14336*pow(Ol, 5)*Or0 - 19008*pow(Ol, 4)*pow(Om0, 2) - 22528*pow(Ol, 4)*Om0*Or0 + 1024*pow(Ol, 4)*pow(Or0, 2) - 11872*pow(Ol, 3)*pow(Om0, 3) - 11072*pow(Ol, 3)*pow(Om0, 2)*Or0 + 7168*pow(Ol, 3)*Om0*pow(Or0, 2) - 2148*pow(Ol, 2)*pow(Om0, 4) - 6560*pow(Ol, 2)*pow(Om0, 3)*Or0 - 3072*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2) + 1530*Ol*pow(Om0, 5) + 2636*Ol*pow(Om0, 4)*Or0 + 1216*Ol*pow(Om0, 3)*pow(Or0, 2) + 97*pow(Om0, 6) + 160*pow(Om0, 5)*Or0 + 64*pow(Om0, 4)*pow(Or0, 2))*pow(log(4 + Om0/Ol), 4))/((24576*pow(Ol, 6) + 36864*pow(Ol, 5)*Om0 + 23040*pow(Ol, 4)*pow(Om0, 2) + 7680*pow(Ol, 3)*pow(Om0, 3) + 1440*pow(Ol, 2)*pow(Om0, 4) + 144*Ol*pow(Om0, 5) + 6*pow(Om0, 6))*pow(log(4 + Om0/Ol), 8))) - 3*b*((-12*Ol*Om0*pow(a, 4) - 12*pow(Om0, 2)*a - 12*Om0*Or0 + (16*pow(Ol, 2)*pow(a, 7) + 8*Ol*Om0*pow(a, 4) + pow(Om0, 2)*a)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + (8*pow(Ol, 2)*pow(a, 7) - 8*Ol*Om0*pow(a, 4) - 8*Ol*Or0*pow(a, 3) - 7*pow(Om0, 2)*a - 8*Om0*Or0)*log(4 + Om0/(Ol*pow(a, 3))))/((48*pow(Ol, 2)*pow(a, 7) + 24*Ol*Om0*pow(a, 4) + 3*pow(Om0, 2)*a)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3)) - (-12*Ol*Om0 - 12*pow(Om0, 2) - 12*Om0*Or0 + (16*pow(Ol, 2) + 8*Ol*Om0 + pow(Om0, 2))*pow(log(4 + Om0/Ol), 2) + (8*pow(Ol, 2) - 8*Ol*Om0 - 8*Ol*Or0 - 7*pow(Om0, 2) - 8*Om0*Or0)*log(4 + Om0/Ol))/((48*pow(Ol, 2) + 24*Ol*Om0 + 3*pow(Om0, 2))*pow(log(4 + Om0/Ol), 3))));	
    		
    		
    	pvecback[pba->index_bg_rho_ds] = rho_fR * pba->Omega0_ds * pow(pba->H0,2);
        pvecback[pba->index_bg_p_ds]   = w_fR *pvecback[pba->index_bg_rho_ds];
    		
    		}if (l == 8){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = ArcTanh (Case 8) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = (1.0/2.0)*(pow(Ol, 2) - pow(3*Om0 + 4*Or0 - 4, 2))*(fR0 - 1)/pow(Ol, 2);

   			w_fR = pow(b, 2)*(262136250*pow(Ol, 13)*Om0*pow(a, 39)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 145071000*pow(Ol, 13)*Om0*pow(a, 39) + 103275000*pow(Ol, 13)*Or0*pow(a, 38)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 50220000*pow(Ol, 13)*Or0*pow(a, 38) + 841833000*pow(Ol, 12)*pow(Om0, 2)*pow(a, 36)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 501778800*pow(Ol, 12)*pow(Om0, 2)*pow(a, 36) + 771120000*pow(Ol, 12)*Om0*Or0*pow(a, 35)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 8113500*pow(Ol, 12)*Om0*Or0*pow(a, 35) - 12150000*pow(Ol, 12)*pow(Or0, 2)*pow(a, 34) + 961295850*pow(Ol, 11)*pow(Om0, 3)*pow(a, 33)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 873155700*pow(Ol, 11)*pow(Om0, 3)*pow(a, 33) + 1250977500*pow(Ol, 11)*pow(Om0, 2)*Or0*pow(a, 32)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 634176000*pow(Ol, 11)*pow(Om0, 2)*Or0*pow(a, 32) - 423144000*pow(Ol, 11)*Om0*pow(Or0, 2)*pow(a, 31) + 458823240*pow(Ol, 10)*pow(Om0, 4)*pow(a, 30)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 564167106*pow(Ol, 10)*pow(Om0, 4)*pow(a, 30) + 887119200*pow(Ol, 10)*pow(Om0, 3)*Or0*pow(a, 29)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 1505579820*pow(Ol, 10)*pow(Om0, 3)*Or0*pow(a, 29) + 400491000*pow(Ol, 10)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 28) - 250440*pow(Ol, 9)*pow(Om0, 5)*pow(a, 27)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 231318204*pow(Ol, 9)*pow(Om0, 5)*pow(a, 27) + 276999420*pow(Ol, 9)*pow(Om0, 4)*Or0*pow(a, 26)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 432230960*pow(Ol, 9)*pow(Om0, 4)*Or0*pow(a, 26) + 403476480*pow(Ol, 9)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 25) - 113325360*pow(Ol, 8)*pow(Om0, 6)*pow(a, 24)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 416500464*pow(Ol, 8)*pow(Om0, 6)*pow(a, 24) - 11216128*pow(Ol, 8)*pow(Om0, 5)*Or0*pow(a, 23)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 436418592*pow(Ol, 8)*pow(Om0, 5)*Or0*pow(a, 23) - 54464016*pow(Ol, 8)*pow(Om0, 4)*pow(Or0, 2)*pow(a, 22) - 63582732*pow(Ol, 7)*pow(Om0, 7)*pow(a, 21)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 182391156*pow(Ol, 7)*pow(Om0, 7)*pow(a, 21) - 41648104*pow(Ol, 7)*pow(Om0, 6)*Or0*pow(a, 20)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 312811472*pow(Ol, 7)*pow(Om0, 6)*Or0*pow(a, 20) - 122107136*pow(Ol, 7)*pow(Om0, 5)*pow(Or0, 2)*pow(a, 19) - 18390960*pow(Ol, 6)*pow(Om0, 8)*pow(a, 18)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 32119548*pow(Ol, 6)*pow(Om0, 8)*pow(a, 18) - 16129088*pow(Ol, 6)*pow(Om0, 7)*Or0*pow(a, 17)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 73810624*pow(Ol, 6)*pow(Om0, 7)*Or0*pow(a, 17) - 38275712*pow(Ol, 6)*pow(Om0, 6)*pow(Or0, 2)*pow(a, 16) - 3138870*pow(Ol, 5)*pow(Om0, 9)*pow(a, 15)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 573420*pow(Ol, 5)*pow(Om0, 9)*pow(a, 15) - 3220560*pow(Ol, 5)*pow(Om0, 8)*Or0*pow(a, 14)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 4981648*pow(Ol, 5)*pow(Om0, 8)*Or0*pow(a, 14) - 3965312*pow(Ol, 5)*pow(Om0, 7)*pow(Or0, 2)*pow(a, 13) - 302328*pow(Ol, 4)*pow(Om0, 10)*pow(a, 12)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 524880*pow(Ol, 4)*pow(Om0, 10)*pow(a, 12) - 358272*pow(Ol, 4)*pow(Om0, 9)*Or0*pow(a, 11)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 609380*pow(Ol, 4)*pow(Om0, 9)*Or0*pow(a, 11) + 108032*pow(Ol, 4)*pow(Om0, 8)*pow(Or0, 2)*pow(a, 10) - 11838*pow(Ol, 3)*pow(Om0, 11)*pow(a, 9)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 52752*pow(Ol, 3)*pow(Om0, 11)*pow(a, 9) - 18868*pow(Ol, 3)*pow(Om0, 10)*Or0*pow(a, 8)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 86864*pow(Ol, 3)*pow(Om0, 10)*Or0*pow(a, 8) + 34496*pow(Ol, 3)*pow(Om0, 9)*pow(Or0, 2)*pow(a, 7) + 360*pow(Ol, 2)*pow(Om0, 12)*pow(a, 6)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 618*pow(Ol, 2)*pow(Om0, 12)*pow(a, 6) - 32*pow(Ol, 2)*pow(Om0, 11)*Or0*pow(a, 5)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 1012*pow(Ol, 2)*pow(Om0, 11)*Or0*pow(a, 5) + 392*pow(Ol, 2)*pow(Om0, 10)*pow(Or0, 2)*pow(a, 4) + 36*Ol*pow(Om0, 13)*pow(a, 3)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 28*Ol*pow(Om0, 12)*Or0*pow(a, 2)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)))/(512578125*pow(Ol, 14)*pow(a, 42) + 1913625000*pow(Ol, 13)*Om0*pow(a, 39) + 3301003125*pow(Ol, 12)*pow(Om0, 2)*pow(a, 36) + 3487050000*pow(Ol, 11)*pow(Om0, 3)*pow(a, 33) + 2519960625*pow(Ol, 10)*pow(Om0, 4)*pow(a, 30) + 1317821400*pow(Ol, 9)*pow(Om0, 5)*pow(a, 27) + 514274985*pow(Ol, 8)*pow(Om0, 6)*pow(a, 24) + 152138976*pow(Ol, 7)*pow(Om0, 7)*pow(a, 21) + 34284999*pow(Ol, 6)*pow(Om0, 8)*pow(a, 18) + 5856984*pow(Ol, 5)*pow(Om0, 9)*pow(a, 15) + 746655*pow(Ol, 4)*pow(Om0, 10)*pow(a, 12) + 68880*pow(Ol, 3)*pow(Om0, 11)*pow(a, 9) + 4347*pow(Ol, 2)*pow(Om0, 12)*pow(a, 6) + 168*Ol*pow(Om0, 13)*pow(a, 3) + 3*pow(Om0, 14)) + b*(4230*pow(Ol, 5)*Om0*pow(a, 15) + 1800*pow(Ol, 5)*Or0*pow(a, 14) + 5892*pow(Ol, 4)*pow(Om0, 2)*pow(a, 12) + 6240*pow(Ol, 4)*Om0*Or0*pow(a, 11) + 1614*pow(Ol, 3)*pow(Om0, 3)*pow(a, 9) + 2372*pow(Ol, 3)*pow(Om0, 2)*Or0*pow(a, 8) - 36*pow(Ol, 2)*pow(Om0, 4)*pow(a, 6) + 128*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 5) - 36*Ol*pow(Om0, 5)*pow(a, 3) - 28*Ol*pow(Om0, 4)*Or0*pow(a, 2))/(10125*pow(Ol, 6)*pow(a, 18) + 16200*pow(Ol, 5)*Om0*pow(a, 15) + 10665*pow(Ol, 4)*pow(Om0, 2)*pow(a, 12) + 3696*pow(Ol, 3)*pow(Om0, 3)*pow(a, 9) + 711*pow(Ol, 2)*pow(Om0, 4)*pow(a, 6) + 72*Ol*pow(Om0, 5)*pow(a, 3) + 3*pow(Om0, 6)) - 1;
   			
   			// dark energy density
   			rho_fR = 1;	
   			
   	pvecback[pba->index_bg_rho_ds] = rho_fR * pba->Omega0_ds * pow(pba->H0,2);
        pvecback[pba->index_bg_p_ds]   = w_fR *pvecback[pba->index_bg_rho_ds];
        
    		}else {	// LCDM
    	double w_ds = pba->w_ds;
   	pvecback[pba->index_bg_rho_ds] = pba->Omega0_ds*pow(pba->H0,2)/pow(a_rel,3.*(1.+w_ds));
    
    	pvecback[pba->index_bg_p_ds] = w_ds*pvecback[pba->index_bg_rho_ds]; 
  	}

	p_tot   += pvecback[pba->index_bg_p_ds];
	rho_tot += pvecback[pba->index_bg_rho_ds]; 

	}
      if (pba->mg_type == fR) {

    double f = pvecback_B[pba->index_bi_f_ds];
    double f_prime = pvecback_B[pba->index_bi_f_prime_ds];
    pvecback[pba->index_bg_f_ds] = f; // value of df/dR
    pvecback[pba->index_bg_f_prime_ds] = f_prime; // value of the df/dR derivative wrt ln(a)

  }



  /* dcdm */
  if (pba->has_dcdm == _TRUE_) {
    /* Pass value of rho_dcdm to output */
    pvecback[pba->index_bg_rho_dcdm] = pvecback_B[pba->index_bi_rho_dcdm];
    rho_tot += pvecback[pba->index_bg_rho_dcdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_dcdm];
  }

  /* dr */
  if (pba->has_dr == _TRUE_) {
    /* Pass value of rho_dr to output */
    pvecback[pba->index_bg_rho_dr] = pvecback_B[pba->index_bi_rho_dr];
    rho_tot += pvecback[pba->index_bg_rho_dr];
    p_tot += (1./3.)*pvecback[pba->index_bg_rho_dr];
    rho_r += pvecback[pba->index_bg_rho_dr];
  }

  /* Scalar field */
  if (pba->has_scf == _TRUE_) {
    phi = pvecback_B[pba->index_bi_phi_scf];
    phi_prime = pvecback_B[pba->index_bi_phi_prime_scf];
    pvecback[pba->index_bg_phi_scf] = phi; // value of the scalar field phi
    pvecback[pba->index_bg_phi_prime_scf] = phi_prime; // value of the scalar field phi derivative wrt conformal time
    pvecback[pba->index_bg_V_scf] = V_scf(pba,phi); //V_scf(pba,phi); //write here potential as function of phi
    pvecback[pba->index_bg_dV_scf] = dV_scf(pba,phi); // dV_scf(pba,phi); //potential' as function of phi
    pvecback[pba->index_bg_ddV_scf] = ddV_scf(pba,phi); // ddV_scf(pba,phi); //potential'' as function of phi
    pvecback[pba->index_bg_rho_scf] = (phi_prime*phi_prime/(2*a*a) + V_scf(pba,phi))/3.; // energy of the scalar field. The field units are set automatically by setting the initial conditions
    pvecback[pba->index_bg_p_scf] =(phi_prime*phi_prime/(2*a*a) - V_scf(pba,phi))/3.; // pressure of the scalar field
    rho_tot += pvecback[pba->index_bg_rho_scf];
    p_tot += pvecback[pba->index_bg_p_scf];
    //divide relativistic & nonrelativistic (not very meaningful for oscillatory models)
    rho_r += 3.*pvecback[pba->index_bg_p_scf]; //field pressure contributes radiation
    rho_m += pvecback[pba->index_bg_rho_scf] - 3.* pvecback[pba->index_bg_p_scf]; //the rest contributes matter
    //printf(" a= %e, Omega_scf = %f, \n ",a_rel, pvecback[pba->index_bg_rho_scf]/rho_tot );
  }

  /* ncdm */
  if (pba->has_ncdm == _TRUE_) {

    /* Loop over species: */
    for(n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++){

      /* function returning background ncdm[n_ncdm] quantities (only
         those for which non-NULL pointers are passed) */
      class_call(background_ncdm_momenta(
                                         pba->q_ncdm_bg[n_ncdm],
                                         pba->w_ncdm_bg[n_ncdm],
                                         pba->q_size_ncdm_bg[n_ncdm],
                                         pba->M_ncdm[n_ncdm],
                                         pba->factor_ncdm[n_ncdm],
                                         1./a_rel-1.,
                                         NULL,
                                         &rho_ncdm,
                                         &p_ncdm,
                                         NULL,
                                         &pseudo_p_ncdm),
                 pba->error_message,
                 pba->error_message);

      pvecback[pba->index_bg_rho_ncdm1+n_ncdm] = rho_ncdm;
      rho_tot += rho_ncdm;
      pvecback[pba->index_bg_p_ncdm1+n_ncdm] = p_ncdm;
      p_tot += p_ncdm;
      pvecback[pba->index_bg_pseudo_p_ncdm1+n_ncdm] = pseudo_p_ncdm;

      /* (3 p_ncdm1) is the "relativistic" contribution to rho_ncdm1 */
      rho_r += 3.* p_ncdm;

      /* (rho_ncdm1 - 3 p_ncdm1) is the "non-relativistic" contribution
         to rho_ncdm1 */
      rho_m += rho_ncdm - 3.* p_ncdm;
    }
  }

  /* Lambda */
  
  /*
  // f(R) models: lcdm, HS(n=1)
double H0 = pba->H0;
double HSb,fR0,om0,Om0,Ol,Or0,w_fR,rho_fR,b,l,N,PI,G;
Om0 = pba->Omega0_b + pba->Omega0_cdm;
om0 = pba->Omega0_b + pba->Omega0_cdm;
Or0 = pba->Omega0_ur + pba->Omega0_g;
Ol = pba->Omega0_lambda;
l = pba->l;
fR0 = 1 - pba->fR0; 
//fR0 = -0.1;

  if (pba->has_lambda == _TRUE_) {
  	if (l == 0){		//lcdm
  		static int printed = 0;
  		if (!printed) {
      		printf("DEBUG: BACKGROUND MODEL = LCDM (Case 0) SELECTED\n");
     		printed = 1;          // non stampare più
   		}
   		
    		pvecback[pba->index_bg_rho_lambda] = pba->Omega0_lambda * pow(pba->H0,2);
    		rho_tot += pvecback[pba->index_bg_rho_lambda];
  		p_tot -= pvecback[pba->index_bg_rho_lambda];
  		
    		}
    	if (l > 0){
    		if (l == 1){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = HS1 (Case 1) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = (2 - 2*fR0)*pow(-3.0/2.0*Om0 - 2*Or0 + 2, 2)/pow(Ol, 2);

   			w_fR = (2.0/3.0)*pow(Ol, 2)*pow(a, 4)*pow(b, 2)*(122880*pow(Ol, 6)*Om0*pow(a, 20) + 49152*pow(Ol, 6)*Or0*pow(a, 19) + 239616*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17) + 53248*pow(Ol, 5)*Om0*Or0*pow(a, 16) - 8192*pow(Ol, 5)*pow(Or0, 2)*pow(a, 15) + 225792*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14) + 363008*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13) - 243712*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12) - 55344*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11) + 417664*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10) + 503296*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9) - 195084*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8) - 447296*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7) - 244352*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6) + 19944*Ol*pow(Om0, 6)*pow(a, 5) + 32960*Ol*pow(Om0, 5)*Or0*pow(a, 4) + 13328*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3) + 327*pow(Om0, 7)*pow(a, 2) + 520*pow(Om0, 6)*Or0*a + 196*pow(Om0, 5)*pow(Or0, 2))/pow(4*Ol*pow(a, 3) + Om0, 9) + (4.0/3.0)*Ol*pow(a, 2)*b*(72*pow(Ol, 2)*Om0*pow(a, 7) + 32*pow(Ol, 2)*Or0*pow(a, 6) + 63*Ol*pow(Om0, 2)*pow(a, 4) + 88*Ol*Om0*Or0*pow(a, 3) - 9*pow(Om0, 3)*a - 7*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 4) - 1
;
   			
   			// dark energy density
   			rho_fR =  exp(-3*pow(b, 2)*((-2621440*pow(Ol, 7)*Om0*pow(a, 21) - 786432*pow(Ol, 7)*Or0*pow(a, 20) - 4849664*pow(Ol, 6)*pow(Om0, 2)*pow(a, 18) - 1048576*pow(Ol, 6)*Om0*Or0*pow(a, 17) + 65536*pow(Ol, 6)*pow(Or0, 2)*pow(a, 16) - 4030464*pow(Ol, 5)*pow(Om0, 3)*pow(a, 15) - 2768896*pow(Ol, 5)*pow(Om0, 2)*Or0*pow(a, 14) + 1441792*pow(Ol, 5)*Om0*pow(Or0, 2)*pow(a, 13) - 964352*pow(Ol, 4)*pow(Om0, 4)*pow(a, 12) - 2801664*pow(Ol, 4)*pow(Om0, 3)*Or0*pow(a, 11) - 1966080*pow(Ol, 4)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 10) + 639488*pow(Ol, 3)*pow(Om0, 5)*pow(a, 9) + 1307648*pow(Ol, 3)*pow(Om0, 4)*Or0*pow(a, 8) + 630784*pow(Ol, 3)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 7) + 9024*pow(Ol, 2)*pow(Om0, 6)*pow(a, 6) + 26624*pow(Ol, 2)*pow(Om0, 5)*Or0*pow(a, 5) + 12544*pow(Ol, 2)*pow(Om0, 4)*pow(Or0, 2)*pow(a, 4) - 352*Ol*pow(Om0, 7)*pow(a, 3) - 11*pow(Om0, 8))/(25165824*pow(Ol, 8)*pow(a, 24) + 50331648*pow(Ol, 7)*Om0*pow(a, 21) + 44040192*pow(Ol, 6)*pow(Om0, 2)*pow(a, 18) + 22020096*pow(Ol, 5)*pow(Om0, 3)*pow(a, 15) + 6881280*pow(Ol, 4)*pow(Om0, 4)*pow(a, 12) + 1376256*pow(Ol, 3)*pow(Om0, 5)*pow(a, 9) + 172032*pow(Ol, 2)*pow(Om0, 6)*pow(a, 6) + 12288*Ol*pow(Om0, 7)*pow(a, 3) + 384*pow(Om0, 8)) - (-2621440*pow(Ol, 7)*Om0 - 786432*pow(Ol, 7)*Or0 - 4849664*pow(Ol, 6)*pow(Om0, 2) - 1048576*pow(Ol, 6)*Om0*Or0 + 65536*pow(Ol, 6)*pow(Or0, 2) - 4030464*pow(Ol, 5)*pow(Om0, 3) - 2768896*pow(Ol, 5)*pow(Om0, 2)*Or0 + 1441792*pow(Ol, 5)*Om0*pow(Or0, 2) - 964352*pow(Ol, 4)*pow(Om0, 4) - 2801664*pow(Ol, 4)*pow(Om0, 3)*Or0 - 1966080*pow(Ol, 4)*pow(Om0, 2)*pow(Or0, 2) + 639488*pow(Ol, 3)*pow(Om0, 5) + 1307648*pow(Ol, 3)*pow(Om0, 4)*Or0 + 630784*pow(Ol, 3)*pow(Om0, 3)*pow(Or0, 2) + 9024*pow(Ol, 2)*pow(Om0, 6) + 26624*pow(Ol, 2)*pow(Om0, 5)*Or0 + 12544*pow(Ol, 2)*pow(Om0, 4)*pow(Or0, 2) - 352*Ol*pow(Om0, 7) - 11*pow(Om0, 8))/(25165824*pow(Ol, 8) + 50331648*pow(Ol, 7)*Om0 + 44040192*pow(Ol, 6)*pow(Om0, 2) + 22020096*pow(Ol, 5)*pow(Om0, 3) + 6881280*pow(Ol, 4)*pow(Om0, 4) + 1376256*pow(Ol, 3)*pow(Om0, 5) + 172032*pow(Ol, 2)*pow(Om0, 6) + 12288*Ol*pow(Om0, 7) + 384*pow(Om0, 8))) - 3*b*((-192*pow(Ol, 2)*Om0*pow(a, 6) - 64*pow(Ol, 2)*Or0*pow(a, 5) - 132*Ol*pow(Om0, 2)*pow(a, 3) - 112*Ol*Om0*Or0*pow(a, 2) - 3*pow(Om0, 3))/(1536*pow(Ol, 3)*pow(a, 9) + 1152*pow(Ol, 2)*Om0*pow(a, 6) + 288*Ol*pow(Om0, 2)*pow(a, 3) + 24*pow(Om0, 3)) - (-192*pow(Ol, 2)*Om0 - 64*pow(Ol, 2)*Or0 - 132*Ol*pow(Om0, 2) - 112*Ol*Om0*Or0 - 3*pow(Om0, 3))/(1536*pow(Ol, 3) + 1152*pow(Ol, 2)*Om0 + 288*Ol*pow(Om0, 2) + 24*pow(Om0, 3))))
;	
    		}
    		
    		if (l == 2){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = HS4 (Case 2) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = (1.0/2.0)*pow(2, 1.0/4.0)*pow(pow(H0, 10)*(fR0 - 1)*pow(3*Om0 + 4*Or0 - 4, 5), 1.0/4.0)/pow(pow(H0, 2)*Ol, 5.0/4.0);

   			w_fR = (4.0/3.0)*pow(Ol, 8)*pow(a, 22)*pow(b, 8)*(3932160*pow(Ol, 6)*Om0*pow(a, 20) + 589824*pow(Ol, 6)*Or0*pow(a, 19) + 12201984*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17) - 2289664*pow(Ol, 5)*Om0*Or0*pow(a, 16) - 65536*pow(Ol, 5)*pow(Or0, 2)*pow(a, 15) + 1301760*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14) + 54906880*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13) - 10833920*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12) - 46567680*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11) - 11008000*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10) + 61829120*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9) - 25864560*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8) - 95334400*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7) - 67425280*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6) + 16514640*Ol*pow(Om0, 6)*pow(a, 5) + 31165456*Ol*pow(Om0, 5)*Or0*pow(a, 4) + 14671360*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3) + 48021*pow(Om0, 7)*pow(a, 2) + 93104*pow(Om0, 6)*Or0*a + 45056*pow(Om0, 5)*pow(Or0, 2))/pow(4*Ol*pow(a, 3) + Om0, 15) + (4.0/3.0)*pow(Ol, 4)*pow(a, 11)*pow(b, 4)*(576*pow(Ol, 2)*Om0*pow(a, 7) + 128*pow(Ol, 2)*Or0*pow(a, 6) + 288*Ol*pow(Om0, 2)*pow(a, 4) + 784*Ol*Om0*Or0*pow(a, 3) - 369*pow(Om0, 3)*a - 352*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 7) - 1;
   			
   			// dark energy density
   			rho_fR = exp(-3*pow(b, 8)*((-171798691840*pow(Ol, 13)*Om0*pow(a, 39) - 19327352832*pow(Ol, 13)*Or0*pow(a, 38) - 545729282048*pow(Ol, 12)*pow(Om0, 2)*pow(a, 36) + 16642998272*pow(Ol, 12)*Om0*Or0*pow(a, 35) + 1073741824*pow(Ol, 12)*pow(Or0, 2)*pow(a, 34) - 564687536128*pow(Ol, 11)*pow(Om0, 3)*pow(a, 33) - 705112834048*pow(Ol, 11)*pow(Om0, 2)*Or0*pow(a, 32) + 129922760704*pow(Ol, 11)*Om0*pow(Or0, 2)*pow(a, 31) + 120420564992*pow(Ol, 10)*pow(Om0, 4)*pow(a, 30) - 322927853568*pow(Ol, 10)*pow(Om0, 3)*Or0*pow(a, 29) - 506940358656*pow(Ol, 10)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 28) + 286218256384*pow(Ol, 9)*pow(Om0, 5)*pow(a, 27) + 634652721152*pow(Ol, 9)*pow(Om0, 4)*Or0*pow(a, 26) + 311116693504*pow(Ol, 9)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 25) - 12924092416*pow(Ol, 8)*pow(Om0, 6)*pow(a, 24) + 2122317824*pow(Ol, 8)*pow(Om0, 5)*Or0*pow(a, 23) + 1073741824*pow(Ol, 8)*pow(Om0, 4)*pow(Or0, 2)*pow(a, 22) - 3992322048*pow(Ol, 7)*pow(Om0, 7)*pow(a, 21) - 873320448*pow(Ol, 6)*pow(Om0, 8)*pow(a, 18) - 145553408*pow(Ol, 5)*pow(Om0, 9)*pow(a, 15) - 18194176*pow(Ol, 4)*pow(Om0, 10)*pow(a, 12) - 1654016*pow(Ol, 3)*pow(Om0, 11)*pow(a, 9) - 103376*pow(Ol, 2)*pow(Om0, 12)*pow(a, 6) - 3976*Ol*pow(Om0, 13)*pow(a, 3) - 71*pow(Om0, 14))/(105553116266496*pow(Ol, 14)*pow(a, 42) + 369435906932736*pow(Ol, 13)*Om0*pow(a, 39) + 600333348765696*pow(Ol, 12)*pow(Om0, 2)*pow(a, 36) + 600333348765696*pow(Ol, 11)*pow(Om0, 3)*pow(a, 33) + 412729177276416*pow(Ol, 10)*pow(Om0, 4)*pow(a, 30) + 206364588638208*pow(Ol, 9)*pow(Om0, 5)*pow(a, 27) + 77386720739328*pow(Ol, 8)*pow(Om0, 6)*pow(a, 24) + 22110491639808*pow(Ol, 7)*pow(Om0, 7)*pow(a, 21) + 4836670046208*pow(Ol, 6)*pow(Om0, 8)*pow(a, 18) + 806111674368*pow(Ol, 5)*pow(Om0, 9)*pow(a, 15) + 100763959296*pow(Ol, 4)*pow(Om0, 10)*pow(a, 12) + 9160359936*pow(Ol, 3)*pow(Om0, 11)*pow(a, 9) + 572522496*pow(Ol, 2)*pow(Om0, 12)*pow(a, 6) + 22020096*Ol*pow(Om0, 13)*pow(a, 3) + 393216*pow(Om0, 14)) - (-171798691840*pow(Ol, 13)*Om0 - 19327352832*pow(Ol, 13)*Or0 - 545729282048*pow(Ol, 12)*pow(Om0, 2) + 16642998272*pow(Ol, 12)*Om0*Or0 + 1073741824*pow(Ol, 12)*pow(Or0, 2) - 564687536128*pow(Ol, 11)*pow(Om0, 3) - 705112834048*pow(Ol, 11)*pow(Om0, 2)*Or0 + 129922760704*pow(Ol, 11)*Om0*pow(Or0, 2) + 120420564992*pow(Ol, 10)*pow(Om0, 4) - 322927853568*pow(Ol, 10)*pow(Om0, 3)*Or0 - 506940358656*pow(Ol, 10)*pow(Om0, 2)*pow(Or0, 2) + 286218256384*pow(Ol, 9)*pow(Om0, 5) + 634652721152*pow(Ol, 9)*pow(Om0, 4)*Or0 + 311116693504*pow(Ol, 9)*pow(Om0, 3)*pow(Or0, 2) - 12924092416*pow(Ol, 8)*pow(Om0, 6) + 2122317824*pow(Ol, 8)*pow(Om0, 5)*Or0 + 1073741824*pow(Ol, 8)*pow(Om0, 4)*pow(Or0, 2) - 3992322048*pow(Ol, 7)*pow(Om0, 7) - 873320448*pow(Ol, 6)*pow(Om0, 8) - 145553408*pow(Ol, 5)*pow(Om0, 9) - 18194176*pow(Ol, 4)*pow(Om0, 10) - 1654016*pow(Ol, 3)*pow(Om0, 11) - 103376*pow(Ol, 2)*pow(Om0, 12) - 3976*Ol*pow(Om0, 13) - 71*pow(Om0, 14))/(105553116266496*pow(Ol, 14) + 369435906932736*pow(Ol, 13)*Om0 + 600333348765696*pow(Ol, 12)*pow(Om0, 2) + 600333348765696*pow(Ol, 11)*pow(Om0, 3) + 412729177276416*pow(Ol, 10)*pow(Om0, 4) + 206364588638208*pow(Ol, 9)*pow(Om0, 5) + 77386720739328*pow(Ol, 8)*pow(Om0, 6) + 22110491639808*pow(Ol, 7)*pow(Om0, 7) + 4836670046208*pow(Ol, 6)*pow(Om0, 8) + 806111674368*pow(Ol, 5)*pow(Om0, 9) + 100763959296*pow(Ol, 4)*pow(Om0, 10) + 9160359936*pow(Ol, 3)*pow(Om0, 11) + 572522496*pow(Ol, 2)*pow(Om0, 12) + 22020096*Ol*pow(Om0, 13) + 393216*pow(Om0, 14))) - 3*pow(b, 4)*((-49152*pow(Ol, 5)*Om0*pow(a, 15) - 8192*pow(Ol, 5)*Or0*pow(a, 14) - 43008*pow(Ol, 4)*pow(Om0, 2)*pow(a, 12) - 32768*pow(Ol, 4)*Om0*Or0*pow(a, 11) - 3840*pow(Ol, 3)*pow(Om0, 3)*pow(a, 9) - 720*pow(Ol, 2)*pow(Om0, 4)*pow(a, 6) - 72*Ol*pow(Om0, 5)*pow(a, 3) - 3*pow(Om0, 6))/(3145728*pow(Ol, 6)*pow(a, 18) + 4718592*pow(Ol, 5)*Om0*pow(a, 15) + 2949120*pow(Ol, 4)*pow(Om0, 2)*pow(a, 12) + 983040*pow(Ol, 3)*pow(Om0, 3)*pow(a, 9) + 184320*pow(Ol, 2)*pow(Om0, 4)*pow(a, 6) + 18432*Ol*pow(Om0, 5)*pow(a, 3) + 768*pow(Om0, 6)) - (-49152*pow(Ol, 5)*Om0 - 8192*pow(Ol, 5)*Or0 - 43008*pow(Ol, 4)*pow(Om0, 2) - 32768*pow(Ol, 4)*Om0*Or0 - 3840*pow(Ol, 3)*pow(Om0, 3) - 720*pow(Ol, 2)*pow(Om0, 4) - 72*Ol*pow(Om0, 5) - 3*pow(Om0, 6))/(3145728*pow(Ol, 6) + 4718592*pow(Ol, 5)*Om0 + 2949120*pow(Ol, 4)*pow(Om0, 2) + 983040*pow(Ol, 3)*pow(Om0, 3) + 184320*pow(Ol, 2)*pow(Om0, 4) + 18432*Ol*pow(Om0, 5) + 768*pow(Om0, 6))));	
    		}
    		
    		if (l == 3){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = Star1 (Case 3) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = (1.0/2.0)*sqrt(pow(H0, 6)*(fR0 - 1)*pow(3*Om0 + 4*Or0 - 4, 3))/pow(pow(H0, 2)*Ol, 3.0/2.0);

   			w_fR = (2.0/3.0)*pow(Ol, 4)*pow(a, 10)*pow(b, 4)*(884736*pow(Ol, 6)*Om0*pow(a, 20) + 229376*pow(Ol, 6)*Or0*pow(a, 19) + 2016768*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17) - 57344*pow(Ol, 5)*Om0*Or0*pow(a, 16) - 32768*pow(Ol, 5)*pow(Or0, 2)*pow(a, 15) + 1936896*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14) + 5477888*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13) - 2050048*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12) - 2307840*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11) + 3363200*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10) + 6776320*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9) - 2905440*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8) - 7610144*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7) - 4563584*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6) + 641490*Ol*pow(Om0, 6)*pow(a, 5) + 1157960*Ol*pow(Om0, 5)*Or0*pow(a, 4) + 519200*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3) + 4845*pow(Om0, 7)*pow(a, 2) + 8844*pow(Om0, 6)*Or0*a + 4000*pow(Om0, 5)*pow(Or0, 2))/pow(4*Ol*pow(a, 3) + Om0, 11) + (2.0/3.0)*pow(Ol, 2)*pow(a, 5)*pow(b, 2)*(384*pow(Ol, 2)*Om0*pow(a, 7) + 128*pow(Ol, 2)*Or0*pow(a, 6) + 300*Ol*pow(Om0, 2)*pow(a, 4) + 496*Ol*Om0*Or0*pow(a, 3) - 111*pow(Om0, 3)*a - 100*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 5) - 1;
   			
   			// dark energy density
   			rho_fR = exp(-3*pow(b, 4)*((-37748736*pow(Ol, 9)*Om0*pow(a, 27) - 7340032*pow(Ol, 9)*Or0*pow(a, 26) - 85491712*pow(Ol, 8)*pow(Om0, 2)*pow(a, 24) - 5767168*pow(Ol, 8)*Om0*Or0*pow(a, 23) + 524288*pow(Ol, 8)*pow(Or0, 2)*pow(a, 22) - 84541440*pow(Ol, 7)*pow(Om0, 3)*pow(a, 21) - 73433088*pow(Ol, 7)*pow(Om0, 2)*Or0*pow(a, 20) + 24117248*pow(Ol, 7)*Om0*pow(Or0, 2)*pow(a, 19) - 12369920*pow(Ol, 6)*pow(Om0, 4)*pow(a, 18) - 61358080*pow(Ol, 6)*pow(Om0, 3)*Or0*pow(a, 17) - 53772288*pow(Ol, 6)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 16) + 21082112*pow(Ol, 5)*pow(Om0, 5)*pow(a, 15) + 44582912*pow(Ol, 5)*pow(Om0, 4)*Or0*pow(a, 14) + 21708800*pow(Ol, 5)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 13) - 169600*pow(Ol, 4)*pow(Om0, 6)*pow(a, 12) + 411648*pow(Ol, 4)*pow(Om0, 5)*Or0*pow(a, 11) + 204800*pow(Ol, 4)*pow(Om0, 4)*pow(Or0, 2)*pow(a, 10) - 53760*pow(Ol, 3)*pow(Om0, 7)*pow(a, 9) - 5040*pow(Ol, 2)*pow(Om0, 8)*pow(a, 6) - 280*Ol*pow(Om0, 9)*pow(a, 3) - 7*pow(Om0, 10))/(805306368*pow(Ol, 10)*pow(a, 30) + 2013265920*pow(Ol, 9)*Om0*pow(a, 27) + 2264924160*pow(Ol, 8)*pow(Om0, 2)*pow(a, 24) + 1509949440*pow(Ol, 7)*pow(Om0, 3)*pow(a, 21) + 660602880*pow(Ol, 6)*pow(Om0, 4)*pow(a, 18) + 198180864*pow(Ol, 5)*pow(Om0, 5)*pow(a, 15) + 41287680*pow(Ol, 4)*pow(Om0, 6)*pow(a, 12) + 5898240*pow(Ol, 3)*pow(Om0, 7)*pow(a, 9) + 552960*pow(Ol, 2)*pow(Om0, 8)*pow(a, 6) + 30720*Ol*pow(Om0, 9)*pow(a, 3) + 768*pow(Om0, 10)) - (-37748736*pow(Ol, 9)*Om0 - 7340032*pow(Ol, 9)*Or0 - 85491712*pow(Ol, 8)*pow(Om0, 2) - 5767168*pow(Ol, 8)*Om0*Or0 + 524288*pow(Ol, 8)*pow(Or0, 2) - 84541440*pow(Ol, 7)*pow(Om0, 3) - 73433088*pow(Ol, 7)*pow(Om0, 2)*Or0 + 24117248*pow(Ol, 7)*Om0*pow(Or0, 2) - 12369920*pow(Ol, 6)*pow(Om0, 4) - 61358080*pow(Ol, 6)*pow(Om0, 3)*Or0 - 53772288*pow(Ol, 6)*pow(Om0, 2)*pow(Or0, 2) + 21082112*pow(Ol, 5)*pow(Om0, 5) + 44582912*pow(Ol, 5)*pow(Om0, 4)*Or0 + 21708800*pow(Ol, 5)*pow(Om0, 3)*pow(Or0, 2) - 169600*pow(Ol, 4)*pow(Om0, 6) + 411648*pow(Ol, 4)*pow(Om0, 5)*Or0 + 204800*pow(Ol, 4)*pow(Om0, 4)*pow(Or0, 2) - 53760*pow(Ol, 3)*pow(Om0, 7) - 5040*pow(Ol, 2)*pow(Om0, 8) - 280*Ol*pow(Om0, 9) - 7*pow(Om0, 10))/(805306368*pow(Ol, 10) + 2013265920*pow(Ol, 9)*Om0 + 2264924160*pow(Ol, 8)*pow(Om0, 2) + 1509949440*pow(Ol, 7)*pow(Om0, 3) + 660602880*pow(Ol, 6)*pow(Om0, 4) + 198180864*pow(Ol, 5)*pow(Om0, 5) + 41287680*pow(Ol, 4)*pow(Om0, 6) + 5898240*pow(Ol, 3)*pow(Om0, 7) + 552960*pow(Ol, 2)*pow(Om0, 8) + 30720*Ol*pow(Om0, 9) + 768*pow(Om0, 10))) - 3*pow(b, 2)*((-512*pow(Ol, 3)*Om0*pow(a, 9) - 128*pow(Ol, 3)*Or0*pow(a, 8) - 392*pow(Ol, 2)*pow(Om0, 2)*pow(a, 6) - 320*pow(Ol, 2)*Om0*Or0*pow(a, 5) - 16*Ol*pow(Om0, 3)*pow(a, 3) - pow(Om0, 4))/(6144*pow(Ol, 4)*pow(a, 12) + 6144*pow(Ol, 3)*Om0*pow(a, 9) + 2304*pow(Ol, 2)*pow(Om0, 2)*pow(a, 6) + 384*Ol*pow(Om0, 3)*pow(a, 3) + 24*pow(Om0, 4)) - (-512*pow(Ol, 3)*Om0 - 128*pow(Ol, 3)*Or0 - 392*pow(Ol, 2)*pow(Om0, 2) - 320*pow(Ol, 2)*Om0*Or0 - 16*Ol*pow(Om0, 3) - pow(Om0, 4))/(6144*pow(Ol, 4) + 6144*pow(Ol, 3)*Om0 + 2304*pow(Ol, 2)*pow(Om0, 2) + 384*Ol*pow(Om0, 3) + 24*pow(Om0, 4))))
;	
    		}
    		
    		if (l == 4){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = Star2 (Case 4) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0) DA MATHEMATICA
   			b = (1.0/2.0)*pow(2, 1.0/4.0)*pow(pow(H0, 10)*(fR0 - 1)*pow(3*Om0 + 4*Or0 - 4, 5), 1.0/4.0)/pow(pow(H0, 2)*Ol, 5.0/4.0);

   			w_fR = -4*pow(Ol, 6)*pow(a, 17)*pow(b, 6)*(768*pow(Ol, 2)*Om0*pow(a, 7) + 128*pow(Ol, 2)*Or0*pow(a, 6) + 132*Ol*pow(Om0, 2)*pow(a, 4) + 1072*Ol*Om0*Or0*pow(a, 3) - 771*pow(Om0, 3)*a - 748*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 9) + (4.0/3.0)*pow(Ol, 4)*pow(a, 11)*pow(b, 4)*(576*pow(Ol, 2)*Om0*pow(a, 7) + 128*pow(Ol, 2)*Or0*pow(a, 6) + 288*Ol*pow(Om0, 2)*pow(a, 4) + 784*Ol*Om0*Or0*pow(a, 3) - 369*pow(Om0, 3)*a - 352*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 7) - 1;
   			
   			// dark energy density
   			rho_fR =exp(-3*pow(b, 6)*((393216*pow(Ol, 7)*Om0*pow(a, 21) + 49152*pow(Ol, 7)*Or0*pow(a, 20) + 377856*pow(Ol, 6)*pow(Om0, 2)*pow(a, 18) + 270336*pow(Ol, 6)*Om0*Or0*pow(a, 17) + 57344*pow(Ol, 5)*pow(Om0, 3)*pow(a, 15) + 17920*pow(Ol, 4)*pow(Om0, 4)*pow(a, 12) + 3584*pow(Ol, 3)*pow(Om0, 5)*pow(a, 9) + 448*pow(Ol, 2)*pow(Om0, 6)*pow(a, 6) + 32*Ol*pow(Om0, 7)*pow(a, 3) + pow(Om0, 8))/(100663296*pow(Ol, 8)*pow(a, 24) + 201326592*pow(Ol, 7)*Om0*pow(a, 21) + 176160768*pow(Ol, 6)*pow(Om0, 2)*pow(a, 18) + 88080384*pow(Ol, 5)*pow(Om0, 3)*pow(a, 15) + 27525120*pow(Ol, 4)*pow(Om0, 4)*pow(a, 12) + 5505024*pow(Ol, 3)*pow(Om0, 5)*pow(a, 9) + 688128*pow(Ol, 2)*pow(Om0, 6)*pow(a, 6) + 49152*Ol*pow(Om0, 7)*pow(a, 3) + 1536*pow(Om0, 8)) - (393216*pow(Ol, 7)*Om0 + 49152*pow(Ol, 7)*Or0 + 377856*pow(Ol, 6)*pow(Om0, 2) + 270336*pow(Ol, 6)*Om0*Or0 + 57344*pow(Ol, 5)*pow(Om0, 3) + 17920*pow(Ol, 4)*pow(Om0, 4) + 3584*pow(Ol, 3)*pow(Om0, 5) + 448*pow(Ol, 2)*pow(Om0, 6) + 32*Ol*pow(Om0, 7) + pow(Om0, 8))/(100663296*pow(Ol, 8) + 201326592*pow(Ol, 7)*Om0 + 176160768*pow(Ol, 6)*pow(Om0, 2) + 88080384*pow(Ol, 5)*pow(Om0, 3) + 27525120*pow(Ol, 4)*pow(Om0, 4) + 5505024*pow(Ol, 3)*pow(Om0, 5) + 688128*pow(Ol, 2)*pow(Om0, 6) + 49152*Ol*pow(Om0, 7) + 1536*pow(Om0, 8))) - 3*pow(b, 4)*((-49152*pow(Ol, 5)*Om0*pow(a, 15) - 8192*pow(Ol, 5)*Or0*pow(a, 14) - 43008*pow(Ol, 4)*pow(Om0, 2)*pow(a, 12) - 32768*pow(Ol, 4)*Om0*Or0*pow(a, 11) - 3840*pow(Ol, 3)*pow(Om0, 3)*pow(a, 9) - 720*pow(Ol, 2)*pow(Om0, 4)*pow(a, 6) - 72*Ol*pow(Om0, 5)*pow(a, 3) - 3*pow(Om0, 6))/(3145728*pow(Ol, 6)*pow(a, 18) + 4718592*pow(Ol, 5)*Om0*pow(a, 15) + 2949120*pow(Ol, 4)*pow(Om0, 2)*pow(a, 12) + 983040*pow(Ol, 3)*pow(Om0, 3)*pow(a, 9) + 184320*pow(Ol, 2)*pow(Om0, 4)*pow(a, 6) + 18432*Ol*pow(Om0, 5)*pow(a, 3) + 768*pow(Om0, 6)) - (-49152*pow(Ol, 5)*Om0 - 8192*pow(Ol, 5)*Or0 - 43008*pow(Ol, 4)*pow(Om0, 2) - 32768*pow(Ol, 4)*Om0*Or0 - 3840*pow(Ol, 3)*pow(Om0, 3) - 720*pow(Ol, 2)*pow(Om0, 4) - 72*Ol*pow(Om0, 5) - 3*pow(Om0, 6))/(3145728*pow(Ol, 6) + 4718592*pow(Ol, 5)*Om0 + 2949120*pow(Ol, 4)*pow(Om0, 2) + 983040*pow(Ol, 3)*pow(Om0, 3) + 184320*pow(Ol, 2)*pow(Om0, 4) + 18432*Ol*pow(Om0, 5) + 768*pow(Om0, 6))));	
    		}
    		
    		if (l == 5){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = Exp (Case 5) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = -1.0/2.0*(fR0 - 1)*pow(3*Om0 + 4*Or0 - 4, 2)/pow(Ol, 2); // GIÀ CAMBIATO

   			w_fR = (1.0/3.0)*pow(Ol, 2)*pow(a, 4)*pow(b, 2)*(344064*pow(Ol, 6)*Om0*pow(a, 20) + 131072*pow(Ol, 6)*Or0*pow(a, 19) + 654336*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17) + 266240*pow(Ol, 5)*Om0*Or0*pow(a, 16) - 16384*pow(Ol, 5)*pow(Or0, 2)*pow(a, 15) + 536832*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14) + 839680*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13) - 487424*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12) - 104160*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11) + 859392*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10) + 1006592*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9) - 395640*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8) - 896128*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7) - 488704*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6) + 38412*Ol*pow(Om0, 6)*pow(a, 5) + 64816*Ol*pow(Om0, 5)*Or0*pow(a, 4) + 26656*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3) + 543*pow(Om0, 7)*pow(a, 2) + 940*pow(Om0, 6)*Or0*a + 392*pow(Om0, 5)*pow(Or0, 2))/pow(4*Ol*pow(a, 3) + Om0, 9) + (4.0/3.0)*Ol*pow(a, 2)*b*(72*pow(Ol, 2)*Om0*pow(a, 7) + 32*pow(Ol, 2)*Or0*pow(a, 6) + 63*Ol*pow(Om0, 2)*pow(a, 4) + 88*Ol*Om0*Or0*pow(a, 3) - 9*pow(Om0, 3)*a - 7*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 4) - 1;
   			
   			// dark energy density
   			rho_fR = exp(-3*pow(b, 2)*((-3670016*pow(Ol, 7)*Om0*pow(a, 21) - 1048576*pow(Ol, 7)*Or0*pow(a, 20) - 6701056*pow(Ol, 6)*pow(Om0, 2)*pow(a, 18) - 1966080*pow(Ol, 6)*Om0*Or0*pow(a, 17) + 65536*pow(Ol, 6)*pow(Or0, 2)*pow(a, 16) - 5259264*pow(Ol, 5)*pow(Om0, 3)*pow(a, 15) - 3522560*pow(Ol, 5)*pow(Om0, 2)*Or0*pow(a, 14) + 1441792*pow(Ol, 5)*Om0*pow(Or0, 2)*pow(a, 13) - 1365760*pow(Ol, 4)*pow(Om0, 4)*pow(a, 12) - 3063808*pow(Ol, 4)*pow(Om0, 3)*Or0*pow(a, 11) - 1966080*pow(Ol, 4)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 10) + 570880*pow(Ol, 3)*pow(Om0, 5)*pow(a, 9) + 1265664*pow(Ol, 3)*pow(Om0, 4)*Or0*pow(a, 8) + 630784*pow(Ol, 3)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 7) + 3072*pow(Ol, 2)*pow(Om0, 6)*pow(a, 6) + 24064*pow(Ol, 2)*pow(Om0, 5)*Or0*pow(a, 5) + 12544*pow(Ol, 2)*pow(Om0, 4)*pow(Or0, 2)*pow(a, 4) - 608*Ol*pow(Om0, 7)*pow(a, 3) - 19*pow(Om0, 8))/(25165824*pow(Ol, 8)*pow(a, 24) + 50331648*pow(Ol, 7)*Om0*pow(a, 21) + 44040192*pow(Ol, 6)*pow(Om0, 2)*pow(a, 18) + 22020096*pow(Ol, 5)*pow(Om0, 3)*pow(a, 15) + 6881280*pow(Ol, 4)*pow(Om0, 4)*pow(a, 12) + 1376256*pow(Ol, 3)*pow(Om0, 5)*pow(a, 9) + 172032*pow(Ol, 2)*pow(Om0, 6)*pow(a, 6) + 12288*Ol*pow(Om0, 7)*pow(a, 3) + 384*pow(Om0, 8)) - (-3670016*pow(Ol, 7)*Om0 - 1048576*pow(Ol, 7)*Or0 - 6701056*pow(Ol, 6)*pow(Om0, 2) - 1966080*pow(Ol, 6)*Om0*Or0 + 65536*pow(Ol, 6)*pow(Or0, 2) - 5259264*pow(Ol, 5)*pow(Om0, 3) - 3522560*pow(Ol, 5)*pow(Om0, 2)*Or0 + 1441792*pow(Ol, 5)*Om0*pow(Or0, 2) - 1365760*pow(Ol, 4)*pow(Om0, 4) - 3063808*pow(Ol, 4)*pow(Om0, 3)*Or0 - 1966080*pow(Ol, 4)*pow(Om0, 2)*pow(Or0, 2) + 570880*pow(Ol, 3)*pow(Om0, 5) + 1265664*pow(Ol, 3)*pow(Om0, 4)*Or0 + 630784*pow(Ol, 3)*pow(Om0, 3)*pow(Or0, 2) + 3072*pow(Ol, 2)*pow(Om0, 6) + 24064*pow(Ol, 2)*pow(Om0, 5)*Or0 + 12544*pow(Ol, 2)*pow(Om0, 4)*pow(Or0, 2) - 608*Ol*pow(Om0, 7) - 19*pow(Om0, 8))/(25165824*pow(Ol, 8) + 50331648*pow(Ol, 7)*Om0 + 44040192*pow(Ol, 6)*pow(Om0, 2) + 22020096*pow(Ol, 5)*pow(Om0, 3) + 6881280*pow(Ol, 4)*pow(Om0, 4) + 1376256*pow(Ol, 3)*pow(Om0, 5) + 172032*pow(Ol, 2)*pow(Om0, 6) + 12288*Ol*pow(Om0, 7) + 384*pow(Om0, 8))) - 3*b*((-192*pow(Ol, 2)*Om0*pow(a, 6) - 64*pow(Ol, 2)*Or0*pow(a, 5) - 132*Ol*pow(Om0, 2)*pow(a, 3) - 112*Ol*Om0*Or0*pow(a, 2) - 3*pow(Om0, 3))/(1536*pow(Ol, 3)*pow(a, 9) + 1152*pow(Ol, 2)*Om0*pow(a, 6) + 288*Ol*pow(Om0, 2)*pow(a, 3) + 24*pow(Om0, 3)) - (-192*pow(Ol, 2)*Om0 - 64*pow(Ol, 2)*Or0 - 132*Ol*pow(Om0, 2) - 112*Ol*Om0*Or0 - 3*pow(Om0, 3))/(1536*pow(Ol, 3) + 1152*pow(Ol, 2)*Om0 + 288*Ol*pow(Om0, 2) + 24*pow(Om0, 3))));	
    		}
    		
    		if (l == 6){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = Tsjk (Case 6) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = -1.0/2.0*sqrt(pow(H0, 6)*(fR0 - 1)*pow(3*Om0 + 4*Or0 - 4, 3))/pow(pow(H0, 2)*Ol, 3.0/2.0);

   			w_fR = (2.0/3.0)*pow(Ol, 4)*pow(a, 10)*pow(b, 3)*(32*pow(Ol, 2)*pow(a, 7) - 32*Ol*Om0*pow(a, 4) - 16*Ol*Or0*pow(a, 3) - 37*pow(Om0, 2)*a - 40*Om0*Or0)*(384*pow(Ol, 2)*Om0*pow(a, 7) + 128*pow(Ol, 2)*Or0*pow(a, 6) + 300*Ol*pow(Om0, 2)*pow(a, 4) + 496*Ol*Om0*Or0*pow(a, 3) - 111*pow(Om0, 3)*a - 100*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 9) + (2.0/3.0)*pow(Ol, 2)*pow(a, 5)*pow(b, 2)*(384*pow(Ol, 2)*Om0*pow(a, 7) + 128*pow(Ol, 2)*Or0*pow(a, 6) + 300*Ol*pow(Om0, 2)*pow(a, 4) + 496*Ol*Om0*Or0*pow(a, 3) - 111*pow(Om0, 3)*a - 100*pow(Om0, 2)*Or0)/pow(4*Ol*pow(a, 3) + Om0, 5) - 1
;
   			
   			// dark energy density
   			rho_fR = 1;	
    		}
    		
    		if (l == 7){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = Log (Case 7) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = (1.0/2.0)*(fR0 - 1)*(3*Om0 + 4*Or0 - 4)*pow(log((-3*Om0 - 4*Or0 + 4)/Ol), 2)/Ol;

   			w_fR = pow(b, 2)*(147456*pow(Ol, 6)*Om0*pow(a, 20)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 344064*pow(Ol, 6)*Om0*pow(a, 20)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 122880*pow(Ol, 6)*Om0*pow(a, 20)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 114688*pow(Ol, 6)*Or0*pow(a, 19)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 32768*pow(Ol, 6)*Or0*pow(a, 19)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 6144*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + 235008*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 664320*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 192768*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 152064*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 28672*pow(Ol, 5)*Om0*Or0*pow(a, 16)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + 215040*pow(Ol, 5)*Om0*Or0*pow(a, 16)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 200704*pow(Ol, 5)*Om0*Or0*pow(a, 16)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) - 159744*pow(Ol, 5)*Om0*Or0*pow(a, 16)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 16384*pow(Ol, 5)*pow(Or0, 2)*pow(a, 15)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) - 15360*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + 101952*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 575808*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 357312*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 17280*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 580608*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*log(4 + Om0/(Ol*pow(a, 3))) - 1024*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + 88320*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 601088*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 1310208*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 1465344*pow(Ol, 4)*pow(Om0, 2)*Or0*pow(a, 13)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 152576*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) - 466944*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) - 405504*pow(Ol, 4)*Om0*pow(Or0, 2)*pow(a, 12)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 14016*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) - 144*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 105264*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) - 398472*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 1181952*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 1143072*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*log(4 + Om0/(Ol*pow(a, 3))) + 217728*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11) - 12032*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + 123968*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 876416*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 1301184*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 155520*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 1669248*pow(Ol, 3)*pow(Om0, 3)*Or0*pow(a, 10)*log(4 + Om0/(Ol*pow(a, 3))) + 111104*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 779520*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 1790208*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 1423872*pow(Ol, 3)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 9)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 3792*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) - 50112*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) - 245472*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) - 707844*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 890568*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 258552*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*log(4 + Om0/(Ol*pow(a, 3))) + 435456*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8) - 4160*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) - 99024*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) - 526208*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) - 1553520*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 2379456*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 1378944*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7)*log(4 + Om0/(Ol*pow(a, 3))) + 435456*pow(Ol, 2)*pow(Om0, 4)*Or0*pow(a, 7) - 47488*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) - 271296*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) - 822528*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) - 1444608*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 1088640*pow(Ol, 2)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 6)*log(4 + Om0/(Ol*pow(a, 3))) - 348*Ol*pow(Om0, 6)*pow(a, 5)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) - 513*Ol*pow(Om0, 6)*pow(a, 5)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 11832*Ol*pow(Om0, 6)*pow(a, 5)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 63876*Ol*pow(Om0, 6)*pow(a, 5)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 192888*Ol*pow(Om0, 6)*pow(a, 5)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 303912*Ol*pow(Om0, 6)*pow(a, 5)*log(4 + Om0/(Ol*pow(a, 3))) + 217728*Ol*pow(Om0, 6)*pow(a, 5) - 448*Ol*pow(Om0, 5)*Or0*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) - 2724*Ol*pow(Om0, 5)*Or0*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 14000*Ol*pow(Om0, 5)*Or0*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 110796*Ol*pow(Om0, 5)*Or0*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 346248*Ol*pow(Om0, 5)*Or0*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 576072*Ol*pow(Om0, 5)*Or0*pow(a, 4)*log(4 + Om0/(Ol*pow(a, 3))) + 435456*Ol*pow(Om0, 5)*Or0*pow(a, 4) - 1952*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 2640*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 47448*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 153792*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 272160*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3)*log(4 + Om0/(Ol*pow(a, 3))) + 217728*Ol*pow(Om0, 4)*pow(Or0, 2)*pow(a, 3) - 3*pow(Om0, 7)*pow(a, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + 63*pow(Om0, 7)*pow(a, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 582*pow(Om0, 7)*pow(a, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 1260*pow(Om0, 7)*pow(a, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 1296*pow(Om0, 7)*pow(a, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 8*pow(Om0, 6)*Or0*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) - 8*pow(Om0, 6)*Or0*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 780*pow(Om0, 6)*Or0*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 2556*pow(Om0, 6)*Or0*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 2592*pow(Om0, 6)*Or0*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 64*pow(Om0, 5)*pow(Or0, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + 192*pow(Om0, 5)*pow(Or0, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 1296*pow(Om0, 5)*pow(Or0, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + 1296*pow(Om0, 5)*pow(Or0, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2))/(49152*pow(Ol, 7)*pow(a, 23)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 86016*pow(Ol, 6)*Om0*pow(a, 20)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 64512*pow(Ol, 5)*pow(Om0, 2)*pow(a, 17)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 26880*pow(Ol, 4)*pow(Om0, 3)*pow(a, 14)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 6720*pow(Ol, 3)*pow(Om0, 4)*pow(a, 11)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 1008*pow(Ol, 2)*pow(Om0, 5)*pow(a, 8)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 84*Ol*pow(Om0, 6)*pow(a, 5)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9) + 3*pow(Om0, 7)*pow(a, 2)*pow(log(4 + Om0/(Ol*pow(a, 3))), 9)) + b*(192*pow(Ol, 2)*Om0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 192*pow(Ol, 2)*Om0*pow(a, 7)*log(4 + Om0/(Ol*pow(a, 3))) + 128*pow(Ol, 2)*Or0*pow(a, 6)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 168*Ol*pow(Om0, 2)*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 204*Ol*pow(Om0, 2)*pow(a, 4)*log(4 + Om0/(Ol*pow(a, 3))) - 108*Ol*pow(Om0, 2)*pow(a, 4) + 208*Ol*Om0*Or0*pow(a, 3)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + 288*Ol*Om0*Or0*pow(a, 3)*log(4 + Om0/(Ol*pow(a, 3))) + 3*pow(Om0, 3)*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 42*pow(Om0, 3)*a*log(4 + Om0/(Ol*pow(a, 3))) - 108*pow(Om0, 3)*a + 8*pow(Om0, 2)*Or0*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) - 36*pow(Om0, 2)*Or0*log(4 + Om0/(Ol*pow(a, 3))) - 108*pow(Om0, 2)*Or0)/(192*pow(Ol, 3)*pow(a, 10)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 144*pow(Ol, 2)*Om0*pow(a, 7)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 36*Ol*pow(Om0, 2)*pow(a, 4)*pow(log(4 + Om0/(Ol*pow(a, 3))), 4) + 3*pow(Om0, 3)*a*pow(log(4 + Om0/(Ol*pow(a, 3))), 4)) - 1;
   			
   			// dark energy density
   			rho_fR = exp(-3*pow(b, 2)*((18144*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) + 36288*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) + 36288*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 7) + 18144*Ol*pow(Om0, 5)*pow(a, 5) + 36288*Ol*pow(Om0, 4)*Or0*pow(a, 4) + 18144*Ol*pow(Om0, 3)*pow(Or0, 2)*pow(a, 3) + (-4096*pow(Ol, 6)*pow(a, 20) - 6144*pow(Ol, 5)*Om0*pow(a, 17) - 3840*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) - 1280*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) - 240*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) - 24*Ol*pow(Om0, 5)*pow(a, 5) - pow(Om0, 6)*pow(a, 2))*pow(log(4 + Om0/(Ol*pow(a, 3))), 6) + (-24192*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) - 33696*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) - 69120*pow(Ol, 3)*pow(Om0, 2)*Or0*pow(a, 10) + 16848*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) - 29376*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 7) - 44928*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 6) + 26352*Ol*pow(Om0, 5)*pow(a, 5) + 51408*Ol*pow(Om0, 4)*Or0*pow(a, 4) + 25056*Ol*pow(Om0, 3)*pow(Or0, 2)*pow(a, 3))*log(4 + Om0/(Ol*pow(a, 3))) + (4096*pow(Ol, 6)*pow(a, 20) + 4096*pow(Ol, 5)*Om0*pow(a, 17) + 512*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) - 2048*pow(Ol, 4)*Om0*Or0*pow(a, 13) + 384*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) - 512*pow(Ol, 3)*pow(Om0, 2)*Or0*pow(a, 10) + 496*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) + 384*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 7) + 152*Ol*pow(Om0, 5)*pow(a, 5) + 160*Ol*pow(Om0, 4)*Or0*pow(a, 4) + 14*pow(Om0, 6)*pow(a, 2) + 16*pow(Om0, 5)*Or0*a)*pow(log(4 + Om0/(Ol*pow(a, 3))), 5) + (-15360*pow(Ol, 5)*Om0*pow(a, 17) - 49344*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) + 9216*pow(Ol, 4)*Om0*Or0*pow(a, 13) - 69552*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) - 71424*pow(Ol, 3)*pow(Om0, 2)*Or0*pow(a, 10) + 18432*pow(Ol, 3)*Om0*pow(Or0, 2)*pow(a, 9) - 16728*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) - 73440*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 7) - 55872*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 6) + 17040*Ol*pow(Om0, 5)*pow(a, 5) + 32760*Ol*pow(Om0, 4)*Or0*pow(a, 4) + 15696*Ol*pow(Om0, 3)*pow(Or0, 2)*pow(a, 3) + 144*pow(Om0, 6)*pow(a, 2) + 288*pow(Om0, 5)*Or0*a + 144*pow(Om0, 4)*pow(Or0, 2))*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + (4096*pow(Ol, 6)*pow(a, 20) - 37888*pow(Ol, 5)*Om0*pow(a, 17) - 4096*pow(Ol, 5)*Or0*pow(a, 16) - 79808*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) - 22528*pow(Ol, 4)*Om0*Or0*pow(a, 13) - 64928*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) - 43008*pow(Ol, 3)*pow(Om0, 2)*Or0*pow(a, 10) + 21504*pow(Ol, 3)*Om0*pow(Or0, 2)*pow(a, 9) - 19528*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) - 42880*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 7) - 23616*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 6) + 5800*Ol*pow(Om0, 5)*pow(a, 5) + 11864*Ol*pow(Om0, 4)*Or0*pow(a, 4) + 6048*Ol*pow(Om0, 3)*pow(Or0, 2)*pow(a, 3) + 168*pow(Om0, 6)*pow(a, 2) + 360*pow(Om0, 5)*Or0*a + 192*pow(Om0, 4)*pow(Or0, 2))*pow(log(4 + Om0/(Ol*pow(a, 3))), 3) + (13312*pow(Ol, 6)*pow(a, 20) - 1536*pow(Ol, 5)*Om0*pow(a, 17) - 14336*pow(Ol, 5)*Or0*pow(a, 16) - 19008*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) - 22528*pow(Ol, 4)*Om0*Or0*pow(a, 13) + 1024*pow(Ol, 4)*pow(Or0, 2)*pow(a, 12) - 11872*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) - 11072*pow(Ol, 3)*pow(Om0, 2)*Or0*pow(a, 10) + 7168*pow(Ol, 3)*Om0*pow(Or0, 2)*pow(a, 9) - 2148*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) - 6560*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 7) - 3072*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 6) + 1530*Ol*pow(Om0, 5)*pow(a, 5) + 2636*Ol*pow(Om0, 4)*Or0*pow(a, 4) + 1216*Ol*pow(Om0, 3)*pow(Or0, 2)*pow(a, 3) + 97*pow(Om0, 6)*pow(a, 2) + 160*pow(Om0, 5)*Or0*a + 64*pow(Om0, 4)*pow(Or0, 2))*pow(log(4 + Om0/(Ol*pow(a, 3))), 4))/((24576*pow(Ol, 6)*pow(a, 20) + 36864*pow(Ol, 5)*Om0*pow(a, 17) + 23040*pow(Ol, 4)*pow(Om0, 2)*pow(a, 14) + 7680*pow(Ol, 3)*pow(Om0, 3)*pow(a, 11) + 1440*pow(Ol, 2)*pow(Om0, 4)*pow(a, 8) + 144*Ol*pow(Om0, 5)*pow(a, 5) + 6*pow(Om0, 6)*pow(a, 2))*pow(log(4 + Om0/(Ol*pow(a, 3))), 8)) - (18144*pow(Ol, 3)*pow(Om0, 3) + 36288*pow(Ol, 2)*pow(Om0, 4) + 36288*pow(Ol, 2)*pow(Om0, 3)*Or0 + 18144*Ol*pow(Om0, 5) + 36288*Ol*pow(Om0, 4)*Or0 + 18144*Ol*pow(Om0, 3)*pow(Or0, 2) + (-4096*pow(Ol, 6) - 6144*pow(Ol, 5)*Om0 - 3840*pow(Ol, 4)*pow(Om0, 2) - 1280*pow(Ol, 3)*pow(Om0, 3) - 240*pow(Ol, 2)*pow(Om0, 4) - 24*Ol*pow(Om0, 5) - pow(Om0, 6))*pow(log(4 + Om0/Ol), 6) + (-24192*pow(Ol, 4)*pow(Om0, 2) - 33696*pow(Ol, 3)*pow(Om0, 3) - 69120*pow(Ol, 3)*pow(Om0, 2)*Or0 + 16848*pow(Ol, 2)*pow(Om0, 4) - 29376*pow(Ol, 2)*pow(Om0, 3)*Or0 - 44928*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2) + 26352*Ol*pow(Om0, 5) + 51408*Ol*pow(Om0, 4)*Or0 + 25056*Ol*pow(Om0, 3)*pow(Or0, 2))*log(4 + Om0/Ol) + (4096*pow(Ol, 6) + 4096*pow(Ol, 5)*Om0 + 512*pow(Ol, 4)*pow(Om0, 2) - 2048*pow(Ol, 4)*Om0*Or0 + 384*pow(Ol, 3)*pow(Om0, 3) - 512*pow(Ol, 3)*pow(Om0, 2)*Or0 + 496*pow(Ol, 2)*pow(Om0, 4) + 384*pow(Ol, 2)*pow(Om0, 3)*Or0 + 152*Ol*pow(Om0, 5) + 160*Ol*pow(Om0, 4)*Or0 + 14*pow(Om0, 6) + 16*pow(Om0, 5)*Or0)*pow(log(4 + Om0/Ol), 5) + (-15360*pow(Ol, 5)*Om0 - 49344*pow(Ol, 4)*pow(Om0, 2) + 9216*pow(Ol, 4)*Om0*Or0 - 69552*pow(Ol, 3)*pow(Om0, 3) - 71424*pow(Ol, 3)*pow(Om0, 2)*Or0 + 18432*pow(Ol, 3)*Om0*pow(Or0, 2) - 16728*pow(Ol, 2)*pow(Om0, 4) - 73440*pow(Ol, 2)*pow(Om0, 3)*Or0 - 55872*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2) + 17040*Ol*pow(Om0, 5) + 32760*Ol*pow(Om0, 4)*Or0 + 15696*Ol*pow(Om0, 3)*pow(Or0, 2) + 144*pow(Om0, 6) + 288*pow(Om0, 5)*Or0 + 144*pow(Om0, 4)*pow(Or0, 2))*pow(log(4 + Om0/Ol), 2) + (4096*pow(Ol, 6) - 37888*pow(Ol, 5)*Om0 - 4096*pow(Ol, 5)*Or0 - 79808*pow(Ol, 4)*pow(Om0, 2) - 22528*pow(Ol, 4)*Om0*Or0 - 64928*pow(Ol, 3)*pow(Om0, 3) - 43008*pow(Ol, 3)*pow(Om0, 2)*Or0 + 21504*pow(Ol, 3)*Om0*pow(Or0, 2) - 19528*pow(Ol, 2)*pow(Om0, 4) - 42880*pow(Ol, 2)*pow(Om0, 3)*Or0 - 23616*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2) + 5800*Ol*pow(Om0, 5) + 11864*Ol*pow(Om0, 4)*Or0 + 6048*Ol*pow(Om0, 3)*pow(Or0, 2) + 168*pow(Om0, 6) + 360*pow(Om0, 5)*Or0 + 192*pow(Om0, 4)*pow(Or0, 2))*pow(log(4 + Om0/Ol), 3) + (13312*pow(Ol, 6) - 1536*pow(Ol, 5)*Om0 - 14336*pow(Ol, 5)*Or0 - 19008*pow(Ol, 4)*pow(Om0, 2) - 22528*pow(Ol, 4)*Om0*Or0 + 1024*pow(Ol, 4)*pow(Or0, 2) - 11872*pow(Ol, 3)*pow(Om0, 3) - 11072*pow(Ol, 3)*pow(Om0, 2)*Or0 + 7168*pow(Ol, 3)*Om0*pow(Or0, 2) - 2148*pow(Ol, 2)*pow(Om0, 4) - 6560*pow(Ol, 2)*pow(Om0, 3)*Or0 - 3072*pow(Ol, 2)*pow(Om0, 2)*pow(Or0, 2) + 1530*Ol*pow(Om0, 5) + 2636*Ol*pow(Om0, 4)*Or0 + 1216*Ol*pow(Om0, 3)*pow(Or0, 2) + 97*pow(Om0, 6) + 160*pow(Om0, 5)*Or0 + 64*pow(Om0, 4)*pow(Or0, 2))*pow(log(4 + Om0/Ol), 4))/((24576*pow(Ol, 6) + 36864*pow(Ol, 5)*Om0 + 23040*pow(Ol, 4)*pow(Om0, 2) + 7680*pow(Ol, 3)*pow(Om0, 3) + 1440*pow(Ol, 2)*pow(Om0, 4) + 144*Ol*pow(Om0, 5) + 6*pow(Om0, 6))*pow(log(4 + Om0/Ol), 8))) - 3*b*((-12*Ol*Om0*pow(a, 4) - 12*pow(Om0, 2)*a - 12*Om0*Or0 + (16*pow(Ol, 2)*pow(a, 7) + 8*Ol*Om0*pow(a, 4) + pow(Om0, 2)*a)*pow(log(4 + Om0/(Ol*pow(a, 3))), 2) + (8*pow(Ol, 2)*pow(a, 7) - 8*Ol*Om0*pow(a, 4) - 8*Ol*Or0*pow(a, 3) - 7*pow(Om0, 2)*a - 8*Om0*Or0)*log(4 + Om0/(Ol*pow(a, 3))))/((48*pow(Ol, 2)*pow(a, 7) + 24*Ol*Om0*pow(a, 4) + 3*pow(Om0, 2)*a)*pow(log(4 + Om0/(Ol*pow(a, 3))), 3)) - (-12*Ol*Om0 - 12*pow(Om0, 2) - 12*Om0*Or0 + (16*pow(Ol, 2) + 8*Ol*Om0 + pow(Om0, 2))*pow(log(4 + Om0/Ol), 2) + (8*pow(Ol, 2) - 8*Ol*Om0 - 8*Ol*Or0 - 7*pow(Om0, 2) - 8*Om0*Or0)*log(4 + Om0/Ol))/((48*pow(Ol, 2) + 24*Ol*Om0 + 3*pow(Om0, 2))*pow(log(4 + Om0/Ol), 3))));	
    		}
    		
    		if (l == 8){
    			static int printed = 0;
    			if (!printed) {
      			printf("DEBUG: BACKGROUND MODEL = ArcTanh (Case 8) SELECTED\n");
     			printed = 1;          // non stampare più
   			}
   			// deviation parameter b(fR0)
   			b = (1.0/2.0)*(pow(Ol, 2) - pow(3*Om0 + 4*Or0 - 4, 2))*(fR0 - 1)/pow(Ol, 2);

   			w_fR = pow(b, 2)*(262136250*pow(Ol, 13)*Om0*pow(a, 39)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 145071000*pow(Ol, 13)*Om0*pow(a, 39) + 103275000*pow(Ol, 13)*Or0*pow(a, 38)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 50220000*pow(Ol, 13)*Or0*pow(a, 38) + 841833000*pow(Ol, 12)*pow(Om0, 2)*pow(a, 36)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 501778800*pow(Ol, 12)*pow(Om0, 2)*pow(a, 36) + 771120000*pow(Ol, 12)*Om0*Or0*pow(a, 35)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 8113500*pow(Ol, 12)*Om0*Or0*pow(a, 35) - 12150000*pow(Ol, 12)*pow(Or0, 2)*pow(a, 34) + 961295850*pow(Ol, 11)*pow(Om0, 3)*pow(a, 33)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 873155700*pow(Ol, 11)*pow(Om0, 3)*pow(a, 33) + 1250977500*pow(Ol, 11)*pow(Om0, 2)*Or0*pow(a, 32)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 634176000*pow(Ol, 11)*pow(Om0, 2)*Or0*pow(a, 32) - 423144000*pow(Ol, 11)*Om0*pow(Or0, 2)*pow(a, 31) + 458823240*pow(Ol, 10)*pow(Om0, 4)*pow(a, 30)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 564167106*pow(Ol, 10)*pow(Om0, 4)*pow(a, 30) + 887119200*pow(Ol, 10)*pow(Om0, 3)*Or0*pow(a, 29)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 1505579820*pow(Ol, 10)*pow(Om0, 3)*Or0*pow(a, 29) + 400491000*pow(Ol, 10)*pow(Om0, 2)*pow(Or0, 2)*pow(a, 28) - 250440*pow(Ol, 9)*pow(Om0, 5)*pow(a, 27)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 231318204*pow(Ol, 9)*pow(Om0, 5)*pow(a, 27) + 276999420*pow(Ol, 9)*pow(Om0, 4)*Or0*pow(a, 26)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 432230960*pow(Ol, 9)*pow(Om0, 4)*Or0*pow(a, 26) + 403476480*pow(Ol, 9)*pow(Om0, 3)*pow(Or0, 2)*pow(a, 25) - 113325360*pow(Ol, 8)*pow(Om0, 6)*pow(a, 24)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 416500464*pow(Ol, 8)*pow(Om0, 6)*pow(a, 24) - 11216128*pow(Ol, 8)*pow(Om0, 5)*Or0*pow(a, 23)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 436418592*pow(Ol, 8)*pow(Om0, 5)*Or0*pow(a, 23) - 54464016*pow(Ol, 8)*pow(Om0, 4)*pow(Or0, 2)*pow(a, 22) - 63582732*pow(Ol, 7)*pow(Om0, 7)*pow(a, 21)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 182391156*pow(Ol, 7)*pow(Om0, 7)*pow(a, 21) - 41648104*pow(Ol, 7)*pow(Om0, 6)*Or0*pow(a, 20)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 312811472*pow(Ol, 7)*pow(Om0, 6)*Or0*pow(a, 20) - 122107136*pow(Ol, 7)*pow(Om0, 5)*pow(Or0, 2)*pow(a, 19) - 18390960*pow(Ol, 6)*pow(Om0, 8)*pow(a, 18)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 32119548*pow(Ol, 6)*pow(Om0, 8)*pow(a, 18) - 16129088*pow(Ol, 6)*pow(Om0, 7)*Or0*pow(a, 17)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 73810624*pow(Ol, 6)*pow(Om0, 7)*Or0*pow(a, 17) - 38275712*pow(Ol, 6)*pow(Om0, 6)*pow(Or0, 2)*pow(a, 16) - 3138870*pow(Ol, 5)*pow(Om0, 9)*pow(a, 15)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 573420*pow(Ol, 5)*pow(Om0, 9)*pow(a, 15) - 3220560*pow(Ol, 5)*pow(Om0, 8)*Or0*pow(a, 14)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) - 4981648*pow(Ol, 5)*pow(Om0, 8)*Or0*pow(a, 14) - 3965312*pow(Ol, 5)*pow(Om0, 7)*pow(Or0, 2)*pow(a, 13) - 302328*pow(Ol, 4)*pow(Om0, 10)*pow(a, 12)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 524880*pow(Ol, 4)*pow(Om0, 10)*pow(a, 12) - 358272*pow(Ol, 4)*pow(Om0, 9)*Or0*pow(a, 11)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 609380*pow(Ol, 4)*pow(Om0, 9)*Or0*pow(a, 11) + 108032*pow(Ol, 4)*pow(Om0, 8)*pow(Or0, 2)*pow(a, 10) - 11838*pow(Ol, 3)*pow(Om0, 11)*pow(a, 9)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 52752*pow(Ol, 3)*pow(Om0, 11)*pow(a, 9) - 18868*pow(Ol, 3)*pow(Om0, 10)*Or0*pow(a, 8)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 86864*pow(Ol, 3)*pow(Om0, 10)*Or0*pow(a, 8) + 34496*pow(Ol, 3)*pow(Om0, 9)*pow(Or0, 2)*pow(a, 7) + 360*pow(Ol, 2)*pow(Om0, 12)*pow(a, 6)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 618*pow(Ol, 2)*pow(Om0, 12)*pow(a, 6) - 32*pow(Ol, 2)*pow(Om0, 11)*Or0*pow(a, 5)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 1012*pow(Ol, 2)*pow(Om0, 11)*Or0*pow(a, 5) + 392*pow(Ol, 2)*pow(Om0, 10)*pow(Or0, 2)*pow(a, 4) + 36*Ol*pow(Om0, 13)*pow(a, 3)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)) + 28*Ol*pow(Om0, 12)*Or0*pow(a, 2)*atanh(Ol*pow(a, 3)/(4*Ol*pow(a, 3) + Om0)))/(512578125*pow(Ol, 14)*pow(a, 42) + 1913625000*pow(Ol, 13)*Om0*pow(a, 39) + 3301003125*pow(Ol, 12)*pow(Om0, 2)*pow(a, 36) + 3487050000*pow(Ol, 11)*pow(Om0, 3)*pow(a, 33) + 2519960625*pow(Ol, 10)*pow(Om0, 4)*pow(a, 30) + 1317821400*pow(Ol, 9)*pow(Om0, 5)*pow(a, 27) + 514274985*pow(Ol, 8)*pow(Om0, 6)*pow(a, 24) + 152138976*pow(Ol, 7)*pow(Om0, 7)*pow(a, 21) + 34284999*pow(Ol, 6)*pow(Om0, 8)*pow(a, 18) + 5856984*pow(Ol, 5)*pow(Om0, 9)*pow(a, 15) + 746655*pow(Ol, 4)*pow(Om0, 10)*pow(a, 12) + 68880*pow(Ol, 3)*pow(Om0, 11)*pow(a, 9) + 4347*pow(Ol, 2)*pow(Om0, 12)*pow(a, 6) + 168*Ol*pow(Om0, 13)*pow(a, 3) + 3*pow(Om0, 14)) + b*(4230*pow(Ol, 5)*Om0*pow(a, 15) + 1800*pow(Ol, 5)*Or0*pow(a, 14) + 5892*pow(Ol, 4)*pow(Om0, 2)*pow(a, 12) + 6240*pow(Ol, 4)*Om0*Or0*pow(a, 11) + 1614*pow(Ol, 3)*pow(Om0, 3)*pow(a, 9) + 2372*pow(Ol, 3)*pow(Om0, 2)*Or0*pow(a, 8) - 36*pow(Ol, 2)*pow(Om0, 4)*pow(a, 6) + 128*pow(Ol, 2)*pow(Om0, 3)*Or0*pow(a, 5) - 36*Ol*pow(Om0, 5)*pow(a, 3) - 28*Ol*pow(Om0, 4)*Or0*pow(a, 2))/(10125*pow(Ol, 6)*pow(a, 18) + 16200*pow(Ol, 5)*Om0*pow(a, 15) + 10665*pow(Ol, 4)*pow(Om0, 2)*pow(a, 12) + 3696*pow(Ol, 3)*pow(Om0, 3)*pow(a, 9) + 711*pow(Ol, 2)*pow(Om0, 4)*pow(a, 6) + 72*Ol*pow(Om0, 5)*pow(a, 3) + 3*pow(Om0, 6)) - 1;
   			
   			// dark energy density
   			rho_fR = 1;	
    		}
    		pvecback[pba->index_bg_rho_lambda] =rho_fR * pba->Omega0_lambda * pow(pba->H0,2);
    		rho_tot += pvecback[pba->index_bg_rho_lambda];
    		p_tot += w_fR * pvecback[pba->index_bg_rho_lambda];
    	}
    	*/
  }

  /* fluid with w=w0+wa(1-a/a0) and constant cs2 */
  if (pba->has_fld == _TRUE_) {
    pvecback[pba->index_bg_rho_fld] = pba->Omega0_fld * pow(pba->H0,2)
      / pow(a_rel,3.*(1.+pba->w0_fld+pba->wa_fld))
      * exp(3.*pba->wa_fld*(a_rel-1.));
    rho_tot += pvecback[pba->index_bg_rho_fld];
    p_tot += (pba->w0_fld+pba->wa_fld*(1.-a_rel)) * pvecback[pba->index_bg_rho_fld];
  }

  /* relativistic neutrinos (and all relativistic relics) */
  if (pba->has_ur == _TRUE_) {
    pvecback[pba->index_bg_rho_ur] = pba->Omega0_ur * pow(pba->H0,2) / pow(a_rel,4);
    rho_tot += pvecback[pba->index_bg_rho_ur];
    p_tot += (1./3.) * pvecback[pba->index_bg_rho_ur];
    rho_r += pvecback[pba->index_bg_rho_ur];
  }

  /** - compute expansion rate H from Friedmann equation: this is the
      only place where the Friedmann equation is assumed. Remember
      that densities are all expressed in units of \f$ [3c^2/8\pi G] \f$, ie
      \f$ \rho_{class} = [8 \pi G \rho_{physical} / 3 c^2]\f$ */
 // pvecback[pba->index_bg_H] = sqrt(rho_tot-pba->K/a/a);
  

//double a = y[pba->index_bi_a];
//H2 = rho_tot-pba->K/a/a;	
//H_prime = - (3./2.) * (rho_tot + p_tot) * a + pba->K/a;

// Hubble parameter and prime derivative
pvecback[pba->index_bg_H] = sqrt( rho_tot-pba->K/a/a );
pvecback[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a + pba->K/a;


  /** - compute relativistic density to total density ratio */
  pvecback[pba->index_bg_Omega_r] = rho_r / rho_tot;



 /**DS - compute  epsilon_H and Ricci  with respect to conformal time */
  pvecback[pba->index_bg_epsilon_H] = - pvecback[pba->index_bg_H_prime]/(a*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]);
  pvecback[pba->index_bg_Ricci] = 12.*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]*(1.-0.5*pvecback[pba->index_bg_epsilon_H]);

  


  
 if (pba->has_ds == _TRUE_) {
    pvecback[pba->index_bg_Omega_ds] = rho_ds / rho_tot;

    double w_ds = pba->w_ds;
    //w_ds = -1.;
    double Omega_m = 1.-pvecback[pba->index_bg_Omega_ds];
    double Omega_ds = pvecback[pba->index_bg_Omega_ds];
    double w_m_Omega_m = -w_ds*pvecback[pba->index_bg_Omega_ds]
                         +(2./3.)*pvecback[pba->index_bg_epsilon_H] -1.;
    double w_m = w_m_Omega_m/(1.-pvecback[pba->index_bg_Omega_ds]);
    pvecback[pba->index_bg_w_m] = w_m;

   
    double Omega_ds_prime =
    Omega_ds*
    (
    2.*pvecback[pba->index_bg_epsilon_H]
    -3.*(1.+w_ds)
     );
  
    double w_c = 0.;
    double w_c_Omega_c = 0.;
    double w_b = 0.;
    double w_b_Omega_b = 0.;
    double w_g = 1./3.;
    double w_g_Omega_g = w_g*(pvecback[pba->index_bg_rho_g]+pvecback[pba->index_bg_rho_ur])/
             (pvecback[pba->index_bg_H]
             *pvecback[pba->index_bg_H]);

    double one_plus_w_g_Omega_g =
      (1.+w_g)
      *(pvecback[pba->index_bg_rho_g]+pvecback[pba->index_bg_rho_ur])
      /(pvecback[pba->index_bg_H]
       *pvecback[pba->index_bg_H]);

    double one_plus_w_b_Omega_b =
      (1.+w_b)
      *(pvecback[pba->index_bg_rho_g]+pvecback[pba->index_bg_rho_ur])
      /(pvecback[pba->index_bg_H]
       *pvecback[pba->index_bg_H]);
    
    double one_plus_w_c_Omega_c =
      (1.+w_c)
      *(pvecback[pba->index_bg_rho_g]+pvecback[pba->index_bg_rho_ur])
      /(pvecback[pba->index_bg_H]
       *pvecback[pba->index_bg_H]);



    double w_m_prime =
      -3.*(
	   (1.+w_b)*w_b_Omega_b
	   +(1.+w_g)*w_g_Omega_g
	   +(1.+w_c)*w_c_Omega_c
	   )/Omega_m
      +3.*(
	   w_b_Omega_b
	   +w_g_Omega_g
	   +w_c_Omega_c
	   )
      *(
	   one_plus_w_b_Omega_b
	   +one_plus_w_g_Omega_g
	   +one_plus_w_c_Omega_c
	)
      /Omega_m
      /Omega_m;

    pvecback[pba->index_bg_ca2_m] =
      w_m - w_m_prime/(3.*(1.+w_m));

    pvecback[pba->index_bg_epsilon_H_bar] =  pvecback[pba->index_bg_epsilon_H]
      -(3./2)*((1.+w_b)*w_b_Omega_b
               +(1.+w_g)*w_g_Omega_g
               +(1.+w_c)*w_c_Omega_c
               +(1.+w_ds)*w_ds*pvecback[pba->index_bg_Omega_ds]);




   pvecback[pba->index_bg_epsilon_H_bar_prime] =
     pvecback[pba->index_bg_epsilon_H_bar]
     -4.*pvecback[pba->index_bg_epsilon_H]
     +2.*pvecback[pba->index_bg_epsilon_H]*pvecback[pba->index_bg_epsilon_H]
     -(3./2)*((1.+w_b)*(2.*pvecback[pba->index_bg_epsilon_H]
                        -3.*(1.+w_b))*w_b_Omega_b
               +(1.+w_g)*(2.*pvecback[pba->index_bg_epsilon_H]
                        -3.*(1.+w_g))*w_g_Omega_g
               +(1.+w_c)*(2.*pvecback[pba->index_bg_epsilon_H]
                        -3.*(1.+w_c))*w_c_Omega_c
               +(1.+w_ds)*(2.*pvecback[pba->index_bg_epsilon_H]
                        -3.*(1.+w_ds))*w_ds*pvecback[pba->index_bg_Omega_ds]);



     
     pvecback[pba->index_bg_Ricci_prime] =
       -6*(pvecback[pba->index_bg_H]
         *pvecback[pba->index_bg_H])
         *pvecback[pba->index_bg_epsilon_H_bar];

if (pba->mg_type == fR) {


 pvecback[pba->index_bg_fR_ds] =
   pvecback[pba->index_bg_f_prime_ds]
   /pvecback[pba->index_bg_Ricci_prime];


pvecback[pba->index_bg_fR_prime_ds] =
  (1.-pvecback[pba->index_bg_epsilon_H])
  *pvecback[pba->index_bg_fR_ds]
  -pvecback[pba->index_bg_f_ds]
  /(6.*(pvecback[pba->index_bg_H]
        *pvecback[pba->index_bg_H]))
  -pvecback[pba->index_bg_Omega_ds];
     
  pvecback[pba->index_bg_B_ds] =
    -pvecback[pba->index_bg_fR_prime_ds]
    /(pvecback[pba->index_bg_epsilon_H]
      *(1+pvecback[pba->index_bg_fR_ds]));


        pba->B0_ds = pvecback[pba->index_bg_B_ds];


	 }
  ////////////ALPHAS////////////
   
  
  pvecback[pba->index_bg_alpha_b_hi] =
     pba->alpha_b_hi
     *Omega_ds;
     

  
    pvecback[pba->index_bg_alpha_m_hi] = 
    pba->alpha_m_hi
      *Omega_ds; 

      
   pvecback[pba->index_bg_alpha_t_hi] =
     pba->alpha_t_hi
     *Omega_ds;
  

   pvecback[pba->index_bg_alpha_k_hi] = 
     pba->alpha_k_hi
     *Omega_ds;


   pvecback[pba->index_bg_alpha_b_hi_prime] =
     pba->alpha_b_hi
     *Omega_ds_prime;

   
   pvecback[pba->index_bg_alpha_k_hi_prime] =
     pba->alpha_k_hi
     *Omega_ds_prime;

   pvecback[pba->index_bg_m2_hi] = pvecback_B[pba->index_bi_m2_hi];
   pvecback[pba->index_bg_m2_minus_1_hi] =
     pvecback_B[pba->index_bi_m2_hi]-1./3.;

   //fR
      if (pba->mg_type == fR){

  pvecback[pba->index_bg_alpha_b_hi] =
   pvecback[pba->index_bg_fR_prime_ds] 
   /(2. 
    *(1+pvecback[pba->index_bg_fR_ds])); 
        
        pvecback[pba->index_bg_alpha_m_hi] = 
       2.*pvecback[pba->index_bg_alpha_b_hi]; 

   pvecback[pba->index_bg_alpha_t_hi] = 0.;
   pvecback[pba->index_bg_alpha_k_hi] = 0.;

        
  pvecback[pba->index_bg_m2_hi] =
  (
   1.
   +
   pvecback[pba->index_bg_fR_ds]
   )/3.;

   pvecback[pba->index_bg_m2_minus_1_hi] =
   pvecback[pba->index_bg_fR_ds]
   /3.;

   pvecback[pba->index_bg_alpha_b_hi_prime] = 0.;
   pvecback[pba->index_bg_alpha_k_hi_prime] = 0.;
   
    }
   
   pvecback[pba->index_bg_alpha_hi] =
     pvecback[pba->index_bg_alpha_k_hi]
     +6.*pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_b_hi];
   //Note: alpha_BBP = D_hiCLASS
   //alpha_b_BBP = alpha_b_hiCLASS/2


   double one_plus_w_m_Omega_m_over_M2 =
      (Omega_m + w_m_Omega_m)/3./pvecback[pba->index_bg_m2_hi];
 
   pvecback[pba->index_bg_gamma_0_hi] = 0.;

          if (pba->mg_type == hi){

   pvecback[pba->index_bg_cs2_hi] =//one_plus_w_m_Omega_m_over_M2 ;
     (1./pvecback[pba->index_bg_alpha_hi])
       *(
           2.*(1.
                +
              pvecback[pba->index_bg_alpha_b_hi]
             )
         *(pvecback[pba->index_bg_epsilon_H]
           +pvecback[pba->index_bg_alpha_m_hi]
           -pvecback[pba->index_bg_alpha_t_hi]
           -pvecback[pba->index_bg_alpha_b_hi]
           *(1.+pvecback[pba->index_bg_alpha_t_hi])
           )
           -2.*pvecback[pba->index_bg_alpha_b_hi_prime]
           -3.*one_plus_w_m_Omega_m_over_M2
           );

            
   pvecback[pba->index_bg_gamma_1_hi] =
     -(3./4.)
     *pvecback[pba->index_bg_alpha_k_hi]
     *one_plus_w_m_Omega_m_over_M2;
     +(1./2.)
     * pvecback[pba->index_bg_alpha_hi]
     *pvecback[pba->index_bg_epsilon_H];
     
     
     //TBC
   pvecback[pba->index_bg_Gamma_hi] =
      (
      pvecback[pba->index_bg_alpha_k_hi]
      *pvecback[pba->index_bg_alpha_b_hi]
      /pvecback[pba->index_bg_gamma_1_hi]
      *3/4.
       ) 
     *(
       (3.+pvecback[pba->index_bg_alpha_m_hi])
       *one_plus_w_m_Omega_m_over_M2
       +3.*w_m*one_plus_w_m_Omega_m_over_M2
       )
     +(
       pvecback[pba->index_bg_alpha_k_hi_prime]
      *pvecback[pba->index_bg_alpha_b_hi]
     -2.*pvecback[pba->index_bg_alpha_k_hi]
      *pvecback[pba->index_bg_alpha_b_hi_prime]
       )
     /pvecback[pba->index_bg_gamma_1_hi]
     *(
       3.*one_plus_w_m_Omega_m_over_M2
       -2.*pvecback[pba->index_bg_epsilon_H]
       )
     -pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_hi]
     /2.
     /pvecback[pba->index_bg_gamma_1_hi]
     *(pvecback[pba->index_bg_epsilon_H_bar]
       +4.*pvecback[pba->index_bg_epsilon_H]
       );

   
   
   pvecback[pba->index_bg_gamma_2_hi] =
      pvecback[pba->index_bg_cs2_hi]
     +pvecback[pba->index_bg_alpha_t_hi]/3.
     -(2./pvecback[pba->index_bg_alpha_hi])
     *(
       pvecback[pba->index_bg_alpha_b_hi]
       *(2.
         +pvecback[pba->index_bg_Gamma_hi])
       +(1.+pvecback[pba->index_bg_alpha_b_hi])
       *(pvecback[pba->index_bg_alpha_m_hi]
         -pvecback[pba->index_bg_alpha_t_hi])
       );

   pvecback[pba->index_bg_gamma_8_hi] =
     pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_t_hi]
     +pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_t_hi]
     -pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_m_hi];
     
   

   pvecback[pba->index_bg_gamma_3_hi] =
     pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_cs2_hi]
     +pvecback[pba->index_bg_gamma_8_hi]
     /3.;
   

   //TBC
   pvecback[pba->index_bg_gamma_4_hi] =
     -1./(
          3.*one_plus_w_m_Omega_m_over_M2
          -2.*pvecback[pba->index_bg_epsilon_H]
          )
     *(
       -9.*pvecback[pba->index_bg_alpha_k_hi]
       *w_m*one_plus_w_m_Omega_m_over_M2
       -pvecback[pba->index_bg_alpha_m_hi]
       *(
         3.*one_plus_w_m_Omega_m_over_M2
         +3.
         -2.*pvecback[pba->index_bg_epsilon_H]
         )
       +6.*pvecback[pba->index_bg_alpha_b_hi]
          *pvecback[pba->index_bg_alpha_b_hi]
          /pvecback[pba->index_bg_alpha_hi]
       *3.*one_plus_w_m_Omega_m_over_M2
       *(3.+pvecback[pba->index_bg_alpha_m_hi])
       +6.*pvecback[pba->index_bg_alpha_b_hi]
          /pvecback[pba->index_bg_alpha_hi]
       *pvecback[pba->index_bg_Gamma_hi]
       *3.*one_plus_w_m_Omega_m_over_M2
       +2.*pvecback[pba->index_bg_epsilon_H_bar]
       +3.*pvecback[pba->index_bg_alpha_m_hi]
       +2.*pvecback[pba->index_bg_epsilon_H]
       *(1.-pvecback[pba->index_bg_alpha_m_hi])
       );

   
   pvecback[pba->index_bg_gamma_5_hi] =
     -pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_b_hi]
     -(
       6.*pvecback[pba->index_bg_alpha_b_hi]
       -pvecback[pba->index_bg_alpha_k_hi]
       )
     *(
       pvecback[pba->index_bg_alpha_t_hi]
       -pvecback[pba->index_bg_alpha_m_hi]
       )
     /6.
     +pvecback[pba->index_bg_alpha_b_hi]
     /pvecback[pba->index_bg_alpha_hi]
     *(
       pvecback[pba->index_bg_alpha_k_hi_prime]
      *pvecback[pba->index_bg_alpha_b_hi]
     -2.*pvecback[pba->index_bg_alpha_k_hi]
      *pvecback[pba->index_bg_alpha_b_hi_prime]
       );


   pvecback[pba->index_bg_gamma_6_hi] =
     -6.
     *pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_b_hi]
     /pvecback[pba->index_bg_alpha_hi]
     *(2.+pvecback[pba->index_bg_Gamma_hi])
     +
     (
      pvecback[pba->index_bg_alpha_k_hi]
     *pvecback[pba->index_bg_alpha_m_hi]
     -6.*pvecback[pba->index_bg_alpha_b_hi]
        *pvecback[pba->index_bg_alpha_b_hi]
      )
     /pvecback[pba->index_bg_alpha_hi];
     

   
   pvecback[pba->index_bg_gamma_7_hi] =
     pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_b_hi]
     /3.
     /pvecback[pba->index_bg_alpha_hi]
     *(
       pvecback[pba->index_bg_alpha_k_hi]
       *pvecback[pba->index_bg_alpha_m_hi]
       -6.*pvecback[pba->index_bg_alpha_b_hi]
       *pvecback[pba->index_bg_alpha_b_hi]
       -(
         6.*pvecback[pba->index_bg_alpha_b_hi]
         -pvecback[pba->index_bg_alpha_k_hi]
         )
       *(
         pvecback[pba->index_bg_alpha_t_hi]
         -pvecback[pba->index_bg_alpha_m_hi]
         )
       );
   
   pvecback[pba->index_bg_gamma_7_not_tilde_hi] =
     1./3.
     /pvecback[pba->index_bg_alpha_hi]
     *(
       pvecback[pba->index_bg_alpha_k_hi]
       *pvecback[pba->index_bg_alpha_m_hi]
       -6.*pvecback[pba->index_bg_alpha_b_hi]
       *pvecback[pba->index_bg_alpha_b_hi]
       -(
         6.*pvecback[pba->index_bg_alpha_b_hi]
         -pvecback[pba->index_bg_alpha_k_hi]
         )
       *(
         pvecback[pba->index_bg_alpha_t_hi]
         -pvecback[pba->index_bg_alpha_m_hi]
         )
       );
   
   
   pvecback[pba->index_bg_gamma_9_hi] =
     pvecback[pba->index_bg_alpha_hi]
     /2.
     *(
       pvecback[pba->index_bg_alpha_t_hi]
       -pvecback[pba->index_bg_alpha_m_hi]
       );
   

      pvecback[pba->index_bg_gamma_10_hi] =
        3.
        *pvecback[pba->index_bg_alpha_b_hi]
        *pvecback[pba->index_bg_alpha_b_hi]
     *(
       pvecback[pba->index_bg_alpha_t_hi]
       -pvecback[pba->index_bg_alpha_m_hi]
       );

        }
    
   //fR
         if (pba->mg_type == fR){

   pvecback[pba->index_bg_cs2_hi] = 1.;


           
   pvecback[pba->index_bg_gamma_1_hi] =
     3.
     *pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_epsilon_H];

   
pvecback[pba->index_bg_Gamma_hi] =
  (
  pvecback[pba->index_bg_epsilon_H_bar]
  /pvecback[pba->index_bg_epsilon_H]
  -4.
    );
           
 pvecback[pba->index_bg_gamma_2_hi] =
   1./3.
     -(4.
       +pvecback[pba->index_bg_Gamma_hi]
       )
     /3./pvecback[pba->index_bg_alpha_b_hi];


   pvecback[pba->index_bg_gamma_8_hi] =
     -2.*pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_b_hi];

           
   pvecback[pba->index_bg_gamma_3_hi] =
     pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_b_hi]/3.;

   pvecback[pba->index_bg_gamma_4_hi] =
       -3.
       -pvecback[pba->index_bg_Gamma_hi];

           
  pvecback[pba->index_bg_gamma_5_hi] =
    pvecback[pba->index_bg_alpha_b_hi]
    *pvecback[pba->index_bg_alpha_b_hi];

           
   pvecback[pba->index_bg_gamma_6_hi] =
      pvecback[pba->index_bg_gamma_4_hi];

   pvecback[pba->index_bg_gamma_7_hi] =
     pvecback[pba->index_bg_gamma_3_hi];

  
   pvecback[pba->index_bg_gamma_7_not_tilde_hi] =
     1./3.;

  
   pvecback[pba->index_bg_gamma_9_hi] =
     -6.
     *pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_b_hi];
     


      pvecback[pba->index_bg_gamma_10_hi] =
     pvecback[pba->index_bg_gamma_9_hi];
   

     }


         if (pba->mg_type == qe){


  pvecback[pba->index_bg_alpha_b_hi] = 0.; 
  pvecback[pba->index_bg_alpha_m_hi] = 0.; 
   pvecback[pba->index_bg_alpha_t_hi] = 0.;
   
   pvecback[pba->index_bg_alpha_k_hi] = 3.*Omega_ds*(1.+w_ds);

        
   pvecback[pba->index_bg_m2_hi] =1./3.;
   pvecback[pba->index_bg_m2_minus_1_hi] =0.;

   pvecback[pba->index_bg_alpha_b_hi_prime] = 0.;
   pvecback[pba->index_bg_alpha_k_hi_prime] = 0.;
   pvecback[pba->index_bg_alpha_hi] =
     pvecback[pba->index_bg_alpha_k_hi]
     +6.*pvecback[pba->index_bg_alpha_b_hi]
     *pvecback[pba->index_bg_alpha_b_hi];
   //Note: alpha_BBP = D_hiCLASS
   //alpha_b_BBP = alpha_b_hiCLASS/2



   pvecback[pba->index_bg_cs2_hi] = pba->cs2_ds;
   pvecback[pba->index_bg_gamma_1_hi] =
     1./4.
     *pvecback[pba->index_bg_alpha_k_hi]
     *pvecback[pba->index_bg_alpha_k_hi];
   pvecback[pba->index_bg_Gamma_hi] = 0.;
   pvecback[pba->index_bg_gamma_2_hi] = pvecback[pba->index_bg_cs2_hi];
   pvecback[pba->index_bg_gamma_8_hi] = 0.;
   pvecback[pba->index_bg_gamma_3_hi] = 0.;
   pvecback[pba->index_bg_gamma_4_hi] = 3.*w_ds; //TBC
   pvecback[pba->index_bg_gamma_5_hi] = 0.;
   pvecback[pba->index_bg_gamma_6_hi] = 0.;
   pvecback[pba->index_bg_gamma_7_hi] = 0.;
   pvecback[pba->index_bg_gamma_7_not_tilde_hi] = 0.;
   pvecback[pba->index_bg_gamma_9_hi] = 0.;
   pvecback[pba->index_bg_gamma_10_hi] = 0.;
         }//qe
         

  
    
  }



  
  /** - compute other quantities in the exhaustive, redundant format */
  if (return_format == pba->long_info) {

    /** - compute critical density */
    pvecback[pba->index_bg_rho_crit] = rho_tot-pba->K/a/a;
    class_test(pvecback[pba->index_bg_rho_crit] <= 0.,
               pba->error_message,
               "rho_crit = %e instead of strictly positive",pvecback[pba->index_bg_rho_crit]);

    /** - compute Omega_m */
    pvecback[pba->index_bg_Omega_m] = rho_m / rho_tot;

    /* one can put other variables here */
    /*  */
    /*  */

  }

  return _SUCCESS_;

}

/**
 * Initialize the background structure, and in particular the
 * background interpolation table.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input/Output: pointer to initialized background structure
 * @return the error status
 */

int background_init(
                    struct precision * ppr,
                    struct background * pba
                    ) {

  /** Summary: */

  /** - define local variables */
  int n_ncdm;
  double rho_ncdm_rel,rho_nu_rel;
  double Neff;
  int filenum=0;
  
  /* ---------F.P. 24-10-2025 ---------*/


  double b;			//for fR models, deviation from lcdm


  /*-------------------------------------*/
  

  /** - in verbose mode, provide some information */
  if (pba->background_verbose > 0) {
    printf("Running CLASS version %s\n",_VERSION_);
    printf("Computing background\n");
    /*   printf("alpha_b_hi=%e\n",pba->alpha_b_hi);
    printf("alpha_k_hi=%e\n",pba->alpha_k_hi);
    printf("alpha_m_hi=%e\n",pba->alpha_m_hi);
    printf("alpha_t_hi=%e\n",pba->alpha_t_hi);
    printf("m2_hi_ini=%e\n",pba->m2_hi_ini);
    */
  printf("Omega0_lambda=%e\n",pba->Omega0_lambda);
  printf("Omega0_m=%e\n", pba->Omega0_b + pba->Omega0_cdm);
  printf("Omega0_r=%e\n", pba->Omega0_g + pba->Omega0_ur);
    
    /* below we want to inform the user about ncdm species*/
    if (pba->N_ncdm > 0) {

      Neff = pba->Omega0_ur/7.*8./pow(4./11.,4./3.)/pba->Omega0_g;

      /* loop over ncdm species */
      for (n_ncdm=0;n_ncdm<pba->N_ncdm; n_ncdm++) {

        /* inform if p-s-d read in files */
        if (pba->got_files[n_ncdm] == _TRUE_) {
          printf(" -> ncdm species i=%d read from file %s\n",n_ncdm+1,pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
          filenum++;
        }

        /* call this function to get rho_ncdm */
        background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                                pba->w_ncdm_bg[n_ncdm],
                                pba->q_size_ncdm_bg[n_ncdm],
                                0.,
                                pba->factor_ncdm[n_ncdm],
                                0.,
                                NULL,
                                &rho_ncdm_rel,
                                NULL,
                                NULL,
                                NULL);

        /* inform user of the contribution of each species to
           radiation density (in relativistic limit): should be
           between 1.01 and 1.02 for each active neutrino species;
           evaluated as rho_ncdm/rho_nu_rel where rho_nu_rel is the
           density of one neutrino in the instantaneous decoupling
           limit, i.e. assuming T_nu=(4/11)^1/3 T_gamma (this comes
           from the definition of N_eff) */
        rho_nu_rel = 56.0/45.0*pow(_PI_,6)*pow(4.0/11.0,4.0/3.0)*_G_/pow(_h_P_,3)/pow(_c_,7)*
          pow(_Mpc_over_m_,2)*pow(pba->T_cmb*_k_B_,4);

        printf(" -> ncdm species i=%d sampled with %d (resp. %d) points for purpose of background (resp. perturbation) integration. In the relativistic limit it gives Delta N_eff = %g\n",
               n_ncdm+1,
               pba->q_size_ncdm_bg[n_ncdm],
               pba->q_size_ncdm[n_ncdm],
               rho_ncdm_rel/rho_nu_rel);

        Neff += rho_ncdm_rel/rho_nu_rel;

      }

      printf(" -> total N_eff = %g (sumed over ultra-relativistic and ncdm species)\n",Neff);

    }
  }

  /** - if shooting failed during input, catch the error here */
  class_test(pba->shooting_failed == _TRUE_,
             pba->error_message,
             "Shooting failed, try optimising input_get_guess(). Error message:\n\n%s",
             pba->shooting_error);

  /** - assign values to all indices in vectors of background quantities with background_indices()*/
  class_call(background_indices(pba),
             pba->error_message,
             pba->error_message);

  /** - control that cosmological parameter values make sense */

  /* H0 in Mpc^{-1} */
  class_test((pba->H0 < _H0_SMALL_)||(pba->H0 > _H0_BIG_),
             pba->error_message,
             "H0=%g out of bounds (%g<H0<%g) \n",pba->H0,_H0_SMALL_,_H0_BIG_);

  class_test(fabs(pba->h * 1.e5 / _c_  / pba->H0 -1.)>ppr->smallest_allowed_variation,
             pba->error_message,
             "inconsistency between Hubble and reduced Hubble parameters: you have H0=%f/Mpc=%fkm/s/Mpc, but h=%f",pba->H0,pba->H0/1.e5* _c_,pba->h);

  /* T_cmb in K */
  class_test((pba->T_cmb < _TCMB_SMALL_)||(pba->T_cmb > _TCMB_BIG_),
             pba->error_message,
             "T_cmb=%g out of bounds (%g<T_cmb<%g)",pba->T_cmb,_TCMB_SMALL_,_TCMB_BIG_);

  /* H0 in Mpc^{-1} */
  class_test((pba->Omega0_k < _OMEGAK_SMALL_)||(pba->Omega0_k > _OMEGAK_BIG_),
             pba->error_message,
             "Omegak = %g out of bounds (%g<Omegak<%g) \n",pba->Omega0_k,_OMEGAK_SMALL_,_OMEGAK_BIG_);

  /* fluid equation of state */
  if (pba->has_fld == _TRUE_) {
    class_test(pba->w0_fld+pba->wa_fld>=1./3.,
               pba->error_message,
               "Your choice for w0_fld+wa_fld=%g is suspicious, there would not be radiation domination at early times\n",
               pba->w0_fld+pba->wa_fld);
  }

  /* in verbose mode, inform the user about the value of the ncdm
     masses in eV and about the ratio [m/omega_ncdm] in eV (the usual
     93 point something)*/
  if ((pba->background_verbose > 0) && (pba->has_ncdm == _TRUE_)) {
    for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
      printf(" -> non-cold dark matter species with i=%d has m_i = %e eV (so m_i / omega_i =%e eV)\n",
             n_ncdm+1,
             pba->m_ncdm_in_eV[n_ncdm],
             pba->m_ncdm_in_eV[n_ncdm]*pba->deg_ncdm[n_ncdm]/pba->Omega0_ncdm[n_ncdm]/pba->h/pba->h);
    }
  }

  /* check other quantities which would lead to segmentation fault if zero */
  class_test(pba->a_today <= 0,
             pba->error_message,
             "input a_today = %e instead of strictly positive",pba->a_today);

  class_test(_Gyr_over_Mpc_ <= 0,
             pba->error_message,
             "_Gyr_over_Mpc = %e instead of strictly positive",_Gyr_over_Mpc_);

  // printf("Omega0_ds=%e\n",pba->Omega0_ds);


if (pba->has_loop_over_b_ds == 1 && pba->has_ds == _TRUE_
    && pba->mg_type == fR) {
  
   class_call(background_loop_over_b_ds(ppr,pba),
             pba->error_message,
             pba->error_message);
   //return _SUCCESS_;
   if (pba->background_verbose>0) printf("B0_ds = %e w_ds =  %e, Omega_m = %e\n",pba->B0_ds,pba->w_ds,1.-pba->Omega0_ds);

 }

// pba->b_ds =  -3.17035e-01;

 
  /** - this function integrates the background over time, allocates
      and fills the background table */
  class_call(background_solve(ppr,pba),
             pba->error_message,
             pba->error_message);

  

  
  return _SUCCESS_;

}

/**
 * Free all memory space allocated by background_init().
 *
 *
 * @param pba Input: pointer to background structure (to be freed)
 * @return the error status
 */

int background_free(
                    struct background *pba
                    ) {
  int err;

  free(pba->tau_table);
  free(pba->z_table);
  free(pba->d2tau_dz2_table);
  free(pba->background_table);
  free(pba->d2background_dtau2_table);

  err = background_free_input(pba);

  return err;
}

/**
 * Free pointers inside background structure which were
 * allocated in input_read_parameters()
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */

int background_free_input(
                          struct background *pba
                          ) {

  int k;
  if (pba->Omega0_ncdm_tot != 0.){
    for(k=0; k<pba->N_ncdm; k++){
      free(pba->q_ncdm[k]);
      free(pba->w_ncdm[k]);
      free(pba->q_ncdm_bg[k]);
      free(pba->w_ncdm_bg[k]);
      free(pba->dlnf0_dlnq_ncdm[k]);
    }

    free(pba->q_ncdm);
    free(pba->w_ncdm);
    free(pba->q_ncdm_bg);
    free(pba->w_ncdm_bg);
    free(pba->dlnf0_dlnq_ncdm);
    free(pba->q_size_ncdm);
    free(pba->q_size_ncdm_bg);
    free(pba->M_ncdm);
    free(pba->T_ncdm);
    free(pba->ksi_ncdm);
    free(pba->deg_ncdm);
    free(pba->Omega0_ncdm);
    free(pba->m_ncdm_in_eV);
    free(pba->factor_ncdm);
    if(pba->got_files!=NULL)
      free(pba->got_files);
    if(pba->ncdm_psd_files!=NULL)
      free(pba->ncdm_psd_files);
    if(pba->ncdm_psd_parameters!=NULL)
      free(pba->ncdm_psd_parameters);
  }

  if (pba->Omega0_scf != 0.){
    if (pba->scf_parameters != NULL)
      free(pba->scf_parameters);
  }
  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of background quantities.
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */

int background_indices(
                       struct background *pba
                       ) {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of background quantities */
  int index_bg;
  /* a running index for the vector of background quantities to be integrated */
  int index_bi;

  /** - initialize all flags: which species are present? */

  pba->has_cdm = _FALSE_;
  pba->has_ncdm = _FALSE_;
  pba->has_dcdm = _FALSE_;
  pba->has_dr = _FALSE_;
  pba->has_scf = _FALSE_;
  pba->has_lambda = _FALSE_;
  pba->has_fld = _FALSE_;
  pba->has_ur = _FALSE_;
  pba->has_curvature = _FALSE_;

  if (pba->Omega0_cdm != 0.)
    pba->has_cdm = _TRUE_;


 if (pba->Omega0_ds != 0.){
   pba->has_ds = _TRUE_;
  }

  
  if (pba->Omega0_ncdm_tot != 0.)
    pba->has_ncdm = _TRUE_;

  if (pba->Omega0_dcdmdr != 0.){
    pba->has_dcdm = _TRUE_;
    if (pba->Gamma_dcdm != 0.)
      pba->has_dr = _TRUE_;
  }

  if (pba->Omega0_scf != 0.)
    pba->has_scf = _TRUE_;

  if (pba->Omega0_lambda != 0.)
    pba->has_lambda = _TRUE_;

  if (pba->Omega0_fld != 0.)
    pba->has_fld = _TRUE_;

  if (pba->Omega0_ur != 0.)
    pba->has_ur = _TRUE_;

  if (pba->sgnK != 0)
    pba->has_curvature = _TRUE_;

  /** - initialize all indices */

  index_bg=0;

  /* index for scale factor */
  class_define_index(pba->index_bg_a,_TRUE_,index_bg,1);

  /* - indices for H and its conformal-time-derivative */
  class_define_index(pba->index_bg_H,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_H_prime,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_epsilon_H,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_epsilon_H_bar,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_epsilon_H_bar_prime,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_w_m,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_ca2_m,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_Ricci,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_Ricci_prime,_TRUE_,index_bg,1);



  
  /* - end of indices in the short vector of background values */
  pba->bg_size_short = index_bg;

  /* - index for rho_g (photon density) */
  class_define_index(pba->index_bg_rho_g,_TRUE_,index_bg,1);

  /* - index for rho_b (baryon density) */
  class_define_index(pba->index_bg_rho_b,_TRUE_,index_bg,1);

  /* - index for rho_cdm */
  class_define_index(pba->index_bg_rho_cdm,pba->has_cdm,index_bg,1);



 /* - index for ds */
  class_define_index(pba->index_bg_Ricci,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_p_ds,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_Omega_ds,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_rho_ds,pba->has_ds,index_bg,1);

  if (pba->mg_type == fR) {
  
  class_define_index(pba->index_bg_f_ds,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_f_prime_ds,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_fR_ds,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_fR_prime_ds,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_B_ds,pba->has_ds,index_bg,1);

  }

  
  /*Functions for Horndeski*/

  class_define_index(pba->index_bg_alpha_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_alpha_m_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_alpha_t_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_alpha_b_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_alpha_b_hi_prime,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_alpha_k_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_alpha_k_hi_prime,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_m2_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_m2_minus_1_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_Gamma_hi,pba->has_ds,index_bg,1);


  class_define_index(pba->index_bg_gamma_0_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_gamma_1_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_gamma_2_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_gamma_3_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_gamma_4_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_gamma_5_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_gamma_6_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_gamma_7_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_gamma_7_not_tilde_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_gamma_8_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_gamma_9_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_gamma_10_hi,pba->has_ds,index_bg,1);
  class_define_index(pba->index_bg_cs2_hi,pba->has_ds,index_bg,1);






  
  /* - indices for ncdm. We only define the indices for ncdm1
     (density, pressure, pseudo-pressure), the other ncdm indices
     are contiguous */
  class_define_index(pba->index_bg_rho_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_p_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_pseudo_p_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);

  /* - index for dcdm */
  class_define_index(pba->index_bg_rho_dcdm,pba->has_dcdm,index_bg,1);

  /* - index for dr */
  class_define_index(pba->index_bg_rho_dr,pba->has_dr,index_bg,1);

  /* - indices for scalar field */
  class_define_index(pba->index_bg_phi_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_phi_prime_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_V_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_dV_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_ddV_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_rho_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_p_scf,pba->has_scf,index_bg,1);

  /* - index for Lambda */
  class_define_index(pba->index_bg_rho_lambda,pba->has_lambda,index_bg,1);

  /* - index for fluid */
  class_define_index(pba->index_bg_rho_fld,pba->has_fld,index_bg,1);

  /* - index for ultra-relativistic neutrinos/species */
  class_define_index(pba->index_bg_rho_ur,pba->has_ur,index_bg,1);

  /* - index for Omega_r (relativistic density fraction) */
  class_define_index(pba->index_bg_Omega_r,_TRUE_,index_bg,1);	

  /* - put here additional ingredients that you want to appear in the
     normal vector */
  /*    */
  /*    */
   /* -> density growth factor in wCDM universe */
  class_define_index(pba->index_bg_D_wCDM,_TRUE_,index_bg,1);

  

  /* - end of indices in the normal vector of background values */
  pba->bg_size_normal = index_bg;

  /* - indices in the long version : */

  /* -> critical density */
  class_define_index(pba->index_bg_rho_crit,_TRUE_,index_bg,1);

  /* - index for Omega_m (non-relativistic density fraction) */
  class_define_index(pba->index_bg_Omega_m,_TRUE_,index_bg,1);

  /* -> conformal distance */
  class_define_index(pba->index_bg_conf_distance,_TRUE_,index_bg,1);

  /* -> angular diameter distance */
  class_define_index(pba->index_bg_ang_distance,_TRUE_,index_bg,1);

  /* -> luminosity distance */
  class_define_index(pba->index_bg_lum_distance,_TRUE_,index_bg,1);

  /* -> proper time (for age of the Universe) */
  class_define_index(pba->index_bg_time,_TRUE_,index_bg,1);

  /* -> conformal sound horizon */
  class_define_index(pba->index_bg_rs,_TRUE_,index_bg,1);

  /* -> density growth factor in dust universe */
  class_define_index(pba->index_bg_D,_TRUE_,index_bg,1);

  /* -> velocity growth factor in dust universe */
  class_define_index(pba->index_bg_f,_TRUE_,index_bg,1);

  /* -> put here additional quantities describing background */
  /*    */
  /*    */

  /* -> end of indices in the long vector of background values */
  pba->bg_size = index_bg;

  /* - now, indices in vector of variables to integrate.
     First {B} variables, then {C} variables. */

  index_bi=0;

  /* -> scale factor */
  class_define_index(pba->index_bi_a,_TRUE_,index_bi,1);

  /* -> energy density in DCDM */
  class_define_index(pba->index_bi_rho_dcdm,pba->has_dcdm,index_bi,1);

  if(pba->mg_type == fR){
  class_define_index(pba->index_bi_f_ds,pba->has_ds,index_bi,1);
  class_define_index(pba->index_bi_f_prime_ds,pba->has_ds,index_bi,1);
 }
  
  /* -> energy density in DR */
  class_define_index(pba->index_bi_rho_dr,pba->has_dr,index_bi,1);

  /* -> scalar field and its derivative wrt conformal time (Zuma) */
  class_define_index(pba->index_bi_phi_scf,pba->has_scf,index_bi,1);
  class_define_index(pba->index_bi_phi_prime_scf,pba->has_scf,index_bi,1);

  /* End of {B} variables, now continue with {C} variables */
  pba->bi_B_size = index_bi;

  /* -> proper time (for age of the Universe) */
  class_define_index(pba->index_bi_time,_TRUE_,index_bi,1);

  /* -> sound horizon */
  class_define_index(pba->index_bi_rs,_TRUE_,index_bi,1);

  /* -> integral for growth factor */
  class_define_index(pba->index_bi_growth,_TRUE_,index_bi,1);
  class_define_index(pba->index_bi_m2_hi,pba->has_ds,index_bi,1);

 /* -> integral for growth factor in wCDM*/
  class_define_index(pba->index_bi_growth_wCDM,_TRUE_,index_bi,1);
  class_define_index(pba->index_bi_growth_prime_wCDM,_TRUE_,index_bi,1);


  
  /* -> index for conformal time in vector of variables to integrate */
  class_define_index(pba->index_bi_tau,_TRUE_,index_bi,1);

  /* -> end of indices in the vector of variables to integrate */
  pba->bi_size = index_bi;

  /* index_bi_tau must be the last index, because tau is part of this vector for the purpose of being stored, */
  /* but it is not a quantity to be integrated (since integration is over tau itself) */
  class_test(pba->index_bi_tau != index_bi-1,
             pba->error_message,
             "background integration requires index_bi_tau to be the last of all index_bi's");

  /* flags for calling the interpolation routine */

  pba->short_info=0;
  pba->normal_info=1;
  pba->long_info=2;

  pba->inter_normal=0;
  pba->inter_closeby=1;

  return _SUCCESS_;

}

/**
 * This is the routine where the distribution function f0(q) of each
 * ncdm species is specified (it is the only place to modify if you
 * need a partlar f0(q))
 *
 * @param pbadist Input:  structure containing all parameters defining f0(q)
 * @param q       Input:  momentum
 * @param f0      Output: phase-space distribution
 */

int background_ncdm_distribution(
                                 void * pbadist,
                                 double q,
                                 double * f0
                                 ) {
  struct background * pba;
  struct background_parameters_for_distributions * pbadist_local;
  int n_ncdm,lastidx;
  double ksi;
  double qlast,dqlast,f0last,df0last;
  double *param;
  /* Variables corresponding to entries in param: */
  //double square_s12,square_s23,square_s13;
  //double mixing_matrix[3][3];
  //int i;

  /** - extract from the input structure pbadist all the relevant information */
  pbadist_local = pbadist;          /* restore actual format of pbadist */
  pba = pbadist_local->pba;         /* extract the background structure from it */
  param = pba->ncdm_psd_parameters; /* extract the optional parameter list from it */
  n_ncdm = pbadist_local->n_ncdm;   /* extract index of ncdm species under consideration */
  ksi = pba->ksi_ncdm[n_ncdm];      /* extract chemical potential */

  /** - shall we interpolate in file, or shall we use analytical formula below? */

  /** - a) deal first with the case of interpolating in files */
  if (pba->got_files[n_ncdm]==_TRUE_) {

    lastidx = pbadist_local->tablesize-1;
    if(q<pbadist_local->q[0]){
      //Handle q->0 case:
      *f0 = pbadist_local->f0[0];
    }
    else if(q>pbadist_local->q[lastidx]){
      //Handle q>qmax case (ensure continuous and derivable function with Boltzmann tail):
      qlast=pbadist_local->q[lastidx];
      f0last=pbadist_local->f0[lastidx];
      dqlast=qlast - pbadist_local->q[lastidx-1];
      df0last=f0last - pbadist_local->f0[lastidx-1];

      *f0 = f0last*exp(-(qlast-q)*df0last/f0last/dqlast);
    }
    else{
      //Do interpolation:
      class_call(array_interpolate_spline(
                                          pbadist_local->q,
                                          pbadist_local->tablesize,
                                          pbadist_local->f0,
                                          pbadist_local->d2f0,
                                          1,
                                          q,
                                          &pbadist_local->last_index,
                                          f0,
                                          1,
                                          pba->error_message),
                 pba->error_message,     pba->error_message);
    }
  }

  /** - b) deal now with case of reading analytical function */
  else{
    /**
       Next enter your analytic expression(s) for the p.s.d.'s. If
       you need different p.s.d.'s for different species, put each
       p.s.d inside a condition, like for instance: if (n_ncdm==2) 
       {*f0=...}.  Remember that n_ncdm = 0 refers to the first
       species.
    */

    /**************************************************/
    /*    FERMI-DIRAC INCLUDING CHEMICAL POTENTIALS   */
    /**************************************************/

    *f0 = 1.0/pow(2*_PI_,3)*(1./(exp(q-ksi)+1.) +1./(exp(q+ksi)+1.));

    /**************************************************/

    /** This form is only appropriate for approximate studies, since in
        reality the chemical potentials are associated with flavor
        eigenstates, not mass eigenstates. It is easy to take this into
        account by introducing the mixing angles. In the later part
        (not read by the code) we illustrate how to do this. */

    if (_FALSE_) {

      /* We must use the list of extra parameters read in input, stored in the
         ncdm_psd_parameter list, extracted above from the structure
         and now called param[..] */

      /* check that this list has been read */
      class_test(param == NULL,
                 pba->error_message,
                 "Analytic expression wants to use 'ncdm_psd_parameters', but they have not been entered!");

      /* extract values from the list (in this example, mixing angles) */
      double square_s12=param[0];
      double square_s23=param[1];
      double square_s13=param[2];

      /* infer mixing matrix */
      double mixing_matrix[3][3];
      int i;

      mixing_matrix[0][0]=pow(fabs(sqrt((1-square_s12)*(1-square_s13))),2);
      mixing_matrix[0][1]=pow(fabs(sqrt(square_s12*(1-square_s13))),2);
      mixing_matrix[0][2]=fabs(square_s13);
      mixing_matrix[1][0]=pow(fabs(sqrt((1-square_s12)*square_s13*square_s23)+sqrt(square_s12*(1-square_s23))),2);
      mixing_matrix[1][1]=pow(fabs(sqrt(square_s12*square_s23*square_s13)-sqrt((1-square_s12)*(1-square_s23))),2);
      mixing_matrix[1][2]=pow(fabs(sqrt(square_s23*(1-square_s13))),2);
      mixing_matrix[2][0]=pow(fabs(sqrt(square_s12*square_s23)-sqrt((1-square_s12)*square_s13*(1-square_s23))),2);
      mixing_matrix[2][1]=pow(sqrt((1-square_s12)*square_s23)+sqrt(square_s12*square_s13*(1-square_s23)),2);
      mixing_matrix[2][2]=pow(fabs(sqrt((1-square_s13)*(1-square_s23))),2);

      /* loop over flavor eigenstates and compute psd of mass eigenstates */
      *f0=0.0;
      for(i=0;i<3;i++){

    	*f0 += mixing_matrix[i][n_ncdm]*1.0/pow(2*_PI_,3)*(1./(exp(q-pba->ksi_ncdm[i])+1.) +1./(exp(q+pba->ksi_ncdm[i])+1.));

      }
    } /* end of region not used, but shown as an example */
  }

  return _SUCCESS_;
}

/**
 * This function is only used for the purpose of finding optimal
 * quadrature weights. The logic is: if we can accurately convolve
 * f0(q) with this function, then we can convolve it accurately with
 * any other relevant function.
 *
 * @param pbadist Input:  structure containing all background parameters
 * @param q       Input:  momentum
 * @param test    Output: value of the test function test(q)
 */

int background_ncdm_test_function(
                                  void * pbadist,
                                  double q,
                                  double * test
                                  ) {

  double c = 2.0/(3.0*_zeta3_);
  double d = 120.0/(7.0*pow(_PI_,4));
  double e = 2.0/(45.0*_zeta5_);

  /** Using a + bq creates problems for otherwise acceptable distributions
      which diverges as \f$ 1/r \f$ or \f$ 1/r^2 \f$ for \f$ r\to 0 \f$*/
  *test = pow(2.0*_PI_,3)/6.0*(c*q*q-d*q*q*q-e*q*q*q*q);

  return _SUCCESS_;
}

/**
 * This function finds optimal quadrature weights for each ncdm
 * species
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int background_ncdm_init(
                         struct precision *ppr,
                         struct background *pba
                         ) {

  int index_q, k,tolexp,row,status,filenum;
  double f0m2,f0m1,f0,f0p1,f0p2,dq,q,df0dq,tmp1,tmp2;
  struct background_parameters_for_distributions pbadist;
  FILE *psdfile;

  pbadist.pba = pba;

  /* Allocate pointer arrays: */
  class_alloc(pba->q_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->dlnf0_dlnq_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);

  /* Allocate pointers: */
  class_alloc(pba->q_size_ncdm,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_size_ncdm_bg,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->factor_ncdm,sizeof(double)*pba->N_ncdm,pba->error_message);

  for(k=0, filenum=0; k<pba->N_ncdm; k++){
    pbadist.n_ncdm = k;
    pbadist.q = NULL;
    pbadist.tablesize = 0;
    /*Do we need to read in a file to interpolate the distribution function? */
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)){
      psdfile = fopen(pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_,"r");
      class_test(psdfile == NULL,pba->error_message,
                 "Could not open file %s!",pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
      // Find size of table:
      for (row=0,status=2; status==2; row++){
        status = fscanf(psdfile,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(psdfile);
      pbadist.tablesize = row-1;

      /*Allocate room for interpolation table: */
      class_alloc(pbadist.q,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.d2f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      for (row=0; row<pbadist.tablesize; row++){
        status = fscanf(psdfile,"%lf %lf",
                        &pbadist.q[row],&pbadist.f0[row]);
        //		printf("(q,f0) = (%g,%g)\n",pbadist.q[row],pbadist.f0[row]);
      }
      fclose(psdfile);
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(pbadist.q,
                                          pbadist.tablesize,
                                          pbadist.f0,
                                          1,
                                          pbadist.d2f0,
                                          _SPLINE_EST_DERIV_,
                                          pba->error_message),
                 pba->error_message,
                 pba->error_message);
      filenum++;
    }

    /* Handle perturbation qsampling: */
    class_alloc(pba->q_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);
    class_alloc(pba->w_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);

    class_call(get_qsampling(pba->q_ncdm[k],
                             pba->w_ncdm[k],
                             &(pba->q_size_ncdm[k]),
                             _QUADRATURE_MAX_,
                             ppr->tol_ncdm,
                             pbadist.q,
                             pbadist.tablesize,
                             background_ncdm_test_function,
                             background_ncdm_distribution,
                             &pbadist,
                             pba->error_message),
               pba->error_message,
               pba->error_message);
    pba->q_ncdm[k]=realloc(pba->q_ncdm[k],pba->q_size_ncdm[k]*sizeof(double));
    pba->w_ncdm[k]=realloc(pba->w_ncdm[k],pba->q_size_ncdm[k]*sizeof(double));


    if (pba->background_verbose > 0)
      printf("ncdm species i=%d sampled with %d points for purpose of perturbation integration\n",
             k+1,
             pba->q_size_ncdm[k]);

    /* Handle background q_sampling: */
    class_alloc(pba->q_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);
    class_alloc(pba->w_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);

    class_call(get_qsampling(pba->q_ncdm_bg[k],
                             pba->w_ncdm_bg[k],
                             &(pba->q_size_ncdm_bg[k]),
                             _QUADRATURE_MAX_BG_,
                             ppr->tol_ncdm_bg,
                             pbadist.q,
                             pbadist.tablesize,
                             background_ncdm_test_function,
                             background_ncdm_distribution,
                             &pbadist,
                             pba->error_message),
               pba->error_message,
               pba->error_message);


    pba->q_ncdm_bg[k]=realloc(pba->q_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double));
    pba->w_ncdm_bg[k]=realloc(pba->w_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double));

    /** - in verbose mode, inform user of number of sampled momenta
        for background quantities */
    if (pba->background_verbose > 0)
      printf("ncdm species i=%d sampled with %d points for purpose of background integration\n",
             k+1,
             pba->q_size_ncdm_bg[k]);

    class_alloc(pba->dlnf0_dlnq_ncdm[k],
                pba->q_size_ncdm[k]*sizeof(double),
                pba->error_message);


    for (index_q=0; index_q<pba->q_size_ncdm[k]; index_q++) {
      q = pba->q_ncdm[k][index_q];
      class_call(background_ncdm_distribution(&pbadist,q,&f0),
                 pba->error_message,pba->error_message);

      //Loop to find appropriate dq:
      for(tolexp=_PSD_DERIVATIVE_EXP_MIN_; tolexp<_PSD_DERIVATIVE_EXP_MAX_; tolexp++){

        if (index_q == 0){
          dq = MIN((0.5-ppr->smallest_allowed_variation)*q,2*exp(tolexp)*(pba->q_ncdm[k][index_q+1]-q));
        }
        else if (index_q == pba->q_size_ncdm[k]-1){
          dq = exp(tolexp)*2.0*(pba->q_ncdm[k][index_q]-pba->q_ncdm[k][index_q-1]);
        }
        else{
          dq = exp(tolexp)*(pba->q_ncdm[k][index_q+1]-pba->q_ncdm[k][index_q-1]);
        }

        class_call(background_ncdm_distribution(&pbadist,q-2*dq,&f0m2),
                   pba->error_message,pba->error_message);
        class_call(background_ncdm_distribution(&pbadist,q+2*dq,&f0p2),
                   pba->error_message,pba->error_message);

        if (fabs((f0p2-f0m2)/f0)>sqrt(ppr->smallest_allowed_variation)) break;
      }

      class_call(background_ncdm_distribution(&pbadist,q-dq,&f0m1),
                 pba->error_message,pba->error_message);
      class_call(background_ncdm_distribution(&pbadist,q+dq,&f0p1),
                 pba->error_message,pba->error_message);
      //5 point estimate of the derivative:
      df0dq = (+f0m2-8*f0m1+8*f0p1-f0p2)/12.0/dq;
      //printf("df0dq[%g] = %g. dlf=%g ?= %g. f0 =%g.\n",q,df0dq,q/f0*df0dq,
      //Avoid underflow in extreme tail:
      if (fabs(f0)==0.)
        pba->dlnf0_dlnq_ncdm[k][index_q] = -q; /* valid for whatever f0 with exponential tail in exp(-q) */
      else
        pba->dlnf0_dlnq_ncdm[k][index_q] = q/f0*df0dq;
    }

    pba->factor_ncdm[k]=pba->deg_ncdm[k]*4*_PI_*pow(pba->T_cmb*pba->T_ncdm[k]*_k_B_,4)*8*_PI_*_G_
      /3./pow(_h_P_/2./_PI_,3)/pow(_c_,7)*_Mpc_over_m_*_Mpc_over_m_;

    /* If allocated, deallocate interpolation table:  */
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)){
      free(pbadist.q);
      free(pbadist.f0);
      free(pbadist.d2f0);
    }
  }


  return _SUCCESS_;
}

/**
 * For a given ncdm species: given the quadrature weights, the mass
 * and the redshift, find background quantities by a quick weighted
 * sum over.  Input parameters passed as NULL pointers are not
 * evaluated for speed-up
 *
 * @param qvec     Input: sampled momenta
 * @param wvec     Input: quadrature weights
 * @param qsize    Input: number of momenta/weights
 * @param M        Input: mass
 * @param factor   Input: normalization factor for the p.s.d.
 * @param z        Input: redshift
 * @param n        Output: number density
 * @param rho      Output: energy density
 * @param p        Output: pressure
 * @param drho_dM  Output: derivative used in next function
 * @param pseudo_p Output: pseudo-pressure used in perturbation module for fluid approx
 *
 */

int background_ncdm_momenta(
                            /* Only calculate for non-NULL pointers: */
                            double * qvec,
                            double * wvec,
                            int qsize,
                            double M,
                            double factor,
                            double z,
                            double * n,
                            double * rho, // density
                            double * p,   // pressure
                            double * drho_dM,  // d rho / d M used in next function
                            double * pseudo_p  // pseudo-p used in ncdm fluid approx
                            ) {

  int index_q;
  double epsilon;
  double q2;
  double factor2;
  /** Summary: */
  
  /** - rescale normalization at given redshift */
  factor2 = factor*pow(1+z,4);

  /** - initialize quantities */
  if (n!=NULL) *n = 0.;
  if (rho!=NULL) *rho = 0.;
  if (p!=NULL) *p = 0.;
  if (drho_dM!=NULL) *drho_dM = 0.;
  if (pseudo_p!=NULL) *pseudo_p = 0.;

  /** - loop over momenta */
  for (index_q=0; index_q<qsize; index_q++) {

    /* squared momentum */
    q2 = qvec[index_q]*qvec[index_q];

    /* energy */
    epsilon = sqrt(q2+M*M/(1.+z)/(1.+z));

    /* integrand of the various quantities */
    if (n!=NULL) *n += q2*wvec[index_q];
    if (rho!=NULL) *rho += q2*epsilon*wvec[index_q];
    if (p!=NULL) *p += q2*q2/3./epsilon*wvec[index_q];
    if (drho_dM!=NULL) *drho_dM += q2*M/(1.+z)/(1.+z)/epsilon*wvec[index_q];
    if (pseudo_p!=NULL) *pseudo_p += pow(q2/epsilon,3)/3.0*wvec[index_q];
  }

  /** - adjust normalization */
  if (n!=NULL) *n *= factor2*(1.+z);
  if (rho!=NULL) *rho *= factor2;
  if (p!=NULL) *p *= factor2;
  if (drho_dM!=NULL) *drho_dM *= factor2;
  if (pseudo_p!=NULL) *pseudo_p *=factor2;

  return _SUCCESS_;
}

/**
 * When the user passed the density fraction Omega_ncdm or
 * omega_ncdm in input but not the mass, infer the mass with Newton iteration method.
 *
 * @param ppr    Input: precision structure
 * @param pba    Input/Output: background structure
 * @param n_ncdm Input: index of ncdm species
 */

int background_ncdm_M_from_Omega(
                                 struct precision *ppr,
                                 struct background *pba,
                                 int n_ncdm
                                 ) {
  double rho0,rho,n,M,deltaM,drhodM;
  int iter,maxiter=50;

  rho0 = pba->H0*pba->H0*pba->Omega0_ncdm[n_ncdm]; /*Remember that rho is defined such that H^2=sum(rho_i) */
  M = 0.0;

  background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                          pba->w_ncdm_bg[n_ncdm],
                          pba->q_size_ncdm_bg[n_ncdm],
                          M,
                          pba->factor_ncdm[n_ncdm],
                          0.,
                          &n,
                          &rho,
                          NULL,
                          NULL,
                          NULL);

  /* Is the value of Omega less than a massless species?*/
  class_test(rho0<rho,pba->error_message,
             "The value of Omega for the %dth species, %g, is less than for a massless species! It should be atleast %g. Check your input.",
             n_ncdm,pba->Omega0_ncdm[n_ncdm],pba->Omega0_ncdm[n_ncdm]*rho/rho0);

  /* In the strict NR limit we have rho = n*(M) today, giving a zeroth order guess: */
  M = rho0/n; /* This is our guess for M. */
  for (iter=1; iter<=maxiter; iter++){

    /* Newton iteration. First get relevant quantities at M: */
    background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                            pba->w_ncdm_bg[n_ncdm],
                            pba->q_size_ncdm_bg[n_ncdm],
                            M,
                            pba->factor_ncdm[n_ncdm],
                            0.,
                            NULL,
                            &rho,
                            NULL,
                            &drhodM,
                            NULL);

    deltaM = (rho0-rho)/drhodM; /* By definition of the derivative */
    if ((M+deltaM)<0.0) deltaM = -M/2.0; /* Avoid overshooting to negative M value. */
    M += deltaM; /* Update value of M.. */
    if (fabs(deltaM/M)<ppr->tol_M_ncdm){
      /* Accuracy reached.. */
      pba->M_ncdm[n_ncdm] = M;
      break;
    }
  }
  class_test(iter>=maxiter,pba->error_message,
             "Newton iteration could not converge on a mass for some reason.");
  return _SUCCESS_;
}

/**
 *  This function integrates the background over time, allocates and
 *  fills the background table
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int background_solve(
                     struct precision *ppr,
                     struct background *pba
                     ) {

  /** Summary: */

  /** - define local variables */

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;
  /* parameters and workspace for the background_derivs function */
  struct background_parameters_and_workspace bpaw;
  /* a growing table (since the number of time steps is not known a priori) */
  growTable gTable;
  /* needed for growing table */
  double * pData;
  /* needed for growing table */
  void * memcopy_result;
  /* initial conformal time */
  double tau_start;
  /* final conformal time */
  double tau_end;
  /* an index running over bi indices */
  int i;
  /* vector of quantities to be integrated */
  double * pvecback_integration;
  /* vector of all background quantities */
  double * pvecback;
  /* necessary for calling array_interpolate(), but never used */
  int last_index=0;
  /* comoving radius coordinate in Mpc (equal to conformal distance in flat case) */
  double comoving_radius=0.;

  bpaw.pba = pba;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  bpaw.pvecback = pvecback;

  /** - allocate vector of quantities to be integrated */
  class_alloc(pvecback_integration,pba->bi_size*sizeof(double),pba->error_message);

  /** - initialize generic integrator with initialize_generic_integrator() */

  /* Size of vector to integrate is (pba->bi_size-1) rather than
   * (pba->bi_size), since tau is not integrated.
   */
  class_call(initialize_generic_integrator((pba->bi_size-1),&gi),
             gi.error_message,
             pba->error_message);

  /** - impose initial conditions with background_initial_conditions() */
  class_call(background_initial_conditions(ppr,pba,pvecback,pvecback_integration),
             pba->error_message,
             pba->error_message);

  /* here tau_end is in fact the initial time (in the next loop
     tau_start = tau_end) */
  tau_end=pvecback_integration[pba->index_bi_tau];

  /** - create a growTable with gt_init() */
  class_call(gt_init(&gTable),
             gTable.error_message,
             pba->error_message);

  /* initialize the counter for the number of steps */
  pba->bt_size=0;

  /** - loop over integration steps: call background_functions(), find step size, save data in growTable with gt_add(), perform one step with generic_integrator(), store new value of tau */

  while (pvecback_integration[pba->index_bi_a] < pba->a_today) {

    tau_start = tau_end;

    /* -> find step size (trying to adjust the last step as close as possible to the one needed to reach a=a_today; need not be exact, difference corrected later) */
    class_call(background_functions(pba,pvecback_integration, pba->short_info, pvecback),
               pba->error_message,
               pba->error_message);

    if ((pvecback_integration[pba->index_bi_a]*(1.+ppr->back_integration_stepsize)) < pba->a_today) {
      tau_end = tau_start + ppr->back_integration_stepsize / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }
    else {
      tau_end = tau_start + (pba->a_today/pvecback_integration[pba->index_bi_a]-1.) / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }

    class_test((tau_end-tau_start)/tau_start < ppr->smallest_allowed_variation,
               pba->error_message,
               "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",(tau_end-tau_start)/tau_start);

    /* -> save data in growTable */
    class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
               gTable.error_message,
               pba->error_message);
    pba->bt_size++;

    /* -> perform one step */
    class_call(generic_integrator(background_derivs,
                                  tau_start,
                                  tau_end,
                                  pvecback_integration,
                                  &bpaw,
                                  ppr->tol_background_integration,
                                  ppr->smallest_allowed_variation,
                                  &gi),
               gi.error_message,
               pba->error_message);

    /* -> store value of tau */
    pvecback_integration[pba->index_bi_tau]=tau_end;

  }

  /** - save last data in growTable with gt_add() */
  class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
             gTable.error_message,
             pba->error_message);
  pba->bt_size++;


  /* integration finished */

  /** - clean up generic integrator with cleanup_generic_integrator() */
  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             pba->error_message);

  /** - retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
             gTable.error_message,
             pba->error_message);

  /** - interpolate to get quantities precisely today with array_interpolate() */
  class_call(array_interpolate(
                               pData,
                               pba->bi_size,
                               pba->bt_size,
                               pba->index_bi_a,
                               pba->a_today,
                               &last_index,
                               pvecback_integration,
                               pba->bi_size,
                               pba->error_message),
             pba->error_message,
             pba->error_message);

  /* substitute last line with quantities today */
  for (i=0; i<pba->bi_size; i++)
    pData[(pba->bt_size-1)*pba->bi_size+i]=pvecback_integration[i];

  /** - deduce age of the Universe */
  /* -> age in Gyears */
  pba->age = pvecback_integration[pba->index_bi_time]/_Gyr_over_Mpc_;
  /* -> conformal age in Mpc */
  pba->conformal_age = pvecback_integration[pba->index_bi_tau];
  /* -> contribution of decaying dark matter and dark radiation to the critical density today: */
  if (pba->has_dcdm == _TRUE_){
    pba->Omega0_dcdm = pvecback_integration[pba->index_bi_rho_dcdm]/pba->H0/pba->H0;
  }
  if (pba->has_dr == _TRUE_){
    pba->Omega0_dr = pvecback_integration[pba->index_bi_rho_dr]/pba->H0/pba->H0;
  }


  /** - allocate background tables */
  class_alloc(pba->tau_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->z_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2tau_dz2_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->background_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2background_dtau2_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  /** - In a loop over lines, fill background table using the result of the integration plus background_functions() */
  for (i=0; i < pba->bt_size; i++) {

    /* -> establish correspondence between the integrated variable and the bg variables */

    pba->tau_table[i] = pData[i*pba->bi_size+pba->index_bi_tau];

    class_test(pData[i*pba->bi_size+pba->index_bi_a] <= 0.,
               pba->error_message,
               "a = %e instead of strictly positiv",pData[i*pba->bi_size+pba->index_bi_a]);

    pba->z_table[i] = pba->a_today/pData[i*pba->bi_size+pba->index_bi_a]-1.;

    pvecback[pba->index_bg_time] = pData[i*pba->bi_size+pba->index_bi_time];
    pvecback[pba->index_bg_conf_distance] = pba->conformal_age - pData[i*pba->bi_size+pba->index_bi_tau];

    if (pba->sgnK == 0) comoving_radius = pvecback[pba->index_bg_conf_distance];
    else if (pba->sgnK == 1) comoving_radius = sin(sqrt(pba->K)*pvecback[pba->index_bg_conf_distance])/sqrt(pba->K);
    else if (pba->sgnK == -1) comoving_radius = sinh(sqrt(-pba->K)*pvecback[pba->index_bg_conf_distance])/sqrt(-pba->K);

    pvecback[pba->index_bg_ang_distance] = pba->a_today*comoving_radius/(1.+pba->z_table[i]);
    pvecback[pba->index_bg_lum_distance] = pba->a_today*comoving_radius*(1.+pba->z_table[i]);
    pvecback[pba->index_bg_rs] = pData[i*pba->bi_size+pba->index_bi_rs];

    /* -> compute all other quantities depending only on {B} variables.
       The value of {B} variables in pData are also copied to pvecback.*/
    class_call(background_functions(pba,pData+i*pba->bi_size, pba->long_info, pvecback),
               pba->error_message,
               pba->error_message);

    /* -> compute growth functions (valid in dust universe) */

    /* D = H \int [da/(aH)^3] = H \int [dtau/(aH^2)] = H * growth */
    pvecback[pba->index_bg_D] = pvecback[pba->index_bg_H]*pData[i*pba->bi_size+pba->index_bi_growth];

    /* f = [dlnD]/[dln a] = 1/(aH) [dlnD]/[dtau] = H'/(aH^2) + 1/(a^2 H^3 growth) */
    pvecback[pba->index_bg_f] = pvecback[pba->index_bg_H_prime]/pvecback[pba->index_bg_a]/pvecback[pba->index_bg_H]/pvecback[pba->index_bg_H] + 1./(pvecback[pba->index_bg_a]*pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]*pData[i*pba->bi_size+pba->index_bi_growth]);

   /* D =  a * growth in wCDM*/
    pvecback[pba->index_bg_D_wCDM] = pvecback[pba->index_bg_a]*pData[i*pba->bi_size+pba->index_bi_growth_wCDM];

    //printf("a = %e, D/DwCDM = %e, DwCDM = %e\n",pvecback[pba->index_bg_a],pvecback[pba->index_bg_D]/pvecback[pba->index_bg_D_wCDM],pvecback[pba->index_bg_D_wCDM]);
    

    /* -> write in the table */
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size,pvecback,pba->bg_size*sizeof(double));

    class_test(memcopy_result != pba->background_table + i*pba->bg_size,
               pba->error_message,
               "cannot copy data back to pba->background_table");
  }

  /** - free the growTable with gt_free() */

  class_call(gt_free(&gTable),
             gTable.error_message,
             pba->error_message);

  /** - fill tables of second derivatives (in view of spline interpolation) */
  class_call(array_spline_table_lines(pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      1,
                                      pba->d2tau_dz2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  class_call(array_spline_table_lines(pba->tau_table,
                                      pba->bt_size,
                                      pba->background_table,
                                      pba->bg_size,
                                      pba->d2background_dtau2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  /** - compute remaining "related parameters" 
   *     - so-called "effective neutrino number", computed at earliest
      time in interpolation table. This should be seen as a
      definition: Neff is the equivalent number of
      instantaneously-decoupled neutrinos accounting for the
      radiation density, beyond photons */
  pba->Neff = (pba->background_table[pba->index_bg_Omega_r]
               *pba->background_table[pba->index_bg_rho_crit]
               -pba->background_table[pba->index_bg_rho_g])
    /(7./8.*pow(4./11.,4./3.)*pba->background_table[pba->index_bg_rho_g]);

  /** - done */
  if (pba->background_verbose > 0) {
    printf(" -> age = %f Gyr\n",pba->age);
    printf(" -> conformal age = %f Mpc\n",pba->conformal_age);
  }

  if (pba->background_verbose > 2) {
    if ((pba->has_dcdm == _TRUE_)&&(pba->has_dr == _TRUE_)){
      printf("    Decaying Cold Dark Matter details: (DCDM --> DR)\n");
      printf("     -> Omega0_dcdm = %f\n",pba->Omega0_dcdm);
      printf("     -> Omega0_dr = %f\n",pba->Omega0_dr);
      printf("     -> Omega0_dr+Omega0_dcdm = %f, input value = %f\n",
             pba->Omega0_dr+pba->Omega0_dcdm,pba->Omega0_dcdmdr);
      printf("     -> Omega_ini_dcdm/Omega_b = %f\n",pba->Omega_ini_dcdm/pba->Omega0_b);
    }
    if (pba->has_scf == _TRUE_){
      printf("    Scalar field details:\n");
      printf("     -> Omega_scf = %g, wished %g\n",
             pvecback[pba->index_bg_rho_scf]/pvecback[pba->index_bg_rho_crit], pba->Omega0_scf);
      if(pba->has_lambda == _TRUE_)
	printf("     -> Omega_Lambda = %g, wished %g\n",
               pvecback[pba->index_bg_rho_lambda]/pvecback[pba->index_bg_rho_crit], pba->Omega0_lambda);
      printf("     -> parameters: [lambda, alpha, A, B] = \n");
      printf("                    [");
      for (i=0; i<pba->scf_parameters_size-1; i++){
        printf("%.3f, ",pba->scf_parameters[i]);
      }
      printf("%.3f]\n",pba->scf_parameters[pba->scf_parameters_size-1]);
    }
  }

  free(pvecback);
  free(pvecback_integration);

  return _SUCCESS_;

}

/**
 * Assign initial values to background integrated variables.
 *
 * @param ppr                  Input: pointer to precision structure
 * @param pba                  Input: pointer to background structure
 * @param pvecback             Input: vector of background quantities used as workspace
 * @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
 * @return the error status
 */

int background_initial_conditions(
                                  struct precision *ppr,
                                  struct background *pba,
                                  double * pvecback, /* vector with argument pvecback[index_bg] (must be already allocated, normal format is sufficient) */
                                  double * pvecback_integration /* vector with argument pvecback_integration[index_bi] (must be already allocated with size pba->bi_size) */
                                  ) {

  /** Summary: */

  /** - define local variables */

  /* scale factor */
  double a;

  double rho_ncdm, p_ncdm, rho_ncdm_rel_tot=0.;
  double f,Omega_rad, rho_rad;
  int counter,is_early_enough,n_ncdm;
  double scf_lambda;

  /** - fix initial value of \f$ a \f$ */
  a = ppr->a_ini_over_a_today_default * pba->a_today;

  /**  If we have ncdm species, perhaps we need to start earlier
      than the standard value for the species to be relativistic.
      This could happen for some WDM models.
  */

  if (pba->has_ncdm == _TRUE_) {

    for (counter=0; counter < _MAX_IT_; counter++) {

      is_early_enough = _TRUE_;
      rho_ncdm_rel_tot = 0.;

      for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++) {

	class_call(background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
					   pba->w_ncdm_bg[n_ncdm],
					   pba->q_size_ncdm_bg[n_ncdm],
					   pba->M_ncdm[n_ncdm],
					   pba->factor_ncdm[n_ncdm],
					   pba->a_today/a-1.0,
					   NULL,
					   &rho_ncdm,
					   &p_ncdm,
					   NULL,
					   NULL),
                   pba->error_message,
                   pba->error_message);
	rho_ncdm_rel_tot += 3.*p_ncdm;
	if (fabs(p_ncdm/rho_ncdm-1./3.)>ppr->tol_ncdm_initial_w)
	  is_early_enough = _FALSE_;
      }
      if (is_early_enough == _TRUE_)
	break;
      else
	a *= _SCALE_BACK_;
    }
    class_test(counter == _MAX_IT_,
	       pba->error_message,
	       "Search for initial scale factor a such that all ncdm species are relativistic failed.");
  }

  pvecback_integration[pba->index_bi_a] = a;

  /* Set initial values of {B} variables: */
  Omega_rad = pba->Omega0_g;
  if (pba->has_ur == _TRUE_)
    Omega_rad += pba->Omega0_ur;
  rho_rad = Omega_rad*pow(pba->H0,2)/pow(a/pba->a_today,4);
  if (pba->has_ncdm == _TRUE_){
    /** - We must add the relativistic contribution from NCDM species */
    rho_rad += rho_ncdm_rel_tot;
  }
  if (pba->has_dcdm == _TRUE_){
    /* Remember that the critical density today in CLASS conventions is H0^2 */
    pvecback_integration[pba->index_bi_rho_dcdm] =
      pba->Omega_ini_dcdm*pba->H0*pba->H0*pow(pba->a_today/a,3);
    if (pba->background_verbose > 3)
      printf("Density is %g. a_today=%g. Omega_ini=%g\n",pvecback_integration[pba->index_bi_rho_dcdm],pba->a_today,pba->Omega_ini_dcdm);
  }

  if (pba->has_dr == _TRUE_){
    if (pba->has_dcdm == _TRUE_){
      /**  - f is the critical density fraction of DR. The exact solution is:
       * 
       * `f = -Omega_rad+pow(pow(Omega_rad,3./2.)+0.5*pow(a/pba->a_today,6)*pvecback_integration[pba->index_bi_rho_dcdm]*pba->Gamma_dcdm/pow(pba->H0,3),2./3.);`
       * 
       * but it is not numerically stable for very small f which is always the case.
       * Instead we use the Taylor expansion of this equation, which is equivalent to
       * ignoring f(a) in the Hubble rate.
       */
      f = 1./3.*pow(a/pba->a_today,6)*pvecback_integration[pba->index_bi_rho_dcdm]*pba->Gamma_dcdm/pow(pba->H0,3)/sqrt(Omega_rad);
      pvecback_integration[pba->index_bi_rho_dr] = f*pba->H0*pba->H0/pow(a/pba->a_today,4);
    }
    else{
      /** There is also a space reserved for a future case where dr is not sourced by dcdm */
      pvecback_integration[pba->index_bi_rho_dr] = 0.0;
    }
  }

  /** - Fix initial value of \f$ \phi, \phi' \f$
   * set directly in the radiation attractor => fixes the units in terms of rho_ur
   * 
   * TODO: 
   * - There seems to be some small oscillation when it starts. 
   * - Check equations and signs. Sign of phi_prime? 
   * - is rho_ur all there is early on?
   */
  if(pba->has_scf == _TRUE_){
    scf_lambda = pba->scf_parameters[0];
    if(pba->attractor_ic_scf == _TRUE_){
      pvecback_integration[pba->index_bi_phi_scf] = -1/scf_lambda*
        log(rho_rad*4./(3*pow(scf_lambda,2)-12))*pba->phi_ini_scf;
      if (3.*pow(scf_lambda,2)-12. < 0){
        /** - --> If there is no attractor solution for scf_lambda, assign some value. Otherwise would give a nan.*/
    	pvecback_integration[pba->index_bi_phi_scf] = 1./scf_lambda;//seems to the work
	if (pba->background_verbose > 0)
	  printf(" No attractor IC for lambda = %.3e ! \n ",scf_lambda);
      }
      pvecback_integration[pba->index_bi_phi_prime_scf] = 2*pvecback_integration[pba->index_bi_a]*
        sqrt(V_scf(pba,pvecback_integration[pba->index_bi_phi_scf]))*pba->phi_prime_ini_scf;
    }
    else{
      printf("Not using attractor initial conditions\n");
      /** - --> If no attractor initial conditions are assigned, gets the provided ones. */
      pvecback_integration[pba->index_bi_phi_scf] = pba->phi_ini_scf;
      pvecback_integration[pba->index_bi_phi_prime_scf] = pba->phi_prime_ini_scf;
    }
    class_test(!isfinite(pvecback_integration[pba->index_bi_phi_scf]) ||
               !isfinite(pvecback_integration[pba->index_bi_phi_scf]),
               pba->error_message,
               "initial phi = %e phi_prime = %e -> check initial conditions",
               pvecback_integration[pba->index_bi_phi_scf],
               pvecback_integration[pba->index_bi_phi_scf]);
  }

  
  if (pba->has_ds == _TRUE_){
  pvecback_integration[pba->index_bi_m2_hi] = pba->m2_hi_ini;

      if (pba->mg_type == fR ){
    pvecback_integration[pba->index_bi_f_ds] = pba->f_ini_ds;
    pvecback_integration[pba->index_bi_f_prime_ds] = pba->f_prime_ini_ds;
 }
 }
  
  /* Infer pvecback from pvecback_integration */
  class_call(background_functions(pba, pvecback_integration, pba->normal_info, pvecback),
	     pba->error_message,
	     pba->error_message);

  /* Just checking that our initial time indeed is deep enough in the radiation
     dominated regime */
  class_test(fabs(pvecback[pba->index_bg_Omega_r]-1.) > ppr->tol_initial_Omega_r,
	     pba->error_message,
	     "Omega_r = %e, not close enough to 1. Decrease a_ini_over_a_today_default in order to start from radiation domination.",
	     pvecback[pba->index_bg_Omega_r]);

  /** - compute initial proper time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ t=1/(2H) \f$ (good
      approximation for most purposes) */

  class_test(pvecback[pba->index_bg_H] <= 0.,
             pba->error_message,
             "H = %e instead of strictly positive",pvecback[pba->index_bg_H]);

  pvecback_integration[pba->index_bi_time] = 1./(2.* pvecback[pba->index_bg_H]);

  /** - compute initial conformal time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ \tau=1/(aH) \f$
      (good approximation for most purposes) */
  pvecback_integration[pba->index_bi_tau] = 1./(a * pvecback[pba->index_bg_H]);

  /** - compute initial sound horizon, assuming \f$ c_s=1/\sqrt{3} \f$ initially */
  pvecback_integration[pba->index_bi_rs] = pvecback_integration[pba->index_bi_tau]/sqrt(3.);

  /** - compute initial value of the integral over \f$ d\tau /(aH^2) \f$,
      assumed to be proportional to \f$ a^4 \f$ during RD, but with arbitrary
      normalization */
  pvecback_integration[pba->index_bi_growth] = 1./(4.*a*a*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]);


pvecback_integration[pba->index_bi_growth_wCDM] = 1.;
    pvecback_integration[pba->index_bi_growth_prime_wCDM] = 0.;

  
  return _SUCCESS_;

}

/**
 * Subroutine for formatting background output
 * 
 */

int background_output_titles(struct background * pba,
                             char titles[_MAXTITLESTRINGLENGTH_]
                             ){

  /** - Length of the column title should be less than _OUTPUTPRECISION_+6
      to be indented correctly, but it can be as long as . */
  int n;
  char tmp[20];
  
  class_store_columntitle(titles,"z",_TRUE_);
  class_store_columntitle(titles,"proper time [Gyr]",_TRUE_);
  class_store_columntitle(titles,"conf. time [Mpc]",_TRUE_);
  class_store_columntitle(titles,"H [1/Mpc]",_TRUE_);
 class_store_columntitle(titles,"gr.fac. D_wCDM",_TRUE_);
  class_store_columntitle(titles,"gr.fac. D",_TRUE_);

  
  class_store_columntitle(titles,"Omega_ds",pba->has_ds);
if (pba->mg_type == fR) {
  class_store_columntitle(titles,"B_ds",pba->has_ds);
  class_store_columntitle(titles,"fR_ds",pba->has_ds);


   }
  
  /*  
  class_store_columntitle(titles,"comov. dist.",_TRUE_);
  class_store_columntitle(titles,"ang.diam.dist.",_TRUE_);
  class_store_columntitle(titles,"lum. dist.",_TRUE_);
  class_store_columntitle(titles,"comov.snd.hrz.",_TRUE_);
  class_store_columntitle(titles,"(.)rho_g",_TRUE_);
  class_store_columntitle(titles,"(.)rho_b",_TRUE_);
  class_store_columntitle(titles,"(.)rho_cdm",pba->has_cdm);
  if (pba->has_ncdm == _TRUE_){
    for (n=0; n<pba->N_ncdm; n++){
      sprintf(tmp,"(.)rho_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
      sprintf(tmp,"(.)p_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
    }
  }
  class_store_columntitle(titles,"(.)rho_lambda",pba->has_lambda);
  class_store_columntitle(titles,"(.)rho_fld",pba->has_fld);
  class_store_columntitle(titles,"(.)rho_ur",pba->has_ur);
  class_store_columntitle(titles,"(.)rho_crit",_TRUE_);
  class_store_columntitle(titles,"(.)rho_dcdm",pba->has_dcdm);
  class_store_columntitle(titles,"(.)rho_dr",pba->has_dr);

  class_store_columntitle(titles,"(.)rho_scf",pba->has_scf);
  class_store_columntitle(titles,"(.)p_scf",pba->has_scf);
  class_store_columntitle(titles,"phi_scf",pba->has_scf);
  class_store_columntitle(titles,"phi'_scf",pba->has_scf);
  class_store_columntitle(titles,"V_scf",pba->has_scf);
  class_store_columntitle(titles,"V'_scf",pba->has_scf);
  class_store_columntitle(titles,"V''_scf",pba->has_scf);

  class_store_columntitle(titles,"gr.fac. D",_TRUE_);
  class_store_columntitle(titles,"gr.fac. f",_TRUE_);
  */
  class_store_columntitle(titles,"alpha_hi",pba->has_ds);
  class_store_columntitle(titles,"alpha_hi",pba->has_ds);
  
  class_store_columntitle(titles,"alpha_hi",pba->has_ds);
  class_store_columntitle(titles,"m2_hi",pba->has_ds);
  class_store_columntitle(titles,"alpha_m_hi",pba->has_ds);
  class_store_columntitle(titles,"alpha_b_hi",pba->has_ds);
  class_store_columntitle(titles,"alpha_t_hi",pba->has_ds);
  class_store_columntitle(titles,"alpha_k_hi",pba->has_ds);
  class_store_columntitle(titles,"alpha_b_hi_prime",pba->has_ds);
  class_store_columntitle(titles,"cs2_hi",pba->has_ds);
  class_store_columntitle(titles,"Gamma_hi",pba->has_ds);
  class_store_columntitle(titles,"gamma_1_hi",pba->has_ds);
  class_store_columntitle(titles,"gamma_2_hi",pba->has_ds);
  class_store_columntitle(titles,"gamma_3_hi",pba->has_ds);
  class_store_columntitle(titles,"gamma_4_hi",pba->has_ds);
  class_store_columntitle(titles,"gamma_5_hi",pba->has_ds);
  class_store_columntitle(titles,"gamma_6_hi",pba->has_ds);
  class_store_columntitle(titles,"gamma_7_hi",pba->has_ds);
  class_store_columntitle(titles,"gamma_7_not_tilde_hi",pba->has_ds);
  class_store_columntitle(titles,"gamma_8_hi",pba->has_ds);
  class_store_columntitle(titles,"gamma_9_hi",pba->has_ds);
  class_store_columntitle(titles,"gamma_10_hi",pba->has_ds);




  return _SUCCESS_;
}

int background_output_data(
                           struct background *pba,
                           int number_of_titles,
                           double *data){
  int index_tau, storeidx, n;
  double *dataptr, *pvecback;

  /** Stores quantities */
  for (index_tau=0; index_tau<pba->bt_size; index_tau++){
    dataptr = data + index_tau*number_of_titles;
    pvecback = pba->background_table + index_tau*pba->bg_size;
    storeidx = 0;

    class_store_double(dataptr,pba->a_today/pvecback[pba->index_bg_a]-1.,_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_time]/_Gyr_over_Mpc_,_TRUE_,storeidx);
    class_store_double(dataptr,pba->conformal_age-pvecback[pba->index_bg_conf_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_H],_TRUE_,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_D_wCDM],_TRUE_,storeidx);
     class_store_double(dataptr,pvecback[pba->index_bg_D],_TRUE_,storeidx);

        class_store_double(dataptr,pvecback[pba->index_bg_Omega_ds],pba->has_ds,storeidx);

 if (pba->mg_type == fR) {
    class_store_double(dataptr,pvecback[pba->index_bg_B_ds],pba->has_ds,storeidx);
      class_store_double(dataptr,pvecback[pba->index_bg_fR_ds],pba->has_ds,storeidx);  }
    
    /*
    class_store_double(dataptr,pvecback[pba->index_bg_conf_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_ang_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_lum_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rs],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_g],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_b],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_cdm],pba->has_cdm,storeidx);
    if (pba->has_ncdm == _TRUE_){
      for (n=0; n<pba->N_ncdm; n++){
        class_store_double(dataptr,pvecback[pba->index_bg_rho_ncdm1+n],_TRUE_,storeidx);
        class_store_double(dataptr,pvecback[pba->index_bg_p_ncdm1+n],_TRUE_,storeidx);
      }
    }
    class_store_double(dataptr,pvecback[pba->index_bg_rho_lambda],pba->has_lambda,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_fld],pba->has_fld,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_ur],pba->has_ur,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_crit],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_dcdm],pba->has_dcdm,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_dr],pba->has_dr,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_rho_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_p_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_phi_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_phi_prime_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_V_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_dV_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_ddV_scf],pba->has_scf,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_f],_TRUE_,storeidx);
    */

    class_store_double(dataptr,pvecback[pba->index_bg_alpha_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_m2_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_alpha_m_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_alpha_b_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_alpha_t_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_alpha_k_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_alpha_b_hi_prime],pba->has_ds,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_cs2_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_Gamma_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_gamma_1_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_gamma_2_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_gamma_3_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_gamma_4_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_gamma_5_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_gamma_6_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_gamma_7_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_gamma_7_not_tilde_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_gamma_8_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_gamma_9_hi],pba->has_ds,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_gamma_10_hi],pba->has_ds,storeidx);


  }

  return _SUCCESS_;
}


/**
 * Subroutine evaluating the derivative with respect to conformal time
 * of quantities which are integrated (a, t, etc).
 *
 * This is one of the few functions in the code which is passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed input parameters and workspaces are passed through a generic
 * pointer. Here, this is just a pointer to the background structure
 * and to a background vector, but generic_integrator() doesn't know
 * its fine structure.
 *
 * - the error management is a bit special: errors are not written as
 * usual to pba->error_message, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param tau                      Input: conformal time
 * @param y                        Input: vector of variable
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output: error message
 */
int background_derivs(
                      double tau,
                      double* y, /* vector with argument y[index_bi] (must be already allocated with size pba->bi_size) */
                      double* dy, /* vector with argument dy[index_bi]
                                     (must be already allocated with
                                     size pba->bi_size) */
                      void * parameters_and_workspace,
                      ErrorMsg error_message
                      ) {

  /** Summary: */

  /** - define local variables */

  struct background_parameters_and_workspace * pbpaw;
  struct background * pba;
  double * pvecback;

  pbpaw = parameters_and_workspace;
  pba =  pbpaw->pba;
  pvecback = pbpaw->pvecback;

  /** - calculate functions of \f$ a \f$ with background_functions() */
  class_call(background_functions(pba, y, pba->normal_info, pvecback),
             pba->error_message,
             error_message);
/* ---------F.P. 26-10-2025 ---------
double a = y[pba->index_bi_a];
double H0 = pba->H0;
double b = pba->b;
int l = pba->l;		//f(R) mod gravity index

// ** RICALCOLA RHO_TOT e P_TOT sommando le componenti che background_functions ha salvato **

double rho_tot = 0.;
double p_tot = 0.;

rho_tot += pvecback[pba->index_bg_rho_g];
rho_tot += pvecback[pba->index_bg_rho_b];
rho_tot += pvecback[pba->index_bg_rho_cdm];
rho_tot += pvecback[pba->index_bg_rho_lambda];

p_tot += (1./3.) * pvecback[pba->index_bg_rho_g]; // Pressione fotoni


double H2; 
double H_prime;
// dichiarazione delle derivate lcdm utile per f(R)
/*double HL = 0;
double H1L = 0;
double H2L = 0;
double H3L = 0;
double H4L = 0;
double H5L = 0; 

 
// switch for choose the f(R) Hubble parameter:

  if (pba->l == 0){		//lcdm
      		printf("DEBUG: BACKGROUND MODEL = LCDM (Case 0) SELECTED\n");
  		H2 = rho_tot-pba->K/a/a;		// lcdm 
  		H_prime = - (3./2.) * (rho_tot + p_tot) * a + pba->K/a;
  		}
  		
  		
  if (pba->l == 1){		//HS1
  		// -------- Hubble parameter lcdm and his 5 derivatives useful for f(R) Hubble parameter---------------
		printf("DEBUG: BACKGROUND MODEL = HS1 (Case 1) SELECTED\n");
 		double HL = H0 * sqrt((pba->Omega0_b + pba->Omega0_cdm) /a/a/a + pba->Omega0_lambda);
  
  		double H1L =-3*H0*(pba->Omega0_b + pba->Omega0_cdm)/a/a/a/(2*sqrt(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a));
  
  		double H2L = -9*H0*(pba->Omega0_b + pba->Omega0_cdm)*((pba->Omega0_b + pba->Omega0_cdm)/a/a/a/(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a) - 2)/a/a/a/(4*sqrt(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a));
  
  		double H3L = -27*H0*(pba->Omega0_b + pba->Omega0_cdm)*(3*(pba->Omega0_b + pba->Omega0_cdm)*(pba->Omega0_b + pba->Omega0_cdm)/a/a/a/a/a/a/((pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)*(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)) - 6*(pba->Omega0_b + pba->Omega0_cdm)/a/a/a/(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a) + 4)/a/a/a/(8*sqrt(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a));
  
 		double H4L = -81*H0*(pba->Omega0_b + pba->Omega0_cdm)*(15*(pba->Omega0_b + pba->Omega0_cdm)*(pba->Omega0_b + pba->Omega0_cdm)*(pba->Omega0_b + pba->Omega0_cdm)/a/a/a/a/a/a/a/a/a/((pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)*(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)*(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)) - 36*(pba->Omega0_b + pba->Omega0_cdm)*(pba->Omega0_b + pba->Omega0_cdm)/a/a/a/a/a/a/((pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)*(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)) + 28*(pba->Omega0_b + pba->Omega0_cdm)/a/a/a/(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a) - 8)/a/a/a/(16*sqrt(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a));
  
  		double H5L = -243*H0*(pba->Omega0_b + pba->Omega0_cdm)*(105*(pba->Omega0_b + pba->Omega0_cdm)*(pba->Omega0_b + pba->Omega0_cdm)*(pba->Omega0_b + pba->Omega0_cdm)*(pba->Omega0_b + pba->Omega0_cdm)/a/a/a/a/a/a/a/a/a/a/a/a/((pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)*(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)*(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)*(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)) - 300*(pba->Omega0_b + pba->Omega0_cdm)*(pba->Omega0_b + pba->Omega0_cdm)*(pba->Omega0_b + pba->Omega0_cdm)/a/a/a/a/a/a/a/a/a/((pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)*(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)*(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a))+ 300*(pba->Omega0_b + pba->Omega0_cdm)*(pba->Omega0_b + pba->Omega0_cdm)/a/a/a/a/a/a/((pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)*(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a)) - 120*(pba->Omega0_b + pba->Omega0_cdm)/a/a/a/(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a) + 16)/a/a/a/(32*sqrt(pba->Omega0_lambda + (pba->Omega0_b + pba->Omega0_cdm)/a/a/a));
  		  //-----------------------------------------------------------------------------------------------------------
  		
  		H2 = 	(1.0/4.0)*pow(H0, 6)*pow(pba->Omega0_lambda, 3)*pow(b, 2)*(32*pow(H0, 2)*pow(H1L, 6)*pba->Omega0_lambda + 278*pow(H0, 2)*pow(H1L, 5)*HL*pba->Omega0_lambda + 44*pow(H0, 2)*pow(H1L, 4)*H2L*HL*pba->Omega0_lambda + 749*pow(H0, 2)*pow(H1L, 4)*pow(HL, 2)*pba->Omega0_lambda + 296*pow(H0, 2)*pow(H1L, 3)*H2L*pow(HL, 2)*pba->Omega0_lambda - 12*pow(H0, 2)*pow(H1L, 3)*H3L*pow(HL, 2)*pba->Omega0_lambda + 168*pow(H0, 2)*pow(H1L, 3)*pow(HL, 3)*pba->Omega0_lambda + 48*pow(H0, 2)*pow(H1L, 2)*pow(H2L, 2)*pow(HL, 2)*pba->Omega0_lambda + 434*pow(H0, 2)*pow(H1L, 2)*H2L*pow(HL, 3)*pba->Omega0_lambda - 84*pow(H0, 2)*pow(H1L, 2)*H3L*pow(HL, 3)*pba->Omega0_lambda + 2*pow(H0, 2)*pow(H1L, 2)*H4L*pow(HL, 3)*pba->Omega0_lambda - 2160*pow(H0, 2)*pow(H1L, 2)*pow(HL, 4)*pba->Omega0_lambda + 246*pow(H0, 2)*H1L*pow(H2L, 2)*pow(HL, 3)*pba->Omega0_lambda - 24*pow(H0, 2)*H1L*H2L*H3L*pow(HL, 3)*pba->Omega0_lambda - 1032*pow(H0, 2)*H1L*H2L*pow(HL, 4)*pba->Omega0_lambda - 96*pow(H0, 2)*H1L*H3L*pow(HL, 4)*pba->Omega0_lambda + 8*pow(H0, 2)*H1L*H4L*pow(HL, 4)*pba->Omega0_lambda - 1152*pow(H0, 2)*H1L*pow(HL, 5)*pba->Omega0_lambda + 42*pow(H0, 2)*pow(H2L, 3)*pow(HL, 3)*pba->Omega0_lambda - 204*pow(H0, 2)*pow(H2L, 2)*pow(HL, 4)*pba->Omega0_lambda - 48*pow(H0, 2)*H2L*H3L*pow(HL, 4)*pba->Omega0_lambda - 120*pow(H0, 2)*H2L*pow(HL, 5)*pba->Omega0_lambda + 48*pow(H0, 2)*H3L*pow(HL, 5)*pba->Omega0_lambda + 8*pow(H0, 2)*H4L*pow(HL, 5)*pba->Omega0_lambda - 144*pow(H0, 2)*pow(HL, 6)*pba->Omega0_lambda + 9*pow(H1L, 6)*pow(HL, 2) + 106*pow(H1L, 5)*pow(HL, 3) + 6*pow(H1L, 4)*H2L*pow(HL, 3) + 496*pow(H1L, 4)*pow(HL, 4) + 48*pow(H1L, 3)*H2L*pow(HL, 4) + 1168*pow(H1L, 3)*pow(HL, 5) + 144*pow(H1L, 2)*H2L*pow(HL, 5) + 1424*pow(H1L, 2)*pow(HL, 6) + 192*H1L*H2L*pow(HL, 6) + 800*H1L*pow(HL, 7) + 96*H2L*pow(HL, 7) + 128*pow(HL, 8))/(pow(HL, 4)*pow(H1L + 2*HL, 8)) + pow(H0, 4)*pow(pba->Omega0_lambda, 2)*b*(-2.0*pow(H1L, 2) - 7.5*H1L*HL - H2L*HL - 3.0*pow(HL, 2))/(HL*(pow(H1L, 3) + 6.0*pow(H1L, 2)*HL + 12.0*H1L*pow(HL, 2) + 8.0*pow(HL, 3))) + pow(HL, 2);
  		
  		H_prime = ( (-1.0/2.0*pow(H0, 6)*H1L*pow(pba->Omega0_lambda, 3)*pow(b, 2)*(32*pow(H0, 2)*pow(H1L, 6)*pba->Omega0_lambda + 278*pow(H0, 2)*pow(H1L, 5)*HL*pba->Omega0_lambda + 44*pow(H0, 2)*pow(H1L, 4)*H2L*HL*pba->Omega0_lambda + 749*pow(H0, 2)*pow(H1L, 4)*pow(HL, 2)*pba->Omega0_lambda + 296*pow(H0, 2)*pow(H1L, 3)*H2L*pow(HL, 2)*pba->Omega0_lambda - 12*pow(H0, 2)*pow(H1L, 3)*H3L*pow(HL, 2)*pba->Omega0_lambda + 168*pow(H0, 2)*pow(H1L, 3)*pow(HL, 3)*pba->Omega0_lambda + 48*pow(H0, 2)*pow(H1L, 2)*pow(H2L, 2)*pow(HL, 2)*pba->Omega0_lambda + 434*pow(H0, 2)*pow(H1L, 2)*H2L*pow(HL, 3)*pba->Omega0_lambda - 84*pow(H0, 2)*pow(H1L, 2)*H3L*pow(HL, 3)*pba->Omega0_lambda + 2*pow(H0, 2)*pow(H1L, 2)*H4L*pow(HL, 3)*pba->Omega0_lambda - 2160*pow(H0, 2)*pow(H1L, 2)*pow(HL, 4)*pba->Omega0_lambda + 246*pow(H0, 2)*H1L*pow(H2L, 2)*pow(HL, 3)*pba->Omega0_lambda - 24*pow(H0, 2)*H1L*H2L*H3L*pow(HL, 3)*pba->Omega0_lambda - 1032*pow(H0, 2)*H1L*H2L*pow(HL, 4)*pba->Omega0_lambda - 96*pow(H0, 2)*H1L*H3L*pow(HL, 4)*pba->Omega0_lambda + 8*pow(H0, 2)*H1L*H4L*pow(HL, 4)*pba->Omega0_lambda - 1152*pow(H0, 2)*H1L*pow(HL, 5)*pba->Omega0_lambda + 42*pow(H0, 2)*pow(H2L, 3)*pow(HL, 3)*pba->Omega0_lambda - 204*pow(H0, 2)*pow(H2L, 2)*pow(HL, 4)*pba->Omega0_lambda - 48*pow(H0, 2)*H2L*H3L*pow(HL, 4)*pba->Omega0_lambda - 120*pow(H0, 2)*H2L*pow(HL, 5)*pba->Omega0_lambda + 48*pow(H0, 2)*H3L*pow(HL, 5)*pba->Omega0_lambda + 8*pow(H0, 2)*H4L*pow(HL, 5)*pba->Omega0_lambda - 144*pow(H0, 2)*pow(HL, 6)*pba->Omega0_lambda + 9*pow(H1L, 6)*pow(HL, 2) + 106*pow(H1L, 5)*pow(HL, 3) + 6*pow(H1L, 4)*H2L*pow(HL, 3) + 496*pow(H1L, 4)*pow(HL, 4) + 48*pow(H1L, 3)*H2L*pow(HL, 4) + 1168*pow(H1L, 3)*pow(HL, 5) + 144*pow(H1L, 2)*H2L*pow(HL, 5) + 1424*pow(H1L, 2)*pow(HL, 6) + 192*H1L*H2L*pow(HL, 6) + 800*H1L*pow(HL, 7) + 96*H2L*pow(HL, 7) + 128*pow(HL, 8))/(pow(HL, 5)*pow(H1L + 2*HL, 8)) + (1.0/8.0)*pow(H0, 6)*pow(pba->Omega0_lambda, 3)*pow(b, 2)*(-16*H1L - 8*H2L)*(32*pow(H0, 2)*pow(H1L, 6)*pba->Omega0_lambda + 278*pow(H0, 2)*pow(H1L, 5)*HL*pba->Omega0_lambda + 44*pow(H0, 2)*pow(H1L, 4)*H2L*HL*pba->Omega0_lambda + 749*pow(H0, 2)*pow(H1L, 4)*pow(HL, 2)*pba->Omega0_lambda + 296*pow(H0, 2)*pow(H1L, 3)*H2L*pow(HL, 2)*pba->Omega0_lambda - 12*pow(H0, 2)*pow(H1L, 3)*H3L*pow(HL, 2)*pba->Omega0_lambda + 168*pow(H0, 2)*pow(H1L, 3)*pow(HL, 3)*pba->Omega0_lambda + 48*pow(H0, 2)*pow(H1L, 2)*pow(H2L, 2)*pow(HL, 2)*pba->Omega0_lambda + 434*pow(H0, 2)*pow(H1L, 2)*H2L*pow(HL, 3)*pba->Omega0_lambda - 84*pow(H0, 2)*pow(H1L, 2)*H3L*pow(HL, 3)*pba->Omega0_lambda + 2*pow(H0, 2)*pow(H1L, 2)*H4L*pow(HL, 3)*pba->Omega0_lambda - 2160*pow(H0, 2)*pow(H1L, 2)*pow(HL, 4)*pba->Omega0_lambda + 246*pow(H0, 2)*H1L*pow(H2L, 2)*pow(HL, 3)*pba->Omega0_lambda - 24*pow(H0, 2)*H1L*H2L*H3L*pow(HL, 3)*pba->Omega0_lambda - 1032*pow(H0, 2)*H1L*H2L*pow(HL, 4)*pba->Omega0_lambda - 96*pow(H0, 2)*H1L*H3L*pow(HL, 4)*pba->Omega0_lambda + 8*pow(H0, 2)*H1L*H4L*pow(HL, 4)*pba->Omega0_lambda - 1152*pow(H0, 2)*H1L*pow(HL, 5)*pba->Omega0_lambda + 42*pow(H0, 2)*pow(H2L, 3)*pow(HL, 3)*pba->Omega0_lambda - 204*pow(H0, 2)*pow(H2L, 2)*pow(HL, 4)*pba->Omega0_lambda - 48*pow(H0, 2)*H2L*H3L*pow(HL, 4)*pba->Omega0_lambda - 120*pow(H0, 2)*H2L*pow(HL, 5)*pba->Omega0_lambda + 48*pow(H0, 2)*H3L*pow(HL, 5)*pba->Omega0_lambda + 8*pow(H0, 2)*H4L*pow(HL, 5)*pba->Omega0_lambda - 144*pow(H0, 2)*pow(HL, 6)*pba->Omega0_lambda + 9*pow(H1L, 6)*pow(HL, 2) + 106*pow(H1L, 5)*pow(HL, 3) + 6*pow(H1L, 4)*H2L*pow(HL, 3) + 496*pow(H1L, 4)*pow(HL, 4) + 48*pow(H1L, 3)*H2L*pow(HL, 4) + 1168*pow(H1L, 3)*pow(HL, 5) + 144*pow(H1L, 2)*H2L*pow(HL, 5) + 1424*pow(H1L, 2)*pow(HL, 6) + 192*H1L*H2L*pow(HL, 6) + 800*H1L*pow(HL, 7) + 96*H2L*pow(HL, 7) + 128*pow(HL, 8))/(pow(HL, 4)*pow(H1L + 2*HL, 9)) + (1.0/8.0)*pow(H0, 6)*pow(pba->Omega0_lambda, 3)*pow(b, 2)*(278*pow(H0, 2)*pow(H1L, 6)*pba->Omega0_lambda + 236*pow(H0, 2)*pow(H1L, 5)*H2L*pba->Omega0_lambda + 1498*pow(H0, 2)*pow(H1L, 5)*HL*pba->Omega0_lambda + 1982*pow(H0, 2)*pow(H1L, 4)*H2L*HL*pba->Omega0_lambda + 20*pow(H0, 2)*pow(H1L, 4)*H3L*HL*pba->Omega0_lambda + 504*pow(H0, 2)*pow(H1L, 4)*pow(HL, 2)*pba->Omega0_lambda + 272*pow(H0, 2)*pow(H1L, 3)*pow(H2L, 2)*HL*pba->Omega0_lambda + 4298*pow(H0, 2)*pow(H1L, 3)*H2L*pow(HL, 2)*pba->Omega0_lambda + 44*pow(H0, 2)*pow(H1L, 3)*H3L*pow(HL, 2)*pba->Omega0_lambda - 6*pow(H0, 2)*pow(H1L, 3)*H4L*pow(HL, 2)*pba->Omega0_lambda - 8640*pow(H0, 2)*pow(H1L, 3)*pow(HL, 3)*pba->Omega0_lambda + 1626*pow(H0, 2)*pow(H1L, 2)*pow(H2L, 2)*pow(HL, 2)*pba->Omega0_lambda - 12*pow(H0, 2)*pow(H1L, 2)*H2L*H3L*pow(HL, 2)*pba->Omega0_lambda - 3624*pow(H0, 2)*pow(H1L, 2)*H2L*pow(HL, 3)*pba->Omega0_lambda + 50*pow(H0, 2)*pow(H1L, 2)*H3L*pow(HL, 3)*pba->Omega0_lambda - 52*pow(H0, 2)*pow(H1L, 2)*H4L*pow(HL, 3)*pba->Omega0_lambda + 2*pow(H0, 2)*pow(H1L, 2)*H5L*pow(HL, 3)*pba->Omega0_lambda - 5760*pow(H0, 2)*pow(H1L, 2)*pow(HL, 4)*pba->Omega0_lambda + 222*pow(H0, 2)*H1L*pow(H2L, 3)*pow(HL, 2)*pba->Omega0_lambda + 52*pow(H0, 2)*H1L*pow(H2L, 2)*pow(HL, 3)*pba->Omega0_lambda + 132*pow(H0, 2)*H1L*H2L*H3L*pow(HL, 3)*pba->Omega0_lambda - 20*pow(H0, 2)*H1L*H2L*H4L*pow(HL, 3)*pba->Omega0_lambda - 4920*pow(H0, 2)*H1L*H2L*pow(HL, 4)*pba->Omega0_lambda - 24*pow(H0, 2)*H1L*pow(H3L, 2)*pow(HL, 3)*pba->Omega0_lambda - 792*pow(H0, 2)*H1L*H3L*pow(HL, 4)*pba->Omega0_lambda - 56*pow(H0, 2)*H1L*H4L*pow(HL, 4)*pba->Omega0_lambda + 8*pow(H0, 2)*H1L*H5L*pow(HL, 4)*pba->Omega0_lambda - 864*pow(H0, 2)*H1L*pow(HL, 5)*pba->Omega0_lambda + 246*pow(H0, 2)*pow(H2L, 3)*pow(HL, 3)*pba->Omega0_lambda + 102*pow(H0, 2)*pow(H2L, 2)*H3L*pow(HL, 3)*pba->Omega0_lambda - 1032*pow(H0, 2)*pow(H2L, 2)*pow(HL, 4)*pba->Omega0_lambda - 504*pow(H0, 2)*H2L*H3L*pow(HL, 4)*pba->Omega0_lambda - 40*pow(H0, 2)*H2L*H4L*pow(HL, 4)*pba->Omega0_lambda - 1152*pow(H0, 2)*H2L*pow(HL, 5)*pba->Omega0_lambda - 48*pow(H0, 2)*pow(H3L, 2)*pow(HL, 4)*pba->Omega0_lambda - 120*pow(H0, 2)*H3L*pow(HL, 5)*pba->Omega0_lambda + 48*pow(H0, 2)*H4L*pow(HL, 5)*pba->Omega0_lambda + 8*pow(H0, 2)*H5L*pow(HL, 5)*pba->Omega0_lambda + 18*pow(H1L, 7)*HL + 318*pow(H1L, 6)*pow(HL, 2) + 72*pow(H1L, 5)*H2L*pow(HL, 2) + 1984*pow(H1L, 5)*pow(HL, 3) + 722*pow(H1L, 4)*H2L*pow(HL, 3) + 6*pow(H1L, 4)*H3L*pow(HL, 3) + 5840*pow(H1L, 4)*pow(HL, 4) + 24*pow(H1L, 3)*pow(H2L, 2)*pow(HL, 3) + 2704*pow(H1L, 3)*H2L*pow(HL, 4) + 48*pow(H1L, 3)*H3L*pow(HL, 4) + 8544*pow(H1L, 3)*pow(HL, 5) + 144*pow(H1L, 2)*pow(H2L, 2)*pow(HL, 4) + 4656*pow(H1L, 2)*H2L*pow(HL, 5) + 144*pow(H1L, 2)*H3L*pow(HL, 5) + 5600*pow(H1L, 2)*pow(HL, 6) + 288*H1L*pow(H2L, 2)*pow(HL, 5) + 3520*H1L*H2L*pow(HL, 6) + 192*H1L*H3L*pow(HL, 6) + 1024*H1L*pow(HL, 7) + 192*pow(H2L, 2)*pow(HL, 6) + 800*H2L*pow(HL, 7) + 96*H3L*pow(HL, 7))/(pow(HL, 4)*pow(H1L + 2*HL, 8)) - pow(H0, 4)*H1L*pow(pba->Omega0_lambda, 2)*b*(-2.0*pow(H1L, 2) - 7.5*H1L*HL - H2L*HL - 3.0*pow(HL, 2))/(pow(HL, 2)*(2*pow(H1L, 3) + 12.0*pow(H1L, 2)*HL + 24.0*H1L*pow(HL, 2) + 16.0*pow(HL, 3))) + 0.0034722222222222199*pow(H0, 4)*pow(pba->Omega0_lambda, 2)*b*(-2.0*pow(H1L, 2) - 7.5*H1L*HL - H2L*HL - 3.0*pow(HL, 2))*(-6.0*pow(H1L, 3) - 3*pow(H1L, 2)*H2L - 24.0*pow(H1L, 2)*HL - 12.0*H1L*H2L*HL - 24.0*H1L*pow(HL, 2) - 12.0*H2L*pow(HL, 2))/(HL*pow(0.083333333333333301*pow(H1L, 3) + 0.5*pow(H1L, 2)*HL + H1L*pow(HL, 2) + 0.66666666666666696*pow(HL, 3), 2)) + pow(H0, 4)*pow(pba->Omega0_lambda, 2)*b*(-7.5*pow(H1L, 2) - 5.0*H1L*H2L - 6.0*H1L*HL - 7.5*H2L*HL - H3L*HL)/(HL*(2*pow(H1L, 3) + 12.0*pow(H1L, 2)*HL + 24.0*H1L*pow(HL, 2) + 16.0*pow(HL, 3))) + H1L*HL)/sqrt((1.0/4.0)*pow(H0, 6)*pow(pba->Omega0_lambda, 3)*pow(b, 2)*(32*pow(H0, 2)*pow(H1L, 6)*pba->Omega0_lambda + 278*pow(H0, 2)*pow(H1L, 5)*HL*pba->Omega0_lambda + 44*pow(H0, 2)*pow(H1L, 4)*H2L*HL*pba->Omega0_lambda + 749*pow(H0, 2)*pow(H1L, 4)*pow(HL, 2)*pba->Omega0_lambda + 296*pow(H0, 2)*pow(H1L, 3)*H2L*pow(HL, 2)*pba->Omega0_lambda - 12*pow(H0, 2)*pow(H1L, 3)*H3L*pow(HL, 2)*pba->Omega0_lambda + 168*pow(H0, 2)*pow(H1L, 3)*pow(HL, 3)*pba->Omega0_lambda + 48*pow(H0, 2)*pow(H1L, 2)*pow(H2L, 2)*pow(HL, 2)*pba->Omega0_lambda + 434*pow(H0, 2)*pow(H1L, 2)*H2L*pow(HL, 3)*pba->Omega0_lambda - 84*pow(H0, 2)*pow(H1L, 2)*H3L*pow(HL, 3)*pba->Omega0_lambda + 2*pow(H0, 2)*pow(H1L, 2)*H4L*pow(HL, 3)*pba->Omega0_lambda - 2160*pow(H0, 2)*pow(H1L, 2)*pow(HL, 4)*pba->Omega0_lambda + 246*pow(H0, 2)*H1L*pow(H2L, 2)*pow(HL, 3)*pba->Omega0_lambda - 24*pow(H0, 2)*H1L*H2L*H3L*pow(HL, 3)*pba->Omega0_lambda - 1032*pow(H0, 2)*H1L*H2L*pow(HL, 4)*pba->Omega0_lambda - 96*pow(H0, 2)*H1L*H3L*pow(HL, 4)*pba->Omega0_lambda + 8*pow(H0, 2)*H1L*H4L*pow(HL, 4)*pba->Omega0_lambda - 1152*pow(H0, 2)*H1L*pow(HL, 5)*pba->Omega0_lambda + 42*pow(H0, 2)*pow(H2L, 3)*pow(HL, 3)*pba->Omega0_lambda - 204*pow(H0, 2)*pow(H2L, 2)*pow(HL, 4)*pba->Omega0_lambda - 48*pow(H0, 2)*H2L*H3L*pow(HL, 4)*pba->Omega0_lambda - 120*pow(H0, 2)*H2L*pow(HL, 5)*pba->Omega0_lambda + 48*pow(H0, 2)*H3L*pow(HL, 5)*pba->Omega0_lambda + 8*pow(H0, 2)*H4L*pow(HL, 5)*pba->Omega0_lambda - 144*pow(H0, 2)*pow(HL, 6)*pba->Omega0_lambda + 9*pow(H1L, 6)*pow(HL, 2) + 106*pow(H1L, 5)*pow(HL, 3) + 6*pow(H1L, 4)*H2L*pow(HL, 3) + 496*pow(H1L, 4)*pow(HL, 4) + 48*pow(H1L, 3)*H2L*pow(HL, 4) + 1168*pow(H1L, 3)*pow(HL, 5) + 144*pow(H1L, 2)*H2L*pow(HL, 5) + 1424*pow(H1L, 2)*pow(HL, 6) + 192*H1L*H2L*pow(HL, 6) + 800*H1L*pow(HL, 7) + 96*H2L*pow(HL, 7) + 128*pow(HL, 8))/(pow(HL, 4)*pow(H1L + 2*HL, 8)) + pow(H0, 4)*pow(pba->Omega0_lambda, 2)*b*(-2.0*pow(H1L, 2) - 7.5*H1L*HL - H2L*HL - 3.0*pow(HL, 2))/(HL*(pow(H1L, 3) + 6.0*pow(H1L, 2)*HL + 12.0*H1L*pow(HL, 2) + 8.0*pow(HL, 3))) + pow(HL, 2)) ); 
  		}
pvecback[pba->index_bg_H] = sqrt(H2);
pvecback[pba->index_bg_H_prime] = H_prime;


-------------------------------------*/



  /** - calculate \f$ a'=a^2 H \f$ */
  dy[pba->index_bi_a] = y[pba->index_bi_a] * y[pba->index_bi_a] * pvecback[pba->index_bg_H];

  /** - calculate \f$ t' = a \f$ */
  dy[pba->index_bi_time] = y[pba->index_bi_a];

  class_test(pvecback[pba->index_bg_rho_g] <= 0.,
             error_message,
             "rho_g = %e instead of strictly positive",pvecback[pba->index_bg_rho_g]);

  /** - calculate \f$ rs' = c_s \f$*/
  dy[pba->index_bi_rs] = 1./sqrt(3.*(1.+3.*pvecback[pba->index_bg_rho_b]/4./pvecback[pba->index_bg_rho_g]))*sqrt(1.-pba->K*y[pba->index_bi_rs]*y[pba->index_bi_rs]); // TBC: curvature correction

  /** - calculate growth' \f$ = 1/(aH^2) \f$ */
  dy[pba->index_bi_growth] = 1./(y[pba->index_bi_a] * pvecback[pba->index_bg_H] * pvecback[pba->index_bg_H]);



  
  /** calculate growth' = 1/(aH^2) */
if (y[pba->index_bi_a]<-1) {
  dy[pba->index_bi_growth_wCDM] = 0.0;
  dy[pba->index_bi_growth_prime_wCDM] = 0.0;

 }
 else  {
   double OmegaDE;
   double wDE;
   double OmegaK;
  
    if (pba->Omega0_lambda != 0.){
      OmegaDE = pvecback[pba->index_bg_rho_lambda]/pvecback[pba->index_bg_H]/pvecback[pba->index_bg_H];
      wDE = -1.;
      OmegaK = 0.;
       }
    else{
      OmegaDE =pvecback[pba->index_bg_Omega_ds];
      wDE = pba->w_ds;
      OmegaK=0.;
      
       }

    //dy[pba->index_bi_growth_wCDM] = 0.0;
    // dy[pba->index_bi_growth_prime_wCDM] = 0.0;
    
    dy[pba->index_bi_growth_wCDM] = (y[pba->index_bi_a] *pvecback[pba->index_bg_H])*y[pba->index_bi_growth_prime_wCDM];
    
    dy[pba->index_bi_growth_prime_wCDM] =-(y[pba->index_bi_a] *pvecback[pba->index_bg_H])*(((5./2.)+(1./2.)*(OmegaK-3.*wDE*OmegaDE))*y[pba->index_bi_growth_prime_wCDM]+(2.*OmegaK+(3./2.)*(1.-wDE)*OmegaDE)*y[pba->index_bi_growth_wCDM]);

   }


  

  if (pba->has_dcdm == _TRUE_){
    /** - compute dcdm density \f$ \rho' = -3aH \rho - a \Gamma \rho \f$*/
    dy[pba->index_bi_rho_dcdm] = -3.*y[pba->index_bi_a]*pvecback[pba->index_bg_H]*y[pba->index_bi_rho_dcdm]-
      y[pba->index_bi_a]*pba->Gamma_dcdm*y[pba->index_bi_rho_dcdm];
  }

  if ((pba->has_dcdm == _TRUE_) && (pba->has_dr == _TRUE_)){
    /** - Compute dr density \f$ \rho' = -4aH \rho - a \Gamma \rho \f$ */
    dy[pba->index_bi_rho_dr] = -4.*y[pba->index_bi_a]*pvecback[pba->index_bg_H]*y[pba->index_bi_rho_dr]+
      y[pba->index_bi_a]*pba->Gamma_dcdm*y[pba->index_bi_rho_dcdm];
  }

  if (pba->has_scf == _TRUE_){
    /** - Scalar field equation: \f$ \phi'' + 2 a H \phi' + a^2 dV = 0 \f$  (note H is wrt cosmic time) */
    dy[pba->index_bi_phi_scf] = y[pba->index_bi_phi_prime_scf];
    dy[pba->index_bi_phi_prime_scf] = - y[pba->index_bi_a]*
      (2*pvecback[pba->index_bg_H]*y[pba->index_bi_phi_prime_scf]
       + y[pba->index_bi_a]*dV_scf(pba,y[pba->index_bi_phi_scf])) ;
  }

  if (pba->has_ds == _TRUE_){
dy[pba->index_bi_m2_hi] =
  y[pba->index_bi_a]
  *pvecback[pba->index_bg_H]
  *pvecback[pba->index_bg_alpha_m_hi];

 
   if (pba->mg_type == fR){

if(y[pba->index_bi_a]>pba->a_start_ds) {
        
    dy[pba->index_bi_f_ds] =
      y[pba->index_bi_a]
      *pvecback[pba->index_bg_H				]
      *y[pba->index_bi_f_prime_ds];

    dy[pba->index_bi_f_prime_ds] =
      y[pba->index_bi_a]*pvecback[pba->index_bg_H]*(
     -(3.*pvecback[pba->index_bg_epsilon_H]
       -1.-pvecback[pba->index_bg_epsilon_H_bar_prime]
       /pvecback[pba->index_bg_epsilon_H_bar])
     *y[pba->index_bi_f_prime_ds]
     +pvecback[pba->index_bg_epsilon_H_bar]*y[pba->index_bi_f_ds]
     +6.*pvecback[pba->index_bg_H]
     *pvecback[pba->index_bg_H]*pvecback[pba->index_bg_epsilon_H_bar]
     *pvecback[pba->index_bg_Omega_ds]
                                                    );
      }


      else {dy[pba->index_bi_f_ds] = 0.;
        dy[pba->index_bi_f_prime_ds] = 0.; }

   }
            }

    
 
  return _SUCCESS_;

}

/**
 * Scalar field potential and its derivatives with respect to the field _scf
 * For Albrecht & Skordis model: 9908085
 * - \f$ V = V_{p_{scf}}*V_{e_{scf}} \f$
 * - \f$ V_e =  \exp(-\lambda \phi) \f$ (exponential) 
 * - \f$ V_p = (\phi - B)^\alpha + A \f$ (polynomial bump) 
 * 
 * TODO: 
 * - Add some functionality to include different models/potentials (tuning would be difficult, though)
 * - Generalize to Kessence/Horndeski/PPF and/or couplings
 * - A default module to numerically compute the derivatives when no analytic functions are given should be added.
 * - Numerical derivatives may further serve as a consistency check.
 * 
 */

/** 
 * 
 * The units of phi, tau in the derivatives and the potential V are the following:
 * - phi is given in units of the reduced Planck mass \f$ m_{pl} = (8 \pi G)^{(-1/2)}\f$
 * - tau in the derivative is given in units of Mpc.
 * - the potential \f$ V(\phi) \f$ is given in units of \f$ m_{pl}^2/Mpc^2 \f$.
 * With this convention, we have
 * \f$ \rho^{class} = (8 \pi G)/3 \rho^{physical} = 1/(3 m_{pl}^2) \rho^{physical} = 1/3 * [ 1/(2a^2) (\phi')^2 + V(\phi) ] \f$
    and \f$ \rho^{class} \f$ has the proper dimension \f$ Mpc^-2 \f$.
 */

double V_e_scf(struct background *pba,
               double phi
               ) {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return  exp(-scf_lambda*phi);
}

double dV_e_scf(struct background *pba,
                double phi
                ) {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return -scf_lambda*V_scf(pba,phi);
}

double ddV_e_scf(struct background *pba,
                 double phi
                 ) {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return pow(-scf_lambda,2)*V_scf(pba,phi);
}


/** parameters and functions for the polynomial coefficient
 * \f$ V_p = (\phi - B)^\alpha + A \f$(polynomial bump) 
 * 
 * double scf_alpha = 2;
 * 
 * double scf_B = 34.8;
 * 
 * double scf_A = 0.01; (values for their Figure 2)
 */

double V_p_scf(
               struct background *pba,
               double phi) {
  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return  pow(phi - scf_B,  scf_alpha) +  scf_A;
}

double dV_p_scf(
                struct background *pba,
                double phi) {

  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return   scf_alpha*pow(phi -  scf_B,  scf_alpha - 1);
}

double ddV_p_scf(
                 struct background *pba,
                 double phi) {
  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return  scf_alpha*(scf_alpha - 1.)*pow(phi -  scf_B,  scf_alpha - 2);
}

/** Fianlly we can obtain the overall potential \f$ V = V_p*V_e \f$
 */

double V_scf(
             struct background *pba,
             double phi) {
  return  V_e_scf(pba,phi)*V_p_scf(pba,phi);
}

double dV_scf(
              struct background *pba,
	      double phi) {
  return dV_e_scf(pba,phi)*V_p_scf(pba,phi) + V_e_scf(pba,phi)*dV_p_scf(pba,phi);
}

double ddV_scf(
               struct background *pba,
               double phi) {
  return ddV_e_scf(pba,phi)*V_p_scf(pba,phi) + 2*dV_e_scf(pba,phi)*dV_p_scf(pba,phi) + V_e_scf(pba,phi)*ddV_p_scf(pba,phi);
}






  int background_loop_over_b_ds(
		       struct precision *ppr,
		       struct background *pba
		       ) {


 
    double b_ds,b_min_ds,b_max_ds;
    b_min_ds = pba->b_min_ds;
    b_max_ds = pba->b_max_ds;
    int i_b_ds;
    int N_b_ds = pba->N_b_ds;
    double error_B0 = 0.0001;
    int verbose_temp = pba->background_verbose;
    pba->background_verbose = 0;
    b_ds = pba->b_ds;


  
  //////////////LOOP OVER VALUES OF b_ds
  /*
FILE * fp;
  
fp=fopen("output/B0.txt", "w");
if(fp == NULL)
    exit(-1); 
fprintf(fp, "#1:b\t\t\t 2:B0\n");
  

for( i_b_ds=0; i_b_ds< N_b_ds; i_b_ds++){
      
      pba->b_ds = b_min_ds + (b_max_ds-b_min_ds)*i_b_ds/(N_b_ds-1.);
      
      class_call(compute_B0_ds(ppr,pba),
             pba->error_message,
		 pba->error_message);
      printf("B0_ds = %e b_ds =  %e w_ds = %e\n ",
	     pba->B0_ds,
	     pba->b_ds,
	     pba->w_ds);

       fprintf(fp,"%e\t\t %e \n",
		    pba->b_ds,
                    pba->B0_ds);
     
     }
		    fclose(fp);

  */

 double B0_right,B0_left;
 double delta_b = 0.1;
      //step 1
      pba->b_ds = pba->b_max_ds;
      
      class_call(compute_B0_ds(ppr,pba),
             pba->error_message,
		 pba->error_message);

      B0_right = pba->B0_ds;

      //printf("B0_asked_ds = %e, B0_ds = %e b_ds =  %e w_ds =  %e\n",pba->B0_asked_ds,pba->B0_ds,pba->b_ds,pba->w_ds);

      
      //step 2
      pba->b_ds = pba->b_max_ds - delta_b;
      
      class_call(compute_B0_ds(ppr,pba),
             pba->error_message,
		 pba->error_message);

      B0_left = pba->B0_ds;

      if(B0_right>pba->B0_asked_ds) {
	if(B0_left<B0_right){
	if(B0_left<pba->B0_asked_ds)
	  pba->b_min_ds = pba->b_ds;

        else {
	while (pba->B0_ds>pba->B0_asked_ds) {
          pba->b_ds = pba->b_ds - delta_b;
          class_call(compute_B0_ds(ppr,pba),
                 pba->error_message,
		 pba->error_message);
        }
	  pba->b_min_ds = pba->b_ds;
	}
      }
	else {
	  pba->b_ds = pba->b_min_ds;
	while (pba->B0_ds>pba->B0_asked_ds) {
          pba->b_ds = pba->b_ds + delta_b;
          class_call(compute_B0_ds(ppr,pba),
                 pba->error_message,
		 pba->error_message);
        }
	  pba->b_max_ds = pba->b_ds;
	     }
      }

      else  {  //B0_right<pba->B0_asked_ds
	if(B0_left>B0_right){
	if(B0_left>pba->B0_asked_ds)
	  pba->b_min_ds = pba->b_ds;

        else {
	while (pba->B0_ds<pba->B0_asked_ds) {
          pba->b_ds = pba->b_ds - delta_b;
          class_call(compute_B0_ds(ppr,pba),
                 pba->error_message,
		 pba->error_message);
        }
	  pba->b_min_ds = pba->b_ds;
	}
      }
	else {
	  pba->b_ds = pba->b_min_ds;
	while (pba->B0_ds<pba->B0_asked_ds) {
          pba->b_ds = pba->b_ds + delta_b;
          class_call(compute_B0_ds(ppr,pba),
                 pba->error_message,
		 pba->error_message);
        }
	  pba->b_max_ds = pba->b_ds;
	     }
      }


      //printf("b_min_ds =  %e b_max_ds = %e\n ",pba->b_min_ds,pba->b_max_ds);

      pba->b_ds = pba->b_min_ds;
          class_call(compute_B0_ds(ppr,pba),
                 pba->error_message,
		 pba->error_message);
	  B0_left = pba->B0_ds;
	  
      pba->b_ds = pba->b_max_ds;
          class_call(compute_B0_ds(ppr,pba),
                 pba->error_message,
		 pba->error_message);
	  B0_right = pba->B0_ds;


     pba->B0_ds = 1.e9;  
     while (fabs(pba->B0_ds/pba->B0_asked_ds-1.)>error_B0){

      pba->b_ds = (pba->b_min_ds+pba->b_max_ds)/2.;
          class_call(compute_B0_ds(ppr,pba),
                 pba->error_message,
		 pba->error_message);
	  if(B0_left<B0_right){
	  if (pba->B0_ds>pba->B0_asked_ds)
	    pba->b_max_ds = pba->b_ds;
          else pba->b_min_ds = pba->b_ds; 
	  }
	  else{
	  if (pba->B0_ds>pba->B0_asked_ds)
	    pba->b_min_ds = pba->b_ds;
          else pba->b_max_ds = pba->b_ds; 
	  }

	  /* printf("B0_ds = %e b_ds =  %e w_ds = %e er = %e\n ",
	     pba->B0_ds,
	     pba->b_ds,
		       pba->w_ds,
		       fabs(pba->B0_ds/pba->B0_asked_ds-1.));*/
     }

 

  pba->background_verbose = verbose_temp;
  return _SUCCESS_;

    }



  int compute_B0_ds(
		       struct precision *ppr,
		       struct background *pba
				) 

 {
 double lambda = 3.*(1.+pba->w_ds);

 double Omega_ds_ini = (pba->Omega0_ds * pow(pba->H0,2)/pow(pba->a_start_ds,3.*(1.+pba->w_ds)))
                          /( (pba->Omega0_cdm * pow(pba->H0,2) / pow(pba->a_start_ds,3))                                                         +(pba->Omega0_b * pow(pba->H0,2) / pow(pba->a_start_ds,3)));

 double Omega_ds_prime_ini= (3.-3.*(1.+pba->w_ds))*Omega_ds_ini; 

 double C_ds = 18.*pba->Omega0_ds*pow(pba->H0,2)
                 /(2.*lambda*lambda
                   -7.*lambda
                   -3);
 double n_plus = (-7.+sqrt(7.*7.-4.*2.*(-3)))/(2.*2.);



    
    double f_part_ini = C_ds*pow(pba->a_start_ds,-lambda);
    double f_part_prime_ini = -lambda*f_part_ini;

    

  double tau_start;
  double tau_end;
  int i;
   
    struct background_parameters_and_workspace bpaw;

  double * pvecback_integration;
  double * pvecback;
  int last_index=0;

  bpaw.pba = pba;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  bpaw.pvecback = pvecback;


  class_alloc(pvecback_integration,pba->bi_size*sizeof(double),pba->error_message);

     growTable gTable;
     struct generic_integrator_workspace gi;
  
     double b_ds = pba->b_ds;
    
    pba->f_ini_ds = b_ds*C_ds*pow(pba->a_start_ds,n_plus)+f_part_ini;
    pba->f_prime_ini_ds = n_plus*b_ds*C_ds*pow(pba->a_start_ds,n_plus)+f_part_prime_ini;
    

  class_call(initialize_generic_integrator((pba->bi_size-1),&gi),
             gi.error_message,
             pba->error_message);

  class_call(background_initial_conditions(ppr,pba,pvecback,pvecback_integration),
             pba->error_message,
             pba->error_message);

  tau_end=pvecback_integration[pba->index_bi_tau];

  class_call(gt_init(&gTable),
             gTable.error_message,
             pba->error_message);
  pba->bt_size=0;

  while (pvecback_integration[pba->index_bi_a] < pba->a_today) {

    tau_start = tau_end;

    class_call(background_functions(pba,pvecback_integration, pba->short_info, pvecback),
               pba->error_message,
               pba->error_message);

    if ((pvecback_integration[pba->index_bi_a]*(1.+ppr->back_integration_stepsize)) < pba->a_today) {
      tau_end = tau_start + ppr->back_integration_stepsize / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);
    }
    else {
      tau_end = tau_start + (pba->a_today/pvecback_integration[pba->index_bi_a]-1.) / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);
    }

    class_test((tau_end-tau_start)/tau_start < ppr->smallest_allowed_variation,
               pba->error_message,
               "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",(tau_end-tau_start)/tau_start);

    class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
               gTable.error_message,
               pba->error_message);
    pba->bt_size++;

    class_call(generic_integrator(background_derivs,
                                  tau_start,
                                  tau_end,
                                  pvecback_integration,
                                  &bpaw,
                                  ppr->tol_background_integration,
                                  ppr->smallest_allowed_variation,
                                  &gi),
               gi.error_message,
               pba->error_message);

    pvecback_integration[pba->index_bi_tau]=tau_end;
  }

  /** - free the growTable with gt_free() */

    class_call(gt_free(&gTable),
             gTable.error_message,
             pba->error_message);

  /** - clean up generic integrator with cleanup_generic_integrator() */
    class_call(cleanup_generic_integrator(&gi),
	       gi.error_message,
	       pba->error_message);


  free(pvecback);
  free(pvecback_integration);


  b_ds = pba->b_ds;
  pba->f_ini_ds = b_ds*C_ds*pow(pba->a_start_ds,n_plus)+f_part_ini;
  pba->f_prime_ini_ds = n_plus*b_ds*C_ds*pow(pba->a_start_ds,n_plus)+f_part_prime_ini;
    
  
    
  return _SUCCESS_;

    }
