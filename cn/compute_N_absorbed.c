/*--------------------------------------------------------------*/
/*                                                              */
/*		compute_N_absorbed				*/
/*                                                              */
/*  NAME                                                        */
/*		compute_N_absorbed				*/
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/*  void compute_N_absorbed(int					*/
/*					double	,		*/
/*					double	,		*/
/*					double	,		*/
/*					double	,		*/
/*					double);		*/
/*                                                              */
/*  OPTIONS                                                     */
/*                                                              */
/*                                                              */
/*  DESCRIPTION                                                 */
/*                                                              */
/*								*/
/*  PROGRAMMER NOTES                                            */
/* 	following the model (and parameters) of 		*/
/*	Kothawala and Moore, Canadian J. For Research, 2009	*/
/*	39:2381-3290						*/
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"
#define  PARTICLE_DENSITY    2.65    /* soil particle density g/cm3 (Dingman) */

double	compute_N_absorbed(int verbose_flag, 
			double z1, //<<--- sat_def_z
			double z2, //<<--- soil depth
			double N_absorption_rate, 
			double p, //porosity_decay in depth (m)
			double n_0) //porosity_0
			
	{
        //z2 must be greater than z1 numerically
	/*------------------------------------------------------*/ 
	/*	Local Function Declarations.						*/ 
	/*------------------------------------------------------*/
   		
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	double nabsorbed;       /*kg*/
	double bulk_density;
    double kg_soil;
    
	//bulk_density = PARTICLE_DENSITY * (1.0 - n_0) * 1000;
	//nabsorbed=n_0*(z2-z1)*N_absorption_rate*bulk_density; // original
    //nabsorbed=n_0* max(z2-z1,0.0) *N_absorption_rate*bulk_density; // upper 6a
        
        //nabsorbed = N_absorption_rate  *  n_0 (rescaling the surface of particule) * (bulk_density*max(z2-z1,0.0)); // upper 6a
        //nabsorbed = N_absorption_rate *  n_0 (where is this come from?) *  (kg soil);
        //nabsorbed = N_absorption_rate *  (kg soil);

        //bulk_density = PARTICLE_DENSITY * (1.0 - porosity) * 1000; <<--- 1000 is converting unit o kg/m3
        //bulk_density = PARTICLE_DENSITY * 1000 * ( z2-z1 -n_0*(exp(-p*z2)-exp(-p*z1))/p ) / (z2-z1);
        // n_0 = patch[0].soil_defaults[0][0].porosity_0
        // p = patch[0].soil_defaults[0][0].porosity_decay
        
        
    //kg_soil = PARTICLE_DENSITY * n_0 * (z2-z1 + n_0*(exp(-z2/p)-exp(-z1/p))*p ) * 1000.0; // correct
    //nabsorbed=  N_absorption_rate * kg_soil ; //integrated over z2-z1  upper 7 testing
    //nabsorbed=max(nabsorbed, 0.0);
        
        // N_absorption_rate * PARTICLE_DENSITY * [porosity * (1-porosity)] | z1,z2
        nabsorbed = (p*n_0*(exp(-z2/p)-exp(-z1/p)) + 0.5*p*n_0*n_0*(exp(-2*z1/p)-exp(-2*z2/p)))*PARTICLE_DENSITY*N_absorption_rate;
        
        
	
	return(nabsorbed);
} /* end compute_N_absorbed */

