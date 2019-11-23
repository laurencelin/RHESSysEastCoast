/*--------------------------------------------------------------*/
/* 								*/
/*		compute_transmissivity_curve			*/
/*								*/
/*	NAME							*/
/*	compute_transmissivity_curve - estimates transmissivity	*/
/*		assuming an exponential decay of transmivivity	*/
/*		with depth - note returned value is 		*/
/*		relative transmissivity (i.e per unit		*/
/*		Ksat at the surface)				*/
/*								*/
/*								*/
/*	SYNOPSIS						*/
/*	compute_transmissivity_curve(				*/
/*				double	,			*/
/*				double	,			*/
/*				double	)			*/
/*								*/
/*	returns:						*/
/*	transmissivity - (unitless) multiplier for Ksat0 	*/
/*		to calculate transmissivity over range of	*/
/*		depths specified				*/
/*								*/
/*	OPTIONS							*/
/*	double	m - Ksat decay parameter			*/
/*	double	z - (m) depth to the water table		*/
/*	double  D - maximum depth				*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	computes transmissivity multiplier over range of 	*/
/*	depths between z and D; assumes an exponential decay	*/
/*	of Ksat with depth (decay given by m) until z_layer1	*/
/*	and then a constant conductivity below that layer	*/
/*	given by Ksat at the surface?				*/
/*	Note that if m is 0, we assume that Ksat is constant    */
/*	with depth						*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

double 	*compute_transmissivity_curve( 
					double  gamma,
					struct patch_object *patch,
					struct command_line_object *command_line
					)
{

	
	/*--------------------------------------------------------------*/
	/*	Local function definition.				*/
	/*--------------------------------------------------------------*/
	double compute_field_capacity(
		int,
		int,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double);


	double	compute_delta_water(
		int,
		double,
		double,
		double,
		double,
		double);


	double	compute_z_final(
		int,
		double,
		double,
		double,
		double,
		double);

	 void    *alloc( size_t, char *, char *);


	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/

	int didx,initial;
    double	lower, depth;
    //double horizontal_ksat_decay;
	double	lower_z, depth_z;
	double start, fclayer, potential_sat;
	double transmissivity_layer;
	double *transmissivity;

    //horizontal_ksat_decay = patch[0].soil_defaults[0][0].m; // * patch[0].horizontal_k_SCALE;

	transmissivity = (double *) alloc((patch[0].num_soil_intervals+1) * sizeof(double),
					"trans","compute_transmissivity_cuve");


	/*--------------------------------------------------------------*/
	/*	for do not include surface overland flow or detention   */
	/*	storage here						*/
	/*--------------------------------------------------------------*/
//                          transmissivity      depth(z) at lower boundary
//    |0                    [1]+something
//    |1                    [2]+something
//    |2                    [3]+something
//                                              soildepth - z(interval_size)
//    |num_soil_intervals   0                   soildepth

    
    
	if (patch[0].soil_defaults[0][0].soil_water_cap > patch[0].soil_defaults[0][0].interval_size){
        //whole column water space usually larger than layor space!
        // soil_water_cap is volumn of water
        
        // looks like looping from bottom to up
		initial = patch[0].num_soil_intervals; // last index
		depth = patch[0].soil_defaults[0][0].soil_water_cap; // full soil vol
		transmissivity[initial]=0.0; // last index; at the bottom; from D to D --> 0 trans
		initial = initial-1;
		for (didx=initial; didx >= 0; didx -= 1) {
            lower = depth;//initially whole column water;
			depth = depth-patch[0].soil_defaults[0][0].interval_size; // set up for next step

			lower_z = compute_z_final(
				command_line[0].verbose_flag,
				patch[0].soil_defaults[0][0].porosity_0,
				patch[0].soil_defaults[0][0].porosity_decay,
				patch[0].soil_defaults[0][0].soil_depth,
				0.0,// initial depth
				-1.0*lower); // the sat_def_z @ vol=lower (bottom)
            
			depth_z = compute_z_final(
				command_line[0].verbose_flag,
				patch[0].soil_defaults[0][0].porosity_0,
				patch[0].soil_defaults[0][0].porosity_decay,
				patch[0].soil_defaults[0][0].soil_depth,
				0.0,
				-1.0*depth); // the sat_def_z @ vol=lower-intervals (top)
            
            
            // field capacity: from water table to surface
            // lower_z and depth_z are the limits on the integration about theta-psi
            fclayer = compute_field_capacity(
				command_line[0].verbose_flag,
				patch[0].soil_defaults[0][0].theta_psi_curve,
				patch[0].soil_defaults[0][0].psi_air_entry,
				patch[0].soil_defaults[0][0].pore_size_index,
				patch[0].soil_defaults[0][0].p3,
				patch[0].soil_defaults[0][0].p4,
				patch[0].soil_defaults[0][0].porosity_0,
				patch[0].soil_defaults[0][0].porosity_decay,
				patch[0].soil_defaults[0][0].soil_depth,
				lower_z, // bottom
				depth_z);// top

            // this block here is doing some field capacity correction, involving gamma.
            // m = patch[0].soil_defaults[0][0].m;
            if (patch[0].soil_defaults[0][0].m > ZERO) transmissivity_layer = gamma *
                ( exp( -1.0 * max(depth, 0.0)/patch[0].soil_defaults[0][0].m ) -
                  exp ( -1.0 * lower/patch[0].soil_defaults[0][0].m) );   // multipled by gamma
            else transmissivity_layer =  gamma * (lower-depth); // when no ksat decay
                // this function is called during initial!

            
            // as if the layer_water is drained horizontially and a fc is left behind.
            fclayer = max( patch[0].soil_defaults[0][0].interval_size - fclayer, 0.0); // volumn space after accounted for field capacity
            transmissivity_layer = min(fclayer, transmissivity_layer/patch[0].area);

            if (gamma > ZERO)
                transmissivity[didx] = transmissivity[didx+1] + transmissivity_layer * patch[0].area / gamma; // divided by gamma
            else
                transmissivity[didx] = transmissivity[didx+1];

		}//for didx
	
	
    }else {
        // whole column water space is less than a layor!!
		initial = 1;
		transmissivity[initial]=0.0;
		lower = patch[0].soil_defaults[0][0].soil_water_cap;
		depth = 0;
		if (patch[0].soil_defaults[0][0].m > ZERO)
			transmissivity[initial-1] =   (exp ( -1.0 * (max(depth, 0.0)*patch[0].soil_defaults[0][0].m)) - exp ( -1.0 * (lower*patch[0].soil_defaults[0][0].m)));
		else
			transmissivity[initial-1] =  (lower-depth);
		
    }//else

	return(transmissivity);

} /*compute_transmissivity_curve*/
