/*--------------------------------------------------------------*/
/* 								*/
/*		compute_varbased_flow			*/
/*								*/
/*	NAME							*/
/*	compute_varbased_flow - estimates subsurface flow	*/
/*	by assuming variance around mean saturation deficit	*/
/*								*/
/*								*/
/*	SYNOPSIS						*/
/*	compute_varbased_flow(				*/
/*				double	,			*/
/*				double	,			*/
/*				double	,			*/
/*				double	,			*/
/*				double	,			*/
/*				struct patch_object *patch)	    	*/
/*								*/
/*	returns:						*/
/*	estimate of subsurface flow from patch			*/
/*								*/
/*	OPTIONS							*/
/*	double	std - standard deviation of normal distrib	*/
/*	double gamma						*/
/*	double	m - Ksat decay parameter			*/
/*	double	z - (m) depth to the water table		*/
/*	double  D - maximum depth				*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	computes subsurface throughflow by computing 		*/
/*	transmissivity based on a normal distribution 		*/
/*	of saturation deficits (around a given mean)		*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

double	compute_varbased_flow(
				int num_soil_intervals, // <-- patch[0].num_soil_intervals
				double std, // <-- patch[0].std * std_scale; std_scale = 0 by default
				int sat_def_pct_index, // <-- patch[0].sat_def_pct_index
                double sat_def_pct_indexM, // <-- patch[0].sat_def_pct_indexM
				double gamma,
				struct patch_object *patch)
{


	/*--------------------------------------------------------------*/
	/*	Local sub	definition				*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/

	double	flow, accum,thre_flow,abovthre_flow;
	double	normal[9], perc[9];
	int i;
	int didx,didthr;
	

	// soil deficit threshold, not the soil moisture threshold, fs_threshold is defined as the soil moiture threshold as
	// percentage of the max soil moisture holding capacity
//	double threshold;
//	double fs_spill;
//	double fs_percolation;


	/*--------------------------------------------------------------*/
	/* calculate or initialize value    				*/
	/*--------------------------------------------------------------*/
    
    // p_0 * p_decay * (1 - exp(-soil_depth/p_decay))
//	threshold = patch[0].soil_defaults[0][0].soil_water_cap * (1 - patch[0].soil_defaults[0][0].fs_threshold);// fs_threshold default is 0.2
//        // threshold = 80% of full column of water
//	fs_spill = patch[0].soil_defaults[0][0].fs_spill;
//	fs_percolation = patch[0].soil_defaults[0][0].fs_percolation;

	normal[0] = 0;
	normal[1] = 0.253;
	normal[2] = 0.524;
	normal[3] = 0.842;
	normal[4] = 1.283;
	normal[5] = -0.253;
	normal[6] = -0.524;
	normal[7] = -0.842;
	normal[8] = -1.283;

	perc[0] = 0.2;
	for (i=1; i<9; i++) perc[i] = 0.1;
    // c(0.2, 0.1, 0.1, 0.1, ..... 0.1); weight of the normal distribution
	
	flow = 0.0;
	//if (s1 < 0.0) s1 = 0.0;	// s1 = sat_def

	if (std > ZERO) {
        
        //std is zero for normal runs  //(trigger by the -stdev flag)
        for (i=0; i <9; i++) {
            // s1 = sat_def
            // interval_size =
            //didx = (int) lround( (s1 + normal[i]*std)/interval_size );
            //if (didx > num_soil_intervals) didx = num_soil_intervals;
            //accum = transmissivity[didx];
            
            didx = sat_def_pct_index + (int)(normal[i]*std*1000);
            if(didx<0) didx = 0;
            if(didx>1000) didx = 1000;
            accum = min( gamma*patch[0].soil_defaults[0][0].transmissivity_dailyflux[didx], patch[0].area * patch[0].soil_defaults[0][0].transmissivity_maxdailyflux[didx]);
            
            /* fill and spill */ // not sure what's this do! basically nothings
//            if ( (patch[0].sat_deficit <= threshold) && (s1+normal[i]*std <= threshold) ){
//                accum=transmissivity[didx];
//            }

            flow += accum * perc[i];
        }// end of for loop
        
	}else{
       //normal runs

        // note: "gamma" has patch[0].area
        flow = min(
                   gamma*(sat_def_pct_indexM * patch[0].soil_defaults[0][0].transmissivity_dailyflux[sat_def_pct_index+1] + (1.0-sat_def_pct_indexM) * patch[0].soil_defaults[0][0].transmissivity_dailyflux[sat_def_pct_index]),
                   patch[0].area*(sat_def_pct_indexM * patch[0].soil_defaults[0][0].transmissivity_maxdailyflux[sat_def_pct_index+1] + (1.0-sat_def_pct_indexM) * patch[0].soil_defaults[0][0].transmissivity_maxdailyflux[sat_def_pct_index]));
        

	}//std or not

	//flow = flow*gamma; // merge this calculation to above
    return(flow); // it's volume, not mm
} /*compute_varbased_flow*/
