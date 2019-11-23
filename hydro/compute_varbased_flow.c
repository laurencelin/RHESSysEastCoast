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
				int num_soil_intervals,
				double std,
				double s1,
				double gamma,	
				double no_used1,
				double no_used2,
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
	double threshold;
	//double p_decay;
	//double p_0;
	//double soil_depth;
	double fs_spill;
	double fs_percolation;


	/*--------------------------------------------------------------*/
	/* calculate or initialize value    				*/
	/*--------------------------------------------------------------*/
//	p_decay = patch[0].soil_defaults[0][0].porosity_decay;
//	p_0 = patch[0].soil_defaults[0][0].porosity_0;
//	soil_depth = patch[0].soil_defaults[0][0].soil_depth;
    
    // p_0 * p_decay * (1 - exp(-soil_depth/p_decay))
	threshold = patch[0].soil_defaults[0][0].soil_water_cap * (1 - patch[0].soil_defaults[0][0].fs_threshold);// fs_threshold default is 0.2
        // threshold = 80% of full column of water
	fs_spill = patch[0].soil_defaults[0][0].fs_spill;
	fs_percolation = patch[0].soil_defaults[0][0].fs_percolation;

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
	if (s1 < 0.0) s1 = 0.0;	// s1 = sat_def

	if (std > ZERO) {
        
        //std is zero for normal runs  //(trigger by the -stdev flag)
        for (i=0; i <9; i++) {
            // s1 = sat_def
            // interval_size =
            //didx = (int) lround( (s1 + normal[i]*std)/interval_size );
            //if (didx > num_soil_intervals) didx = num_soil_intervals;
            //accum = transmissivity[didx];
            
            didx = patch[0].sat_def_pct_index + (int)(normal[i]*std*1000);
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
        
//        double    compute_varbased_flow(
//                                        int num_soil_intervals,
//                                        double std,
//                                        double s1,  <<---- sat-deficit
//                                        double gamma, <<---- total gamma
//                                        double interval_size,
//                                        double *transmissivity,
//                                        struct patch_object *patch)
        
        
//		didx = (int) lround(s1/interval_size); //sat_def
//		didthr = (int) lround(threshold/interval_size);///<<---------- this is new!!
//		if (didx > num_soil_intervals) didx = num_soil_intervals;

		/* default lateral flow (seepage) below the threshold is 1/3 of the original value, the multiplier is arbitary. Xiaoli */
        // why not calculate transimissivity at real time? performance issue?
        // to have a pre-calculated table to look up because gamma is transimissivity at surface and it decay with depth

		/* if sat_deficit > threshold */
		if(patch[0].sat_deficit > threshold){
            // water table is very low, less than 80% of the full capacity
		    //flow = transmissivity[didx] * fs_percolation; // fs_percolation defaults = 1
            //fs_percolation = 1 means no loss flow to GW?
            
            //didx=patch[0].sat_def_pct_index
            //if(didx<0) didx = 0;
            
            flow = min( gamma*(patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].transmissivity_dailyflux[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].transmissivity_dailyflux[patch[0].sat_def_pct_index]), patch[0].area*(patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].transmissivity_maxdailyflux[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].transmissivity_maxdailyflux[patch[0].sat_def_pct_index]));
            flow *= fs_percolation;
    
		}else{
            // if water level exceed moisture threshold (or sat_deficit <= soil deficit threshold)
            // threshold = patch[0].soil_defaults[0][0].sat_zZ[0] * (1 - patch[0].soil_defaults[0][0].fs_threshold);
            
            //  didx = (int) lround(s1/interval_size); //sat_def
            //  didthr = (int) lround(threshold/interval_size);///<<---------- this is new!!
            
            didthr = (int)(threshold*1000*patch[0].soil_defaults[0][0].max_sat_def_1); // volumn
            thre_flow = min(gamma*patch[0].soil_defaults[0][0].transmissivity_dailyflux[didthr], patch[0].area*patch[0].soil_defaults[0][0].transmissivity_maxdailyflux[didthr]);
            
            abovthre_flow = min(gamma*(patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].transmissivity_dailyflux[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].transmissivity_dailyflux[patch[0].sat_def_pct_index]), patch[0].area *(patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].transmissivity_maxdailyflux[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].transmissivity_maxdailyflux[patch[0].sat_def_pct_index]));
            abovthre_flow -= thre_flow;
            abovthre_flow *= (abovthre_flow>0.0? fs_spill : 0.0); // fs_spill default = 1
            
            //thre_flow=transmissivity[didthr];
		    //abovthre_flow = (transmissivity[didx]-thre_flow) * fs_spill; // fs_spill default value is 1
            
            flow = abovthre_flow + thre_flow*fs_percolation;  // fs_percolation defaults = 1
		}

	
	}//std or not

	//flow = flow*gamma; // merge this calculation to above
    return(flow); // it's volume, not mm
} /*compute_varbased_flow*/
