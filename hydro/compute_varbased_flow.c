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
				double interval_size,
				double *transmissivity,
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
	double p;
	double n_0;
	double soil_depth;
	double fs_spill;
	double fs_percolation;


	/*--------------------------------------------------------------*/
	/* calculate or initialize value    				*/
	/*--------------------------------------------------------------*/
	p = patch[0].soil_defaults[0][0].porosity_decay;
	n_0 = patch[0].soil_defaults[0][0].porosity_0;
	soil_depth = patch[0].soil_defaults[0][0].soil_depth;
    
	threshold = n_0 * p * (1 - exp(-soil_depth/p)) * (1 - patch[0].soil_defaults[0][0].fs_threshold);// fs_threshold default is 0.2
        // threshold = 80% of full column of water
        // if porosity_decay is negative, then (1 - exp(-soil_depth/p)) is negative --> threshold <0 
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

	
	flow = 0.0;
	if (s1 < 0.0) s1 = 0.0;	// s1 = sat_def

	if (std > ZERO) {
        
        //std is zero for normal runs  //(trigger by the -stdev flag)
        
        for (i=0; i <9; i++) {
            didx = (int) lround((s1 + normal[i]*std)/interval_size);
            if (didx > num_soil_intervals) didx = num_soil_intervals;

            accum = transmissivity[didx] * 1;
            /* fill and spill */
            if ((patch[0].sat_deficit <= threshold) && ((s1 + normal[i]*std) <= threshold)){
                accum=transmissivity[didx] * 1;
            }

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
        
        
		didx = (int) lround(s1/interval_size); //sat_def
		didthr = (int) lround(threshold/interval_size);///<<---------- this is new!!
		if (didx > num_soil_intervals) didx = num_soil_intervals;

		/* default lateral flow (seepage) below the threshold is 1/3 of the original value, the multiplier is arbitary. Xiaoli */
        // why not calculate transimissivity at real time? performance issue?
        // to have a pre-calculated table to look up because gamma is transimissivity at surface and it decay with depth

		/* if sat_deficit > threshold */
		if(patch[0].sat_deficit > threshold){
            // water table is very low, less than 80% of the full capacity
		    flow = transmissivity[didx] * fs_percolation; // fs_percolation defaults = 1 
            //fs_percolation = 1 mean no loss flow to GW?
		}else{
            // if water level exceed moisture threshold (or sat_deficit <= soil deficit threshold)
		    thre_flow=transmissivity[didthr];
		    abovthre_flow = (transmissivity[didx]-thre_flow) * fs_spill; // fs_spill default value is 1 
            flow = abovthre_flow + thre_flow * fs_percolation;  // fs_percolation defaults = 1
		}

	
	}//std

	flow = flow*gamma;

	  return(flow);
} /*compute_varbased_flow*/
