/*--------------------------------------------------------------*/
/* 								*/
/*		update_gw_drainage					*/
/*								*/
/*								*/
/*	NAME							*/
/*	update_gw_drainage -  					*/
/* 		drainage shallow subsurface saturation zone	*/
/*		from each patch to a regional (hillslope scale)	*/
/*		groundwater store				*/
/*		nitrogen is also drained using assumption 	*/
/*		of an exponential decay of N with depth		*/
/*								*/
/*	SYNOPSIS						*/
/*	int update_gw_drainage(					*/
/*			struct patch_object *			*/
/*			struct hillslope_object *		*/
/*			struct command_line_object *		*/
/*			struct date,				*/
/*			)					*/
/*								*/
/*	returns:						*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*	preset code just uses a user assigned loading rate	*/
/*	and all of it is nitrate				*/
/*								*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

int update_gw_drainage(
				  struct  patch_object   *patch,
				  struct  hillslope_object *hillslope,
				  struct  zone_object *zone,
				  struct  command_line_object *command_line,
				  struct	date	current_date)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.			*/
	/*------------------------------------------------------*/
	
	double  compute_z_final(
		int,
		double,
		double,
		double,
		double,
		double); // not in use


	/*------------------------------------------------------*/
	/*	Local Variable Definition. 			*/
	/*------------------------------------------------------*/
	int ok = 1;
	double drainage,sat_store,N_loss;
	double preday_sat_deficit_z, add_field_capacity;
	double sat_to_gw_coeff;
    double active_zone_z;
    double decay_rate;
    double perc_inroot;
    double sat_def, sat_defz;
    double unsat;
	/*------------------------------------------------------*/
	/*		assume percent of incoming precip	*/
	/*------------------------------------------------------*/
	// Note: Mult by Ksat_vertical so that precip onto impervious surfaces does not contribute to GW
	if (zone[0].hourly_rain_flag==1){
	  sat_to_gw_coeff = patch[0].soil_defaults[0][0].sat_to_gw_coeff / 24.0;
	}
	else sat_to_gw_coeff = patch[0].soil_defaults[0][0].sat_to_gw_coeff; 

    drainage = sat_to_gw_coeff * patch[0].detention_store * patch[0].Ksat_vertical; // correcting with % imprevious surface
	patch[0].detention_store -= drainage;
	patch[0].gw_drainage = drainage;
	hillslope[0].gw.storage += (drainage * patch[0].area / hillslope[0].area);

	/*------------------------------------------------------*/
	/*	determine associated N leached			*/
	/*------------------------------------------------------*/
	if (patch[0].surface_DON > ZERO && command_line[0].grow_flag > 0) {
		N_loss = sat_to_gw_coeff * patch[0].surface_DON * command_line[0].Rsolute2gw;
		hillslope[0].gw.DON += (N_loss * patch[0].area / hillslope[0].area);
		patch[0].ndf.DON_to_gw = N_loss;
		patch[0].surface_DON -= N_loss;
		}
	if (patch[0].surface_DOC > ZERO && command_line[0].grow_flag > 0) {
		N_loss = sat_to_gw_coeff * patch[0].surface_DOC * command_line[0].Rsolute2gw;
		hillslope[0].gw.DOC += (N_loss * patch[0].area / hillslope[0].area);
		patch[0].cdf.DOC_to_gw = N_loss;
		patch[0].surface_DOC -= N_loss;
		}
    if (patch[0].surface_NH4 > ZERO && command_line[0].grow_flag > 0) {
        N_loss = sat_to_gw_coeff * patch[0].surface_NH4 * command_line[0].Rsolute2gw;
        hillslope[0].gw.NH4 += (N_loss * patch[0].area / hillslope[0].area);
        patch[0].ndf.N_to_gw += N_loss;
        patch[0].surface_NH4 -= N_loss;
    }
    
    if (patch[0].surface_NO3 > ZERO && command_line[0].grow_flag > 0) {
        N_loss = sat_to_gw_coeff * patch[0].surface_NO3 * command_line[0].Rsolute2gw;
        hillslope[0].gw.NO3 += (N_loss * patch[0].area / hillslope[0].area);
        patch[0].ndf.N_to_gw += N_loss;
        patch[0].surface_NO3 -= N_loss;
    }
    
    ///// ------ newly added (try this for NO3)
    
    
//    active_zone_z = (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z));
//    decay_rate = (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NO3decayRate : patch[0].soil_defaults[0][0].N_decay_rate);
//    perc_inroot = max(1 - (1.0-exp(-decay_rate * patch[0].sat_deficit_z)) / (1.0 - exp(-decay_rate * active_zone_z)),0.0);
    
   // ---- is it the problem?
//    perc_inroot = command_line[0].soluteLoss2GW;
//    if (patch[0].soil_ns.DON > ZERO) {
//        N_loss = sat_to_gw_coeff * patch[0].soil_ns.DON*perc_inroot * command_line[0].Rsolute2gw;
//        hillslope[0].gw.DON += (N_loss * patch[0].area / hillslope[0].area); // does this area correction necessary?
//        patch[0].ndf.DON_to_gw += N_loss;
//        patch[0].soil_ns.DON -= N_loss;
//    }
//    if (patch[0].soil_cs.DOC > ZERO) {
//        N_loss = sat_to_gw_coeff * patch[0].soil_cs.DOC*perc_inroot * command_line[0].Rsolute2gw;
//        hillslope[0].gw.DOC += (N_loss * patch[0].area / hillslope[0].area);
//        patch[0].cdf.DOC_to_gw += N_loss;
//        patch[0].soil_cs.DOC -= N_loss;
//    }
//    if (patch[0].soil_ns.nitrate > ZERO) {
//        N_loss = sat_to_gw_coeff * patch[0].soil_ns.nitrate*perc_inroot * command_line[0].Rsolute2gw;
//        hillslope[0].gw.NO3 += (N_loss * patch[0].area / hillslope[0].area);
//        patch[0].ndf.N_to_gw += N_loss;
//        patch[0].soil_ns.nitrate -= N_loss;
//    }
    
    
//    active_zone_z = (command_line[0].NH4root2active>0.0? patch[0].rootzone.depth * command_line[0].NH4root2active : (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active>0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)));
//    decay_rate = (command_line[0].NH4root2active>0.0? patch[0].soil_defaults[0][0].N_decay_rate : (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NH4decayRate : patch[0].soil_defaults[0][0].N_decay_rate));
    
//    perc_inroot = (1.0-exp(-decay_rate * patch[0].rootzone.depth)) / (1.0 - exp(-decay_rate * active_zone_z));
//    perc_inroot = max(0.1,min(perc_inroot,1.0));
//    sum_avail += max(0.0,perc_inroot * ns_soil->sminn); //sum_avail = perc_inroot * sum_avail;
	
	
	

	return (!ok);
} /* end update_gw_drainage.c */
