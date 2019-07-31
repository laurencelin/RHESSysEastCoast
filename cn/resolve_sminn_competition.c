/*--------------------------------------------------------------*/
/* 								*/
/*								*/
/*	resolve_sminn_competition				*/
/*								*/
/*	NAME							*/
/*		resolve_sminn_competition			*/
/*								*/
/*	SYNOPSIS						*/
/*	int resolve_sminn_competition(				*/
/*   			    struct  soil_n_object   *ns_soil,	*/
/*     			    struct ndayflux_patch_struct *ndf)  */
/*								*/
/*	returns:						*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*								*/
/*              modified from Peter Thornton (1998)             */
/*                      dynamic - 1d-bgc ver4.0                 */
/*		plant uptake under limiting nitrogen is 	*/
/*		now a function of root density			*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

int resolve_sminn_competition(
							  struct  soil_n_object   *ns_soil,
							  struct patch_object *patch,
							  struct command_line_object *command_line,
							  struct ndayflux_patch_struct *ndf)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	double sum_ndemand, sum_avail, litterndemand, totalDecompDemand;
	double actual_immob, actual_uptake, perc_inroot;
    double active_zone_z, decay_rate;
    double p_decayRate = 1.0 / patch[0].soil_defaults[0][0].porosity_decay;
    double z1 = patch[0].sat_deficit_z>0? patch[0].sat_deficit_z : 0.0; //<<--- count for negative
    double active_zone_z_sat = max(z1 + 0.33, patch[0].rootzone.depth);
    double constantHold1, constantHold2, constantHold3;
    constantHold1 = exp(-z1*p_decayRate);
    constantHold2 = exp(-patch[0].rootzone.depth*p_decayRate);
    constantHold3 = exp(-active_zone_z_sat*p_decayRate);
    
    //double rtzNO3, rtzSatNO3, rtzNH4, rtzSatNH4;
    
	/*--------------------------------------------------------------*/
	/* compare the combined decomposition immobilization and plant*/
	/*growth N demands against the available soil mineral N pool. */
	/*--------------------------------------------------------------*/
    litterndemand = max(0.0, ndf->potential_immoblitter - patch[0].surface_NO3 + patch[0].litter.NO3_stored - patch[0].surface_NH4);
    totalDecompDemand = ndf->potential_immob - ndf->potential_immoblitter + litterndemand;//N demand for the subsurface N pools
	sum_ndemand = ndf->plant_potential_ndemand + totalDecompDemand;
	//sum_avail = max(ns_soil->sminn + ns_soil->nitrate + ndf->mineralized, 0.0);
    sum_avail = max(ndf->mineralized, 0.0);
	/*--------------------------------------------------------------*/
	/* limit available N for plants by rooting depth		*/
	/* for really small rooting depths this can be problematic	*/
	/* for now provide a minimum access				*/
	/*--------------------------------------------------------------*/
    active_zone_z = (command_line[0].NH4root2active>0.0? patch[0].rootzone.depth * command_line[0].NH4root2active : (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active>0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)));
    decay_rate = (command_line[0].NH4root2active>0.0? patch[0].soil_defaults[0][0].N_decay_rate : (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NH4decayRate : patch[0].soil_defaults[0][0].N_decay_rate));
    perc_inroot = (1.0-exp(-decay_rate * patch[0].rootzone.depth)) / (1.0 - exp(-decay_rate * active_zone_z));
    perc_inroot = min(perc_inroot,1.0);
    patch[0].rtzNH4 = max(0.0,perc_inroot * ns_soil->sminn);
    patch[0].rtzSatNH4 = z1 < patch[0].rootzone.depth? (constantHold1-constantHold2)/(constantHold1-constantHold3)*patch[0].sat_NH4 : 0.0;
    
    active_zone_z = (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z));
    decay_rate = (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NO3decayRate : patch[0].soil_defaults[0][0].N_decay_rate);
    perc_inroot = (1.0-exp(-decay_rate * (patch[0].rootzone.depth+patch[0].ndistz))) / (1.0 - exp(-decay_rate * (active_zone_z+patch[0].ndistz)));
    perc_inroot = min(perc_inroot,1.0);
    patch[0].rtzNO3 = max(0.0,perc_inroot * ns_soil->nitrate); //sum_avail = perc_inroot * sum_avail;
    patch[0].rtzSatNO3 = z1 < patch[0].rootzone.depth? (constantHold1-constantHold2)/(constantHold1-constantHold3)*patch[0].sat_NO3 : 0.0;
    
    sum_avail = patch[0].rtzNO3 + patch[0].rtzSatNO3 + patch[0].rtzNH4 + patch[0].rtzSatNH4;
    
	if (sum_ndemand <= sum_avail){
	/* N availability is not limiting immobilization or plant
		uptake, and both can proceed at their potential rates */
		actual_immob = ndf->potential_immob;
		ns_soil->nlimit = 0;
		ns_soil->fract_potential_immob = 1.0;
		ns_soil->fract_potential_uptake = 1.0;
		ndf->plant_avail_uptake = ndf->plant_potential_ndemand;
	}
	else{
	/* N availability can not satisfy the sum of immobiliation and
	plant growth demands, so these two demands compete for available
		soil mineral N */
        // sum_avail,
        
        ns_soil->nlimit = 1;
		actual_immob = sum_avail * (ndf->potential_immob - ndf->potential_immoblitter + litterndemand)/sum_ndemand;// microbe first
		actual_uptake = max(sum_avail - actual_immob,0.0); // <<--- could be small negative!!
		
        // set limit for microbe
        if (ndf->potential_immob == 0 && totalDecompDemand == 0){
			ns_soil->fract_potential_immob = 0.0; // absolutely zero immobilization
        }else if (ndf->potential_immob > 0 && totalDecompDemand == 0){
            ns_soil->fract_potential_immob = 1.0; // soilOM mineralizing while surface litter immobilizing from the surface N sources.
        }else{
            ns_soil->fract_potential_immob = actual_immob/totalDecompDemand;
        }
        
        // set limit for plant uptake
		if (ndf->plant_potential_ndemand == 0) {
			ns_soil->fract_potential_uptake = 0.0;
			ndf->plant_avail_uptake = 0.0; 
		}else {
			ns_soil->fract_potential_uptake = actual_uptake / ndf->plant_potential_ndemand;
			ndf->plant_avail_uptake = actual_uptake; // <<--- could be small negative!!
		}
        
        // fract_potential_immob   is very important;  potential_immoblitter
//        if(ndf->plant_potential_ndemand>0.0) printf(" avail %e = (%e + %e + %e)=%e   demand %e = (%e + %e, %e) -> reduction (immob %f, uptake %f, %d) perc_inroot = %f, nflux(%e,%e,%e,%e, %e,%e,%e,%e)\n",
//               sum_avail, ns_soil->sminn, ns_soil->nitrate, ndf->mineralized, (ns_soil->sminn+ ns_soil->nitrate+ ndf->mineralized),
//               sum_ndemand, ndf->plant_potential_ndemand, totalDecompDemand, ndf->potential_immob,
//               ns_soil->fract_potential_immob, ns_soil->fract_potential_uptake,ns_soil->nlimit,
//               perc_inroot,
//               ndf->pmnf_l1s1, ndf->pmnf_l2s2, ndf->pmnf_l3l2, ndf->pmnf_l4s3, ndf->pmnf_s1s2, ndf->pmnf_s2s3, ndf->pmnf_s3s4, ndf->pmnf_s4);
        
	}// nlimit

	return(0);
} /* end resolve_sminn_competition.c */

