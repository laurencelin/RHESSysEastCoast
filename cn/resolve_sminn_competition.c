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
	double sum_ndemand, litterndemand, totalDecompDemand;
	double actual_immob, actual_uptake, perc_inroot;
    //double active_zone_z, decay_rate;
    double p_decayRate = 1.0 / patch[0].soil_defaults[0][0].porosity_decay;
    double z1 = patch[0].sat_deficit_z>0? patch[0].sat_deficit_z : 0.0; //<<--- count for negative
    
    //double rtzNO3, rtzSatNO3, rtzNH4, rtzSatNH4;
    
	/*--------------------------------------------------------------*/
	/* compare the combined decomposition immobilization and plant*/
	/*growth N demands against the available soil mineral N pool. */
	/*--------------------------------------------------------------*/
    
    // demands:
    // 1) ndf->plant_potential_ndemand
    // 2) ndf->potential_immob-ndf->potential_immoblitter
    // 3) ndf->potential_immoblitter
    litterndemand = max(0.0, ndf->potential_immoblitter - patch[0].surface_NO3 + patch[0].litter.NO3_stored + patch[0].surface_NH4);
    totalDecompDemand = (ndf->potential_immob - ndf->potential_immoblitter) + litterndemand;//N demand for the subsurface N pools
	sum_ndemand = ndf->plant_potential_ndemand + totalDecompDemand;
	
    
    // resources
    // 1) patch[0].surface_NO3 + patch[0].litter.NO3_stored + patch[0].surface_NH4 :: ndf->potential_immoblitter
    // 2) max(ndf->mineralized, 0.0) :: (ndf->potential_immob-ndf->potential_immoblitter)
    // 3) patch[0].rtzNO3 + patch[0].rtzSatNO3 + patch[0].rtzNH4 + patch[0].rtzSatNH4 :: ndf->plant_potential_ndemand
    
    // we have 1) patch[0].soil_defaults[0][0].active_zone_index; 2) patch[0].soil_defaults[0][0].maxrootdepth_index; and 3) patch[0].rtz2_index
    // surface layer: litter only
    // 0-rtz: plant uptake first
    // mineralized: decomp first
    
    // important note:: "ndf->plant_avail_uptake" will be updated to "ndf->sminn_to_npool"
    // important note:: "update_decomp.c" is the code does the N substraction to soilN based on immobilization and uptake; it also does the mineralization.
    
    
    double avail_rtzN = 0.0;
    patch[0].rtzNH4 = max(0.0, ns_soil->sminn*patch[0].soil_defaults[0][0].rtz2NH4prop[patch[0].rtz2_index]);
    patch[0].rtzSatNH4 = 0.0;
    if(patch[0].available_soil_water>0){
        patch[0].rtzSatNH4 = patch[0].sat_NH4 * max(patch[0].rootzone.potential_sat-patch[0].sat_deficit,0.0)/patch[0].available_soil_water;
    }//if
    patch[0].rtzNO3 = max(0.0, ns_soil->nitrate*patch[0].soil_defaults[0][0].rtz2NO3prop[patch[0].rtz2_index]);
    patch[0].rtzSatNO3 = 0.0;
    if(patch[0].available_soil_water>0){
        patch[0].rtzSatNO3 = patch[0].sat_NO3 * max(patch[0].rootzone.potential_sat-patch[0].sat_deficit,0.0)/patch[0].available_soil_water;
    }//if
    avail_rtzN = patch[0].rtzNO3 + patch[0].rtzSatNO3 + patch[0].rtzNH4 + patch[0].rtzSatNH4;
    
    // assume 95% soil OM is within 0 - maxrootz
    double avail_activeN = 0.0; // no overlapping with the "avail_rtzN"
    patch[0].rtzNH4 = max(0.0, ns_soil->sminn*patch[0].soil_defaults[0][0].rtz2NH4prop[patch[0].soil_defaults[0][0].maxrootdepth_index]);
    patch[0].rtzSatNH4 = 0.0;
    if(patch[0].available_soil_water>0){
        patch[0].rtzSatNH4 = patch[0].sat_NH4 * max(patch[0].soil_defaults[0][0].maxrootdepth_sat_0z-patch[0].sat_deficit,0.0)/patch[0].available_soil_water;
    }//if
    patch[0].rtzNO3 = max(0.0, ns_soil->nitrate*patch[0].soil_defaults[0][0].rtz2NO3prop[patch[0].soil_defaults[0][0].maxrootdepth_index]);
    patch[0].rtzSatNO3 = 0.0;
    if(patch[0].available_soil_water>0){
        patch[0].rtzSatNO3 = patch[0].sat_NO3 * max(patch[0].soil_defaults[0][0].maxrootdepth_sat_0z-patch[0].sat_deficit,0.0)/patch[0].available_soil_water;
    }//if
    avail_activeN = patch[0].rtzNO3 + patch[0].rtzSatNO3 + patch[0].rtzNH4 + patch[0].rtzSatNH4;
    avail_activeN -= avail_rtzN;
//    if(patch[0].soil_defaults[0][0].active_zone_index > patch[0].rtz2_index){
//        patch[0].rtzNH4 = max(0.0, ns_soil->sminn*patch[0].soil_defaults[0][0].rtz2NH4prop[patch[0].soil_defaults[0][0].active_zone_index]);
//        patch[0].rtzSatNH4 = 0.0;
//        if(patch[0].available_soil_water>0){
//            patch[0].rtzSatNH4 = patch[0].sat_NH4 * max(patch[0].rootzone.potential_sat-patch[0].sat_deficit,0.0)/patch[0].available_soil_water;
//        }//if
//        patch[0].rtzNO3 = max(0.0, ns_soil->nitrate*patch[0].soil_defaults[0][0].rtz2NO3prop[patch[0].soil_defaults[0][0].active_zone_index]);
//        patch[0].rtzSatNO3 = 0.0;
//        if(patch[0].available_soil_water>0){
//            patch[0].rtzSatNO3 = patch[0].sat_NO3 * max(patch[0].rootzone.potential_sat-patch[0].sat_deficit,0.0)/patch[0].available_soil_water;
//        }//if
//        avail_activeN = patch[0].rtzNO3 + patch[0].rtzSatNO3 + patch[0].rtzNH4 + patch[0].rtzSatNH4;
//        avail_activeN -= avail_rtzN;
//    }// end of if
    
    
    double avail_mineral = max(ndf->mineralized, 0.0);
    

    
    // this is wrong below!!!
    if(ndf->plant_potential_ndemand > avail_rtzN){
        // plant needs additional N
        if(totalDecompDemand > avail_activeN + avail_mineral){
            // both are n limited
            ns_soil->nlimit = 1;
            ns_soil->fract_potential_immob = (avail_activeN + avail_mineral)/totalDecompDemand;
            ndf->plant_avail_uptake = avail_rtzN;
            ns_soil->fract_potential_uptake = ndf->plant_avail_uptake/ndf->plant_potential_ndemand;
        }else{
            // no n limit for decomp
            ns_soil->fract_potential_immob = 1.0;
            
            ndf->plant_avail_uptake = avail_rtzN + (avail_mineral-max(0.0,totalDecompDemand-avail_activeN))*patch[0].soil_defaults[0][0].rtz2NH4prop[patch[0].rtz2_index];
            if( ndf->plant_avail_uptake >= ndf->plant_potential_ndemand){
                // no n limit
                ns_soil->nlimit = 0;
                ndf->plant_avail_uptake = ndf->plant_potential_ndemand;
                ns_soil->fract_potential_uptake = 1.0;
            }else{
                // yes plant n limit
                ns_soil->nlimit = 1;
                ns_soil->fract_potential_uptake = ndf->plant_avail_uptake/ndf->plant_potential_ndemand;
            }
        }
    }else{
        // no n limit for plant uptake
        ns_soil->fract_potential_uptake = 1.0;
        ndf->plant_avail_uptake = ndf->plant_potential_ndemand;
        
        // decomp
        double remainN = avail_activeN + (avail_rtzN - ndf->plant_potential_ndemand) + avail_mineral;
        if(totalDecompDemand > remainN){
            //decomp needs additional N
            ns_soil->nlimit = 1;
            ns_soil->fract_potential_immob = remainN/totalDecompDemand;
        }else{
            // no n limit for decomp
            ns_soil->nlimit = 0;
            ns_soil->fract_potential_immob = 1.0;
        }
    }//end of if
    
    if( (ns_soil->fract_potential_immob!=ns_soil->fract_potential_immob) || ns_soil->fract_potential_immob<0 ||
    (ns_soil->fract_potential_uptake!=ns_soil->fract_potential_uptake) || ns_soil->fract_potential_uptake<0 )
    printf("error resolve_sminn_competition:%d, %e %e, %d: %e,%e,%e, (%e - %e)\n",
           patch[0].ID,
           ns_soil->fract_potential_immob, ns_soil->fract_potential_uptake, ns_soil->nlimit,
           //-----------
           litterndemand, totalDecompDemand, sum_ndemand,
           ndf->potential_immob, ndf->potential_immoblitter
           );
    
    
    
    
    
//	if (sum_ndemand <= sum_avail){
//	/* N availability is not limiting immobilization or plant
//		uptake, and both can proceed at their potential rates */
//		actual_immob = ndf->potential_immob;
//		ns_soil->nlimit = 0;
//		ns_soil->fract_potential_immob = 1.0;
//		ns_soil->fract_potential_uptake = 1.0;
//		ndf->plant_avail_uptake = ndf->plant_potential_ndemand;
//	}else{
//	/* competition */
//
//        ns_soil->nlimit = 1;
//        // litterndemand has been taken care of
//        // totalDecompDemand = (ndf->potential_immob - ndf->potential_immoblitter) + litterndemand;
//
//        actual_uptake = sum_avail * max(1.0 - totalDecompDemand/sum_ndemand, 0.0); //microbe first
//        if (ndf->plant_potential_ndemand > 0) {
//            ns_soil->fract_potential_uptake = actual_uptake / ndf->plant_potential_ndemand;
//            ndf->plant_avail_uptake = actual_uptake;
//        }else{
//            ns_soil->fract_potential_uptake = 0.0;
//            ndf->plant_avail_uptake = 0.0;
//        }
//
//
//		actual_immob = max(0.0, sum_avail - actual_uptake);// microbe first
//        if(totalDecompDemand > 0){
//            ns_soil->fract_potential_immob = actual_immob/totalDecompDemand;
//            // affect only the immobilizing fluxes; no effect on mineralizing fluxes
//        }else{
//            // totalDecompDemand = 0 --> ndf->potential_immob = ndf->potential_immoblitter && litterndemand = 0
//            ns_soil->fract_potential_immob = 1.0;
//        }
//
//
//        if( (ns_soil->fract_potential_immob!=ns_soil->fract_potential_immob) || ns_soil->fract_potential_immob<0 ||
//            (ns_soil->fract_potential_uptake!=ns_soil->fract_potential_uptake) || ns_soil->fract_potential_uptake<0 )
//            printf("error resolve_sminn_competition:%d, %e %e, %d: %e,%e,%e,%e (%e - %e)\n",
//                   patch[0].ID,
//                   ns_soil->fract_potential_immob, ns_soil->fract_potential_uptake, ns_soil->nlimit,
//                   //-----------
//                   litterndemand, totalDecompDemand, sum_ndemand, sum_avail,
//                   ndf->potential_immob, ndf->potential_immoblitter
//                   );
//

        // fract_potential_immob   is very important;  potential_immoblitter
//        if(ndf->plant_potential_ndemand>0.0) printf(" avail %e = (%e + %e + %e)=%e   demand %e = (%e + %e, %e) -> reduction (immob %f, uptake %f, %d) perc_inroot = %f, nflux(%e,%e,%e,%e, %e,%e,%e,%e)\n",
//               sum_avail, ns_soil->sminn, ns_soil->nitrate, ndf->mineralized, (ns_soil->sminn+ ns_soil->nitrate+ ndf->mineralized),
//               sum_ndemand, ndf->plant_potential_ndemand, totalDecompDemand, ndf->potential_immob,
//               ns_soil->fract_potential_immob, ns_soil->fract_potential_uptake,ns_soil->nlimit,
//               perc_inroot,
//               ndf->pmnf_l1s1, ndf->pmnf_l2s2, ndf->pmnf_l3l2, ndf->pmnf_l4s3, ndf->pmnf_s1s2, ndf->pmnf_s2s3, ndf->pmnf_s3s4, ndf->pmnf_s4);
        
//	}// nlimit

	return(0);
} /* end resolve_sminn_competition.c */

