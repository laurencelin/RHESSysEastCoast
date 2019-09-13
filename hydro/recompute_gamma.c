
/*--------------------------------------------------------------*/
/* 								*/
/*		recompute_gamma			*/
/*								*/
/*	NAME							*/
/*	recompute_gamma - recomputes total and individual 	*/
/*	gamma values for subsurface routing to replace		*/
/*	tographic gradients by water table gradients		*/ 
/*								*/
/*								*/
/*	SYNOPSIS						*/
/*	recompute_gamma(					*/
/*			 struct  patch_object *,		*/
/*				double)				*/
/*								*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

double recompute_gamma(
                       struct patch_object *patch,
                       double total_gamma)
{
	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/ 
	/*--------------------------------------------------------------*/ 
	int i, d;
    double adjustment;
    double revised_total_gamma;
    double totaledge;
	double z1, z2, water_table_z1, water_table_z2;
	/*--------------------------------------------------------------*/ 
	/*	for now, if water table is above the surface we		*/
	/*	we set saturation deficit to zero, since return flow	*/
	/*	is modelled separately and should not be taken into	*/
	/*	account in modelling surface gradients			*/
	/*--------------------------------------------------------------*/ 

	adjustment = 0.0;
    totaledge = 0.0;
	z1 = patch[0].z;  
    if(patch[0].sat_deficit_z > ZERO){
		water_table_z1	 = (z1 - patch[0].sat_deficit_z);
    }else{
		water_table_z1 = z1;
    }
    
	d = 0;
    if (patch[0].innundation_list[d].num_neighbours > 0){
        
        // have neighbor
        revised_total_gamma = 0.0;
		for (i =0; i < patch[0].innundation_list[d].num_neighbours; i++) {
            
			z2 = patch[0].innundation_list[d].neighbours[i].patch[0].z;
            if(patch[0].innundation_list[d].neighbours[i].patch[0].sat_deficit_z > 0){
				water_table_z2	 = (z2 - patch[0].innundation_list[d].neighbours[i].patch[0].sat_deficit_z);
            }else{
				water_table_z2 = z2;
            }
            
            
            // old code
//            if (fabs(z1-z2) > ZERO) {
//                // z1 = current patch DEM
//                // z2 = neighbor DEM
//                // if there is a difference in elevation
//                adjustment += max(
//                           ((water_table_z1 - water_table_z2) / (z1 - z2) * patch[0].innundation_list[d].neighbours[i].gamma),
//                           0.0);
//            }else{
//                adjustment += 0.0;
//            }
            
            
            
            // new code, Sept 19, 2018
            if( (water_table_z1 - water_table_z2 - patch[0].innundation_list[d].neighbours[i].patch[0].constraintWaterTableTopDepth )>0 ){
                // updating the gamma fractions too
                patch[0].innundation_list[d].neighbours[i].gamma = (water_table_z1 - water_table_z2) * patch[0].innundation_list[d].neighbours[i].edgedistance; // (edge/distance)
                revised_total_gamma += patch[0].innundation_list[d].neighbours[i].gamma;
                totaledge += patch[0].innundation_list[d].neighbours[i].edge;
            }else{
                patch[0].innundation_list[d].neighbours[i].gamma = 0.0;
            }// end of if
            
		}//end of for neighbour i loop
        
        
        if(revised_total_gamma>0){
            
            for (i =0; i < patch[0].innundation_list[d].num_neighbours; i++) {
                patch[0].innundation_list[d].neighbours[i].gamma /= revised_total_gamma; // gamma fraction

                //patch[0].innundation_list[d].neighbours[i].gamma *= 0.5;
                //patch[0].innundation_list[d].neighbours[i].gamma += 0.5* patch[0].innundation_list[d].neighbours[i].gammaCONST;

            }//end of for neighbour i loop
           
            if(patch[0].drainage_type != STREAM){
                revised_total_gamma /= totaledge;
                revised_total_gamma *= patch[0].area;
                revised_total_gamma *= patch[0].soil_defaults[0][0].Ksat_0 * (patch[0].soil_defaults[0][0].m>0? patch[0].soil_defaults[0][0].m : 1.0);
                //revised_total_gamma *= patch[0].horizontal_k_SCALE;
                //revised_total_gamma *= patch[0].horizontal_k_SCALE;
            }else{
                //stream; no change at all
                revised_total_gamma = patch[0].tanSlope * patch[0].area;
                revised_total_gamma *= patch[0].soil_defaults[0][0].Ksat_0 * (patch[0].soil_defaults[0][0].m>0? patch[0].soil_defaults[0][0].m : 1.0);
                //revised_total_gamma *= patch[0].horizontal_k_SCALE;
                //revised_total_gamma *= patch[0].horizontal_k_SCALE;
            }
        }else{
            revised_total_gamma = 0.0;
        }//
        
        //revised_total_gamma *= 0.5;
        //revised_total_gamma += 0.5*adjustment * total_gamma; // average with the old code to bring directional flow prefer.
        
    }else{
        // no neighbor; num_neighbours = 0
		// adjustment = 1.0;
        revised_total_gamma = total_gamma;
    }// end of if
    
	//revised_total_gamma = adjustment * total_gamma;
		

	return(revised_total_gamma);
} /*recompute_gamma*/
