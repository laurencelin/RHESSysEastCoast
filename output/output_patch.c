/*--------------------------------------------------------------*/
/* 																*/
/*					output_patch						*/
/*																*/
/*	output_patch - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_patch - outputs current contents of a patch.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_patch(										*/
/*					struct	patch_object	*patch,				*/
/*					struct	date	date,  						*/
/*					FILE 	*outfile)							*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	outputs spatial structure according to commandline			*/
/*	specifications to specific files							*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	We only permit one fileset per spatial modelling level.     */
/*	Each fileset has one file for each timestep.  				*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

void	output_patch(
					 int basinID, int hillID, int zoneID,
					 struct	patch_object	*patch,
					 struct	zone_object	*zone,
					 struct	date	current_date,
					 FILE *outfile)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	int check, c, layer;
    
    double coverf;
    double alai = 0.0;
    double treeLAI = 0.0;
    double nontreeLAI = 0.0;
	for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
		for ( c=0 ; c<patch[0].layers[layer].count; c++ ){
            
            coverf = patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cover_fraction;
			alai += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].epv.proj_lai;
            
            if(patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].epc.veg_type == TREE && patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].ID!=802){
                treeLAI += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].epv.proj_lai;
            }//tree
            if(patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].epc.veg_type == GRASS || patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].ID==802){
                nontreeLAI += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].epv.proj_lai;
            }//non-tree
            
        }// for c
	}//for layer
    
    double wilting = exp(-1.0*log(-100.0*patch[0].psi_max_veg/patch[0].soil_defaults[0][0].psi_air_entry) * patch[0].soil_defaults[0][0].pore_size_index);
    double wilting_mm = min(patch[0].rootzone.potential_sat,patch[0].sat_deficit) * wilting;
    double totalfc = patch[0].rootzone.field_capacity + patch[0].field_capacity;
    double vksat0 = patch[0].soil_defaults[0][0].Ksat_0_v;
    double vksat_decay = patch[0].soil_defaults[0][0].mz_v;
    double vksat_decay_1 = -1.0/patch[0].soil_defaults[0][0].mz_v;
    
    double top12cm_storage;
    double top12cm_potential_sat = patch[0].soil_defaults[0][0].rtz2sat_def_0z[120];
    double top30cm_storage;
    double top30cm_potential_sat = patch[0].soil_defaults[0][0].rtz2sat_def_0z[300];
    double top60cm_storage;
    double top60cm_potential_sat = patch[0].soil_defaults[0][0].rtz2sat_def_0z[600];
    
    
    if(patch[0].rootzone.potential_sat>0){
        // plants and root
        if(patch[0].rz_storage > patch[0].rootzone.field_capacity){
            // drainage pattern
            // vksat_decay* vksat0*(1.0-exp(vksat_decay_1*0.12))
            // vksat_decay* vksat0*(1.0-exp(vksat_decay_1*patch[0].rootzone.depth))
            top12cm_storage = totalfc * patch[0].soil_defaults[0][0].fc1_00012r[patch[0].sat_def_pct_index];
            top12cm_storage += (patch[0].rz_storage-patch[0].rootzone.field_capacity) * (1.0 - (1.0-exp(vksat_decay_1*0.12))/(1.0-exp(vksat_decay_1*patch[0].rootzone.depth))); // approximation
            
            top30cm_storage = totalfc * patch[0].soil_defaults[0][0].fc1_003r[patch[0].sat_def_pct_index];
            top30cm_storage += (patch[0].rz_storage-patch[0].rootzone.field_capacity) * (1.0 - (1.0-exp(vksat_decay_1*0.30))/(1.0-exp(vksat_decay_1*patch[0].rootzone.depth))); // approximation
            
            top60cm_storage = totalfc * patch[0].soil_defaults[0][0].fc1_006r[patch[0].sat_def_pct_index];
            top60cm_storage += (patch[0].rz_storage-patch[0].rootzone.field_capacity) * (1.0 - (1.0-exp(vksat_decay_1*0.60))/(1.0-exp(vksat_decay_1*patch[0].rootzone.depth))); // approximation
            
        }else if(patch[0].rz_storage < patch[0].rootzone.field_capacity){
            //bounded by wilting + (patch[0].rz_storage - wilting) follows ?
            top12cm_storage = patch[0].sat_deficit>0.0? wilting_mm * top12cm_potential_sat/min(patch[0].rootzone.potential_sat,patch[0].sat_deficit) : 0.0;
            if( patch[0].rz_storage > wilting_mm && patch[0].sat_deficit>0.0){
                top12cm_storage += (patch[0].rz_storage-wilting_mm) * totalfc * patch[0].soil_defaults[0][0].fc1_00012r[patch[0].sat_def_pct_index] / patch[0].rootzone.field_capacity;
            }
            
            top30cm_storage = patch[0].sat_deficit>0.0? wilting_mm * top30cm_potential_sat/min(patch[0].rootzone.potential_sat,patch[0].sat_deficit) : 0.0;
            if( patch[0].rz_storage > wilting_mm && patch[0].sat_deficit>0.0){
                top30cm_storage += (patch[0].rz_storage-wilting_mm) * totalfc * patch[0].soil_defaults[0][0].fc1_003r[patch[0].sat_def_pct_index] / patch[0].rootzone.field_capacity;
            }
            
            top60cm_storage = patch[0].sat_deficit>0.0? wilting_mm * top60cm_potential_sat/min(patch[0].rootzone.potential_sat,patch[0].sat_deficit) : 0.0;
            if( patch[0].rz_storage > wilting_mm && patch[0].sat_deficit>0.0){
                top60cm_storage += (patch[0].rz_storage-wilting_mm) * totalfc * patch[0].soil_defaults[0][0].fc1_006r[patch[0].sat_def_pct_index] / patch[0].rootzone.field_capacity;
            }
            
        }else{
            // patch[0].rz_storage = patch[0].rootzone.field_capacity
            top12cm_storage = totalfc * patch[0].soil_defaults[0][0].fc1_00012r[patch[0].sat_def_pct_index];
            top30cm_storage = totalfc * patch[0].soil_defaults[0][0].fc1_003r[patch[0].sat_def_pct_index];
            top60cm_storage = totalfc * patch[0].soil_defaults[0][0].fc1_006r[patch[0].sat_def_pct_index];
        }
    }else{
        // no root / veg
        if(patch[0].unsat_storage > patch[0].field_capacity){
            top12cm_storage = totalfc * patch[0].soil_defaults[0][0].fc1_00012r[patch[0].sat_def_pct_index];
            top12cm_storage += (patch[0].unsat_storage-patch[0].field_capacity) * (1.0 - (1.0-exp(vksat_decay_1*0.12))/(1.0-exp(vksat_decay_1*patch[0].sat_deficit_z))); // approximation
            
            top30cm_storage = totalfc * patch[0].soil_defaults[0][0].fc1_003r[patch[0].sat_def_pct_index];
            top30cm_storage += (patch[0].unsat_storage-patch[0].field_capacity) * (1.0 - (1.0-exp(vksat_decay_1*0.30))/(1.0-exp(vksat_decay_1*patch[0].sat_deficit_z))); // approximation
            
            top60cm_storage = totalfc * patch[0].soil_defaults[0][0].fc1_006r[patch[0].sat_def_pct_index];
            top60cm_storage += (patch[0].unsat_storage-patch[0].field_capacity) * (1.0 - (1.0-exp(vksat_decay_1*0.60))/(1.0-exp(vksat_decay_1*patch[0].sat_deficit_z))); // approximation
        }else{
            top12cm_storage = totalfc * patch[0].soil_defaults[0][0].fc1_00012r[patch[0].sat_def_pct_index];
            top30cm_storage = totalfc * patch[0].soil_defaults[0][0].fc1_003r[patch[0].sat_def_pct_index];
            top60cm_storage = totalfc * patch[0].soil_defaults[0][0].fc1_006r[patch[0].sat_def_pct_index];
        }
    }//if else
    

    
	check = fprintf(outfile,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                    
					current_date.year, current_date.month, current_date.day, //1,2,3,
					patch[0].ID, //4
                    (patch[0].Qout_total - patch[0].Qin_total) * 1000.0, //5
                    (patch[0].surface_Qout_total - patch[0].surface_Qin_total) * 1000.0, //6
                    patch[0].detention_store*1000.0, //7
                    patch[0].stormdrainYield*1000.0, //8
                    patch[0].overland_flow * 1000.0, //9 <-- locally yielded returnflow (not all made to the stream); not include surface Qin and Qout
                    patch[0].rain_throughfall*1000.0, //10
					(patch[0].rain_throughfall - patch[0].recharge)*1000.0,//11
                    (patch[0].cap_rise - patch[0].unsat_drainage)*1000.0,//12
                    patch[0].sat_deficit_z*1000.0,//13 wtz
                    patch[0].sat_deficit*1000.0, //14 total subS -> sat_def
                    patch[0].rz_storage*1000.0, //15 rtS -> rtz storage
                        
                    (patch[0].transpiration_sat_zone + patch[0].transpiration_unsat_zone + patch[0].evaporation + patch[0].evaporation_surf  + patch[0].exfiltration_sat_zone + patch[0].exfiltration_unsat_zone)*1000.0, //16 ET mm
                    
                    treeLAI, //17
                    nontreeLAI, //18
                    patch[0].grassIrrigation_m,
                    
                    patch[0].rootzone.potential_sat*1000.0,
                    patch[0].field_capacity*1000.0,
                    patch[0].rootzone.field_capacity*1000.0,
                    patch[0].unsat_storage*1000.0,
                    top12cm_storage * 1000.0,
                    top12cm_potential_sat * 1000.0,
                    patch[0].rootzone.depth * 1000.0,
                    patch[0].soil_defaults[0][0].soil_depth * 1000.0,
                    top30cm_storage * 1000.0,
                    top30cm_potential_sat * 1000.0,
                    top60cm_storage * 1000.0,
                    top60cm_potential_sat * 1000.0
                    );

	if (check <= 0) {
		fprintf(stdout, "\nWARNING: output error has occured in output_patch, file");
	}
	return;
} /*end output_patch*/
