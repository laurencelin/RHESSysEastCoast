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
    double alai; //, asub, apsn, litterS, aheight, agsi;
    double coverf;
//    double mean_gl, gl_scalar;
//    double coverf, totalcoverf;
//    struct mult_conduct_struct *mc;
//    double mc_scalar;
    
//    if (patch[0].litter.rain_capacity > ZERO)
//        litterS = patch[0].litter.rain_stored / patch[0].litter.rain_capacity;
//    else
//        litterS = 1.0;

//    apsn = 0.0;
//    asub = 0.0;
	alai = 0.0;
//    aheight = 0.0;
//    agsi = 0.0;
//    mean_gl = 0.0;
//    gl_scalar= 0.0;
//    totalcoverf = 0.0;
    double treeLAI = 0.0;
    double nontreeLAI = 0.0;
	for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
		for ( c=0 ; c<patch[0].layers[layer].count; c++ ){
            
            coverf = patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cover_fraction;
//            mc = &patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].mult_conductance;
//            mc_scalar = mc->APAR * mc->tavg * mc->LWP * mc->CO2 * mc->tmin * mc->vpd;
            
//            apsn += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cs.net_psn ;
//            asub += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].sublimation;
//            aheight += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].epv.height;
			alai += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].epv.proj_lai;
            
            if(patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].epc.veg_type == TREE && patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].ID!=802){
                treeLAI += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].epv.proj_lai;
            }//tree
            if(patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].epc.veg_type == GRASS || patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].ID==802){
                nontreeLAI += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].epv.proj_lai;
            }//non-tree
            
            

            
//            agsi += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].phen.gsi;
//            mean_gl += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].epc.gl_smax * mc_scalar;
//            gl_scalar += coverf * mc_scalar;
//            totalcoverf += coverf;
//
		}
	}
//    if(totalcoverf>1.0){ mean_gl/=totalcoverf; gl_scalar/=totalcoverf;  }
    
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
    
//    default_object_list[i].vksat_0zm[ii] = default_object_list[i].sat_def_z[ii]>0? vksat_decay* vksat0*(1-exp(vksat_decay_1*default_object_list[i].sat_def_z[ii]))/default_object_list[i].sat_def_z[ii] : vksat0;
//    patch[0].sat_deficit
//    patch[0].rootzone.potential_sat
//    patch[0].soil_defaults[0][0].rtz2sat_def_0z[120] / patch[0].rootzone.potential_sat; //wilting prop
//    patch[0].psi_max_veg
//    patch[0].sat_def_pct_index
//    patch[0].rtz2_index = (int)(round(patch[0].rootzone.depth*1000));
//    patch[0].rootzone.potential_sat = patch[0].soil_defaults[0][0].rtz2sat_def_0z[patch[0].rtz2_index];
//    patch[0].rootzone.depth
//    patch[0].soil_defaults[0][0].soil_depth
    
	check = fprintf(outfile,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", // \
                             //%lf %lf %lf %lf %lf \
                             //%lf %lf %lf %lf %lf \
                             //%lf %lf %lf\n",
                                // added 4 extra to track fc
					current_date.year, current_date.month, current_date.day, //1,2,3,
					patch[0].ID, //4
                    (patch[0].Qout_total - patch[0].Qin_total) * 1000.0, //5
                    (patch[0].surface_Qout_total - patch[0].surface_Qin_total) * 1000.0, //6
                    patch[0].detention_store*1000.0, //7
                    patch[0].stormdrainYield*1000.0, //8
                    patch[0].overland_flow * 1000.0, //9 <-- local yielded returnflow (not all made to the stream); not include surface Qin and Qout
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
 
    
                    //patch[0].sat_deficit*1000.0, //patch[0].sat_DOC*1000.0,
                    //patch[0].rootzone.potential_sat*1000.0, //patch[0].soil_ns.nitrate*1000.0,
                    
                    //---extra
//                    patch[0].sat_NO3*1000.0,
//                    patch[0].sat_NH4*1000.0,
//                    patch[0].soil_ns.sminn*1000.0,
//                    patch[0].soil_cs.DOC,
//                    patch[0].z,
//                    patch[0].ndf.denitrif*1000.0,
//                    patch[0].ndf.sminn_to_nitrate*1000.0,
//                    patch[0].ndf.plant_avail_uptake*1000.0, // plant uptake?
//                    patch[0].ndf.net_mineralized*1000.0, // net from decay
//                    (patch[0].ndf.net_mineralized - patch[0].ndf.mineralized)*1000.0, // decay uptake?
//                    patch[0].soil_ns.NO3_Qout_total*1000.0,
//                    patch[0].soil_ns.NO3_Qin_total*1000.0,
//                    patch[0].sat_deficit_z>0? patch[0].rootzone.depth/patch[0].sat_deficit_z : 1000.0
//                    );
	

	if (check <= 0) {
		fprintf(stdout, "\nWARNING: output error has occured in output_patch, file");
	}
	return;
} /*end output_patch*/
