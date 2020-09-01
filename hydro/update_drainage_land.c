/*--------------------------------------------------------------*/
/* 											*/
/*					update_drainage_land			*/
/*											*/
/*	update_drainage_land.c - creates a patch object				*/
/*											*/
/*	NAME										*/
/*	update_drainage_land.c - creates a patch object				*/
/*											*/
/*	SYNOPSIS									*/
/*	void update_drainage_land( 							*/
/*					struct patch_object *patch			*/
/*				 			double,			 	*/
/*				 			double,			 	*/
/*				 			double,			 	*/
/*							int,				*/
/*							int)				*/
/*											*/
/* 											*/
/*											*/
/*	OPTIONS										*/
/*											*/
/*											*/
/*	DESCRIPTION									*/
/*											*/
/*											*/
/*											*/
/*											*/
/*	PROGRAMMER NOTES								*/
/*											*/
/*											*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"


void  update_drainage_land(
					struct patch_object *patch,
					 struct command_line_object *command_line,
					 double time_int,
					 int kk)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.				*/
	/*--------------------------------------------------------------*/
    int verbose_flag = 0;
	double  compute_delta_water(
		int,
		double,
		double,
		double,
		double,
		double);
	
	double compute_varbased_returnflow(
		double,
		double,
		double,
		struct litter_object *);


	double compute_varbased_flow(
		int,
		double,
		int,
		double,
		double,
		struct patch_object *);


	double compute_N_leached(
         int verbose_flag,
         double total_nitrate,
         double Qout,
         double N_decay_rate,
         double activedepthz,
         double N_absorption_rate,
         int signal,
         struct patch_object *patch);
		
	
	double recompute_gamma(	
		struct patch_object *,
		double);


	double compute_infiltration( int,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double);
	
	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
	int j, d;
	//double m, Ksat, std_scale;
    double std_scale;
	double NH4_leached_to_patch;
	double NO3_leached_to_patch;
	double DON_leached_to_patch;
	double DOC_leached_to_patch;
	double NO3_leached_to_surface; /* kg/m2 */
	double NH4_leached_to_surface; /* kg/m2 */
	double DON_leached_to_surface; /* kg/m2 */
	double DOC_leached_to_surface; /* kg/m2 */

	double route_to_surface;  /* m3 */
	double return_flow, route_to_patch ;  /* m3 */
	double available_sat_water; /* m3 */
	double Qin, Qout;  /* m */
	double innundation_depth, infiltration; /* m */
	double total_gamma;
	double Nin, Nout; /* kg/m2 */ 
    double extrawater;

	struct patch_object *neigh;
	route_to_patch = 0.0;
	route_to_surface = 0.0;
	return_flow=0.0;
    extrawater = 0.0;
    
	DON_leached_to_patch = 0.0;
	DOC_leached_to_patch = 0.0;
	NH4_leached_to_patch = 0.0;
	NO3_leached_to_patch = 0.0;
	NO3_leached_to_surface = 0.0;
	NH4_leached_to_surface = 0.0;
	DOC_leached_to_surface = 0.0;
	DON_leached_to_surface = 0.0;
    
    if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
       patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
       patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
       patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
       patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
       patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
       patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
       patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
       patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0){
        printf("update_drainage_land0 (%d,%d) [%e %e %e] [%e %e %e] [%e %e %e] %e %e[%d] [%e %e %e %e]\n",
               patch[0].ID, kk,
               patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
               patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
               patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC,
               patch[0].sat_deficit_z, patch[0].sat_deficit, patch[0].sat_def_pct_index,
               patch[0].soil_ns.NO3_Qin,patch[0].soil_ns.NH4_Qin,patch[0].soil_cs.DOC_Qin,patch[0].Qin);
    }//debug
    
    
    int sat_def_z_index = (int)(max(0.0, patch[0].sat_deficit_z*1000));
    double rate_;
    
    // soil solute transfer to sat solute because of change in water table depth; this is different from drainage or cap_rise
    rate_ = patch[0].soil_defaults[0][0].rtz2NO3prop[sat_def_z_index];
    patch[0].sat_NO3 += patch[0].soil_ns.nitrate * (1.0-rate_); // o_z / o_Z
    patch[0].soil_ns.nitrate *=  rate_; // this rate_ is negative?
    
    rate_ = patch[0].soil_defaults[0][0].rtz2NH4prop[sat_def_z_index];
    patch[0].sat_NH4 += patch[0].soil_ns.sminn * (1.0-rate_);
    patch[0].soil_ns.sminn *= rate_;
    
    rate_ = patch[0].soil_defaults[0][0].rtz2DOMprop[sat_def_z_index];
    patch[0].sat_DOC += patch[0].soil_cs.DOC * (1.0-rate_);
    patch[0].sat_DON += patch[0].soil_ns.DON * (1.0-rate_);
    patch[0].soil_cs.DOC *= rate_;
    patch[0].soil_ns.DON *= rate_;
    
    //    if(patch[0].soil_ns.nitrate<0 || rate_<0){
    //        printf("update_drainage_land0 [%d,%d,%e]: %e,%e, %d,%d,%d\n",
    //               patch[0].ID, kk, rate_,
    //               patch[0].sat_deficit_z, patch[0].sat_deficit,
    //               sat_def_z_index,patch[0].sat_def_pct_index,patch[0].soil_defaults[0][0].soildepthLen);
    //    }//debug
        
    
    /*--------------------------------------------------------------*/
    /*    available sat water  */
    /*--------------------------------------------------------------*/
    std_scale = command_line[0].std_scale;
    available_sat_water = max(
                patch[0].area*( patch[0].soil_defaults[0][0].soil_water_cap - max(patch[0].sat_deficit, 0.0) ),
                0.0);
        
    //    if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity || patch[0].sat_def_pct_index!=patch[0].sat_def_pct_index || patch[0].sat_def_pct_indexM!=patch[0].sat_def_pct_indexM){
    //        printf("drainage_land(1,%d): %lf %lf %lf %lf (%lf, %lf)\n", patch[0].ID,
    //               patch[0].sat_deficit, patch[0].sat_deficit_z,
    //               patch[0].rootzone.field_capacity, patch[0].field_capacity,
    //               patch[0].sat_def_pct_index, patch[0].sat_def_pct_indexM);
    //    }//debug
    
    
    
	/*--------------------------------------------------------------*/
	/*	calculate drain-to fluxes (i.e., drain out from current patch)
    /*  move gamma caluclation (function call) to here, Feb 11, 2020, Lin  */
    /*  move varbased_flow caluclation (function call) to here, Feb 11, 2020, Lin  */
	/*--------------------------------------------------------------*/
    int i; //j and d are used
    double adjustment = 0.0;
    double revised_total_gamma = 0.0;
    double totaledge = 0.0;
    double z1, z2, water_table_z1, water_table_z2;
    double soil_depth_z1, soil_depth_z2;
    int z1_sat_def_pct_index_z2;
    double z1_sat_def_pct_indexM_z2;
    double z1_soildepth_z2;
    double sat_to_gw_coeff = patch[0].soil_defaults[0][0].sat_to_gw_coeff;
    
    z1 = patch[0].z;
    if(patch[0].sat_deficit_z > 0){
        water_table_z1 = (z1 - patch[0].sat_deficit_z);
    }else{
        water_table_z1 = z1;
    }
    soil_depth_z1 = z1 - patch[0].soil_defaults[0][0].soil_depth; // elevation of soil bottom
//    if(patch[0].sat_deficit_z>patch[0].soil_defaults[0][0].soil_depth) printf("warning0 patch %d(%f,%f)->(%f,%f)\n", patch[0].ID, patch[0].sat_deficit_z, patch[0].soil_defaults[0][0].soil_depth, water_table_z1, soil_depth_z1);

    d = 0;
    if (patch[0].innundation_list[d].num_neighbours > 0){
        
        // have neighbor
        for (i =0; i < patch[0].innundation_list[d].num_neighbours; i++) {
            
            z2 = patch[0].innundation_list[d].neighbours[i].patch[0].z;
            if(patch[0].innundation_list[d].neighbours[i].patch[0].sat_deficit_z > 0){
                water_table_z2 = (z2 - patch[0].innundation_list[d].neighbours[i].patch[0].sat_deficit_z);
            }else{
                water_table_z2 = z2;
            }
            
            // new code, Sept 19, 2018
            // - patch[0].innundation_list[d].neighbours[i].patch[0].constraintWaterTableTopDepth
            if( (water_table_z1 - water_table_z2 )>0 ){
                // updating the gamma fractions: step 1 comparing water table depths
                patch[0].innundation_list[d].neighbours[i].gamma = (water_table_z1 - water_table_z2) * patch[0].innundation_list[d].neighbours[i].edgedistance; // (edge/distance)
                revised_total_gamma += patch[0].innundation_list[d].neighbours[i].gamma;
                totaledge += patch[0].innundation_list[d].neighbours[i].edge;
            }else{
                patch[0].innundation_list[d].neighbours[i].gamma = 0.0;
            }// end of if
            
            // how to handle cliff effect?
            
            
        }//end of for neighbour i loop
        
        
        if(revised_total_gamma>0){
            
            for (i =0; i < patch[0].innundation_list[d].num_neighbours; i++) {
                
                z2 = patch[0].innundation_list[d].neighbours[i].patch[0].z;
                soil_depth_z2 = z2 - patch[0].innundation_list[d].neighbours[i].patch[0].soil_defaults[0][0].soil_depth;
                if(patch[0].innundation_list[d].neighbours[i].patch[0].sat_deficit_z > 0){
                    water_table_z2 = (z2 - patch[0].innundation_list[d].neighbours[i].patch[0].sat_deficit_z);
                }else{
                    water_table_z2 = z2;
                }
                // new code, Feb 11, 2020
                if( patch[0].innundation_list[d].neighbours[i].gamma>0.0 ){
                    
                    // updating the gamma fractions: step 2 elevation of soil bottoms
                    patch[0].innundation_list[d].neighbours[i].transmissivity_flux2neighbour = time_int * compute_varbased_flow(
                        patch[0].num_soil_intervals,
                        patch[0].std * std_scale,
                        patch[0].sat_def_pct_index,
                        patch[0].sat_def_pct_indexM,
                        patch[0].innundation_list[d].neighbours[i].gamma*patch[0].area/totaledge,
                        patch);
//                    if(patch[0].innundation_list[d].neighbours[i].transmissivity_flux2neighbour<0.0) printf("warning1 negative routing from patch %d(%f,%f) -> %d(%f,%f) with gamma (fraction) %lf\n",patch[0].ID, water_table_z1, soil_depth_z1, patch[0].innundation_list[d].neighbours[i].patch[0].ID, water_table_z2, soil_depth_z2, patch[0].innundation_list[d].neighbours[i].transmissivity_flux2neighbour);
                    
                    
                    
                    // count for soil depth blocking
                    if( soil_depth_z2 > soil_depth_z1){
                        z1_soildepth_z2 = patch[0].soil_defaults[0][0].soil_depth-(soil_depth_z2-soil_depth_z1);
                        z1_sat_def_pct_index_z2 = patch[0].soil_defaults[0][0].rtz2sat_def_pct_index[(int)(z1_soildepth_z2*1000)];
                        z1_sat_def_pct_indexM_z2 = 1000*(patch[0].soil_defaults[0][0].rtz2sat_def_0z[(int)(z1_soildepth_z2*1000)]*patch[0].soil_defaults[0][0].max_sat_def_1 - z1_sat_def_pct_index_z2*0.001);
                        
//                        if(z1_sat_def_pct_indexM_z2<0.0 || z1_sat_def_pct_indexM_z2>1.0) printf("warning1.5 patch %d %f<%f %f %f\n", patch[0].ID, z1_soildepth_z2, patch[0].soil_defaults[0][0].soil_depth, patch[0].soil_defaults[0][0].rtz2sat_def_0z[(int)(z1_soildepth_z2*1000)], z1_sat_def_pct_index_z2*0.001);
                        
                        patch[0].innundation_list[d].neighbours[i].transmissivity_flux2neighbour -= time_int * compute_varbased_flow(
                            patch[0].num_soil_intervals,
                            patch[0].std * std_scale,
                            z1_sat_def_pct_index_z2,//<---- the function is reading "patch[0].sat_def_pct_index"
                            z1_sat_def_pct_indexM_z2,
                            patch[0].innundation_list[d].neighbours[i].gamma*patch[0].area/totaledge,
                            patch);
                        
                        if(patch[0].innundation_list[d].neighbours[i].transmissivity_flux2neighbour<0.0 && patch[0].innundation_list[d].neighbours[i].transmissivity_flux2neighbour > -0.0001){
                            
                            patch[0].innundation_list[d].neighbours[i].transmissivity_flux2neighbour = 0.0;
                        }else if(patch[0].innundation_list[d].neighbours[i].transmissivity_flux2neighbour<0.0) printf("warning2 negative routing from patch %d(%f,%f) -> %d(%f,%f) with gamma (fraction) %lf, %d(%f,%f), %d(%f,%f)\n",patch[0].ID, water_table_z1, soil_depth_z1, patch[0].innundation_list[d].neighbours[i].patch[0].ID, water_table_z2, soil_depth_z2, patch[0].innundation_list[d].neighbours[i].transmissivity_flux2neighbour, patch[0].sat_def_pct_index, patch[0].sat_def_pct_indexM, patch[0].sat_deficit_z, z1_sat_def_pct_index_z2, z1_sat_def_pct_indexM_z2, z1_soildepth_z2 );
                    }//if
                    
                }else{
                    patch[0].innundation_list[d].neighbours[i].transmissivity_flux2neighbour = 0.0;
                }// end of if
                
                route_to_patch += patch[0].innundation_list[d].neighbours[i].transmissivity_flux2neighbour;
                
            }//end of for neighbour i loop
            
            if(route_to_patch!=route_to_patch){
                printf("drainage_land(2,%d): %lf %lf %lf %lf (%lf,%d, %lf)\n", patch[0].ID,
                       patch[0].sat_deficit, patch[0].sat_deficit_z,
                       patch[0].rootzone.field_capacity, patch[0].field_capacity,
                       patch[0].sat_def_pct, patch[0].sat_def_pct_index, patch[0].sat_def_pct_indexM);
            }//debug
            if (route_to_patch <=0.0){
                route_to_patch = 0.0;
                for (i =0; i < patch[0].innundation_list[d].num_neighbours; i++) {
                    patch[0].innundation_list[d].neighbours[i].gamma = 0.0;
                }//end of for neighbour i loop
                revised_total_gamma = 0.0;
            }else{
                
                for (i =0; i < patch[0].innundation_list[d].num_neighbours; i++) {
                    patch[0].innundation_list[d].neighbours[i].gamma = patch[0].innundation_list[d].neighbours[i].transmissivity_flux2neighbour/route_to_patch;
                    // neighbours' gamma (fractions) are still in use in the code below
                }//end of for neighbour i loop
            }// end of if else
            
            // ---- just keep the total gamma here, but no use in further down
            revised_total_gamma /= totaledge;
            revised_total_gamma *= patch[0].area;
        }else{
            revised_total_gamma = 0.0;
            route_to_patch = 0.0;
            //default: patch[0].innundation_list[d].neighbours[i].gamma = 0;
        }//
        total_gamma = revised_total_gamma;
        
    }else{
        // no neighbor; num_neighbours = 0
        // adjustment = 1.0;
        total_gamma = patch[0].innundation_list[d].gamma;
        route_to_patch = time_int * compute_varbased_flow(
            patch[0].num_soil_intervals,
            patch[0].std * std_scale,
            patch[0].sat_def_pct_index,
            patch[0].sat_def_pct_indexM,
            total_gamma,
            patch);
    }// end of if

    if (route_to_patch < 0.0) route_to_patch = 0.0;
    if (route_to_patch > 0.0 && route_to_patch > available_sat_water) route_to_patch=available_sat_water; //volumn
    
    /*--------------------------------------------------------------*/
    /*  calculate drain-in fluxes (i.e., drain into current patch from remote patches through "irrigration")
    /*  this is different from pipe / sewer drainage that any "left-over" water will be drain from the current patch, */
    /*  and there are no concerns of how much destination patches will receive the water */
    /*  "irrigration" is a control-transport of water, and its controls are on the destination patches, not the source ones. */
    /*  For example, when a lawn does need irrigration, then there is no transfer of the water from the source even if there */
    /*  is "left-over" water.  */
    /*  the workflow of this drain-in process starts from patch_dailyF(). */
    /*  in patch_dailyF(), the drainIN transfer flux is calculated and distributed to the source patches */
    /*  note that "maxDailyDrain" and "DrainFrac" in the drainIN object control the how much water can be withdraw from source */
    /*--------------------------------------------------------------*/
//    // 1) calculating transfers; this is calculated in patch_dailyF()
//    for (i =0; i < patch[0].innundation_list[d].num_drainIN; i++){
//        patch[0].innundation_list[d].drainIN[i].transfer_flux;
//    }// end of for loop i
    
//    // 2) adding transfers to the patch; this is calculated in patch_dailyF()
//    patch[0].grassIrrigation_m; // added in patch_dailyF(), patch[0].detention_store @LINE 1687
//    patch[0].landuse_defaults[0][0].septic_water_load // added in patch_dailyF() @LINE 920
//    patch[0].landuse_defaults[0][0].septic_NO3_load // added in patch_dailyF() @LINE 920
    
    
    //Sept 18, 2019
    extrawater = patch[0].rz_storage+patch[0].unsat_storage - patch[0].sat_deficit - route_to_patch/patch[0].area;// + patch[0].constraintWaterTableTopDepth_def;
    if(extrawater>0){
        // bounded by impervious surface
        // trigger returnflow and fully saturation
        double imp_unsat_storage = (patch[0].rz_storage+patch[0].unsat_storage)*patch[0].Ksat_vertical;
        double imp_sat_def = (patch[0].sat_deficit+route_to_patch/patch[0].area)*patch[0].Ksat_vertical;
        return_flow = max( (imp_unsat_storage-imp_sat_def),0.0); // no counting litter storage
        route_to_patch += (extrawater>return_flow? (extrawater-return_flow)*patch[0].area : 0.0);  // volumn
        
        // counting litter storage for returnflow
        return_flow = compute_varbased_returnflow(
            patch[0].std * std_scale,
            imp_unsat_storage,
            imp_sat_def,
            &(patch[0].litter));// mm
        
        patch[0].detention_store += return_flow;
    }// extra water
    
	/*--------------------------------------------------------------*/
	/* compute Nitrogen leaching amount		(subsurface)		*/
	/*--------------------------------------------------------------*/
	if(command_line[0].grow_flag > 0) {
        
        Nout = compute_N_leached(
             verbose_flag,
             patch[0].sat_NO3, //total solute <--- projected total_sat_solute
             route_to_patch / patch[0].area,//Qout (already time_int adjusted)
             patch[0].soil_defaults[0][0].NO3decayRate,
             patch[0].soil_defaults[0][0].active_zone_z, //activedepthz
             patch[0].soil_defaults[0][0].NO3_adsorption_rate, // N_absorption_rate
             5,patch); // signal, patch
		NO3_leached_to_patch = Nout * patch[0].area;
		patch[0].soil_ns.NO3_Qout += Nout;//<<-----------*********
        if(Nout<0 || Nout!=Nout || (Nout<=0&&patch[0].sat_NO3>ZERO&&route_to_patch>ZERO) ||
           patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
           patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
           patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
           patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
           patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
           patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
           patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
           patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
           patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0 ) printf("update_drainage_land[%d,%d, %e,%e]:[%e %e %e] [%e %e %e] [%e %e %e]\n", patch[0].ID, kk, route_to_patch, Nout,
               patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
               patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
               patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);


		Nout = compute_N_leached(
			verbose_flag,
			patch[0].sat_NH4,
			route_to_patch / patch[0].area,
			patch[0].soil_defaults[0][0].NH4decayRate,
            patch[0].soil_defaults[0][0].active_zone_z,
			patch[0].soil_defaults[0][0].NH4_adsorption_rate,
            8,patch);
		NH4_leached_to_patch = Nout * patch[0].area;
		patch[0].soil_ns.NH4_Qout += Nout;
        if(Nout<0 || Nout!=Nout ) printf("update_drainage_land[%d,%e]: soil NH4 (%e) flux (%e)\n", patch[0].ID,route_to_patch,patch[0].soil_ns.sminn, Nout);

		Nout = compute_N_leached(
			verbose_flag,
            patch[0].sat_DON,
			route_to_patch / patch[0].area,
			patch[0].soil_defaults[0][0].DOMdecayRate,
			patch[0].soil_defaults[0][0].active_zone_z,
			patch[0].soil_defaults[0][0].DON_adsorption_rate,
			11,patch);
		DON_leached_to_patch = Nout * patch[0].area;
		patch[0].soil_ns.DON_Qout += Nout;
        if(Nout<0 || Nout!=Nout) printf("update_drainage_land[%d,%e]: soil DON (%e) flux (%e)\n", patch[0].ID,route_to_patch,patch[0].soil_ns.DON, Nout);

		Nout = compute_N_leached(
			verbose_flag,
			patch[0].sat_DOC,
			route_to_patch / patch[0].area,
			patch[0].soil_defaults[0][0].DOMdecayRate,
			patch[0].soil_defaults[0][0].active_zone_z,
			patch[0].soil_defaults[0][0].DOC_adsorption_rate,
			14,patch);
		DOC_leached_to_patch = Nout * patch[0].area;
		patch[0].soil_cs.DOC_Qout += Nout;
        if(Nout<0 || Nout!=Nout) printf("update_drainage_land[%d,%e]: soil DOC (%e) flux (%e)\n", patch[0].ID,route_to_patch,patch[0].soil_cs.DOC, Nout);


	}//leaching from subsurface growth flag

//	if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity || patch[0].sat_def_pct_index!=patch[0].sat_def_pct_index || patch[0].sat_def_pct_indexM!=patch[0].sat_def_pct_indexM){
//        printf("drainage_land(3,%d): %lf %lf %lf %lf (%lf, %lf)\n", patch[0].ID,
//               patch[0].sat_deficit, patch[0].sat_deficit_z,
//               patch[0].rootzone.field_capacity, patch[0].field_capacity,
//               patch[0].sat_def_pct_index, patch[0].sat_def_pct_indexM);
//    }//debug
    
	patch[0].Qout += (route_to_patch / patch[0].area);
    if(extrawater>0){
        // @line 260 extrawater = patch[0].rz_storage+patch[0].unsat_storage - patch[0].sat_deficit - route_to_patch/patch[0].area + patch[0].constraintWaterTableTopDepth_def;
        // @line 271 route_to_patch += (1-extrawater)*patch[0].Ksat_vertical
        //
        // set this to be  "-route_to_patch/patch[0].area" and later += Qout = "route_to_patch/patch[0].area"
        //patch[0].sat_deficit = -route_to_patch/patch[0].area + patch[0].constraintWaterTableTopDepth_def; // -- problem: what of patch[0].sat_deficit<0 from start?
        // extrawater*patch[0].Ksat_vertical = amount of water travels vertically up to surface; while the rest of it goes subsurface
        //
        patch[0].sat_deficit -= patch[0].unsat_storage + patch[0].rz_storage - extrawater*patch[0].Ksat_vertical; // i don't know the answer
        patch[0].unsat_storage = 0.0; // converted to be part of sat
        patch[0].rz_storage = 0.0; // converted to be part of sat
    }// extra water

//    if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity || patch[0].sat_def_pct_index!=patch[0].sat_def_pct_index || patch[0].sat_def_pct_indexM!=patch[0].sat_def_pct_indexM){
//        printf("drainage_land(4,%d): %lf %lf %lf %lf (%lf, %lf)\n", patch[0].ID,
//               patch[0].sat_deficit, patch[0].sat_deficit_z,
//               patch[0].rootzone.field_capacity, patch[0].field_capacity,
//               patch[0].sat_def_pct_index, patch[0].sat_def_pct_indexM);
//    }//debug
    
   
    // below is the returnflow associated solute processes
	/*--------------------------------------------------------------*/
	/*	calculated any N-transport associated with return flow  */
	/*	-note available N reduced by what has already been 	*/
	/*	we assume that only nitrate follows return flow		*/
	/*	lost in subsurface flow routing				*/
	/*--------------------------------------------------------------*/
    // leaching from surface;
    // when return_flow>0, then extrawater = patch[0].rz_storage+patch[0].unsat_storage - patch[0].sat_deficit - route_to_patch/patch[0].area + patch[0].constraintWaterTableTopDepth_def >0
		if (command_line[0].grow_flag > 0 && return_flow>0) {
			Nout = compute_N_leached(
				verbose_flag,
                /// problem re-project!! && "soil_ns -= going2sat_NO3"
				patch[0].sat_NO3 - patch[0].soil_ns.NO3_Qout,
				return_flow,
				patch[0].soil_defaults[0][0].NO3decayRate,
				patch[0].soil_defaults[0][0].active_zone_z,
				patch[0].soil_defaults[0][0].NO3_adsorption_rate,
                17,patch);
			patch[0].surface_NO3 += Nout;
			patch[0].soil_ns.NO3_Qout += Nout;
            if(Nout<0 || Nout!=Nout || patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 || NO3_leached_to_patch<0 || NO3_leached_to_patch!=NO3_leached_to_patch) printf("update_drainage_land[%d,%e]: return NO3 (%e,%e) flux (%e)\n", patch[0].ID,return_flow,patch[0].soil_ns.nitrate - (NO3_leached_to_patch/patch[0].area),NO3_leached_to_patch, Nout);

			Nout = compute_N_leached(
				verbose_flag,
				patch[0].sat_NH4 - patch[0].soil_ns.NH4_Qout,
				return_flow,
				patch[0].soil_defaults[0][0].NH4decayRate,
                patch[0].soil_defaults[0][0].active_zone_z,
				patch[0].soil_defaults[0][0].NH4_adsorption_rate,
				20,patch);
			patch[0].surface_NH4 += Nout;
			patch[0].soil_ns.NH4_Qout += Nout;
            if(Nout<0 || Nout!=Nout) printf("update_drainage_land[%d,%e]: return NH4 (%e) flux (%e)\n", patch[0].ID,return_flow,patch[0].soil_ns.sminn - (NH4_leached_to_patch/patch[0].area), Nout);

			Nout = compute_N_leached(
				verbose_flag,
				patch[0].sat_DON - patch[0].soil_ns.DON_Qout,
				return_flow,
				patch[0].soil_defaults[0][0].DOMdecayRate,
				patch[0].soil_defaults[0][0].active_zone_z,
				patch[0].soil_defaults[0][0].DON_adsorption_rate,
				23,patch);
			patch[0].surface_DON += Nout;
			patch[0].soil_ns.DON_Qout += Nout;
            if(Nout<0 || Nout!=Nout) printf("update_drainage_land[%d,%e]: return DON (%e) flux (%e)\n", patch[0].ID,return_flow,patch[0].soil_ns.DON - (DON_leached_to_patch/patch[0].area), Nout);
            
			Nout = compute_N_leached(
				verbose_flag,
				patch[0].sat_DOC - patch[0].soil_cs.DOC_Qout,
				return_flow,
				patch[0].soil_defaults[0][0].DOMdecayRate,
				patch[0].soil_defaults[0][0].active_zone_z,
				patch[0].soil_defaults[0][0].DOC_adsorption_rate,
				26,patch);
			patch[0].surface_DOC += Nout;
			patch[0].soil_cs.DOC_Qout += Nout;
            if(Nout<0 || Nout!=Nout) printf("update_drainage_land[%d,%e]: return DOC (%e) flux (%e)\n", patch[0].ID,return_flow,patch[0].soil_cs.DOC - (DOC_leached_to_patch/patch[0].area), Nout);
		}//
    
    patch[0].overland_flow += return_flow;
    
	/*--------------------------------------------------------------*/
	/*	route water and nitrogen lossed due to infiltration excess */
	/*--------------------------------------------------------------*/
	if ( (patch[0].detention_store > patch[0].landuse_defaults[0][0].detention_store_size* (1.0 - patch[0].Ksat_vertical)+patch[0].landuse_defaults[0][0].pond_size * patch[0].waterFrac) && (patch[0].detention_store > ZERO) ){

		Qout = (patch[0].detention_store - patch[0].landuse_defaults[0][0].detention_store_size* (1.0 - patch[0].Ksat_vertical)-patch[0].landuse_defaults[0][0].pond_size * patch[0].waterFrac);
		if (command_line[0].grow_flag > 0) {
                Nout = (min(1.0, (Qout/ patch[0].detention_store))) * patch[0].surface_DOC;
                DOC_leached_to_surface = Nout * patch[0].area;
                patch[0].surface_DOC -= Nout;
            
                Nout = (min(1.0, (Qout/ patch[0].detention_store))) * patch[0].surface_DON;
                DON_leached_to_surface = Nout * patch[0].area;
                patch[0].surface_DON -= Nout;
            
                Nout = (min(1.0, (Qout/ patch[0].detention_store))) * patch[0].surface_NO3;
                NO3_leached_to_surface = Nout * patch[0].area;
                patch[0].surface_NO3 -= Nout;
            
                Nout = (min(1.0, (Qout/ patch[0].detention_store))) * patch[0].surface_NH4;
                NH4_leached_to_surface = Nout * patch[0].area;
                patch[0].surface_NH4 -= Nout;
        }// end of grow_flag if statement
        
        route_to_surface = (Qout *  patch[0].area);
        patch[0].detention_store -= Qout;
        patch[0].surface_Qout += Qout;

    }// end of if statement
			

	if (NO3_leached_to_surface < 0.0)
		printf("WARNING %d %lf",patch[0].ID, NO3_leached_to_surface);

    
    
    
    
    
    
    
    
	/*--------------------------------------------------------------*/
	/*	route flow to neighbours				*/
	/*	route n_leaching if grow flag specfied			*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/* regular downslope routing */
	/*--------------------------------------------------------------*/
	if (command_line[0].noredist_flag == 0) {
        d=0;
        for (j = 0; j < patch[0].innundation_list[d].num_neighbours; j++) {
            neigh = patch[0].innundation_list[d].neighbours[j].patch;
            /*--------------------------------------------------------------*/
            /* first transfer subsurface water and nitrogen */  // --------- subsurface
            /*--------------------------------------------------------------*/
            /* some "transmissivity_flux2neighbour" is loss to GW_storage  */
            if(command_line[0].gw_flag > 0 && (patch[0].drainage_type>0 && patch[0].drainage_type % actionGWDRAIN==0)){
                //how do we know how fill is the GW?
                patch[0].gw_drainage += route_to_patch * sat_to_gw_coeff;// has multiplied patch[0].area // reset to zero every patch_daily_I()
                route_to_patch *= 1.0 - sat_to_gw_coeff;
            }//end of if
            Qin = (patch[0].innundation_list[d].neighbours[j].gamma * route_to_patch) / neigh[0].area;
            if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){neigh[0].fromLAND_Q+=Qin; }//reset @subsurface_routing
            if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){neigh[0].fromRIPARIAN_Q+=Qin; }
            if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){neigh[0].fromSTREAM_Q+=Qin; }
            
            
            if (Qin < 0) printf("warning negative routing from patch %d with gamma %lf\n", patch[0].ID, total_gamma);
            if (command_line[0].grow_flag > 0) {
                
                double coef;
                if(command_line[0].gw_flag > 0 && (patch[0].drainage_type>0 && patch[0].drainage_type % actionGWDRAIN==0)){
                    coef = 1.0 - sat_to_gw_coeff;
                    coef /= neigh[0].area;
                    patch[0].gw_drainage_DON += DON_leached_to_patch * sat_to_gw_coeff;// has multiplied patch[0].area
                    patch[0].gw_drainage_DOC += DOC_leached_to_patch * sat_to_gw_coeff;// has multiplied patch[0].area
                    patch[0].gw_drainage_NO3 += NO3_leached_to_patch * sat_to_gw_coeff;// has multiplied patch[0].area
                    patch[0].gw_drainage_NH4 += NH4_leached_to_patch * sat_to_gw_coeff;// has multiplied patch[0].area
                }else{
                    coef = 1.0/neigh[0].area;
                }

                //Note that: "xxx_leached_to_patch" is flux, i.e., not areal
                //Note that: "Nin" is an areal adjustment to the neigh due to influx
                Nin = (patch[0].innundation_list[d].neighbours[j].gamma * DON_leached_to_patch) * coef;
                neigh[0].soil_ns.DON_Qin += Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){neigh[0].fromLAND_DON+=Nin; }//reset @subsurface_routing
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){neigh[0].fromRIPARIAN_DON+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){neigh[0].fromSTREAM_DON+=Nin; }
                
                Nin = (patch[0].innundation_list[d].neighbours[j].gamma * DOC_leached_to_patch) * coef;
                neigh[0].soil_cs.DOC_Qin += Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){neigh[0].fromLAND_DOC+=Nin; }//reset @subsurface_routing
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){neigh[0].fromRIPARIAN_DOC+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){neigh[0].fromSTREAM_DOC+=Nin; }
                
                Nin = (patch[0].innundation_list[d].neighbours[j].gamma * NO3_leached_to_patch) * coef;
                neigh[0].soil_ns.NO3_Qin += Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){neigh[0].fromLAND_NO3+=Nin; }//reset @subsurface_routing
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){neigh[0].fromRIPARIAN_NO3+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){neigh[0].fromSTREAM_NO3+=Nin; }
                
                Nin = (patch[0].innundation_list[d].neighbours[j].gamma * NH4_leached_to_patch) * coef;
                neigh[0].soil_ns.NH4_Qin += Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){neigh[0].fromLAND_NH4+=Nin; }//reset @subsurface_routing
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){neigh[0].fromRIPARIAN_NH4+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){neigh[0].fromSTREAM_NH4+=Nin; }
                
                }//if
            neigh[0].Qin += Qin; // subsurface --> neighbour sat_def
        }// end of for subsurface routing loop

        
        
        /*--------------------------------------------------------------*/
        /* surface downslope routing */
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        /* determine which innundation depth to consider		*/
        /*--------------------------------------------------------------*/
        if (patch[0].num_innundation_depths > 0) {
            innundation_depth = patch[0].detention_store + route_to_surface/patch[0].area;
            d=0;
            while ((innundation_depth > patch[0].innundation_list[d].critical_depth) && (d < patch[0].num_innundation_depths-1)) {
                d++;}// while
            }// if
        else d=0;
        
        for (j = 0; j < patch[0].surface_innundation_list[d].num_neighbours; j++) {

            neigh = patch[0].surface_innundation_list[d].neighbours[j].patch;

            /*--------------------------------------------------------------*/
            /* now transfer surface water and nitrogen */ // -------------- surface (updated Spet 12)
            /*	- first nitrogen					*/
            /*--------------------------------------------------------------*/
            // route_to_surface = patch[0].area * (patch[0].detention_store - patch[0].landuse_defaults[0][0].detention_store_size);
            if (command_line[0].grow_flag > 0) {
                if(neigh[0].drainage_type==STREAM){
                    
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * NO3_leached_to_surface) / neigh[0].area;
                    neigh[0].streamflow_NO3 += Nin;
                    if(neigh[0].ID==command_line[0].outletPatchID && patch[0].drainage_type!=STREAM){neigh[0].stormdrained_NO3 += Nin;}
                    
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * NH4_leached_to_surface) / neigh[0].area;
                    neigh[0].streamflow_NH4 += Nin;
                    if(neigh[0].ID==command_line[0].outletPatchID && patch[0].drainage_type!=STREAM){neigh[0].stormdrained_NH4 += Nin;}
                    
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * DON_leached_to_surface) / neigh[0].area;
                    neigh[0].streamflow_DON += Nin;
                    if(neigh[0].ID==command_line[0].outletPatchID && patch[0].drainage_type!=STREAM){neigh[0].stormdrained_DON += Nin;}
                    
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * DOC_leached_to_surface) / neigh[0].area;
                    neigh[0].streamflow_DOC += Nin;
                    if(neigh[0].ID==command_line[0].outletPatchID && patch[0].drainage_type!=STREAM){neigh[0].stormdrained_DOC += Nin;}
                }else{
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * NO3_leached_to_surface) / neigh[0].area;
                    neigh[0].surface_NO3 += Nin;
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * NH4_leached_to_surface) / neigh[0].area;
                    neigh[0].surface_NH4 += Nin;
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * DON_leached_to_surface) / neigh[0].area;
                    neigh[0].surface_DON += Nin;
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * DOC_leached_to_surface) / neigh[0].area;
                    neigh[0].surface_DOC += Nin;
                }
            }//if growth
            
            /*--------------------------------------------------------------*/
            /*	- now surface water 					*/
            /*	surface stores should be updated to facilitate transfer */
            /* added net surface water transfer to detention store		*/
            /*--------------------------------------------------------------*/

            // re-work here!!
            Qin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * route_to_surface) / neigh[0].area;
            if( neigh[0].drainage_type==STREAM && neigh[0].ID==command_line[0].outletPatchID && patch[0].drainage_type>0 && patch[0].drainage_type % actionSTORMDRAIN==0 ){
                
                //patch[0].drainage_type!=STREAM
                patch[0].stormdrainYield += Qin; //Spet 17 tracking how much is storm drain yielded at "local patch"
                neigh[0].stormdrained += Qin; //Spet 15 tracking how much storm drain to the "outlet" (basin)
                
                neigh[0].surface_Qin += Qin;
                neigh[0].streamflow += Qin; // short-cut by passing the detention and infilration at stream patch; but trying to seperate return and stormdrain
                // perviously was if neigh[0].drainage_type==STREAM, then all surface water becomes streamflow
                
            }else{
                neigh[0].detention_store += Qin;// need fix this ****
                neigh[0].surface_Qin += Qin;
            }
            
            /*--------------------------------------------------------------*/
            /* try to infiltrate this water					*/
            /* use time_int as duration */
            /*--------------------------------------------------------------*/
            if (neigh[0].detention_store > ZERO) {
                
                infiltration = compute_infiltration(
                    command_line[0].verbose_flag,
                    neigh[0].sat_deficit_z,
                    0.0, //neigh[0].aboveWT_SatPct, // initiated in daily_I()
                    neigh[0].Ksat_vertical, // 1- impervious
                    neigh[0].sat_def_pct_indexM * neigh[0].soil_defaults[0][0].vksat_0zm[neigh[0].sat_def_pct_index+1] + (1.0-neigh[0].sat_def_pct_indexM) * neigh[0].soil_defaults[0][0].vksat_0zm[neigh[0].sat_def_pct_index],
                    neigh[0].rz_storage+neigh[0].unsat_storage,
                    neigh[0].sat_def_pct_indexM * neigh[0].soil_defaults[0][0].sat_def_0zm[neigh[0].sat_def_pct_index+1] + (1.0-neigh[0].sat_def_pct_indexM) * neigh[0].soil_defaults[0][0].sat_def_0zm[neigh[0].sat_def_pct_index],
                    neigh[0].sat_deficit,
                    neigh[0].detention_store,
                    time_int,
                    neigh[0].soil_defaults[0][0].psi_air_entry);
                

            } else infiltration = 0.0;
            /*--------------------------------------------------------------*/
            /* added an surface N flux to surface N pool	and		*/
            /* allow infiltration of surface N				*/
            /*--------------------------------------------------------------*/
            if ((command_line[0].grow_flag > 0 ) && (infiltration > ZERO)) {
                Nin = ((infiltration / neigh[0].detention_store) * neigh[0].surface_DOC);
                neigh[0].soil_cs.DOC_Qin += Nin;
                neigh[0].surface_DOC -= Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){patch[0].fromLAND_surfsubDOC+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){patch[0].fromRIPARIAN_surfsubDOC+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){patch[0].fromSTREAM_surfsubDOC+=Nin; }
                
                Nin = ((infiltration / neigh[0].detention_store) * neigh[0].surface_DON);
                neigh[0].soil_ns.DON_Qin += Nin;
                neigh[0].surface_DON -= Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){patch[0].fromLAND_surfsubDON+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){patch[0].fromRIPARIAN_surfsubDON+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){patch[0].fromSTREAM_surfsubDON+=Nin; }
                
                Nin = ((infiltration / neigh[0].detention_store) * neigh[0].surface_NO3);
                neigh[0].soil_ns.NO3_Qin += Nin;
                neigh[0].surface_NO3 -= Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){patch[0].fromLAND_surfsubNO3+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){patch[0].fromRIPARIAN_surfsubNO3+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){patch[0].fromSTREAM_surfsubNO3+=Nin; }
                
                Nin = ((infiltration / neigh[0].detention_store) * neigh[0].surface_NH4);
                neigh[0].soil_ns.NH4_Qin += Nin;
                neigh[0].surface_NH4 -= Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){patch[0].fromLAND_surfsubNH4+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){patch[0].fromRIPARIAN_surfsubNH4+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){patch[0].fromSTREAM_surfsubNH4+=Nin; }
            }

            if (infiltration > neigh[0].sat_deficit - neigh[0].unsat_storage - neigh[0].rz_storage) {
                neigh[0].sat_deficit -= (infiltration + neigh[0].unsat_storage + neigh[0].rz_storage);
                neigh[0].unsat_storage = 0.0;
                neigh[0].rz_storage = 0.0;
                neigh[0].field_capacity = 0.0;
                neigh[0].rootzone.field_capacity = 0.0;
            }

            else if ((neigh[0].sat_deficit > neigh[0].rootzone.potential_sat) &&
                (infiltration > neigh[0].rootzone.potential_sat - neigh[0].rz_storage)) {
            /*------------------------------------------------------------------------------*/
            /*		Just add the infiltration to the rz_storage and unsat_storage	*/
            /*------------------------------------------------------------------------------*/
                neigh[0].unsat_storage += infiltration - (neigh[0].rootzone.potential_sat - neigh[0].rz_storage);
                neigh[0].rz_storage = neigh[0].rootzone.potential_sat;
            }
            /* Only rootzone layer saturated - perched water table case */
            else if ((neigh[0].sat_deficit > neigh[0].rootzone.potential_sat) &&
                (infiltration <= neigh[0].rootzone.potential_sat - neigh[0].rz_storage)) {
                /*--------------------------------------------------------------*/
                /*		Just add the infiltration to the rz_storage	*/
                /*--------------------------------------------------------------*/
                neigh[0].rz_storage += infiltration;
            }
            else if ((neigh[0].sat_deficit <= neigh[0].rootzone.potential_sat) &&
                (infiltration <= neigh[0].sat_deficit - neigh[0].rz_storage - neigh[0].unsat_storage)) {
                neigh[0].rz_storage += neigh[0].unsat_storage;
                /* transfer left water in unsat storage to rootzone layer */
                neigh[0].unsat_storage = 0;
                neigh[0].rz_storage += infiltration;
                neigh[0].field_capacity = 0;
            }

            neigh[0].detention_store -= infiltration;

        }// surface routing loop
        
        

        
        
        

	} /* end if noredistribution flag */

	return;

} /*end update_drainage_land.c*/

