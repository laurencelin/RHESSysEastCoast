/*--------------------------------------------------------------*/
/* 											*/
/*					compute_subsurface_routing			*/
/*											*/
/*	compute_subsurface_routing.c - creates a patch object				*/
/*											*/
/*	NAME										*/
/*	compute_subsurface_routing.c - creates a patch object				*/
/*											*/
/*	SYNOPSIS									*/
/*	struct routing_list_object compute_subsurface_routing( 				*/
/*							struct command_line_object command */
/*							struct basin_object *basinn)	*/
/*				 			int,			 	*/
/*							struct date *current_date)	*/
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
/*	June 16, 98 C.Tague								*/
/*	limit drainage to maximum saturation deficit defined by soil depth		*/
/*											*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

void compute_subsurface_routing(struct command_line_object *command_line,
		struct basin_object *basin, int n_timesteps, struct date current_date) {
	/*--------------------------------------------------------------*/
	/*	Local function definition.				*/
	/*--------------------------------------------------------------*/

	void update_drainage_stream(struct patch_object *,
			struct command_line_object *, double, int);

	void update_drainage_road(struct patch_object *,
			struct command_line_object *, double, int);

	void update_drainage_land(struct patch_object *,
			struct command_line_object *, double, int);

	double compute_infiltration(int, double, double, double, double, double,
			double, double, double, double, double);

	double compute_z_final(int, double, double, double, double, double);

	double compute_N_leached(
             int verbose_flag,
             double total_nitrate,
             double Qout,
             double N_decay_rate,
             double activedepthz,
             double N_absorption_rate,
             int signal,
             struct patch_object *patch);

	double compute_layer_field_capacity(int, int, double, double, double,
			double, double, double, double, double, double);

	double compute_unsat_zone_drainage(int, int, double, double, double, double,
			double, double, double);

	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
	int i, d;
	int j, k;
	int grow_flag, verbose_flag;
	double time_int, tmp;
	double theta, m, Ksat, Nout;
	double NO3_out, NH4_out, DON_out, DOC_out;
	double return_flow;
	double water_balance, infiltration;
	double innundation_depth;
	double basin_outflow;
	double basin_rz_storage;
	double basin_unsat_storage;
	double basin_sat_deficit;
	double basin_return_flow;
	double basin_detention_store;
	double basin_area;
	double preday_basin_unsat_storage;
	double preday_basin_rz_storage;
	double preday_basin_sat_deficit;
	//double preday_sat_deficit;
	double preday_basin_return_flow;
	double preday_basin_detention_store;
	double add_field_capacity, rz_drainage, unsat_drainage;
	double streamflow, Qout, Qin_total, Qstr_total;
	struct patch_object *patch;
	struct hillslope_object *hillslope;
	struct patch_object *neigh;
    double temp, totalfc;
    double excess;
    double areaRatio;
    double detention_store_1;
	/*--------------------------------------------------------------*/
	/*	initializations						*/
	/*--------------------------------------------------------------*/
	grow_flag = command_line[0].grow_flag;
	verbose_flag = command_line[0].verbose_flag;

	time_int = 1.0 / (1.0 * n_timesteps);
	basin_outflow = 0.0;
	basin_area = 0.0;
	basin_unsat_storage = 0.0;
	basin_rz_storage = 0.0;
	basin_sat_deficit = 0.0;
	basin_return_flow = 0.0;
	basin_detention_store = 0.0;
	preday_basin_rz_storage = 0.0;
	preday_basin_unsat_storage = 0.0;
	preday_basin_sat_deficit = 0.0;
	preday_basin_return_flow = 0.0;
	preday_basin_detention_store = 0.0;
	streamflow = 0.0;
	Qin_total = 0.0;
	Qstr_total = 0.0;
	d = 0;
	basin[0].basin_outflow = 0.0;       ////<<---------------------repeated? maybe for multple basin
	basin[0].basin_area = 0.0;
	basin[0].basin_unsat_storage = 0.0;
	basin[0].basin_rz_storage = 0.0;
	basin[0].basin_sat_deficit = 0.0;
	basin[0].basin_return_flow = 0.0;
	basin[0].basin_detention_store = 0.0;
	basin[0].preday_basin_rz_storage = 0.0;
	basin[0].preday_basin_unsat_storage = 0.0;
	basin[0].preday_basin_sat_deficit = 0.0;
	basin[0].preday_basin_return_flow = 0.0;
	basin[0].preday_basin_detention_store = 0.0;	
	
	// Note: this assumes that the set of patches in the surface routing table is identical to
	//       the set of patches in the subsurface flow table
	for (i = 0; i < basin->route_list->num_patches; i++) {
		patch = basin->route_list->list[i];
		patch[0].streamflow = 0.0;
		patch[0].return_flow = 0.0;
		patch[0].base_flow = 0.0;
        
        patch[0].stormdrained = 0.0;//outlet
        patch[0].stormdrained_NO3 = 0.0;//outlet
        patch[0].stormdrained_NH4 = 0.0;//outlet
        patch[0].stormdrained_DON = 0.0;//outlet
        patch[0].stormdrained_DOC = 0.0;//outlet
        patch[0].stormdrainYield = 0.0; // being calculated in update_drainage
        //patch[0].pipedrainYield = 0.0; // being calculated in patch_daily_F
        
        patch[0].fromSTREAM_Q = 0.0;
        patch[0].fromSTREAM_NO3 = 0.0;
        patch[0].fromSTREAM_NH4 = 0.0;
        patch[0].fromSTREAM_DON = 0.0;
        patch[0].fromSTREAM_DOC = 0.0;
        patch[0].fromLAND_Q = 0.0;
        patch[0].fromLAND_NO3 = 0.0;
        patch[0].fromLAND_NH4 = 0.0;
        patch[0].fromLAND_DON = 0.0;
        patch[0].fromLAND_DOC = 0.0;
        patch[0].fromRIPARIAN_Q = 0.0;
        patch[0].fromRIPARIAN_NO3 = 0.0;
        patch[0].fromRIPARIAN_NH4 = 0.0;
        patch[0].fromRIPARIAN_DON = 0.0;
        patch[0].fromRIPARIAN_DOC = 0.0;
        
        patch[0].fromSTREAM_surfsubQ = 0.0;
        patch[0].fromSTREAM_surfsubNO3 = 0.0;
        patch[0].fromSTREAM_surfsubNH4 = 0.0;
        patch[0].fromSTREAM_surfsubDON = 0.0;
        patch[0].fromSTREAM_surfsubDOC = 0.0;
        patch[0].fromLAND_surfsubQ = 0.0;
        patch[0].fromLAND_surfsubNO3 = 0.0;
        patch[0].fromLAND_surfsubNH4 = 0.0;
        patch[0].fromLAND_surfsubDON = 0.0;
        patch[0].fromLAND_surfsubDOC = 0.0;
        patch[0].fromRIPARIAN_surfsubQ = 0.0;
        patch[0].fromRIPARIAN_surfsubNO3 = 0.0;
        patch[0].fromRIPARIAN_surfsubNH4 = 0.0;
        patch[0].fromRIPARIAN_surfsubDON = 0.0;
        patch[0].fromRIPARIAN_surfsubDOC = 0.0;
        
        
		patch[0].infiltration_excess = 0.0;
		basin[0].preday_basin_rz_storage += patch[0].rz_storage * patch[0].area;
		basin[0].preday_basin_unsat_storage += patch[0].unsat_storage * patch[0].area;
		basin[0].preday_basin_sat_deficit += patch[0].sat_deficit * patch[0].area;
		basin[0].preday_basin_return_flow += patch[0].return_flow * patch[0].area;
		basin[0].preday_basin_detention_store += patch[0].detention_store
				* patch[0].area;
		basin[0].basin_area += patch[0].area;
		patch[0].Qin_total = 0.0;
		patch[0].Qout_total = 0.0;
        patch[0].surface_Qin_total = 0.0;
        patch[0].surface_Qout_total = 0.0;
		patch[0].Qin = 0.0;
		patch[0].Qout = 0.0;
		patch[0].surface_Qin = 0.0;
		patch[0].surface_Qout = 0.0;
		
		patch[0].overland_flow = 0.0; ///<<--------------- being customized to track return_flow, Sept 7
        
        
		patch[0].interim_sat = patch[0].sat_deficit - patch[0].unsat_storage; // ok not sure what this does
        
		if (grow_flag > 0) {
			patch[0].soil_ns.NO3_Qin = 0.0;
			patch[0].soil_ns.NO3_Qout = 0.0;
			patch[0].soil_ns.NH4_Qin = 0.0;
			patch[0].soil_ns.NH4_Qout = 0.0;
			patch[0].soil_ns.NO3_Qin_total = 0.0;
			patch[0].soil_ns.NO3_Qout_total = 0.0;
			patch[0].soil_ns.NH4_Qin_total = 0.0;
			patch[0].soil_ns.NH4_Qout_total = 0.0;
			patch[0].streamflow_DON = 0.0;
			patch[0].streamflow_DOC = 0.0;
			patch[0].streamflow_NO3 = 0.0;
			patch[0].streamflow_NH4 = 0.0;
			patch[0].soil_ns.DON_Qin_total = 0.0;
			patch[0].soil_ns.DON_Qout_total = 0.0;
			patch[0].soil_cs.DOC_Qin_total = 0.0;
			patch[0].soil_cs.DOC_Qout_total = 0.0;
			patch[0].surface_DON_Qin_total = 0.0;
			patch[0].surface_DON_Qout_total = 0.0;
			patch[0].surface_DOC_Qin_total = 0.0;
			patch[0].surface_DOC_Qout_total = 0.0;
			patch[0].soil_ns.leach = 0.0;
			patch[0].surface_ns_leach = 0.0;
			patch[0].soil_ns.DON_Qout = 0.0;
			patch[0].soil_ns.DON_Qin = 0.0;
			patch[0].soil_cs.DOC_Qout = 0.0;
			patch[0].soil_cs.DOC_Qin = 0.0;
			patch[0].surface_DON_Qout = 0.0;
			patch[0].surface_DON_Qin = 0.0;
			patch[0].surface_DOC_Qout = 0.0;
			patch[0].surface_DOC_Qin = 0.0;
			
			patch[0].streamNO3_from_surface	= 0.0;
			patch[0].streamNO3_from_sub = 0.0;
		}// growth flag
        
    }// reset values in all patches; outside of the k hourly loop


    
	/*--------------------------------------------------------------*/
	/*	calculate Qout for each patch and add appropriate	*/
	/*	proportion of subsurface outflow to each neighbour	*/
	/*--------------------------------------------------------------*/
	for (k = 0; k < n_timesteps; k++) {
        
		//within hourly loop, looping through the patches to reset fluxes and calling update_drainage_xxx()
		for (i = 0; i < basin->route_list->num_patches; i++) {
			patch = basin->route_list->list[i];
            
            patch[0].hourly_subsur2stream_flow = 0;
			patch[0].hourly_sur2stream_flow = 0;
			patch[0].hourly_stream_flow = 0;
			patch[0].hourly[0].streamflow_NO3 = 0;
			patch[0].hourly[0].streamflow_NO3_from_sub = 0;
			patch[0].hourly[0].streamflow_NO3_from_surface = 0;			
			/*--------------------------------------------------------------*/
			/*	for roads, saturated throughflow beneath road cut	*/
			/*	is routed to downslope patches; saturated throughflow	*/
			/*	above the cut and overland flow is routed to the stream	*/
			/*								*/
			/*	for streams, no routing - all exported from basin	*/
			/*								*/
			/*	regular land patches - route to downslope neighbours    */
			/*--------------------------------------------------------------*/
            
            if(patch[0].soil_defaults[0][0].soil_water_cap+ZERO < patch[0].sat_deficit && k==0){
                printf("sub routing (%d,%d) %f %f %f\n",
                       patch[0].ID, k,
                       patch[0].soil_defaults[0][0].soil_water_cap, patch[0].sat_deficit,
                       patch[0].constraintWaterTableTopDepth_def);
            }//debug
            if(patch[0].soil_defaults[0][0].soil_water_cap+ZERO < patch[0].sat_deficit || patch[0].sat_deficit!=patch[0].sat_deficit){
                printf("sub routing (%d,%d) %f %f %f\n",
                       patch[0].ID, k,
                       patch[0].soil_defaults[0][0].soil_water_cap, patch[0].sat_deficit,
                       patch[0].constraintWaterTableTopDepth_def);
            }//debug
            if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
               patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
               patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
               patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
               patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
               patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
               patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
               patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
               patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0){
                printf("sub routing1 (%d,%d) [%e %e %e] [%e %e %e] [%e %e %e]\n",
                       patch[0].ID, k,
                       patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
                       patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
                       patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);
            }//debug
//            if(patch[0].ID==193917 || patch[0].ID==131724 || patch[0].ID==182853 || patch[0].ID==167273){
//                printf("sub routing1 (forced) (%d,%d) [%e %e %e] [%e %e %e] [%e %e %e]\n",
//                       patch[0].ID, k,
//                       patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
//                       patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
//                       patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);
//            }//debug
            
            
            // need to be careful here: patch[0].sat_deficit could be negative.
            if(patch[0].sat_deficit >= 0){
                patch[0].available_soil_water = patch[0].soil_defaults[0][0].soil_water_cap - patch[0].sat_deficit;
                patch[0].sat_def_pct = patch[0].sat_deficit * patch[0].soil_defaults[0][0].max_sat_def_1;
                patch[0].sat_def_pct_index = (int)(patch[0].sat_def_pct*1000);
                patch[0].sat_def_pct_indexM = 1000*(patch[0].sat_def_pct - patch[0].sat_def_pct_index*0.001);
                
                patch[0].sat_deficit_z = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index];
                
                
                patch[0].preday_sat_deficit = patch[0].sat_deficit;
                patch[0].preday_sat_deficit_z = patch[0].sat_deficit_z;
                patch[0].preday_totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index]; // total fc before ET water consumption to SAT
            }else{
                // surface
                patch[0].available_soil_water = patch[0].soil_defaults[0][0].soil_water_cap;
                patch[0].sat_deficit_z = patch[0].sat_deficit;
                patch[0].sat_def_pct = 0.0;
                patch[0].sat_def_pct_index = 0;
                patch[0].sat_def_pct_indexM = 0;
                
                patch[0].preday_sat_deficit = patch[0].sat_deficit;
                patch[0].preday_sat_deficit_z = patch[0].sat_deficit_z;
                patch[0].preday_totalfc = patch[0].soil_defaults[0][0].fc1_0z[0];
            }//ifelse
            
            // --- now going into the flow extraction
			if ((patch[0].drainage_type == ROAD) && (command_line[0].road_flag == 1)) {
                // road_flag is always on from commandline
				update_drainage_road(patch, command_line, time_int, verbose_flag);
			} else if (patch[0].drainage_type == STREAM) {
				update_drainage_stream(patch, command_line, time_int, verbose_flag);
			} else {
                // any other codes, e.g.,PAVEDROAD, IMP, RIPARIAN // Sept 7
				//update_drainage_land(patch, command_line, time_int, verbose_flag);
                update_drainage_land(patch, command_line, time_int, k);
			}
            
            if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
               patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
               patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
               patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
               patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
               patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
               patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
               patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
               patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0){
                printf("sub routing2 (%d,%d) [%e %e %e] [%e %e %e] [%e %e %e]\n",
                       patch[0].ID, k,
                       patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
                       patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
                       patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);
            }//debug
            
//            if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity){
//                printf("subsurface_routing(1,%d): %d %lf %lf %lf %lf\n", patch[0].ID,
//                       k, patch[0].sat_deficit, patch[0].sat_deficit_z,
//                       patch[0].rootzone.field_capacity, patch[0].field_capacity);
//            }//debug
		} // end of loop i, the patch loop within an hour.

        /*--------------------------------------------------------------*/
        // within a hourly loop; applying the lateral IN/OUT on each patch (below)
        /*--------------------------------------------------------------*/
        for (i = 0; i < basin->route_list->num_patches; i++) {
            patch = basin->route_list->list[i];

            // problem:
            // 1. route_to_patch = time_int * compute_varbased_flow()
            // 2. route_to_patch += (extrawater>return_flow? (extrawater-return_flow)*patch[0].area : 0.0);
            // the extra water is supposed to be surface water
            //patch[0].satzZ_balance = min(0.0, patch[0].available_soil_water - (patch[0].Qout - patch[0].Qin)); // [negative-0] only
            patch[0].sat_deficit += patch[0].Qout - patch[0].Qin + patch[0].satzZ_balance; // hourly; subsurface only
                // got patch[0].sat_deficit = NaN problem
            if(patch[0].sat_deficit>patch[0].soil_defaults[0][0].soil_water_cap){
                if(patch[0].sat_deficit>patch[0].soil_defaults[0][0].soil_water_cap+ZERO)
                    printf("patch routing trouble (%d,%d) %f %f (%f %f %f, %f, %f)\n",
                           patch[0].ID, k,
                           patch[0].soil_defaults[0][0].soil_water_cap,
                           patch[0].sat_deficit,
                           patch[0].Qout, patch[0].Qin, patch[0].satzZ_balance,
                           patch[0].available_soil_water,
                           patch[0].constraintWaterTableTopDepth_def
                           );
                patch[0].sat_deficit=patch[0].soil_defaults[0][0].soil_water_cap;
            }//if
            // how to fix this?
            // (soil_defaults[0][0].soil_water_cap-patch[0].sat_deficit);
            // potential solution: satzZ_balance; [negative-0];
            // patch[0].Qout+satzZ_balance = pot. Qout
            //
            
            // need to be careful here: patch[0].sat_deficit could be negative.
            if(patch[0].sat_deficit >= 0){
                patch[0].available_soil_water = patch[0].soil_defaults[0][0].soil_water_cap - patch[0].sat_deficit;
                patch[0].sat_def_pct = patch[0].sat_deficit * patch[0].soil_defaults[0][0].max_sat_def_1;
                patch[0].sat_def_pct_index = (int)(patch[0].sat_def_pct*1000);
                patch[0].sat_def_pct_indexM = 1000*(patch[0].sat_def_pct - patch[0].sat_def_pct_index*0.001);
                
                patch[0].sat_deficit_z = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index];
            }else{
                // surface
                patch[0].available_soil_water = patch[0].soil_defaults[0][0].soil_water_cap;
                patch[0].sat_deficit_z = patch[0].sat_deficit;
                patch[0].sat_def_pct = 0.0;
                patch[0].sat_def_pct_index = 0;
                patch[0].sat_def_pct_indexM = 0;
            }
            // fc & SatPct
            totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index];
            totalfc *= (1.0-patch[0].basementFrac); // <---- second thought on this, Oct 8, 2019; basement is 3m at most
            
            if (patch[0].sat_deficit < ZERO) {
                //patch[0].aboveWT_SatPct = 1.0;
                //patch[0].rootzone.SatPct = 1.0;
                patch[0].rootzone.field_capacity = 0.0;
                patch[0].field_capacity = 0.0;
            } else {
                patch[0].rootzone.field_capacity = totalfc * patch[0].zeroRootCoef * (patch[0].rootzone_scale_ref*patch[0].rootzone_end_reffc[patch[0].sat_def_pct_index] + (1.0-patch[0].rootzone_scale_ref)*patch[0].rootzone_start_reffc[patch[0].sat_def_pct_index]);
                patch[0].rootzone.field_capacity = min(patch[0].rootzone.field_capacity,patch[0].rootzone.potential_sat);
                patch[0].field_capacity = max(0.0,min(patch[0].sat_deficit-patch[0].rootzone.potential_sat, totalfc - patch[0].rootzone.field_capacity));
            }//if else

            
            if (grow_flag > 0) {

                if(patch[0].sat_NO3*1.001 < patch[0].soil_ns.NO3_Qout)
                    printf("sub_routing _Qout[%d]{%e,%e}\n",
                           patch[0].ID, patch[0].sat_NO3, patch[0].soil_ns.NO3_Qout);
                
                patch[0].soil_ns.NO3_Qout = min(patch[0].sat_NO3, patch[0].soil_ns.NO3_Qout);
                patch[0].soil_ns.NH4_Qout = min(patch[0].sat_NH4, patch[0].soil_ns.NH4_Qout);
                patch[0].soil_cs.DOC_Qout = min(patch[0].sat_DOC, patch[0].soil_cs.DOC_Qout);
                patch[0].soil_ns.DON_Qout = min(patch[0].sat_DON, patch[0].soil_ns.DON_Qout);
                patch[0].sat_NO3 += (patch[0].soil_ns.NO3_Qin - patch[0].soil_ns.NO3_Qout);
                patch[0].sat_NH4 += (patch[0].soil_ns.NH4_Qin - patch[0].soil_ns.NH4_Qout);
                patch[0].sat_DOC += (patch[0].soil_cs.DOC_Qin - patch[0].soil_cs.DOC_Qout);
                patch[0].sat_DON += (patch[0].soil_ns.DON_Qin - patch[0].soil_ns.DON_Qout);
                
                if( patch[0].sat_deficit_z > patch[0].preday_sat_deficit_z){
                    // water table drops
                    //double sat_leftbehind_frac = (exp(-patch[0].preday_sat_deficit_z/patch[0].soil_defaults[0][0].porosity_decay)-exp(-patch[0].sat_deficit_z/patch[0].soil_defaults[0][0].porosity_decay)) / (exp(-patch[0].preday_sat_deficit_z/patch[0].soil_defaults[0][0].porosity_decay)-exp(-patch[0].soil_defaults[0][0].soil_depth/patch[0].soil_defaults[0][0].porosity_decay));
                        // 0 < sat_leftbehind_frac < 1
                    
                    double sat_leftbehind_frac = (patch[0].sat_deficit - patch[0].preday_sat_deficit) / (patch[0].soil_defaults[0][0].soil_water_cap - patch[0].preday_sat_deficit);// wrong here?
                    
                    if(sat_leftbehind_frac<0 || sat_leftbehind_frac>1.0 || patch[0].sat_NO3<0 || patch[0].sat_NO3!=patch[0].sat_NO3) printf("sub_routing (%d) [%e,%e,%e], %f %f %f \n", patch[0].ID, sat_leftbehind_frac, patch[0].soil_ns.nitrate, patch[0].sat_NO3,
                        patch[0].sat_deficit,  patch[0].preday_sat_deficit, patch[0].soil_defaults[0][0].soil_water_cap);
                    
                    // for negative "soil_ns.nitrate" problem
                    patch[0].soil_ns.nitrate += patch[0].sat_NO3 * sat_leftbehind_frac;
                    patch[0].soil_ns.sminn += patch[0].sat_NH4 * sat_leftbehind_frac;
                    patch[0].soil_cs.DOC += patch[0].sat_DOC * sat_leftbehind_frac;
                    patch[0].soil_ns.DON += patch[0].sat_DON * sat_leftbehind_frac;
                    
                    patch[0].sat_NO3 *= 1.0 - sat_leftbehind_frac;
                    patch[0].sat_NH4 *= 1.0 - sat_leftbehind_frac;
                    patch[0].sat_DOC *= 1.0 - sat_leftbehind_frac;
                    patch[0].sat_DON *= 1.0 - sat_leftbehind_frac;
                    
                }// water table drops
            }//grow_flag
            if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
               patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
               patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
               patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
               patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
               patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
               patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
               patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
               patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0){
                printf("sub routing3 (%d,%d) [%e %e %e] [%e %e %e] [%e %e %e]\n",
                       patch[0].ID, k,
                       patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
                       patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
                       patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);
            }//debug

            /*--------------------------------------------------------------*/
            /*    reset iterative  patch fluxes to zero & do some acummulation */
            /*--------------------------------------------------------------*/
            patch[0].soil_ns.leach += (patch[0].soil_ns.DON_Qout
                    + patch[0].soil_ns.NH4_Qout + patch[0].soil_ns.NO3_Qout
                    - patch[0].soil_ns.NH4_Qin - patch[0].soil_ns.NO3_Qin
                    - patch[0].soil_ns.DON_Qin);
            patch[0].surface_ns_leach += ((patch[0].surface_NO3_Qout
                    - patch[0].surface_NO3_Qin)
                    + (patch[0].surface_NH4_Qout - patch[0].surface_NH4_Qin)
                    + (patch[0].surface_DON_Qout - patch[0].surface_DON_Qin));
            patch[0].Qin_total += patch[0].Qin;
            patch[0].Qout_total += patch[0].Qout;
            patch[0].surface_Qin_total += patch[0].surface_Qin;
            patch[0].surface_Qout_total += patch[0].surface_Qout;

            // these are being reset every hour
            patch[0].surface_Qin = 0.0;
            patch[0].surface_Qout = 0.0;
            patch[0].Qin = 0.0;
            patch[0].Qout = 0.0;
            if (grow_flag > 0) {
                patch[0].soil_cs.DOC_Qin_total += patch[0].soil_cs.DOC_Qin;
                patch[0].soil_cs.DOC_Qout_total += patch[0].soil_cs.DOC_Qout;
                patch[0].soil_ns.NH4_Qin_total += patch[0].soil_ns.NH4_Qin;
                patch[0].soil_ns.NH4_Qout_total += patch[0].soil_ns.NH4_Qout;
                patch[0].soil_ns.NO3_Qin_total += patch[0].soil_ns.NO3_Qin;
                patch[0].soil_ns.NO3_Qout_total += patch[0].soil_ns.NO3_Qout;
                patch[0].soil_ns.DON_Qin_total += patch[0].soil_ns.DON_Qin;
                patch[0].soil_ns.DON_Qout_total += patch[0].soil_ns.DON_Qout;
                patch[0].surface_DON_Qin_total += patch[0].surface_DON_Qin;
                patch[0].surface_DON_Qout_total += patch[0].surface_DON_Qout;
                patch[0].surface_DOC_Qin_total += patch[0].surface_DOC_Qin;
                patch[0].surface_DOC_Qout_total += patch[0].surface_DOC_Qout;

                patch[0].soil_ns.NH4_Qin = 0.0;
                patch[0].soil_ns.NH4_Qout = 0.0;
                patch[0].soil_ns.NO3_Qin = 0.0;
                patch[0].soil_ns.NO3_Qout = 0.0;
                patch[0].soil_ns.DON_Qout = 0.0;
                patch[0].soil_ns.DON_Qin = 0.0;
                patch[0].soil_cs.DOC_Qout = 0.0;
                patch[0].soil_cs.DOC_Qin = 0.0;
                patch[0].surface_NH4_Qout = 0.0;
                patch[0].surface_NH4_Qin = 0.0;
                patch[0].surface_NO3_Qout = 0.0;
                patch[0].surface_NO3_Qin = 0.0;
                patch[0].surface_DON_Qout = 0.0;
                patch[0].surface_DON_Qin = 0.0;
                patch[0].surface_DOC_Qout = 0.0;
                patch[0].surface_DOC_Qin = 0.0;

            }//reset done
                    
            /*--------------------------------------------------------------*/
            /*     leave behind field capacity            */
            /*    if sat deficit has been lowered            */
            /*    this should be an interactive process, we will use     */
            /*    0th order approximation                    */
            /*     we do not do this once sat def is below 0.9 soil depth    */
            /*     we use 0.9 to prevent numerical instability        */
            /*--------------------------------------------------------------*/
            if ((patch[0].sat_deficit_z > patch[0].preday_sat_deficit_z) && (patch[0].sat_deficit_z < patch[0].soil_defaults[0][0].soil_depth * 0.9)) {
            
                add_field_capacity = min( patch[0].available_soil_water, max(0.0, totalfc - patch[0].preday_totalfc));
                patch[0].sat_deficit += add_field_capacity;// could go deeper than soil depth
                
                if ((patch[0].sat_deficit_z > patch[0].rootzone.depth) && (patch[0].preday_sat_deficit_z > patch[0].rootzone.depth))
                    patch[0].unsat_storage += add_field_capacity;
                else
                    patch[0].rz_storage += add_field_capacity*(1.0-patch[0].basementFrac);

                
                // need to be careful here: patch[0].sat_deficit could be negative.
                if(patch[0].sat_deficit >= 0){
                    patch[0].available_soil_water = patch[0].soil_defaults[0][0].soil_water_cap - patch[0].sat_deficit;
                    patch[0].sat_def_pct = patch[0].sat_deficit * patch[0].soil_defaults[0][0].max_sat_def_1;
                    patch[0].sat_def_pct_index = (int)(patch[0].sat_def_pct*1000);
                    patch[0].sat_def_pct_indexM = 1000*(patch[0].sat_def_pct - patch[0].sat_def_pct_index*0.001);

                    patch[0].sat_deficit_z = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index];
                }else{
                    // surface
                    patch[0].available_soil_water = patch[0].soil_defaults[0][0].soil_water_cap;
                    patch[0].sat_deficit_z = patch[0].sat_deficit;
                    patch[0].sat_def_pct = 0.0;
                    patch[0].sat_def_pct_index = 0;
                    patch[0].sat_def_pct_indexM = 0;
                }

                totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index]; // total fc after ET water consumption to SAT
                totalfc *= (1.0-patch[0].basementFrac); // <---- second thought on this, Oct 8, 2019; basement is 3m at most

                if (patch[0].sat_deficit < ZERO) {
                    //patch[0].aboveWT_SatPct = 1.0;
                    //patch[0].rootzone.SatPct = 1.0;
                    patch[0].rootzone.field_capacity = 0.0;
                    patch[0].field_capacity = 0.0;
                } else {
                    patch[0].rootzone.field_capacity = totalfc * patch[0].zeroRootCoef * (patch[0].rootzone_scale_ref*patch[0].rootzone_end_reffc[patch[0].sat_def_pct_index] + (1.0-patch[0].rootzone_scale_ref)*patch[0].rootzone_start_reffc[patch[0].sat_def_pct_index]);
                    patch[0].rootzone.field_capacity = min(patch[0].rootzone.field_capacity,patch[0].rootzone.potential_sat);
                    patch[0].field_capacity = max(0.0,min(patch[0].sat_deficit-patch[0].rootzone.potential_sat, totalfc - patch[0].rootzone.field_capacity));
                }//if else
            }// if left-behind water
            
            
            
            //---- adding FC back to rtz; rtz was previously set to zero due to sat_def_z <= 0
            if (patch[0].rootzone.depth > ZERO) {
                if ((patch[0].sat_deficit > ZERO) && (patch[0].rz_storage == 0.0)) {
                    
                    // solution is to track "patch[0].soil_defaults[0][0].sat_zZ[ii]" = patch[0].available_soil_water
                    patch[0].sat_deficit += patch[0].rootzone.field_capacity; // patch[0].sat_deficit<patch[0].soil_defaults[0][0].soil_water_cap
                    patch[0].rz_storage += patch[0].rootzone.field_capacity;
                    
                    patch[0].sat_deficit += patch[0].field_capacity;
                    patch[0].unsat_storage += patch[0].field_capacity;
                }
            } else {
                // ---- adding FC back to unsat; unsat was previously set to zero due to sat_def_z <= 0
                // ---- however, why is it "else"? does not seem water balance to me; (Jan 2, 2019)
                if ((patch[0].sat_deficit > ZERO) && (patch[0].unsat_storage == 0.0)) {

                    patch[0].sat_deficit += patch[0].field_capacity;
                    patch[0].unsat_storage += patch[0].field_capacity;
                }
            }//--------------------------- field capacity
           
    //        // need to be careful here: patch[0].sat_deficit could be negative.
    //        if(patch[0].sat_deficit >= 0){
    //            patch[0].sat_def_pct = patch[0].sat_deficit * patch[0].soil_defaults[0][0].max_sat_def_1;
    //            patch[0].sat_def_pct_index = (int)(patch[0].sat_def_pct*1000);
    //            patch[0].sat_def_pct_indexM = 1000*(patch[0].sat_def_pct - patch[0].sat_def_pct_index*0.001);
    //
    //            patch[0].sat_deficit_z = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index];
    //        }else{
    //            // surface
    //            patch[0].sat_deficit_z = patch[0].sat_deficit;
    //            patch[0].sat_def_pct = 0.0;
    //            patch[0].sat_def_pct_index = 0;
    //            patch[0].sat_def_pct_indexM = 0;
    //        }
    //
    //        // fc & SatPct
    //        totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index];
    //        totalfc *= (1.0-patch[0].basementFrac); // <---- second thought on this, Oct 8, 2019; basement is 3m at most
    //
    //        if (patch[0].sat_deficit < ZERO) {
    //            //patch[0].aboveWT_SatPct = 1.0;
    //            //patch[0].rootzone.SatPct = 1.0;
    //            patch[0].rootzone.field_capacity = 0.0;
    //            patch[0].field_capacity = 0.0;
    //        } else {
    //            patch[0].rootzone.field_capacity = totalfc * patch[0].zeroRootCoef * (patch[0].rootzone_scale_ref*patch[0].rootzone_end_reffc[patch[0].sat_def_pct_index] + (1.0-patch[0].rootzone_scale_ref)*patch[0].rootzone_start_reffc[patch[0].sat_def_pct_index]);
    //            patch[0].rootzone.field_capacity = min(patch[0].rootzone.field_capacity,patch[0].rootzone.potential_sat);
    //            patch[0].field_capacity = max(0.0,min(patch[0].sat_deficit-patch[0].rootzone.potential_sat, totalfc - patch[0].rootzone.field_capacity));
    //        }//if else
            
            /*--------------------------------------------------------------*/
            /*    final overland flow routing (vertially)           */
            /*--------------------------------------------------------------*/
            excess = patch[0].rz_storage+patch[0].unsat_storage - patch[0].sat_deficit + patch[0].constraintWaterTableTopDepth_def;
            if(excess > 0) {
                  patch[0].detention_store += excess;
                  patch[0].sat_deficit -= patch[0].unsat_storage + patch[0].rz_storage - excess;
                  patch[0].unsat_storage = 0.0; // converted to be part of sat
                  patch[0].rz_storage = 0.0; // converted to be part of sat
                  
                if(patch[0].sat_deficit <= -1 || patch[0].sat_deficit >= patch[0].soil_defaults[0][0].soil_depth){
                    printf("subsurface_routing(end1,%d): %lf %lf %lf %lf\n",patch[0].ID,
                    patch[0].sat_deficit, patch[0].sat_deficit_z,
                    patch[0].rootzone.field_capacity, patch[0].field_capacity);
                }//debug
                
                  // need to be careful here: patch[0].sat_deficit could be negative.
                  if(patch[0].sat_deficit >= 0){
                      patch[0].sat_def_pct = patch[0].sat_deficit * patch[0].soil_defaults[0][0].max_sat_def_1;
                      patch[0].sat_def_pct_index = (int)(patch[0].sat_def_pct*1000);
                      patch[0].sat_def_pct_indexM = 1000*(patch[0].sat_def_pct - patch[0].sat_def_pct_index*0.001);

                      patch[0].sat_deficit_z = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index];
                  }else{
                      // surface
                      patch[0].sat_deficit_z = patch[0].sat_deficit;
                      patch[0].sat_def_pct = 0.0;
                      patch[0].sat_def_pct_index = 0;
                      patch[0].sat_def_pct_indexM = 0;
                  }
                  
                  totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index]; // total fc after ET water consumption to SAT
                  totalfc *= (1.0-patch[0].basementFrac); // <---- second thought on this, Oct 8, 2019; basement is 3m at most
                  
                  if (patch[0].sat_deficit < ZERO) {
                      //patch[0].aboveWT_SatPct = 1.0;
                      //patch[0].rootzone.SatPct = 1.0;
                      patch[0].rootzone.field_capacity = 0.0;
                      patch[0].field_capacity = 0.0;
                  } else {
                      patch[0].rootzone.field_capacity = totalfc * patch[0].zeroRootCoef * (patch[0].rootzone_scale_ref*patch[0].rootzone_end_reffc[patch[0].sat_def_pct_index] + (1.0-patch[0].rootzone_scale_ref)*patch[0].rootzone_start_reffc[patch[0].sat_def_pct_index]);
                      patch[0].rootzone.field_capacity = min(patch[0].rootzone.field_capacity,patch[0].rootzone.potential_sat);
                      patch[0].field_capacity = max(0.0,min(patch[0].sat_deficit-patch[0].rootzone.potential_sat, totalfc - patch[0].rootzone.field_capacity));
                  }//if else
                
                  if (grow_flag > 0) {
                      
                      Nout =
                              compute_N_leached(
                                      verbose_flag,
                                      patch[0].sat_DOC,
                                      excess,
                                      patch[0].soil_defaults[0][0].DOMdecayRate,
                                      patch[0].soil_defaults[0][0].active_zone_z,
                                      patch[0].soil_defaults[0][0].DOC_adsorption_rate,
                                      26,patch);
                      patch[0].surface_DOC += Nout;
                      patch[0].sat_DOC -= Nout; //patch[0].soil_cs.DOC -= Nout;
                 
                      Nout =
                              compute_N_leached(verbose_flag,
                                      patch[0].sat_DON,
                                      excess,
                                      patch[0].soil_defaults[0][0].DOMdecayRate,
                                      patch[0].soil_defaults[0][0].active_zone_z,
                                      patch[0].soil_defaults[0][0].DON_adsorption_rate,
                                      23,patch);
                      patch[0].surface_DON += Nout;
                      patch[0].sat_DON -= Nout; //patch[0].soil_ns.DON -= Nout;
                  
                      
                      Nout =
                              compute_N_leached(verbose_flag,
                                      patch[0].sat_NO3,
                                      excess,
                                      patch[0].soil_defaults[0][0].NO3decayRate,
                                      patch[0].soil_defaults[0][0].active_zone_z,
                                      patch[0].soil_defaults[0][0].NO3_adsorption_rate,
                                      17,patch);
                      patch[0].surface_NO3 += Nout;
                      patch[0].sat_NO3 -= Nout; //patch[0].soil_ns.nitrate -= Nout;
                  
                      
                      Nout =
                              compute_N_leached(verbose_flag,
                                      patch[0].sat_NH4,
                                      excess,
                                      patch[0].soil_defaults[0][0].NH4decayRate,
                                      patch[0].soil_defaults[0][0].active_zone_z,
                                      patch[0].soil_defaults[0][0].NH4_adsorption_rate,
                                      20,patch);
                      patch[0].surface_NH4 += Nout;
                      patch[0].sat_NH4 -= Nout;//patch[0].soil_ns.sminn -= Nout;
               }//grow_flag
            }//if return flow
            
            /*--------------------------------------------------------------*/
            /*    final overland flow routing (horizontally)              */
            /*--------------------------------------------------------------*/
            excess = patch[0].detention_store - patch[0].landuse_defaults[0][0].detention_store_size* (1.0 - patch[0].Ksat_vertical);
            if ( excess > ZERO && patch[0].detention_store > ZERO ) {
                detention_store_1 = 1.0 / patch[0].detention_store;
                
                if (patch[0].drainage_type == STREAM) {
                    // current patch is STREAM grid
                    patch[0].return_flow += excess; // directly add into return_flow at current patch
                    patch[0].detention_store -= excess;
                    patch[0].Qout_total += excess;
                    patch[0].hourly_sur2stream_flow += excess;
                    excess *= detention_store_1;
                    if (grow_flag > 0) {
                        patch[0].streamflow_DON += excess* patch[0].surface_DON;
                        patch[0].streamflow_DOC += excess* patch[0].surface_DOC;

                        patch[0].streamflow_NO3 += excess* patch[0].surface_NO3;
                        patch[0].streamNO3_from_surface += excess* patch[0].surface_NO3;
                        patch[0].hourly[0].streamflow_NO3 += excess* patch[0].surface_NO3;
                        patch[0].hourly[0].streamflow_NO3_from_surface += excess* patch[0].surface_NO3;

                        patch[0].streamflow_NH4 += excess* patch[0].surface_NH4;
                        patch[0].surface_DON -= excess* patch[0].surface_DON;
                        patch[0].surface_DOC -= excess* patch[0].surface_DOC;
                        patch[0].surface_NO3 -= excess* patch[0].surface_NO3;
                        patch[0].surface_NH4 -= excess* patch[0].surface_NH4;
                    }// growth flag
                    
                } else {
                    // current patch is not a stream grid
                    patch[0].detention_store -= excess;
                    patch[0].Qout_total += excess;
                    /*--------------------------------------------------------------*/
                    /* determine which innundation depth to consider        */
                    /*--------------------------------------------------------------*/
                    if (patch[0].num_innundation_depths > 0) {
                        innundation_depth = patch[0].detention_store;
                        d = 0;
                        while ( (innundation_depth > patch[0].innundation_list[d].critical_depth) &&
                                (d < patch[0].num_innundation_depths - 1)) {
                            d++;
                        }//while
                    } else {
                        d = 0;
                    }// if else

                    for (j = 0; j < patch->surface_innundation_list[d].num_neighbours; j++) {
                        neigh = patch->surface_innundation_list[d].neighbours[j].patch;
                        Qout = excess * patch->surface_innundation_list[d].neighbours[j].gamma;
                        if (grow_flag > 0) {
                            NO3_out = Qout * detention_store_1 * patch[0].surface_NO3;
                            NH4_out = Qout * detention_store_1 * patch[0].surface_NH4;
                            DON_out = Qout * detention_store_1 * patch[0].surface_DON;
                            DOC_out = Qout * detention_store_1 * patch[0].surface_DOC;
                            Nout = NO3_out + NH4_out + DON_out;
                        }//growth flag
                        areaRatio = patch[0].area / neigh[0].area;
                        if (neigh[0].drainage_type == STREAM) {
                            neigh[0].Qin_total += Qout * areaRatio;
                            neigh[0].return_flow += Qout * areaRatio;
                            if (grow_flag > 0) {
                                neigh[0].streamflow_DOC += (DOC_out * areaRatio);
                                neigh[0].streamflow_DON += (DON_out * areaRatio);

                                neigh[0].streamflow_NO3 += (NO3_out * areaRatio);
                                neigh[0].streamNO3_from_surface +=(NO3_out * areaRatio);
                                neigh[0].hourly[0].streamflow_NO3 += (NO3_out * areaRatio);
                                neigh[0].hourly[0].streamflow_NO3_from_sub +=(NO3_out * areaRatio);

                                neigh[0].streamflow_NH4 += (NH4_out * areaRatio);
                                neigh[0].surface_ns_leach += (Nout * areaRatio);
                            }//growth_flag
                        } else {
                            neigh[0].Qin_total += Qout * areaRatio;
                            neigh[0].detention_store += Qout * areaRatio;
                            if (grow_flag > 0) {
                                neigh[0].surface_DOC += (DOC_out * areaRatio);
                                neigh[0].surface_DON += (DON_out * areaRatio);
                                neigh[0].surface_NO3 += (NO3_out * areaRatio);
                                neigh[0].surface_ns_leach -= (Nout * areaRatio);
                                neigh[0].surface_NH4 += (NH4_out * areaRatio);
                            }//growth flag
                        }//if else
                    }// end of for loop j
                    if (grow_flag > 0) {
                        patch[0].surface_DOC -= excess * detention_store_1 * patch[0].surface_DOC;
                        patch[0].surface_DON -= excess * detention_store_1 * patch[0].surface_DON;
                        patch[0].surface_NO3 -= excess * detention_store_1 * patch[0].surface_NO3;

                        patch[0].surface_NH4 -= excess * detention_store_1 * patch[0].surface_NH4;
                        patch[0].surface_ns_leach += excess * detention_store_1 * patch[0].surface_NO3;
                    }// growth flag
                }// if else
            }//return_flow
            
      

        }// for loop i for patches
	} /* end k hourly step */


    // outside the hourly loop, the end of day
    for (i = 0; i < basin->route_list->num_patches; i++) {
        patch = basin->route_list->list[i];
        // this block should be done once at the end of day because these processes have been done in drainage_xxx(). July 12, 2019
        // this is necessarily because the aggregation process is outside of drainage_xxx() and balancing processes are in drainage_xxx().
        // this block is acting on current patch only.
        //
        // return_flow calculation one last time at the end of day (for lateral IN/OUT balance) -- indeed we do this once more at the end (but not hourly)
        // start to calculate the new "field capacity" (not done in update_drain_land)
        // patch[0].sat_deficit_z is updated above after lateral IN/OUT
        // this is not smooth by daily step; that the "extra water" need to save for the next next day (not going to fix it)

        /*------------------------------------------------------------------------------*/
        /*        Just add the infiltration to the rz_storage and unsat_storage    */
        /*------------------------------------------------------------------------------*/
        if (patch[0].detention_store > ZERO){
            
            infiltration = compute_infiltration(
                command_line[0].verbose_flag,
                patch[0].sat_deficit_z,
                0.0, //patch[0].aboveWT_SatPct, // initiated in daily_I()
                patch[0].Ksat_vertical, // 1- impervious
                patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index],
                patch[0].rz_storage+patch[0].unsat_storage,
                patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].sat_def_0zm[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].sat_def_0zm[patch[0].sat_def_pct_index],
                patch[0].sat_deficit,
                patch[0].detention_store,
                time_int,
                patch[0].soil_defaults[0][0].psi_air_entry);
            
        } else infiltration = 0.0;
        
        if ((grow_flag > 0) && (infiltration > ZERO)) {
            patch[0].soil_ns.DON += ((infiltration / patch[0].detention_store) * patch[0].surface_DON);
            patch[0].surface_DON -= ((infiltration / patch[0].detention_store) * patch[0].surface_DON);
            
            patch[0].soil_cs.DOC += ((infiltration / patch[0].detention_store) * patch[0].surface_DOC);
            patch[0].surface_DOC -= ((infiltration / patch[0].detention_store) * patch[0].surface_DOC);
            
            patch[0].soil_ns.nitrate += ((infiltration / patch[0].detention_store) * patch[0].surface_NO3);
            patch[0].surface_NO3 -= ((infiltration / patch[0].detention_store) * patch[0].surface_NO3);
            
            patch[0].soil_ns.sminn += ((infiltration / patch[0].detention_store) * patch[0].surface_NH4);
            patch[0].surface_NH4 -= ((infiltration / patch[0].detention_store) * patch[0].surface_NH4);
        }
        
        if (infiltration > patch[0].sat_deficit - patch[0].unsat_storage - patch[0].rz_storage) {
            /*--------------------------------------------------------------*/
            /*        Yes the unsat zone will be filled so we may    */
            /*        as well treat the unsat_storage and infiltration*/
            /*        as water added to the water table.        */
            /*--------------------------------------------------------------*/
            patch[0].sat_deficit -= (infiltration + patch[0].unsat_storage + patch[0].rz_storage);
            patch[0].unsat_storage = 0.0;
            patch[0].rz_storage = 0.0;
            patch[0].field_capacity = 0.0;
            patch[0].rootzone.field_capacity = 0.0;
        } else if ((patch[0].sat_deficit > patch[0].rootzone.potential_sat) && (infiltration > patch[0].rootzone.potential_sat - patch[0].rz_storage)) {
            patch[0].unsat_storage += infiltration - (patch[0].rootzone.potential_sat - patch[0].rz_storage);
            patch[0].rz_storage = patch[0].rootzone.potential_sat;
        } else if ((patch[0].sat_deficit > patch[0].rootzone.potential_sat) && (infiltration <= patch[0].rootzone.potential_sat - patch[0].rz_storage)) {
            /* Only rootzone layer saturated - perched water table case */
            patch[0].rz_storage += infiltration;
        } else if ((patch[0].sat_deficit <= patch[0].rootzone.potential_sat) && (infiltration <= patch[0].sat_deficit - patch[0].rz_storage - patch[0].unsat_storage)) {
            patch[0].rz_storage += patch[0].unsat_storage;
            /* transfer left water in unsat storage to rootzone layer */
            patch[0].unsat_storage = 0;
            patch[0].rz_storage += infiltration;
            patch[0].field_capacity = 0;
        }

        if (patch[0].sat_deficit < 0.0) {
            patch[0].detention_store -= (patch[0].sat_deficit - patch[0].unsat_storage);
            patch[0].sat_deficit = 0.0;
            patch[0].unsat_storage = 0.0;
        }

        patch[0].detention_store -= infiltration;
        
        
        
        /*------------------------------------------------------------------------------*/
        /*        drainge    */
        /*------------------------------------------------------------------------------*/
        // need to be careful here: patch[0].sat_deficit could be negative.
        if(patch[0].sat_deficit >= 0){
            patch[0].sat_def_pct = patch[0].sat_deficit * patch[0].soil_defaults[0][0].max_sat_def_1;
            patch[0].sat_def_pct_index = (int)(patch[0].sat_def_pct*1000);
            patch[0].sat_def_pct_indexM = 1000*(patch[0].sat_def_pct - patch[0].sat_def_pct_index*0.001);
            
            patch[0].sat_deficit_z = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index];
        }else{
            // surface
            patch[0].sat_deficit_z = patch[0].sat_deficit;
            patch[0].sat_def_pct = 0.0;
            patch[0].sat_def_pct_index = 0;
            patch[0].sat_def_pct_indexM = 0;
        }
        
        
        // fc & SatPct
        totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index];
        totalfc *= (1.0-patch[0].basementFrac); // <---- second thought on this, Oct 8, 2019; basement is 3m at most
        
        if (patch[0].sat_deficit < ZERO) {
             //patch[0].aboveWT_SatPct = 1.0;
             //patch[0].rootzone.SatPct = 1.0;
            patch[0].rootzone.field_capacity = 0.0;
            patch[0].field_capacity = 0.0;
            unsat_drainage = 0.0;
            rz_drainage = 0.0;
        } else {
             patch[0].rootzone.field_capacity = totalfc * patch[0].zeroRootCoef * (patch[0].rootzone_scale_ref*patch[0].rootzone_end_reffc[patch[0].sat_def_pct_index] + (1.0-patch[0].rootzone_scale_ref)*patch[0].rootzone_start_reffc[patch[0].sat_def_pct_index]);
             patch[0].rootzone.field_capacity = min(patch[0].rootzone.field_capacity,patch[0].rootzone.potential_sat);
             patch[0].field_capacity = max(0.0,min(patch[0].sat_deficit-patch[0].rootzone.potential_sat, totalfc - patch[0].rootzone.field_capacity));

             /// drainage
             rz_drainage = compute_unsat_zone_drainage(
                     command_line[0].verbose_flag,
                     patch[0].soil_defaults[0][0].theta_psi_curve,
                     patch[0].soil_defaults[0][0].pore_size_index,
                     patch[0].rootzone.potential_sat, //patch[0].rootzone.SatPct,
                     patch[0].rootdepth_indexM * patch[0].soil_defaults[0][0].vksat_0zm[patch[0].rootdepth_index+1] + (1.0-patch[0].rootdepth_indexM)* patch[0].soil_defaults[0][0].vksat_0zm[patch[0].rootdepth_index],
                     patch[0].rootdepth_indexM * patch[0].soil_defaults[0][0].vksat_z[patch[0].rootdepth_index+1] + (1.0-patch[0].rootdepth_indexM) * patch[0].soil_defaults[0][0].vksat_z[patch[0].rootdepth_index],
                     patch[0].rz_storage,
                     patch[0].rootzone.field_capacity,
                     patch[0].sat_deficit);
             patch[0].rz_storage -=  rz_drainage;
             patch[0].unsat_storage +=  rz_drainage;
             
            
             unsat_drainage = compute_unsat_zone_drainage(
                     command_line[0].verbose_flag,
                     patch[0].soil_defaults[0][0].theta_psi_curve,
                     patch[0].soil_defaults[0][0].pore_size_index,
                     patch[0].sat_deficit - patch[0].rootzone.potential_sat,//patch[0].aboveWT_SatPct,
                     patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM)* patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index],
                     patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].vksat_z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].vksat_z[patch[0].sat_def_pct_index],
                     patch[0].unsat_storage,
                     patch[0].field_capacity,
                     patch[0].sat_deficit);
             patch[0].unsat_storage -=  unsat_drainage;
             patch[0].sat_deficit -=  unsat_drainage;
             
         }//if else
        // tracking fluxes
        patch[0].unsat_drainage += unsat_drainage;
        patch[0].rz_drainage += rz_drainage;
        
        
        
        // need to be careful here: patch[0].sat_deficit could be negative.
        if(patch[0].sat_deficit >= 0){
            patch[0].sat_def_pct = patch[0].sat_deficit * patch[0].soil_defaults[0][0].max_sat_def_1;
            patch[0].sat_def_pct_index = (int)(patch[0].sat_def_pct*1000);
            patch[0].sat_def_pct_indexM = 1000*(patch[0].sat_def_pct - patch[0].sat_def_pct_index*0.001);
            
            patch[0].sat_deficit_z = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].sat_def_z[patch[0].sat_def_pct_index];
        }else{
            // surface
            patch[0].sat_deficit_z = patch[0].sat_deficit;
            patch[0].sat_def_pct = 0.0;
            patch[0].sat_def_pct_index = 0;
            patch[0].sat_def_pct_indexM = 0;
        }
        
        // ----------------------------------------- finalize for printing output
        // fc & SatPct
        totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index];
        totalfc *= (1.0-patch[0].basementFrac); // <---- second thought on this, Oct 8, 2019; basement is 3m at most
        
        if (patch[0].sat_deficit < ZERO) {
            //patch[0].aboveWT_SatPct = 1.0;
            //patch[0].rootzone.SatPct = 1.0;
            patch[0].rootzone.field_capacity = 0.0;
            patch[0].field_capacity = 0.0;
        } else {
            patch[0].rootzone.field_capacity = totalfc * patch[0].zeroRootCoef * (patch[0].rootzone_scale_ref*patch[0].rootzone_end_reffc[patch[0].sat_def_pct_index] + (1.0-patch[0].rootzone_scale_ref)*patch[0].rootzone_start_reffc[patch[0].sat_def_pct_index]);
            patch[0].rootzone.field_capacity = min(patch[0].rootzone.field_capacity,patch[0].rootzone.potential_sat);
            patch[0].field_capacity = max(0.0,min(patch[0].sat_deficit-patch[0].rootzone.potential_sat, totalfc - patch[0].rootzone.field_capacity));
        }//if else
        
        if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity){
            printf("patch_daily_F(9): (%d,%d,%d) %lf %lf %lf %lf\n",
                   patch[0].ID, patch[0].soil_defaults[0][0].ID, patch[0].sat_def_pct_index,
                   patch[0].sat_deficit, patch[0].sat_deficit_z,
                   patch[0].rootzone.field_capacity, patch[0].field_capacity);
        }//debug
        
        // ----------------------------------------- aggregate streamflow from return and base flow
        if (patch[0].drainage_type == STREAM) {
           patch[0].streamflow += patch[0].return_flow
                   + patch[0].base_flow; //something is wrong here
        }
        
        
        
//        /* ******************************** this is done by each hour*/
//        patch[0].hourly_stream_flow += patch[0].hourly_subsur2stream_flow
//                      + patch[0].hourly_sur2stream_flow;
//
//        basin[0].basin_return_flow += (patch[0].return_flow) * patch[0].area;
//        /*--------------------------------------------------------------*/
//        /* final stream flow calculations                */
//        /*--------------------------------------------------------------*/
//
//        basin[0].basin_outflow += (patch[0].streamflow) * patch[0].area;
//        basin[0].basin_unsat_storage += patch[0].unsat_storage * patch[0].area;
//        basin[0].basin_sat_deficit += patch[0].sat_deficit * patch[0].area;
//        basin[0].basin_rz_storage += patch[0].rz_storage * patch[0].area;
//        basin[0].basin_detention_store += patch[0].detention_store
//                * patch[0].area;

        
        

    }// for loop i for patches
    
//    basin[0].basin_outflow /= basin_area;
//    basin[0].preday_basin_rz_storage /= basin_area;
//    basin[0].preday_basin_unsat_storage /= basin_area;
//    basin[0].preday_basin_detention_store /= basin_area;
//    basin[0].preday_basin_sat_deficit /= basin_area;
//    basin[0].basin_rz_storage /= basin_area;
//    basin[0].basin_unsat_storage /= basin_area;
//    basin[0].basin_detention_store /= basin_area;
//    basin[0].basin_sat_deficit /= basin_area;
//    water_balance = basin[0].preday_basin_rz_storage + basin[0].preday_basin_unsat_storage
//            + basin[0].preday_basin_detention_store - basin[0].preday_basin_sat_deficit
//            - (basin[0].basin_rz_storage + basin[0].basin_unsat_storage + basin[0].basin_detention_store
//                    - basin[0].basin_sat_deficit) - basin[0].basin_outflow;

    
	return;

} /*end compute_subsurface_routing.c*/

