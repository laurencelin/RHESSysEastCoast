/*--------------------------------------------------------------*/
/* 																*/
/*						patch_hourly							*/
/*																*/
/*	NAME														*/
/*	patch_hourly 												*/
/*				 - performs cycling and output of a patch		*/
/*																*/
/*																*/
/*	SYNOPSIS													*/
/*	void patch_hourly(											*/
/*						struct	world_object *,   				*/	
/*						struct	basin_object *,   				*/	
/*						struct	hillslope_object *,				*/	
/*						struct	zone_object *,   				*/	
/*						struct 	patch_object *,					*/
/*						struct 	command_line_object *,			*/
/*						struct  tec_entry   *,					*/
/*						struct  date );							*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	This routine performs simulation cycles on an identified	*/
/*	canopy_stata in the patch.									*/ 
/*																*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*																*/
/*																*/
/*--------------------------------------------------------------*/
#include "rhessys.h"

void		patch_hourly(
						 struct world_object *world,
						 struct basin_object *basin,
						 struct hillslope_object *hillslope,
						 struct zone_object *zone,
						 struct patch_object *patch,
						 struct command_line_object *command_line,
						 struct	tec_entry	*event,
						 struct	date current_date)
{
    //printf("now you are in a patch hourly %d\n",patch->ID);
	/*--------------------------------------------------------------*/
	/*  Local Function Declarations.                                */
	/*--------------------------------------------------------------*/
	void   canopy_stratum_hourly (
		struct world_object *,
		struct basin_object *,
		struct hillslope_object *,
		struct zone_object *,
		struct patch_object *,
		struct canopy_strata_object *,
		struct command_line_object *,
		struct tec_entry *,
		struct date);
	
	double	compute_delta_water(
		int,
		double,
		double,
		double,
		double,
		double);

	double	compute_infiltration(
		int,
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

	double compute_layer_field_capacity(
		int,
		int,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double);
	
	double  compute_unsat_zone_drainage(
		int,
		int,
		double,
		double,
		double,
		double,
		double,
		double, double);
 
	double  compute_z_final(
		int,
		double,
		double,
		double,
		double,
		double);

	void 	surface_hourly(
		struct world_object *,
		struct basin_object *,
		struct hillslope_object *,
		struct zone_object *,
		struct patch_object *,
		struct command_line_object *,
		struct tec_entry *,
		struct date);

	void	update_soil_moisture(
		int	verbose_flag,
		double	infiltration,
		double	net_inflow,
		struct	patch_object	*patch,
		struct 	command_line_object *command_line,
		struct	date 			current_date);

	int 	update_gw_drainage(
			struct patch_object *,
			struct hillslope_object *,
			struct zone_object *,
			struct command_line_object *,
			struct date); 
	/*--------------------------------------------------------------*/
	/*	Local Variable Declarations.								*/
	/*--------------------------------------------------------------*/
	int	stratum;
	int	layer;
	double  net_inflow, duration, infiltration;
	double 	rz_drainage, unsat_drainage;
	double  theta;
	struct 	litter_object *litter;
	/*--------------------------------------------------------------*/
	/*	process any hourly rainfall				*/
	/*--------------------------------------------------------------*/

	if ( zone[0].hourly_rain_flag == 1) {
		patch[0].hourly[0].rain_throughfall = zone[0].hourly[0].rain;
		patch[0].precip_with_assim += zone[0].hourly[0].rain;
		}
	else
		patch[0].hourly[0].rain_throughfall = 0.0;

	patch[0].hourly[0].NO3_throughfall = (1-command_line[0].fracDirectNdep) * zone[0].ndep_NO3/24.0;// in stdzone file, the n_deposition is in kg/m2/yr

	/*--------------------------------------------------------------*/
	/*	Cycle through the canopy strata								*/
	/*	above the snowpack					*/
	/*--------------------------------------------------------------*/
	for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
		if ( (patch[0].layers[layer].height > patch[0].snowpack.height) ){
			patch[0].rain_throughfall_final = 0.0;
			/* NO3_throughfall_final collects NO3_throughfall first from null_cover area */
			/* Then use it to sum up the NO3_thoughfall from canopy cover */
	    		patch[0].hourly[0].NO3_throughfall_final = patch[0].layers[layer].null_cover * patch[0].hourly[0].NO3_throughfall;
			for (stratum=0 ;stratum<patch[0].layers[layer].count; stratum++ ){
				canopy_stratum_hourly(
					world,
					basin,
					hillslope,
					zone,
					patch,
					patch[0].canopy_strata[stratum],
					command_line,
					event,
					current_date );
			}
			/*--------------------------------------------------------------*/
			/*	process any hourly throughfallthat falls on a snowpack */
			/*--------------------------------------------------------------*/
			patch[0].hourly[0].rain_throughfall = patch[0].rain_throughfall_final;
			patch[0].hourly[0].NO3_throughfall = patch[0].hourly[0].NO3_throughfall_final;

		}
	}
      	
	if (patch[0].snowpack.water_equivalent_depth > 0.0) {
		patch[0].snowpack.water_equivalent_depth
			+= patch[0].hourly[0].rain_throughfall;
		patch[0].hourly[0].rain_throughfall = 0.0;
	}
	/*--------------------------------------------------------------*/
	/*	Cycle through the canopy strata				*/
	/*	below the snowpack					*/
	/*--------------------------------------------------------------*/
	for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
		if ( (patch[0].layers[layer].height <= patch[0].snowpack.height) ){
			patch[0].rain_throughfall_final = 0.0;
			patch[0].hourly[0].NO3_throughfall_final = patch[0].layers[layer].null_cover * patch[0].hourly[0].NO3_throughfall;

			for ( stratum=0;stratum<patch[0].layers[layer].count; stratum++ ){
				canopy_stratum_hourly(
					world,
					basin,
					hillslope,
					zone,
					patch,
					patch[0].canopy_strata[stratum],
					command_line,
					event,
					current_date );
			}
		patch[0].hourly[0].rain_throughfall = patch[0].rain_throughfall_final;
		patch[0].hourly[0].NO3_throughfall = patch[0].hourly[0].NO3_throughfall_final;	
		}

	}


	patch[0].surface_NO3 += patch[0].hourly[0].NO3_throughfall;// on a rain day, n does not add to surfaceN in daily step because @115 rain_throughfall = 0.0;
    if ( zone[0].hourly_rain_flag == 1) {patch[0].surface_NO3 += (command_line[0].fracDirectNdep) * zone[0].ndep_NO3/24.0;}
	patch[0].detention_store += patch[0].hourly[0].rain_throughfall;	

	/*--------------------------------------------------------------*/
	/*	include any detention storage as throughfall		*/
	/*--------------------------------------------------------------*/
	if (zone[0].hourly_rain_flag == 1) {
		/*--------------------------------------------------------------*/
		/*	calculate the litter interception			*/
		/*--------------------------------------------------------------*/
		surface_hourly(
					world,
					basin,
					hillslope,
					zone,
					patch,
					command_line,
					event,
					current_date);

        /*--------------------------------------------------------------*/
        /*     Above ground Hydrologic Processes            */
        /*     compute infiltration into the soil            */
        /*    from snowmelt or rain_throughfall            */
        /*    for now assume that all water infilatrates        */
        /*--------------------------------------------------------------*/
        if (patch[0].detention_store > 0.0) {
            /*------------------------------------------------------------------------*/
            /*    drainage to a deeper groundwater store                  */
            /*    move both nitrogen and water                           */
            /*------------------------------------------------------------------------*/
            if (command_line[0].gw_flag > 0 ){

                // it constraints itself to hourly by reading the hourly data.
                if ( update_gw_drainage(patch,
                    hillslope,
                    zone,
                    command_line,
                    current_date) != 0) {
                    fprintf(stderr,"fATAL ERROR: in update_decomp() ... Exiting\n");
                    exit(EXIT_FAILURE);
                }
            }
	  
        
            net_inflow=patch[0].detention_store;
            /*--------------------------------------------------------------*/
            /*      - if rain duration is zero, then input is from snow     */
            /*      melt  assume full daytime duration                      */
            /*--------------------------------------------------------------*/
            if (zone[0].hourly[0].rain_duration <= ZERO)
                duration = 60*60/(86400);
            else
                duration = zone[0].hourly[0].rain_duration/(86400);
            
            
            infiltration = compute_infiltration(
                command_line[0].verbose_flag,
                patch[0].sat_deficit_z,
                0.0,//patch[0].aboveWT_SatPct, // initiated in daily_I()
                patch[0].Ksat_vertical, // 1- impervious
                patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index],
                patch[0].rz_storage+patch[0].unsat_storage,
                patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].sat_def_0zm[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].sat_def_0zm[patch[0].sat_def_pct_index],
                patch[0].sat_deficit,
                net_inflow,
                duration,// constrain the infiltration by hour
                patch[0].soil_defaults[0][0].psi_air_entry);
            } else infiltration = 0.0;

        if (infiltration < 0.0)
            printf("\nInfiltration %lf < 0 for %d on %d",
                infiltration,
                patch[0].ID, current_date.day);
        /*--------------------------------------------------------------*/
        /* determine fate of hold infiltration excess in detention store */
        /* infiltration excess will removed during routing portion    */
        /*--------------------------------------------------------------*/
        infiltration=min(infiltration,patch[0].detention_store);

        patch[0].detention_store -= infiltration;
        
        if (infiltration>ZERO) {
            /*--------------------------------------------------------------*/
            /*    Update patch level soil moisture with final infiltration.    */
            /*--------------------------------------------------------------*/
            update_soil_moisture(
                command_line[0].verbose_flag,
                infiltration,
                net_inflow,
                patch,
                command_line,
                current_date );
        } /* end if infiltration > ZERO */

        /* aggregate the hourly recharge */
        patch[0].recharge += infiltration;
	
	} /* end if rain throughfall */
	/*--------------------------------------------------------------*/
	/*	Destroy the patch hourly object.							*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/*	use rain_throughfall_24hours to collect the accumulative rain_throughfall		*/
	/*--------------------------------------------------------------*/

	patch[0].rain_throughfall_24hours+=patch[0].hourly[0].rain_throughfall;


	/*-------------------------------------------------------------------------*/
	/*	Compute current actual depth to water table				*/
	/*------------------------------------------------------------------------*/
    // reason to update "sat_deficit" here because
    //canopy_stratum_hourly -- throughfall and no phenology call --> no change to rootzone.depth
    //surface_hourly -- detention and litter
    //update_gw_drainage -- detention and gw
    //compute_infiltration -- calculate infiltration
    //update_soil_moisture -- does the infiltration add water to soil and modifies "sat_deficit" 

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
//	patch[0].sat_deficit_z = compute_z_final(
//		command_line[0].verbose_flag,
//		patch[0].soil_defaults[0][0].porosity_0,
//		patch[0].soil_defaults[0][0].porosity_decay,
//		patch[0].soil_defaults[0][0].soil_depth,
//		0.0,
//		-1.0 * patch[0].sat_deficit);

	/*--------------------------------------------------------------*/
	/*	compute new field capacity				*/
	/*--------------------------------------------------------------*/


    double totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index];
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

	/*-------------------------------------------------------------------------*/
	/*	Compute current actual depth to water table				*/
	/*------------------------------------------------------------------------*/
    // just updated "sat_deficit_z" above before FC, and FC does not move water.
//	patch[0].sat_deficit_z = compute_z_final(
//		command_line[0].verbose_flag,
//		patch[0].soil_defaults[0][0].porosity_0,
//		patch[0].soil_defaults[0][0].porosity_decay,
//		patch[0].soil_defaults[0][0].soil_depth,
//		0.0,
//		-1.0 * patch[0].sat_deficit);


	/*--------------------------------------------------------------*/
	/*      Recompute patch soil moisture storage                   */
	/*--------------------------------------------------------------*/
	if (patch[0].sat_deficit < ZERO) {
		//patch[0].aboveWT_SatPct = 1.0;
		//patch[0].rootzone.SatPct = 1.0;
		rz_drainage = 0.0;
		unsat_drainage = 0.0;
	} else {
        /* Constant vertical profile of soil porosity */
		/*-------------------------------------------------------*/
		/*	soil drainage and storage update	     	 */
		/*-------------------------------------------------------*/
		
		//patch[0].rootzone.SatPct = min(patch[0].rz_storage/patch[0].rootzone.potential_sat/(1.0-patch[0].basementFrac), 1.0);
        //patch[0].aboveWT_SatPct = (patch[0].rz_storage+patch[0].unsat_storage)/patch[0].sat_deficit;
        
		rz_drainage = compute_unsat_zone_drainage(
			command_line[0].verbose_flag,
			patch[0].soil_defaults[0][0].theta_psi_curve,
			patch[0].soil_defaults[0][0].pore_size_index,
			patch[0].rootzone.potential_sat, //patch[0].rootzone.SatPct,
			0.04166667*(patch[0].rootdepth_indexM * patch[0].soil_defaults[0][0].vksat_0zm[patch[0].rootdepth_index+1] + (1.0-patch[0].rootdepth_indexM)* patch[0].soil_defaults[0][0].vksat_0zm[patch[0].rootdepth_index]),
			0.04166667*(patch[0].rootdepth_indexM * patch[0].soil_defaults[0][0].vksat_z[patch[0].rootdepth_index+1] + (1.0-patch[0].rootdepth_indexM) * patch[0].soil_defaults[0][0].vksat_z[patch[0].rootdepth_index]),
			patch[0].rz_storage,
            patch[0].rootzone.field_capacity,
            patch[0].sat_deficit);
	
		patch[0].rz_storage -=  rz_drainage;
		patch[0].unsat_storage +=  rz_drainage;
		
		unsat_drainage = compute_unsat_zone_drainage(
			command_line[0].verbose_flag,
			patch[0].soil_defaults[0][0].theta_psi_curve,
			patch[0].soil_defaults[0][0].pore_size_index,
			patch[0].sat_deficit - patch[0].rootzone.potential_sat, //patch[0].aboveWT_SatPct,
			0.04166667*(patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM)* patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index]),
			0.04166667*(patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].vksat_z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].vksat_z[patch[0].sat_def_pct_index]),
			patch[0].unsat_storage,
            patch[0].field_capacity,
            patch[0].sat_deficit);
	
		patch[0].unsat_storage -=  unsat_drainage;
		patch[0].sat_deficit -=  unsat_drainage;
        
	}//if else
	
	patch[0].unsat_drainage += unsat_drainage;
	patch[0].rz_drainage += rz_drainage;
	patch[0].hourly_unsat_drainage = unsat_drainage;
	patch[0].hourly_rz_drainage = rz_drainage;
	
	/* ---------------------------------------------- */
	/*     Final rootzone saturation calculation      */
	/* ---------------------------------------------- */
	
// no change to the rootzone.depth here
//	/*-----------------------------------------------------*/
//	/*  re-Compute potential saturation for rootzone layer   */
//	/*-----------------------------------------------------*/
//	if (patch[0].rootzone.depth > ZERO)
//		patch[0].rootzone.potential_sat = compute_delta_water(
//		command_line[0].verbose_flag,
//		patch[0].soil_defaults[0][0].porosity_0,
//		patch[0].soil_defaults[0][0].porosity_decay,
//		patch[0].soil_defaults[0][0].soil_depth,
//		patch[0].rootzone.depth, 0.0);

	/*------------------------------------------------------------------------*/
	/*	Compute current actual depth to water table				*/
	/*------------------------------------------------------------------------*/
    // above drainage makes changes to the sat_deficit;

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
//	patch[0].sat_deficit_z = compute_z_final(
//		command_line[0].verbose_flag,
//		patch[0].soil_defaults[0][0].porosity_0,
//		patch[0].soil_defaults[0][0].porosity_decay,
//		patch[0].soil_defaults[0][0].soil_depth,
//		0.0,
//		-1.0 * patch[0].sat_deficit);

    
    if(patch[0].rootzone.potential_sat>ZERO){
        if (patch[0].sat_deficit > patch[0].rootzone.potential_sat) theta = min(patch[0].rz_storage/patch[0].rootzone.potential_sat/(1.0-patch[0].basementFrac), 1.0);
        else theta = min((patch[0].rz_storage + patch[0].rootzone.potential_sat - patch[0].sat_deficit)/patch[0].rootzone.potential_sat/(1.0-patch[0].basementFrac),1.0);
    }else{ theta = 0.0; }
	patch[0].theta_std = (patch[0].soil_defaults[0][0].theta_mean_std_p2*theta*theta +
				patch[0].soil_defaults[0][0].theta_mean_std_p1*theta);
	


} /*end patch_hourly.c*/
