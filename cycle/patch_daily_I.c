/*--------------------------------------------------------------*/
/* 								*/
/*		patch_daily_I					*/
/*								*/
/*	NAME							*/
/*	patch_daily_I 						*/
/*		 - performs cycling of patch state vars		*/
/*			at the START OF THE DAY			*/ 
/*								*/
/*	SYNOPSIS 						*/
/*	void patch_daily_I(					*/
/*			 struct patch_object	,		*/
/*			 struct command_line_object ,		*/
/*			 struct tec_entry,			*/
/*			 struct date)				*/
/*								*/
/*	OPTIONS							*/
/*	struct	world_object *world,				*/
/*	struct	basin_object *basin,				*/
/*	struct 	hillslope_object *hillslope,			*/
/*	struct  zone_object *zone,				*/
/*	struct patch_object *patch,				*/
/*	struct command_line_object *command_line,		*/
/*	struct command_line_object *command_line,		*/
/*	struct	tec_entry	*event,				*/
/*	struct	date current_date - local time (?)		*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	This routine performs simulation cycles on an identified*/
/*	canopy_stata in the patch.				*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*	March 12, 1997	C.Tague					*/
/*	- added calculation for patch effective lai		*/
/*								*/
/*								*/
/*	Sept 15 1997	RAF					*/
/*	Substantially modified accounting of current water	*/
/*	equivalent depth to sat zone and unsat_storage		*/
/*								*/
/*	Sept 29 1997 CT						*/
/*	switched above to an implementation using sat_deficit	*/
/*	as the TOP_model control volume - see discussion 	*/
/*	there							*/
/*								*/
/*	Oct 22 1997 CT						*/
/*								*/
/*	unsat storage now subtracted from sat_deficit after	*/
/*	return flow calculated - in previous version this	*/
/*	step was missing which  results in a 			*/
/*	serious over-estimation of sat_deficit after		*/
/*	return flow events					*/	
/*								*/
/*	Feb 2 1998 RAF						*/
/*	Included potential exfiltration module.			*/
/*								*/
/*	April 28 1998 RAF					*/
/*	Excluded stratum of height 0 from computation of	*/
/*	effective LAI of site.					*/
/*--------------------------------------------------------------*/
#include <stdlib.h>
#include "rhessys.h"
#include "functions.h"

void		patch_daily_I(
						  struct	world_object *world,
						  struct	basin_object *basin,
						  struct 	hillslope_object *hillslope,
						  struct  zone_object *zone,
						  struct patch_object *patch,
						  struct command_line_object *command_line,
						  struct	tec_entry	*event,
						  struct	date current_date)
{
	/*--------------------------------------------------------------*/
	/*  Local Function Declarations.                                */
	/*--------------------------------------------------------------*/
	void   canopy_stratum_daily_I(
		struct	world_object *,
		struct	basin_object *,
		struct 	hillslope_object *,
		struct  zone_object *,
		struct	patch_object *,
		struct canopy_strata_object *,
		struct command_line_object *,
		struct tec_entry *,
		struct date);
	
	double	compute_layer_field_capacity(
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
	
	
	double	compute_delta_water(
		int,
		double,
		double,
		double,
		double,
		double);

	double	compute_z_final(
		int,
		double,
		double,
		double,
		double,
		double);

	int	update_rootzone_moist(
		struct patch_object	*,
		struct	rooting_zone_object	*,
		struct command_line_object *);
	
	double	compute_capillary_rise(
		int,
		double,
		double,
		double,
		double,
		double);
	
	double  compute_soil_water_potential(
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
		double,
		double);
	
	int  compute_potential_decomp(
		double,
		double,
		double,
		double,
		double,
		struct  soil_c_object   *,
		struct  soil_n_object   *,
		struct  litter_c_object *,
		struct  litter_n_object *,
		struct  cdayflux_patch_struct *,
		struct  ndayflux_patch_struct *,
        struct    patch_object *,
        double, int, int);
	
	void    sort_patch_layers(struct patch_object *);

		
	void	update_litter_interception_capacity (double, 
		double,
		struct litter_c_object *,
		struct litter_object *);

	int	zero_patch_daily_flux(
		struct	patch_object *,
		struct  cdayflux_patch_struct *,
		struct  ndayflux_patch_struct *);
	
	
	long julday( struct date);
	/*--------------------------------------------------------------*/
	/*  Local variable definition.                                  */
	/*--------------------------------------------------------------*/
	int	layer, inx;
	int	stratum;
	double	cnt, count, theta;
    int vegtype;
	
	double  edible_leafc, grazing_mean_nc, grazing_Closs;
	struct  canopy_strata_object *strata;
	struct  dated_sequence	clim_event;

	/*--------------------------------------------------------------*/
	/*	zero out daily fluxes					*/
	/*--------------------------------------------------------------*/
	if (zero_patch_daily_flux(patch, &(patch[0].cdf), &(patch[0].ndf))){
		fprintf(stderr,
			"Error in zero_patch_daily_flux() from patch_daily_I.c... Exiting\n");
		exit(EXIT_FAILURE);
	}


	patch[0].precip_with_assim = 0.0;
    
    // -- update sat_def related variables
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
    
	/*-----------------------------------------------------*/
	/*  Compute potential saturation for rootzone layer   */
	/*-----------------------------------------------------*/
    if( patch[0].rootzone.depth > 0.0 && patch[0].rootzone.depth<0.01) patch[0].rootzone.depth= 0.01; // 10 mmm
    patch[0].rtz2_index = (int)(round(patch[0].rootzone.depth*1000));
    patch[0].rootzone.potential_sat = patch[0].soil_defaults[0][0].rtz2sat_def_0z[patch[0].rtz2_index];
    patch[0].rootdepth_index = patch[0].soil_defaults[0][0].rtz2sat_def_pct_index[patch[0].rtz2_index];
    //patch[0].rootdepth_indexM = patch[0].soil_defaults[0][0].rtz2sat_def_pct_indexM[patch[0].rtz2_index];
    //1000*(patch[0].rootzone.potential_sat*patch[0].soil_defaults[0][0].max_sat_def_1 - patch[0].rootdepth_index*0.001);
    if(patch[0].rootdepth_indexM<0.0 || patch[0].rootdepth_indexM>1.0)
        printf("rootdepth_indexM at (%d,%d,%f, %f %f)\n",patch[0].ID, patch[0].rootdepth_index, patch[0].rootdepth_indexM,
               patch[0].rootzone.potential_sat,patch[0].rootzone.depth);
    // how to handle no rootzone
    // quick solution is to set pointers to 002 and make an (double) zeroRoot[1 or 0]
    if(patch[0].rootzone.depth < 0.012){
        patch[0].rootzone_start_refcap = patch[0].soil_defaults[0][0].pot_caprise_00012r;
        patch[0].rootzone_end_refcap = patch[0].soil_defaults[0][0].pot_caprise_00012r;
        patch[0].rootzone_start_reffc = patch[0].soil_defaults[0][0].fc1_00012r;
        patch[0].rootzone_end_reffc = patch[0].soil_defaults[0][0].fc1_00012r;
        patch[0].rootzone_scale_ref = 1.0;// don't matter what value in this case
        patch[0].zeroRootCoef = 0.0; // treat it as no root; no real root depth be less than 1 cm.
    }else if(patch[0].rootzone.depth>0.012 && patch[0].rootzone.depth<=0.3){
        patch[0].rootzone_start_refcap = patch[0].soil_defaults[0][0].pot_caprise_00012r;
        patch[0].rootzone_end_refcap = patch[0].soil_defaults[0][0].pot_caprise_003r;
        patch[0].rootzone_start_reffc = patch[0].soil_defaults[0][0].fc1_00012r;
        patch[0].rootzone_end_reffc = patch[0].soil_defaults[0][0].fc1_003r;
        patch[0].rootzone_scale_ref = (patch[0].rootzone.depth-0.012)/(0.3-0.012);
        patch[0].zeroRootCoef = 1.0;
    }else if(patch[0].rootzone.depth>0.3 && patch[0].rootzone.depth<=0.6){
        patch[0].rootzone_start_refcap = patch[0].soil_defaults[0][0].pot_caprise_003r;
        patch[0].rootzone_end_refcap = patch[0].soil_defaults[0][0].pot_caprise_006r;
        patch[0].rootzone_start_reffc = patch[0].soil_defaults[0][0].fc1_003r;
        patch[0].rootzone_end_reffc = patch[0].soil_defaults[0][0].fc1_006r;
        patch[0].rootzone_scale_ref = (patch[0].rootzone.depth-0.3)/(0.6-0.3);
        patch[0].zeroRootCoef = 1.0;
    }else if(patch[0].rootzone.depth>0.6 && patch[0].rootzone.depth<=1.0){
        patch[0].rootzone_start_refcap = patch[0].soil_defaults[0][0].pot_caprise_006r;
        patch[0].rootzone_end_refcap = patch[0].soil_defaults[0][0].pot_caprise_010r;
        patch[0].rootzone_start_reffc = patch[0].soil_defaults[0][0].fc1_006r;
        patch[0].rootzone_end_reffc = patch[0].soil_defaults[0][0].fc1_010r;
        patch[0].rootzone_scale_ref = (patch[0].rootzone.depth-0.6)/(1.0-0.6);
        patch[0].zeroRootCoef = 1.0;
    }else if(patch[0].rootzone.depth>1.0 && patch[0].rootzone.depth<=1.5){
        patch[0].rootzone_start_refcap = patch[0].soil_defaults[0][0].pot_caprise_010r;
        patch[0].rootzone_end_refcap = patch[0].soil_defaults[0][0].pot_caprise_015r;
        patch[0].rootzone_start_reffc = patch[0].soil_defaults[0][0].fc1_010r;
        patch[0].rootzone_end_reffc = patch[0].soil_defaults[0][0].fc1_015r;
        patch[0].rootzone_scale_ref = (patch[0].rootzone.depth-1.0)/(1.5-1.0);
        patch[0].zeroRootCoef = 1.0;
    }else if(patch[0].rootzone.depth>1.5 && patch[0].rootzone.depth<=2.0){
        patch[0].rootzone_start_refcap = patch[0].soil_defaults[0][0].pot_caprise_015r;
        patch[0].rootzone_end_refcap = patch[0].soil_defaults[0][0].pot_caprise_020r;
        patch[0].rootzone_start_reffc = patch[0].soil_defaults[0][0].fc1_015r;
        patch[0].rootzone_end_reffc = patch[0].soil_defaults[0][0].fc1_020r;
        patch[0].rootzone_scale_ref = (patch[0].rootzone.depth-1.5)/(2.0-1.5);
        patch[0].zeroRootCoef = 1.0;
    }else if(patch[0].rootzone.depth>2.0 && patch[0].rootzone.depth<=2.5){
        patch[0].rootzone_start_refcap = patch[0].soil_defaults[0][0].pot_caprise_020r;
        patch[0].rootzone_end_refcap = patch[0].soil_defaults[0][0].pot_caprise_025r;
        patch[0].rootzone_start_reffc = patch[0].soil_defaults[0][0].fc1_020r;
        patch[0].rootzone_end_reffc = patch[0].soil_defaults[0][0].fc1_025r;
        patch[0].rootzone_scale_ref = (patch[0].rootzone.depth-2.0)/(2.5-2.0);
        patch[0].zeroRootCoef = 1.0;
    }else if(patch[0].rootzone.depth>2.5 && patch[0].rootzone.depth<=3.0){
        patch[0].rootzone_start_refcap = patch[0].soil_defaults[0][0].pot_caprise_025r;
        patch[0].rootzone_end_refcap = patch[0].soil_defaults[0][0].pot_caprise_030r;
        patch[0].rootzone_start_reffc = patch[0].soil_defaults[0][0].fc1_025r;
        patch[0].rootzone_end_reffc = patch[0].soil_defaults[0][0].fc1_030r;
        patch[0].rootzone_scale_ref = (patch[0].rootzone.depth-2.5)/(3.0-2.5);
        patch[0].zeroRootCoef = 1.0;
    }//if else
    
    
    
//	if (patch[0].rootzone.depth > ZERO)  {
//        // how is lulc frac affecting this part? Aug 8, 2019
//        // daily updated "patch[0].rootzone.depth" is done by this daily_I @ LINE 472 (below)
//        // stratum[0].rootzone.depth is first updated via "update_phenology", then aggregated to here
//
////        patch[0].rootzone.potential_sat = compute_delta_water(
////            command_line[0].verbose_flag,
////            patch[0].soil_defaults[0][0].porosity_0,
////            patch[0].soil_defaults[0][0].porosity_decay,
////            patch[0].soil_defaults[0][0].soil_depth,
////            patch[0].rootzone.depth, 0.0);
//        if (patch[0].sat_deficit_z > patch[0].rootzone.depth)
//            patch[0].rootzone.SatPct = min(patch[0].rz_storage/patch[0].rootzone.potential_sat/(1.0-patch[0].basementFrac), 1.0);
//        else
//            patch[0].rootzone.SatPct = min((patch[0].rz_storage + patch[0].rootzone.potential_sat - patch[0].sat_deficit)/patch[0].rootzone.potential_sat/(1.0-patch[0].basementFrac),1.0);
//
//	}else{
//        // no rootzone
//		patch[0].rootzone.potential_sat = 0.0;
//		patch[0].rootzone.SatPct = 0.0;
//    }
//
//	if (patch[0].sat_deficit < ZERO)
//		patch[0].aboveWT_SatPct = 1.0;
//	else
//		patch[0].aboveWT_SatPct = (patch[0].rz_storage+patch[0].unsat_storage)/patch[0].sat_deficit;
//
	/*--------------------------------------------------------------*/
	/*  compute standard deviation of theta based on soil parameters */
	/* assume no decay of porosity here 				*/
	/*--------------------------------------------------------------*/
//	theta = patch[0].aboveWT_SatPct * patch[0].soil_defaults[0][0].porosity_0;
//	patch[0].theta_std = (patch[0].soil_defaults[0][0].theta_mean_std_p2*theta*theta +
//				patch[0].soil_defaults[0][0].theta_mean_std_p1*theta);
    
    if(patch[0].rootzone.potential_sat>ZERO){
        if (patch[0].sat_deficit > patch[0].rootzone.potential_sat) theta = min(patch[0].rz_storage/patch[0].rootzone.potential_sat/(1.0-patch[0].basementFrac), 1.0);
        else theta = min((patch[0].rz_storage + patch[0].rootzone.potential_sat - patch[0].sat_deficit)/patch[0].rootzone.potential_sat/(1.0-patch[0].basementFrac),1.0);
    }else{ theta = 0.0; }
    //theta *= patch[0].soil_defaults[0][0].porosity_0; // really?
    patch[0].theta_std = (patch[0].soil_defaults[0][0].theta_mean_std_p2*theta*theta +
                patch[0].soil_defaults[0][0].theta_mean_std_p1*theta);
    
    
    
    //for biochemical mode accounting for basement;
    // constraintWaterTableTopDepth is patch-level average basement depth; basement depth is 3m --> basementFrac = constraintWaterTableTopDepth/3


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
        //(1*patch[0].rootzone_end_ref[187] + (1-1)*patch[0].rootzone_start_ref[187]);
        
        //patch[0].rootzone.field_capacity = totalfc * patch[0].zeroRootCoef * ( patch[0].rootdepth_indexM*(patch[0].rootzone_scale_ref*patch[0].rootzone_end_ref[patch[0].sat_def_pct_index+1] + (1.0-patch[0].rootzone_scale_ref)*patch[0].rootzone_start_ref[patch[0].sat_def_pct_index+1]) + (1.0-patch[0].rootdepth_indexM) * (patch[0].rootzone_scale_ref*patch[0].rootzone_end_ref[patch[0].sat_def_pct_index] + (1.0-patch[0].rootzone_scale_ref)*patch[0].rootzone_start_ref[patch[0].sat_def_pct_index]) );
        // rootzone.field_capacity is too little?
        patch[0].field_capacity = max(0.0,min(patch[0].sat_deficit-patch[0].rootzone.potential_sat, totalfc - patch[0].rootzone.field_capacity));
    }//if else
    
    if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity || (patch[0].sat_deficit>patch[0].rootzone.potential_sat && patch[0].field_capacity >= patch[0].sat_deficit-patch[0].rootzone.potential_sat+ZERO) || (patch[0].rootzone.field_capacity > patch[0].rootzone.potential_sat) ){
        printf("patch_daily_I(1): (%d,%d,%d) sat(%lf %lf) rtz(%lf %lf %lf %d,%lf) unsat(%lf %lf)\n",
               patch[0].ID, patch[0].soil_defaults[0][0].ID, patch[0].sat_def_pct_index,
               patch[0].sat_deficit, patch[0].sat_deficit_z,
               patch[0].rootzone.depth, patch[0].rootzone.potential_sat, patch[0].rootzone.field_capacity, patch[0].rootdepth_index,patch[0].rootzone_scale_ref,
               patch[0].sat_deficit-patch[0].rootzone.potential_sat, patch[0].field_capacity);
    }//debug


	/*--------------------------------------------------------------*/
	/*	Estimate potential cap rise.				*/
	/*	limited to water in sat zone.				*/
	/*--------------------------------------------------------------*/
    patch[0].potential_cap_rise = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].pot_caprise_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) *  patch[0].soil_defaults[0][0].pot_caprise_0z[patch[0].sat_def_pct_index];
    patch[0].potential_cap_rise = min(patch[0].potential_cap_rise, patch[0].available_soil_water);
    patch[0].potential_cap_rise = min(patch[0].potential_cap_rise, max(0.0,totalfc - patch[0].unsat_storage - patch[0].rz_storage));
    
//    patch[0].potential_cap_rise = compute_capillary_rise(
//		command_line[0].verbose_flag,
//		patch[0].sat_deficit_z,
//		patch[0].soil_defaults[0][0].psi_air_entry,
//		patch[0].soil_defaults[0][0].pore_size_index,
//		patch[0].soil_defaults[0][0].mz_v,
//		patch[0].soil_defaults[0][0].Ksat_0_v );
    
    
//    if(command_line[0].capreduction < 1.0) {patch[0].potential_cap_rise *=command_line[0].capreduction; }
//    if(command_line[0].capMax>0.0){
//        if(patch[0].sat_deficit_z<patch[0].newcapZ0){patch[0].potential_cap_rise = patch[0].newcapSlope*patch[0].sat_deficit_z + patch[0].newcapInter;}
//    }//pre-calculated in patch ini.
    
//if(patch[0].potential_cap_rise >= command_line[0].capMax){patch[0].potential_cap_rise *= command_line[0].capreduction;}
    
    



	/*--------------------------------------------------------------*/
	/*	Compute the max exfiltration rate.			*/
	/*								*/
	/*	First check if we are saturated.  If so the 		*/
	/*	potential exfiltration rate is the capillary rise rate.	*/
	/*								*/
	/*	If we are within the active unsat zone then we assume	*/
	/*	that the potential exfiltration rate is the minimum of	*/
	/*	the computed exfiltration rate and the potential cap	*/
	/*	rise - i.e. hydrologic connectivity between surface and	*/
	/*	water table.						*/
	/*								*/
	/*--------------------------------------------------------------*/
	if ( patch[0].sat_deficit_z <= patch[0].soil_defaults[0][0].psi_air_entry ){
		patch[0].potential_exfiltration = patch[0].potential_cap_rise;
	}else{
        patch[0].potential_exfiltration = compute_potential_exfiltration(
            command_line[0].verbose_flag,
            0.0,//patch[0].aboveWT_SatPct,
            patch[0].soil_defaults[0][0].exfiltration_wilting_point,
            patch[0].sat_def_pct_indexM *patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index+1] + (1.0 - patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index],
            patch[0].sat_def_pct_indexM *patch[0].soil_defaults[0][0].exfiltration_coef[patch[0].sat_def_pct_index+1] + (1.0 - patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].exfiltration_coef[patch[0].sat_def_pct_index],
            patch[0].rz_storage+patch[0].unsat_storage, //patch[0].soil_defaults[0][0].psi_air_entry,
            patch[0].sat_deficit, //patch[0].soil_defaults[0][0].pore_size_index,
            patch[0].sat_def_pct_indexM *patch[0].soil_defaults[0][0].sat_def_0zm[patch[0].sat_def_pct_index+1] + (1.0 - patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].sat_def_0zm[patch[0].sat_def_pct_index],
            patch[0].soil_defaults[0][0].exfiltration_S_pow);
	}//if else
    // debug
    if( patch[0].potential_exfiltration<0 || patch[0].potential_cap_rise<0){
        printf("patch dailyI patch[0].potential_exfiltration & patch[0].potential_cap_rise(%d) %f %f %f\n",
               patch[0].ID, patch[0].potential_exfiltration, patch[0].potential_cap_rise, patch[0].sat_deficit_z);
    }// if debug


	/*-----------------------------------------------------*/
	/* 	Check for any grazing activity from a land use default file			*/
	/*-----------------------------------------------------*/
	if (patch[0].base_stations != NULL) {
		inx = patch[0].base_stations[0][0].dated_input[0].grazing_Closs.inx;
		if (inx > -999) {
			clim_event = patch[0].base_stations[0][0].dated_input[0].grazing_Closs.seq[inx];
			while (julday(clim_event.edate) < julday(current_date)) {
				patch[0].base_stations[0][0].dated_input[0].grazing_Closs.inx += 1;
				inx = patch[0].base_stations[0][0].dated_input[0].grazing_Closs.inx;
				clim_event = patch[0].base_stations[0][0].dated_input[0].grazing_Closs.seq[inx];
				}
			if ((clim_event.edate.year != 0) && ( julday(clim_event.edate) == julday(current_date)) ) {
				grazing_Closs = clim_event.value;
				}
			else grazing_Closs = 0.0;
			} 
		else grazing_Closs = patch[0].landuse_defaults[0][0].grazing_Closs;
		}
	else grazing_Closs = patch[0].landuse_defaults[0][0].grazing_Closs;
	patch[0].grazing_Closs = grazing_Closs;

	/*--------------------------------------------------------------*/
	/*	Cycle through the canopy layers.			*/
	/*--------------------------------------------------------------*/
	edible_leafc = 0.0;
	grazing_mean_nc = 0.0;
	cnt = 0;
    vegtype=0;
	for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
		/*--------------------------------------------------------------*/
		/*	Cycle through the canopy strata				*/
		/*--------------------------------------------------------------*/
		for ( stratum=0 ; stratum<patch[0].layers[layer].count; stratum++ ){

			strata = patch[0].canopy_strata[(patch[0].layers[layer].strata[stratum])];
            
			patch[0].preday_rain_stored += strata->cover_fraction * strata->rain_stored;
			patch[0].preday_snow_stored += strata->cover_fraction * strata->snow_stored;
			if ((strata[0].defaults[0][0].epc.edible == 1) && (strata[0].cs.leafc > ZERO)) {
				edible_leafc += strata->cs.leafc * strata->cover_fraction;
				cnt += 1;
				grazing_mean_nc += strata->ns.leafn/strata->cs.leafc * strata->cover_fraction;
            }//edible
            
            if ( strata[0].defaults[0][0].epc.veg_type != NON_VEG ){
                vegtype = 1;
            }
            
			canopy_stratum_daily_I(
				world,
				basin,
				hillslope,
				zone,
				patch,
				patch[0].canopy_strata[(patch[0].layers[layer].strata[stratum])],
				command_line,
				event,
				current_date );
		}//stratum
	}//layer
	patch[0].grazing_Closs = min(edible_leafc, patch[0].grazing_Closs);
	if (cnt > 0)
		patch[0].grazing_mean_nc = grazing_mean_nc / cnt;

	/*--------------------------------------------------------------*/
	/*	Calculate effective patch lai from stratum					*/
	/*	- for later use by zone_daily_F								*/
	/*      Accumulate root biomass for patch soil -		*/
	/*      required for N updake from soil                         */
	/*	also determine total plant carbon			*/
	/*	- if grow option is specified				*/
	/*--------------------------------------------------------------*/
	patch[0].effective_lai = 0.0;
	patch[0].soil_cs.frootc = 0.0;
	patch[0].rootzone.depth = 0.0;
    count = 0.0;
	for ( stratum=0 ; stratum<patch[0].num_canopy_strata; stratum++){
		patch[0].effective_lai += patch[0].canopy_strata[stratum][0].epv.proj_lai;
		if (command_line[0].grow_flag > 0) {
			patch[0].soil_cs.frootc
				+= patch[0].canopy_strata[stratum][0].cover_fraction
				* patch[0].canopy_strata[stratum][0].cs.frootc;
			patch[0].preday_totalc
				+= patch[0].canopy_strata[stratum][0].cover_fraction
				* patch[0].canopy_strata[stratum][0].cs.preday_totalc;
			patch[0].preday_totaln
				+= patch[0].canopy_strata[stratum][0].cover_fraction
				* patch[0].canopy_strata[stratum][0].ns.preday_totaln;
		}//if
		patch[0].rootzone.depth = max(patch[0].rootzone.depth, 
			 patch[0].canopy_strata[stratum][0].rootzone.depth);
        
	}//stratum
	patch[0].effective_lai = patch[0].effective_lai / patch[0].num_canopy_strata;
	/*--------------------------------------------------------------*/
	/*	re-sort patch layers to account for any changes in 	*/
	/*	height							*/
	/*------------------------------------------------------------------------*/
	sort_patch_layers(patch);


	/*------------------------------------------------------------------------*/
	/*	compute current soil moisture potential					*/
	/*	do this before nitrogen updake occurs later in the day			*/
	/*------------------------------------------------------------------------*/
	patch[0].psi = compute_soil_water_potential(
		command_line[0].verbose_flag,
		patch[0].soil_defaults[0][0].theta_psi_curve,
		patch[0].Tsoil,
		-1.0*patch[0].soil_defaults[0][0].psi_max,
		-10.0,
		patch[0].soil_defaults[0][0].psi_air_entry,
		patch[0].soil_defaults[0][0].pore_size_index,
		patch[0].soil_defaults[0][0].p3,
		patch[0].soil_defaults[0][0].p4,
		patch[0].rz_storage,//patch[0].soil_defaults[0][0].porosity_0,
		patch[0].rootzone.potential_sat,//patch[0].soil_defaults[0][0].porosity_decay,
		patch[0].sat_deficit //patch[0].rootzone.SatPct //patch[0].aboveWT_SatPct // which one?
        );


	if (command_line[0].grow_flag > 0) {

		/*--------------------------------------------------------------*/
		/*	update litter interception capacity			*/
		/*--------------------------------------------------------------*/
		update_litter_interception_capacity(
			patch[0].litter.moist_coef,
			patch[0].litter.density,
			&(patch[0].litter_cs),
			&(patch[0].litter));

	
        if (compute_potential_decomp(
			patch[0].Tsoil,
			patch[0].soil_defaults[0][0].psi_max,
			patch[0].soil_defaults[0][0].psi_air_entry,
			theta,///
			patch[0].theta_std,
			&(patch[0].soil_cs),
			&(patch[0].soil_ns),
			&(patch[0].litter_cs),
			&(patch[0].litter_ns),
			&(patch[0].cdf),
			&(patch[0].ndf),
            patch,
            command_line[0].soilDecayScalar,
            command_line[0].soilCNadaptation_flag,
            vegtype
			) != 0){
			fprintf(stderr,"fATAL ERROR: in compute_potential_decomp() ... Exiting\n");
			exit(EXIT_FAILURE);
        }//if
	}//grow_flag

		/*--------------------------------------------------------------*/
		/*	zeros the accumulative rain_throughfall for 24 hours	*/
		/*--------------------------------------------------------------*/
	patch[0].rain_throughfall_24hours=0.0;
	patch[0].recharge=0;
	patch[0].rz_drainage=0;
	patch[0].unsat_drainage=0;

    
    
	return;
}/*end patch_daily_I.c*/
