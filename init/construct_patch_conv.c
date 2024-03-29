/*--------------------------------------------------------------*/
/* 																*/
/*					construct_patch								*/
/*																*/
/*	construct_patch.c - creates a patch object					*/
/*																*/
/*	NAME														*/
/*	construct_patch.c - creates a patch object					*/
/*																*/
/*	SYNOPSIS					 								*/
/*	struct patch_object  construct_patch( 						*/
/*					FILE	*world_file,						*/
/*					struct	command_line_object	*command_line,	*/
/*					struct	default_object)						*/
/* 																*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	Allocates memory and reads in parameters from the hillslope	*/
/*	file to create a patch object.  Invokes construction		*/
/*	of canopy_stratum objects.									*/
/*																*/
/*	Refer to cnostruct_basin.c for a specification of the		*/
/*	hillslopes file.											*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	We assume that the FILE pointers to the 					*/
/*	hillslope file are positioned properly.						*/
/*	 															*/
/*	We assume that the basin and hillslope files have correct	*/
/*	syntax.														*/
/*																*/
/*	Original code, January 16, 1996.							*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"
#include "phys_constants.h"

struct patch_object *construct_patch(
									 struct	command_line_object	*command_line,
									 FILE	*world_file,
									 int     num_world_base_stations,
									 struct  base_station_object **world_base_stations,
									 struct	default_object	*defaults)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	struct base_station_object *assign_base_station(
		int ,
		int ,
		struct base_station_object **);
	struct 	canopy_strata_object *construct_canopy_strata(
		struct command_line_object *,
		FILE	*,
		struct	patch_object *,
		int     num_world_base_stations,
		struct  base_station_object **world_base_stations,
		struct	default_object	*defaults);
	double	compute_z_final( 	int,
		double,
		double,
		double,
		double,
		double);
	void	update_litter_interception_capacity (
		double, 
		struct litter_c_object *,
		struct litter_object *);
	
	void	sort_patch_layers(struct patch_object *);
	void	*alloc(	size_t, char *, char *);
	
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int		base_stationID;
	int		i;
	int		soil_default_object_ID;
	int		landuse_default_object_ID;
	char		record[MAXSTR];
	struct patch_object *patch;
	double	mpar;
	
	/*--------------------------------------------------------------*/
	/*  Allocate a patch object.                                */
	/*--------------------------------------------------------------*/
	patch = (struct patch_object *) alloc( 1 *
		sizeof( struct patch_object ),"patch","construct_patch");
	
	/*--------------------------------------------------------------*/
	/*	Read in the next patch record for this hillslope.			*/
	/*--------------------------------------------------------------*/
	fscanf(world_file,"%d",&(patch[0].ID));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].x));

	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].y));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].z));
	read_record(world_file, record);
	fscanf(world_file,"%d",&(soil_default_object_ID));
	read_record(world_file, record);
	fscanf(world_file,"%d",&(landuse_default_object_ID));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].area));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].slope));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].lna));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].Ksat_vertical));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(mpar));
	read_record(world_file, record);
	if (command_line[0].std_flag == 1) {
		fscanf(world_file,"%lf",&(patch[0].std));
		read_record(world_file, record);
		}
	else patch[0].std = 0.0;
	fscanf(world_file,"%lf",&(patch[0].unsat_storage));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].sat_deficit));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].snowpack.water_equivalent_depth));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].snowpack.water_depth));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].snowpack.T));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].snowpack.surface_age));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].snowpack.energy_deficit));
	read_record(world_file, record);
	if (command_line[0].snow_scale_flag == 1) {
		fscanf(world_file,"%lf",&(patch[0].snow_redist_scale));
		read_record(world_file, record);
		}

	patch[0].slope = patch[0].slope * DtoR;
	patch[0].surface_Tday = -999.9;
	patch[0].surface_Tnight = -999.9;
	patch[0].preday_sat_deficit  = patch[0].sat_deficit;
	patch[0].deltaS = 0.0;
	patch[0].streamflow = 0.0;
	patch[0].return_flow = 0.0;
	patch[0].gw_drainage = 0.0;
	patch[0].infiltration_excess = 0.0;
	patch[0].streamflow_N = 0.0;
	patch[0].snowpack.height = patch[0].snowpack.water_equivalent_depth *10.0;
	patch[0].tmp = 0.0;
	patch[0].detention_store = 0.0;	
	patch[0].rz_storage = 0.2;
	
	/*--------------------------------------------------------------*/
	/*      initialize accumulator variables for this patch         */
	/*--------------------------------------------------------------*/
	patch[0].acc_month.snowpack = 0.0;
	patch[0].acc_month.et = 0.0;
	patch[0].acc_month.psn = 0.0;
	patch[0].acc_month.denitrif = 0.0;
	patch[0].acc_month.sm_deficit = 0.0;
	patch[0].acc_month.DOC_loss = 0.0;
	patch[0].acc_month.DON_loss = 0.0;
	patch[0].acc_month.leach = 0.0;	
	patch[0].acc_month.length = 0;
	patch[0].acc_year.length = 0;
	patch[0].acc_year.num_threshold = 0;
	patch[0].acc_year.et = 0.0;	
	patch[0].acc_year.psn = 0.0;	
	patch[0].acc_year.denitrif = 0.0;	
	patch[0].acc_year.leach = 0.0;	
	
	/*--------------------------------------------------------------*/
	/*	Variables for the dynamic version are included here     */
	/*--------------------------------------------------------------*/
	fscanf(world_file,"%lf",&(patch[0].litter.rain_stored));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].litter_cs.litr1c));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].litter_ns.litr1n));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].litter_cs.litr2c));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].litter_cs.litr3c));
	read_record(world_file, record); 
	fscanf(world_file,"%lf",&(patch[0].litter_cs.litr4c));
	read_record(world_file, record);
	
	patch[0].litter_ns.litr2n = patch[0].litter_cs.litr2c / CEL_CN;
	patch[0].litter_ns.litr3n = patch[0].litter_cs.litr3c / CEL_CN;
	patch[0].litter_ns.litr4n = patch[0].litter_cs.litr4c / LIG_CN;
	
	fscanf(world_file,"%lf",&(patch[0].soil_cs.soil1c));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].soil_ns.sminn));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].soil_ns.nitrate));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].soil_cs.soil2c));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].soil_cs.soil3c));
	read_record(world_file, record);
	fscanf(world_file,"%lf",&(patch[0].soil_cs.soil4c));
	read_record(world_file, record);


	patch[0].soil_ns.soil1n = patch[0].soil_cs.soil1c / SOIL1_CN;
	patch[0].soil_ns.soil2n = patch[0].soil_cs.soil2c / SOIL2_CN;
	patch[0].soil_ns.soil3n = patch[0].soil_cs.soil3c / SOIL3_CN;
	patch[0].soil_ns.soil4n = patch[0].soil_cs.soil4c / SOIL4_CN;
	
	/*--------------------------------------------------------------*/
	/*	initialize sinks					*/
	/*--------------------------------------------------------------*/
	
	patch[0].litter_cs.litr1c_hr_snk = 0.0;
	patch[0].litter_cs.litr2c_hr_snk = 0.0;
	patch[0].litter_cs.litr4c_hr_snk = 0.0;
	
	patch[0].soil_cs.soil1c_hr_snk = 0.0;
	patch[0].soil_cs.soil2c_hr_snk = 0.0;
	patch[0].soil_cs.soil4c_hr_snk = 0.0;
	
	patch[0].soil_ns.nfix_src = 0.0;
	patch[0].soil_ns.ndep_src = 0.0;
	patch[0].soil_ns.nleached_snk = 0.0;
	patch[0].soil_ns.nvolatilized_snk = 0.0;

	patch[0].surface_NO3 = 0.0;
	patch[0].surface_NH4 = 0.0;
	patch[0].fertilizer_NO3 = 0.0;
	patch[0].fertilizer_NH4 = 0.0;
	
	/*--------------------------------------------------------------*/
	/*	Assign	defaults for this patch								*/
	/*--------------------------------------------------------------*/
	patch[0].soil_defaults = (struct soil_default **)
		alloc( sizeof(struct soil_default *),"defaults",
		"construct_patch" );
	
	i = 0;
	while (defaults[0].soil[i].ID != soil_default_object_ID) {
		i++;
		/*--------------------------------------------------------------*/
		/*  Report an error if no match was found.  Otherwise assign    */
		/*  the default to point to this patch.						    */
		/*--------------------------------------------------------------*/
		if ( i>= defaults[0].num_soil_default_files ){
			fprintf(stderr,
				"\nFATAL ERROR: in construct_patch, soil default ID %d not found for patch %d\n" ,
				soil_default_object_ID, patch[0].ID);
			exit(EXIT_FAILURE);
		}
	} /* end-while */
	patch[0].soil_defaults[0] = &defaults[0].soil[i];

	patch[0].landuse_defaults = (struct landuse_default **)
		alloc( sizeof(struct landuse_default *),"defaults",
		"construct_patch" );
	
	i = 0;
	while (defaults[0].landuse[i].ID != landuse_default_object_ID) {
		i++;
		/*--------------------------------------------------------------*/
		/*  Report an error if no match was found.  Otherwise assign    */
		/*  the default to point to this patch.						    */
		/*--------------------------------------------------------------*/
		if ( i>= defaults[0].num_landuse_default_files ){
			fprintf(stderr,
				"\nFATAL ERROR: in construct_patch, landuse default ID %d not found for patch %d\n" ,
				landuse_default_object_ID, patch[0].ID);
			exit(EXIT_FAILURE);
		}
	} /* end-while */
	patch[0].landuse_defaults[0] = &defaults[0].landuse[i];
	/*--------------------------------------------------------------*/
	/* FOR NOW, assign DON loss rate from command line		*/
	/* this should be moved to a patch default parameter		*/
	/*--------------------------------------------------------------*/
	patch[0].soil_defaults[0][0].DON_loss_rate = command_line[0].don_value;
	
	/*--------------------------------------------------------------*/
	/* FOR now substitute worldfile m (if > 0) in defaults			*/
	/*--------------------------------------------------------------*/
	if (mpar > ZERO) {
		patch[0].soil_defaults[0][0].original_m = mpar;
		patch[0].soil_defaults[0][0].m = mpar * command_line[0].sen[M];
		patch[0].soil_defaults[0][0].m_z = mpar * command_line[0].sen[M] / 
				patch[0].soil_defaults[0][0].porosity_0;
	}


	/*--------------------------------------------------------------*/
	/* detention store size can vary with both soil and landuse		*/
	/*	use the maximum of the two									*/
	/*--------------------------------------------------------------*/
	patch[0].soil_defaults[0][0].detention_store_size = 
				max(patch[0].landuse_defaults[0][0].detention_store_size,
				patch[0].soil_defaults[0][0].detention_store_size);
	/*--------------------------------------------------------------*/
	/*	Read in the number of  patch base stations 					*/
	/*--------------------------------------------------------------*/
	fscanf(world_file,"%d",&(patch[0].num_base_stations));
	read_record(world_file, record);
	/*--------------------------------------------------------------*/
	/*    Allocate a list of base stations for this patch.			*/
	/*--------------------------------------------------------------*/
	patch[0].base_stations = (struct base_station_object **)
		alloc(patch[0].num_base_stations *
		sizeof(struct base_station_object *),
		"base_stations","construct_patch" );
	/*--------------------------------------------------------------*/
	/*      Read each base_station ID and then point to that base_statio*/
	/*--------------------------------------------------------------*/
	for (i=0 ; i<patch[0].num_base_stations; i++){
		fscanf(world_file,"%d",&(base_stationID));
		read_record(world_file, record);
		/*--------------------------------------------------------------*/
		/*	Point to the appropriate base station in the base       	*/
		/*              station list for this world.					*/
		/*																*/
		/*--------------------------------------------------------------*/
		patch[0].base_stations[i] = assign_base_station(
			base_stationID,
			num_world_base_stations,
			world_base_stations);
	} /*end for*/
	/*--------------------------------------------------------------*/
	/*	Read in number of canopy strata objects in this patch		*/
	/*--------------------------------------------------------------*/
	fscanf(world_file,"%d",&(patch[0].num_canopy_strata));
	read_record(world_file, record);
	
	/*--------------------------------------------------------------*/
	/*	Allocate list of pointers to stratum objects .				*/
	/*--------------------------------------------------------------*/
	patch[0].canopy_strata = ( struct canopy_strata_object ** )
		alloc( patch[0].num_canopy_strata *
		sizeof( struct canopy_strata_object *),
		"canopy_strata","construct_patch");
	
	/*--------------------------------------------------------------*/
	/*      Allocate the patch hourly object.	  */
	/*--------------------------------------------------------------*/
	if ((patch[0].hourly = (struct patch_hourly_object *) calloc(1,
		sizeof(struct patch_hourly_object))) == NULL ){
		fprintf(stderr,"FATAL ERROR: in patch_hourly\n");
		exit(EXIT_FAILURE);
	}
	
	/*--------------------------------------------------------------*/
	/*      Initialize patch level rainand snow stored              */
	/*--------------------------------------------------------------*/
	patch[0].rain_stored = 0.0;
	patch[0].snow_stored = 0.0;
	patch[0].soil_defaults[0][0].daily_fire_litter_turnover = 0.0;
	patch[0].litter.gl_c = 0.0;
	patch[0].litter.gsurf_slope = 0.0;
	patch[0].litter.gsurf_intercept = 0.0;
	
	/*--------------------------------------------------------------*/
	/*	Construct the strata in this patch.						*/
	/*--------------------------------------------------------------*/
	for ( i=0 ; i<patch[0].num_canopy_strata ; i++ ){
		patch[0].canopy_strata[i] = construct_canopy_strata(
			command_line,
			world_file,
			patch,
			num_world_base_stations,
			world_base_stations,defaults);
		/*--------------------------------------------------------------*/
		/*      Aggregate rain and snow stored already for water balance*/
		/*--------------------------------------------------------------*/
		patch[0].rain_stored += patch[0].canopy_strata[i][0].rain_stored
			* patch[0].canopy_strata[i][0].cover_fraction;
		patch[0].snow_stored += patch[0].canopy_strata[i][0].snow_stored
			* patch[0].canopy_strata[i][0].cover_fraction;
		patch[0].soil_defaults[0][0].daily_fire_litter_turnover +=
			patch[0].canopy_strata[i][0].defaults[0][0].epc.daily_fire_turnover
			* patch[0].canopy_strata[i][0].cover_fraction;
		patch[0].litter.gl_c +=
			patch[0].canopy_strata[i][0].defaults[0][0].epc.gl_c
			* patch[0].canopy_strata[i][0].cover_fraction;
		patch[0].litter.gsurf_slope +=
			patch[0].canopy_strata[i][0].defaults[0][0].epc.litter_gsurf_slope
			* patch[0].canopy_strata[i][0].cover_fraction;
		patch[0].litter.gsurf_intercept +=
			patch[0].canopy_strata[i][0].defaults[0][0].epc.litter_gsurf_intercept
			* patch[0].canopy_strata[i][0].cover_fraction;
		patch[0].litter.moist_coef +=
			patch[0].canopy_strata[i][0].defaults[0][0].epc.litter_moist_coef
			* patch[0].canopy_strata[i][0].cover_fraction;
	} /*end for*/
	/*--------------------------------------------------------------*/
	/*	initialize litter capacity				*/
	/*--------------------------------------------------------------*/
	update_litter_interception_capacity(
		patch[0].litter.moist_coef,
		&(patch[0].litter_cs),
		&(patch[0].litter));
	
	/*--------------------------------------------------------------*/
	/*	Define a list of canopy strata layers that can at least	*/
	/*	fit all of the canopy strata.				*/
	/*--------------------------------------------------------------*/
	patch[0].layers = (struct layer_object *) alloc( patch[0].num_canopy_strata *
		sizeof( struct layer_object ),"layers","construct_patch");
	patch[0].num_layers = 0;
	sort_patch_layers(patch);
	
	/*--------------------------------------------------------------*/
	/*	compute actual depth to water tablke			*/
	/*--------------------------------------------------------------*/
	patch[0].unsat_zone_volume = patch[0].sat_deficit + patch[0].unsat_storage;
	patch[0].sat_deficit_z = compute_z_final(
		command_line[0].verbose_flag,
		patch[0].soil_defaults[0][0].porosity_0,
		patch[0].soil_defaults[0][0].porosity_decay,
		patch[0].soil_defaults[0][0].soil_depth,
		0,
		-1*patch[0].sat_deficit);
	patch[0].preday_sat_deficit_z = patch[0].sat_deficit_z;
	return(patch);
} /*end construct_patch.c*/

