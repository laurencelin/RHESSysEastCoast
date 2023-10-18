/*--------------------------------------------------------------*/
/* 																*/
/*				 		patch_daily_F					*/
/*																*/
/*	NAME														*/
/*	patch_daily_F 										*/
/*				 - performs cycling and output of a patch		*/
/*																*/
/*																*/
/*	SYNOPSIS													*/
/*	void patch_daily_F(								*/
/*						struct	world_object	*,				*/
/*						struct	basin_object	*,				*/
/*						struct	hillslope_object	*,			*/
/*						struct	zone_object		*,				*/
/*						struct patch_object	,					*/
/*						struct command_line_object ,			*/
/*						struct tec_entry,						*/
/*						struct date)							*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	This routine performs simulation cycles on an identified	*/
/*	canopy_stata in the patch. The routine also prints out results*/
/*	where specified by current tec events files.				*/
/*																*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	March 4, 1997	  C.Tague 									*/
/*	baseflow moved to hillslope level 							*/
/*																*/
/*	March 9, 1997	  C.Tague 									*/
/*	moss routines only called if moss depth > 0 				*/
/*																*/
/*																*/
/*	March 13, 19987 - RAF			*/
/*	Took comments offurface_daily_F call	*/
/*	since it was under an if condition 	*/
/*	Put the comments around the if condition*/  
/*																*/
/*	July 28, 1997 - C.Tague			*/
/*	removed capillary rise routines 	*/
/*											*/
/*	Sep 2 1997								*/
/*	Removed references to moss or humic layer	*/
/*	Everything is now either a stratum or soil	*/
/*							*/
/*	Sep 3 1997					*/
/*	Took out condition that strata under snowpack	*/
/*	must be at z>0					*/
/*							*/
/*	Sept 29 1997 CT   				*/
/*	unsat drainage and cap rise are now based 	*/
/*	on actual depth to water table which is		*/
/*	determined by assuming a porosity decay		*/
/*	with depth;  transpiration is now sep.		*/
/*	into sat_zone and unsat_zone trans		*/
/*							*/
/*	Oct 22 1997 CT   				*/
/*	included calculation of water balance		*/
/*							*/
/*	Oct 22 1997 CT   				*/
/*	canopy_strata sublimation now added to 		*/
/*	patch level sum of canopy evaporation		*/
/*							*/
/*	added check_zero_stores				*/
/*							*/
/*	Feb 2010 AD											*/
/*	Added detention store evaporation, including more	*/
/*	substantial updates to surface daily F				*/
/*														*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rhessys.h"

void		patch_daily_F(
						  struct	world_object	*world,
						  struct	basin_object	*basin,
						  struct	hillslope_object	*hillslope,
						  struct	zone_object		*zone,
						  struct 	patch_object 	*patch,
						  struct 	command_line_object *command_line,
						  struct	tec_entry		*event,
						  struct	date 			current_date)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	double compute_delta_water(int, double, double,	double, double, double);
	
	
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
	
	void canopy_stratum_daily_F(
		struct world_object *,
		struct basin_object *,
		struct hillslope_object *,
		struct zone_object *,
		struct patch_object *,
		struct layer_object *,
		struct canopy_strata_object *,
		struct canopy_strata_object *,
		struct command_line_object *,
		struct tec_entry *,
		struct date);
	
	void   surface_daily_F(
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

	
	double	snowpack_daily_F (
		struct date,
		int,
		struct zone_object *,
		struct patch_object *,
		struct	snowpack_object	*,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double *,
		double *,
		double *,
		double *,
		double *,
		double *,
		double,
		double,
		double,
		double,
		double,
		double,
		int);
	
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

	double  compute_surface_heat_flux(
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
	
	double	compute_unsat_zone_drainage(
		int,
		int,
		double,
		double,
		double,
		double,
		double,
		double, double);
	
	double  compute_radiative_fluxes(
		int,
		double *,
		double *,
		double,
		double,
		double);
	
	
		/*-------------------------------------
		double  compute_stability_correction(
		int ,
		double,
		double,
		double,
		double,
		double,
		double);
	----------------------------------------*/
	
	double  compute_z_final(
		int,
		double,
		double,
		double,
		double,
		double);

	int 	check_zero_stores(
		struct  soil_c_object   *,
		struct  soil_n_object   *,
		struct  litter_c_object *,
		struct  litter_n_object *
		);
	
	int	update_decomp(
		struct	date,
		struct  soil_c_object   *,
		struct  soil_n_object   *,
		struct  litter_c_object *,
		struct  litter_n_object *,
		struct cdayflux_patch_struct *,
		struct ndayflux_patch_struct *,
		struct patch_object *,
        struct command_line_object *);
	
	int	update_dissolved_organic_losses(
		struct	date,
		double,
		struct  soil_c_object   *,
		struct  soil_n_object   *,
		struct  litter_c_object *,
		struct  litter_n_object *,
		struct cdayflux_patch_struct *,
		struct ndayflux_patch_struct *,
        struct patch_object *, int);
	
	int	update_septic(
		struct	date,
		struct  patch_object   *);

	int	update_nitrif(
		struct  soil_c_object   *,
		struct  soil_n_object   *,
		struct cdayflux_patch_struct *,
		struct ndayflux_patch_struct *,
		struct patch_object *,
        struct command_line_object *,
		double);
	
	int	update_denitrif(
		struct  soil_c_object   *,
		struct  soil_n_object   *,
		struct cdayflux_patch_struct *,
		struct ndayflux_patch_struct *,
        struct patch_object *,
        struct command_line_object *,
		double);

	
	int	resolve_sminn_competition(
		struct  soil_n_object   *,
		struct patch_object *,
        struct command_line_object *,
		struct ndayflux_patch_struct *);
	
	void   canopy_stratum_growth(
		struct world_object *,
		struct basin_object *,
		struct hillslope_object *,
		struct zone_object *,
		struct patch_object *,
		struct canopy_strata_object *,
		struct command_line_object *,
		struct tec_entry *,
		struct date);

	int update_gw_drainage(
			struct patch_object *,
			struct hillslope_object *,
			struct zone_object *,
			struct command_line_object *,
			struct date);
			
	double	penman_monteith(
		int,
		double,
		double,
		double,
		double,
		double,
		double,
		int);
		
	double	compute_diffuse_radiative_fluxes(
		int,
		double *,
		double *,
		double,
		double,
		double,
		double,
		double,
		double);

	
	double	compute_direct_radiative_fluxes(
		int,
		double *,
		double *,
		double,
		double,
		double,
		double,
		double,
		double);
	
	double compute_toc_wind(int,
							double,
							double,
							double,
							double);
	
	double  compute_ra_overstory(
								 int     ,
								 double  ,
								 double  ,
								 double *,
								 double *,
								 double *,
								 double *,
								 double  ,
								 double  ,
								 double	,
								 double *,
								 double *);
	
    double    compute_N_leached(
                                int verbose_flag,
                                double total_nitrate,
                                double Qout,
                                double N_decay_rate,
                                double activedepthz,
                                double N_absorption_rate,
                                int signal,
                                struct patch_object *patch);

	long julday( struct date);
	int	compute_Lstar(int	verbose_flag,
					  struct	basin_object	*basin,
					  struct	zone_object	*zone,
					  struct	patch_object	*patch);
	/*--------------------------------------------------------------*/
	/*  Local variable definition.                                  */
	/*--------------------------------------------------------------*/
	int	layer;
	int stratum, ch, inx;
	int	vegtype;
	int dum;
    int j,i,d;
	double	cap_rise, cap_rise_unsat, tmp, wilting_point;
	double	delta_unsat_zone_storage;
	double  infiltration, lhvap;
	double	infiltration_ini;
	double	infiltration_fin;
	double	net_inflow, theta;
	double	preday_snowpack_height;
	double	sat_zone_patch_demand;
	double	sat_zone_patch_demand_initial;
	double	available_sat_water;
	double	temp,scale;
	double	unsat_zone_patch_demand;
	double	unsat_zone_patch_demand_initial;
	double  add_field_capacity;
	double	water_above_field_cap;
	double	water_below_field_cap;
	double 	duration, irrigation;
	double	snow_melt_input;
	double  fertilizer_NO3, fertilizer_NH4;
	double	resp, transpiration_reduction_percent;
	double 	surfaceN_to_soil;
	double	FERT_TO_SOIL;
	double	pond_height;
	double tmpra, tmpga, tmpgasnow, tmpwind, tmpwindcan, tmpwindsnow, tmpustar;
	double Kup_direct_snow, Kup_diffuse_snow;
	double Kdown_direct_covered, Kdown_diffuse_covered, Kdown_direct_exposed, Kdown_diffuse_exposed;
	double Kup_direct_snow_covered, Kup_diffuse_snow_covered, Kup_direct_snow_exposed, Kup_diffuse_snow_exposed;
	double PAR_direct_covered, PAR_diffuse_covered, PAR_direct_exposed, PAR_diffuse_exposed;
	double snow_melt_covered, snow_melt_exposed;
	double 	rz_drainage,unsat_drainage;
	struct	canopy_strata_object	*strata;
	struct	litter_object	*litter;
	struct  dated_sequence	clim_event;
    double totalfc;
	/*--------------------------------------------------------------*/
	/*	We assume the zone soil temp applies to the patch as well.	*/
	/* 	unless we are using the surface energy iteration code 	in which */
	/* 	case we use the temperature from the previous day		*/
	/*	alos for the Kdowns and PAR (for now Ldown can be kept )	*/
	/*--------------------------------------------------------------*/

    if (command_line[0].surface_energy_flag == 0) patch[0].Tsoil = zone[0].metv.tsoil;// +1-patch[0].Ksat_vertical;

    d=0;
	patch[0].Kdown_direct = zone[0].Kdown_direct;
	patch[0].Kdown_diffuse = zone[0].Kdown_diffuse;
	patch[0].PAR_direct = zone[0].PAR_direct;
	patch[0].PAR_diffuse = zone[0].PAR_diffuse;
	patch[0].evaporation_surf = 0.0;
	patch[0].potential_evaporation = 0.0;
	patch[0].Ldown = zone[0].Ldown;
	patch[0].Ldown_night = zone[0].Ldown_night;
	patch[0].Ldown_day = zone[0].Ldown_day;
	patch[0].Ldown_final = 0.0;
	patch[0].Ldown_final_night = 0.0;
	patch[0].Ldown_final_day = 0.0;
	
	patch[0].Kstar_canopy = 0.0;
	patch[0].Kstar_canopy_final = 0.0;
	patch[0].LE_canopy = 0.0;
	patch[0].LE_canopy_final = 0.0;
	patch[0].Kdown_direct_subcanopy = 0.0;
	patch[0].Kdown_diffuse_subcanopy = 0.0;
	patch[0].Ldown_subcanopy = 0.0;
	
	patch[0].wind_final = 0.0;
	patch[0].windsnow_final = 0.0;
	patch[0].ustar = 0.0;
	patch[0].ustar_final = 0.0;
	
	patch[0].snowpack.overstory_fraction = 0.0;
	patch[0].overstory_fraction = 0.0;
	
	patch[0].ga = 0.0;
	patch[0].gasnow = 0.0;
	
	tmpwind = zone[0].wind;
	
	dum = 0;
	patch[0].Kup_direct = 0.0;
	patch[0].Kup_diffuse = 0.0;
	patch[0].Kup_direct_final = 0.0;
	patch[0].Kup_diffuse_final = 0.0;
	Kup_direct_snow = 0.0;
	Kup_diffuse_snow = 0.0;
	
	patch[0].Kdown_direct_bare = 0.0;
	
	patch[0].snowpack.sublimation = 0.0;
	patch[0].snowpack.Rnet = 0.0;
	patch[0].snowpack.Q_H = 0.0;
	patch[0].snowpack.Q_LE = 0.0;
	patch[0].snowpack.Q_rain = 0.0;
	patch[0].snowpack.Q_melt = 0.0;
	
	patch[0].LE_soil = 0.0;
	
	
	patch[0].snowpack.K_reflectance = 0.0;
	patch[0].snowpack.K_absorptance = 0.0;
	patch[0].snowpack.PAR_reflectance = 0.0;
	patch[0].snowpack.PAR_absorptance = 0.0;
	patch[0].snowpack.Kstar_direct = 0.0;
	patch[0].snowpack.Kstar_diffuse = 0.0;
	patch[0].snowpack.APAR_direct = 0.0;
	patch[0].snowpack.APAR_diffuse = 0.0;
	
	Kdown_direct_covered = 0.0;
	Kdown_diffuse_covered = 0.0;
	Kdown_direct_exposed = 0.0;
	Kdown_diffuse_exposed = 0.0;
	Kup_direct_snow_covered = 0.0;
	Kup_diffuse_snow_covered = 0.0;
	Kup_direct_snow_exposed = 0.0;
	Kup_diffuse_snow_exposed = 0.0;	
	PAR_direct_covered = 0.0;
	PAR_diffuse_covered = 0.0;
	PAR_direct_exposed = 0.0;
	PAR_diffuse_exposed = 0.0;
	snow_melt_covered = 0.0;
	snow_melt_exposed = 0.0;
	
	patch[0].exfiltration_unsat_zone = 0.0;
	patch[0].exfiltration_sat_zone = 0.0;
	
	patch[0].T_canopy = zone[0].metv.tavg;
	patch[0].T_canopy_final = 0.0;

	patch[0].ex_inundation_depth = 0.0;
	patch[0].ex_inundation_dur = 0.0;

	if ( command_line[0].verbose_flag == -5 ){
        printf("\nPATCH DAILY F:");
	}

	/*--------------------------------------------------------------*/
	/*	INUNDATION	*/
	/*--------------------------------------------------------------*/
	#include <string.h>

	// Declare variables for the column names
	double inundation_PatchID[50];
	double inundation_date[50];
	double inundation_duration[50];
	double inundation_depth[50];
	FILE *file;

    	// Open the CSV file for reading
  	// char url[255];
   	//USER INPUT URL
   	// printf("https://github.com/hanneborstlap/RHESSysEastCoast_orig/blob/inundation/CobbMill_output_edited.csv");

	// scanf("%s", );
   	
    	// file = fopen("https://github.com/hanneborstlap/RHESSysEastCoast_orig/blob/inundation/CobbMill_output_edited.csv", "r");

	struct date createDateFromDateString(const char* dateString) {
    	struct date result;
    	char* token;
    	char* copy = strdup(dateString); // Make a copy of the input string
    	// Tokenize the string using "-" as the delimiter
    	token = strtok(copy, "-");
    	result.month = atol(token);
    	token = strtok(NULL, "-");
    	result.day = atoi(token);
    	token = strtok(NULL, "-");
    	result.year = atoi(token);
    	free(copy); 
    	return result;
	}

	file = fopen("/scratch/tpv4jw/RHESSys/5_INUNDATION/CobbMill_output_edited.csv", "r");
	int count = 0;
	while (fscanf(file, "%lf, %s, %lf, %lf", &inundation_PatchID[count], &inundation_date[count], &inundation_duration[count], &inundation_depth[count]) == 4) {
		// printf("PatchID: %s\n", inundation_patchID);
        	// printf("Date: %s\n", inundation_date);
        	// printf("Depth: %s\n", inundation_depth);
        	// printf("Duration: %s\n", inundation_duration);
		count++;
    	}

	struct date inundation_date_f = createDateFromDateString(inundation_date);

	for (int i = 0; i < count; i++) {
		if (patch[0].ID == inundation_PatchID[i]) {
		if (julday(inundation_date_f[i]) != julday(current_date)) {
			patch[0].ex_inundation_depth = 0.0; 
			patch[0].ex_inundation_dur = 0.0; 
		}
		if (julday(inundation_date_f[i]) == julday(current_date)) {
			patch[0].ex_inundation_depth = inundation_depth[i]; 
			patch[0].ex_inundation_depth = duration[i]; 
		}
	} 

    }

    // Close the file when you're done
    fclose(file);
    return 0;

	/*--------------------------------------------------------------*/
	/*	Set the patch rain and snow throughfall equivalent to the	*/
	/*	rain and snow coming down over the zone.					*/
	/* check to see if there are base station inputs 		*/
	/*--------------------------------------------------------------*/

	if (patch[0].base_stations != NULL) {
		inx = patch[0].base_stations[0][0].dated_input[0].irrigation.inx;
		if (inx > -999) {
			clim_event = patch[0].base_stations[0][0].dated_input[0].irrigation.seq[inx];
			while (julday(clim_event.edate) < julday(current_date)) {
				patch[0].base_stations[0][0].dated_input[0].irrigation.inx += 1;
				inx = patch[0].base_stations[0][0].dated_input[0].irrigation.inx;
				clim_event = patch[0].base_stations[0][0].dated_input[0].irrigation.seq[inx];
            }//while
			if ((clim_event.edate.year != 0) && ( julday(clim_event.edate) == julday(current_date)) ) {
				irrigation = clim_event.value;
            }else irrigation = 0.0;
        }else irrigation = 0.0;
    }else irrigation = 0.0;
    //patch[0].landuse_defaults[0][0].irrigation; format change;
    //patch[0].landuse_defaults[0][0].irrigation is used below
    
    
    patch[0].grassIrrigation_m = 0.0;
    if(command_line[0].grassIrrigation_flag==1 && patch[0].drainage_type>0 && patch[0].drainage_type % actionIRRIGRATION==0){
        //---------------------------- lawn or crop
        // first calculate water demand by lawn or crop
        double grassET, grassPET;
        for ( j=0 ; j<patch[0].num_canopy_strata ; j++ ){
            if(patch[0].canopy_strata[j][0].defaults[0][0].epc.veg_type == GRASS && patch[0].canopy_strata[j][0].phen.gwseasonday>0){
                //gwseasonday > 0 when plant is growing; its -1 otherwise.
                grassET = patch[0].canopy_strata[j][0].evaporation + patch[0].canopy_strata[j][0].sublimation;
                grassPET = grassET;
                grassET += patch[0].canopy_strata[j][0].transpiration_sat_zone + patch[0].canopy_strata[j][0].transpiration_unsat_zone;
                grassPET += patch[0].canopy_strata[j][0].PET;
                
                if( grassPET>0 ){
                      // 4mm/day/m2 irrigation -> 0.004 m/day/m2 (max rate)
                      // patch[0].landuse_defaults[0][0].irrigation be the daily max. irrigation rate
                    patch[0].grassIrrigation_m += (1-grassET/grassPET)*patch[0].landuse_defaults[0][0].irrigation * patch[0].canopy_strata[j][0].cover_fraction;
                }//if
            }// if
        }// for loop j
        
        // second calculate water supplies from known sources
        double sourceTransferWater = 0.0;
        double totalTransferWater = 0.0;
        double tmp_fraction = 0.0;
        for (i =0; i < patch[0].innundation_list[d].num_drainIN_irrigation; i++){
        
            // from surface to surface
            sourceTransferWater = min(
                              patch[0].innundation_list[d].drainIN_irrigation[i].maxDailyDrain,
                              patch[0].innundation_list[d].drainIN_irrigation[i].DrainFrac * patch[0].grassIrrigation_m); // a fraction of water demand from this konwn source
            
            patch[0].innundation_list[d].drainIN_irrigation[i].transfer_flux_surf = max(0.0, min(sourceTransferWater * patch[0].innundation_list[d].drainIN_irrigation[i].propDrainFrmSurf,
                patch[0].innundation_list[d].drainIN_irrigation[i].patch[0].detention_store));// actual available water from known source
            
            // perform water and solute transfer
            patch[0].detention_store += patch[0].innundation_list[d].drainIN_irrigation[i].transfer_flux_surf;
            tmp_fraction = min(1.0,(patch[0].innundation_list[d].drainIN_irrigation[i].patch[0].detention_store>0? patch[0].innundation_list[d].drainIN_irrigation[i].transfer_flux_surf/patch[0].innundation_list[d].drainIN_irrigation[i].patch[0].detention_store : 0.0));
            patch[0].surface_NO3 += tmp_fraction * patch[0].innundation_list[d].drainIN_irrigation[i].patch[0].surface_NO3;
            patch[0].surface_NH4 += tmp_fraction * patch[0].innundation_list[d].drainIN_irrigation[i].patch[0].surface_NH4;
            patch[0].surface_DON += tmp_fraction * patch[0].innundation_list[d].drainIN_irrigation[i].patch[0].surface_DON;
            patch[0].surface_DOC += tmp_fraction * patch[0].innundation_list[d].drainIN_irrigation[i].patch[0].surface_DOC;
            // substrate the source water and solute
            tmp_fraction = min(1.0,max(0.0, 1.0 - tmp_fraction));
            patch[0].innundation_list[d].drainIN_irrigation[i].patch[0].detention_store *= tmp_fraction;
            patch[0].innundation_list[d].drainIN_irrigation[i].patch[0].surface_NO3 *= tmp_fraction;
            patch[0].innundation_list[d].drainIN_irrigation[i].patch[0].surface_NH4 *= tmp_fraction;
            patch[0].innundation_list[d].drainIN_irrigation[i].patch[0].surface_DON *= tmp_fraction;
            patch[0].innundation_list[d].drainIN_irrigation[i].patch[0].surface_DOC *= tmp_fraction;
            
            // drawing from deep GW
            patch[0].innundation_list[d].drainIN_irrigation[i].transfer_flux_sub = max(0.0,min(sourceTransferWater * (1.0 - patch[0].innundation_list[d].drainIN_irrigation[i].propDrainFrmSurf), hillslope[0].gw.storage*hillslope[0].area/patch[0].area)); // actual available water from known source
            
            // perform water and solute transfer
            patch[0].detention_store += patch[0].innundation_list[d].drainIN_irrigation[i].transfer_flux_sub;
            tmp_fraction = min(1.0,(hillslope[0].gw.storage>0? patch[0].innundation_list[d].drainIN_irrigation[i].transfer_flux_sub*patch[0].area/hillslope[0].area/hillslope[0].gw.storage : 0.0));
            patch[0].surface_NO3 += tmp_fraction * hillslope[0].gw.NO3;
            patch[0].surface_NH4 += tmp_fraction * hillslope[0].gw.NH4;
            patch[0].surface_DON += tmp_fraction * hillslope[0].gw.DOC;
            patch[0].surface_DOC += tmp_fraction * hillslope[0].gw.DON;
            // substrate the source water and solute
            tmp_fraction = min(1.0,max(0.0, 1.0 - tmp_fraction));
            hillslope[0].gw.storage *= tmp_fraction;
            hillslope[0].gw.NO3 *= tmp_fraction;
            hillslope[0].gw.NH4 *= tmp_fraction;
            hillslope[0].gw.DOC *= tmp_fraction;
            hillslope[0].gw.DON *= tmp_fraction;
            
            // calculate total water gain from the source to this patch
            totalTransferWater += patch[0].innundation_list[d].drainIN_irrigation[i].transfer_flux_sub + patch[0].innundation_list[d].drainIN_irrigation[i].transfer_flux_surf;
        }// end of for loop i; number of sources
        if(patch[0].innundation_list[d].num_drainIN_irrigation>0){ patch[0].grassIrrigation_m = totalTransferWater; }
        // note:: solutes and water are transferred by the processes above
        // note:: solutes at sources are substrated by the processes above
        
    }// irrigation
    
    patch[0].sewerdrained = 0.0;
    patch[0].sewerdrained_NO3 = 0.0;
    patch[0].sewerdrained_NH4 = 0.0;
    patch[0].sewerdrained_DOC = 0.0;
    patch[0].sewerdrained_DON = 0.0;
    if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate ||
    patch[0].soil_ns.nitrate<0 ||
    patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn ||
    patch[0].soil_ns.sminn<0 ||
    patch[0].soil_ns.DON!=patch[0].soil_ns.DON ||
    patch[0].soil_ns.DON<0 ||
    patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC ||
    patch[0].soil_cs.DOC<0 ||
    patch[0].sat_NO3!=patch[0].sat_NO3 ||
    patch[0].sat_NO3<0 || patch[0].sat_NH4!=patch[0].sat_NH4 ||
    patch[0].sat_NH4<0 || patch[0].sat_DOC!=patch[0].sat_DOC ||
    patch[0].sat_DOC<0) printf("patch daily F0N %d-%d-%d [%d,%d,%d,%d]{%e,%e,%e,%e}[%e,%e,%e,%e]\n",
       current_date.year, current_date.month, current_date.day,
       patch[0].ID, patch[0].drainage_type, actionPIPEDRAIN, actionSEWER,
       patch[0].soil_ns.nitrate,
       patch[0].soil_ns.sminn,
       patch[0].soil_cs.DOC,
       patch[0].soil_ns.DON,
       patch[0].sat_NO3, patch[0].sat_NH4, patch[0].sat_DOC, patch[0].sat_DON);
    
    if(command_line[0].sewer_flag==1 && patch[0].drainage_type>0 && patch[0].drainage_type % actionSEWER==0){
        //---------------------------- lawn / impervious drainage (sewerage)
        
        if(patch[0].sat_deficit_z < patch[0].landuse_defaults[0][0].sewer_infiltrationSatDefZThreshold){
            // infiltration into the sewer
            patch[0].sewerdrained += min(min((patch[0].landuse_defaults[0][0].sewer_infiltrationSatDefZThreshold-patch[0].sat_deficit_z)*patch[0].landuse_defaults[0][0].sewer_infiltrationSatDefZHeadSpace,1.0)*patch[0].landuse_defaults[0][0].sewer_infiltrationRate, patch[0].landuse_defaults[0][0].sewer_infiltrationSatDefZThreshold-patch[0].sat_deficit_z);
            // Moore et al. 2015; when water is above the empty header, no solute leaching out
            if(patch[0].sewerdrained>0 && command_line[0].grow_flag > 0){
                // how? compute_N_leached? here?
                double tmpNout = 0.0;
                tmpNout = compute_N_leached(
                                         command_line[0].verbose_flag,
                                         patch[0].sat_NO3, //patch[0].soil_ns.nitrate,
                                         patch[0].sewerdrained,
                                         patch[0].soil_defaults[0][0].NO3decayRate,
                                         patch[0].soil_defaults[0][0].active_zone_z,
                                         patch[0].soil_defaults[0][0].NO3_adsorption_rate,
                                         29, patch);
                //printf("patch [%d] sewerdrained(%e) NO3{%e - > %e}\n", patch[0].ID, patch[0].sewerdrained, patch[0].soil_ns.nitrate, tmpNout);
                patch[0].sewerdrained_NO3 += tmpNout;
                patch[0].sat_NO3 -= tmpNout; //patch[0].soil_ns.nitrate -= tmpNout;
    //            if(patch[0].sewerdrained > 0) printf("patch[%d, %d,%d,%d] sewer yields (%e) and %e[%e, %e,%e]\n",
    //                                                   patch[0].ID, current_date.year,current_date.month,current_date.day,
    //                                                   patch[0].sewerdrained, patch[0].sewerdrained_NO3, tmpNout,
    //                                                   patch[0].sat_deficit_z, patch[0].soil_ns.nitrate);
                
                tmpNout = compute_N_leached(
                                         command_line[0].verbose_flag,
                                         patch[0].sat_NH4, //patch[0].soil_ns.sminn,
                                         patch[0].sewerdrained,
                                         patch[0].soil_defaults[0][0].NH4decayRate,
                                         patch[0].soil_defaults[0][0].active_zone_z,
                                         patch[0].soil_defaults[0][0].NH4_adsorption_rate,
                                         32,patch);
                //printf("patch [%d] sewerdrained(%e) NH4{%e - > %e}\n", patch[0].ID,patch[0].sewerdrained,patch[0].soil_ns.sminn,tmpNout);
                patch[0].sewerdrained_NH4 += tmpNout;
                patch[0].sat_NH4 -= tmpNout; //patch[0].soil_ns.sminn -= tmpNout;

                tmpNout = compute_N_leached(
                                         command_line[0].verbose_flag,
                                         patch[0].sat_DON, //patch[0].soil_ns.DON,
                                         patch[0].sewerdrained,
                                         patch[0].soil_defaults[0][0].DOMdecayRate,
                                         patch[0].soil_defaults[0][0].active_zone_z,
                                         patch[0].soil_defaults[0][0].DON_adsorption_rate,
                                         35,patch);
                //printf("patch [%d] sewerdrained(%e) DON{%e - > %e}\n", patch[0].ID,patch[0].sewerdrained,patch[0].soil_ns.DON,tmpNout);
                patch[0].sewerdrained_DON += tmpNout;
                patch[0].sat_DON -= tmpNout; //patch[0].soil_ns.DON -= tmpNout;

                tmpNout = compute_N_leached(
                                         command_line[0].verbose_flag,
                                         patch[0].sat_DOC,//patch[0].soil_cs.DOC,
                                         patch[0].sewerdrained,
                                         patch[0].soil_defaults[0][0].DOMdecayRate,
                                         patch[0].soil_defaults[0][0].active_zone_z,
                                         patch[0].soil_defaults[0][0].DOC_adsorption_rate,
                                         38,patch);
                //printf("patch [%d] sewerdrained(%e) DOC{%e - > %e}\n", patch[0].ID,patch[0].sewerdrained,patch[0].soil_cs.DOC,tmpNout);
                patch[0].sewerdrained_DOC += tmpNout;
                patch[0].sat_DOC -= tmpNout; //patch[0].soil_cs.DOC -= tmpNout;
                
            }// growth_flag
        }else{
            // exfiltration out from sewer
            patch[0].sewerdrained -= min(min((patch[0].sat_deficit_z-patch[0].landuse_defaults[0][0].sewer_infiltrationSatDefZThreshold)*patch[0].landuse_defaults[0][0].sewer_exfiltrationDepthVol,1.0)*patch[0].landuse_defaults[0][0].sewer_exfiltrationRate,
                patch[0].sat_deficit_z-patch[0].landuse_defaults[0][0].sewer_infiltrationSatDefZThreshold);
            // use negative value to indicate exfiltration from the pipe.
            // count for the filling up process
            
            patch[0].sewerdrained_NH4 -= fabs(patch[0].sewerdrained) * patch[0].landuse_defaults[0][0].sewerNH4c;
            patch[0].sewerdrained_NO3 -= fabs(patch[0].sewerdrained) * patch[0].landuse_defaults[0][0].sewerNO3c; // (kgN/m3) * ["m3/m2" * area m2] / area m2 = kgN/area m2
            patch[0].sewerdrained_DON -= fabs(patch[0].sewerdrained) * patch[0].landuse_defaults[0][0].sewerDONc;
            patch[0].sewerdrained_DOC -= fabs(patch[0].sewerdrained) * patch[0].landuse_defaults[0][0].sewerDOCc;
        }//end of if
            
        
        if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
        patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
        patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
        patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
        patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
        patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
        patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
        patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
        patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0) printf("patch daily F1N %d-%d-%d [%d,%d,%d,%d] [%e %e %e] [%e %e %e] [%e %e %e]\n",
           current_date.year, current_date.month, current_date.day,
           patch[0].ID, patch[0].drainage_type, actionPIPEDRAIN, actionSEWER,
           patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
           patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
           patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);
        
        patch[0].sat_deficit += patch[0].sewerdrained;
        
        if(patch[0].sat_deficit >= 0){
            patch[0].sat_deficit = min(patch[0].sat_deficit,patch[0].soil_defaults[0][0].soil_water_cap);
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
        
        totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index];
        //totalfc *= (1.0-patch[0].basementFrac); // <---- second thought on this, Oct 8, 2019; basement is 3m at most
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
            printf("patch_daily_F(1): (%d,%d,%d) %lf %lf %lf %lf\n",
                   patch[0].ID, patch[0].soil_defaults[0][0].ID, patch[0].sat_def_pct_index,
                   patch[0].sat_deficit, patch[0].sat_deficit_z,
                   patch[0].rootzone.field_capacity, patch[0].field_capacity);
        }//debug
    }// Laurence: this is for SLB
    
    patch[0].pipedrainYield = 0.0;
    patch[0].pipedrainYield_NO3 = 0.0;
    patch[0].pipedrainYield_NH4 = 0.0;
    patch[0].pipedrainYield_DON = 0.0;
    patch[0].pipedrainYield_DOC = 0.0;
    if(patch[0].drainage_type>0 && patch[0].drainage_type % actionPIPEDRAIN==0){
        
        if(patch[0].sat_deficit_z < 1.0){
            patch[0].pipedrainYield += min(patch[0].available_soil_water, 0.0005*max(0.0, 1.0-patch[0].sat_deficit_z*0.143));
        }
        if(patch[0].pipedrainYield>0 && command_line[0].grow_flag > 0){
            // how? compute_N_leached?
            double tmpNout = 0.0;
            tmpNout = compute_N_leached(
                                     command_line[0].verbose_flag,
                                     patch[0].sat_NO3, //patch[0].soil_ns.nitrate,
                                     patch[0].pipedrainYield,
                                     patch[0].soil_defaults[0][0].NO3decayRate,
                                     patch[0].soil_defaults[0][0].active_zone_z,
                                     patch[0].soil_defaults[0][0].NO3_adsorption_rate,
                                     29,patch);
            //printf("patch [%d] pipedrainYield(%e) NO3{%e - > %e}\n", patch[0].ID,patch[0].pipedrainYield,patch[0].soil_ns.nitrate,tmpNout);
            patch[0].pipedrainYield_NO3 += tmpNout;
            patch[0].sat_NO3 -= tmpNout; //patch[0].soil_ns.nitrate -= tmpNout;
//            if(patch[0].pipedrainYield > 0) printf("patch[%d, %d,%d,%d] pipe yields (%e) and %e[%e, %e,%e]\n",
//                                                   patch[0].ID, current_date.year,current_date.month,current_date.day,
//                                                   patch[0].pipedrainYield, patch[0].pipedrainYield_NO3, tmpNout,
//                                                   patch[0].sat_deficit_z, patch[0].soil_ns.nitrate);

            tmpNout = compute_N_leached(
                                     command_line[0].verbose_flag,
                                     patch[0].sat_NH4, //patch[0].soil_ns.sminn,
                                     patch[0].pipedrainYield,
                                     patch[0].soil_defaults[0][0].NH4decayRate,
                                     patch[0].soil_defaults[0][0].active_zone_z,
                                     patch[0].soil_defaults[0][0].NH4_adsorption_rate,
                                     32,patch);
            //printf("patch [%d] pipedrainYield(%e) NH4{%e - > %e}\n", patch[0].ID,patch[0].pipedrainYield,patch[0].soil_ns.sminn,tmpNout);
            patch[0].pipedrainYield_NH4 += tmpNout;
            patch[0].sat_NH4 -= tmpNout; //patch[0].soil_ns.sminn -= tmpNout;

            tmpNout = compute_N_leached(
                                     command_line[0].verbose_flag,
                                     patch[0].sat_DON, //patch[0].soil_ns.DON,
                                     patch[0].pipedrainYield,
                                     patch[0].soil_defaults[0][0].DOMdecayRate,
                                     patch[0].soil_defaults[0][0].active_zone_z,
                                     patch[0].soil_defaults[0][0].DON_adsorption_rate,
                                     35,patch);
            //printf("patch [%d] pipedrainYield(%e) DON{%e - > %e}\n", patch[0].ID,patch[0].pipedrainYield,patch[0].soil_ns.DON,tmpNout);
            patch[0].pipedrainYield_DON += tmpNout;
            patch[0].sat_DON -= tmpNout; //patch[0].soil_ns.DON -= tmpNout;

            tmpNout = compute_N_leached(
                                     command_line[0].verbose_flag,
                                     patch[0].sat_DOC, //patch[0].soil_cs.DOC,
                                     patch[0].pipedrainYield,
                                     patch[0].soil_defaults[0][0].DOMdecayRate,
                                     patch[0].soil_defaults[0][0].active_zone_z,
                                     patch[0].soil_defaults[0][0].DOC_adsorption_rate,
                                     38,patch);
            //printf("patch [%d] pipedrainYield(%e) DOC{%e - > %e}\n", patch[0].ID,patch[0].pipedrainYield,patch[0].soil_cs.DOC,tmpNout);
            patch[0].pipedrainYield_DOC += tmpNout;
            patch[0].sat_DOC -= tmpNout; //patch[0].soil_cs.DOC -= tmpNout;

        }// growth flag

        patch[0].sat_deficit += patch[0].pipedrainYield;
        
        if(patch[0].sat_deficit >= 0){
             patch[0].sat_deficit = min(patch[0].sat_deficit,patch[0].soil_defaults[0][0].soil_water_cap);
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
         
        totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index];
         //totalfc *= (1.0-patch[0].basementFrac); // <---- second thought on this, Oct 8, 2019; basement is 3m at most
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
             printf("patch_daily_I(2): (%d,%d,%d) %lf %lf %lf %lf\n", patch[0].ID,
                    patch[0].sat_deficit, patch[0].sat_deficit_z,
                    patch[0].rootzone.field_capacity, patch[0].field_capacity);
         }//debug
         
    }// Laurence: this is for SLB
    if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
    patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
    patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
    patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
    patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
    patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
    patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
    patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
    patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0) printf("patch daily F2N %d-%d-%d [%d,%d,%d,%d] [%e %e %e] [%e %e %e] [%e %e %e]\n",
       current_date.year, current_date.month, current_date.day,
       patch[0].ID, patch[0].drainage_type, actionPIPEDRAIN, actionSEWER,
       patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
       patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
       patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);
    
	/*--------------------------------------------------------------*/
	/*	process any daily rainfall				*/
	/*--------------------------------------------------------------*/
    patch[0].rain_throughfall = zone[0].rain;
    // problem: irrigation should be adding to rain_throughfall
    // irrigation should be added to patch[0].detention_store @ LINE 1583

	/* the N_depo is add in patch_hourly.c in hourly */
	/* it could be washed away hourly or daily, depending on whether the precipitation data is hourly or daily */
    // patch_hourly.c runs even when no hourly; zone[0].ndep is divided by hourly and by the directNdep fraction in patch_hourly.c
    // then passed into canopy_stratum_hourly.c and stored on the stratum[0].NO3_stored
    // finally is released from stratum[0].NO3_stored in patch_dailyF.c and canopy_stratum_dailyF.c
	patch[0].NO3_throughfall = 0; // this is updated by the hourly accumulated stratum[0].NO3_stored in the canopy_stratum_dailyF.c


	if (command_line[0].snow_scale_flag == 1)
		patch[0].snow_throughfall = zone[0].snow * patch[0].snow_redist_scale;
    else
        patch[0].snow_throughfall = zone[0].snow;

	patch[0].wind = zone[0].wind;
	patch[0].windsnow = zone[0].wind;

	patch[0].precip_with_assim += patch[0].rain_throughfall + patch[0].snow_throughfall;

    
    patch[0].septicReleaseQ_m = 0.0;
	if (command_line[0].septicProcess_flag==1 && patch[0].drainage_type>0 && patch[0].drainage_type % actionSeptic==0) {
        double sourceTransferWater = 0.0;
        double totalTransferWater = 0.0;
        double tmp_fraction = 0.0;
        patch[0].septicReleaseQ_m = patch[0].landuse_defaults[0][0].septic_water_load/patch[0].area; // release to patch
        
        for (i =0; i < patch[0].innundation_list[d].num_drainIN_septic; i++){
            // from surface to surface
            sourceTransferWater = min(
                              patch[0].innundation_list[d].drainIN_septic[i].maxDailyDrain,
                              patch[0].innundation_list[d].drainIN_septic[i].DrainFrac *
                                patch[0].landuse_defaults[0][0].septic_water_load/patch[0].area);// water depth release to patch
            
            
            patch[0].innundation_list[d].drainIN_septic[i].transfer_flux_surf = max(0.0, min(sourceTransferWater * patch[0].innundation_list[d].drainIN_septic[i].propDrainFrmSurf,
                                          patch[0].innundation_list[d].drainIN_septic[i].patch[0].detention_store));
            
            patch[0].detention_store += patch[0].innundation_list[d].drainIN_septic[i].transfer_flux_surf;
            
            tmp_fraction = min(1.0,(patch[0].innundation_list[d].drainIN_septic[i].patch[0].detention_store>0? patch[0].innundation_list[d].drainIN_septic[i].transfer_flux_surf/patch[0].innundation_list[d].drainIN_septic[i].patch[0].detention_store : 0.0));
            patch[0].surface_NO3 += tmp_fraction * patch[0].innundation_list[d].drainIN_septic[i].patch[0].surface_NO3;
            patch[0].surface_NH4 += tmp_fraction * patch[0].innundation_list[d].drainIN_septic[i].patch[0].surface_NH4;
            patch[0].surface_DON += tmp_fraction * patch[0].innundation_list[d].drainIN_septic[i].patch[0].surface_DON;
            patch[0].surface_DOC += tmp_fraction * patch[0].innundation_list[d].drainIN_septic[i].patch[0].surface_DOC;
            tmp_fraction = min(1.0,max(0.0, 1.0 - tmp_fraction));
            patch[0].innundation_list[d].drainIN_septic[i].patch[0].detention_store *= tmp_fraction;
            patch[0].innundation_list[d].drainIN_septic[i].patch[0].surface_NO3 *= tmp_fraction;
            patch[0].innundation_list[d].drainIN_septic[i].patch[0].surface_NH4 *= tmp_fraction;
            patch[0].innundation_list[d].drainIN_septic[i].patch[0].surface_DON *= tmp_fraction;
            patch[0].innundation_list[d].drainIN_septic[i].patch[0].surface_DOC *= tmp_fraction;
            
            // drawing from deep GW
            patch[0].innundation_list[d].drainIN_septic[i].transfer_flux_sub = max(0.0,min(sourceTransferWater * (1.0 - patch[0].innundation_list[d].drainIN_septic[i].propDrainFrmSurf), hillslope[0].gw.storage*hillslope[0].area/patch[0].area)); // actual available water from known source
            
            // perform water and solute transfer
            patch[0].detention_store += patch[0].innundation_list[d].drainIN_septic[i].transfer_flux_sub;
            tmp_fraction = min(1.0,(hillslope[0].gw.storage>0? patch[0].innundation_list[d].drainIN_septic[i].transfer_flux_sub*patch[0].area/hillslope[0].area/hillslope[0].gw.storage : 0.0));
            patch[0].surface_NO3 += tmp_fraction * hillslope[0].gw.NO3;
            patch[0].surface_NH4 += tmp_fraction * hillslope[0].gw.NH4;
            patch[0].surface_DON += tmp_fraction * hillslope[0].gw.DOC;
            patch[0].surface_DOC += tmp_fraction * hillslope[0].gw.DON;
            // substrate the source water and solute
            tmp_fraction = min(1.0,max(0.0, 1.0 - tmp_fraction));
            hillslope[0].gw.storage *= tmp_fraction;
            hillslope[0].gw.NO3 *= tmp_fraction;
            hillslope[0].gw.NH4 *= tmp_fraction;
            hillslope[0].gw.DOC *= tmp_fraction;
            hillslope[0].gw.DON *= tmp_fraction;
            
            totalTransferWater += patch[0].innundation_list[d].drainIN_septic[i].transfer_flux_sub + patch[0].innundation_list[d].drainIN_septic[i].transfer_flux_surf;
        }// for loop of sources
        if(patch[0].innundation_list[d].num_drainIN_septic>0){
            patch[0].septicReleaseQ_m = totalTransferWater;
            // one source per septic output
            patch[0].surface_NO3 += patch[0].landuse_defaults[0][0].septic_NO3_load/patch[0].area * patch[0].innundation_list[d].num_drainIN_septic;
        }
    }//if



	/*--------------------------------------------------------------*/
	/* if snowmelt is from another model (and input rather than computed */
	/* get that value and set it up to substitute for rhessys internal snowmelt */
	/*--------------------------------------------------------------*/
	snow_melt_input=-999.0;
	if (patch[0].base_stations != NULL) {
		inx = patch[0].base_stations[0][0].dated_input[0].snow_melt_input.inx;
		if (inx > -999) {
			clim_event = patch[0].base_stations[0][0].dated_input[0].snow_melt_input.seq[inx];
			while (julday(clim_event.edate) < julday(current_date)) {
				patch[0].base_stations[0][0].dated_input[0].snow_melt_input.inx += 1;
				inx = patch[0].base_stations[0][0].dated_input[0].snow_melt_input.inx;
				clim_event = patch[0].base_stations[0][0].dated_input[0].snow_melt_input.seq[inx];
				}
			if ((clim_event.edate.year != 0) && ( julday(clim_event.edate) == julday(current_date)) ) {
				snow_melt_input = clim_event.value;
				}
			else snow_melt_input = 0.0;
			} 
		else snow_melt_input=-999.0;
	}



	/*--------------------------------------------------------------*/
	/*	Compute the stability correction factor for aero cond	*/
	/*--------------------------------------------------------------*/
	patch[0].stability_correction = 1.0;
	
	
	/*--------------------------------------------------------------*/
	/*      Determine patch SOIL heat flux.                         */
	/*      (This is ignored if there is a 0 height stratum.        */
	/*--------------------------------------------------------------*/

	patch[0].surface_heat_flux = -1 * compute_surface_heat_flux(
		command_line[0].verbose_flag,
		patch[0].snow_stored,
		patch[0].unsat_storage,
		patch[0].sat_deficit,
		zone[0].metv.tavg,
		zone[0].metv.tnightmax,
		zone[0].metv.tsoil,
		patch[0].soil_defaults[0][0].deltaz,
		patch[0].soil_defaults[0][0].min_heat_capacity,
		patch[0].soil_defaults[0][0].max_heat_capacity);

	/*----------------------------------------------------------------------*/
	/*	Compute the no-canopy aerodynamic resistance at 2m as a baseline	*/
	/*	for null covers	*/
	/*----------------------------------------------------------------------*/
	tmpra = compute_ra_overstory(
								 command_line[0].verbose_flag,
								 0.0,
								 0.4,
								 &(tmpwind),
								 &(tmpwindcan),
								 &(tmpwindsnow),
								 &(tmpustar),
								 zone[0].base_stations[0][0].screen_height,
								 0.0,
								 0.0,
								 &(tmpga),
								 &(tmpgasnow));
	
	/* Set values for no stratum case. These will be overwritten if veg present. */
	patch[0].ga = tmpga;
	patch[0].gasnow = tmpgasnow;
	patch[0].wind = tmpwind;
	patch[0].windsnow = tmpwindsnow;
	patch[0].ustar = tmpustar;
	
	/*--------------------------------------------------------------*/
	/*	Cycle through patch layers with height greater than the	*/
	/*	snowpack.						*/
	/*--------------------------------------------------------------*/
	
	/*	Calculate initial pond height		*/
	pond_height = max(0.0, -patch[0].sat_deficit_z + patch[0].detention_store);
	
	/*--------------------------------------------------------------*/
	/* Layers above snowpack and pond */
	/*--------------------------------------------------------------*/
	for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
		patch[0].snowpack.overstory_height = zone[0].base_stations[0][0].screen_height;
		if ( (patch[0].layers[layer].height > patch[0].snowpack.height) && (patch[0].layers[layer].height > pond_height) ){
            
            // find a tall layer within a layer looping
                if ( command_line[0].verbose_flag == -5 ){ printf("\n     ABOVE SNOWPACK AND POND"); }
			patch[0].snowpack.overstory_fraction = max(patch[0].snowpack.overstory_fraction, (1.0 - patch[0].layers[layer].null_cover));
			patch[0].snowpack.overstory_height = max(patch[0].snowpack.overstory_height,  patch[0].layers[layer].height);
			patch[0].overstory_fraction = max(patch[0].overstory_fraction, (1.0 - patch[0].layers[layer].null_cover));
			
            patch[0].Kdown_direct_final = patch[0].layers[layer].null_cover * patch[0].Kdown_direct; //modified by canopy_stratum_dailyF() //<<--- zone[0].
            patch[0].Kdown_diffuse_final = patch[0].layers[layer].null_cover * patch[0].Kdown_diffuse; //modified by canopy_stratum_dailyF() //<<--- zone[0].
            // Kup_direct_final is 0 initially; calculated in canopy_stratum_dailyF() and accumulated
            // Kup_diffuse_final is 0 initially; calculated in canopy_stratum_dailyF() and accumulated
			patch[0].PAR_direct_final = patch[0].layers[layer].null_cover * patch[0].PAR_direct; //modified by canopy_stratum_dailyF() //<<--- zone[0].
			patch[0].PAR_diffuse_final = patch[0].layers[layer].null_cover * patch[0].PAR_diffuse; //modified by canopy_stratum_dailyF() //<<--- zone[0].
			
            patch[0].Ldown_final = patch[0].layers[layer].null_cover * patch[0].Ldown; //<<--- zone[0].Ldown
			patch[0].Ldown_final_night = patch[0].layers[layer].null_cover * patch[0].Ldown_night; //<<--- zone[0].Ldown_night
			patch[0].Ldown_final_day = patch[0].layers[layer].null_cover * patch[0].Ldown_day; //<<--- zone[0].Ldown_day
            
			patch[0].Kstar_canopy_final = patch[0].Kstar_canopy;
			patch[0].LE_canopy_final = patch[0].LE_canopy;
            
			patch[0].rain_throughfall_final = patch[0].layers[layer].null_cover * patch[0].rain_throughfall;
			patch[0].snow_throughfall_final = patch[0].layers[layer].null_cover * patch[0].snow_throughfall;
			patch[0].NO3_throughfall_final = patch[0].layers[layer].null_cover * patch[0].NO3_throughfall;
			patch[0].T_canopy_final = patch[0].layers[layer].null_cover * patch[0].T_canopy;
            
			if (dum == 0) {
				patch[0].ga_final = patch[0].layers[layer].null_cover * tmpga;
				patch[0].gasnow_final = patch[0].layers[layer].null_cover * tmpgasnow;
				patch[0].wind_final = patch[0].layers[layer].null_cover * tmpwind;
				patch[0].windsnow_final = patch[0].layers[layer].null_cover * tmpwindsnow;
				patch[0].ustar_final = patch[0].layers[layer].null_cover * tmpustar;
				if ( command_line[0].verbose_flag == -5 ){
					printf("\n     ***TOP: ga=%lf gasnow=%lf wind=%lf windsnow=%lf",patch[0].ga_final, patch[0].gasnow_final, patch[0].wind_final, patch[0].windsnow_final);
				}
			}
			else {				
				patch[0].ga_final = patch[0].layers[layer].null_cover * patch[0].ga;
				patch[0].gasnow_final = patch[0].layers[layer].null_cover * patch[0].gasnow;
				patch[0].wind_final = patch[0].layers[layer].null_cover * patch[0].wind;
				patch[0].windsnow_final = patch[0].layers[layer].null_cover * patch[0].windsnow;
				patch[0].ustar_final = patch[0].layers[layer].null_cover * patch[0].ustar;
				if ( command_line[0].verbose_flag == -5 ){
					printf("\n     ***NOT TOP: ga=%lf gasnow=%lf wind=%lf windsnow=%lf",patch[0].ga_final, patch[0].gasnow_final, patch[0].wind_final, patch[0].windsnow_final);
				}
			}

			/*--------------------------------------------------------------*/
			/*		Cycle through the canopy strata in this layer	*/
			/*--------------------------------------------------------------*/
            // plants within the same layer does not affect each other!
			for ( stratum=0 ; stratum<patch[0].layers[layer].count; stratum++ ){
					canopy_stratum_daily_F(
						world,
						basin,
						hillslope,
						zone,
						patch,
						&(patch[0].layers[layer]),
						patch[0].canopy_strata[(patch[0].layers[layer].strata[stratum])],
                        patch[0].shadow_strata[(patch[0].layers[layer].strata[stratum])],
						command_line,
						event,
						current_date );
				dum += 1;
			}
            
            if(patch[0].PAR_direct_final + patch[0].PAR_diffuse_final < 0.0){
                printf("patch dailyF abovesnow & pond (%d, %d, %d, %d) %lf,%lf -> %lf,%lf  %lf,%lf -> %lf,%lf  %lf,%lf -> %lf,%lf\n",
                       patch[0].ID, layer, layer, stratum,
                       patch[0].Kdown_direct, patch[0].Kup_direct, patch[0].Kdown_direct_final, patch[0].Kup_direct_final,
                       patch[0].Kdown_diffuse, patch[0].Kup_diffuse, patch[0].Kdown_diffuse_final, patch[0].Kup_diffuse_final,
                       patch[0].PAR_direct, patch[0].PAR_diffuse, patch[0].PAR_direct_final, patch[0].PAR_diffuse_final
                       );
            }//debug
            
			patch[0].Kdown_direct = patch[0].Kdown_direct_final; //modified by += canopy_stratum_dailyF() and update to patch
			patch[0].Kup_direct = patch[0].Kup_direct_final; //modified by += canopy_stratum_dailyF() and update to patch
			patch[0].Kdown_diffuse = patch[0].Kdown_diffuse_final; //modified by += canopy_stratum_dailyF() and update to patch
			patch[0].Kup_diffuse = patch[0].Kup_diffuse_final; //modified by += canopy_stratum_dailyF() and update to patch
			patch[0].PAR_direct = patch[0].PAR_direct_final; //modified by += canopy_stratum_dailyF() and update to patch
			patch[0].PAR_diffuse = patch[0].PAR_diffuse_final; //modified by += canopy_stratum_dailyF() and update to patch
            
			patch[0].Ldown = patch[0].Ldown_final;
			patch[0].Ldown_night = patch[0].Ldown_final_night;
			patch[0].Ldown_day = patch[0].Ldown_final_day;
			patch[0].Kstar_canopy = patch[0].Kstar_canopy_final;
			patch[0].LE_canopy = patch[0].LE_canopy_final;
			patch[0].rain_throughfall = patch[0].rain_throughfall_final;
			patch[0].snow_throughfall = patch[0].snow_throughfall_final;
			patch[0].NO3_throughfall = patch[0].NO3_throughfall_final;
			patch[0].ga = patch[0].ga_final;
			patch[0].gasnow = patch[0].gasnow_final;
			patch[0].wind = patch[0].wind_final;
			patch[0].windsnow = patch[0].windsnow_final;
			patch[0].ustar = patch[0].ustar_final;
			patch[0].T_canopy = patch[0].T_canopy_final;
		}
	}//layers
	
	/*--------------------------------------------------------------*/
	/*	Compute patch level long wave radiation processes.			*/
	/*--------------------------------------------------------------*/
	if (command_line[0].evap_use_longwave_flag) {
		compute_Lstar(command_line[0].verbose_flag,
					  basin,
					  zone,
					  patch); // what?
	}
	
	
	/*--------------------------------------------------------------*/
	/*	We assume the snowpack is conceptually over the		*/
	/*	current ponded water.					*/
	/*--------------------------------------------------------------*/
	/*--------------------------------------------------------------*/
	/*	Now add the throughfall of snow	to the snowpack	 	 		*/
	/*		rain is added to the snowpack if it exists				*/
	/*		and snowpack melt allowed to occur		 				*/
	/*		this means that in rain on snow - rain is included		*/
	/*		as snowmelt												*/
	/*--------------------------------------------------------------*/
	preday_snowpack_height = patch[0].snowpack.height;
	patch[0].snowpack.water_equivalent_depth += patch[0].snow_throughfall;
	
	patch[0].Kdown_direct_subcanopy = patch[0].Kdown_direct;
	patch[0].Kdown_diffuse_subcanopy = patch[0].Kdown_diffuse;
	
	if ( command_line[0].verbose_flag == -5 ){
	printf("\n     wind=%lf windfin=%lf windsnow=%lf SWE=%lf Kstarcan=%lf Kdowndirpch=%lf Kdowndifpch=%lf detstore=%lf T_canopy=%lf", 
			patch[0].wind, 
			patch[0].wind_final, 
			patch[0].windsnow, 
			patch[0].snowpack.water_equivalent_depth, 
			patch[0].Kstar_canopy/86.4, 
			patch[0].Kdown_direct/86.4, 
			patch[0].Kdown_diffuse/86.4,
		   patch[0].detention_store,
		   patch[0].T_canopy);
	}

	
	patch[0].Kdown_direct_bare = patch[0].Kdown_direct;
	patch[0].Kdown_diffuse_bare = patch[0].Kdown_diffuse;
	

	if ( patch[0].snowpack.water_equivalent_depth > ZERO ) {
		
		/*--------------------------------------------------------------*/
		/*	Calculate snowmelt 											*/
		/*--------------------------------------------------------------*/
		/*	Check to see if snowpack is above pond. If so, proceed 	*/
		/*	with snowpack daily F.  Otherwise, melt all snow and add 	*/
		/*	to rain throughfall.										*/
		/*--------------------------------------------------------------*/
		if ( patch[0].snowpack.water_equivalent_depth > pond_height ) {
			
			/* COVER FRACTION */
			if ((patch[0].snowpack.overstory_fraction < 1) && (patch[0].snowpack.overstory_fraction > 0)) {
				if ( command_line[0].verbose_flag == -5 ){
					printf("\nSNOWPACK WITH COVER FRACTION %lf", 
						   patch[0].snowpack.overstory_fraction);
				}				
				/* Separate Kdown under canopy from patch-average Kdown below canopy layers */
				/* Using zone Kdown for exposed portion */
				Kdown_direct_covered = ( patch[0].Kdown_direct - zone[0].Kdown_direct 
											* (1 - patch[0].snowpack.overstory_fraction) )
											/ patch[0].snowpack.overstory_fraction;
				Kdown_diffuse_covered = ( patch[0].Kdown_diffuse - zone[0].Kdown_diffuse 
											* (1 - patch[0].snowpack.overstory_fraction) )
											/ patch[0].snowpack.overstory_fraction;
				Kdown_direct_exposed = zone[0].Kdown_direct;
				Kdown_diffuse_exposed = zone[0].Kdown_diffuse;
				PAR_direct_covered = ( patch[0].PAR_direct - zone[0].PAR_direct 
											* (1 - patch[0].snowpack.overstory_fraction) )
											/ patch[0].snowpack.overstory_fraction;
				PAR_diffuse_covered = ( patch[0].PAR_diffuse - zone[0].PAR_diffuse 
											* (1 - patch[0].snowpack.overstory_fraction) )
											/ patch[0].snowpack.overstory_fraction;
				PAR_direct_exposed = zone[0].PAR_direct;
				PAR_diffuse_exposed = zone[0].PAR_diffuse;
				
				/* Lundberg 1994 reduce conductance for snow vs. rain by factor of 10 */
				snow_melt_covered = snowpack_daily_F(
						current_date,
						command_line[0].verbose_flag,
						zone,
						patch,
						&patch[0].snowpack,
						basin[0].theta_noon,
						zone[0].metv.tavg,
						zone[0].e_dewpoint,
						patch[0].gasnow/10.0,
						zone[0].metv.pa,
						zone[0].cloud_fraction,
						patch[0].rain_throughfall,
						patch[0].snow_throughfall,
						&Kdown_direct_covered,
						&Kup_direct_snow_covered,
						&Kdown_diffuse_covered,
						&Kup_diffuse_snow_covered,
						&PAR_direct_covered,
						&PAR_diffuse_covered,
						patch[0].soil_defaults[0][0].maximum_snow_energy_deficit,
						patch[0].soil_defaults[0][0].snow_water_capacity,
						patch[0].soil_defaults[0][0].snow_light_ext_coef,
						patch[0].soil_defaults[0][0].snow_melt_Tcoef,
						1.0,
						patch[0].snowpack.overstory_fraction,
						0);
				snow_melt_exposed = snowpack_daily_F(
						current_date,
						command_line[0].verbose_flag,
						zone,
						patch,
						&patch[0].snowpack,
						basin[0].theta_noon,
						zone[0].metv.tavg,
						zone[0].e_dewpoint,
						patch[0].gasnow/10.0,
						zone[0].metv.pa,
						zone[0].cloud_fraction,
						patch[0].rain_throughfall,
						patch[0].snow_throughfall,
						&Kdown_direct_exposed,
						&Kup_direct_snow_exposed,
						&Kdown_diffuse_exposed,
						&Kup_diffuse_snow_exposed,
						&PAR_direct_exposed,
						&PAR_diffuse_exposed,
						patch[0].soil_defaults[0][0].maximum_snow_energy_deficit,
						patch[0].soil_defaults[0][0].snow_water_capacity,
						patch[0].soil_defaults[0][0].snow_light_ext_coef,
						patch[0].soil_defaults[0][0].snow_melt_Tcoef,
						0.0,
						(1-patch[0].snowpack.overstory_fraction),
						1);
				
				patch[0].snow_melt = (snow_melt_covered * patch[0].snowpack.overstory_fraction)
						+ (snow_melt_exposed * (1-patch[0].snowpack.overstory_fraction));
				patch[0].Kdown_direct = (Kdown_direct_covered * patch[0].snowpack.overstory_fraction)
						+ (Kdown_direct_exposed * (1-patch[0].snowpack.overstory_fraction));
				patch[0].Kdown_diffuse = (Kdown_diffuse_covered * patch[0].snowpack.overstory_fraction)
						+ (Kdown_diffuse_exposed * (1-patch[0].snowpack.overstory_fraction));
				patch[0].PAR_direct = (PAR_direct_covered * patch[0].snowpack.overstory_fraction)
						+ (PAR_direct_exposed * (1-patch[0].snowpack.overstory_fraction));
				patch[0].PAR_diffuse = (PAR_diffuse_covered * patch[0].snowpack.overstory_fraction)
						+ (PAR_diffuse_exposed * (1-patch[0].snowpack.overstory_fraction));
				
			}
			/* NO COVER FRACTION */
			else {
				if ( command_line[0].verbose_flag == -5 ){
					printf("\nSNOWPACK WITHOUT COVER FRACTION %lf", 
						   patch[0].snowpack.overstory_fraction);
				}								
				patch[0].snow_melt = snowpack_daily_F(
					current_date,
					command_line[0].verbose_flag,
					zone,
					patch,
					&patch[0].snowpack,
					basin[0].theta_noon,
					zone[0].metv.tavg,
					zone[0].e_dewpoint,
					patch[0].gasnow/10.0,
					zone[0].metv.pa,
					zone[0].cloud_fraction,
					patch[0].rain_throughfall,
					patch[0].snow_throughfall,
					&patch[0].Kdown_direct,
					&Kup_direct_snow,
					&patch[0].Kdown_diffuse,
					&Kup_diffuse_snow,
					&patch[0].PAR_direct,
					&patch[0].PAR_diffuse,
					patch[0].soil_defaults[0][0].maximum_snow_energy_deficit,
					patch[0].soil_defaults[0][0].snow_water_capacity,
					patch[0].soil_defaults[0][0].snow_light_ext_coef,
					patch[0].soil_defaults[0][0].snow_melt_Tcoef,
					patch[0].snowpack.overstory_fraction,
					1.0,
					1);
					if (patch[0].snowpack.overstory_fraction == 0) {
						Kup_direct_snow_exposed = Kup_direct_snow;
						Kup_diffuse_snow_exposed = Kup_diffuse_snow;
						}
				}
				
			/* FOR ALL COVER FRACTIONS */
			patch[0].Kup_direct += Kup_direct_snow_exposed * (1 - patch[0].snowpack.overstory_fraction);
			patch[0].Kup_diffuse += Kup_diffuse_snow_exposed * (1 - patch[0].snowpack.overstory_fraction);
		
			patch[0].snowpack.water_equivalent_depth -= patch[0].snow_melt;
            patch[0].snowpack.water_equivalent_depth = max(0.0, patch[0].snowpack.water_equivalent_depth);
			patch[0].snowpack.sublimation = min(patch[0].snowpack.sublimation, patch[0].snowpack.water_equivalent_depth);
			patch[0].snowpack.height = patch[0].snowpack.water_equivalent_depth / 0.1; /* snow density ~ 0.1 */

			if (snow_melt_input == -999.0) 
				patch[0].rain_throughfall += patch[0].snow_melt;
			else {
				patch[0].rain_throughfall += snow_melt_input;
				patch[0].snow_melt = snow_melt_input;
			}
			patch[0].snow_throughfall = 0.0;
			patch[0].snowpack.water_equivalent_depth -= patch[0].snowpack.sublimation;
            patch[0].snowpack.water_equivalent_depth = max(0.0, patch[0].snowpack.water_equivalent_depth);
			/* Force turbulent fluxes to 0 under snowpack */
			patch[0].ga = 0.0;
			patch[0].wind = 0.0;
		} else {
			patch[0].rain_throughfall += patch[0].snowpack.water_equivalent_depth;
			patch[0].snow_throughfall = 0.0;
			patch[0].snowpack.water_equivalent_depth = 0.0;
			patch[0].snowpack.height = 0.0;
		}
	}else{
		/*--------------------------------------------------------------*/
		/*	Just to create symmetrical output for snow and no snow	*/
		/*	days we do some fake calls which snowpack_daily_F does	*/
		/*--------------------------------------------------------------*/
		temp = 1;
		temp = compute_radiative_fluxes(0,&temp,&temp,1,1,1);
		temp = compute_radiative_fluxes(0,&temp,&temp,1,1,1);
		temp = compute_radiative_fluxes(0,&temp,&temp,1,1,1);
		temp = compute_radiative_fluxes(0,&temp,&temp,1,1,1);
		patch[0].snow_melt = 0.0;
		patch[0].snowpack.energy_deficit = 0.001;
		patch[0].snowpack.Kstar_direct = 0.0;
		patch[0].snowpack.Kstar_diffuse = 0.0;
		patch[0].snowpack.APAR_direct = 0.0;
		patch[0].snowpack.APAR_diffuse = 0.0;
		patch[0].snowpack.water_equivalent_depth = 0.0;
	}
	
	if (patch[0].snowpack.water_equivalent_depth < 0.0001) {
		patch[0].rain_throughfall += patch[0].snowpack.water_equivalent_depth;
		patch[0].snowpack.water_equivalent_depth = 0.0;
		patch[0].snowpack.energy_deficit = 0.001;
		patch[0].snowpack.surface_age = 0.0;
		patch[0].snowpack.T = 0.0;
		patch[0].snowpack.height = 0.0;
		}
	
	if ( command_line[0].verbose_flag == -5 ){
		printf("\n     AFTER SNOWPACK: Kup_direct=%lf Kup_diffuse=%lf", 
			   patch[0].Kup_direct/86.4, 
			   patch[0].Kup_diffuse/86.4);
	}
	
	
	/*--------------------------------------------------------------*/
	/*	Cycle through patch layers with height less than the	*/
	/*	snowpack but more than  0				*/
	/*	Note that the rain throughfall should equal the*/
	/*	rain and snow melt getting through the snowpack.	*/
	/*	There should be no snow_throughfall if there is a snow	*/
	/*	pack (i.e. only moss layers will have any incident	*/
	/*	snow throughfall in the loop below)			*/
	/*	We then conceptually look at it as a snowpack 		*/
	/*	"HOVERING" above any patch layers lower than its current*/
	/*	maximum height.  This is fine since if we assume that	*/
	/*	The snowpack does not transmit shortwave or par		*/
	/*	no vapour fluxes will take place anyways and only	*/
	/*	respiration will occurr which we want .			*/
	/*								*/
	/*	Patches under the snowpack but over the pond.		*/
	/*	need to use previous day (or beginning of the day)	*/
	/*	snowpack height to avoid processing some strata		*/
	/*	twice							*/
	/*--------------------------------------------------------------*/
	/* Layers below snowpack and above pond */
	/*--------------------------------------------------------------*/
	for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
		if ( (preday_snowpack_height > 0.0) && (patch[0].layers[layer].height <= preday_snowpack_height) &&
			(patch[0].layers[layer].height > pond_height) ){
			if ( command_line[0].verbose_flag == -5 ){
				printf("\n     BELOW SNOWPACK AND ABOVE POND");
			}
			patch[0].Tday_surface_offset_final = 0.0;
			patch[0].Kdown_direct_final = patch[0].layers[layer].null_cover * patch[0].Kdown_direct;
			patch[0].Kdown_diffuse_final = patch[0].layers[layer].null_cover * patch[0].Kdown_diffuse;
			patch[0].PAR_direct_final = patch[0].layers[layer].null_cover * patch[0].PAR_direct;
			patch[0].PAR_diffuse_final = patch[0].layers[layer].null_cover * patch[0].PAR_diffuse;
			patch[0].rain_throughfall_final = patch[0].layers[layer].null_cover * patch[0].rain_throughfall;
			patch[0].snow_throughfall_final = patch[0].layers[layer].null_cover * patch[0].snow_throughfall;
			patch[0].NO3_throughfall_final = patch[0].layers[layer].null_cover * patch[0].NO3_throughfall;
			patch[0].ga_final = patch[0].layers[layer].null_cover * patch[0].ga;
			patch[0].wind_final = patch[0].layers[layer].null_cover * patch[0].wind;
			patch[0].T_canopy_final = patch[0].layers[layer].null_cover * patch[0].T_canopy;
			for ( stratum=0 ;stratum<patch[0].layers[layer].count; stratum++ ){
					canopy_stratum_daily_F(
						world,
						basin,
						hillslope,
						zone,
						patch,
						&(patch[0].layers[layer]),
						patch[0].canopy_strata[(patch[0].layers[layer].strata[stratum])],
						patch[0].shadow_strata[(patch[0].layers[layer].strata[stratum])],
						command_line,
						event,
						current_date );
			}
            
//            if(patch[0].PAR_direct_final<=ZERO || patch[0].PAR_diffuse_final<=ZERO){
//                printf("patch dailyF under snow but above pond (%d, %d, %d) %lf,%lf -> %lf,%lf  %lf,%lf -> %lf,%lf  %lf,%lf -> %lf,%lf\n",
//                       patch[0].ID, layer, patch[0].canopy_strata[(patch[0].layers[layer].strata[stratum])]->ID,
//                       patch[0].Kdown_direct, patch[0].Kup_direct, patch[0].Kdown_direct_final, patch[0].Kup_direct_final,
//                       patch[0].Kdown_diffuse, patch[0].Kup_diffuse, patch[0].Kdown_diffuse_final, patch[0].Kup_diffuse_final,
//                       patch[0].PAR_direct, patch[0].PAR_diffuse, patch[0].PAR_direct_final, patch[0].PAR_diffuse_final
//                       );
//            }//debug
            
			patch[0].Kdown_direct = patch[0].Kdown_direct_final;
			patch[0].Kdown_diffuse = patch[0].Kdown_diffuse_final;
			patch[0].PAR_direct = patch[0].PAR_direct_final;
			patch[0].PAR_diffuse = patch[0].PAR_diffuse_final;
			patch[0].rain_throughfall = patch[0].rain_throughfall_final;
			patch[0].snow_throughfall = patch[0].snow_throughfall_final;
			patch[0].NO3_throughfall = patch[0].NO3_throughfall_final;
			patch[0].Tday_surface_offset = patch[0].Tday_surface_offset_final;
			patch[0].ga = patch[0].ga_final;
			patch[0].wind = patch[0].wind_final;
			patch[0].ustar = patch[0].ustar_final;
			patch[0].T_canopy = patch[0].T_canopy_final;
		}// height if
	}//layer loop
		
	/*--------------------------------------------------------------*/
	/*	Layers below the pond.					*/
	/*--------------------------------------------------------------*/
	for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
		if (patch[0].layers[layer].height <= pond_height){
			if ( command_line[0].verbose_flag == -5 ){
				printf("\n     BELOW POND (INCLUDING SURFACE)");
			}
			patch[0].Tday_surface_offset_final = 0.0;
			patch[0].Kdown_direct_final = patch[0].layers[layer].null_cover * patch[0].Kdown_direct;
			patch[0].Kdown_diffuse_final = patch[0].layers[layer].null_cover * patch[0].Kdown_diffuse;
			patch[0].PAR_direct_final = patch[0].layers[layer].null_cover * patch[0].PAR_direct;
			patch[0].PAR_diffuse_final = patch[0].layers[layer].null_cover * patch[0].PAR_diffuse;
			patch[0].rain_throughfall_final = patch[0].layers[layer].null_cover * patch[0].rain_throughfall;
			patch[0].snow_throughfall_final = patch[0].layers[layer].null_cover * patch[0].snow_throughfall;
			patch[0].NO3_throughfall_final = patch[0].layers[layer].null_cover * patch[0].NO3_throughfall;
			patch[0].ga_final = patch[0].layers[layer].null_cover * patch[0].ga;
			patch[0].wind_final = patch[0].layers[layer].null_cover * patch[0].wind;
			patch[0].T_canopy_final = patch[0].layers[layer].null_cover * patch[0].T_canopy;
			for ( stratum=0 ; stratum<patch[0].layers[layer].count; stratum++ ){
					canopy_stratum_daily_F(
						world,
						basin,
						hillslope,
						zone,
						patch,
						&(patch[0].layers[layer]),
						patch[0].canopy_strata[(patch[0].layers[layer].strata[stratum])],
						patch[0].shadow_strata[(patch[0].layers[layer].strata[stratum])],
						command_line,
						event,
						current_date );
			}
            
//            if(patch[0].PAR_direct_final<=ZERO || patch[0].PAR_diffuse_final<=ZERO){
//                printf("patch dailyF under snow and under pond (%d, %d, %d) %lf,%lf -> %lf,%lf  %lf,%lf -> %lf,%lf  %lf,%lf -> %lf,%lf\n",
//                       patch[0].ID, layer, patch[0].canopy_strata[(patch[0].layers[layer].strata[stratum])]->ID,
//                       patch[0].Kdown_direct, patch[0].Kup_direct, patch[0].Kdown_direct_final, patch[0].Kup_direct_final,
//                       patch[0].Kdown_diffuse, patch[0].Kup_diffuse, patch[0].Kdown_diffuse_final, patch[0].Kup_diffuse_final,
//                       patch[0].PAR_direct, patch[0].PAR_diffuse, patch[0].PAR_direct_final, patch[0].PAR_diffuse_final
//                       );
//            }//debug
            
			patch[0].Kdown_direct = patch[0].Kdown_direct_final;
			patch[0].Kdown_diffuse = patch[0].Kdown_diffuse_final;
			patch[0].PAR_direct = patch[0].PAR_direct_final;
			patch[0].PAR_diffuse = patch[0].PAR_diffuse_final;
			patch[0].rain_throughfall = patch[0].rain_throughfall_final;
			patch[0].snow_throughfall = patch[0].snow_throughfall_final;
			patch[0].NO3_throughfall = patch[0].NO3_throughfall_final;
			patch[0].Tday_surface_offset = patch[0].Tday_surface_offset_final;
			patch[0].ga = patch[0].ga_final;
			patch[0].wind = patch[0].wind_final;
			patch[0].ustar = patch[0].ustar_final;
			patch[0].T_canopy = patch[0].T_canopy_final;
		}//heigh if
	}//layer loop
	
	if ( command_line[0].verbose_flag == -5 ){
		printf("\n     PATCH DAILY POST LAYERS: ga=%lf Kdowndirpch=%lf Kdowndiffpch=%lf rainthru=%lf snowthru=%lf wind=%lf ustar=%lf Tcan=%lf", 
			   patch[0].ga,patch[0].Kdown_direct/86.4, patch[0].Kdown_diffuse/86.4, 
			   patch[0].rain_throughfall, patch[0].snow_throughfall,
			   patch[0].wind, patch[0].ustar,
			   patch[0].T_canopy);
	}

	/*--------------------------------------------------------------*/
	/*	Need to account for snow "throughfall" that is a result	*/
	/*	of leaves dropping, reducing potential interception, 	*/
	/*	and thus reducing snow storage.	This should be small	*/
	/*	enough (and generally in the fall) to ignore melting	*/
	/*	of this snow.						*/
	/*--------------------------------------------------------------*/
	patch[0].snowpack.water_equivalent_depth += patch[0].snow_throughfall;

	/*--------------------------------------------------------------*/
	/*	Do the pond energy balance now.				*/
	/*	Currently we just throw the pond in as rain.		*/
	/* 	before adding in detention store, determine input	*/
	/*	due to atm N wet dep concentration			*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/*      added in nitrogen deposition                            */
	/*	will consider N deposition as a concentration		*/
	/*	in throughfall - at moment				*/
	/*	simply added N deposition as a total flux		*/
	/*--------------------------------------------------------------*/


	if (patch[0].base_stations != NULL) {
		inx = patch[0].base_stations[0][0].dated_input[0].fertilizer_NO3.inx;
		if (inx > -999) {
			clim_event = patch[0].base_stations[0][0].dated_input[0].fertilizer_NO3.seq[inx];
			while (julday(clim_event.edate) < julday(current_date)) {
				patch[0].base_stations[0][0].dated_input[0].fertilizer_NO3.inx += 1;
				inx = patch[0].base_stations[0][0].dated_input[0].fertilizer_NO3.inx;
				clim_event = patch[0].base_stations[0][0].dated_input[0].fertilizer_NO3.seq[inx];
				}
			if ((clim_event.edate.year != 0) && ( julday(clim_event.edate) == julday(current_date)) ) {
				fertilizer_NO3 = clim_event.value;
				}
			else fertilizer_NO3 = 0.0;
			} 
		else fertilizer_NO3 = patch[0].landuse_defaults[0][0].fertilizer_NO3;
		}
	else fertilizer_NO3 = patch[0].landuse_defaults[0][0].fertilizer_NO3;

	if (patch[0].base_stations != NULL) {
		inx = patch[0].base_stations[0][0].dated_input[0].fertilizer_NH4.inx;
		if (inx > -999) {
			clim_event = patch[0].base_stations[0][0].dated_input[0].fertilizer_NH4.seq[inx];
			while (julday(clim_event.edate) < julday(current_date)) {
				patch[0].base_stations[0][0].dated_input[0].fertilizer_NH4.inx += 1;
				inx = patch[0].base_stations[0][0].dated_input[0].fertilizer_NH4.inx;
				clim_event = patch[0].base_stations[0][0].dated_input[0].fertilizer_NH4.seq[inx];
				}
			if ((clim_event.edate.year != 0) && ( julday(clim_event.edate) == julday(current_date)) ) {
				fertilizer_NH4 = clim_event.value;
				}
			else fertilizer_NH4 = 0.0;
			} 
		else fertilizer_NH4 = patch[0].landuse_defaults[0][0].fertilizer_NH4;
		}
	else fertilizer_NH4 = patch[0].landuse_defaults[0][0].fertilizer_NH4;

	/*
	if (patch[0].drainage_type == STREAM) {
		fertilizer_NO3 = 0.0;
		fertilizer_NH4 = 0.0;
		}
	*/

    // patch[0].fertilizer_NO3 is accumulating fertilizer_NO3 that has been added to the patch. nothing consumes it!
//    patch[0].fertilizer_NO3 += fertilizer_NO3;
//    patch[0].fertilizer_NH4 += fertilizer_NH4;
	//patch[0].surface_NO3 += zone[0].ndep_NO3; // old code
    patch[0].surface_NO3 += 0.5 * patch[0].NO3_throughfall; // stratum[0].NO3_stored is one in patch_hour() and stratum_hourly()
	patch[0].surface_NH4 += zone[0].ndep_NH4;

	/*--------------------------------------------------------------*/
	/*	a certain amount of surface_N is incorporated into the */
	/*	soil each day - we used 66% based on fertilizer experiments 	*/
	/*	Agronomy Guide 1989-1990 pg 17, Penn State web site fact sheet */
	/*--------------------------------------------------------------*/
//    FERT_TO_SOIL = 100.0;// this is 100, not 0.66 or 60%
//    additionally, patch[0].fertilizer_NO3 and patch[0].fertilizer_NH4 are accumulating the inputs and not even a bit of reduce!
//    it's also adding to soil Npool but not to surface N pool and not count for cover fraction
//    the code makes no sense!
//    if (patch[0].fertilizer_NH4 > ZERO) {
//        surfaceN_to_soil = FERT_TO_SOIL * patch[0].fertilizer_NH4;
//        patch[0].fertilizer_NH4 -= surfaceN_to_soil;
//        patch[0].soil_ns.sminn += surfaceN_to_soil;
//        }
//
//    if (patch[0].fertilizer_NO3 > ZERO) {
//        surfaceN_to_soil = FERT_TO_SOIL * patch[0].fertilizer_NO3;
//        patch[0].fertilizer_NO3 -= surfaceN_to_soil;
//        patch[0].soil_ns.nitrate += surfaceN_to_soil;
//        }
    // new fertilizer code
    if(command_line[0].fertilizer_flag==1 && patch[0].drainage_type>0 && patch[0].drainage_type % actionFERTILIZE==0){
        
        int fertilizerAdded = 0;
        for ( j=0 ; j<patch[0].num_canopy_strata ; j++ ){
            
            if(patch[0].canopy_strata[j][0].defaults[0][0].epc.veg_type == GRASS && patch[0].canopy_strata[j][0].phen.gwseasonday>0){
                //gwseasonday > 0 when plant is growing; its -1 otherwise.
                // fertilizer_NO3/fertilizer_NH4 is areal rate (gN/m2) set in LULC def, not from time series.
                if(patch[0].fertilizerDaysCount % patch[0].landuse_defaults[0][0].fertilizer_freq == 0){
                    patch[0].stored_fertilizer_NO3 += fertilizer_NO3 * patch[0].canopy_strata[j][0].cover_fraction;
                    patch[0].stored_fertilizer_NH4 += fertilizer_NH4 * patch[0].canopy_strata[j][0].cover_fraction;
                }//fertilizer
                fertilizerAdded++;
            }
        }// for loop j
        if(fertilizerAdded>0) patch[0].fertilizerDaysCount++;
        else patch[0].fertilizerDaysCount = 0;
        
        patch[0].soil_ns.nitrate +=  0.2 * patch[0].stored_fertilizer_NO3 * patch[0].landuse_defaults[0][0].fertilizer_decay_rate;
        patch[0].surface_NO3 += 0.8 * patch[0].stored_fertilizer_NO3 * patch[0].landuse_defaults[0][0].fertilizer_decay_rate;
        patch[0].stored_fertilizer_NO3 *= 1.0 - patch[0].landuse_defaults[0][0].fertilizer_decay_rate; // remaining
        patch[0].soil_ns.sminn += 0.2 * patch[0].stored_fertilizer_NH4 * patch[0].landuse_defaults[0][0].fertilizer_decay_rate;
        patch[0].surface_NH4 += 0.8 * patch[0].stored_fertilizer_NH4 * patch[0].landuse_defaults[0][0].fertilizer_decay_rate;
        patch[0].stored_fertilizer_NH4 *= 1.0 - patch[0].landuse_defaults[0][0].fertilizer_decay_rate; // remaining
    }// fertilizer_flag
		
	/*--------------------------------------------------------------*/
	/* adjust PH using data patch level inputs			*/
	/*--------------------------------------------------------------*/
	if (patch[0].base_stations != NULL) {
		inx = patch[0].base_stations[0][0].dated_input[0].PH.inx;
		if (inx > -999) {
			clim_event = patch[0].base_stations[0][0].dated_input[0].PH.seq[inx];
			while (julday(clim_event.edate) < julday(current_date)) {
				patch[0].base_stations[0][0].dated_input[0].PH.inx += 1;
				inx = patch[0].base_stations[0][0].dated_input[0].PH.inx;
				clim_event = patch[0].base_stations[0][0].dated_input[0].PH.seq[inx];
            }//while
			if ((clim_event.edate.year != 0) && ( julday(clim_event.edate) == julday(current_date)) ) {
				patch[0].PH = clim_event.value;
            }//if
        }//if
    }//if

	
	/*	Add rain throughfall to detention store for infiltration	*/
	/*	and evaporation routines.									*/
	
	patch[0].detention_store += 0.5 * patch[0].rain_throughfall;
	patch[0].surface_NO3 += 0.5 * patch[0].NO3_throughfall; // stratum[0].NO3_stored is one in patch_hour() and stratum_hourly()
    if (zone[0].hourly_rain_flag!=1) {patch[0].surface_NO3 += (command_line[0].fracDirectNdep) * zone[0].ndep_NO3;}

	/* Calculate det store, litter, and bare soil evap first */
	
	surface_daily_F(
					world,
					basin,
					hillslope,
					zone,
					patch,
					command_line,
					event,
					current_date );	
	
	if ( command_line[0].verbose_flag == -5 ){
		printf("\n     AFTER SURFACE: Kup_direct=%lf Kup_diffuse=%lf", 
			   patch[0].Kup_direct/86.4, 
			   patch[0].Kup_diffuse/86.4);
	}
	
	patch[0].detention_store += 0.5*patch[0].rain_throughfall + irrigation + patch[0].grassIrrigation_m;
	
    if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate ||
    patch[0].soil_ns.nitrate<0 ||
    patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn ||
    patch[0].soil_ns.sminn<0 ||
    patch[0].soil_ns.DON!=patch[0].soil_ns.DON ||
    patch[0].soil_ns.DON<0 ||
    patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC ||
    patch[0].soil_cs.DOC<0 ||
    patch[0].sat_NO3!=patch[0].sat_NO3 ||
    patch[0].sat_NO3<0 || patch[0].sat_NH4!=patch[0].sat_NH4 ||
    patch[0].sat_NH4<0 ) printf("patch daily F3N %d-%d-%d [%d,%d,%d,%d]{%e,%e,%e,%e}[%e,%e,%e,%e]\n",
       current_date.year, current_date.month, current_date.day,
       patch[0].ID, patch[0].drainage_type, actionPIPEDRAIN, actionSEWER,
       patch[0].soil_ns.nitrate,
       patch[0].soil_ns.sminn,
       patch[0].soil_cs.DOC,
       patch[0].soil_ns.DON,
       patch[0].sat_NO3, patch[0].sat_NH4, patch[0].sat_DOC, patch[0].sat_DON);
	/*--------------------------------------------------------------*/
	/* if there is hourly rain input, don't run the daily infiltration	*/
	/*--------------------------------------------------------------*/
	if (zone[0].hourly_rain_flag!=1) {
		/*--------------------------------------------------------------*/
		/* 	Above ground Hydrologic Processes			*/
		/* 	compute infiltration into the soil			*/
		/*	from snowmelt or rain_throughfall			*/
		/*	for now assume that all water infilatrates		*/
		/*	and some may be lost to deeper groundwater		*/
		/*--------------------------------------------------------------*/
		net_inflow = 0.0;
		duration = 0.0;
		if (patch[0].detention_store > ZERO) {
            
			/*------------------------------------------------------------------------*/
			/*	drainage to a deeper groundwater store				  */
			/*	move both nitrogen and water				    	*/
			/*------------------------------------------------------------------------*/
			if (command_line[0].gw_flag > 0 && ((patch[0].drainage_type>0 && patch[0].drainage_type % actionGWDRAIN==0) || patch[0].drainage_type==ROAD)) {
                
                // these patch[0].gw_drainage.XX are calculated hourly in subsurface_routing -> land_drainage ()
                hillslope[0].gw.storage += patch[0].gw_drainage / hillslope[0].area;
                hillslope[0].gw.DON += patch[0].gw_drainage_DON / hillslope[0].area;
                hillslope[0].gw.DOC += patch[0].gw_drainage_DOC / hillslope[0].area;
                hillslope[0].gw.NH4 += patch[0].gw_drainage_NH4 / hillslope[0].area;
                hillslope[0].gw.NO3 += patch[0].gw_drainage_NO3 / hillslope[0].area;
                patch[0].gw_drainage /= patch[0].area; // for water balance below;
                // patch[0].gw_drainage and patch[0].gw_drainage_DON are from yesterday, right
                
                // diffusion
//                if(patch[0].soil_defaults[0][0].sat_to_gw_coeff>0 && patch[0].available_soil_water>0 && hillslope[0].gw.storage>0){
//                    // hillslope[0].gw.NO3*hillslope[0].gw.soluteConc0coef is the [N0]
//                    patch[0].gw_diffuse = hillslope[0].defaults[0][0].gw_diffusion_coef * (patch[0].sat_NO3 - hillslope[0].gw.NO3*hillslope[0].gw.soluteConc0coef * patch[0].available_soil_water); // a flux
//                    patch[0].gw_diffuse += hillslope[0].defaults[0][0].gw_diffusion_coef * (patch[0].sat_NO3 - hillslope[0].gw.NO3*patch[0].available_soil_water/hillslope[0].gw.storage); // a flux
//                    patch[0].gw_diffuse *= 0.5;
//
//                    patch[0].gw_diffuse = max(min(patch[0].gw_diffuse, patch[0].sat_NO3), -hillslope[0].gw.NO3*hillslope[0].area/patch[0].area);
//                    //if(solute2gw>0) // bounded by patch[0].sat_NO3
//                    //if(solute2gw<0) // bounded by hillslope[0].gw.NO3*hillslope[0].area/patch[0].area
//
//                    patch[0].sat_NO3 -= patch[0].gw_diffuse;
//                    hillslope[0].gw.NO3 += patch[0].gw_diffuse * patch[0].area/hillslope[0].area;
//                    if(patch[0].sat_NO3<0){
//                        if(patch[0].sat_NO3 < -1e-8) printf("sat-gw diffusion@%d %f %f %f %f\n",patch[0].ID, patch[0].sat_NO3, hillslope[0].gw.NO3, patch[0].gw_diffuse, patch[0].area/hillslope[0].area);
//                        patch[0].sat_NO3 = 0.0;
//                    }
//                    if(hillslope[0].gw.NO3<0){
//                        if(hillslope[0].gw.NO3 < -1e-8) printf("sat-gw diffusion@%d %f %f %f %f\n",patch[0].ID, patch[0].sat_NO3, hillslope[0].gw.NO3, patch[0].gw_diffuse, patch[0].area/hillslope[0].area);
//                        hillslope[0].gw.NO3 = 0.0;
//                    }
//                }// end of if
                
                if ( update_gw_drainage(patch,
					hillslope,
					zone,
					command_line,
					current_date) != 0) {
					fprintf(stderr,"fATAL ERROR: in update_gw_drainage() ... Exiting\n");
					exit(EXIT_FAILURE);
				}
			}
			net_inflow = patch[0].detention_store;
			/*--------------------------------------------------------------*/
			/*      - if rain duration is zero, then input is from snow     */
			/*      melt  assume full daytime duration                      */
			/*--------------------------------------------------------------*/
            // this is the first patch[0].S call in patch_daily_F
            // in patch_daily_I.c, patch[0].S = (patch[0].rz_storage+patch[0].unsat_storage)/patch[0].sat_deficit;
            // 1/86400 = 1.15741e-05
            duration = 0.00001157407 * (zone[0].rain_duration <= ZERO? zone[0].metv.dayl : zone[0].rain_duration ); //makes rain all daytime
            infiltration = compute_infiltration(
                command_line[0].verbose_flag,
                patch[0].sat_deficit_z,
                0.0, //patch[0].aboveWT_SatPct,
                patch[0].Ksat_vertical,
                patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index],
                patch[0].rz_storage+patch[0].unsat_storage,
                patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].sat_def_0zm[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].sat_def_0zm[patch[0].sat_def_pct_index],
                patch[0].sat_deficit,
                net_inflow,
                duration,
                patch[0].soil_defaults[0][0].psi_air_entry);
            
        } else infiltration = 0.0;

		if (infiltration < 0.0) {
			printf("\nInfiltration %lf < 0 for %d on %ld",
				infiltration,
				patch[0].ID, current_date.day);
		}//if
		/*--------------------------------------------------------------*/
		/* determine fate of hold infiltration excess in detention store */
		/* infiltration excess will removed during routing portion	*/
		/*--------------------------------------------------------------*/
		
		infiltration = min(infiltration,patch[0].detention_store);

		/*--------------------------------------------------------------*/
		/* now take infiltration out of detention store 	*/
		/*--------------------------------------------------------------*/
		
		patch[0].detention_store -= infiltration;
			/*--------------------------------------------------------------*/
			/*	Determine if the infifltration will fill up the unsat	*/
			/*	zone or not.						*/
			/*	We use the strict assumption that sat deficit is the	*/
			/*	amount of water needed to saturate the soil.		*/
			/*--------------------------------------------------------------*/
		
		if (infiltration > ZERO) {
			/*--------------------------------------------------------------*/
			/*	Update patch level soil moisture with final infiltration.	*/
			/*--------------------------------------------------------------*/
			update_soil_moisture(
				command_line[0].verbose_flag,
				infiltration,
				net_inflow,
				patch,
				command_line,
				current_date );
		}
		
		patch[0].recharge = infiltration;
	} // end if hourly rain flag
    if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity){
        printf("patch_daily_F(3): (%d,%d,%d) %lf %lf %lf %lf\n",
               patch[0].ID, patch[0].soil_defaults[0][0].ID, patch[0].sat_def_pct_index,
               patch[0].sat_deficit, patch[0].sat_deficit_z,
               patch[0].rootzone.field_capacity, patch[0].field_capacity);
    }//debug
    
	/*--------------------------------------------------------------*/
	/*	Calculate patch level transpiration			*/
	/*--------------------------------------------------------------*/
	patch[0].transpiration_sat_zone = 0.0;
	patch[0].transpiration_unsat_zone = 0.0;
	patch[0].evaporation = patch[0].snowpack.sublimation;
	patch[0].PET = 0.0;
	patch[0].rain_stored = patch[0].litter.rain_stored;
	patch[0].snow_stored = 0.0;
	patch[0].ndf.plant_potential_ndemand = 0.0;
	patch[0].net_plant_psn = 0.0;
	patch[0].totalc = 0.0;
	patch[0].totaln = 0.0;
	patch[0].lai = 0.0;
	unsat_zone_patch_demand = patch[0].exfiltration_unsat_zone;
	sat_zone_patch_demand = patch[0].exfiltration_sat_zone;
	for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
		/*--------------------------------------------------------------*/
		/*	Cycle through the canopy strata in this layer		*/
		/*								*/
		/*	We assume that the stratum already has computed a	*/
		/*	unsat_zone and sat_zone transpiration demand based on	*/
		/*	its rooting depth profile.				*/
		/*--------------------------------------------------------------*/
		for ( stratum=0 ; stratum<patch[0].layers[layer].count; stratum++ ){
			strata =
				patch[0].canopy_strata[(patch[0].layers[layer].strata[stratum])];
			/*--------------------------------------------------------------*/
			/*	Add up nitrogen demand					*/
			/*--------------------------------------------------------------*/
			//patch[0].ndf.plant_potential_ndemand += strata[0].cover_fraction * (strata[0].ndf.potential_N_uptake);
			/*--------------------------------------------------------------*/
			/*	Add up evaporation.					*/
			/*--------------------------------------------------------------*/
			patch[0].evaporation +=	strata[0].cover_fraction
				* (strata[0].evaporation +	strata[0].sublimation) ;
            if(strata[0].evaporation<0 || strata[0].sublimation<0){
                printf("patch dailyF stratum [%d,%d] %e(%e,%e,%e,%e),%e(%e),\n",
                       patch[0].ID, strata[0].ID,
                       strata[0].evaporation,strata[0].potential_evaporation, strata[0].PET, strata[0].transpiration_sat_zone, strata[0].transpiration_unsat_zone,
                       strata[0].sublimation,strata[0].snow_stored);
            }//debug
			/*--------------------------------------------------------------*/
			/*      Add up canopy snow and rain stored for water balance.   */
			/*--------------------------------------------------------------*/
			patch[0].rain_stored += strata->cover_fraction * strata->rain_stored ;
			patch[0].snow_stored += strata->cover_fraction * strata->snow_stored ;
			/*--------------------------------------------------------------*/
			/*	Add uptranspiration demand 				*/
			/*--------------------------------------------------------------*/
			unsat_zone_patch_demand += strata->cover_fraction
				* strata->transpiration_unsat_zone;
			sat_zone_patch_demand += strata->cover_fraction
				* strata->transpiration_sat_zone;
			if ( command_line[0].verbose_flag > 1 ) {
				printf("\n%ld %ld %ld  -334.1 ",
					current_date.year, current_date.month, current_date.day);
				printf("\n %d %f %f %f %f %f %f %f", strata->ID,
					strata->cover_fraction,	strata->evaporation,
					strata->sublimation,strata->snow_stored,strata->rain_stored,
					strata->transpiration_unsat_zone,
					strata->transpiration_sat_zone);
			}
		}
	}

			/*--------------------------------------------------------------*/
		/* add canopy evaporation and snow sublimation to PET						*/
			/*--------------------------------------------------------------*/
		patch[0].PET = patch[0].evaporation+patch[0].evaporation_surf;

	if ( command_line[0].verbose_flag > 1 ) {
		printf("\n%ld %ld %ld  -335.1 ",
			current_date.year, current_date.month, current_date.day);
		printf("\n %d %f %f %f %f %f %f",
			patch[0].ID, patch[0].evaporation, patch[0].evaporation_surf,
			patch[0].rain_stored, patch[0].snow_stored, unsat_zone_patch_demand,
			sat_zone_patch_demand);
	}
	
	if ( command_line[0].verbose_flag == -5 ){
		printf("\n***ET DEMANDS START: rzdepth=%lf rzstor=%lf rzS=%lf rzFC=%lf rzpotsat=%lf unsatstor=%lf FC=%lf WP=%lf\n***                  S=%lf satdefz_preday=%lf satdefz=%lf satdef=%lf exfil_unsat=%lf exfil_sat=%lf unsatdemand=%lf satdemand=%lf",
			   patch[0].rootzone.depth,
			   patch[0].rz_storage,
			   0.0,//patch[0].rootzone.SatPct,
			   patch[0].rootzone.field_capacity,
			   patch[0].rootzone.potential_sat,
			   patch[0].unsat_storage,
			   patch[0].field_capacity,
			   patch[0].wilting_point,
			   0.0,//patch[0].aboveWT_SatPct,
			   patch[0].preday_sat_deficit_z,
			   patch[0].sat_deficit_z,
			   patch[0].sat_deficit,
			   patch[0].exfiltration_unsat_zone,
			   patch[0].exfiltration_sat_zone,
			   unsat_zone_patch_demand,
			   sat_zone_patch_demand);
	}
	
	
	/*--------------------------------------------------------------*/
	/* 	Fulfill transpiration demands				*/
	/*--------------------------------------------------------------*/
	/*--------------------------------------------------------------*/
	/*	Remember the transpiration demands.			*/
	/*--------------------------------------------------------------*/
	unsat_zone_patch_demand_initial = unsat_zone_patch_demand;
	sat_zone_patch_demand_initial = sat_zone_patch_demand;
    //debug
    if(unsat_zone_patch_demand<0 || sat_zone_patch_demand<0){
        // first
        // unsat_zone_patch_demand = patch[0].exfiltration_unsat_zone;
        // sat_zone_patch_demand = patch[0].exfiltration_sat_zone;
        // then, they += veg water demand
        printf("patch dailyF negative water demand (%d) %f %f, %f %f\n",
               patch[0].ID, unsat_zone_patch_demand, sat_zone_patch_demand,
               patch[0].exfiltration_unsat_zone, patch[0].exfiltration_sat_zone);
    }//if
	/*--------------------------------------------------------------*/
	/*	Figure out how much of the demand for unsat and sat zone*/
	/*	is met by supply from the zone.				*/
	/*	Update the storages while we are at it.			*/
	/*	We do not allow excess unsat storage demand beyond 	*/
	/*	current storage and what is available through cap rise. */
	/*	  Sat zone demand is not strictly limited here,		*/
	/*	although it is limited previously in 			*/
	/*	previously by rooting depth/ sat_deficit_z) 		*/
	/*								*/
	/*	this follows because					*/
	/*	ii) if the sat_zone demand is high it is because the	*/
	/*		water table is high - which means there will 	*/
	/*		is no problem meeting the sat_zone demand	*/
	/*	iii) if the unsat_zone demand is high the end of day	*/
	/*		cap rise will fill up the the unsat zone anyways*/
	/*		(note however that demands over field capacity	*/
	/*		of the unsat zone will not even be satisfied 	*/
	/*		by cap rise but thats too bad).			*/
	/*								*/
	/*	Note that the demands were decided based on 		*/
	/*	water table depths BEFORE infiltration and may bias	*/
	/*	demands on unsat zone.  However given that infiltration	*/
	/*	is big during a rain and demand small during rain it	*/
	/*	 is likely not a big deal.				*/
	/*--------------------------------------------------------------*/

	/*-------------------------------------------------------------------------*/
	/*	Compute current actual depth to water table				*/
	/*-------------------------------------------------------------------------*/
//	patch[0].sat_deficit_z = compute_z_final(
//		command_line[0].verbose_flag,
//		patch[0].soil_defaults[0][0].porosity_0,
//		patch[0].soil_defaults[0][0].porosity_decay,
//		patch[0].soil_defaults[0][0].soil_depth,
//		0.0,
//		-1.0 * patch[0].sat_deficit);
    
    // need to be careful here: patch[0].sat_deficit could be negative.
    if(patch[0].sat_deficit >= 0){
        patch[0].sat_deficit = min(patch[0].sat_deficit,patch[0].soil_defaults[0][0].soil_water_cap);
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
	//temp = patch[0].sat_deficit_z;
    temp = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index]; // total fc before ET water consumption to SAT
	
    //patch[0].rootzone.potential_sat //0-rtz
    //double sat_def = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].sat_def[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].sat_def[patch[0].sat_def_pct_index];
    available_sat_water = max(0.0, min(sat_zone_patch_demand, patch[0].rootzone.potential_sat-patch[0].sat_deficit));
    available_sat_water = min(available_sat_water, patch[0].available_soil_water);
//	available_sat_water = min(
//                              sat_zone_patch_demand,
//                              compute_delta_water(
//                                    0,
//                                    patch[0].soil_defaults[0][0].porosity_0,
//                                    patch[0].soil_defaults[0][0].porosity_decay,
//                                    patch[0].soil_defaults[0][0].soil_depth,
//                                    patch[0].rootzone.depth,
//                                    patch[0].sat_deficit_z)
//                              );

    
	patch[0].sat_deficit += available_sat_water;
	sat_zone_patch_demand -= available_sat_water;        

    // need to be careful here: patch[0].sat_deficit could be negative.
    if(patch[0].sat_deficit >= 0){
        patch[0].sat_deficit = min(patch[0].sat_deficit,patch[0].soil_defaults[0][0].soil_water_cap);
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
    
    if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity){
        printf("patch_daily_F(4): (%d,%d,%d) %lf %lf %lf %lf\n",
               patch[0].ID, patch[0].soil_defaults[0][0].ID, patch[0].sat_def_pct_index,
               patch[0].sat_deficit, patch[0].sat_deficit_z,
               patch[0].rootzone.field_capacity, patch[0].field_capacity);
    }//debug
//	patch[0].sat_deficit_z = compute_z_final(
//		command_line[0].verbose_flag,
//		patch[0].soil_defaults[0][0].porosity_0,
//		patch[0].soil_defaults[0][0].porosity_decay,
//		patch[0].soil_defaults[0][0].soil_depth,
//		0.0,
//		-1.0 * patch[0].sat_deficit);


    /*--------------------------------------------------------------*/
    /* 	leave behind field capacity			*/
    /*	if sat deficit has been lowered			*/
    /*	this should be an interactive process, we will use 	*/
    /*	0th order approximation					*/
    /* 	we do not do this once sat def is below 0.9 soil depth	*/
    /*     we use 0.9 to prevent numerical instability		*/
    /*--------------------------------------------------------------*/
    if (available_sat_water > ZERO) {
        // if watertable today is not in rootzone, sat_zone_patch_demand=0 & available_sat_water=0
        // if available_sat_water>0 then watertable today must be in the rootzone
//        add_field_capacity = compute_layer_field_capacity(
//            command_line[0].verbose_flag,
//            patch[0].soil_defaults[0][0].theta_psi_curve,
//            patch[0].soil_defaults[0][0].psi_air_entry,
//            patch[0].soil_defaults[0][0].pore_size_index,
//            patch[0].soil_defaults[0][0].p3,
//            patch[0].soil_defaults[0][0].p4,
//            patch[0].soil_defaults[0][0].porosity_0,
//            patch[0].soil_defaults[0][0].porosity_decay,
//            patch[0].sat_deficit_z,
//            patch[0].sat_deficit_z,// after ET consumption (bottom)
//            temp); // before ET consumption (top) // this layer calculation is not correct for what we need.
//        add_field_capacity = max(add_field_capacity, 0.0);
        
        totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index]; // total fc after ET water consumption to SAT
        add_field_capacity = min( patch[0].available_soil_water, max(0.0, totalfc - temp));// only acount for leave-behind water *****
        
        patch[0].sat_deficit += add_field_capacity;
        
        if ((patch[0].sat_deficit_z > patch[0].rootzone.depth) && (patch[0].preday_sat_deficit_z > patch[0].rootzone.depth))
            patch[0].unsat_storage += add_field_capacity;
            // preday_sat_deficit_z is calculated by subsurface routing from previous day
        else if ((patch[0].sat_deficit_z <= patch[0].rootzone.depth) && (patch[0].preday_sat_deficit_z <= patch[0].rootzone.depth))
            patch[0].rz_storage += add_field_capacity;
        else {
            patch[0].rz_storage += add_field_capacity * (patch[0].rootzone.depth -patch[0].preday_sat_deficit_z)
                / (patch[0].sat_deficit_z -patch[0].preday_sat_deficit_z);
            patch[0].unsat_storage += add_field_capacity * (patch[0].sat_deficit_z - patch[0].rootzone.depth)
                / (patch[0].sat_deficit_z -patch[0].preday_sat_deficit_z);
        }//
        
        if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity){
            printf("patch_daily_F(5): (%d,%d,%d) %lf %lf %lf %lf\n",
                   patch[0].ID, patch[0].soil_defaults[0][0].ID, patch[0].sat_def_pct_index,
                   patch[0].sat_deficit, patch[0].sat_deficit_z,
                   patch[0].rootzone.field_capacity, patch[0].field_capacity);
        }//debug
        
    }//if
    
	/*--------------------------------------------------------------*/
	/*	See how much of  unsat zone demand can be 		*/
	/*	met and still field capacity. 				*/
	/*--------------------------------------------------------------*/
    if (patch[0].rootzone.depth > ZERO ) { /* VEG CASE */
		water_above_field_cap = max((patch[0].rz_storage - patch[0].rootzone.field_capacity), 0);
		water_above_field_cap = min(unsat_zone_patch_demand, water_above_field_cap);
		patch[0].rz_storage -= water_above_field_cap;
		unsat_zone_patch_demand -= water_above_field_cap; 
	}
	else { /* NO VEG CASE (NEED TO CHECK THIS) */
		water_above_field_cap = max((patch[0].unsat_storage - patch[0].field_capacity), 0);  
		water_above_field_cap = min(unsat_zone_patch_demand, water_above_field_cap);
		patch[0].unsat_storage -= water_above_field_cap;
		unsat_zone_patch_demand -= water_above_field_cap; 
	}

	/*--------------------------------------------------------------*/
	/*	compute new field capacity				*/
	/*--------------------------------------------------------------*/
    
    // need to be careful here: patch[0].sat_deficit could be negative.
    if(patch[0].sat_deficit >= 0){
        patch[0].sat_deficit = min(patch[0].sat_deficit,patch[0].soil_defaults[0][0].soil_water_cap);
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
    // update fc
    totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index];
    //totalfc *= (1.0-patch[0].basementFrac); // <---- second thought on this, Oct 8, 2019; basement is 3m at most
    
    
	if (patch[0].rootzone.depth > ZERO ) { /* VEG CASE */
        
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
        
//		if (patch[0].sat_deficit_z < patch[0].rootzone.depth)  {
//			patch[0].rootzone.field_capacity = compute_layer_field_capacity(
//				command_line[0].verbose_flag,
//				patch[0].soil_defaults[0][0].theta_psi_curve,
//				patch[0].soil_defaults[0][0].psi_air_entry,
//				patch[0].soil_defaults[0][0].pore_size_index,
//				patch[0].soil_defaults[0][0].p3,
//				patch[0].soil_defaults[0][0].p4,
//				patch[0].soil_defaults[0][0].porosity_0,
//				patch[0].soil_defaults[0][0].porosity_decay,
//				patch[0].sat_deficit_z,
//				patch[0].rootzone.depth, 0.0)*(1.0-patch[0].basementFrac);
//
//			patch[0].field_capacity = 0.0;
//			}
//		else  {
//			patch[0].rootzone.field_capacity = compute_layer_field_capacity(
//				command_line[0].verbose_flag,
//				patch[0].soil_defaults[0][0].theta_psi_curve,
//				patch[0].soil_defaults[0][0].psi_air_entry,
//				patch[0].soil_defaults[0][0].pore_size_index,
//				patch[0].soil_defaults[0][0].p3,
//				patch[0].soil_defaults[0][0].p4,
//				patch[0].soil_defaults[0][0].porosity_0,
//				patch[0].soil_defaults[0][0].porosity_decay,
//				patch[0].sat_deficit_z,
//				patch[0].rootzone.depth, 0.0);
//
//			patch[0].field_capacity = compute_layer_field_capacity(
//				command_line[0].verbose_flag,
//				patch[0].soil_defaults[0][0].theta_psi_curve,
//				patch[0].soil_defaults[0][0].psi_air_entry,
//				patch[0].soil_defaults[0][0].pore_size_index,
//				patch[0].soil_defaults[0][0].p3,
//				patch[0].soil_defaults[0][0].p4,
//				patch[0].soil_defaults[0][0].porosity_0,
//				patch[0].soil_defaults[0][0].porosity_decay,
//				patch[0].sat_deficit_z,
//				patch[0].sat_deficit_z, 0.0) - patch[0].rootzone.field_capacity;
//
//            patch[0].rootzone.field_capacity *= (1.0-patch[0].basementFrac);
//			}
        // patch[0].rootzone.depth must be no zero
		if (patch[0].sat_deficit_z > patch[0].rootzone.depth) 
				water_below_field_cap = patch[0].rootzone.field_capacity - patch[0].rz_storage;
		else
			water_below_field_cap = patch[0].rootzone.field_capacity - patch[0].rz_storage - (patch[0].rootzone.potential_sat-patch[0].sat_deficit);
    } else  {
        /* NO VEG CASE (NEED TO CHECK THIS) */
        patch[0].field_capacity = totalfc;
        //patch[0].rootzone.SatPct = 0.0;
        //patch[0].aboveWT_SatPct = (patch[0].rz_storage+patch[0].unsat_storage)/patch[0].sat_deficit;
//        patch[0].field_capacity = compute_layer_field_capacity(
//           command_line[0].verbose_flag,
//           patch[0].soil_defaults[0][0].theta_psi_curve,
//           patch[0].soil_defaults[0][0].psi_air_entry,
//           patch[0].soil_defaults[0][0].pore_size_index,
//           patch[0].soil_defaults[0][0].p3,
//           patch[0].soil_defaults[0][0].p4,
//           patch[0].soil_defaults[0][0].porosity_0,
//           patch[0].soil_defaults[0][0].porosity_decay,
//           patch[0].sat_deficit_z,
//           patch[0].sat_deficit_z, 0.0);
    
        water_below_field_cap = patch[0].field_capacity - patch[0].unsat_storage;
    } /* END NO VEG CASE */
	
    if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity){
        printf("patch_daily_F(6): (%d,%d,%d) %lf %lf %lf %lf\n",
               patch[0].ID, patch[0].soil_defaults[0][0].ID, patch[0].sat_def_pct_index,
               patch[0].sat_deficit, patch[0].sat_deficit_z,
               patch[0].rootzone.field_capacity, patch[0].field_capacity);
    }//debug
    
	/*--------------------------------------------------------------*/
	/*	fill the leftover demand with cap rise.			*/
	/*--------------------------------------------------------------*/
    // not sure why compute this block again at here; should not have changed.
    if(patch[0].sat_deficit_z >= patch[0].rootzone.depth ){
        patch[0].unsat_deficit = patch[0].sat_deficit - patch[0].unsat_storage - patch[0].rootzone.potential_sat;
        if(patch[0].rootzone.potential_sat > patch[0].rz_storage){
            patch[0].rz_deficit = patch[0].rootzone.potential_sat - patch[0].rz_storage;
        }else{
            patch[0].rz_deficit=0;
        }
    }else{
        patch[0].unsat_deficit = 0.0;
        if(patch[0].sat_deficit > patch[0].rz_storage){
            patch[0].rz_deficit = patch[0].sat_deficit - patch[0].rz_storage;
        }else{
            patch[0].rz_deficit=0;
        }
    }//if else
    //patch[0].UNSAT_fieldcap = patch[0].field_capacity + patch[0].rootzone.field_capacity;
    
   
    // new cap rise code
    // patch[0].potential_cap_rise is a daily flux without consideration of existing soil moisture;
    // "unsat_zone_patch_demand" at this point has been -= water_above_field_cap; it's time to determine if cap_rise can satisfy it.
    // cap_rise = (patch[0].potential_cap_rise)@rtz vs unsat_zone_patch_demand

    // because cap_rise is very sensitive to the sat_def. the potential_cap_rise was calculated @ patch dailyI; here we calculate potential_cap_rise once more and take the average of the two (testing).
    // rtz vs rtz fc
    double plant_drive_cap_rise = 1.0;
    if(patch[0].rz_storage < patch[0].rootzone.field_capacity){
        plant_drive_cap_rise = min(3.0, patch[0].rootzone.field_capacity/(patch[0].rz_storage+0.001));
    }
    
    cap_rise = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].pot_caprise_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) *  patch[0].soil_defaults[0][0].pot_caprise_0z[patch[0].sat_def_pct_index];
    cap_rise = min(cap_rise, patch[0].available_soil_water);
    cap_rise = min(cap_rise, max(0.0,totalfc - patch[0].unsat_storage - patch[0].rz_storage));
    patch[0].potential_cap_rise += cap_rise;
    patch[0].potential_cap_rise *= 0.5;
    //patch[0].potential_cap_rise *= plant_drive_cap_rise;
    
    cap_rise = patch[0].potential_cap_rise * patch[0].zeroRootCoef * (patch[0].rootzone_scale_ref*patch[0].rootzone_end_refcap[patch[0].sat_def_pct_index] + (1.0-patch[0].rootzone_scale_ref)*patch[0].rootzone_start_refcap[patch[0].sat_def_pct_index]);// cap_rise that can reach rtz
    //cap_rise *= plant_drive_cap_rise;
    cap_rise_unsat = max(0.0, min(patch[0].potential_cap_rise-cap_rise, patch[0].field_capacity-patch[0].unsat_storage));
    patch[0].unsat_storage += cap_rise_unsat;
    // run into problem here. FC_0z is great when satz is deep; and most of that FC_0z is un unsat! at the same time, cap_rise becomes less as satz is deep
    // so how to solve this? it make all satwater into unsat as satz drops
    
    if( cap_rise >= unsat_zone_patch_demand ){
        cap_rise = max(0.0,min(cap_rise - unsat_zone_patch_demand, patch[0].rootzone.field_capacity-patch[0].rz_storage));
        patch[0].rz_storage += cap_rise;
        cap_rise += unsat_zone_patch_demand + cap_rise_unsat; // becomes total cap_rise
        unsat_zone_patch_demand=0.0;
    }else{
        // cap_rise <= unsat_zone_patch_demand;
        unsat_zone_patch_demand -= cap_rise;
        cap_rise += cap_rise_unsat;
        // done
    }// if else
    patch[0].sat_deficit += cap_rise;
    patch[0].cap_rise = cap_rise;
    if(cap_rise > 0.0 && patch[0].available_soil_water > 0.0 && command_line[0].grow_flag > 0){
        // cap_rise also carries satSolute up to unsat
        double tmp_ratio = min(1.0, cap_rise/patch[0].available_soil_water);
        //tmp_ratio /= patch[0].soil_defaults[0][0].porosity_0 * patch[0].soil_defaults[0][0].porosity_decay * (exp(-(patch[0].sat_deficit_z>0? patch[0].sat_deficit_z : 0.0)/patch[0].soil_defaults[0][0].porosity_decay) - exp(-patch[0].soil_defaults[0][0].soil_depth/patch[0].soil_defaults[0][0].porosity_decay));// integration of sat_z and soildepth
        patch[0].soil_ns.nitrate += patch[0].sat_NO3 * tmp_ratio;
        patch[0].soil_ns.sminn += patch[0].sat_NH4 * tmp_ratio;
        patch[0].soil_cs.DOC += patch[0].sat_DOC * tmp_ratio;
        patch[0].soil_ns.DON += patch[0].sat_DON * tmp_ratio;
        patch[0].sat_NO3 *= 1.0 - tmp_ratio;
        patch[0].sat_NH4 *= 1.0 - tmp_ratio;
        patch[0].sat_DOC *= 1.0 - tmp_ratio;
        patch[0].sat_DON *= 1.0 - tmp_ratio;
    }// if
    if( unsat_zone_patch_demand > 0.0){
        /*--------------------------------------------------------------*/
        /*    Now supply the remaining demand with water left in    */
        /*    the unsat zone.  We are going below field cap now!!    */
        /*    First guess at change in sat storage to meet demand.    */
        /*--------------------------------------------------------------*/
        delta_unsat_zone_storage = min(unsat_zone_patch_demand, patch[0].rz_storage); // the ET removed water starts for RZ first, but it was previously instant reduced by CAP [related to above]!
        
        if ((patch[0].rz_storage > ZERO) && (patch[0].sat_deficit > ZERO)) {
            patch[0].wilting_point = exp(-patch[0].soil_defaults[0][0].pore_size_index * log(-100.0*patch[0].psi_max_veg/patch[0].soil_defaults[0][0].psi_air_entry));
            patch[0].wilting_point *= min(patch[0].sat_deficit, patch[0].rootzone.potential_sat);
            //patch[0].wilting_point *= (1.0-patch[0].basementFrac); // needed here? probably no; no 'SatPct' involved here. basement effect is still massed up when interact with veg uptake. however basement effect on routing is fine.
            delta_unsat_zone_storage = min(patch[0].rz_storage-patch[0].wilting_point, delta_unsat_zone_storage); // make sure RZ does not go below wilting point!! [hard correction]
            delta_unsat_zone_storage = max(delta_unsat_zone_storage, 0.0);
        }
        else {
            patch[0].wilting_point = 0;
        }
        
        patch[0].rz_storage -= delta_unsat_zone_storage;
        unsat_zone_patch_demand -= delta_unsat_zone_storage;
    }// if
    
    if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity){
        printf("patch_daily_F(7): (%d,%d,%d) %lf %lf %lf %lf\n",
               patch[0].ID, patch[0].soil_defaults[0][0].ID, patch[0].sat_def_pct_index,
               patch[0].sat_deficit, patch[0].sat_deficit_z,
               patch[0].rootzone.field_capacity, patch[0].field_capacity);
    }//debug
	/*--------------------------------------------------------------*/
    // problem could be from cap_rise and gw_dainage
    if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
    patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
    patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
    patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
    patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
    patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
    patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
    patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
    patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0) printf("patch daily F4N %d-%d-%d [%d,%d,%d,%d] [%e %e %e] [%e %e %e] [%e %e %e]\n",
       current_date.year, current_date.month, current_date.day,
       patch[0].ID, patch[0].drainage_type, actionPIPEDRAIN, actionSEWER,
       patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
       patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
       patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);
//	/*--------------------------------------------------------------*/
//	/* 	Resolve plant uptake and soil microbial N demands	*/
//	/*--------------------------------------------------------------*/
//	if (command_line[0].grow_flag > 0)  {
//        resolve_sminn_competition(
//                                  &(patch[0].soil_ns),
//                                  patch,
//                                  command_line,
//                                  &(patch[0].ndf));
//	}//growth_flag
	/*--------------------------------------------------------------*/
	/*	Reduce the stratum actual transpiration and compute 	*/
	/*	final N-uptake and daily allocation as a function 	*/
	/*	of available C and N (if growth flag is on)		*/
	/* 	we reduce carbon flux by current water use efficiency 	*/
	/*--------------------------------------------------------------*/

	patch[0].trans_reduc_perc = 1.0;
	transpiration_reduction_percent = 1.0; // this ratio is used to "correct" plant growth!!
    // transpiration_reduction_percent = 1 means no reduction
	
	if (patch[0].rootzone.depth > ZERO ) {
        if ( unsat_zone_patch_demand_initial+sat_zone_patch_demand_initial > ZERO)
            transpiration_reduction_percent = max(0.0, 1.0 - (unsat_zone_patch_demand + sat_zone_patch_demand) / (unsat_zone_patch_demand_initial + sat_zone_patch_demand_initial));
        else
            transpiration_reduction_percent = 1.0;
	}
    
    //debug
    if(transpiration_reduction_percent<0){
        printf("patch_daily transpiration_reduction_percent (%d) %f = %f + %f over their inits %f + %f, %f %f %f %d, %f %f %f %f\n",
               patch[0].ID, transpiration_reduction_percent,
               unsat_zone_patch_demand, sat_zone_patch_demand,
               unsat_zone_patch_demand_initial, sat_zone_patch_demand_initial,
               patch[0].cap_rise, patch[0].potential_cap_rise, cap_rise_unsat, patch[0].sat_def_pct_index,
               patch[0].field_capacity, patch[0].unsat_storage,
               patch[0].rootzone.field_capacity, patch[0].rz_storage);
    }//if debug

	if ( command_line[0].verbose_flag == -5 ){
		printf("\n***START: exfil_unsat=%lf exfil_sat=%lf unsatdemand_ini=%lf unsatdemand=%lf satdemand_ini=%lf satdemand=%lf",
			   patch[0].exfiltration_unsat_zone,
			   patch[0].exfiltration_sat_zone,
			   unsat_zone_patch_demand_initial,
			   unsat_zone_patch_demand,
			   sat_zone_patch_demand_initial,
			   sat_zone_patch_demand);
	}
	
	
	if ( unsat_zone_patch_demand_initial > 0 ){
		patch[0].exfiltration_unsat_zone = patch[0].exfiltration_unsat_zone
			* (1 - unsat_zone_patch_demand / unsat_zone_patch_demand_initial );
        
		patch[0].transpiration_unsat_zone = patch[0].transpiration_unsat_zone
			* (1 - unsat_zone_patch_demand / unsat_zone_patch_demand_initial );
        
		if ( command_line[0].verbose_flag == -5 ){
			printf("\n***CASE1 TRIGGERED: exfil_unsat=%lf demand_ini=%lf demand=%lf",patch[0].exfiltration_unsat_zone,unsat_zone_patch_demand_initial,unsat_zone_patch_demand);
			}
    }
	if ( sat_zone_patch_demand_initial > 0 ) {
		patch[0].exfiltration_sat_zone = patch[0].exfiltration_sat_zone
			* (1 - sat_zone_patch_demand /  sat_zone_patch_demand_initial );
		patch[0].transpiration_sat_zone = patch[0].transpiration_sat_zone
			* (1 - sat_zone_patch_demand /  sat_zone_patch_demand_initial );
		if ( command_line[0].verbose_flag == -5 ){
			printf("\n***CASE2 TRIGGERED: exfil_sat=%lf demand_ini=%lf demand=%lf",patch[0].exfiltration_sat_zone,sat_zone_patch_demand_initial,sat_zone_patch_demand);
		}		
    }
	
	patch[0].trans_reduc_perc = transpiration_reduction_percent;

	/*--------------------------------------------------------------*/
	/* add soil evap to PET																					*/
	/*--------------------------------------------------------------*/

	patch[0].PET += (patch[0].exfiltration_sat_zone + patch[0].exfiltration_unsat_zone);

	/*--------------------------------------------------------------*/
	/* in order to restrict denitri/nitrific on non-veg patches type */
	/* 	tag vegtype							*/	
	/*--------------------------------------------------------------*/
  vegtype = 0;
  patch[0].target_status = 1;

	for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
		for ( stratum=0 ; stratum < patch[0].layers[layer].count; stratum++ ){
            
			strata = patch[0].canopy_strata[(patch[0].layers[layer].strata[stratum])];
   
            if(command_line[0].vegspinup_flag > 0){
                if (strata->target.met == 0) patch[0].target_status = 0;
            }
			  
            if ( strata[0].defaults[0][0].epc.veg_type != NON_VEG ){

			   	if (transpiration_reduction_percent < 1.0) {
                    // strata->cdf.psn_to_cpool is updated by water availability
                    // strata->cs.availc = cdf.psn_to_cpool - cdf.total_mr + ifelse(cpool<0, cpool, 0); and "cpool" has be modified
                    double correction = strata->cdf.psn_to_cpool * (1.0-transpiration_reduction_percent);
                    double correction2 = (strata->cs.availc>0? max(strata->cs.availc-correction,0.0)/strata->cs.availc : 0.0);
                    strata->ndf.potential_N_uptake *= correction2;
                    strata->cdf.psn_to_cpool = strata->cdf.psn_to_cpool * transpiration_reduction_percent;
                    strata->cs.availc -= correction;
                    if(strata->cs.availc<0) strata->cs.availc = 0.0;
                    
                    // does not seem to affect a thing for these variables below.
                    strata->gs_sunlit *= transpiration_reduction_percent;//gs_sunlit is calculated in canopy_stratum_daily_F @1383 in this patchDailyF
                    strata->gs_shade *= transpiration_reduction_percent;//gs_shade is calculated in canopy_stratum_daily_F @ in this patchDailyF
                    //strata->mult_conductance.LWP *= transpiration_reduction_percent;//useful?
                    //strata->ndf.potential_N_uptake *= transpiration_reduction_percent; wrong to do this
				}//if transpiration_reduction_percent<1
                
                
                if(strata[0].phen.gwseasonday>0){
                    // transpiration_reduction_percent = 1 when available water satistfy water demand; =0 when not
                    strata[0].wFactor += transpiration_reduction_percent;
                }
                    

            } else {
                //non_VEG
				if ( strata->phen.annual_allocation == 1){
					strata->cdf.leafc_store_to_leafc_transfer = strata->cs.leafc_store;
					strata->cs.leafc_transfer += strata->cdf.leafc_store_to_leafc_transfer;
					strata->cs.leafc_store -= strata->cdf.leafc_store_to_leafc_transfer;
                    
					strata->ndf.leafn_store_to_leafn_transfer = strata->ns.leafn_store;
					strata->ns.leafn_transfer += strata->ndf.leafn_store_to_leafn_transfer;
					strata->ns.leafn_store -= strata->ndf.leafn_store_to_leafn_transfer;
				}
			}//if
                 
			if ( unsat_zone_patch_demand_initial > 0.0 )
				strata->transpiration_unsat_zone = strata->transpiration_unsat_zone *(1-unsat_zone_patch_demand/unsat_zone_patch_demand_initial);
			if ( sat_zone_patch_demand_initial > 0.0 )
				strata->transpiration_sat_zone = strata->transpiration_sat_zone *(1 - sat_zone_patch_demand / sat_zone_patch_demand_initial );
			
            patch[0].ndf.plant_potential_ndemand += strata->cover_fraction * strata->ndf.potential_N_uptake; // stratum_dailyF() from far above
            
            patch[0].transpiration_unsat_zone += strata->cover_fraction * strata->transpiration_unsat_zone; // incorrect? for fire model
			patch[0].transpiration_sat_zone += strata->cover_fraction * strata->transpiration_sat_zone; // incorrect? for fire model
            
			patch[0].PET += strata->cover_fraction * strata->PET; // correct
			patch[0].totalc += strata->cover_fraction * strata->cs.totalc;
			patch[0].totaln += strata->cover_fraction * strata->ns.totaln;
			patch[0].net_plant_psn += strata->cover_fraction *	strata->cs.net_psn; // correct
			patch[0].lai += strata->cover_fraction * strata->epv.proj_lai; // correct
		}//strata loop
	}// layer loop
    
    /*--------------------------------------------------------------*/
    /*     Resolve plant uptake and soil microbial N demands    */
    /*--------------------------------------------------------------*/
    if (command_line[0].grow_flag > 0)  {
        resolve_sminn_competition(
                                  &(patch[0].soil_ns),
                                  patch,
                                  command_line,
                                  &(patch[0].ndf));
    }//growth_flag
    
    for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
         for ( stratum=0 ; stratum < patch[0].layers[layer].count; stratum++ ){
             strata = patch[0].canopy_strata[(patch[0].layers[layer].strata[stratum])];
             if ( strata[0].defaults[0][0].epc.veg_type != NON_VEG ){
                 // problem "canopy_stratum_growth()" need nlimit and adds "sminn_to_npool" to patch[0].ndf->sminn_to_npool
                 vegtype=1;
                 canopy_stratum_growth(
                     world,
                     basin,
                     hillslope,
                     zone,
                     patch,
                     strata,
                     command_line,
                     event,
                     current_date );
             }// end of if
         }//strata loop
     }// layer loop
    
    
	/*-------------------------------------------------------------------------*/
	/*	Compute current actual depth to water table				*/
	/*------------------------------------------------------------------------*/

    
    // need to be careful here: patch[0].sat_deficit could be negative.
    if(patch[0].sat_deficit >= 0){
        patch[0].sat_deficit = min(patch[0].sat_deficit,patch[0].soil_defaults[0][0].soil_water_cap);
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
    // update fc
    totalfc = patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].fc1_0z[patch[0].sat_def_pct_index];
    //totalfc *= (1.0-patch[0].basementFrac); // <---- second thought on this, Oct 8, 2019; basement is 3m at most
    
    if (patch[0].sat_deficit < ZERO) {
        //patch[0].aboveWT_SatPct = 1.0;
        //patch[0].rootzone.SatPct = 1.0;
        patch[0].rootzone.field_capacity = 0.0;
        patch[0].field_capacity = 0.0;
        rz_drainage = 0.0;
        unsat_drainage = 0.0;
    } else {
        patch[0].rootzone.field_capacity = totalfc * patch[0].zeroRootCoef * (patch[0].rootzone_scale_ref*patch[0].rootzone_end_reffc[patch[0].sat_def_pct_index] + (1.0-patch[0].rootzone_scale_ref)*patch[0].rootzone_start_reffc[patch[0].sat_def_pct_index]);
        patch[0].rootzone.field_capacity = min(patch[0].rootzone.field_capacity,patch[0].rootzone.potential_sat);
        patch[0].field_capacity = max(0.0,min(patch[0].sat_deficit-patch[0].rootzone.potential_sat, totalfc - patch[0].rootzone.field_capacity));
        
        // vertical drainage from top soil layer to SAT
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
                patch[0].sat_deficit - patch[0].rootzone.potential_sat, //patch[0].aboveWT_SatPct,
                patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM)* patch[0].soil_defaults[0][0].vksat_0zm[patch[0].sat_def_pct_index],
                patch[0].sat_def_pct_indexM * patch[0].soil_defaults[0][0].vksat_z[patch[0].sat_def_pct_index+1] + (1.0-patch[0].sat_def_pct_indexM) * patch[0].soil_defaults[0][0].vksat_z[patch[0].sat_def_pct_index],
                patch[0].unsat_storage,
                patch[0].field_capacity,
                patch[0].sat_deficit);
        patch[0].unsat_storage -=  unsat_drainage;
        patch[0].sat_deficit -=  unsat_drainage;
   
        // need to move the solutes as well
        if(unsat_drainage > 0.0 && command_line[0].grow_flag > 0){
            // cap_rise also carries satSolute up to unsat
            if(patch[0].rz_storage>0){
                double tmp_ratio = max(0.0,min(1.0,unsat_drainage/patch[0].rz_storage)); // rz_drainage/patch[0].rz_storage * unsat_drainage/rz_drainage
                patch[0].sat_NO3 += patch[0].soil_ns.nitrate * patch[0].soil_defaults[0][0].rtz2NO3prop[patch[0].soil_defaults[0][0].active_zone_index]*tmp_ratio;
                patch[0].sat_NH4 += patch[0].soil_ns.sminn * patch[0].soil_defaults[0][0].rtz2NH4prop[patch[0].soil_defaults[0][0].active_zone_index]*tmp_ratio;
                patch[0].sat_DOC += patch[0].soil_cs.DOC * patch[0].soil_defaults[0][0].rtz2DOMprop[patch[0].soil_defaults[0][0].active_zone_index]*tmp_ratio;
                patch[0].sat_DON += patch[0].soil_ns.DON * patch[0].soil_defaults[0][0].rtz2DOMprop[patch[0].soil_defaults[0][0].active_zone_index]*tmp_ratio;
                
                tmp_ratio = max(0.0, 1.0 - tmp_ratio);
                patch[0].soil_ns.nitrate *= tmp_ratio;
                patch[0].soil_ns.sminn *= tmp_ratio;
                patch[0].soil_cs.DOC *= tmp_ratio;
                patch[0].soil_ns.DON *= tmp_ratio;
            }else if(patch[0].unsat_storage>0){
                double tmp_ratio = max(0.0,min(1.0,unsat_drainage/patch[0].unsat_storage));
                patch[0].sat_NO3 += patch[0].soil_ns.nitrate * patch[0].soil_defaults[0][0].rtz2NO3prop[patch[0].soil_defaults[0][0].active_zone_index]*tmp_ratio;
                patch[0].sat_NH4 += patch[0].soil_ns.sminn * patch[0].soil_defaults[0][0].rtz2NH4prop[patch[0].soil_defaults[0][0].active_zone_index]*tmp_ratio;
                patch[0].sat_DOC += patch[0].soil_cs.DOC * patch[0].soil_defaults[0][0].rtz2DOMprop[patch[0].soil_defaults[0][0].active_zone_index]*tmp_ratio;
                patch[0].sat_DON += patch[0].soil_ns.DON * patch[0].soil_defaults[0][0].rtz2DOMprop[patch[0].soil_defaults[0][0].active_zone_index]*tmp_ratio;
                
                tmp_ratio = max(0.0, 1.0 - tmp_ratio);
                patch[0].soil_ns.nitrate *= tmp_ratio;
                patch[0].soil_ns.sminn *= tmp_ratio;
                patch[0].soil_cs.DOC *= tmp_ratio;
                patch[0].soil_ns.DON *= tmp_ratio;
            }
        }// if
        
    }//if else
    patch[0].unsat_drainage += unsat_drainage;
    patch[0].rz_drainage += rz_drainage;
    patch[0].hourly_unsat_drainage += unsat_drainage;
    patch[0].hourly_rz_drainage += rz_drainage;
    
    if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity){
        printf("patch_daily_F(8): (%d,%d,%d) %lf %lf %lf %lf\n",
               patch[0].ID, patch[0].soil_defaults[0][0].ID, patch[0].sat_def_pct_index,
               patch[0].sat_deficit, patch[0].sat_deficit_z,
               patch[0].rootzone.field_capacity, patch[0].field_capacity);
    }//debug

    
    // need to be careful here: patch[0].sat_deficit could be negative.
    if(patch[0].sat_deficit >= 0){
        patch[0].sat_deficit = min(patch[0].sat_deficit,patch[0].soil_defaults[0][0].soil_water_cap);
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
    //totalfc *= (1.0-patch[0].basementFrac); // <---- second thought on this, Oct 8, 2019; basement is 3m at most
    
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
    
    if(patch[0].sat_deficit!=patch[0].sat_deficit || patch[0].sat_deficit_z!=patch[0].sat_deficit_z || patch[0].rootzone.field_capacity!=patch[0].rootzone.field_capacity || patch[0].field_capacity!=patch[0].field_capacity || patch[0].sat_deficit>patch[0].soil_defaults[0][0].soil_water_cap+ZERO){
        printf("patch_daily_F(9): (%d,%d,%d) %lf %lf %lf %lf, %lf\n",
               patch[0].ID, patch[0].soil_defaults[0][0].ID, patch[0].sat_def_pct_index,
               patch[0].sat_deficit, patch[0].sat_deficit_z,
               patch[0].rootzone.field_capacity, patch[0].field_capacity, patch[0].soil_defaults[0][0].soil_water_cap);
    }//debug
    
    /*-----------------------------------------------------*/
    /*  re-Compute potential saturation for rootzone layer   */
    /*-----------------------------------------------------*/
//    if (patch[0].rootzone.depth > ZERO)
//        patch[0].rootzone.potential_sat = compute_delta_water(
//            command_line[0].verbose_flag,
//            patch[0].soil_defaults[0][0].porosity_0,
//            patch[0].soil_defaults[0][0].porosity_decay,
//            patch[0].soil_defaults[0][0].soil_depth,
//            patch[0].rootzone.depth, 0.0);
    
     // doing this in daily_I()
//    if (patch[0].rootzone.depth > ZERO)  {
//        // how is lulc frac affecting this part? Aug 8, 2019
//        // daily updated "patch[0].rootzone.depth" is done by this daily_I @ LINE 472 (below)
//        // stratum[0].rootzone.depth is first updated via "update_phenology", then aggregated to here
//

    
    if(patch[0].rootzone.potential_sat>ZERO){
        if (patch[0].sat_deficit > patch[0].rootzone.potential_sat) theta = min(patch[0].rz_storage/patch[0].rootzone.potential_sat, 1.0);//(1.0-patch[0].basementFrac)
        else theta = min((patch[0].rz_storage + patch[0].rootzone.potential_sat - patch[0].sat_deficit)/patch[0].rootzone.potential_sat,1.0);//(1.0-patch[0].basementFrac)
    }else{ theta = 0.0; }
    patch[0].theta_std = patch[0].soil_defaults[0][0].active_zone_sat_0z*theta;
    patch[0].theta_std *= patch[0].soil_defaults[0][0].theta_mean_std_p2 * patch[0].theta_std;
    patch[0].theta_std += patch[0].soil_defaults[0][0].theta_mean_std_p1 * (patch[0].soil_defaults[0][0].active_zone_sat_0z*theta);
    patch[0].theta_std = max(0.0, patch[0].theta_std);
    
    
    
    /*-----------------------------------------------------*/
    /*  snow related   */
    /*-----------------------------------------------------*/
    patch[0].delta_snowpack = patch[0].snowpack.water_depth
        + patch[0].snowpack.water_equivalent_depth - patch[0].preday_snowpack;
    patch[0].delta_rain_stored = patch[0].rain_stored
        - patch[0].preday_rain_stored;
    patch[0].delta_snow_stored = patch[0].snow_stored
        - patch[0].preday_snow_stored;
    
	/*-------------------------------------------------------------------------*/
	/*	finalized soil and litter decomposition					*/
	/* 	and any septic losses							*/
	/*------------------------------------------------------------------------*/
    if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
    patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
    patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
    patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
    patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
    patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
    patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
    patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
    patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0) printf("patch daily F5N %d-%d-%d [%d,%d,%d,%d] [%e %e %e] [%e %e %e] [%e %e %e]\n",
       current_date.year, current_date.month, current_date.day,
       patch[0].ID, patch[0].drainage_type, actionPIPEDRAIN, actionSEWER,
       patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
       patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
       patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);
	
    if ((command_line[0].grow_flag > 0) && (vegtype == 1) ) {
		
        if ( update_decomp(
            current_date,
            &(patch[0].soil_cs),
            &(patch[0].soil_ns),
            &(patch[0].litter_cs),
            &(patch[0].litter_ns),
            &(patch[0].cdf),
            &(patch[0].ndf),
            patch,
            command_line
            ) != 0){
            fprintf(stderr,"fATAL ERROR: in update_decomp() ... Exiting\n");
            exit(EXIT_FAILURE);
        }

        if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
        patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
        patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
        patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
        patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
        patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
        patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
        patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
        patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0) printf("patch daily F6N %d-%d-%d after update decomp [%d,%d,%d,%d] [%e %e %e] [%e %e %e] [%e %e %e]\n",
           current_date.year, current_date.month, current_date.day,
           patch[0].ID, patch[0].drainage_type, actionPIPEDRAIN, actionSEWER,
           patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
           patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
           patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);
        
        if (patch[0].soil_defaults[0][0].DON_production_rate > ZERO) {
            if ( update_dissolved_organic_losses(
                current_date,
                patch[0].soil_defaults[0][0].DON_production_rate,
                &(patch[0].soil_cs),
                &(patch[0].soil_ns),
                &(patch[0].litter_cs),
                &(patch[0].litter_ns),
                &(patch[0].cdf),
                &(patch[0].ndf),
                patch,
                command_line[0].soilCNadaptation_flag) != 0){
                fprintf(stderr,"fATAL ERROR: in update_dissolved_organic_losses() ... Exiting\n");
                exit(EXIT_FAILURE);
            }
        patch[0].surface_DOC += (patch[0].cdf.do_litr1c_loss +
                patch[0].cdf.do_litr2c_loss + patch[0].cdf.do_litr3c_loss + patch[0].cdf.do_litr4c_loss);
        patch[0].surface_DON += (patch[0].ndf.do_litr1n_loss + patch[0].ndf.do_litr2n_loss + patch[0].ndf.do_litr3n_loss +
                 patch[0].ndf.do_litr4n_loss);
        }

        if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
        patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
        patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
        patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
        patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
        patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
        patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
        patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
        patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0) printf("patch daily F7N after DOloss %d-%d-%d [%d,%d,%d,%d,%e] [%e %e %e] [%e %e %e] [%e %e %e]\n",
           current_date.year, current_date.month, current_date.day,
           patch[0].ID, patch[0].drainage_type, actionPIPEDRAIN, actionSEWER, patch[0].soil_ns.fract_potential_immob,
           patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
           patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
           patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);
        
		if ( update_nitrif(
			&(patch[0].soil_cs),
			&(patch[0].soil_ns),
			&(patch[0].cdf),
			&(patch[0].ndf),
            patch,
            command_line,
            patch[0].theta_std) != 0){
			fprintf(stderr,"fATAL ERROR: in update_nitrific() ... Exiting\n");
			exit(EXIT_FAILURE);
		}
        if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
        patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
        patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
        patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
        patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
        patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
        patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
        patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
        patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0) printf("patch daily F8N after nitrif %d-%d-%d [%d,%d,%d,%d] [%e %e %e] [%e %e %e] [%e %e %e]\n",
           current_date.year, current_date.month, current_date.day,
           patch[0].ID, patch[0].drainage_type, actionPIPEDRAIN, actionSEWER,
           patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
           patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
           patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);
        
		if ( update_denitrif(
			&(patch[0].soil_cs),
			&(patch[0].soil_ns),
			&(patch[0].cdf),
			&(patch[0].ndf),
			patch,
            command_line,
            patch[0].theta_std) != 0){
			fprintf(stderr,"fATAL ERROR: in update_denitrif() ... Exiting\n");
			exit(EXIT_FAILURE);
		}
        if(patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 ||
        patch[0].soil_ns.sminn!=patch[0].soil_ns.sminn || patch[0].soil_ns.sminn<0 ||
        patch[0].soil_cs.DOC!=patch[0].soil_cs.DOC || patch[0].soil_cs.DOC<0 ||
        patch[0].sat_NO3!=patch[0].sat_NO3 || patch[0].sat_NO3<0 ||
        patch[0].sat_NH4!=patch[0].sat_NH4 || patch[0].sat_NH4<0 ||
        patch[0].sat_DOC!=patch[0].sat_DOC || patch[0].sat_DOC<0 ||
        patch[0].surface_NO3!=patch[0].surface_NO3 || patch[0].surface_NO3<0 ||
        patch[0].surface_NH4!=patch[0].surface_NH4 || patch[0].surface_NH4<0 ||
        patch[0].surface_DOC!=patch[0].surface_DOC || patch[0].surface_DOC<0) printf("patch daily F9N after denitrif %d-%d-%d [%d,%d,%d,%d] [%e %e %e] [%e %e %e] [%e %e %e]\n",
              current_date.year, current_date.month, current_date.day,
              patch[0].ID, patch[0].drainage_type, actionPIPEDRAIN, actionSEWER,
              patch[0].soil_ns.nitrate,patch[0].soil_ns.sminn,patch[0].soil_cs.DOC,
              patch[0].sat_NO3,patch[0].sat_NH4,patch[0].sat_DOC,
              patch[0].surface_NO3,patch[0].surface_NH4,patch[0].surface_DOC);
        
        

        
    }//growth flag: decomp, nitrification, denitrification, sat_solute
    
    
    


	/* track variables for snow assimilation  */
	if (patch[0].snowpack.water_equivalent_depth > ZERO) {
		basin[0].snowpack.energy_deficit += patch[0].snowpack.energy_deficit * patch[0].area;
		basin[0].snowpack.surface_age += patch[0].snowpack.surface_age * patch[0].area;
		basin[0].snowpack.T += patch[0].snowpack.T * patch[0].area;
		basin[0].area_withsnow += patch[0].area;
		}

	/* track variables for fire spread */
	if (command_line[0].firespread_flag == 1) {
		patch[0].fire.et = (patch[0].fire_defaults[0][0].ndays_average*patch[0].fire.et  +  
		(patch[0].transpiration_sat_zone + patch[0].transpiration_unsat_zone
		+ patch[0].evaporation + patch[0].evaporation_surf 
		+ patch[0].exfiltration_unsat_zone + patch[0].exfiltration_sat_zone))/
		(patch[0].fire_defaults[0][0].ndays_average + 1);

		patch[0].fire.pet = (patch[0].fire_defaults[0][0].ndays_average*patch[0].fire.pet    
				+ patch[0].PET) / 
		(patch[0].fire_defaults[0][0].ndays_average + 1);
		}
	


	patch[0].soil_cs.totalc = ((patch[0].soil_cs.soil1c)
		+ (patch[0].soil_cs.soil2c) +	(patch[0].soil_cs.soil3c)
		+ (patch[0].soil_cs.soil4c));
	patch[0].totalc += ((patch[0].soil_cs.totalc) + (patch[0].litter_cs.litr1c)
		+ (patch[0].litter_cs.litr2c) + (patch[0].litter_cs.litr3c)
		+ (patch[0].litter_cs.litr4c));
	patch[0].soil_ns.totaln = ((patch[0].soil_ns.soil1n)
		+ (patch[0].soil_ns.soil2n) + (patch[0].soil_ns.soil3n)
		+ (patch[0].soil_ns.soil4n) + (patch[0].soil_ns.nitrate)
		+ (patch[0].soil_ns.sminn));
	patch[0].totaln += (patch[0].soil_ns.totaln + (patch[0].litter_ns.litr1n)
		+ (patch[0].litter_ns.litr2n) + (patch[0].litter_ns.litr3n)
		+ (patch[0].litter_ns.litr4n));
	patch[0].nitrogen_balance = patch[0].preday_totaln - patch[0].totaln - patch[0].ndf.N_to_gw
		+ zone[0].ndep_NO3 + zone[0].ndep_NH4 - patch[0].ndf.denitrif + fertilizer_NO3 + fertilizer_NH4;

	resp =  (patch[0].cdf.litr1c_hr + patch[0].cdf.litr2c_hr
		+ patch[0].cdf.litr4c_hr + patch[0].cdf.soil1c_hr
		+ patch[0].cdf.soil2c_hr + patch[0].cdf.soil3c_hr
		+ patch[0].cdf.soil4c_hr);
	patch[0].carbon_balance = patch[0].preday_totalc + patch[0].net_plant_psn
		- patch[0].totalc - (patch[0].cdf.litr1c_hr + patch[0].cdf.litr2c_hr
		+ patch[0].cdf.litr4c_hr + patch[0].cdf.soil1c_hr
		+ patch[0].cdf.soil2c_hr + patch[0].cdf.soil3c_hr
		+ patch[0].cdf.soil4c_hr);

	if (command_line[0].snow_scale_flag == 1)
	  patch[0].water_balance = zone[0].rain + zone[0].snow*patch[0].snow_redist_scale 
		+ patch[0].preday_detention_store +
		+ irrigation 
		+ patch[0].landuse_defaults[0][0].septic_water_load/patch[0].area
		+ zone[0].rain_hourly_total - ( patch[0].gw_drainage
		+ patch[0].transpiration_sat_zone + patch[0].transpiration_unsat_zone
		+ patch[0].evaporation + patch[0].evaporation_surf 
		+ patch[0].exfiltration_unsat_zone + patch[0].exfiltration_sat_zone)
		- (patch[0].rz_storage - patch[0].preday_rz_storage)		
		- (patch[0].unsat_storage - patch[0].preday_unsat_storage)
		- (patch[0].preday_sat_deficit - patch[0].sat_deficit)
		- patch[0].delta_snowpack - patch[0].delta_rain_stored
		- patch[0].delta_snow_stored - patch[0].detention_store;
	else	
	  patch[0].water_balance = zone[0].rain + zone[0].snow 
		+ patch[0].preday_detention_store +
		+ irrigation 
		+ patch[0].landuse_defaults[0][0].septic_water_load/patch[0].area
		+ zone[0].rain_hourly_total - ( patch[0].gw_drainage
		+ patch[0].transpiration_sat_zone + patch[0].transpiration_unsat_zone
		+ patch[0].evaporation + patch[0].evaporation_surf 
		+ patch[0].exfiltration_unsat_zone + patch[0].exfiltration_sat_zone)
		- (patch[0].rz_storage - patch[0].preday_rz_storage)			
		- (patch[0].unsat_storage - patch[0].preday_unsat_storage)
		- (patch[0].preday_sat_deficit - patch[0].sat_deficit)
		- patch[0].delta_snowpack - patch[0].delta_rain_stored
		- patch[0].delta_snow_stored - patch[0].detention_store;

	/*
	if ((patch[0].water_balance > 0.00000001)||
		(patch[0].water_balance < -0.00000001)){
		printf("\n Water Balance is %12.8f on %ld %ld %ld for patch %d of type %d",
			patch[0].water_balance,
			current_date.day,
			current_date.month,
			current_date.year,
			patch[0].ID,
			patch[0].drainage_type);
		printf("\nRain %lf %lf, dt %lf, rh %lf, T %lf %lf, E %lf ES %lf, Ex %lf %lf, RZ %lf %lf, US %lf %lf, SD %lf %lf, S %lf, RS %lf (%lf %lf), SS %lf, DT %lf",
		 zone[0].rain , zone[0].snow , patch[0].preday_detention_store ,
		 zone[0].rain_hourly_total ,  
		 patch[0].transpiration_sat_zone , patch[0].transpiration_unsat_zone
		, patch[0].evaporation , patch[0].evaporation_surf
		, patch[0].exfiltration_unsat_zone , patch[0].exfiltration_sat_zone
		, patch[0].rz_storage , patch[0].preday_rz_storage
		, patch[0].unsat_storage , patch[0].preday_unsat_storage
		, patch[0].preday_sat_deficit , patch[0].sat_deficit
		, patch[0].delta_snowpack , patch[0].delta_rain_stored
		, patch[0].preday_rain_stored, patch[0].rain_stored
		, patch[0].delta_snow_stored , patch[0].detention_store);	
	
	}
	 
	*/
	
	/* Calculate LE for surface evap */
	/* soil&litter&detstore evap x latent heat vaporization x water density */
	patch[0].LE_soil = (patch[0].evaporation_surf + patch[0].exfiltration_sat_zone 
						+ patch[0].exfiltration_unsat_zone) 
						* (2.5023e6 - 2430.54 * zone[0].metv.tday) / 1000 * 1000;
	
	
	/*---------------------------------------------------------------------*/
	/*	get rid of any negative soil or litter stores			*/
	/*---------------------------------------------------------------------*/

	if (command_line[0].grow_flag > 0)
		ch = check_zero_stores(
			&(patch[0].soil_cs),
			&(patch[0].soil_ns),
			&(patch[0].litter_cs),
			&(patch[0].litter_ns));

	if ( command_line[0].verbose_flag > 1 ) {
		printf("\n%ld %ld %ld  -335.2 ",
			current_date.year, current_date.month, current_date.day);
		printf("\n   %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f ",
			infiltration,
			patch[0].snowpack.water_equivalent_depth,
			patch[0].infiltration_excess,
			patch[0].transpiration_sat_zone + patch[0].transpiration_unsat_zone,
			patch[0].unsat_drainage,
			patch[0].unsat_storage,
			patch[0].cap_rise,
			patch[0].sat_deficit);
		printf("\n%ld %ld %ld  -335.3 ",
			current_date.year, current_date.month, current_date.day);
		printf("\n   %8.5f, %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f",
			zone[0].rain + zone[0].snow,
			patch[0].infiltration_excess,
			patch[0].transpiration_sat_zone + patch[0].transpiration_unsat_zone,
			patch[0].evaporation + patch[0].exfiltration_sat_zone
			+ patch[0].exfiltration_unsat_zone,
			(patch[0].unsat_storage - patch[0].preday_unsat_storage),
			(patch[0].preday_sat_deficit - patch[0].sat_deficit),
			patch[0].delta_snowpack,
			patch[0].delta_rain_stored + patch[0].delta_snow_stored);
	}
if ( command_line[0].verbose_flag == -5 ){
	printf("\n***END PATCH DAILY: exfil_unsat=%lf",patch[0].exfiltration_unsat_zone);
}

	return;
} /*end patch_daily_F.c*/
