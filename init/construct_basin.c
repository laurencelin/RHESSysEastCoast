/*--------------------------------------------------------------*/
/* 																*/
/*					construct_basin								*/
/*																*/
/*	construct_basin.c - creates a basin object					*/
/*																*/
/*	NAME														*/
/*	construct_basin.c - creates a basin object					*/
/*																*/
/*	SYNOPSIS													*/
/*	void construct_basin(										*/
/*			struct	command_line_object	*command_line,			*/
/*			FILE	*world_file									*/
/*			int		num_world_base_stations,					*/
/*			struct base_station_object	**world_base_stations,	*/
/*			struct basin_object	**basin_list,					*/
/*			struct default_object *defaults)					*/
/* 																*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	Constructs the basin object which consists of:				*/
/*		- basin specific parameters and identification			*/
/*		- a possible extension to a grow object					*/
/*		- a list of hillslopes in the basin.					*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	Basins dont own climate files since all of their hillslopes	*/
/*	own them instead.  I guess this means we have to compute	*/
/*	a lot of local climate info but perhaps this is more		*/
/*	correct anyways.  As usual, it is up to the user to         */
/*	aggregate model  output from each hillslope if they want	*/
/*	basin values.												*/
/*																*/
/*	We use a list of pointers to hillslope objects rather than	*/
/*	a contiguous array of hillslope objects with a pointer to	*/
/*	the head of the array.  The list of pointers is a bit less	*/
/*	efficient since the pointer must be placed in the heap (RAM)*/
/*	at the start of each object BUT								*/
/*																*/
/*		1.  We can dynamically add and remove hillslopes.		*/
/*		2.  Most of the processing time is required on the 		*/
/*			sub-hillslope basis so the repositioning of pointers*/
/*			will not be too drastic if it is limited to the 	*/
/*			hillslope level or up.								*/
/*		3.  We will be able to make use of smaller chunks of 	*/
/*			RAM.  												*/
/*	Original code, January 16, 1996.							*/
/*	May 7, 1997	C.Tague											*/
/* 		- added a routine to sort hierarchy by elevation 		*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rhessys.h"
#include "functions.h"
#include "params.h"

struct basin_object *construct_basin(
									 struct	command_line_object	*command_line,
									 FILE	*world_file,
									 int	*num_world_base_stations,
									 struct base_station_object	**world_base_stations,
									 struct	default_object	*defaults,
									 struct base_station_ncheader_object *base_station_ncheader,
									 struct world_object *world)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	struct base_station_object *assign_base_station(
								int,
								int,
								struct base_station_object **);
	
	struct hillslope_object *construct_hillslope(
		struct	command_line_object *,
		FILE    *,
		int		*,
		struct base_station_object **,
		struct	default_object *,
		struct base_station_ncheader_object *,
		struct world_object *);
	
	void	*alloc( 	size_t, char *, char *);
	
	void	sort_by_elevation( struct basin_object *);
	
//	struct routing_list_object construct_ddn_routing_topology(
//		char *,
//		struct basin_object *);

//	struct routing_list_object construct_routing_topology(
//		char *,
//		struct basin_object *,
//		struct	command_line_object *);
	
	struct stream_list_object construct_stream_routing_topology(
		char *,
		struct basin_object *, 
		struct	command_line_object *);
	
    long  julday(struct date);
    
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int	base_stationID;
	int		i,j,z;
	int		default_object_ID;
	double		check_snow_scale;
	double		n_routing_timesteps;
	char		record[MAXSTR];
	struct basin_object	*basin;
    param    *paramPtr=NULL;
    int    paramCnt=0;
    
	/*--------------------------------------------------------------*/
	/*	Allocate a basin object.								*/
	/*--------------------------------------------------------------*/
    basin = (struct basin_object *) alloc( 1 *
		sizeof( struct basin_object ),"basin","construct_basin");
	
	/*--------------------------------------------------------------*/
	/*	Read in the basinID.									*/
	/*--------------------------------------------------------------*/
//	fscanf(world_file,"%d",&(basin[0].ID));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(basin[0].x));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(basin[0].y));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(basin[0].z));
//	read_record(world_file, record);
//	fscanf(world_file,"%d",&(default_object_ID));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(basin[0].latitude));
//	read_record(world_file, record);
//	fscanf(world_file,"%d",&(basin[0].num_base_stations));
//	read_record(world_file, record);
	
    paramPtr=readtag_worldfile(&paramCnt,world_file,"Basin");
    
    basin[0].ID = getIntWorldfile(&paramCnt,&paramPtr,"basin_ID","%d",1,1);
    basin[0].x = getDoubleWorldfile(&paramCnt,&paramPtr,"x","%lf",0.0,1);
    basin[0].y = getDoubleWorldfile(&paramCnt,&paramPtr,"y","%lf",0.0,1);
    basin[0].z = getDoubleWorldfile(&paramCnt,&paramPtr,"z","%lf",-9999,0);
    default_object_ID = getIntWorldfile(&paramCnt,&paramPtr,"basin_parm_ID","%d",1,1);
    basin[0].latitude = getDoubleWorldfile(&paramCnt,&paramPtr,"latitude","%lf",-9999,0);
    basin[0].num_base_stations = getIntWorldfile(&paramCnt,&paramPtr,"basin_n_basestations","%d",0,1);
    
    
	/*--------------------------------------------------------------*/
	/*	Create cosine of latitude to save future computations.		*/
	/*--------------------------------------------------------------*/
	basin[0].cos_latitude = cos(basin[0].latitude*DtoR);
	basin[0].sin_latitude = sin(basin[0].latitude*DtoR);
	
	/*--------------------------------------------------------------*/
	/*    Allocate a list of base stations for this basin.			*/
	/*--------------------------------------------------------------*/
    basin[0].base_stations = NULL;
//    basin[0].base_stations = (struct base_station_object **) alloc(basin[0].num_base_stations * sizeof(struct base_station_object *),"base_stations","construct_basin");
    
    /*--------------------------------------------------------------*/
    /*    calculate 21 day running average for this basin; assuming only one basin for the worldfile	*/
    /*--------------------------------------------------------------*/
    for (i=0; i<world[0].num_base_stations; i++ ) {
        // duration is the size of the time series vectors/arrays
        //world[0].start_date,
        //world[0].duration
        double tmpP;
        int ii;
        int currentJ = julday(world[0].start_date);
        int currentJy = world[0].start_date.year;
        
        // world_base_stations = world[0].base_stations
        world_base_stations[i]->daily_clim[0].pretmin21ravg = (double *) alloc(world[0].duration.day*sizeof(double), "sequence","construct_clim_sequence");
        world_base_stations[i]->daily_clim[0].predaylestimate = (double *) alloc(world[0].duration.day*sizeof(double), "sequence","construct_clim_sequence");
        world_base_stations[i]->daily_clim[0].predayl21ravg = (double *) alloc(world[0].duration.day*sizeof(double), "sequence","construct_clim_sequence");
        
        double JulianCentury;
        double GeomMean_long_sun_deg;
        double GeomMean_Anom_sun_deg;
        double Sun_eq_Ctr;
        double Sun_True_long_deg;
        double Sun_app_long_deg;
        double mean_Obliq_Ecliptic_deg;
        double sixtyth = 1.0/60.0;
        double Obliq_Corr_deg;
        double Sun_delin_deg;
        double HA_sunrise_deg;
        double sunlight_min;
        
        for(ii=0; ii<world[0].duration.day; ii++){
            if(ii<10 || ii+10>=world[0].duration.day){world_base_stations[i]->daily_clim[0].pretmin21ravg[ii] = world_base_stations[i]->daily_clim[0].tmin[ii]; }
            else{
                world_base_stations[i]->daily_clim[0].pretmin21ravg[ii] = (
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii-10]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii-9]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii-8]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii-7]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii-6]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii-5]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii-4]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii-3]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii-2]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii-1]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii+1]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii+2]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii+3]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii+4]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii+5]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii+6]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii+7]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii+8]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii+9]+
                                                                           world_base_stations[i]->daily_clim[0].tmin[ii+10])*0.04761905;
            }
            
            /*
             http://mathforum.org/library/drmath/view/56478.html
            _Ecological Modeling_, volume 80 (1995) pp. 87-95, called "A Model
            Comparison for Daylength as a Function of Latitude and Day of the
                Year."
             
             This article presented a model that apparently does a very good
             job of estimating the daylight - the error is less than one minute
             within 40 degrees of the equator, and less than seven minutes within
             60 degrees and usually within two minutes for these latitudes.
             
             */
//            tmpP = asin(.39795*cos(.2163108 + 2*atan(.9671396*tan(.00860*(currentJ-186)))));
//            world_base_stations[i]->daily_clim[0].predaylestimate[ii] = 24 - 7.639437 * acos( (0.01454332 + basin[0].sin_latitude*sin(tmpP))/(basin[0].cos_latitude*cos(tmpP)) );
//            world_base_stations[i]->daily_clim[0].predaylestimate[ii] *= 3600.0;
            
            /*
             NOAA calculation https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
             century starts from 2000 Jan 1st
             */
            
            JulianCentury =(currentJ + (currentJy-1)*365.25-2451545)*2.737851e-05;
            GeomMean_long_sun_deg = fmod( (280.46646+JulianCentury*(36000.76983+JulianCentury*0.0003032)), 360.0);
            GeomMean_Anom_sun_deg = 357.52911+ JulianCentury*(35999.05029 - 0.0001537* JulianCentury);
            Sun_eq_Ctr = sin(GeomMean_Anom_sun_deg*DtoR)*(1.914602-JulianCentury*(0.004817+0.000014* JulianCentury))+sin((2* GeomMean_Anom_sun_deg)*DtoR)*(0.019993-0.000101* JulianCentury)+sin((3* GeomMean_Anom_sun_deg)*DtoR)*0.000289;
            Sun_True_long_deg = GeomMean_long_sun_deg + Sun_eq_Ctr;
            Sun_app_long_deg = Sun_True_long_deg-0.00569-0.00478*sin( (125.04-1934.136*JulianCentury)*DtoR );
            mean_Obliq_Ecliptic_deg = 23.0+(26.0+((21.448-JulianCentury*(46.815+ JulianCentury*(0.00059-JulianCentury*0.001813))))*sixtyth)*sixtyth;
            Obliq_Corr_deg = mean_Obliq_Ecliptic_deg + 0.00256*cos((125.04-1934.136* JulianCentury)*DtoR);
            Sun_delin_deg = asin( sin(Obliq_Corr_deg*DtoR)*sin(Sun_app_long_deg*DtoR) )*RtoD;
            HA_sunrise_deg= acos( cos(90.833*DtoR)/(basin[0].cos_latitude*cos(Sun_delin_deg*DtoR)) - basin[0].sin_latitude/basin[0].cos_latitude*tan(Sun_delin_deg*DtoR) )*RtoD;
            sunlight_min = (HA_sunrise_deg - 2.076*sqrt(basin[0].z)*sixtyth )*8;
            world_base_stations[i]->daily_clim[0].predaylestimate[ii] = sunlight_min*60.0; // second
            

            currentJ ++;
            if(currentJ>365 && currentJy%4!=0){currentJ=1; currentJy++;}
            if(currentJ>366 && currentJy%4==0){currentJ=1; currentJy++;}
        }//ii
        
        for(ii=0; ii<world[0].duration.day; ii++){
            if(ii<10 || ii+10>=world[0].duration.day){world_base_stations[i]->daily_clim[0].predayl21ravg[ii] = world_base_stations[i]->daily_clim[0].predaylestimate[ii]; }
            else{
                world_base_stations[i]->daily_clim[0].predayl21ravg[ii] = (
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii-10]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii-9]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii-8]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii-7]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii-6]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii-5]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii-4]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii-3]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii-2]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii-1]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii+1]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii+2]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii+3]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii+4]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii+5]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii+6]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii+7]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii+8]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii+9]+
                                                                           world_base_stations[i]->daily_clim[0].predaylestimate[ii+10])*0.04761905;
            }
        }//ii
    } /*end for*/
    
    
	/*--------------------------------------------------------------*/
	/*      Read each base_station ID and then point to that base_statio*/
	/*--------------------------------------------------------------*/
//	for (i=0 ; i<basin[0].num_base_stations; i++) {
//
//
//		fscanf(world_file,"%d",&(base_stationID));
//		read_record(world_file, record);
//
//		/*--------------------------------------------------------------*/
//		/*	Point to the appropriate base station in the base       	*/
//		/*              station list for this world.					*/
//		/*--------------------------------------------------------------*/
//		basin[0].base_stations[i] = assign_base_station(
//			base_stationID,
//			*num_world_base_stations,
//			world_base_stations);
//
//	} /*end for*/
    
	/*--------------------------------------------------------------*/
	/*	Create the grow subobject if needed.						*/
	/*--------------------------------------------------------------*/
	if ( command_line[0].grow_flag == 1 ){
		/*--------------------------------------------------------------*/
		/*		Allocate memory for the grow subobject.					*/
		/*--------------------------------------------------------------*/
		basin[0].grow = (struct grow_basin_object *)
			alloc(1 * sizeof(struct grow_basin_object),
			"grow","construct_basin");
		/*--------------------------------------------------------------*/
		/*	NOTE:  PUT READS FOR GROW SUBOBJECT HERE.					*/
		/*--------------------------------------------------------------*/
	} /*end if*/
	/*--------------------------------------------------------------*/
	/*  Assign  defaults for this basin                             */
	/*--------------------------------------------------------------*/
	basin[0].defaults = (struct basin_default **)
		alloc( sizeof(struct basin_default *),"defaults","construct_basin" );
	i = 0;
	while (defaults[0].basin[i].ID != default_object_ID) {
		i++;
		/*--------------------------------------------------------------*/
		/*  Report an error if no match was found.  Otherwise assign    */
		/*  the default to point to this basin.                         */
		/*--------------------------------------------------------------*/
		if ( i>= defaults[0].num_basin_default_files ){
			fprintf(stderr,
				"\nFATAL ERROR: in construct_basin,basin default ID %d not found.\n",
				default_object_ID);
			exit(EXIT_FAILURE);
		}
	} /* end-while */
	basin[0].defaults[0] = &defaults[0].basin[i];

	/*--------------------------------------------------------------*/
	/*	Read in the number of hillslopes.						*/
	/*--------------------------------------------------------------*/
//	fscanf(world_file,"%d",&(basin[0].num_hillslopes));
//	read_record(world_file, record);
	basin[0].num_hillslopes = getIntWorldfile(&paramCnt,&paramPtr,"NUM_of_","%d",0,0);
	/*--------------------------------------------------------------*/
	/*	Allocate a list of pointers to hillslope objects.			*/
	/*--------------------------------------------------------------*/
	basin[0].hillslopes = (struct hillslope_object **)
		alloc(basin[0].num_hillslopes * sizeof(struct hillslope_object *),
		"hillslopes","construct_basin");

	basin[0].area = 0.0;
	basin[0].max_slope = 0.0;
	n_routing_timesteps = 0.0;
	check_snow_scale = 0.0;
	/*--------------------------------------------------------------*/
	/*	Construct the hillslopes for this basin.					*/
	/*--------------------------------------------------------------*/
	for (i=0; i<basin[0].num_hillslopes; i++){
        // hillslope before drainage type is read in
		basin[0].hillslopes[i] = construct_hillslope(
			command_line, world_file, num_world_base_stations,
			world_base_stations, defaults, base_station_ncheader, world);
		basin[0].area += basin[0].hillslopes[i][0].area;
		n_routing_timesteps += basin[0].hillslopes[i][0].area *
			basin[0].hillslopes[i][0].defaults[0][0].n_routing_timesteps;
		if (basin[0].max_slope < basin[0].hillslopes[i][0].slope)
			basin[0].max_slope = basin[0].hillslopes[i][0].slope;
		if (command_line[0].snow_scale_flag == 1) {
			for (z = 0; z < basin[0].hillslopes[i][0].num_zones; z++) {
				for (j=0; j < basin[0].hillslopes[i][0].zones[z][0].num_patches; j++) { 
				check_snow_scale +=
				basin[0].hillslopes[i][0].zones[z][0].patches[j][0].snow_redist_scale *
				basin[0].hillslopes[i][0].zones[z][0].patches[j][0].area;
				}
			}	
		}
	};
	
	basin[0].defaults[0][0].n_routing_timesteps = 
			(int) (n_routing_timesteps / basin[0].area);

	if (basin[0].defaults[0][0].n_routing_timesteps < 1)
		basin[0].defaults[0][0].n_routing_timesteps = 1;
 
	if (command_line[0].snow_scale_flag == 1) {
		check_snow_scale /= basin[0].area;
		if (fabs(check_snow_scale - 1.0) > ZERO	) {
			printf("\n *******  WARNING  ********** ");
			printf("\n Basin-wide  average snow scale is %lf", check_snow_scale);
			printf("\n Snow rescaling will alter net precip input by this scale factor\n\n");
			}
		if (command_line[0].snow_scale_tol > ZERO) {
			if ((check_snow_scale > command_line[0].snow_scale_tol) || 
				(check_snow_scale < 1/command_line[0].snow_scale_tol)) {
				printf("Basin-wide  average snow scale %lf is outside tolerance %lf", 
				check_snow_scale, command_line[0].snow_scale_tol);
				printf("\n Exiting\n");
				exit(EXIT_FAILURE);
				}
			}
	}

	

	/*--------------------------------------------------------------*/
	/*	Sort sub-hierarchy in the basin by elevation				*/
	/*--------------------------------------------------------------*/
	sort_by_elevation(basin);

	/*--------------------------------------------------------------*/
	/*	Read in flow routing topology for routing option	*/
	/*--------------------------------------------------------------*/
	if ( command_line[0].routing_flag == 1 ) {
		basin[0].outside_region = (struct patch_object *) alloc (1 *
			sizeof(struct patch_object) , "patch",
			"construct_basin");
		basin[0].outside_region[0].sat_deficit = 0.0;
		basin[0].outside_region[0].ID = 0;
		if ( command_line[0].ddn_routing_flag == 1 ) {
			basin->route_list = construct_ddn_routing_topology( command_line[0].routing_filename, basin);
		} else {
			basin->route_list = construct_routing_topology( command_line[0].routing_filename, basin, command_line, false);
			
            if ( command_line->surface_routing_flag ) {
                
				printf("\tReading surface routing table\n");
				basin->surface_route_list = construct_routing_topology( command_line->surface_routing_filename, basin, command_line, true);
                
				if ( basin->surface_route_list->num_patches != basin->route_list->num_patches ) {
					fprintf(stderr,
							"\nFATAL ERROR: in construct_basin, surface routing table has %d patches, but subsurface routing table has %d patches. The number of patches must be identical.\n",
							basin->surface_route_list->num_patches, basin->route_list->num_patches);
					exit(EXIT_FAILURE);
				}
			} else {
				// No surface routing table specified, use sub-surface for surface
				basin->surface_route_list =
						construct_routing_topology( command_line->routing_filename, basin,
								command_line, true);
			}
		}// done with reading flow table
        
        for (i=0; i<basin[0].num_hillslopes; i++){
            // fix riparian zone in hillslopes
            for (z = 0; z < basin[0].hillslopes[i][0].num_zones; z++) {
                for (j=0; j < basin[0].hillslopes[i][0].zones[z][0].num_patches; j++) {

                    if (basin[0].hillslopes[i][0].zones[z][0].patches[j][0].drainage_type>0 && basin[0].hillslopes[i][0].zones[z][0].patches[j][0].drainage_type % actionRIPARIAN==0)
                    basin[0].hillslopes[i][0].riparian_area += basin[0].hillslopes[i][0].zones[z][0].patches[j][0].area;

                }// end of for loop j
            }// end of for loop z
            printf("hillslope riparian: %d %lf\n",basin[0].hillslopes[i][0].ID,basin[0].hillslopes[i][0].riparian_area);
        }// end of for loop i
        
        
	} else { // command_line[0].routing_flag != 1
		// For TOPMODEL mode, make a dummy route list consisting of all patches
		// in the basin, in no particular order.
		basin->route_list = construct_topmodel_patchlist(basin);
	}
	
	/*--------------------------------------------------------------*/
	/*	Read in stream routing topology if needed	*/
	/*--------------------------------------------------------------*/
	if ( command_line[0].stream_routing_flag == 1) {
			basin[0].stream_list = construct_stream_routing_topology( command_line[0].stream_routing_filename, basin,
						command_line);
	}
		else { 
			basin[0].stream_list.stream_network = NULL;
			basin[0].stream_list.streamflow = 0.0;
		}

	return(basin);
} /*end construct_basin.c*/
