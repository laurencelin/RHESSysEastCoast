/*--------------------------------------------------------------*/
/* 																*/
/*					construct_routing_topology					*/
/*																*/
/*	construct_routing_topology.c - creates a patch object		*/
/*																*/
/*	NAME														*/
/*	construct_routing_topology.c - creates a patch object		*/
/*																*/
/*	SYNOPSIS													*/
/*	struct routing_list_object construct_routing_topology( 		*/
/*							struct basin_object *basin)			*/
/*																*/
/* 																*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*  reads routing topology from input file						*/
/*	creates neighbourhood structure for each patch in the basin */
/*	returns a list giving order for patch-level routing			*/
/*																*/
/*																*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rhessys.h"

struct routing_list_object *construct_routing_topology(char *routing_filename,
		  struct basin_object *basin,
		  struct command_line_object *command_line,
		  bool surface)
													  
{
    
    double    compute_delta_water(
                                  int,
                                  double,
                                  double,
                                  double,
                                  double,
                                  double);
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	struct patch_object *find_patch(int, int, int, struct basin_object *);
	
	int assign_neighbours (struct neighbour_object *,
		int,
		struct basin_object *,
		FILE *);
	
	void *alloc(size_t, char *, char *);

	double * compute_transmissivity_curve( double, struct patch_object *, struct command_line_object *);
	
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int		i,d,j;
	int		num_patches, num_neighbours;
	int		patch_ID, zone_ID, hill_ID, aggregate_ID, aggregate_index;
	int		drainage_type;
	double	x,y,z, area, gamma, width, wttd;// --> constraintWaterTableTopDepth
	FILE	*routing_file;
	struct routing_list_object	*rlist;
	struct	patch_object	*patch;
	struct	patch_object	*stream;
	struct	innundation_object *innundation_list;
    int aggregateLength;
    struct soil_default *soilD;
   
	
	rlist = (struct routing_list_object	*)alloc( sizeof(struct routing_list_object), "rlist", "construct_routing_topology");
	
	/*--------------------------------------------------------------*/
	/*  Try to open the routing file in read mode.                    */
	/*--------------------------------------------------------------*/
	if ( (routing_file = fopen(routing_filename,"r")) == NULL ){
		fprintf(stderr,"FATAL ERROR:  Cannot open routing file %s\n",
			routing_filename);
		exit(EXIT_FAILURE);
	} /*end if*/
	fscanf(routing_file,"%d",&num_patches);
	rlist->num_patches = num_patches;
	rlist->list = (struct patch_object **)alloc(
		num_patches * sizeof(struct patch_object *), "patch list",
		"construct_routing_topography");
	
	/*--------------------------------------------------------------*/
	/*	Read in  each patch record and find it		.				*/
	/*	if it is a stream add it to the basin level routing list	*/
	/*	otherwise add it to the hillslope level routing list		*/
	/*--------------------------------------------------------------*/
    aggregateLength = 0;
	for (i=0; i< num_patches; ++i) {
        // under patch loop
        if(command_line[0].aggregate_flag>0){
            fscanf(routing_file,"%d %d %d %lf %lf %lf %lf %lf %d %lf %d %d %d",
                &patch_ID,
                &zone_ID,
                &hill_ID,
                &x,&y,&z,
                &area,
                &wttd,
                &drainage_type,
                &gamma,
                &num_neighbours,
                &aggregate_ID,  //[1 - 4]
                &aggregate_index); //[0, 3]
            //printf("now it reads %d: %d %d %d %lf %lf %lf %lf %lf %d %lf %d %d\n", i, patch_ID, zone_ID, hill_ID, x, y, z, area, area, drainage_type, gamma, num_neighbours, aggregate_ID);
            //fscanf does not read a line. it reads the "format string" as a block
            aggregateLength = max(aggregateLength,aggregate_index);
        }else{
            fscanf(routing_file,"%d %d %d %lf %lf %lf %lf %lf %d %lf %d",
                   &patch_ID,
                   &zone_ID,
                   &hill_ID,
                   &x,&y,&z,
                   &area,
                   &wttd,
                   &drainage_type,
                   &gamma,
                   &num_neighbours);
            aggregate_ID=0;
            aggregate_index=0;
        }
        
        // under patch loop
		if  ( (patch_ID != 0) && (zone_ID != 0) && (hill_ID != 0) )
			patch = find_patch(patch_ID, zone_ID, hill_ID, basin);
		else
			patch = basin[0].outside_region;
        
        if ( !surface ){
            patch[0].aggregate_ID = aggregate_ID;
            patch[0].aggregate_index = aggregate_index;
            patch[0].constraintWaterTableTopDepth = wttd;
            // 3m basement depth then 3*baseFrac + 0*(1-baseFrac) = wttd --> baseFrac = wttd/3.0;
            patch[0].basementFrac = wttd/3.0;
            patch[0].constraintWaterTableTopDepth_def = patch[0].basementFrac * compute_delta_water(
                                                                              command_line[0].verbose_flag,
                                                                              patch[0].soil_defaults[0][0].porosity_0,
                                                                              patch[0].soil_defaults[0][0].porosity_decay,
                                                                              patch[0].soil_defaults[0][0].soil_depth,
                                                                              3.0, 0.0);
            
//            patch[0].constraintWaterTableTopDepth_def = compute_delta_water(
//                                                                            command_line[0].verbose_flag,
//                                                                            patch[0].soil_defaults[0][0].porosity_0,
//                                                                            patch[0].soil_defaults[0][0].porosity_decay,
//                                                                            patch[0].soil_defaults[0][0].soil_depth,
//                                                                            patch[0].constraintWaterTableTopDepth, 0.0);
            
            
            
            patch[0].horizontal_k_SCALE = 1.0; //area;
            
//            patch[0].rootzone.NO3decayRate = 0.12; //k*0.01;
//            patch[0].rootzone.NH4decayRate = 0.35; //k*0.01;
//            patch[0].rootzone.DOMdecayRate = 0.12 ; //k*0.01; 0.35, 0.12
//            patch[0].surface_NO3 = 0.0;
//            patch[0].surface_NH4 = 0.0;
//            patch[0].surface_DOC = 0.0;
//            patch[0].surface_DON = 0.0;

            // assume soil OM is 0.35 decay; (patch averaged) basement; disturbed soil lost 60% OM is set in worldfile.
            double remainPerc = (exp(-0.12*wttd)-exp(-0.12*patch[0].soil_defaults[0][0].soil_depth))/(1.0-exp(-0.12*patch[0].soil_defaults[0][0].soil_depth));
            
            patch[0].soil_ns.sminn *= remainPerc;
            patch[0].soil_ns.nitrate *= remainPerc;
            
            patch[0].soil_cs.soil1c *= remainPerc;
            patch[0].soil_cs.soil2c *= remainPerc;
            patch[0].soil_cs.soil3c *= remainPerc;
            patch[0].soil_cs.soil4c *= remainPerc;
            
            patch[0].soil_ns.soil1n *= remainPerc;
            patch[0].soil_ns.soil2n *= remainPerc;
            patch[0].soil_ns.soil3n *= remainPerc;
            patch[0].soil_ns.soil4n *= remainPerc;
            
            // veg root for water (does not work!!)
            //patch[0].rootzone.depth += wttd*10; // it affects nitrif, denitrif, and decomp processes
            
            
        }// not surface
		rlist->list[i] = patch;

        // under patch loop
        soilD = &patch[0].soil_defaults[0][0];
		if ((soilD->Ksat_0 < ZERO))
			printf("\n WARNING lateral Ksat (%lf) are close to zero for patch %d",
				soilD->Ksat_0, patch[0].ID);
		
        
		if (soilD->m < ZERO)
		 	gamma = gamma * soilD->Ksat_0;
		else	
		 	gamma = gamma * soilD->m * soilD->Ksat_0;

		/*--------------------------------------------------------------*/
		/*  Allocate innundation list array				*/
		/*	note for this routing there is only one innundation depth 	*/
		/*	however it is need to be compatablability 		*/
		/*--------------------------------------------------------------*/
        // under patch loop
		d=0;
		if ( surface ) {
            //surface
			patch->surface_innundation_list = (struct innundation_object *)alloc( 1 *
								sizeof(struct innundation_object), "surface_innundation_list", "construct_routing_topology");
			innundation_list = patch->surface_innundation_list;
		} else {
            // subsurface
			patch->innundation_list = (struct innundation_object *)alloc( 1 *
					sizeof(struct innundation_object), "innundation_list", "construct_routing_topology");
			innundation_list = patch->innundation_list;
		}

		if ( surface ) {
			patch[0].num_innundation_depths = 1;
		}

		innundation_list->num_neighbours = num_neighbours;
		innundation_list->gamma = gamma;
		// TODO: what should critical depth be for a surface flow table?
		innundation_list->critical_depth = NULLVAL;

		if ( !surface ) {
            // <<----------------- very important: patch[0].drainage_type is defined by subsurface flowtable
            // "stream_gamma" is no use in the model so far
			patch[0].stream_gamma = 0.0;
			patch[0].drainage_type = drainage_type;
			if ( (patch[0].drainage_type != STREAM) && (patch[0].innundation_list[d].gamma < ZERO) ) {
				printf(
						"\n non-stream patches with zero gamma %d switched to stream for now",
						patch[0].ID);
				patch[0].drainage_type = STREAM;
			}
		}

        // under patch loop
		/*--------------------------------------------------------------*/
		/*  Allocate neighbour array									*/
		/*--------------------------------------------------------------*/
		innundation_list->neighbours = (struct neighbour_object *)alloc(num_neighbours *
				sizeof(struct neighbour_object), "neighbours", "construct_routing_topology");
		num_neighbours = assign_neighbours(innundation_list->neighbours, num_neighbours, basin, routing_file);
		if ((num_neighbours == -9999) && (patch[0].drainage_type != STREAM)) {
			printf("\n WARNING sum of patch %d neigh gamma is not equal to 1.0", patch[0].ID); 
		} else {
			innundation_list->num_neighbours = num_neighbours;
		}

		if ( drainage_type == ROAD) {
            // road patch has additional neighbour row in the flow table (not included in the neighbour count
			fscanf(routing_file,"%d %d %d %lf",
				&patch_ID,
				&zone_ID,
				&hill_ID,
				&width);
			// TODO: Decide if we need separate stream_gamma, road_cut_depth, and next_stream values for surface flow table
			if ( !surface ) {
                // subsurface only!
				patch[0].stream_gamma = gamma;
				patch[0].road_cut_depth = width * tan(patch[0].slope);
				stream = find_patch(patch_ID, zone_ID, hill_ID, basin);
				patch[0].next_stream = stream;
			}
		}

        // under patch loop
		if ( !surface ) {
			/*--------------------------------------------------------------*/
			/*	create a vector of transmssivities 			*/
			/*--------------------------------------------------------------*/
			patch[0].num_soil_intervals = (int) lround(patch[0].soil_defaults[0][0].soil_water_cap / patch[0].soil_defaults[0][0].interval_size);
			if (patch[0].num_soil_intervals > MAX_NUM_INTERVAL) {
				patch[0].num_soil_intervals = MAX_NUM_INTERVAL;
				patch[0].soil_defaults[0][0].interval_size = patch[0].soil_defaults[0][0].soil_water_cap / MAX_NUM_INTERVAL;
				}
			patch[0].transmissivity_profile = compute_transmissivity_curve(gamma, patch, command_line);
        }

        command_line[0].outletPatchID = patch[0].ID;// shortcut to save this info on commandline obj; ideally it would be basin.
	}// end of i for loop (scanning flow table)
    aggregateLength++; //<<---------------
    basin->aggregateLength = aggregateLength;//<<---------------

	fclose(routing_file);

	return(rlist);
} /*end construct_routing_topology.c*/

