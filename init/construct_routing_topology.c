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
    int assign_drainIN (struct drainIN_object *,
        int,
        struct basin_object *,
        FILE *,
        struct patch_object *);
        
	void *alloc(size_t, char *, char *);

	
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int		i,d,j;
    int		num_patches, num_neighbours, num_drainIN_septic, num_drainIN_irrigation;
    double  fnum_drainIN_septic, fnum_drainIN_irrigation;
	int		patch_ID, zone_ID, hill_ID, aggregate_ID, aggregate_index;
	int		drainage_type;
	double	gamma, width, wttd;// --> constraintWaterTableTopDepth
    double waterFrac, sewerFrac;
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
                &waterFrac,
                &sewerFrac,
                &fnum_drainIN_septic,
                &fnum_drainIN_irrigation,
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
                   &waterFrac,
                   &sewerFrac,
                   &fnum_drainIN_septic,
                   &fnum_drainIN_irrigation,
                   &wttd,
                   &drainage_type,
                   &gamma,
                   &num_neighbours);
            aggregate_ID=0;
            aggregate_index=0;
        }//end of if
        
        // under patch loop
		if  ( (patch_ID != 0) && (zone_ID != 0) && (hill_ID != 0) )
			patch = find_patch(patch_ID, zone_ID, hill_ID, basin);
		else
			patch = basin[0].outside_region;
        
        if ( !surface ){
            patch[0].aggregate_ID = aggregate_ID;
            patch[0].aggregate_index = aggregate_index;
            patch[0].waterFrac = (waterFrac>0 && waterFrac<=1? waterFrac : 0.0);
            patch[0].sewerFrac = (sewerFrac>0 && sewerFrac<=1? sewerFrac : 0.0);
            // basement implementation is disabled
            // 3m basement depth then 3*baseFrac + 0*(1-baseFrac) = wttd --> baseFrac = wttd/3.0;
            //patch[0].constraintWaterTableTopDepth = wttd;
            //patch[0].basementFrac = 0.0; //min(1.0,wttd/3.0);
            //patch[0].constraintWaterTableTopDepth_def = 0.0;
            
            if(patch[0].sat_def_head>0){
                //do nothing
            }else if(wttd>0){
                //need to convert depth to sat_def
                int rtz2_index = (int)(round(wttd*1000));
                patch[0].sat_def_head = patch[0].soil_defaults[0][0].rtz2sat_def_0z[rtz2_index];
                //patch[0].rootdepth_index = patch[0].soil_defaults[0][0].rtz2sat_def_pct_index[patch[0].rtz2_index]; // do i need this?
            }
        }// not surface
		rlist->list[i] = patch;

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
        
        fnum_drainIN_septic = (fnum_drainIN_septic<0? -fnum_drainIN_septic : 0.0);
        fnum_drainIN_irrigation = (fnum_drainIN_irrigation<0? -fnum_drainIN_irrigation : 0.0);
        num_drainIN_septic = (int) fnum_drainIN_septic;
        num_drainIN_irrigation = (int) fnum_drainIN_irrigation;
        innundation_list->num_drainIN_septic = num_drainIN_septic;
        innundation_list->num_drainIN_irrigation = num_drainIN_irrigation;
		
//        if ( !surface ) printf("(%d %d %d %lf %lf %lf %lf %lf %d %lf %d)\n",
//                               patch_ID,
//                               zone_ID,
//                               hill_ID,
//                               x,y,
//                               fnum_drainIN_septic,
//                               fnum_drainIN_irrigation,
//                               wttd,
//                               drainage_type,
//                               gamma,
//                               num_neighbours,
//                               innundation_list->num_drainIN_septic,
//                               innundation_list->num_drainIN_irrigation);

//        // not sure what is below.
//        if ( !surface ) {
//            // <<----------------- very important: patch[0].drainage_type is defined by subsurface flowtable
//            // "stream_gamma" is no use in the model so far
//			patch[0].stream_gamma = 0.0;
//			patch[0].drainage_type = drainage_type;
//			if ( (patch[0].drainage_type != STREAM) && (patch[0].innundation_list[d].gamma < ZERO) ) {
//				printf(
//						"\n non-stream patches with zero gamma %d switched to stream for now (%d %d %d %lf %lf %lf %lf %lf %d %lf %d)",
//						patch[0].ID,
//                       patch_ID,
//                       zone_ID,
//                       hill_ID,
//                       waterFrac,
//                       sewerFrac,
//                       fnum_drainIN_septic,
//                       fnum_drainIN_irrigation,
//                       wttd,
//                       drainage_type,
//                       gamma,
//                       num_neighbours);
//				patch[0].drainage_type = STREAM;
//			}//end of if
//		}// end of if

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
			// note: Decide if we need separate stream_gamma, road_cut_depth, and next_stream values for surface flow table
			if ( !surface ) {
                // subsurface only!
				patch[0].stream_gamma = gamma;
				patch[0].road_cut_depth = width * tan(patch[0].slope);
                patch[0].road_cut_depth_def = patch[0].soil_defaults[0][0].rtz2sat_def_0z[(int)(patch[0].road_cut_depth*1000)];
				stream = find_patch(patch_ID, zone_ID, hill_ID, basin);
				patch[0].next_stream = stream;
			}
		}//end of "road" if

        // 3) reading the drainIN rows information
        if(num_drainIN_septic>0){
            // septic
            innundation_list->drainIN_septic = (struct drainIN_object *)alloc(num_drainIN_septic *
                    sizeof(struct drainIN_object), "drainIN_septic", "construct_routing_topology");
            num_drainIN_septic = assign_drainIN(innundation_list->drainIN_septic, num_drainIN_septic, basin, routing_file,patch);
            if (num_drainIN_septic == -9999) { printf("WARNING drainIN_septic error", patch[0].ID); }
        }else{
            innundation_list->drainIN_septic = NULL;
        }// end of if "num_drainIN_septic"

        if(num_drainIN_irrigation>0){
            // septic
            innundation_list->drainIN_irrigation = (struct drainIN_object *)alloc(num_drainIN_irrigation *
                    sizeof(struct drainIN_object), "drainIN_irrigation", "construct_routing_topology");
            num_drainIN_irrigation = assign_drainIN(innundation_list->drainIN_irrigation, num_drainIN_irrigation, basin, routing_file,patch);
            if (num_drainIN_irrigation == -9999) { printf("WARNING drainIN_irrigation error", patch[0].ID); }
        }else{
            innundation_list->drainIN_irrigation = NULL;
        }// end of if "num_drainIN_septic"
        

        command_line[0].outletPatchID = patch[0].ID;// shortcut to save this info on commandline obj; ideally it would be basin.
	}// end of i for loop (scanning flow table)
    aggregateLength++; //<<---------------
    basin->aggregateLength = aggregateLength;//<<---------------

	fclose(routing_file);

	return(rlist);
} /*end construct_routing_topology.c*/

