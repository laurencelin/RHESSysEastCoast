/*--------------------------------------------------------------*/
/*                                                              */ 
/*		assign_neighbours									*/
/*                                                              */
/*  NAME                                                        */
/*		assign_neighbours									*/
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/*  assign_neighbours( struct patch_object *patch)				*/
/*						int num_neighbours,						*/
/*						FILE *routing_file)						*/
/*                                                              */
/*  OPTIONS                                                     */
/*                                                              */
/*  DESCRIPTION                                                 */
/*                                                              */
/*                                                              */
/*	assigns pointers to neighbours of each patch				*/
/*	as given in topology input file								*/
/*                                                              */
/*  PROGRAMMER NOTES                                            */
/*                                                              */
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"
int assign_neighbours( struct neighbour_object *neighbours,
					   int num_neighbours,
					   struct basin_object *basin,
					   FILE *routing_file)
{
	/*--------------------------------------------------------------*/
	/*  Local function declaration                                  */
	/*--------------------------------------------------------------*/
	void *alloc (size_t, char *, char *);
	struct patch_object *find_patch( int, int, int,
		struct basin_object *);
	
	/*--------------------------------------------------------------*/
	/*  Local variable definition.                                  */
	/*--------------------------------------------------------------*/
	int i, inx, new_num_neighbours;
	int patch_ID, zone_ID, hill_ID;
	double gamma;
    double edgedistance, edge;
	double sum_gamma;
	struct patch_object *neigh;
	
	/*--------------------------------------------------------------*/
	/*  find and assign each neighbour to array						*/
	/*--------------------------------------------------------------*/
	inx = 0;
	sum_gamma = 0.0;
    edgedistance = 0.0;
    edge = 0.0;
	new_num_neighbours = num_neighbours;
	for (i=0; i< num_neighbours; i++) {
        
		fscanf(routing_file,"%d %d %d %lf %lf %lf",
			&patch_ID,
			&zone_ID,
			&hill_ID,
			&gamma,
            &edgedistance,
            &edge);
		/*----------------------------------------------------------------------*/
		/*	only attach neighbours which have a gamma > 0 			*/
		/* patches should point only to other patches in the same basin		*/
		/*	- ie. gamma to adjacent areas outside the basin 		*/
		/* 	should be zero							*/
		/* 	excepth in the case of a stream patch 				*/
		/*	 - this however is not strictly enforced			*/
		/*----------------------------------------------------------------------*/

//		printf("\t\tassign_neigh(%d): patch: %d, zone: %d, hill: %d\n", i, patch_ID, zone_ID, hill_ID);

        if (edgedistance > 0.0) { //gamma > 0.0
            
            if( (patch_ID != 0) && (zone_ID != 0) && (hill_ID != 0) ){
				neigh = find_patch(patch_ID, zone_ID, hill_ID, basin);
            }else{
                neigh = basin[0].outside_region; // ok
            }
            
			sum_gamma += gamma; // check if all add to ONE.
			neighbours[inx].gamma = gamma;
            neighbours[inx].gammaCONST = gamma;
            neighbours[inx].edgedistance = edgedistance;
            neighbours[inx].edge = edge;
			neighbours[inx].patch = neigh; // patch_object
            neighbours[inx].transmissivity_flux2neighbour = 0.0; // new variable Feb 11, 2020, Lin
			inx += 1;
            
        }else{
            new_num_neighbours -= 1;
        }// if gamma
        
	}// end of for loop


	return(new_num_neighbours);
}/*end assign_neighbours.c*/
