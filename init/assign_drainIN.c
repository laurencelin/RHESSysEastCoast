/*--------------------------------------------------------------*/
/*                                                              */ 
/*		assign_drainIN									*/
/*                                                              */
/*  NAME                                                        */
/*		assign_drainIN									*/
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/*  assign_drainIN( struct patch_object *patch)				*/
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
int assign_drainIN( struct drainIN_object *drainIN,
					   int num_drainIN,
					   struct basin_object *basin,
					   FILE *routing_file,
                       struct patch_object *currentPatch)
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
	int i, inx, new_num_drainIN;
	int patch_ID, zone_ID, hill_ID;
	double maxDailyDrain;
    double propDrainFrmSurf;
	double DrainFrac;
	struct patch_object *neigh;
	
	/*--------------------------------------------------------------*/
	/*  find and assign each neighbour to array						*/
	/*--------------------------------------------------------------*/
	inx = 0;
	maxDailyDrain = 0.0;
    propDrainFrmSurf = 0.0;
    DrainFrac = 0.0;
	new_num_drainIN = num_drainIN;
	for (i=0; i< num_drainIN; i++) {
        
		fscanf(routing_file,"%d %d %d %lf %lf %lf",
			&patch_ID,
			&zone_ID,
			&hill_ID,
			&maxDailyDrain,
            &propDrainFrmSurf,
            &DrainFrac);
		

        if (maxDailyDrain > 0.0) {
            
            if( (patch_ID != 0) && (zone_ID != 0) && (hill_ID != 0) ){
				neigh = find_patch(patch_ID, zone_ID, hill_ID, basin);
            }else{
                neigh = basin[0].outside_region; // ok
            }
            
			drainIN[inx].maxDailyDrain = maxDailyDrain;
            drainIN[inx].propDrainFrmSurf = propDrainFrmSurf;
            drainIN[inx].DrainFrac = DrainFrac;
            drainIN[inx].transfer_flux_surf = 0.0;
            drainIN[inx].transfer_flux_sub = 0.0;
			drainIN[inx].patch = neigh; // patch_object
			inx += 1;
            
        }else{
            new_num_drainIN -= 1;
        }// if gamma
        
	}// end of for loop


	return(new_num_drainIN);
}/*end assign_neighbours.c*/
