/*--------------------------------------------------------------*/
/* 																*/
/*					execute_yearly_output_event					*/
/*																*/
/*	execute_yearly_output_event - outputs yearly data			*/
/*																*/
/*	NAME														*/
/*	execute_yearly_output_event - outputs yearly data 	.		*/
/*																*/
/*	SYNOPSIS													*/
/*	void	execute_yearly_output_event(						*/
/*					struct	world_object	*world,				*/
/*					struct	command_line_object *command_line,	*/
/*					struct	date	date,  						*/
/*					struct	world_output_file_object 	*outfile)*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	outputs spatial structure according to commandline			*/
/*	specifications to specific files, for yearly info			*/
/*	a spatial structure is output only if its appropriate		*/
/*	option flag has been set and its ID matches a 				*/
/*	specified ID, or -999 which indicates all					*/
/*	units at that spatial scale are to be output				*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	We only permit one fileset per spatial modelling level.     */
/*	Each fileset has one file for each timestep.  				*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

void	execute_yearly_output_event(
									struct	world_object	*world,
									struct	command_line_object *command_line,
									struct	date	date,
									struct	world_output_file_object	*outfile)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	void output_yearly_basin(
		struct	basin_object *,
		struct	date,
		FILE	*);
	
	void output_yearly_hillslope( int,
		struct	basin_object *,
        int hID,
		struct	date,
		FILE	*);
	
	void output_yearly_zone( 	int, int,
		struct	zone_object *,
		struct	date,
		FILE	*);
	
	void output_yearly_patch(	int, int, int,
		struct	patch_object *,
		struct	date,
		FILE	*);
	
	void output_yearly_canopy_stratum(
		int, int, int, int,
		struct	canopy_strata_object *,
		struct	date,
		FILE	*);
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int	basinID, hillID, patchID, zoneID, stratumID;
	int b,h,p,z,c;
	
	
    /*--------------------------------------------------------------*/
    /*	output basins												*/
    /*--------------------------------------------------------------*/
    for (b=0; b < world[0].num_basin_files; ++ b ) {
        /*--------------------------------------------------------------*/
        /*	Construct the basin output files.							*/
        /*--------------------------------------------------------------*/
        if ( command_line[0].b != NULL ){
            basinID = command_line[0].b->basinID;
            if (( world[0].basins[b][0].ID == basinID) || (basinID == -999))
                output_yearly_basin(
                world[0].basins[b],
                date,
                outfile->basin->yearly);
        }// end of if
        /*------------------------------------------------------------*/
        /*	output hillslopes 											*/
        /*-----------------------------------------------------------*/
        for (h=0; h < world[0].basins[b][0].num_hillslopes; ++h) {
            /*--------------------------------------------------------*/
            /*	Construct the hillslope output files.						*/
            /*--------------------------------------------------------*/
            if ( command_line[0].h != NULL ){
                basinID = command_line[0].h->basinID;
                hillID = command_line[0].h->hillID;
                if (( world[0].basins[b][0].ID == basinID)
                    || (basinID == -999))
                    if (( world[0].basins[b][0].hillslopes[h][0].ID == hillID)
                        || (hillID == -999))
                        output_yearly_hillslope(
                            world[0].basins[b][0].ID,
                            world[0].basins[b],
                            world[0].basins[b]->hillslopes[h]->ID,
                            date,
                            outfile->hillslope->yearly);
            }
            /*-----------------------------------------------*/
            /*	output zones		  									*/
            /*---------------------------------------------*/
            for (z=0; z < world[0].basins[b][0].hillslopes[h][0].num_zones; ++z){
                /*--------------------------------------------*/
                /*	Construct the zone output files.					*/
                /*---------------------------------------------*/
                if ( command_line[0].z != NULL ){
                    basinID = command_line[0].z->basinID;
                    hillID = command_line[0].z->hillID;
                    zoneID = command_line[0].z->zoneID;
                    if (( world[0].basins[b][0].ID == basinID)
                        || (basinID == -999))
                        if (( world[0].basins[b][0].hillslopes[h][0].ID == hillID)
                            || (hillID == -999))
                            if (( world[0].basins[b][0].hillslopes[h][0].zones[z][0].ID == zoneID)
                                || (zoneID == -999))
                                output_yearly_zone(
                                world[0].basins[b][0].ID,
                                world[0].basins[b][0].hillslopes[h][0].ID,
                                world[0].basins[b]->hillslopes[h]->zones[z],
                                date,
                                outfile->zone->yearly);
                }
                /*---------------------------------------------*/
                /*	output patches 										*/
                /*---------------------------------------------*/
                for(p=0; p < world[0].basins[b][0].hillslopes[h][0].zones[z][0].num_patches; ++p) {
                    /*------------------------------------------*/
                    /*	Construct the patch output files.			*/
                    /*------------------------------------------*/
                    if ( command_line[0].p != NULL ){
                        basinID = command_line[0].p->basinID;
                        hillID = command_line[0].p->hillID;
                        zoneID = command_line[0].p->zoneID;
                        patchID = command_line[0].p->patchID;
                        if (( world[0].basins[b][0].ID == basinID)
                            || (basinID == -999))
                            if (( world[0].basins[b][0].hillslopes[h][0].ID == hillID)
                                || (hillID == -999))
                                if (( world[0].basins[b][0].hillslopes[h][0].zones[z][0].ID == zoneID)
                                    || (zoneID == -999))
                                    if (( world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].ID == patchID)
                                        || (patchID == -999))
                                        output_yearly_patch(
                                        world[0].basins[b]->ID,
                                        world[0].basins[b]->hillslopes[h]->ID,
                                        world[0].basins[b]->hillslopes[h]->zones[z]->ID,
                                        world[0].basins[b]->hillslopes[h]->zones[z]->patches[p],
                                        date,
                                        outfile->patch->yearly);
                    } //end of if
                    /*------------------------------------------*/
                    /*	Construct the canopy_stratum output files*/
                    /*------------------------------------------*/
                    if ( command_line[0].c != NULL ){
                        /*---------------------------------------*/
                        /*	output canopy stratum 					 */
                        /*---------------------------------------*/
                        for(c=0; c < world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].num_canopy_strata; ++c) {
                            basinID = command_line[0].c->basinID;
                            hillID = command_line[0].c->hillID;
                            zoneID = command_line[0].c->zoneID;
                            patchID = command_line[0].c->patchID;
                            stratumID = command_line[0].c->stratumID;
                            if(( world[0].basins[b][0].ID == basinID)
                                || (basinID == -999))
                                if(( world[0].basins[b][0].hillslopes[h][0].ID == hillID)
                                    || (hillID == -999))
                                    if(( world[0].basins[b][0].hillslopes[h][0].zones[z][0].ID == zoneID)
                                        || (zoneID == -999))
                                        if(( world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].ID == patchID)
                                            || (patchID == -999))
                                            if (( world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].canopy_strata[c][0].ID == stratumID)
                                                || (stratumID == -999))
                                                output_yearly_canopy_stratum(
                                                world[0].basins[b][0].ID,
                                                world[0].basins[b][0].hillslopes[h][0].ID,
                                                world[0].basins[b][0].hillslopes[h][0].zones[z][0].ID,
                                                world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].ID,
                                                world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->canopy_strata[c],
                                                date,
                                                outfile->canopy_stratum->yearly);
                        } /* end stratum (c) for loop */
                    } /* end if options */
                                
                    // reset accumulated values at each patch
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.subQnet = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.surfQnet = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.subQvnet = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.precip = 0.0;
                    //world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.recharge = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.PET = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.ET = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.sat_deficit_z = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.peakLAI = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.meanLAI = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.psn = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.days = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.satChance = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.plantlimitN = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.plantlimitQ = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.activeS = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.denitrifaspQs = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.denitrifspQs = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.nitrifaspQs = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.nitrifspQs = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.denitrif = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.mineralization = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.uptake = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.subNO3net = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.subNO3vnet = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.subDOCnet = 0.0;
                    world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->acc_year.no3drain2gw = 0.0;
                    
                    
                    
                    } /* end patch (p) for loop */
                } /* end zone (z) for  loop*/
            } /* end hillslope (h) for loop */
        } /* end basin (b) for loop */
		
		return;
} /*end execute_yearly_output_event*/
