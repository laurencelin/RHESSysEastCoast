/*--------------------------------------------------------------*/
/* 																*/
/*					execute_daily_output_event					*/
/*																*/
/*	execute_daily_output_event - outputs daily data			*/
/*																*/
/*	NAME														*/
/*	execute_daily_output_event - outputs daily data 	.		*/
/*																*/
/*	SYNOPSIS													*/
/*	void	execute_daily_output_event(						*/
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
/*	specifications to specific files, for daily info			*/
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
/*	March 14, 1997	- 	RAF				*/ 
/*	Allowed output patch to also output the moss strata if	*/
/*		moss is present.				*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

void	execute_daily_output_event(
								   struct	world_object	*world,
								   struct	command_line_object *command_line,
								   struct	date	date,
								   struct	world_output_file_object	*outfile)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	void output_24hours_basin(
		int,
		struct	basin_object *,
		struct	date,
		FILE	*);
	
	void output_hillslope(	int,
		struct	hillslope_object *,
		struct	date,
		FILE	*);
	
	void output_zone(	int, int,
		struct	zone_object *,
		struct	date,
		FILE	*);
	
	void output_patch(
		int, int, int,
		struct	patch_object *,
		struct	zone_object *,
		struct	date,
		FILE	*);
	
    void output_patch_waterstress(
          int, int, int,
          struct    patch_object *,
          struct    zone_object *,
          struct    date,
          FILE    *);
    
	void output_canopy_stratum(
		int, int, int, int,
		struct	canopy_strata_object *,
		struct	date,
		FILE	*);
	void output_shadow_strata(
		int, int, int, int,
		struct	canopy_strata_object *,
		struct	date,
		FILE	*);
        void output_stream_routing(
		struct	stream_network_object *,
		struct	date,
		FILE	*);
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int	basinID, hillID, patchID, zoneID, stratumID,reachID;
	int b,h,p,z,c,s;
	/*--------------------------------------------------------------*/
	/*	check to see if there are any print options					*/
	/*--------------------------------------------------------------*/
    
    
    /*--------------------------------------------------------------*/
    /*	output stream_routing												*/
    /*--------------------------------------------------------------*/
     for (b=0; b < world[0].num_basin_files; ++ b ) {
            for (s=0; s < world[0].basins[b][0].stream_list.num_reaches; ++s) {
        /*--------------------------------------------------------------*/
        /*	Construct the stream output files.							*/
        /*--------------------------------------------------------------*/
            if ( command_line[0].stro != NULL ){
                        reachID = command_line[0].stro->reachID;
                        if (( world[0].basins[b][0].stream_list.stream_network[s].reach_ID == reachID) || (reachID == -999))
                        {
                            output_stream_routing(
                            &(world[0].basins[b]->stream_list.stream_network[s]),
                            date,
                            outfile->stream_routing->daily);}
            }//end of if
       }//end of for s
     }//end of for b
    
    /*--------------------------------------------------------------*/
    /* spatial aggregation */
    /*--------------------------------------------------------------*/
    int i,layer;
    int tmp_ID = 0;
    int *aggregate_ID;
    float *AETa, *PETa, *AREAa;
    float *AsoilETa, *AsurfETa, *aPSN, *aroots, *aLAI, *arootPorosity;
    float area_1;
    struct patch_object* MyPatch;
    if(command_line[0].aggregate_flag>0){
        aggregate_ID = (int*)calloc(world[0].basins[0][0].aggregateLength, sizeof(int));
        // total ET, PET, soilET, PSN, root.S, LAI, TDR
        // size (num of grids)
        AETa = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        PETa = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        AREAa = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        AsoilETa = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        AsurfETa = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        aPSN = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        aroots = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        aLAI = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        arootPorosity = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        for(i = 0; i < world[0].basins[0][0].aggregateLength; i++){
            aggregate_ID[i] = 0;
            AETa[i] = 0.0;
            PETa[i] = 0.0;
            AREAa[i] = 0.0;
            
            AsoilETa[i] = 0.0;
            AsurfETa[i] = 0.0;
            aPSN[i] = 0.0;
            aroots[i] = 0.0;
            aLAI[i] = 0.0;
            arootPorosity[i] = 0.0;
        }//i
        
        for(i = 0; i < world[0].basins[0][0].route_list->num_patches; i++){

            MyPatch = world[0].basins[0][0].route_list->list[i];
            
            aggregate_ID[MyPatch->aggregate_index] = MyPatch->aggregate_ID;
            
            AREAa[MyPatch->aggregate_index] += MyPatch[0].area;
            PETa[MyPatch->aggregate_index] += MyPatch[0].PET * MyPatch[0].area;
            AETa[MyPatch->aggregate_index] += (MyPatch[0].evaporation + MyPatch[0].evaporation_surf + MyPatch[0].exfiltration_sat_zone + MyPatch[0].exfiltration_unsat_zone + MyPatch[0].transpiration_sat_zone + MyPatch[0].transpiration_unsat_zone)*MyPatch[0].area;
            
            AsoilETa[MyPatch->aggregate_index] += (MyPatch[0].exfiltration_sat_zone + MyPatch[0].exfiltration_unsat_zone)*MyPatch[0].area;
            AsurfETa[MyPatch->aggregate_index] += (MyPatch[0].evaporation + MyPatch[0].evaporation_surf)*MyPatch[0].area;
            aroots[MyPatch->aggregate_index] += MyPatch[0].rootzone.S *MyPatch[0].area;
            arootPorosity[MyPatch->aggregate_index] += MyPatch[0].rootzone.potential_sat / MyPatch[0].rootzone.depth *MyPatch[0].area;
            for ( layer=0 ; layer<MyPatch[0].num_layers; layer++ ){
                for ( c=0 ; c<MyPatch[0].layers[layer].count; c++ ){
                    aPSN[MyPatch->aggregate_index] += MyPatch[0].canopy_strata[(MyPatch[0].layers[layer].strata[c])][0].cover_fraction
                        * MyPatch[0].canopy_strata[(MyPatch[0].layers[layer].strata[c])][0].cs.net_psn
                        * MyPatch[0].area;
                    aLAI[MyPatch->aggregate_index] += MyPatch[0].canopy_strata[(MyPatch[0].layers[layer].strata[c])][0].cover_fraction
                        * MyPatch[0].canopy_strata[(MyPatch[0].layers[layer].strata[c])][0].epv.proj_lai
                        * MyPatch[0].area;
                }//c
            }//layer
            //MyPatch[0].canopy_strata[0][0].ID != 4;  not correct. it refers to StratumID
            //MyPatch[0].canopy_strata[0][0].defaults[0][0].epc.veg_type != NON_VEG
            //if(MyPatch[0].canopy_strata[0][0].defaults[0][0].epc.veg_type != NON_VEG){}// vegetation only
     
           
        }//i for
        
        // <----- here print out
        for(i = 0; i < world[0].basins[0][0].aggregateLength; i++){
            area_1 = 1.0/AREAa[i];
            fprintf(outfile->aggregate->daily,"%d,%d,%d,%d, %f,%f,%f,%f,%f,%f,%f,%f,%f\n",date.day,date.month,date.year,
                    aggregate_ID[i],
                    AETa[i]*1000.0*area_1,
                    PETa[i]*1000.0*area_1,
                    AREAa[i],
                    AsoilETa[i]*1000.0*area_1,
                    AsurfETa[i]*1000.0*area_1,
                    aPSN[i]*area_1,
                    aroots[i]*area_1,
                    aLAI[i]*area_1,
                    arootPorosity[i]*area_1
                );
        }//i
        
    }//spatial re-aggregate
    
    
    
	if ((command_line[0].b != NULL) || (command_line[0].h != NULL) ||
		(command_line[0].z != NULL) || (command_line[0].p != NULL) ||
		(command_line[0].c != NULL)){
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
					output_24hours_basin(
					command_line[0].routing_flag,
					world[0].basins[b],
					date,
					outfile->basin->daily);
			}
			/*--------------------------------------------------------------*/
			/*	check to see if there are any lower print options			*/
			/*--------------------------------------------------------------*/
			if ((command_line[0].h != NULL) || (command_line[0].z != NULL) ||
				(command_line[0].p != NULL) || (command_line[0].c != NULL)){
				/*--------------------------------------------------------------*/
				/*	output hillslopes 											*/
				/*--------------------------------------------------------------*/
				for (h=0; h < world[0].basins[b][0].num_hillslopes; ++h) {
					/*-----------------------------------------------------------*/
					/*	Construct the hillslope output files.						*/
					/*-----------------------------------------------------------*/
					if ( command_line[0].h != NULL ){
						basinID = command_line[0].h->basinID;
						hillID = command_line[0].h->hillID;
						if (( world[0].basins[b][0].ID == basinID)
							|| (basinID == -999))
							if (( world[0].basins[b][0].hillslopes[h][0].ID == hillID)
								|| (hillID == -999))
								output_hillslope(
								world[0].basins[b][0].ID,
								world[0].basins[b]->hillslopes[h],
								date,
								outfile->hillslope->daily);
					}
					/*-------------------------------------------------------------*/
					/*	check to see if there are any lower print options			*/
					/*-------------------------------------------------------------*/
					if ((command_line[0].z != NULL) || (command_line[0].p != NULL)
						|| (command_line[0].c != NULL)){
						/*---------------------------------------------------------*/
						/*	output zones												*/
						/*---------------------------------------------------------*/
						for(z=0;
						z < world[0].basins[b][0].hillslopes[h][0].num_zones;
						++z){
							/*------------------------------------------------------*/
							/*	Construct the zone output files.						  */
							/*-------------------------------------------------------*/
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
											output_zone(
											world[0].basins[b][0].ID,
											world[0].basins[b][0].hillslopes[h][0].ID,
											world[0].basins[b]->hillslopes[h]->zones[z],
											date, outfile->zone->daily);
							}
							/*-------------------------------------------------------*/
							/*	check to see if there are any lower print options		*/
							/*-------------------------------------------------------*/
							if ((command_line[0].p != NULL)
								|| (command_line[0].c != NULL)){
								/*----------------------------------------------------*/
								/*	output patches 												*/
								/*---------------------------------------------------*/
								for(p=0;
								p < world[0].basins[b][0].hillslopes[h][0].zones[z][0].num_patches;
								++p){									
									/*-------------------------------------------------*/
									/*	Construct the patch output files.					*/
									/*-------------------------------------------------*/
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
														|| (patchID == -999 && world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->landuse_defaults[0]->ID >= command_line[0].patchPrintTh)){
                                                            output_patch(
                                                                     world[0].basins[b]->ID,
                                                                     world[0].basins[b]->hillslopes[h]->ID,
                                                                     world[0].basins[b]->hillslopes[h]->zones[z]->ID,
                                                                     world[0].basins[b]->hillslopes[h]->zones[z]->patches[p],
                                                                     world[0].basins[b]->hillslopes[h]->zones[z],
                                                                     date,
                                                                     outfile->patch->daily);
                                                        
                                                        }//if
									}//if
									/*------------------------------------------------*/
									/*	Construct the canopy_stratum output files		  */
									/*------------------------------------------------*/
									if ( command_line[0].c != NULL ){
										/*----------------------------------------------*/
										/*	output canopy stratum 								*/
										/*----------------------------------------------*/
										for(c=0;
										c < world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].num_canopy_strata;
										++c){
											basinID = command_line[0].c->basinID;
											hillID = command_line[0].c->hillID;
											zoneID = command_line[0].c->zoneID;
											patchID = command_line[0].c->patchID;
											stratumID = command_line[0].c->stratumID;
											if (( world[0].basins[b][0].ID == basinID)
												|| (basinID == -999))
												if (( world[0].basins[b][0].hillslopes[h][0].ID == hillID)
													|| (hillID == -999))
													if (( world[0].basins[b][0].hillslopes[h][0].zones[z][0].ID == zoneID)
														|| (zoneID == -999))
														if (( world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].ID == patchID)
															||	(patchID == -999 && world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->landuse_defaults[0]->ID >= command_line[0].patchPrintTh))
															if (( world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].canopy_strata[c][0].ID == stratumID)
																|| (stratumID == -999)) {
																output_canopy_stratum(
																world[0].basins[b][0].ID,
																world[0].basins[b][0].hillslopes[h][0].ID,
																world[0].basins[b][0].hillslopes[h][0].zones[z][0].ID,
																world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].ID,
																world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->canopy_strata[c],
																date, outfile->canopy_stratum->daily);
															}
										} /* end stratum (c) for loop */
									} /* end if options */
								} /* end patch (p) for loop */
							} /* end if options */
						} /* end zone (z) for  loop*/
					} /* end if options */
					} /* end hillslope (h) for loop */
				} /* end if options */
			} /* end basin (b) for loop */
		} /* end if options */
		return;
} /*end execute_daily_output_event*/
