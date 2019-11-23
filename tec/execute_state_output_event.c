/*--------------------------------------------------------------*/
/* 																*/
/*					execute_state_output_event					*/
/*																*/
/*	execute_state_output_event - outputs state data			*/
/*																*/
/*	NAME														*/
/*	execute_state_output_event - outputs state data 	.		*/
/*																*/
/*	SYNOPSIS													*/
/*	void	execute_state_output_event(						*/
/*					struct	world_object	*world,				*/
/*					struct	date	current_date,				*/
/*					struct	date	end_date,				*/
/*					struct	command_line_object *command_line)	*/
/*																*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	outputs current world state - in worldfile format			*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"

void	execute_state_output_event(
								   struct	world_object	*world,
								   struct	date	current_date,
								   struct	date	end_date,
								   struct	command_line_object *command_line)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	void output_basin_state(
		struct	basin_object *,
		struct	date,
		struct	command_line_object *,
		FILE	*);
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int b,i;
	FILE	*outfile;
	char	filename[MAXSTR+100];
	char	ext[20];
	/*--------------------------------------------------------------*/
	/*	Try to open the world file in read mode.					*/
	/*--------------------------------------------------------------*/
	sprintf(ext,".Y%4dM%dD%dH%d",
            current_date.year,
            current_date.month,
            current_date.day,
            current_date.hour);

	strcpy(filename, command_line[0].world_filename);
	strcat(filename, ext);
	strcat(filename, ".csv");

	/*--------------------------------------------------------------*/
	/*	open output file											*/
	/*--------------------------------------------------------------*/
	if ( ( outfile = fopen(filename, "w")) == NULL ){
		fprintf(stderr,"FATAL ERROR: in execute_state_output_event.\n");
		exit(EXIT_FAILURE);
	}
    
    //fprintf(outfile, "%d %s\n", world[0].ID, "world_id");
    //fprintf(outfile, "%d %s", world[0].num_basin_files, "NUM_of_");
    
    //full title
    struct basin_object *basin;
    struct hillslope_object *hillslope;
    struct zone_object *zone;
    struct patch_object *patch;
    struct canopy_strata_object *stratum;
    char world_line [50];
    char basin_line [1000];
    char hill_line [2000];
    char zone_line [2000];
    char patch_line [3000];
    char stratum_line [4000];
    
    fprintf(outfile, "world_ID,basin_ID,x,y,z,basin_parm_ID,latitude,hillslope_ID,x,y,z,hill_parm_ID,gw.storage,gw.NO3,gw.NH4,gw.DOC,gw.DON,zone_ID,x,y,z,zone_parm_ID,area,slope,aspect,precip_lapse_rate,e_horizon,w_horizon,base_station_ID,patch_ID,x,y,z,soil_parm_ID,landuse_parm_ID,area,slope,lna,Ksat_vertical,detention_store,surface_DOC,surface_DON,surface_NO3,surface_NH4,rz_storage,unsat_storage,sat_deficit,satzZ_balance,snowpack.water_equivalent_depth,snowpack.water_depth,snowpack.T,snowpack.surface_age,snowpack.energy_deficit,litter.cover_fraction,litter.rain_stored,litter.NO3_stored,litter_cs.litr1c,litter_cs.litr2c,litter_cs.litr3c,litter_cs.litr4c,litter_ns.litr1n,litter_ns.litr2n,litter_ns.litr3n,litter_ns.litr4n,soil_cs.soil1c,soil_cs.soil2c,soil_cs.soil3c,soil_cs.soil4c,soil_ns.soil1n,soil_ns.soil2n,soil_ns.soil3n,soil_ns.soil4n,patch_SOIL1_CN,patch_SOIL2_CN,patch_SOIL3_CN,patch_SOIL4_CN,soil_ns.sminn,soil_ns.nitrate,soil_cs.DOC,soil_ns.DON,sat_DOC,sat_DON,sat_NO3,sat_NH4,stored_fertilizer_NO3,stored_fertilizer_NH4,fertilizerDaysCount,canopy_strata_ID,veg_parm_ID,cover_fraction,gap_fraction,rootzone.depth,cs.age,snow_stored,rain_stored,NO3_stored,cs.cpool,cs.leafc,cs.dead_leafc,cs.leafc_store,cs.leafc_transfer,cs.live_stemc,cs.livestemc_store,cs.livestemc_transfer,cs.dead_stemc,cs.deadstemc_store,cs.deadstemc_transfer,cs.live_crootc,cs.livecrootc_store,cs.livecrootc_transfer,cs.dead_crootc,cs.deadcrootc_store,cs.deadcrootc_transfer,cs.frootc,cs.frootc_store,cs.frootc_transfer,cs.cwdc,epv.prev_leafcalloc,ns.npool,ns.leafn,ns.dead_leafn,ns.leafn_store,ns.leafn_transfer,ns.live_stemn,ns.livestemn_store,ns.livestemn_transfer,ns.dead_stemn,ns.deadstemn_store,ns.deadstemn_transfer,ns.live_crootn,ns.livecrootn_store,ns.livecrootn_transfer,ns.dead_crootn,ns.deadcrootn_store,ns.deadcrootn_transfer,ns.frootn,ns.frootn_store,ns.frootn_transfer,ns.cwdn,ns.retransn,ns.cwdN_stored,epv.wstress_days,epv.max_fparabs,epv.min_vwc\n");
    

    sprintf(world_line, "%d",world[0].ID);
	for (int b=0; b < world[0].num_basin_files; ++ b ) {
//		output_basin_state(
//            world[0].basins[b],
//            current_date,
//            command_line,
//            outfile,
//            title,
//            line);
        basin = world[0].basins[b];
        sprintf(basin_line, "%d,%e,%e,%e,%d,%e",
                basin[0].ID,
                basin[0].x,
                basin[0].y,
                basin[0].z,
                basin[0].defaults[0][0].ID,
                basin[0].latitude);
        
        for (int h=0; h < basin[0].num_hillslopes; ++ h ) {
            hillslope = basin[0].hillslopes[h];
            sprintf(hill_line,"%d,%e,%e,%e,%d,%e,%e,%e,%e,%e",
                    hillslope[0].ID,
                    hillslope[0].x,
                    hillslope[0].y,
                    hillslope[0].z,
                    hillslope[0].defaults[0][0].ID,
                    hillslope[0].gw.storage,
                    hillslope[0].gw.NO3,
                    hillslope[0].gw.NH4,
                    hillslope[0].gw.DOC,
                    hillslope[0].gw.DON);
             for (int z=0; z < hillslope[0].num_zones; ++ z ) {
                 zone = hillslope[0].zones[z];
                 sprintf(zone_line,"%d,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%d",
                         zone[0].ID,
                         zone[0].x,
                         zone[0].y,
                         zone[0].z,
                         zone[0].defaults[0][0].ID,
                         zone[0].area,
                         zone[0].slope,
                         zone[0].aspect,
                         zone[0].precip_lapse_rate,
                         zone[0].e_horizon,
                         zone[0].w_horizon,
                         zone[0].base_stations[0][0].ID); // one station per zone
                 for (int p=0; p < zone[0].num_patches; ++ p ) {
                     patch = zone[0].patches[p];
                     sprintf(patch_line,"%d,%e,%e,%e,%d,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d",
                             patch[0].ID,
                             patch[0].x,
                             patch[0].y,
                             patch[0].z,
                             patch[0].soil_defaults[0][0].ID,
                             patch[0].landuse_defaults[0][0].ID,
                             patch[0].area,//1
                             patch[0].slope,//2
                             patch[0].lna,//3
                             patch[0].Ksat_vertical,//4
                             patch[0].detention_store,//41
                             patch[0].surface_DOC,//42
                             patch[0].surface_DON,//43
                             patch[0].surface_NO3,//44
                             patch[0].surface_NH4,//45
                             patch[0].rz_storage,//5
                             patch[0].unsat_storage,//6
                             patch[0].sat_deficit,//7
                             patch[0].satzZ_balance,//48
                             patch[0].snowpack.water_equivalent_depth,//8
                             patch[0].snowpack.water_depth,//9
                             patch[0].snowpack.T,//10
                             patch[0].snowpack.surface_age,//11
                             patch[0].snowpack.energy_deficit,//12
                             patch[0].litter.cover_fraction,//13
                             patch[0].litter.rain_stored,//14
                             patch[0].litter.NO3_stored,//46
                             patch[0].litter_cs.litr1c,//15
                             patch[0].litter_cs.litr2c,//16
                             patch[0].litter_cs.litr3c,//17
                             patch[0].litter_cs.litr4c,//18
                             patch[0].litter_ns.litr1n,//19
                             patch[0].litter_ns.litr2n,//20
                             patch[0].litter_ns.litr3n,//21
                             patch[0].litter_ns.litr4n,//22
                             patch[0].soil_cs.soil1c,//23
                             patch[0].soil_cs.soil2c,//24
                             patch[0].soil_cs.soil3c,//25
                             patch[0].soil_cs.soil4c,//26
                             patch[0].soil_ns.soil1n,//27
                             patch[0].soil_ns.soil2n,//28
                             patch[0].soil_ns.soil3n,//29
                             patch[0].soil_ns.soil4n,//30
                             patch[0].patch_SOIL1_CN,//49
                             patch[0].patch_SOIL2_CN,//50
                             patch[0].patch_SOIL3_CN,//51
                             patch[0].patch_SOIL4_CN,//52
                             patch[0].soil_ns.sminn,//31
                             patch[0].soil_ns.nitrate,//32
                             patch[0].soil_cs.DOC,//33
                             patch[0].soil_ns.DON,//34
                             patch[0].sat_DOC,//35
                             patch[0].sat_DON,//36
                             patch[0].sat_NO3,//37
                             patch[0].sat_NH4,//38
                             patch[0].stored_fertilizer_NO3,//39
                             patch[0].stored_fertilizer_NH4, //40
                             patch[0].fertilizerDaysCount);//47
                     for (int s=0; s < patch[0].num_canopy_strata; ++ s ) {
                         stratum = patch[0].canopy_strata[s];
//                         printf("state output stratum (%d,%d,%d,%f).Y%4dM%dD%dH%d: %f or %lf or %e\n",
//                                patch[0].ID, stratum[0].ID, stratum[0].defaults[0][0].ID, stratum[0].cover_fraction,
//                                current_date.year,
//                                current_date.month,
//                                current_date.day,
//                                current_date.hour,
//                                stratum[0].rain_stored, stratum[0].rain_stored, stratum[0].rain_stored);
                        // somehow the string below is not correct
                         sprintf(stratum_line,"%d,%d,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%e",
                                 stratum[0].ID,//1
                                 stratum[0].defaults[0][0].ID,//2
                                 stratum[0].cover_fraction, //1
                                 stratum[0].gap_fraction, //2
                                 stratum[0].rootzone.depth, //3
                                 stratum[0].cs.age, //4
                                 stratum[0].snow_stored, //5
                                 stratum[0].rain_stored, //6
                                 stratum[0].NO3_stored, // from deposition
                                 stratum[0].cs.cpool, //7
                                 stratum[0].cs.leafc, //8
                                 stratum[0].cs.dead_leafc, //9
                                 stratum[0].cs.leafc_store, //10
                                 stratum[0].cs.leafc_transfer, //11
                                 stratum[0].cs.live_stemc, //12
                                 stratum[0].cs.livestemc_store, //13
                                 stratum[0].cs.livestemc_transfer, //14
                                 stratum[0].cs.dead_stemc, //15
                                 stratum[0].cs.deadstemc_store, //16
                                 stratum[0].cs.deadstemc_transfer, //17
                                 stratum[0].cs.live_crootc, //18
                                 stratum[0].cs.livecrootc_store, //19
                                 stratum[0].cs.livecrootc_transfer, //20
                                 stratum[0].cs.dead_crootc, //21
                                 stratum[0].cs.deadcrootc_store, //22
                                 stratum[0].cs.deadcrootc_transfer, //23
                                 stratum[0].cs.frootc, //24
                                 stratum[0].cs.frootc_store, //25
                                 stratum[0].cs.frootc_transfer, //26
                                 stratum[0].cs.cwdc, //27
                                 stratum[0].epv.prev_leafcalloc, //28
                                 stratum[0].ns.npool, //29
                                 stratum[0].ns.leafn, //30
                                 stratum[0].ns.dead_leafn, //31
                                 stratum[0].ns.leafn_store, //32
                                 stratum[0].ns.leafn_transfer, //33
                                 stratum[0].ns.live_stemn, //34
                                 stratum[0].ns.livestemn_store, //35
                                 stratum[0].ns.livestemn_transfer, //36
                                 stratum[0].ns.dead_stemn, //37
                                 stratum[0].ns.deadstemn_store, //38
                                 stratum[0].ns.deadstemn_transfer, //39
                                 stratum[0].ns.live_crootn, //40
                                 stratum[0].ns.livecrootn_store, //41
                                 stratum[0].ns.livecrootn_transfer, //42
                                 stratum[0].ns.dead_crootn, //43
                                 stratum[0].ns.deadcrootn_store, //44
                                 stratum[0].ns.deadcrootn_transfer, //45
                                 stratum[0].ns.frootn, //46
                                 stratum[0].ns.frootn_store, //47
                                 stratum[0].ns.frootn_transfer, //48
                                 stratum[0].ns.cwdn, //49
                                 stratum[0].ns.retransn, //50
                                 stratum[0].ns.cwdN_stored, //51
                                 stratum[0].epv.wstress_days, //52
                                 stratum[0].epv.max_fparabs, //53
                                 stratum[0].epv.min_vwc); //54
                         
                         fprintf(outfile,"%s,%s,%s,%s,%s,%s\n",
                                 world_line,
                                 basin_line,
                                 hill_line,
                                 zone_line,
                                 patch_line,
                                 stratum_line);
                         
                     }//stratum
                 }//patch
             }//zone
        }// hill
	}// for basin
    
	fclose(outfile);
	return;
} /*end execute_state_output_event*/
