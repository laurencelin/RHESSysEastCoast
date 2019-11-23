/*--------------------------------------------------------------*/
/* 																*/
/*					construct_canopy_strata						*/
/*																*/
/*	construct_canopy_strata.c - creates a canopy_strata object	*/
/*																*/
/*	NAME														*/
/*	construct_canopy_strata.c - creates a canopy_strata object	*/
/*																*/
/*	SYNOPSIS													*/
/*	struct canopy_strata_object  construct_canopy_strata( 		*/
/*				struct	command_line_object	*command_line,		*/
/*				FILE	*world_file)						*/
/*																*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	Allocates memory for a canopy strata object and reads in 	*/
/*	parameters for the object from an open hillslopes file.		*/
/*																*/
/*	Refer to construct_basin.c for a specification of the 		*/
/*	hillslopes file.											*/
/*																*/
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
/*	Sep 2 1997	RAF												*/
/*	Removed all references to grow or extended state variables.	*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"
#include "params.h"

struct canopy_strata_object *construct_canopy_strata(
													 struct	command_line_object	*command_line,
													 FILE	*world_file,
													 struct	patch_object	*patch,
													 int		num_world_base_stations,
													 struct base_station_object **world_base_stations,
													 struct	default_object	*defaults)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	struct base_station_object *assign_base_station(
		int ,
		int ,
		struct base_station_object **);
	
	
	int compute_annual_turnover(struct epconst_struct,
		struct epvar_struct *,
		struct cstate_struct *);


	int	compute_annual_litfall(
		struct epconst_struct,
		struct phenology_struct *,
		struct cstate_struct *,
		int);
	
	int	update_rooting_depth(
		struct rooting_zone_object *,
		double,
		double,
		double,
		double);

	void	*alloc(size_t, char *, char *);
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int	base_stationID;
	int	i;
	double	sai, rootc;
	int	default_object_ID; 
	int     spinup_default_object_ID; 
	char	record[MAXSTR];
	struct	canopy_strata_object	*canopy_strata;
	int    paramCnt=0;
    param    *paramPtr=NULL;
    
	/*--------------------------------------------------------------*/
	/*  Allocate a canopy_strata object.                                */
	/*--------------------------------------------------------------*/
	canopy_strata = (struct canopy_strata_object *) alloc( 1 *
		sizeof( struct canopy_strata_object ),"canopy_strata",
		"construct_canopy_strata" );
	
	/*--------------------------------------------------------------*/
	/*	Read in the next canopy strata record for this patch.	*/
	/*--------------------------------------------------------------*/
    paramPtr = readtag_worldfile(&paramCnt,world_file,"Canopy_Strata");
    canopy_strata[0].ID = getIntWorldfile(&paramCnt, &paramPtr, "canopy_strata_ID","%d",-9999,0);//1
    default_object_ID = getIntWorldfile(&paramCnt, &paramPtr, "veg_parm_ID","%d",-9999,0);//2
    
    if (command_line[0].vegspinup_flag > 0){
        spinup_default_object_ID = getIntWorldfile(&paramCnt,&paramPtr,"spinup_object_ID","%d",-9999,0);//3
    }//if
    
    canopy_strata[0].cover_fraction = getDoubleWorldfile(&paramCnt,&paramPtr,"cover_fraction","%lf",1.0,1);//4
    canopy_strata[0].gap_fraction = getDoubleWorldfile(&paramCnt,&paramPtr,"gap_fraction","%lf",0.0,1);//5
    canopy_strata[0].rootzone.depth = getDoubleWorldfile(&paramCnt,&paramPtr,"rootzone.depth","%lf",-9999,0);//6
    canopy_strata[0].rootzone.depth *= command_line[0].rootdepthz;
    
//	fscanf(world_file,"%d",&(canopy_strata[0].ID));
//	read_record(world_file, record);
//	fscanf(world_file,"%d",&(default_object_ID));
//	read_record(world_file, record);
//	if (command_line[0].vegspinup_flag > 0){
//	    fscanf(world_file,"%d",&(spinup_default_object_ID));
//	    read_record(world_file, record);
//        }
//  fscanf(world_file,"%lf",&(canopy_strata[0].cover_fraction));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].gap_fraction));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].rootzone.depth));
//	read_record(world_file, record);
//		if (command_line[0].tmp_value > ZERO)
//			canopy_strata[0].rootzone.depth *= command_line[0].tmp_value;
    
    // something new
//    canopy_strata[0].cs.stem_density = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.stem_density","%lf",-9999,1);
//    canopy_strata[0].cs.leafc_age2 = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.leafc_age2","%lf",canopy_strata[0].cs.leafc/2.0,1);
//    canopy_strata[0].cs.leafc_age1 = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.leafc_age1","%lf",canopy_strata[0].cs.leafc/2.0,1);
    
    canopy_strata[0].cs.age = getIntWorldfile(&paramCnt,&paramPtr,"cs.age","%d",0,1);//7
    canopy_strata[0].snow_stored = getDoubleWorldfile(&paramCnt,&paramPtr,"snow_stored","%lf",0.0,1);//8
    canopy_strata[0].rain_stored = getDoubleWorldfile(&paramCnt,&paramPtr,"rain_stored","%lf",0.0,1);//9
    canopy_strata[0].NO3_stored = getDoubleWorldfile(&paramCnt,&paramPtr,"NO3_stored","%lf",0.0,1);//9
    
    canopy_strata[0].cs.cpool = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.cpool","%lf",0.0,1);//10
    canopy_strata[0].cs.leafc = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.leafc","%lf",0.0,1);//11
    canopy_strata[0].cs.dead_leafc = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.dead_leafc","%lf",0.0,1);//12
    canopy_strata[0].cs.leafc_store = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.leafc_store","%lf",0.0,1);//13
    canopy_strata[0].cs.leafc_transfer = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.leafc_transfer","%lf",0.0,1);//14
    canopy_strata[0].cs.live_stemc = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.live_stemc","%lf",0.0,1);//15
    canopy_strata[0].cs.livestemc_store = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.livestemc_store","%lf",0.0,1);//16
    canopy_strata[0].cs.livestemc_transfer = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.livestemc_transfer","%lf",0.0,1);//17
    canopy_strata[0].cs.dead_stemc = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.dead_stemc","%lf",0.0,1);//18
    canopy_strata[0].cs.deadstemc_store = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.deadstemc_store","%lf",0.0,1);//19
    canopy_strata[0].cs.deadstemc_transfer = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.deadstemc_transfer","%lf",0.0,1);//20
    canopy_strata[0].cs.live_crootc = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.live_crootc","%lf",0.0,1);//21
    canopy_strata[0].cs.livecrootc_store = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.livecrootc_store","%lf",0.0,1);//22
    canopy_strata[0].cs.livecrootc_transfer = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.livecrootc_transfer","%lf",0.0,1);//23
    canopy_strata[0].cs.dead_crootc = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.dead_crootc","%lf",0.0,1);//24
    canopy_strata[0].cs.deadcrootc_store = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.deadcrootc_store","%lf",0.0,1);//25
    canopy_strata[0].cs.deadcrootc_transfer = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.deadcrootc_transfer","%lf",0.0,1);//26
    canopy_strata[0].cs.frootc = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.frootc","%lf",0.0,1);//27
    canopy_strata[0].cs.frootc_store = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.frootc_store","%lf",0.0,1);//28
    canopy_strata[0].cs.frootc_transfer = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.frootc_transfer","%lf",0.0,1);//29
    canopy_strata[0].cs.cwdc = getDoubleWorldfile(&paramCnt,&paramPtr,"cs.cwdc","%lf",0.0,1);//30
    canopy_strata[0].epv.prev_leafcalloc = getDoubleWorldfile(&paramCnt,&paramPtr,"epv.prev_leafcalloc","%lf",0.0,1);//31
    
    canopy_strata[0].ns.npool = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.npool","%lf",0.0,1);//32
    canopy_strata[0].ns.leafn = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.leafn","%lf",0.0,1);//33
    canopy_strata[0].ns.dead_leafn = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.dead_leafn","%lf",0.0,1);//34
    canopy_strata[0].ns.leafn_store = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.leafn_store","%lf",0.0,1);//35
    canopy_strata[0].ns.leafn_transfer = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.leafn_transfer","%lf",0.0,1);//36
    canopy_strata[0].ns.live_stemn = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.live_stemn","%lf",0.0,1);//37
    canopy_strata[0].ns.livestemn_store = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.livestemn_store","%lf",0.0,1);//38
    canopy_strata[0].ns.livestemn_transfer = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.livestemn_transfer","%lf",0.0,1);//39
    canopy_strata[0].ns.dead_stemn = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.dead_stemn","%lf",0.0,1);//40
    canopy_strata[0].ns.deadstemn_store = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.deadstemn_store","%lf",0.0,1);//41
    canopy_strata[0].ns.deadstemn_transfer = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.deadstemn_transfer","%lf",0.0,1);//42
    canopy_strata[0].ns.live_crootn = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.live_crootn","%lf",0.0,1);//43
    canopy_strata[0].ns.livecrootn_store = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.livecrootn_store","%lf",0.0,1);//44
    canopy_strata[0].ns.livecrootn_transfer = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.livecrootn_transfer","%lf",0.0,1);//45
    canopy_strata[0].ns.dead_crootn = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.dead_crootn","%lf",0.0,1);//46
    canopy_strata[0].ns.deadcrootn_store = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.deadcrootn_store","%lf",0.0,1);//47
    canopy_strata[0].ns.deadcrootn_transfer = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.deadcrootn_transfer","%lf",0.0,1);//48
    canopy_strata[0].ns.frootn = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.frootn","%lf",0.0,1);//49
    canopy_strata[0].ns.frootn_store = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.frootn_store","%lf",0.0,1);//50
    canopy_strata[0].ns.frootn_transfer = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.frootn_transfer","%lf",0.0,1);//51
    canopy_strata[0].ns.cwdn = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.cwdn","%lf",0.0,1);//52
    canopy_strata[0].ns.retransn = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.retransn","%lf",0.0,1);//53
    canopy_strata[0].ns.cwdN_stored = getDoubleWorldfile(&paramCnt,&paramPtr,"ns.cwdN_stored","%lf",0.0,1);//58
    
//	fscanf(world_file,"%lf",&(canopy_strata[0].snow_stored));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].rain_stored));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.cpool));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.leafc));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.dead_leafc));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.leafc_store));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.leafc_transfer));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.live_stemc));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.livestemc_store));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.livestemc_transfer));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.dead_stemc));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.deadstemc_store));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.deadstemc_transfer));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.live_crootc));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.livecrootc_store));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.livecrootc_transfer));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.dead_crootc));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.deadcrootc_store));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.deadcrootc_transfer));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.frootc));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.frootc_store));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.frootc_transfer));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].cs.cwdc));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].epv.prev_leafcalloc));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.npool));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.leafn));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.dead_leafn));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.leafn_store));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.leafn_transfer));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.live_stemn));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.livestemn_store));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.livestemn_transfer));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.dead_stemn));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.deadstemn_store));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.deadstemn_transfer));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.live_crootn));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.livecrootn_store));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.livecrootn_transfer));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.dead_crootn));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.deadcrootn_store));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.deadcrootn_transfer));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.frootn));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.frootn_store));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.frootn_transfer));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.cwdn));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].ns.retransn));
//	read_record(world_file, record);


	if (command_line[0].vegspinup_flag > 0){
        canopy_strata[0].target.lai = NULLVAL;
        canopy_strata[0].target.total_stemc = NULLVAL;
        canopy_strata[0].target.met = 2;
        //canopy_strata[0].target.height = NULLVAL;
        //canopy_strata[0].target.age = NULLVAL;
   }//if
    
	/*--------------------------------------------------------------*/
	/*	intialized annual flux variables			*/
	/*--------------------------------------------------------------*/
//	fscanf(world_file,"%d",&(canopy_strata[0].epv.wstress_days));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].epv.max_fparabs));
//	read_record(world_file, record);
//	fscanf(world_file,"%lf",&(canopy_strata[0].epv.min_vwc));
//	read_record(world_file, record);
	
    canopy_strata[0].epv.wstress_days  = getIntWorldfile(&paramCnt,&paramPtr,"epv.wstress_days","%d",0,1);//54
    canopy_strata[0].epv.max_fparabs = getDoubleWorldfile(&paramCnt,&paramPtr,"epv.max_fparabs","%lf",0.0,1);//55
    canopy_strata[0].epv.min_vwc = getDoubleWorldfile(&paramCnt,&paramPtr,"epv.min_vwc","%lf",0.0,1);//56

	/*--------------------------------------------------------------*/
	/*	Assign	defaults for this canopy_strata								*/
	/*--------------------------------------------------------------*/
	canopy_strata[0].defaults = (struct stratum_default **)
		alloc( sizeof(struct stratum_default *),"defaults",
		"construct_canopy_strata" );
	i = 0;
	while (defaults[0].stratum[i].ID != default_object_ID) {
		i++;
		/*--------------------------------------------------------------*/
		/*  Report an error if no match was found.  Otherwise assign    */
		/*  the default to point to this canopy_strata.						    */
		/*--------------------------------------------------------------*/
		if ( i>= defaults[0].num_stratum_default_files ){
			fprintf(stderr,
				"\nFATAL ERROR: in construct_canopy_strata, canopy_strata default ID %d not found.\n" ,
				default_object_ID);
			exit(EXIT_FAILURE);
		}
	} /* end-while */
	canopy_strata[0].defaults[0] = &defaults[0].stratum[i];
    
    canopy_strata[0].local_max_lai = canopy_strata[0].defaults[0]->epc.max_lai;
//    if(canopy_strata[0].defaults[0]->ID == 802 || canopy_strata[0].defaults[0]->epc.veg_type == GRASS || canopy_strata[0].defaults[0]->epc.veg_type == C4GRASS){
//        // understory
//        //canopy_strata[0].local_max_lai = min( canopy_strata[0].defaults[0]->epc.max_lai, (canopy_strata[0].cs.leafc+canopy_strata[0].cs.leafc_transfer+canopy_strata[0].cs.leafc_store)* canopy_strata[0].defaults[0]->epc.proj_sla * 1.1);
//        canopy_strata[0].local_max_lai = canopy_strata[0].defaults[0]->epc.max_lai;
//    }else{
//        //canopy_strata[0].local_max_lai = min( canopy_strata[0].defaults[0]->epc.max_lai, (canopy_strata[0].cs.leafc+canopy_strata[0].cs.leafc_transfer+canopy_strata[0].cs.leafc_store)* canopy_strata[0].defaults[0]->epc.proj_sla * 1.5);
//        canopy_strata[0].local_max_lai = canopy_strata[0].defaults[0]->epc.max_lai;
//    }
    
//    printf("[%d, %d] {%e, %e, %e}\n", patch[0].ID, canopy_strata[0].defaults[0]->ID,
//           canopy_strata[0].local_max_lai, canopy_strata[0].defaults[0]->epc.proj_sla,
//           canopy_strata[0].cs.leafc+canopy_strata[0].cs.leafc_transfer+canopy_strata[0].cs.leafc_store);

	/*--------------------------------------------------------------*/
	/* if spinup module is called assign spinup defaults            */
	/*--------------------------------------------------------------*/
	if (command_line[0].vegspinup_flag > 0) {
        canopy_strata[0].spinup_defaults = (struct spinup_default **)
            alloc( sizeof(struct spinup_default *),"defaults",
            "construct_stratum" );
        i = 0;
        while (defaults[0].spinup[i].ID != spinup_default_object_ID) {
            i++;
            /*--------------------------------------------------------------*/
            /*  Report an error if no match was found.  Otherwise assign    */
            /*  the default to point to this patch.						    */
            /*--------------------------------------------------------------*/
            if ( i>= defaults[0].num_spinup_default_files ){
                fprintf(stderr,
                    "\nFATAL ERROR: in construct_stratum, spinup default ID %d not found for patch %d\n" ,
                    spinup_default_object_ID, patch[0].ID);
                exit(EXIT_FAILURE);
            }
        } /* end-while */
        canopy_strata[0].spinup_defaults[0] = &defaults[0].spinup[i];
	}//if

	/*--------------------------------------------------------------*/
	/* zero all non tree stem and wood variables			*/
	/*--------------------------------------------------------------*/
	if (canopy_strata[0].defaults[0][0].epc.veg_type != TREE) {
		canopy_strata[0].cs.live_stemc = 0.0;
		canopy_strata[0].cs.dead_stemc = 0.0;
		canopy_strata[0].cs.live_crootc = 0.0;
		canopy_strata[0].cs.dead_crootc = 0.0;
		canopy_strata[0].cs.livestemc_store = 0.0;
		canopy_strata[0].cs.deadstemc_store = 0.0;
		canopy_strata[0].cs.livestemc_transfer = 0.0;
		canopy_strata[0].cs.deadstemc_transfer = 0.0;
		canopy_strata[0].cs.livecrootc_store = 0.0;
		canopy_strata[0].cs.deadcrootc_store = 0.0;
		canopy_strata[0].cs.livecrootc_transfer = 0.0;
		canopy_strata[0].cs.deadcrootc_transfer = 0.0;
		canopy_strata[0].cs.cwdc = 0.0;

		canopy_strata[0].ns.live_stemn = 0.0;
		canopy_strata[0].ns.dead_stemn = 0.0;
		canopy_strata[0].ns.live_crootn = 0.0;
		canopy_strata[0].ns.dead_crootn = 0.0;
		canopy_strata[0].ns.livestemn_store = 0.0;
		canopy_strata[0].ns.deadstemn_store = 0.0;
		canopy_strata[0].ns.livestemn_transfer = 0.0;
		canopy_strata[0].ns.deadstemn_transfer = 0.0;
		canopy_strata[0].ns.livecrootn_store = 0.0;
		canopy_strata[0].ns.deadcrootn_store = 0.0;
		canopy_strata[0].ns.livecrootn_transfer = 0.0;
		canopy_strata[0].ns.deadcrootn_transfer = 0.0;
		canopy_strata[0].ns.cwdn = 0.0;
	}
	/*--------------------------------------------------------------*/
	/* 	zero other carbon stores for non veg			*/
	/*--------------------------------------------------------------*/
	if (canopy_strata[0].defaults[0][0].epc.veg_type == NON_VEG) {
		canopy_strata[0].cs.cpool = 0.0;
		canopy_strata[0].cs.leafc = 0.0;
		canopy_strata[0].cs.dead_leafc = 0.0;
		canopy_strata[0].cs.leafc_store = 0.0;
		canopy_strata[0].cs.leafc_transfer = 0.0;
		canopy_strata[0].cs.frootc = 0.0;
		canopy_strata[0].cs.frootc_store = 0.0;
		canopy_strata[0].cs.frootc_transfer = 0.0;

		canopy_strata[0].ns.npool = 0.0;
		canopy_strata[0].ns.leafn = 0.0;
		canopy_strata[0].ns.dead_leafn = 0.0;
		canopy_strata[0].ns.leafn_store = 0.0;
		canopy_strata[0].ns.leafn_transfer = 0.0;
		canopy_strata[0].ns.frootn = 0.0;
		canopy_strata[0].ns.frootn_store = 0.0;
		canopy_strata[0].ns.frootn_transfer = 0.0;
	}

	/*--------------------------------------------------------------*/
	/*	zero all long term sinks				*/
	/*--------------------------------------------------------------*/
	
	canopy_strata[0].cs.gpsn_src = 0.0;
	canopy_strata[0].cs.leaf_mr_snk = 0.0;
	canopy_strata[0].cs.leaf_gr_snk = 0.0;
	canopy_strata[0].cs.livestem_mr_snk = 0.0;
	canopy_strata[0].cs.livestem_gr_snk = 0.0;
	canopy_strata[0].cs.deadstem_gr_snk = 0.0;
	canopy_strata[0].cs.livecroot_mr_snk = 0.0;
	canopy_strata[0].cs.livecroot_gr_snk = 0.0;
	canopy_strata[0].cs.deadcroot_gr_snk = 0.0;
	canopy_strata[0].cs.froot_mr_snk = 0.0;
	canopy_strata[0].cs.froot_gr_snk = 0.0;
	canopy_strata[0].NO3_stored = 0.0; // this is for the NO3 deposition on leaves
	
	

    canopy_strata[0].cs.Tacc = 20.0;	
	/*--------------------------------------------------------------*/
	/*	determine current lai and height  based on current leaf carbon	*/
	/* 	we need to initialize the sunlit/shaded proportions of LAI here */
	/*	(these will later be updated in update_phenology	*/
	/*	using Chen;s method					*/
	/*--------------------------------------------------------------*/
	canopy_strata[0].epv.proj_sla_sunlit = canopy_strata[0].defaults[0][0].epc.proj_sla;
	canopy_strata[0].epv.proj_sla_shade = canopy_strata[0].defaults[0][0].epc.proj_sla *
			canopy_strata[0].defaults[0][0].epc.shade_sla_mult;



	if ( canopy_strata[0].cs.leafc <= 1.0/canopy_strata[0].epv.proj_sla_sunlit) {
		canopy_strata[0].epv.proj_lai = canopy_strata[0].cs.leafc *
				canopy_strata[0].epv.proj_sla_sunlit;
		canopy_strata[0].epv.proj_lai_sunlit = canopy_strata[0].epv.proj_lai;
		canopy_strata[0].epv.proj_lai_shade = 0.0;
    } else {
		canopy_strata[0].epv.proj_lai = 1.0 + ( canopy_strata[0].cs.leafc -
				1.0/canopy_strata[0].epv.proj_sla_sunlit) *
				canopy_strata[0].epv.proj_sla_shade;
		canopy_strata[0].epv.proj_lai_sunlit = 1.0;
		canopy_strata[0].epv.proj_lai_shade = canopy_strata[0].epv.proj_lai - 1.0;
    }


	canopy_strata[0].epv.all_lai = canopy_strata[0].epv.proj_lai *
		canopy_strata[0].defaults[0][0].epc.lai_ratio;
	canopy_strata[0].epv.max_proj_lai =  canopy_strata[0].epv.proj_lai;
	
	if (canopy_strata[0].defaults[0][0].epc.veg_type == TREE)
		canopy_strata[0].epv.height =
		canopy_strata[0].defaults[0][0].epc.height_to_stem_coef
		* pow((canopy_strata[0].cs.live_stemc+canopy_strata[0].cs.dead_stemc),
		canopy_strata[0].defaults[0][0].epc.height_to_stem_exp);
	else
		canopy_strata[0].epv.height =
		canopy_strata[0].defaults[0][0].epc.height_to_stem_coef
		* pow((canopy_strata[0].cs.leafc + canopy_strata[0].cs.dead_leafc),
		canopy_strata[0].defaults[0][0].epc.height_to_stem_exp);

	/*--------------------------------------------------------------*/
	/*	calculate all sided  and project pai from max projected lai	*/
	/*--------------------------------------------------------------*/

	if (canopy_strata[0].defaults[0][0].epc.veg_type == TREE) {
		sai = 0.55*(1.0-exp(-0.175*(canopy_strata[0].cs.live_stemc+canopy_strata[0].cs.dead_stemc)));
		canopy_strata[0].epv.proj_pai = max(canopy_strata[0].epv.proj_lai + sai, 0.0);
		canopy_strata[0].epv.all_pai = max(canopy_strata[0].epv.all_lai + sai, 0.0);
	}
	else {
		canopy_strata[0].epv.proj_pai =  canopy_strata[0].epv.proj_lai;
		canopy_strata[0].epv.all_pai =  canopy_strata[0].epv.all_lai;
	}

	/*--------------------------------------------------------------*/
	/*	initializae turnovers and litterfall 			*/
	/*--------------------------------------------------------------*/
	if (compute_annual_turnover(canopy_strata[0].defaults[0][0].epc,
		&(canopy_strata[0].epv),
		&(canopy_strata[0].cs)) ){
		fprintf(stderr,"FATAL ERROR: in compute_annual_turnover() ... Exiting\n");
		exit(EXIT_FAILURE);
	}

	if (compute_annual_litfall(canopy_strata[0].defaults[0][0].epc,
		&(canopy_strata[0].phen),
		&(canopy_strata[0].cs), command_line[0].grow_flag) ){
		fprintf(stderr,"FATAL ERROR: in compute_annual_litfall() ... Exiting\n");
		exit(EXIT_FAILURE);
	}

	/*--------------------------------------------------------------*/
	/*	compute new rooting depth based on current root carbon (for growth mode)  */
	/* 	for static non-grow version use the worldfile rooting depth */
	/*	as read in above					*/
	/*--------------------------------------------------------------*/
    // ((command_line[0].grow_flag > 0) && (command_line[0].dynRtZoff_flag == 0))
	if (  command_line[0].iniBioRZ > 0 ){
		rootc = canopy_strata[0].cs.frootc+canopy_strata[0].cs.live_crootc+canopy_strata[0].cs.dead_crootc;
		if (rootc > ZERO){
			if (update_rooting_depth(
				&(canopy_strata[0].rootzone), 
				rootc, 
				canopy_strata[0].defaults[0][0].epc.root_growth_direction, 
				canopy_strata[0].defaults[0][0].epc.root_distrib_parm,
				patch[0].soil_defaults[0][0].maxrootdepth)){
				fprintf(stderr,
					"FATAL ERROR: in compute_rooting_depth() from construct_canopy_strata()\n");
				exit(EXIT_FAILURE);
            }//if
		}//if
	}//if

	/*--------------------------------------------------------------*/
	/*	initialize leaf out for non-grow version		*/
	/*--------------------------------------------------------------*/
//    if (( command_line[0].grow_flag == 0) &&
//        (canopy_strata[0].defaults[0][0].epc.veg_type != NON_VEG) ){
//        /*
//        canopy_strata[0].cs.leafc_transfer = canopy_strata[0].phen.leaflitfallc;
//        canopy_strata[0].ns.leafn_transfer = canopy_strata[0].phen.leaflitfallc
//            * canopy_strata[0].ns.leafn
//            / canopy_strata[0].cs.leafc;
//        canopy_strata[0].cs.leafc_store = 0.0;
//        canopy_strata[0].ns.leafn_store = 0.0;
//        */
//    }
	/*--------------------------------------------------------------*/
	/*	set phenology timing if static allocation		*/
	/*  and initialize for dynamic runs				*/
	/*--------------------------------------------------------------*/
		canopy_strata[0].phen.expand_startday = canopy_strata[0].defaults[0][0].epc.day_leafon;
		canopy_strata[0].phen.expand_stopday = canopy_strata[0].phen.expand_startday + canopy_strata[0].defaults[0][0].epc.ndays_expand;
		canopy_strata[0].phen.litfall_startday = canopy_strata[0].defaults[0][0].epc.day_leafoff;
		canopy_strata[0].phen.litfall_stopday = canopy_strata[0].phen.litfall_startday + canopy_strata[0].defaults[0][0].epc.ndays_litfall;
		if (canopy_strata[0].phen.expand_stopday > 365)
			canopy_strata[0].phen.expand_stopday -= 365;
		if (canopy_strata[0].phen.litfall_stopday > 365)
			canopy_strata[0].phen.litfall_stopday -= 365;
		/*---------------------------------------------------------------*/
		/* assume this is 365 for now since we don't know when next      */
		/* year's growing season will start                              */
		/*---------------------------------------------------------------*/
		canopy_strata[0].phen.nretdays = 365;
		canopy_strata[0].phen.gwseasonday = -1;
		canopy_strata[0].phen.lfseasonday = -1;
    
        canopy_strata[0].phen.future_gsi = malloc(sizeof(double)*NUM_fday_Pred);
        memset(canopy_strata[0].phen.future_gsi, 0, sizeof(double)*NUM_fday_Pred);
        canopy_strata[0].phen.GSI_vpd=0;
        canopy_strata[0].phen.GSI_tmin=0;
        canopy_strata[0].phen.GSI_dlen=0;
        canopy_strata[0].phen.gsi=0;
	/*--------------------------------------------------------------*/
	/*	set critical soil moisture (at stomatal closure)	*/
	/*      psi_close is converted to m water tension from MPa using     */
	/*      1m water tension = 10000 Pa                             */
	/*                      = 0.01 Mpa                              */
	/*--------------------------------------------------------------*/
	canopy_strata[0].epv.crit_vwc = exp(patch[0].soil_defaults[0][0].pore_size_index * (log(-100.0 * canopy_strata[0].defaults[0][0].epc.psi_close)-log(patch[0].soil_defaults[0][0].psi_air_entry)));
    //pow(( (-100.0 * canopy_strata[0].defaults[0][0].epc.psi_close) / patch[0].soil_defaults[0][0].psi_air_entry),patch[0].soil_defaults[0][0].pore_size_index );
        
    
	
	/*--------------------------------------------------------------*/
	/* initialize runnning average of psi **** should actually  calc */
	/* current day psi						*/
	/*--------------------------------------------------------------*/
	canopy_strata[0].epv.psi_ravg = canopy_strata[0].defaults[0][0].epc.psi_open;

	/*--------------------------------------------------------------*/
	/*	for now initialize these accumuling variables		*/
	/*	note that age really should be a state variable 	*/
	/*	and initialized in the worldfile			*/
	/*--------------------------------------------------------------*/
	canopy_strata[0].cs.num_resprout = 0;
	canopy_strata[0].cs.age = 0;
	canopy_strata[0].epv.wstress_days = 0;
	canopy_strata[0].epv.max_fparabs = 0.0;
	canopy_strata[0].epv.min_vwc = 1.0;
    
    canopy_strata[0].gDayCount = 0;
    canopy_strata[0].nFactor = 0.0;
    canopy_strata[0].wFactor = 0.0;
    canopy_strata[0].lFactor = 0.0;
    canopy_strata[0].gFactor = 0.0;
    canopy_strata[0].gwPSN = 0.0;
    canopy_strata[0].gwMResp = 0.0;
    canopy_strata[0].gwAPAR = 0.0;
    canopy_strata[0].gwLWP = 0.0;
    canopy_strata[0].gwVPD = 0.0;
	/*--------------------------------------------------------------*/
	/*	Read in the number of  strata base stations 					*/
	/*--------------------------------------------------------------*/
//	fscanf(world_file,"%d",&(canopy_strata[0].num_base_stations));
//	read_record(world_file, record);
    canopy_strata[0].num_base_stations = getIntWorldfile(&paramCnt,&paramPtr,"canopy_strata_n_basestations","%d",0,1);//57
	/*--------------------------------------------------------------*/
	/*    Allocate a list of base stations for this strata.			*/
	/*--------------------------------------------------------------*/
    canopy_strata[0].base_stations = NULL;
//    canopy_strata[0].base_stations = (struct base_station_object **)
//		alloc(canopy_strata[0].num_base_stations *
//		sizeof(struct base_station_object *),"base_stations",
//		"construct_canopy_strata");
//	/*--------------------------------------------------------------*/
//	/* Read each base_station ID and then point to that base_station*/
//	/*--------------------------------------------------------------*/
//	for (i=0 ; i<canopy_strata[0].num_base_stations; i++){
//		fscanf(world_file,"%d",&(base_stationID));
//		read_record(world_file, record);
//		/*--------------------------------------------------------------*/
//		/*	Point to the appropriate base station in the base       	*/
//		/*              station list for this world.					*/
//		/*--------------------------------------------------------------*/
//		canopy_strata[0].base_stations[i] = assign_base_station(
//			base_stationID,
//			num_world_base_stations,
//			world_base_stations);
//	} /*end for*/
	return(canopy_strata);
} /*end construct_canopy_strata.c*/
