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
		struct	basin_object *,
        int hID,
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
    struct patch_object* MyPatch;
    int *aggregate_ID;
    
    float *denitri, *nitrif, *Nuptake, *Nnmineral, *soilNO3, *soilNH4, *infiltrate, *wtz;
    //float *frmStrQ, *frmStrNO3, *frmStrNH4; //, *frmStrDON;
    //float *frmRipQ, *frmRipNO3, *frmRipNH4; //, *frmRipDON;
    //float *frmLndQ, *frmLndNO3, *frmLndNH4; //, *frmLndDON;
    float *Pot_denitrif_SS, *Pot_denitrif_CO2, *exfiltration;
    float *sat_q, *satNO3, *satNH4;
    //float *satDON, *satDOC;
    float *ZONEdirectPARa, *ZONEdiffusePARa;
    float *AETa, *PETa, *AREAa, *aPSN, *aroots, *aLAI, *directPARa, *diffusePARa, area_1;
    float *nlimit, *qlimit;
    float *activeS;
    float *litterLigninC, *litterCelluloseC, *litterLabileC;
    float *litterLigninN, *litterCelluloseN, *litterLabileN;
    if(command_line[0].aggregate_flag>0){
        aggregate_ID = (int*)calloc(world[0].basins[0][0].aggregateLength, sizeof(int));
        
        AETa = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        PETa = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        AREAa = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        aPSN = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        aroots = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        aLAI = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        nlimit = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        qlimit = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        directPARa = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        diffusePARa = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        ZONEdirectPARa = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        ZONEdiffusePARa = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        
        litterLigninN = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        litterCelluloseN = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        litterLabileN = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        litterLigninC = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        litterCelluloseC = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        litterLabileC = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        
        activeS = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        
        denitri = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        Pot_denitrif_SS = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        Pot_denitrif_CO2 = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        nitrif = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        Nuptake = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        Nnmineral = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        soilNO3 = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        soilNH4 = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        infiltrate = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        exfiltration = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        wtz = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        sat_q = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        satNO3 = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        satNH4 = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        satDON = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        satDOC = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        
//        frmStrQ = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        frmStrNO3 = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        frmStrNH4 = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        frmRipQ = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        frmRipNO3 = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        frmRipNH4 = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        frmLndQ = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        frmLndNO3 = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
//        frmLndNH4 = (float*)calloc(world[0].basins[0][0].aggregateLength, sizeof(float));
        
        
        for(i = 0; i < world[0].basins[0][0].aggregateLength; i++){
            aggregate_ID[i] = 0;
            AETa[i] = 0.0;
            PETa[i] = 0.0;
            AREAa[i] = 0.0;
            aPSN[i] = 0.0;
            aroots[i] = 0.0;
            aLAI[i] = 0.0;
            nlimit[i] = 0.0;
            qlimit[i] = 0.0;
//            directPARa[i] = 0.0;
//            diffusePARa[i] = 0.0;
//            ZONEdirectPARa[i] = 0.0;
//            ZONEdiffusePARa[i] = 0.0;
            
            litterLigninN[i] = 0.0;
            litterCelluloseN[i] = 0.0;
            litterLabileN[i] = 0.0;
            litterLigninC[i] = 0.0;
            litterCelluloseC[i] = 0.0;
            litterLabileC[i] = 0.0;
            
            activeS[i] = 0.0;
            
            denitri[i] = 0.0;
            Pot_denitrif_SS[i] = 0.0;
            Pot_denitrif_CO2[i] = 0.0;
            nitrif[i] = 0.0;
            Nuptake[i] = 0.0;
            Nnmineral[i] = 0.0;
            soilNO3[i] = 0.0;
            soilNH4[i] = 0.0;
            infiltrate[i] = 0.0;
            exfiltration[i] = 0.0;
            wtz[i] = 0.0;
            sat_q[i] = 0.0;
            satNO3[i] = 0.0;
            satNH4[i] = 0.0;
//            satDON[i] = 0.0;
//            satDOC[i] = 0.0;
//
//            frmStrQ[i] = 0.0;
//            frmStrNO3[i] = 0.0;
//            frmStrNH4[i] = 0.0;
//            frmRipQ[i] = 0.0;
//            frmRipNO3[i] = 0.0;
//            frmRipNH4[i] = 0.0;
//            frmLndQ[i] = 0.0;
//            frmLndNO3[i] = 0.0;
//            frmLndNH4[i] = 0.0;
        }//i
        
        for(i = 0; i < world[0].basins[0][0].route_list->num_patches; i++){

            MyPatch = world[0].basins[0][0].route_list->list[i];
            
            aggregate_ID[MyPatch->aggregate_index] = MyPatch->aggregate_ID;
            
            AREAa[MyPatch->aggregate_index] += MyPatch[0].area;
            PETa[MyPatch->aggregate_index] += MyPatch[0].PET * MyPatch[0].area;
            AETa[MyPatch->aggregate_index] += (MyPatch[0].evaporation + MyPatch[0].evaporation_surf + MyPatch[0].exfiltration_sat_zone + MyPatch[0].exfiltration_unsat_zone + MyPatch[0].transpiration_sat_zone + MyPatch[0].transpiration_unsat_zone)*MyPatch[0].area;
            
//            directPARa[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].PAR_direct;
//            diffusePARa[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].PAR_diffuse;
//
//            ZONEdirectPARa[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].zone[0].PAR_direct;
//            ZONEdiffusePARa[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].zone[0].PAR_diffuse;
            
            // standing litter
            litterLigninN[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].litter_ns.litr4n;
            litterCelluloseN[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].litter_ns.litr2n + MyPatch[0].litter_ns.litr3n);
            litterLabileN[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].litter_ns.litr1n;

            litterLigninC[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].litter_cs.litr4c;
            litterCelluloseC[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].litter_cs.litr2c + MyPatch[0].litter_cs.litr3c);
            litterLabileC[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].litter_cs.litr1c;

            // leaffall litter
//            litterLigninN[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].ndf.leafn_to_litr4n + MyPatch[0].ndf.frootn_to_litr4n);
//            litterCelluloseN[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].ndf.leafn_to_litr2n  + MyPatch[0].ndf.leafn_to_litr3n + MyPatch[0].ndf.frootn_to_litr2n  + MyPatch[0].ndf.frootn_to_litr3n);
//            litterLabileN[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].ndf.leafn_to_litr1n + MyPatch[0].ndf.frootn_to_litr1n);
//
//            litterLigninC[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].cdf.leafc_to_litr4c + MyPatch[0].cdf.frootc_to_litr4c);
//            litterCelluloseC[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].cdf.leafc_to_litr2c  + MyPatch[0].cdf.leafc_to_litr3c + MyPatch[0].cdf.frootc_to_litr2c  + MyPatch[0].cdf.frootc_to_litr3c);
//            litterLabileC[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].cdf.leafc_to_litr1c + MyPatch[0].cdf.frootc_to_litr1c);

          
            
            aroots[MyPatch->aggregate_index] += MyPatch[0].rootzone.SatPct *MyPatch[0].area;
            nlimit[MyPatch->aggregate_index] += MyPatch[0].soil_ns.nlimit *MyPatch[0].area;
            qlimit[MyPatch->aggregate_index] += MyPatch[0].trans_reduc_perc *MyPatch[0].area;
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
     
            denitri[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].ndf.denitrif;
            Pot_denitrif_SS[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].ndf.Pot_denitrif_SS;
            Pot_denitrif_CO2[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].ndf.Pot_denitrif_SS;
            nitrif[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].ndf.sminn_to_nitrate;
            Nuptake[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].ndf.plant_avail_uptake;
            Nnmineral[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].ndf.net_mineralized;
            soilNO3[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].soil_ns.nitrate;
            soilNH4[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].soil_ns.sminn;
            infiltrate[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].recharge; // overland_flow;
            exfiltration[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].overland_flow;
            wtz[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].sat_deficit_z;
            
            
            //MyPatch[0].rootzone.potential_sat*1000
            //MyPatch[0].zone[0].PAR_direct+MyPatch[0].zone[0].PAR_diffuse;
            
            // estimate of hyperhric Q
//            double wtz = max(MyPatch[0].sat_deficit_z, 0.0);
//            double wtz2m = wtz+2.0;
//            double p0 = MyPatch[0].soil_defaults[0][0].porosity_0;
//            double qq = MyPatch[0].soil_defaults[0][0].porosity_decay
//            sat_q[MyPatch->aggregate_index] += MyPatch[0].area * p0 * qq * (exp(-wtz/qq) - exp(-wtz2m/qq));//m just sample water 2.0 m down from wtz
            
            // how much water in sat zone
             sat_q[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].available_soil_water;//m
            
            
            
            satNO3[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].sat_NO3;//kgN
            satNH4[MyPatch->aggregate_index] += MyPatch[0].area * MyPatch[0].sat_NH4;//kgN
            
            // S in active zone
            if( MyPatch[0].soil_defaults[0][0].active_zone_z > MyPatch[0].sat_deficit_z){
                activeS[MyPatch->aggregate_index] += ((MyPatch[0].rz_storage + MyPatch[0].unsat_storage + MyPatch[0].soil_defaults[0][0].active_zone_sat_0z - MyPatch[0].sat_deficit) * MyPatch[0].soil_defaults[0][0].active_zone_sat_0z_1) * MyPatch[0].area;
            }else if(MyPatch[0].soil_defaults[0][0].active_zone_z > MyPatch[0].rootzone.depth){
                activeS[MyPatch->aggregate_index] += ((MyPatch[0].rz_storage+MyPatch[0].unsat_storage) / MyPatch[0].sat_deficit) * MyPatch[0].area; // approximate
            }else{
                activeS[MyPatch->aggregate_index] += (MyPatch[0].rz_storage/MyPatch[0].rootzone.potential_sat) * MyPatch[0].area;
            }
           
            
            
//            frmStrQ[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].fromSTREAM_Q + MyPatch[0].fromSTREAM_surfsubQ);
//            frmStrNO3[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].fromSTREAM_NO3 + MyPatch[0].fromSTREAM_surfsubNO3);
//            frmStrNH4[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].fromSTREAM_NH4 + MyPatch[0].fromSTREAM_surfsubNH4);
//            frmRipQ[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].fromRIPARIAN_Q + MyPatch[0].fromRIPARIAN_surfsubQ);
//            frmRipNO3[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].fromRIPARIAN_NO3 + MyPatch[0].fromRIPARIAN_surfsubNO3);
//            frmRipNH4[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].fromRIPARIAN_NH4 + MyPatch[0].fromRIPARIAN_surfsubNH4);
//            frmLndQ[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].fromLAND_Q + MyPatch[0].fromRIPARIAN_surfsubQ);
//            frmLndNO3[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].fromRIPARIAN_NO3 + MyPatch[0].fromRIPARIAN_surfsubNO3);
//            frmLndNH4[MyPatch->aggregate_index] += MyPatch[0].area * (MyPatch[0].fromRIPARIAN_NH4 + MyPatch[0].fromRIPARIAN_surfsubNH4);
            
        }//i for
        
        // <----- here print out
        for(i = 0; i < world[0].basins[0][0].aggregateLength; i++){
            area_1 = 1.0/AREAa[i];
            fprintf(outfile->aggregate->daily,"%d,%d,%d,%d, %e,%e,%e,%e,%e,%e,%e,%e, %e,%e,%e,%e,%e,%e, %e,%e,%e,%e,%e, %e,%e,%e,%e,%e, %e,%e,%e,%e,%e\n",
                    date.day, date.month, date.year, aggregate_ID[i],
                    
                    AREAa[i], //m2
                    AETa[i]*1000.0*area_1, // mm
                    PETa[i]*1000.0*area_1,
                    aPSN[i]*area_1,
                    aroots[i]*area_1,
                    aLAI[i]*area_1,
                    nlimit[i]*area_1,
                    qlimit[i]*area_1,
//                    directPARa[i]*area_1,
//                    diffusePARa[i]*area_1,
//                    ZONEdirectPARa[i]*area_1,
//                    ZONEdiffusePARa[i]*area_1,//[11]
                    
                    litterLigninN[i]*1000.0*area_1, // gC/m2
                    litterCelluloseN[i]*1000.0*area_1, // gN/m2
                    litterLabileN[i]*1000.0*area_1,
                    litterLigninC[i]*1000.0*area_1,
                    litterCelluloseC[i]*1000.0*area_1,
                    litterLabileC[i]*1000.0*area_1,//[6]
                    
                    activeS[i]*area_1,
                    denitri[i]*1000.0*area_1, // gN/m2/day
                    Pot_denitrif_SS[i]*1000.0*area_1,
                    Pot_denitrif_CO2[i]*1000.0*area_1,
                    nitrif[i]*1000.0*area_1,
                    Nuptake[i]*1000.0*area_1,
                    Nnmineral[i]*1000.0*area_1,
                    soilNO3[i]*1000.0*area_1,
                    soilNH4[i]*1000.0*area_1,
                    infiltrate[i]*1000.0*area_1, // mm/day
                    exfiltration[i]*1000.0*area_1,
                    wtz[i]*1000.0*area_1, // mm
                    sat_q[i]*1000.0*area_1,//mm/2
                    satNO3[i]*1000.0*area_1,//gN/m2 //[13]
                    satNH4[i]*1000.0*area_1//gN/m2 //[13]
//                    satDON[i]*1000.0*area_1,//gN/m2 //[13]
//                    satDOC[i]*1000.0*area_1,//gN/m2 //[13]
                    
//                    frmStrQ[i]*1000.0*area_1, // mm
//                    frmStrNO3[i]*1000.0*area_1,//gN/m2/day
//                    frmStrNH4[i]*1000.0*area_1,
//                    frmRipQ[i]*1000.0*area_1, // mm
//                    frmRipNO3[i]*1000.0*area_1,
//                    frmRipNH4[i]*1000.0*area_1,
//                    frmLndQ[i]*1000.0*area_1, // mm
//                    frmLndNO3[i]*1000.0*area_1,
//                    frmLndNH4[i]*1000.0*area_1
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
								world[0].basins[b],
                                world[0].basins[b]->hillslopes[h]->ID,
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
													if ( (world[0].basins[b][0].hillslopes[h][0].zones[z][0].patches[p][0].ID == patchID)
														|| (patchID == -999 && world[0].basins[b]->hillslopes[h]->zones[z]->patches[p]->landuse_defaults[0]->ID >= command_line[0].patchPrintTh) ){
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
