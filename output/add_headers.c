/*--------------------------------------------------------------*/
/* 																*/
/*					add_headers					*/
/*																*/
/*	add_headers - 
															*/
/*	NAME														*/
/*	add_headers 
																*/
/*	SYNOPSIS													*/
/*	void add_headers(struct world output_file_object *,				*/
/*			struct command_line_object *)					*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	Adds headers for yearly, monthly, daily and	*/
/*	hourly basin output 					*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*																*/
/*--------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include "rhessys.h"


void add_headers(struct world_output_file_object *world_output_files, 
			struct command_line_object *command_line)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	FILE *outfile;
	int check;
    
    if (command_line[0].aggregate_flag>0) {
        outfile = world_output_files[0].aggregate[0].daily;
        fprintf(outfile,"%s,%s,%s,%s, %s, %s,%s,%s,%s,%s,%s,%s, %s,%s,%s,%s,%s,%s, %s,%s,%s,%s,%s, %s,%s,%s,%s,%s, %s,%s,%s,%s,%s\n",
                "day", "month", "year", "aggregateID",//[4]
                
                "area",
                "ET",
                "PET",
                "psn",
                "rootS",
                "LAI",
                "nlimit",
                "qlimit",
//                "PARdirect",
//                "PARdiffuse",
//                "ZONEPARdirect",
//                "ZONEPARdiffuse",//[11]
                
                "litterLigninN",
                "litterCelluloseN",
                "litterLabileN",
                "litterLigninC",
                "litterCelluloseC",
                "litterLabileC",//[6]
                
                "activeS",
                "denitri",
                "Pot_denitrif_SS",
                "Pot_denitrif_CO2",
                "nitrif",
                "Nuptake",
                "Nnmineral",
                "soilNO3",
                "soilNH4",
                "infiltration",
                "exfiltration",
                "wtz",
                "satq",
                "satNO3",//[13]
                "satNH4"
//                "satDON",
//                "satDOC",
//
//                "frmStrQ",
//                "frmStrNO3",
//                "frmStrNH4",
//                "frmRipQ",
//                "frmRipNO3",
//                "frmRipNH4",
//                "frmLndQ",
//                "frmLndNO3",
//                "frmLndNH4"//[9]
                );
                
    }//
    
	/*--------------------------------------------------------------*/
	/*	Basin file headers					*/
	/*--------------------------------------------------------------*/

	if (command_line[0].b != NULL) {
	outfile = world_output_files[0].basin[0].hourly;
	fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \n",
	// the unit is based on mm and day
		"hour",		
		"day",
		"month",
		"year",
		"basinID",
		"pot_surface_infil",
		//asnow_throughfall * 1000.0,
		"sat_def_z",
		"sat_def",
		"rz_stor",
		"unsat_stor",
		"rz_drainage",
		"unsat_drainage",
		//acap_rise * 1000.0,
		//aevaporation * 1000.0,
		//asnowpack * 1000.0,
		//atranspiration * 1000.0,
		"subsur2stream_flow",
		"sur2stream_flow",
		"streamflow",
		//apsn,
		//alai,
		"gw.Qout",
		"gw.storage",
		"detention_store",
		"%sat_area",
		"litter_store",
		"canopy_store", 
		//aperc_snow *100,
		//asublimation * 1000.0,
		//var_trans,
		//aacctrans*1000,
		//var_acctrans,
		//aPET*1000,
		//adC13, 
		"precip", 
		//amortality_fract*100,
	  	//atmax, 
		//atmin, 
		//asnow*1000.0 ,
		"routedstreamflow");
		
	









	/*--------------------------------------------------------------*/
	/*	Daily 							*/
	/*--------------------------------------------------------------*/
    ///<<<---------- here basin daily
	outfile = world_output_files[0].basin[0].daily;
	fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" , // added 3 extra
		"day",
		"month",
		"year",
		"basinID",
		"rain_thr",
		"snow_thr",
		"sat_def_z",
		"sat_def",
		"rz_storage",
		"unsat_stor",
		"rz_drainage",
		"unsat_drain",
		"cap",
		"evap",
		"snowpack",
		"trans",
		"baseflow",
		"return",
		"streamflow",
		"psn",
		"laiTREE",
		"gw.Qout",
		"gw.storage",
		"detention_store",
		"%sat_area",
		"litter_store",
		"canopy_store",
		"%snow_cover",
		"snow_subl",
		"trans_var",
		"acc_trans",
		"acctransv_var",
		"pet",
		"dC13",
		"precip",
		"pcp_assim",
		"mortf",
		"tmax",
		"tmin",
		"tavg",
		"vpd",
		"snowfall",
		"recharge",
		"gpsn",
		"resp",
		"gs",
		"rootdepth",
		"plantc",
		"snowmelt",
		"canopysubl",
		"routedstreamflow",
		"canopy_snow",
		"height",
		"evap_can","evap_lit","evap_soil",
		"litrc",
		"Kdown","Ldown","Kup","Lup",
		"Kstar_can","Kstar_soil","Kstar_snow",
		"Lstar_can","Lstar_soil","Lstar_snow",
        "LE_canopy","LE_soil","LE_snow","Lstar_strat","canopydrip","ga","srad","dlen",
        "gsi","gsi_vpd","gsi_dlen","gsi_tmin","cwdc","leafc","stemc","frootc",
        "stormdrain",
            "stormdrainNO3",
            "stormdrainNH4",
            "stormdrainDON",
            "stormdrainDOC",
        "sewerdrain",
            "asewerdrainNO3",
            "asewerdrainNH4",
            "asewerdrainDON",
            "asewerdrainDOC",
        "pipedrain",
            "apipedrainNO3",
            "apipedrainNH4",
            "apipedrainDON",
            "apipedrainDOC",
        "lawnIrrigated",
        "septicQ",
        "laiNontree",
        "PAR",
        "unsat_cap",
        "unsat_fc",
        "rtz_fc"); // basin daily

	/*--------------------------------------------------------------*/
	/*	Monthly							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].basin[0].monthly;
	check = fprintf(outfile,
		"%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", 
		"month",
		"year",
		"basinID",
		"streamflow",
		"streamflow_NO3",
		"denitrif",
		"DOC",
		"DON",
		"et",
		"psn",
		"lai",
		"nitrif",
		"mineralized",
		"uptake");
	/*--------------------------------------------------------------*/
	/*	Yearly 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].basin[0].yearly;
	check = fprintf(outfile, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
        "year",
        "subQnet",
        "surfQnet",
        "subQvnet",
        "precip",
        "pet",
        "et",
        "sat_deficit_z",
        "peakLAI",
        "meanLAI",
        "psn",
        "denitrif",
        "mineralization","uptake","subNO3net", "subNO3vnet","subDOCnet",
        "no3drain2gw"
        );
    }// basin
    
	/*--------------------------------------------------------------*/
	/*	Hillslope file headers					*/
	/*--------------------------------------------------------------*/
	if (command_line[0].h != NULL) {
	/*--------------------------------------------------------------*/
	/*	Daily 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].hillslope[0].daily;
	//fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n,"
    fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
		"day", //1
		"month", //2
		"year", //3
		"hillID", //5
        "area", //6
		"sat_def_z", //7
		"sat_def", //8
        "detention_store", //9
        "sat_area",//10
		"rz_storage",//11
		"cap", //12
        "drainage",//13
        "baseflow",//14
		"return",//15
		"streamflow",//16
        "gw.Qout",//17
        "gw.storage",//18
        "snowmelt",//19
		"psn",//20
        "evap",//21
        "trans",//22
		"laiTREE",//23
        "laiGRASS",//24
        "lawnirrigated",//25
        "precip", // 26
        "infiltration" //27
		);

	/*--------------------------------------------------------------*/
	/*	Monthly							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].hillslope[0].monthly;
	check = fprintf(outfile,
		"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", 
		"month",
		"year",
		"basinID",
		"hillslopeID",
		"streamflow",
		"streamflow_NO3",
		"snowpack",
		"denitrif",
		"DOC",
		"DON",
		"et",
		"psn",
		"lai",
		"nitrif",
		"mineralized",
		"uptake","area");
	/*--------------------------------------------------------------*/
	/*	Yearly 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].hillslope[0].yearly;
        
        fprintf(outfile, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
        "year",
        "patchID",
        "subQnet",
        "surfQnet",
        "subQvnet", //5
        "precip",
        "pet",
        "et",
        "sat_deficit_z",
        "peakLAI", //10
        "meanLAI",
        "psn",
        "denitrif",
        "mineralization","uptake","subNO3net", "subNO3vnet","subDOCnet"
        );
	}
	/*--------------------------------------------------------------*/
	/*	Zone file headers					*/
	/*--------------------------------------------------------------*/
	if (command_line[0].z != NULL) {
	/*--------------------------------------------------------------*/
	/*	Daily 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].zone[0].daily;
	fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n " ,
		"day",
		"month",
		"year",
		"basinID",
		"hillID",
		"ID",
		"rain",
		"snow",
		"tmax",
		"tmin",
		"vpd",
		"Kdown_direct",
		"Kdown_diffuse",
		"PAR_direct",
		"PAR_diffuse",
		"Ldown",
		"relH","aspect","z","slope","ehr","whr",
		"tdew","edew",
		"transmis",
		"wind",
		"deltaT","clearskytransmis","tcoeff1","cloudfrac");

	/*--------------------------------------------------------------*/
	/*	Monthly							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].zone[0].monthly;
	check = fprintf(outfile,
		"%s %s %s %s %s %s %s %s %s %s\n" ,
		"month",
		"year",
		"basinID",
		"hillID",
		"zoneID",
		"precip",
		"K_direct",
		"K_diffuse",
		"tmax", "tmin");

	/*--------------------------------------------------------------*/
	/*	Hourly 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].zone[0].hourly;
	fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n " ,
		"day",
		"month",
		"year",
		"hour",
		"basinID",
		"hillID",
		"ID",
		"rain",
		"snow",
		"tday",
		"tavg",
		"vpd",
		"Kdown_direct",
		"Kdown_diffuse",
		"PAR_direct",
		"PAR_diffuse");



	}

	/*--------------------------------------------------------------*/
	/*	Patch file headers					*/
	/*--------------------------------------------------------------*/
	if (command_line[0].p != NULL) {
	/*--------------------------------------------------------------*/
	/*	Daily 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].patch[0].daily;
		check = fprintf(outfile,
//                        "%s-%s-%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
                        "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", // added 2 extra
                        
                        // we are looking for hydrology behavoirs
						"year", //1
						"month", //2
						"day", //3
						"patchID", //4
                        "subsurfaceQnet", //5 Qout-Qin  [+ = source; - = sink]
                        "surfaceQnet", //6 surface_Qout(sum of locally return_flow and locally rain on surface) - surface_Qin
                        "detention",//7
                        "stormdrainYield",//8
                        "return", //9 overland_flow(customized)
                        "rain_thr", //10
                        "thr_recharge", //11 rain_thr - recharge -  [+ = source; - = sink]
                        "cap_drain", //12 patch[0].cap_rise - patch[0].unsat_drainage  [+ = rise; - = down]
                        "sat_def_z", //13
                        "sat_def", //14 (sat_def>0)? (rz_storage+unsat_stor)/sat_def : -1
                        "rtzStorage", //15 (sat_def>0)? rz_storage/potential_rz_store : -1
                        "ET",
                        "treeLAI",//17
                        "nontreeLAI",//18
                        "SmartIrrigation", //19
                        "rtz_totalvol", //20
                        "unsat_fc", //21
                        "rtz_fc", //22
                        "unsat_storage", //23
                        "top12cm_storage", //24
                        "top12cm_potential_sat", //25
                        "rootdepth", //26
                        "soildepth", //27
                        "top30cm_storage", //28
                        "top30cm_potential_sat", //29
                        "top60cm_storage", //30
                        "top60cm_potential_sat", //31
			"ex_inundation_depth", //32
			"ex_inundation_dur", //33
						); 
                        // patch daily
	/*--------------------------------------------------------------*/
	/*	Monthly							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].patch[0].monthly;
	check = fprintf(outfile,
        "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
		"year",
        "month",
		"patchID",
		"subQnet",
		"surfQnet",//5
        "subQvnet",
		"precip",
		"pet",
        "et",
		"sat_deficit_z",//10
		"peakLAI",
        "meanLAI",
        "psn",
        "denitrif",
        "mineralization",//15
        "uptake",
        "subNO3net",
        "subNO3vnet",
        "subDOCnet",
        "no3drain2gw",//20
        "satChance",
        "plantlimitN",
        "plantlimitQ");
	/*--------------------------------------------------------------*/
	/*	Yearly							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].patch[0].yearly;
	fprintf(outfile, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
            "year",
            "patchID",
            "subQnet",
            "surfQnet",
            "subQvnet",//5
            "precip",
            "pet",
            "et",
            "sat_deficit_z",
            "peakLAI",//10
            "meanLAI",
            "psn",
            "denitrif",
            "mineralization",
            "uptake",//15
            "subNO3net",
            "subNO3vnet",
            "subDOCnet",
            "no3drain2gw",
            "satChance",//20
            "plantlimitN",
            "plantlimitQ");
	}

	/*--------------------------------------------------------------*/
	/*	Stratum file headers					*/
	/*--------------------------------------------------------------*/
	if (command_line[0].c != NULL) {
	/*--------------------------------------------------------------*/
	/*	Daily 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].canopy_stratum[0].daily;
	fprintf(outfile,
		"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
		"day",
		"month",
		"year",
		"basinID",
		"hillID",
		"zoneID",
		"patchID",
		"stratumID",
		"lai",
		"evap",
		"K_direct",
		"K_diffuse",
		"sublim",
		"trans",
		"ga",
		"gsurf",
		"gs",
		"psi",
		"leaf_day_mr",
		"psn_to_cpool",
		"rain_stored",
		"snow_stored",
		"rootzone.S",
		"m_APAR","m_tavg","m_LWP","m_CO2","m_tmin","m_vpd","dC13",
		"Lstar","surf_heat",
		"height","covfrac","vegID","PAR_direct","PAR_diffuse","PET","PPSN","ppfd_sunlit","ppfd_shade");
	/*--------------------------------------------------------------*/
	/*	Monthly							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].canopy_stratum[0].monthly;
	fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s \n", 
		"month",
		"year",
		"basinID",
		"hillID",
		"zoneID",
		"patchID",
		"stratumID",
		"lai",
		"psn",
		"lwp");
	/*--------------------------------------------------------------*/
	/*	Yearly							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].canopy_stratum[0].yearly;
	fprintf(outfile,"%s %s %s %s %s %s %s %s %s\n",
		"year",
		"basinID",
		"hillID",
		"zoneID",
		"patchID",
		"stratumID",
		"psn",
		"lwp","root_depth");
	}
	/*--------------------------------------------------------------*/
	/*	Stream routing file headers					*/
	/*--------------------------------------------------------------*/
	if (command_line[0].stro != NULL) {
		/*--------------------------------------------------------------*/
		/*	Daily 							*/
		/*--------------------------------------------------------------*/
		
        outfile = world_output_files[0].stream_routing[0].daily;
		fprintf(outfile,
				"%s %s %s %s %s\n" ,
				"day",
				"month",
				"year",
				"reachID",
				"routedstreamflow");
	}	
	return;
} /*end add_headers*/
