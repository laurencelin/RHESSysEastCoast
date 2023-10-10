/*--------------------------------------------------------------*/
/* 																*/
/*					add_growth_headers					*/
/*																*/
/*	add_growth_headers - 												    	*/
/*	NAME														*/
/*	add_growth_headers    													*/
/*	SYNOPSIS													*/
/*	void add_growth_headers(struct world output_file_object *,				*/
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


void add_growth_headers(struct world_output_file_object *world_output_files, 
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
	/*--------------------------------------------------------------*/
	/*	Basin file headers					*/
	/*--------------------------------------------------------------*/

	if (command_line[0].b != NULL) {

	/*--------------------------------------------------------------*/
	/*	Hourly							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].basin[0].hourly;
	fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
		"hour",
		"day",
		"month",
		"year",
		"basinID",
		//"lai",
		//"gpsn",
		//"plant_resp",
		//"soil_resp",
		//"nitrate",
		//"sminn",
		//"surfaceN",
		//"plantc",
		//"plantn",
		//"npool",
		//"litrc",
		//"litrn",
		//"soilc",
		//"soiln",
		"gwNO3",
		"gwNH4",
		"gwDON",
		"gwDOC",
		"streamflow_NO3",
		"streamflow_NH4",
		"streamflow_DON",
		"streamflow_DOC",
		"gwNO3out",
		"gwNH4out",
		"gwDONout",
		"gwDOCout",
		//"denitrif",
		//"nitrif",
		//"DOC",
		//"DON",
		//"root_depth",
		//"nfix",
		//"nuptake",
		//"grazingC",
		"StreamNO3_from_surface",
		"StreamNO3_from_sub");	  
	/*--------------------------------------------------------------*/
	/*	Daily 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].basin[0].daily;

	fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s, %s %s %s %s, %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
		"day",
		"month",
		"year",
		"basinID",
		"lai",
		"gpsn",
		"plant_resp",
		"soil_resp",
		"nitrate",
		"sminn",
		"surfaceN",
		"plantc",
		"plantn",
		"cpool",
		"npool",
		"litrc",
		"litrn",
		"soilc",
		"soiln",
		"gwNO3",
		"gwNH4",
		"gwDON",
		"gwDOC",
		"streamflow_NO3",
		"streamflow_NH4",
		"streamflow_DON",
		"streamflow_DOC",
		"gwNO3out",
		"gwNH4out",
		"gwDONout",
		"gwDOCout",
		"denitrif",
		"nitrif",
		"DOC",
        "surfaceDOC",
		"DON",
		"root_depth",
		"nfix",
		"nuptake",
        "netMineral",
        "immob",
        "mineralized",
		"grazingC",
		"StreamNO3_from_surface",
		"StreamNO3_from_sub",
        "pmnfl1s1",
        "pmnfl2s2",
        "pmnfl3l2",
        "pmnfl4s3",
        "pmnfs1s2",
        "pmnfs2s3",
        "pmnfs3s4",
        "pmnfs4",
        "soil3c",
        "soil4c",
        "nlimit",
        "litterNO3stored",
        "stratumNO3stored",
        "surfaceDIN",
        "rain_throughfall",
        "NO3_throughfall",
        "ndeposition"
        );
	/*--------------------------------------------------------------*/
	/*	Yearly 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].basin[0].yearly;
	fprintf(outfile, "%s %s %s %s %s %s %s %s %s %s \n",
		"year",
		"basinID",
		"gpsn",
		"plantresp",
		"newC",
		"soilhr",
		"strN",
		"denitrif","root_depth","mortf");
		
	}

	/*--------------------------------------------------------------*/
	/*	Hillslope file headers					*/
	/*--------------------------------------------------------------*/
	if (command_line[0].h != NULL) {
	/*--------------------------------------------------------------*/
	/*	Daily 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].hillslope[0].daily;

	fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
		"day", //1
		"month", //2
		"year", //3
		"hillID", //4
		"gpsn", //5
		"plant_resp", //6
		"soil_resp", //7
		"nitrate", //8
		"sminn", //9
		"surfaceN", //10
		"plantc", //11
		"plantn", //12
		"litrc", //13
		"litrn", //14
		"soilc", //15
		"soiln", //16
		"gwNO3", //17
		"gwNH4", //18
		"gwDON", //19
		"gwDOC", //20
		"streamflow_NO3", //21
		"streamflow_NH4", //22
		"streamflow_DON", //23
		"streamflow_DOC", //24
		"gwNO3out", //25
		"gwNH4out", //26
		"gwDONout", //27
		"gwDOCout", //28
		"denitrif", //29
		"nitrif", //30
		"DOC", //31
		"DON", //32
		"root_depth", //33
		"nfix", //34
		"nuptake"); //35
	}

	/*--------------------------------------------------------------*/
	/*	Zone file headers					*/
	/*--------------------------------------------------------------*/
	if (command_line[0].z != NULL) {
	/*--------------------------------------------------------------*/
	/*	Daily 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].zone[0].daily;
	fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n ", 
		"day",
		"month",
		"year",
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
                    //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
                    "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
                    "year", //1
                    "month", //2
                    "day", //3
                    "patchID", //4
                    
                    "soilNO3",//5
                    "satNO3",//6
                    "subNO3outin",//7
                    
                    "soilNH4",//5
                    "satNH4",//6
                    "subNH4outin",//7
                    
                    "soilDOC",//8
                    "satDOC",//9
                    "subDOCoutin",//10
                    "soilc",//11
                    "soiln",//12
                    "denitrif",//13
                    "nitrifi",//14
                    "uptake",//15
                    "immob",//16
                    "mineral",//17
                    "psn",//18
                    "plant_resp",//19
                    "soil_resp",//20
                    "cFrac",//21
                    "gDayCount",//22
                    "nFactor",//23
                    "wFactor",//24
                    "lFactor",//25
                    "gFactor",//26
                    "gwAPAR",//27
                    "gwLWP",//28
                    "gwVPD"//29
                
                    
//                    "day",
//                    "month",
//                    "year",
//                    "basinID",
//                    "hillID",
//                    "zoneID",
//                    "patchID",
//                    "lai",
//                    "plantc",
//                    "plantn",
//                    "net_psn",
//                    "plant_resp",
//                    "soil_resp",
//                    "litr1c",
//                    "litr2c",
//                    "litr3c",
//                    "litr4c",
//                    "litr1n",
//                    "litr2n",
//                    "litr3n",
//                    "litr4n",
//                    "lit.rain_cap",
//                    "soil1c",
//                    "soil2c",
//                    "soil3c",
//                    "soil4c",
//                    "soil1n",
//                    "soil2n",
//                    "soil3n",
//                    "soil4n",
//                    "soilDON",
//                    "soilDOC",
//                    "denitrif",
//                    "netleach",
//                    "DON_loss",
//                    "DOC_loss",
//                    "soilNO3",
//                    "soilNH4",
//                    "streamNO3",
//                    "streamNH4",
//                    "streamDON",
//                    "streamDOC",
//                    "surfaceNO3",
//                    "surfaceNH4",
//                    "surfaceDOC",
//                    "surfaceDON",
//                    "height",
//                    "nuptake",
//                    "root_depth",
//                    "nfix",
//                    "grazingC",
//                    "area",
//                    "cwdc","cwdn",
//                    "gpsn_apar",
//                    "gpsn_tavg",
//                    "gpsn_lwp",
//                    "gpsn_co2",
//                    "gpsn_tmin",
//                    "gpsn_vpd"
                    );
	/*--------------------------------------------------------------*/
	/*	Yearly 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].patch[0].yearly;
	fprintf(outfile, "%s %s %s %s %s %s %s %s %s %s %s %s\n",
		"year",
		"basinID",
		"hillID",
		"zoneID",
		"patchID",
		"litter_c",
		"soil_c",
		"litter_n",
		"soil_n",
		"nitrate",
		"sminn","root_depth");
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
		"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
		"day",
		"month",
		"year",
		"basinID",
		"hillID",
		"zoneID",
		"patchID",
		"stratumID",
		"proj_lai",
		"leafc",
		"leafn",
		"cpool",
		"npool",
		"dead_leafc",
		"frootc",
		"frootn",
		"live_stemc",
		"live_stemn",
		"leafc_store",
		"leafn_store",
		"dead_stemc",
		"dead_stemn",
		"live_crootc",
		"live_crootn",
		"dead_crootc",
		"dead_crootn",
		"cwdc",
		"mresp",
		"gresp",
		"psn_to_cpool","age","root_depth","gwseasonday","lfseasonday","gsi", "nlimit",
		"fleaf","froot","fwood","Nuptake","smin2pl","retrans2pl","mort_fract");

  /*--------------------------------------------------------------*/
	/* Shadow	Daily 			                                   				*/
	/*--------------------------------------------------------------*/
	if (command_line[0].vegspinup_flag > 0){
    outfile = world_output_files[0].shadow_strata[0].daily;
	  fprintf(outfile,
	  	"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
	  	"day",
	  	"month",
	  	"year",
	  	"basinID",
		  "hillID  ",
		  "zoneID  ",
		  "patchID  ",
		  "stratumID  ",
		  "proj_lai  ",
		  "leafc  ",
		  "leafn  ",
		  "cpool  ",
		  "npool  ",
		  "dead_leafc  ",
		  "frootc  ",
		  "frootn  ",
		  "live_stemc  ",
		  "live_stemn  ",
		  "leafc_store  ",
		  "leafn_store  ",
		  "dead_stemc  ",
		  "dead_stemn  ",
		  "live_crootc  ",
		  "live_crootn  ",
		  "dead_crootc  ",
		  "dead_crootn  ",
		  "cwdc  ",
		  "mresp  ",
		  "gresp  ",
		  "psn_to_cpool","age","root_depth","gwseasonday","lfseasonday","gsi", "nlimit",
		  "fleaf","froot","fwood","Nuptake","smin2pl","retrans2pl","mort_fract");
  }

	/*--------------------------------------------------------------*/
	/*	Yearly 							*/
	/*--------------------------------------------------------------*/
	outfile = world_output_files[0].canopy_stratum[0].yearly;
	fprintf(outfile, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
		"year",
		"basinID",
		"hillID",
		"zoneID",
		"patchID",
		"stratumID",
		"proj_lai",
		"leafc",
		"leafn",
		"frootc",
		"frootn",
		"stemc",
		"stemn",
		"cwdc",
		"cwdn",
		"psn","cpool", "mortfract");
	}


	return;
} /*end add_growth_headers*/
