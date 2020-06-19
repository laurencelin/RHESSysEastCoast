/*--------------------------------------------------------------*/
/* 								*/
/*	construct_landuse_defaults				*/
/*								*/
/*	construct_landuse_defaults.c - makes patch default	*/
/*			objects.				*/
/*								*/
/*	NAME							*/
/*	construct_landuse_defaults.c - makes patch default	*/
/*			objects.				*/
/*								*/
/*	SYNOPSIS						*/
/*	struct patch_default *construct_landuse_defaults(       */
/*			num_default_files,			*/
/*			default_files,				*/
/*			grow_flag,				*/
/*			default_object_list )			*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*	Original code, January 15, 1996.			*/
/*	July 28, 1997	C.Tague					*/
/*	removed capillary rise landuse variables 		*/
/*	i.e rooting depth 					*/		
/*	pore size index and suction				*/
/*								*/
/*	Sep 15 97 RAF						*/
/*	added cap rise variables back as we now do a cap rise	*/
/*	routine but removed wilting point which was never used.	*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"
#include "params.h"

struct landuse_default *construct_landuse_defaults(
	int	num_default_files,
	char	**default_files,
	struct command_line_object *command_line) 
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.				*/
	/*--------------------------------------------------------------*/
	void	*alloc(	size_t,
		char	*,
		char	*);
	
	
	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
	int	i;
        int strbufLen = 256;
        int filenameLen = 1024;
	double 	landuse, ftmp;
	FILE	*default_file;
        char	strbuf[strbufLen];
        char	outFilename[filenameLen];
	char	record[MAXSTR];
	char	*newrecord;
	struct 	landuse_default *default_object_list;
        param *paramPtr = NULL;
        int paramCnt = 0;
	
	/*--------------------------------------------------------------*/
	/*	Allocate an array of default objects.			*/
	/*--------------------------------------------------------------*/
	default_object_list = (struct landuse_default *)
		alloc(num_default_files *
		sizeof(struct landuse_default),"default_object_list",
		"construct_landuse_defaults");
	
	/*--------------------------------------------------------------*/
	/*	Loop through the default files list.			*/
	/*--------------------------------------------------------------*/
	for (i=0 ; i<num_default_files; i++){
		/*--------------------------------------------------------------*/
		/*		Try to open the ith default file.		*/
		/*--------------------------------------------------------------*/
		printf("Reading %s\n", default_files[i]);
                paramCnt = 0;
                if (paramPtr != NULL)
                    free(paramPtr);

                paramPtr = readParamFile(&paramCnt, default_files[i]);

		/*--------------------------------------------------------------*/
		/*		read the ith default file into the ith object.	*/
		/*--------------------------------------------------------------*/

		default_object_list[i].ID = 			getIntParam(&paramCnt, &paramPtr, "landuse_default_ID", "%d", 1, 0);
        default_object_list[i].irrigation = 		getDoubleParam(&paramCnt, &paramPtr, "irrigation", "%lf", 0.0, 1)*0.001;// daily max rate (mm) --> (m);
		default_object_list[i].fertilizer_NO3 = 	getDoubleParam(&paramCnt, &paramPtr, "fertilizer_NO3", "%lf", 0.0, 1);//kgN/m2/each time
		default_object_list[i].fertilizer_NH4 = 	getDoubleParam(&paramCnt, &paramPtr, "fertilizer_NH4", "%lf", 0.0, 1);//kgN/m2/each time
        default_object_list[i].fertilizer_freq =     getIntParam(&paramCnt, &paramPtr, "fertilizer_freq", "%d", 30, 1); // # of days between each time
        default_object_list[i].fertilizer_decay_rate = -log(0.1)/(1.0*default_object_list[i].fertilizer_freq);
        // 1.0 - exp(-10/(1.0*default_object_list[i].fertilizer_freq));
        // for period "default_object_list[i].fertilizer_freq", the stored fertilizer should have been gone 10% N0)
        double numPPLinHouse =     getDoubleParam(&paramCnt, &paramPtr, "numPPLinHouse", "%lf", 1.0, 1);
		default_object_list[i].septic_NO3_load = 	getDoubleParam(&paramCnt, &paramPtr, "septic_NO3_load", "%lf", 0.0, 1) / 365.0 * numPPLinHouse; // convert annual kgN/yr to daily
		default_object_list[i].septic_water_load = 	getDoubleParam(&paramCnt, &paramPtr, "septic_water_load", "%lf", 0.0, 1) / 365.0 * numPPLinHouse; // convert annual m3/yr to daily
		default_object_list[i].detention_store_size = 	getDoubleParam(&paramCnt, &paramPtr, "detention_store_size", "%lf", 0.0, 1);
        default_object_list[i].pond_size = getDoubleParam(&paramCnt, &paramPtr, "pondDepth", "%lf", 0.0, 1);
        // ------------ non sense below
		default_object_list[i].PH = 			getDoubleParam(&paramCnt, &paramPtr, "PH", "%lf", 7.0, 1);
		//default_object_list[i].percent_impervious = 	getDoubleParam(&paramCnt, &paramPtr, "landuse.percent_impervious", "%lf", 0.0, 1);
		default_object_list[i].grazing_Closs = 	getDoubleParam(&paramCnt, &paramPtr, "grazing_Closs", "%lf", 0.0, 1) / 365;
        
        // sanitary sewer
        default_object_list[i].sewerDiameter =     getDoubleParam(&paramCnt, &paramPtr, "sewerDiamInch", "%lf", 8.0, 1)*0.0254; //convert to m
        default_object_list[i].sewerDensity =     getDoubleParam(&paramCnt, &paramPtr, "sewerDensityMileAcre", "%lf", 0.0, 1)*1609.344/4046.85642; // convert to m/m2; 0.3976776; related to population density
        default_object_list[i].sewerDepth =     getDoubleParam(&paramCnt, &paramPtr, "sewerDepthft", "%lf", 11, 1)*0.3048;
            // 15ft deepest = 4.572 m; 11ft = 3.3528 m; 4ft = 1.2192 m // assume this is the depth of the bottom of the pipe.
        default_object_list[i].sewer_exfiltrationRate =  getDoubleParam(&paramCnt, &paramPtr, "MaxExf_gpimd", "%lf", 9061.0, 1)* 0.00378541178/0.0254/1609.344*default_object_list[i].sewerDiameter*default_object_list[i].sewerDensity;// m3/m2/d = m/d
        default_object_list[i].sewer_exfiltrationPercent = getDoubleParam(&paramCnt, &paramPtr, "MaxExf_Percent", "%lf", 0.491, 1);
        default_object_list[i].sewer_emptyPercent = getDoubleParam(&paramCnt, &paramPtr, "emptyPercent", "%lf", 0.25, 1);
        default_object_list[i].sewer_infiltrationRate = default_object_list[i].sewer_exfiltrationRate / default_object_list[i].sewer_exfiltrationPercent/(1.0-default_object_list[i].sewer_emptyPercent)*default_object_list[i].sewer_emptyPercent;
        
        //lookup table
        //double upee[] = {1.0,0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05,0.0};
        int jj;
        double emptyPercent[] = {0.0,0.01869304,0.05204402,0.0940602,0.14237849,0.19550111,0.25231579,0.31191883,0.37353004,0.43644429,0.5,0.51869304,0.55204402,0.5940602,0.64237849,0.69550111,0.75231579,0.81191883,0.87353004,0.93644429,1.0};
        for(jj=0; jj<21; jj++){
            if(emptyPercent[jj] >= default_object_list[i].sewer_emptyPercent) break;
        }// end of for jj
        double upee = -(default_object_list[i].sewer_emptyPercent-emptyPercent[jj-1])/(emptyPercent[jj]-emptyPercent[jj-1])*0.05+(1-0.05*(jj-1));
        default_object_list[i].sewer_infiltrationSatDefZThreshold = default_object_list[i].sewerDepth - default_object_list[i].sewerDiameter*upee;// this is the depth of the empty space of the pipe, direct comparable to sat_def_z
        default_object_list[i].sewer_infiltrationSatDefZHeadSpace = 1.0/(default_object_list[i].sewer_infiltrationSatDefZThreshold - (default_object_list[i].sewerDepth-default_object_list[i].sewerDiameter)); // 1/"thickness of empty space"
        // default_object_list[i].sewerDepth-default_object_list[i].sewerDiameter = top depth of sewer pipe
        default_object_list[i].sewer_exfiltrationDepthVol = 1.0/(default_object_list[i].sewerDepth-default_object_list[i].sewer_infiltrationSatDefZThreshold);
        
        
        default_object_list[i].sewerNO3c = getDoubleParam(&paramCnt, &paramPtr, "sewerNO3mgNL", "%lf", 2.5, 1)*0.001; // convert to kgN/m3
        default_object_list[i].sewerNH4c = getDoubleParam(&paramCnt, &paramPtr, "sewerNH4mgNL", "%lf", 2.5, 1)*0.001;
        default_object_list[i].sewerDONc = getDoubleParam(&paramCnt, &paramPtr, "sewerDONmgNL", "%lf", 2.5, 1)*0.001;
        default_object_list[i].sewerDOCc = getDoubleParam(&paramCnt, &paramPtr, "sewerDOCmgCL", "%lf", 2.5, 1)*0.001;
        /*--------------------------------------------------------------*/
		/*		Close the ith default file.								*/
		/*--------------------------------------------------------------*/

                memset(strbuf, '\0', strbufLen);
                strcpy(strbuf, default_files[i]);
                char *s = strbuf;
                char *y = NULL;
                char *token = NULL;
                char filename[256];
    
                // Store filename portion of path in 't'
                while ((token = strtok(s, "/")) != NULL) {
                    // Save the latest component of the filename
                    strcpy(filename, token);
                    s = NULL;
                } 
    
                // Remove the file extension, if one exists
                memset(strbuf, '\0', strbufLen);
                strcpy(strbuf, filename);
                free(s);
                s = strbuf;
                token = strtok(s, ".");
                if (token != NULL) {
                    strcpy(filename, token);
                }
        
                memset(outFilename, '\0', filenameLen);


        
    
            // Concatenate the output prefix with the filename of the input .def file
            // and "_landuse.params"
            if (command_line[0].output_prefix != NULL) {
                strcat(outFilename, command_line[0].output_prefix);
                if (filename != NULL) {
                    strcat(outFilename, "_");
                    strcat(outFilename, filename);
                }
                strcat(outFilename, "_landuse.params");
            } 
            else {
                if (filename != NULL) {
                    strcat(outFilename, "_");
                    strcat(outFilename, filename);
                }
                strcat(outFilename, "landuse.params");
            }
    
                printParams(paramCnt, paramPtr, outFilename);
	} /*end for*/
        free(paramPtr);
	return(default_object_list);
} /*end construct_landuse_defaults*/
