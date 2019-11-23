/*--------------------------------------------------------------*/
/* 																*/
/*					construct_command_line						*/
/*																*/
/*	construct_command_line.c - creates command line object		*/
/*																*/
/*	NAME														*/
/*	construct_command_line.c - creates command line object		*/
/*																*/
/*	SYNOPSIS													*/
/*	struct	command_line_object	*construct_command_line(		*/
/*								argc, argv, command_line)		*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	Creates a command_line object which consists of flags		*/
/*	entered on the command line during execution of rhessys.	*/
/*	Some error checking is performed but error checking must	*/
/*	wait until the world object has been specified.				*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	Original code, January 15, 1996.							*/
/*	valid_option to be written still - determines if the next	*/
/*			arguement is a valid option.						*/
/*	added routing option - May 7, 1997, C.Tague					*/
/*																*/
/*																*/
/*	Sep 1997	RAF												*/
/*	Removed extended output option flag as all output is		*/
/*	now of a single format specified by output routines.		*/
/*								*/
/*	Sept 1998	C.Tague					*/
/* 	added comma separated output option			*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rhessys.h"
struct	command_line_object	*construct_command_line(
													int main_argc,
													char **main_argv)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	int	valid_option( char * );
	void	*alloc( size_t, char *, char * );
	void	output_template_structure();
	
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int		i;
	struct	command_line_object	*command_line;
	
	/*--------------------------------------------------------------*/
	/*	Allocate a command line object.								*/
	/*--------------------------------------------------------------*/
	command_line = (struct command_line_object *)
		alloc(1 * sizeof(struct command_line_object),
		"command_line","construct_command_line");
	
	/*--------------------------------------------------------------*/
	/*	Initialize the options as null				*/
	/*--------------------------------------------------------------*/
	command_line[0].gridded_ascii_flag = 0;
	command_line[0].gridded_netcdf_flag = 0;
	command_line[0].grow_flag = 0;	
	command_line[0].std_scale = 0;
	command_line[0].stdev_flag = 0;
	command_line[0].road_flag = 1;
	command_line[0].prefix_flag = 0;
	command_line[0].verbose_flag = 0;
	command_line[0].routing_flag = 0;
	command_line[0].surface_routing_flag = 0;
	command_line[0].stream_routing_flag = 0;
	command_line[0].reservoir_operation_flag = 0;
	command_line[0].dclim_flag = 0;
	command_line[0].ddn_routing_flag = 0;
	command_line[0].tec_flag = 0;
	command_line[0].world_flag = 0;
	command_line[0].world_header_flag = 0;
	command_line[0].start_flag = 0;
	command_line[0].end_flag = 0;
	command_line[0].sen_flag = 0;
	command_line[0].vsen_flag = 0;
	command_line[0].vsen_alt_flag = 0;
	command_line[0].precip_scale_flag = 0;
	command_line[0].snow_scale_flag = 0;
	command_line[0].noredist_flag = 0;
	command_line[0].surface_energy_flag = 0;
	command_line[0].firespread_flag = 0;
	command_line[0].vegspinup_flag = 0;		
	command_line[0].vgsen_flag = 0;
	command_line[0].FillSpill_flag=0;	
	command_line[0].evap_use_longwave_flag = 0;
	command_line[0].veg_sen1 = 1.0;
	command_line[0].veg_sen2 = 1.0;
	command_line[0].veg_sen3 = 1.0;
	command_line[0].vmort_flag = 0;
	command_line[0].version_flag = 0;
	command_line[0].vsen[M] = 1.0;
	command_line[0].vsen[K] = 1.0;
    command_line[0].psen[M] = 1.0;
    command_line[0].psen[K] = 1.0;
	command_line[0].sen[M] = 1.0;
	command_line[0].sen[K] = 1.0;
	command_line[0].sen[SOIL_DEPTH] = 1.0;
	command_line[0].prev_flag = 0;
	command_line[0].gw_flag = 0;
	command_line[0].tchange_flag = 0;
	command_line[0].tmax_add = 0.0;
	command_line[0].tmin_add = 0.0;
	command_line[0].output_prefix = NULL;
	command_line[0].output_flags.yearly = 0;
	command_line[0].output_flags.monthly = 0;
	command_line[0].output_flags.daily = 0;
	command_line[0].output_flags.hourly = 0;
	command_line[0].stro = NULL;
	command_line[0].b = NULL;
	command_line[0].h = NULL;
	command_line[0].z = NULL;
	command_line[0].p = NULL;
	command_line[0].c = NULL;
	command_line[0].tmp_value = 1.0;
	command_line[0].thresholds[SATDEF] = 0.0;
	command_line[0].thresholds[STREAMFLOW] = 0.0;
	command_line[0].snow_scale_tol = 999999999;
    
    command_line[0].gwtoriparian_flag = 0;
    command_line[0].gwtoriparian_ID = 400; // any LU ID between [gwtoriparian_ID, gwtoriparian_ID+100) = ripiarn
    
    // mostly added by Laurence
	command_line[0].capreduction = 1.0; // scaler to cap rise
    command_line[0].caprsplit = 1.0; // % to rootzone and unsat below root zont
    command_line[0].capMax = 0.0;// disable at 0.0; upper bound fpr cap rise
    command_line[0].newcaprise_flag = 0; // new cap rise mechanism, passing unsat first then root
    command_line[0].slowDrain_flag = 0; // reduce drainage
    
    command_line[0].snowT_scaler = 1.0;// snowT scaler
    command_line[0].snowE_scaler = 1.0;
    command_line[0].toptoff_flag = 0; // switch off "topt scaling" for plant gl calculation
    command_line[0].rootdepthz = 1.0; // scaler for root depth
    command_line[0].dynRtZoff_flag =0;// dwitch off dynRT depth calcuation in runtime
    command_line[0].iniBioRZ = 0; // set initial root depth using "update_rooting_depth.c"
    
    command_line[0].fracDirectNdep = 0.0; // % directly hit dep N hit the ground, passing canopy
    command_line[0].soilDecayScalar = 1.0;// scale to litter/soilc decay rate base values
    command_line[0].BGC_flag = 0; // BGC decomposition flag
    command_line[0].soilCNadaptation_flag = 0; // current soil4CN is fixed <<---- still confusing!
    
    command_line[0].rootNdecayRate = 0;// (flag) <<---- still confusing!
    command_line[0].root2active = -1.0; // scale (x) <<---- still confusing!
    command_line[0].NH4root2active = -1.0; // scale (x); <<---- still confusing!
        /*
                            root2active(x)      NH4root2active(x)        rootNdecayRate
            nitrif_decay    N_decay_rate        N_decay_rate             NH4decayRate       N_decay_rate
            nitrif_depth    x*rootdepth         x*rootdepth              soildepth          active_zone_z
         
            denitrif_decay  N_decay_rate        N_decay_rate             NO3decayRate       N_decay_rate
            denitrif_depth  x*rootdepth         x*rootdepth              soildepth          active_zone_z
         
                                                                         DOMdecayRate
            leaching_depth  x*rootdepth         x*rootdepth              soildepth          active_zone_z
            decomp/immobilization/uptake
         */
    command_line[0].Rsolute2gw = 1.0; // reduce solute directly entering GW from surface
    command_line[0].soluteLoss2GW = 0; // subsurface solute to GW without return to surface
    
    command_line[0].leafDarkRespScalar = 1.0; // scaler
    command_line[0].frootRespScalar = 1.0; // scaler
    command_line[0].StemWoodRespScalar = 1.0; // scaler
    command_line[0].aggregate_flag = 0; // make aggreated result based on a special flow table.
    command_line[0].patchPrintTh = 0; // any LU ID between [patchPrintTh, patchPrintTh+100) = print out (patch output)
    
    command_line[0].grassIrrigation_flag=0;
    command_line[0].fertilizer_flag=0;
    command_line[0].sewer_flag=0;
    command_line[0].stormDrainFrac = 0.0; //<------ becomes no use Spet 12
    command_line[0].readinWFdoc_flag=0;
	/*-------------------------------------------------*/
	/* Loop through each arguement in the command line.*/
	/*-------------------------------------------------*/
	i = 1;
	while  ( i < main_argc){
		/*------------------------------------------*/
		/* Check for the print version flag         */
		/*------------------------------------------*/
		if ( strcmp(main_argv[i], "-version") == 0) {
			command_line[0].version_flag = 1;
			++i;
		}
		/*------------------------------------------*/
		/*	Check if the verbose flag is next.  */
		/*------------------------------------------*/
		if ( i< main_argc ){
			if ( strcmp(main_argv[i],"-v") == 0 ){
				/*-----------------------------------------------*/
				/*			Check if "-v" was the last agruement.   */
				/*-----------------------------------------------*/
				i++;
				if ( i == main_argc ){
					/*------------------------------------------*/
					/* assume we want verbose level 1			  */
					/*------------------------------------------*/
					command_line[0].verbose_flag= 1;
				}
				else if ( valid_option(main_argv[i]) == 1 ){
					/*----------------------------------------------*/
					/*	check if the next arguement is an option.		*/
					/*----------------------------------------------*/
					command_line[0].verbose_flag= 1;
				}
				else{
					/*-------------------------------------------------*/
					/*	read in the value of the verbose level.			*/
					/*-------------------------------------------------*/
					command_line[0].verbose_flag = (int)atoi(main_argv[i]);
					i++;
				}/*end if*/
			}/*end if*/
			/*------------------------------------------*/
			/* output_template file - this will terminate RHESSys RUN */
			/*------------------------------------------*/
			else if (strcmp(main_argv[i],"-template") == 0) {
				printf("\n Outputting template file structure and Exiting\n");
				output_template_structure();
				exit(EXIT_FAILURE);
			}
			/*------------------------------------------*/
			/*Check if the no redistribution flag is next.           */
			/*------------------------------------------*/
			else if ( strcmp(main_argv[i],"-noredist") == 0 ){
				printf("\n Running with no lateral redistribution - water balance not maintained ");
				command_line[0].noredist_flag = 1;
				i++;
			}
			/*------------------------------------------*/
			/*Check if the variable mortality flag is next.           */
			/*------------------------------------------*/
			else if ( strcmp(main_argv[i],"-vmort") == 0 ){
				command_line[0].vmort_flag = 1;
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Value for vmort flag not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*-------------------------------*/
				/*Read in the tmp value		*/
				/*-------------------------------*/
				command_line[0].cpool_mort_fract = (double)atof(main_argv[i]);
				i++;
			}
			/*------------------------------------------*/
			/*Check if the distributed climate flag is next.           */
			/*------------------------------------------*/
			else if ( strcmp(main_argv[i],"-dclim") == 0 ){
				command_line[0].dclim_flag = 1;
				i++;
			}
			/*------------------------------------------*/
			/*Check if the grow flag is next.           */
			/*------------------------------------------*/
			else if ( strcmp(main_argv[i],"-g") == 0 ){
				command_line[0].grow_flag = 1;
				i++;
			}
			/*-------------------------------------------------*/
			/*Check tmp value option is next.				*/
			/* Currently tmp is used for sensitivity analysis of rooting_depth  */ 
			/*-------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-tmp") == 0 ){
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Value for Tmp variable not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*-------------------------------*/
				/*Read in the tmp value		*/
				/*-------------------------------*/
				command_line[0].tmp_value = (double)atof(main_argv[i]);
				i++;
			}/* end if */


			/*-------------------------------------------------------*/
			/*Check if the snow distribution flag is next.           */ //<<-----------------------------
			/*-------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-snowdistb") == 0 ){
				printf("\n Running wiith snow redistribution");
				command_line[0].snow_scale_flag = 1;
				i++;
				/*--------------------------------------------------------------*/
				/*	check to see if there is a tolerance parameter 		*/
				/*--------------------------------------------------------------*/
				if (  (i != main_argc) && (valid_option(main_argv[i])==0) ){
					command_line[0].snow_scale_tol = (double)atof(main_argv[i]);
					i++;
				}/*end if*/
			}/* end if */
			/*-------------------------------------------------*/
			/*	fire spread option and coeffcients	  */                //<<-----------------------------
			/*-------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-firespread") == 0 ){
				printf("\n Running with FIRE SPREAD turned on");
				command_line[0].firespread_flag = 1;
				i++;
				command_line[0].fire_grid_res = 30;
				/*-------------------------------*/
				/*Read in the fire spread grid parameters		*/
				/*-------------------------------*/
				if (  (i != main_argc) && (valid_option(main_argv[i])==0) ){
					command_line[0].fire_grid_res = (double)atof(main_argv[i]);
					i++;
				}/*end if*/
			}/* end if */

			/*-------------------------------------------------*/
			/*	surface energy option */
			/*-------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-surfaceenergy") == 0 ){
				i++;
				printf("\n Running with SURFACE ENERGY CALC turned on");
				command_line[0].surface_energy_flag = 1;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Values for fire grid parameters not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
			}/* end if */
            /*------------------------------------------*/
            /*Check if the aggregate_flag flag is next.           */
            /*------------------------------------------*/
            else if ( strcmp(main_argv[i],"-aggregate_flag") == 0 ){
                printf("active aggregate_flag\n");
                command_line[0].aggregate_flag = 1;
                i++;
            }
            /*------------------------------------------*/
            /*Check if the grassIrrigation_flag flag is next.           */
            /*------------------------------------------*/
            else if ( strcmp(main_argv[i],"-grassIrrigation_flag") == 0 ){
                printf("active grassIrrigation_flag\n");
                command_line[0].grassIrrigation_flag = 1;
                i++;
            }
            /*------------------------------------------*/
            /*Check if the fertilizer_flag flag is next.           */
            /*------------------------------------------*/
            else if ( strcmp(main_argv[i],"-fertilizer_flag") == 0 ){
                printf("active fertilizer_flag\n");
                command_line[0].fertilizer_flag = 1;
                i++;
            }
            /*------------------------------------------*/
            /*Check if the sewer_flag flag is next.           */
            /*------------------------------------------*/
            else if ( strcmp(main_argv[i],"-sewer_flag") == 0 ){
                printf("active sewer_flag\n");
                command_line[0].sewer_flag = 1;
                i++;
            }
            /*------------------------------------------*/
            /*Check if the readinWFdoc_flag flag is next.           */
            /*------------------------------------------*/
            else if ( strcmp(main_argv[i],"-readinWFdoc_flag") == 0 ){
                printf("active readinWFdoc_flag\n");
                command_line[0].readinWFdoc_flag = 1;
                i++;
            }
            /*-------------------------------------------------*/
            /*	new cap rise option */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-newcaprise") == 0 ){
                printf("\n Running with newcaprise turned on\n");
                command_line[0].newcaprise_flag = 1;
                i++;
            }/* end if */
            
            /*-------------------------------------------------*/
            /*    rootNdecayRate option */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-rootNdecayRate") == 0 ){
                printf("\n Running with rootNdecayRate turned on\n");
                command_line[0].rootNdecayRate = 1;
                i++;
            }/* end if */
           
            /*-------------------------------------------------*/
            /*    BGC_flag option */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-BGC_flag") == 0 ){
                printf("\n Running with BGC_flag turned on\n");
                command_line[0].BGC_flag = 1;
                i++;
            }/* end if */
            
            /*-------------------------------------------------*/
            /*    soilCNadaptation_flag option */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-soilCNadaptation_flag") == 0 ){
                printf("\n Running with soilCNadaptation_flag turned on\n");
                command_line[0].soilCNadaptation_flag = 1;
                i++;
            }/* end if */
            
            /*-------------------------------------------------*/
            /*    new unsat drain option */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-slowDrain") == 0 ){
                printf("\n Running with slowdrain turned on\n");
                command_line[0].slowDrain_flag = 1;
                i++;
            }/* end if */
            
            /*-------------------------------------------------*/
            /*    iniBioRZ option */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-iniBioRZ") == 0 ){
                printf("\n Running with iniBioRZ\n");
                command_line[0].iniBioRZ = 1;
                i++;
            }/* end if */
            
			/*------------------------------------------*/
			/*Check if spinup flag next.                */
			/*------------------------------------------*/
			else if ( strcmp(main_argv[i],"-vegspinup") == 0 ){
				printf("\n Running with SPINUP turned on \n");
				command_line[0].vegspinup_flag = 1;
				i++;
     
      	/*--------------------------------------------------------------*/
				/*			Read in the vegspinup file name.						          	*/
				/*--------------------------------------------------------------*/
				strncpy(command_line[0].vegspinup_filename, main_argv[i], FILEPATH_LEN);
				i++;
      }

			/*-------------------------------------------------*/
			/*	routing gw to riparian option */
			/*-------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-gwtoriparian") == 0 ){
				i++;
				printf("\n Running with hillslope gw routed to riparian areas\n ");
				command_line[0].gwtoriparian_flag = 1;
			}/* end if */
            /*-------------------------------------------------*/
            /*    routing gw to riparian IDs */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-gwtoriparian_ID") == 0 ){
                i++;
                command_line[0].gwtoriparian_flag = 1;
                if ( (i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for gwtoriparian_ID\n");
                    exit(EXIT_FAILURE);
                } /*end if*/
                
                /*-------------------------------*/
                /*Read in the gwtoriparian_ID values        */
                /*-------------------------------*/
                command_line[0].gwtoriparian_ID = (int)atoi(main_argv[i]);
                i++;
            }/* end if */
			/*-------------------------------------------------*/
			/*	groundwater flag and coeffcients	  */
			/*-------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-gw") == 0 ){
				i++;
				command_line[0].gw_flag = 1;
				if ( (i == main_argc) || (i == main_argc-1) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Values for gw coefficients not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*-------------------------------*/
				/*Read in the loss to gw rate multiplier values		*/
				/*-------------------------------*/
				command_line[0].sat_to_gw_coeff_mult = (double)atof(main_argv[i]);
				i++;
				command_line[0].gw_loss_coeff_mult = (double)atof(main_argv[i]);
				i++;
			}/* end if */
            /*-------------------------------------------------*/
            /*    scaler to adjust fracDirectNdep      */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-fracDirectNdep") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for fracDirectNdep\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values        */
                /*-------------------------------*/
                command_line[0].fracDirectNdep = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*    scaler to adjust stormDrainFrac      */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-stormDrainFrac") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for stormDrainFrac\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values        */
                /*-------------------------------*/
                command_line[0].stormDrainFrac = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*    scaler to adjust leafDarkRespScalar      */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-leafDarkRespScalar") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for leafDarkRespScalar\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values        */
                /*-------------------------------*/
                command_line[0].leafDarkRespScalar = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*    scaler to adjust frootRespScalar      */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-frootRespScalar") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for frootRespScalar\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values        */
                /*-------------------------------*/
                command_line[0].frootRespScalar = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*    scaler to adjust frootRespScalar      */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-StemWoodRespScalar") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for StemWoodRespScalar\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values        */
                /*-------------------------------*/
                command_line[0].StemWoodRespScalar = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*    scaler to adjust NH4root2active      */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-NH4root2active") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for NH4root2active\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values        */
                /*-------------------------------*/
                command_line[0].NH4root2active = (double)atof(main_argv[i]);
                i++;
                if(command_line[0].NH4root2active>0.0) printf("NH4root2active is triggered. \n");
                
            }/* end if */
            /*-------------------------------------------------*/
            /*    scaler to adjust soluteLoss2GW      */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-soluteLoss2GW") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for NH4root2active\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values        */
                /*-------------------------------*/
                command_line[0].soluteLoss2GW = (double)atof(main_argv[i]);
                i++;
                if(command_line[0].soluteLoss2GW>0.0) printf("soluteLoss2GW is triggered. \n");
                
            }/* end if */
            /*-------------------------------------------------*/
            /*    scaler to adjust root2active      */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-root2active") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for root2active\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values        */
                /*-------------------------------*/
                command_line[0].root2active = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*    scaler to adjust soilDecayScalar      */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-soilDecayScalar") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for soilDecayScalar\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values        */
                /*-------------------------------*/
                command_line[0].soilDecayScalar = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*	scaler to adjust snow melt due to temperature	  */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-snowTs") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for snowTs not specified\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values		*/
                /*-------------------------------*/
                command_line[0].snowT_scaler = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*    scaler to adjust snow melt due to energy      */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-snowEs") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for snowEs not specified\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values        */
                /*-------------------------------*/
                command_line[0].snowE_scaler = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*    scaler to Rsolute2gw      */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-Rsolute2gw") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for Rsolute2gw not specified\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                printf("Rsolute2gw is triggered.\n");
                /*-------------------------------*/
                /*Read in the reduction multiplier values        */
                /*-------------------------------*/
                command_line[0].Rsolute2gw = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*	patchPrintTh	  */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-patchPrintTh") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for print patch threshold not specified\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values		*/
                /*-------------------------------*/
                command_line[0].patchPrintTh = (int)atoi(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*	cap rise reduction flag and coeffcients	  */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-capr") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for cap rise reduction coefficients not specified\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values		*/
                /*-------------------------------*/
                command_line[0].capreduction = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*	cap rise split flag and coeffcients	  */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-caprsplit") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for cap rise split not specified\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values		*/
                /*-------------------------------*/
                command_line[0].caprsplit = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*	set capMax	  */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-capMax") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for capMax not specified\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values		*/
                /*-------------------------------*/
                command_line[0].capMax = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*	root depth	  */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-rtz") == 0 ){
                
                i++;
                if ((i == main_argc) || (valid_option(main_argv[i])==1)){
                    fprintf(stderr,"FATAL ERROR: Values for rtz not specified\n");
                    exit(EXIT_FAILURE);
                    printf("Here\n");
                } /*end if*/
                /*-------------------------------*/
                /*Read in the reduction multiplier values		*/
                /*-------------------------------*/
                command_line[0].rootdepthz = (double)atof(main_argv[i]);
                i++;
                
            }/* end if */
            /*-------------------------------------------------*/
            /*	turn off glmax response to tavg */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-dynRtZoff") == 0 ){
                i++;
                printf("\n Running without dynamic root depth\n");
                command_line[0].dynRtZoff_flag = 1;
            }/* end if */
            /*-------------------------------------------------*/
            /*	turn off glmax response to tavg */
            /*-------------------------------------------------*/
            else if ( strcmp(main_argv[i],"-toptoff") == 0 ){
                i++;
                printf("\n Running without topt control on stomatal closure\n");
                command_line[0].toptoff_flag = 1;
            }/* end if */
			/*-------------------------------------------------*/
			/* simple addition of temperature increases 		*/
			/*-------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-tchange") == 0 ){
				i++;
				command_line[0].tchange_flag = 1;
				if ( (i == main_argc) || (i == main_argc-1) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Values for tchange not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*-------------------------------*/
				/*Read in the loss to gw rate multiplier values		*/
				/*-------------------------------*/
				command_line[0].tmax_add = (double)atof(main_argv[i]);
				i++;
				command_line[0].tmin_add = (double)atof(main_argv[i]);
				i++;
			}/* end if */
			/*-------------------------------------------------*/
			/* 	soil moisture standard deviation	-std */
			/*	if this flag is set there must be an extra variable in the worldfile */
			/*	at the patch level which inputs a std for that patch 	*/	
			/*-------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-stdev") == 0 ){
				i++;
				command_line[0].stdev_flag = 1;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Value for soil moisture std not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				command_line[0].std_scale = (double)atof(main_argv[i]);
				i++;
			}/* end if */
			/*-------------------------------------------------*/
			/*Check if the threshold option is next.				*/
			/*-------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-th") == 0 ){
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Thresholds not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*-------------------------------*/
				/*Read in the thresholds 	*/
				/* 	at present there are up to two values - threshold */
				/*	sat_deficit - used for outputting water stress days 	*/
				/*	and streamflow - used for outputting low flow days	*/
				/*-------------------------------*/
				command_line[0].thresholds[SATDEF] = (double)atof(main_argv[i]);
				i++;
				/*--------------------------------------------------------------*/
				/*	check to see if there is a second threshold parameter 	*/
				/*--------------------------------------------------------------*/
				if (  (i != main_argc) && (valid_option(main_argv[i])==0) ){
					command_line[0].thresholds[STREAMFLOW] = (double)atof(main_argv[i]);
					i++;
				}/*end if*/
			}/* end if */
			/*-----------------------------------------------------------*/
			/*Check if the sensitivity analysis option is next.				*/
			/*-----------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-s") == 0 ){
				i++;
				command_line[0].sen_flag = 1;
				if (  (i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,
						"FATAL ERROR: Sensitivity perturbation not specified\n");
					exit(EXIT_FAILURE);
				}/*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the sensitivity parameter values. */
				/*--------------------------------------------------------------*/
				command_line[0].sen[M] = (double)atof(main_argv[i]);
				i++;
				if (  (i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,
						"FATAL ERROR: Sensitivity perturbation not specified\n");
					exit(EXIT_FAILURE);
				}/*end if*/
				command_line[0].sen[K] = (double)atof(main_argv[i]);
				i++;
				/*--------------------------------------------------------------*/
				/*	check to see if there is a 3rd sensitivity parameter 	*/
				/*	if not set to 1.0					*/
				/*--------------------------------------------------------------*/
				if (  (i != main_argc) && (valid_option(main_argv[i]) == 0) ){
					command_line[0].sen[SOIL_DEPTH] = (double)atof(main_argv[i]);
					i++;
				}  /*end if*/
			} /* end if */


			/*-----------------------------------------------------------*/
			/*Check if the vertical  sensitivity analysis option is next.				*/
			/*-----------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-sv") == 0 ){
				i++;
				command_line[0].vsen_flag = 1;
				if (  (i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,
						"FATAL ERROR: Sensitivity perturbation not specified\n");
					exit(EXIT_FAILURE);
				}/*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the sensitivity parameter values. */
				/*--------------------------------------------------------------*/
				command_line[0].vsen[M] = (double)atof(main_argv[i]);
				i++;
				if (  (i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,
						"FATAL ERROR: Sensitivity perturbation not specified\n");
					exit(EXIT_FAILURE);
				}/*end if*/
				command_line[0].vsen[K] = (double)atof(main_argv[i]);
				i++;
			} /* end if */
            
            else if ( strcmp(main_argv[i],"-spor") == 0 ){
                i++;
                if (  (i == main_argc) || (valid_option(main_argv[i])==1) ){
                    fprintf(stderr,
                        "FATAL ERROR: Sensitivity perturbation not specified\n");
                    exit(EXIT_FAILURE);
                }/*end if*/
                /*--------------------------------------------------------------*/
                /*            Read in the sensitivity parameter values. */
                /*--------------------------------------------------------------*/
                command_line[0].psen[M] = (double)atof(main_argv[i]);
                i++;
                if (  (i == main_argc) || (valid_option(main_argv[i])==1) ){
                    fprintf(stderr,
                        "FATAL ERROR: Sensitivity perturbation not specified\n");
                    exit(EXIT_FAILURE);
                }/*end if*/
                command_line[0].psen[K] = (double)atof(main_argv[i]);
                i++;
            } /* end if */

			/*-------------------------------------------------------*/
			/*Check if the precip scaling (using random dist) flag is next.           */
			/*-------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-precip") == 0 ){
				printf("\n Running wiith stochastic precipitation scaling ");
				command_line[0].precip_scale_flag = 1;
				i++;

			} /* end if */

			/*-----------------------------------------------------------*/
			/*alternatively use pore size inidex and psi air entry				*/
			/*-----------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-svalt") == 0 ){
				i++;
				command_line[0].vsen_alt_flag = 1;
				if (  (i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,
						"FATAL ERROR: Sensitivity perturbation not specified\n");
					exit(EXIT_FAILURE);
				}/*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the sensitivity parameter values. */
				/*--------------------------------------------------------------*/
				command_line[0].vsen_alt[PA] = (double)atof(main_argv[i]);
				i++;
				if (  (i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,
						"FATAL ERROR: Sensitivity perturbation not specified\n");
					exit(EXIT_FAILURE);
				}/*end if*/
				command_line[0].vsen_alt[PO] = (double)atof(main_argv[i]);
				i++;
			} /* end if */

			/*-----------------------------------------------------------*/
			/*Check if the vegetation  sensitivity analysis option is next.				*/
			/*-----------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-vgsen") == 0 ){
				i++;
				if (  (i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,
						"FATAL ERROR: Sensitivity perturbation not specified\n");
					exit(EXIT_FAILURE);
				}/*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the sensitivity parameter name. */
				/*--------------------------------------------------------------*/
				command_line[0].vgsen_flag = 1;
				command_line[0].veg_sen1 = (double)atof(main_argv[i]);
				i++;
				if (  (i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,
						"FATAL ERROR: Sensitivity perturbation not specified\n");
					exit(EXIT_FAILURE);
				}/*end if*/
				command_line[0].veg_sen2 = (double)atof(main_argv[i]);
				i++;
				if (  (i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,
						"FATAL ERROR: 3rd Vegetation Sensitivity perturbation not specified\n");
					exit(EXIT_FAILURE);
				}/*end if*/
				command_line[0].veg_sen3 = (double)atof(main_argv[i]);
				i++;
			} /* end if */
			/*--------------------------------------------------------------*/
			/* check for start and end dates				*/
			/*--------------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-st") == 0 ){
				command_line[0].start_flag = 1;
				/*-------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*-------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Start date year not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*-------------------------------------------------*/
				/*			Read in the start year.							*/
				/*-------------------------------------------------*/
				command_line[0].start_date.year = (int)atoi(main_argv[i]);
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Start date month not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the start month.							*/
				/*--------------------------------------------------------------*/
				command_line[0].start_date.month = (int)atoi(main_argv[i]);
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if (  (i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Start date day not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the start day.							*/
				/*--------------------------------------------------------------*/
				command_line[0].start_date.day = (int)atoi(main_argv[i]);
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Start date hour not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the hour year.							*/
				/*--------------------------------------------------------------*/
				command_line[0].start_date.hour = (int)atoi(main_argv[i]);
				i++;
			}/*end if*/
			/*--------------------------------------------------------------*/
			/* check for start and end dates				*/
			/*--------------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-ed") == 0 ){
				command_line[0].end_flag = 1;
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Start date year not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the end year.							*/
				/*--------------------------------------------------------------*/
				command_line[0].end_date.year = (int)atoi(main_argv[i]);
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Start date month not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the end month.							*/
				/*--------------------------------------------------------------*/
				command_line[0].end_date.month = (int)atoi(main_argv[i]);
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Start date day not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the end day.							*/
				/*--------------------------------------------------------------*/
				command_line[0].end_date.day = (int)atoi(main_argv[i]);
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Start date hour not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the hour year.							*/
				/*--------------------------------------------------------------*/
				command_line[0].end_date.hour = (int)atoi(main_argv[i]);
				i++;
			}/*end if*/
			/*--------------------------------------------------------------*/
			/*		Check if the routing option file is next.				*/
			/*--------------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-r") == 0 ){
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,"FATAL ERROR: Routing file name not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the routing file name.						*/
				/*--------------------------------------------------------------*/
				command_line[0].routing_flag = 1;
				strncpy(command_line[0].routing_filename, main_argv[i], FILEPATH_LEN);
				i++;
				/*--------------------------------------------------------------*/
				/*	Attempt to read in surface routing file name.				*/
				/*--------------------------------------------------------------*/
				if ( (i < main_argc) && !valid_option(main_argv[i]) ) {
					command_line[0].surface_routing_flag = 1;
					strncpy(command_line[0].surface_routing_filename, main_argv[i], FILEPATH_LEN);
					i++;
				}
			} /*end if*/


			/*--------------------------------------------------------------*/
			/*		Check if the stream routing option file is next.				*/
			/*--------------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-str") == 0 ){
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,"FATAL ERROR: Routing file name not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the routing file name.						*/
				/*--------------------------------------------------------------*/
				command_line[0].stream_routing_flag = 1;
				strcpy(command_line[0].stream_routing_filename,main_argv[i]);
				i++;
			} /*end if*/


			/*--------------------------------------------------------------*/
			/*		Check if the reservoir option file is next.				*/
			/*--------------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-res") == 0 ){
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,"FATAL ERROR: Routing file name not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the reservoir file name.						*/
				/*--------------------------------------------------------------*/
				command_line[0].reservoir_operation_flag = 1;
				strcpy(command_line[0].reservoir_operation_filename,main_argv[i]);
				i++;
			} /*end if*/


			/*--------------------------------------------------------------*/
			/*		Check if the ddn routing option file is next.				*/
			/*--------------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-rddn") == 0 ){
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,"FATAL ERROR: Routing file name not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the routing file name.						*/
				/*--------------------------------------------------------------*/
				command_line[0].ddn_routing_flag = 1;
				command_line[0].routing_flag = 1;
				strcpy(command_line[0].routing_filename,main_argv[i]);
				i++;
			} /*end if*/

			/*--------------------------------------------------------------*/
			/*		Check if the world file is next.						*/
			/*--------------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-w") == 0 ){
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,"FATAL ERROR: World file name not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the world file name.						*/
				/*--------------------------------------------------------------*/
				command_line[0].world_flag = 1;
				strcpy(command_line[0].world_filename,main_argv[i]);
				i++;
			} /*end if*/
			/*--------------------------------------------------------------*/
			/*		Check if the world header file is next.						*/
			/*--------------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-whdr") == 0 ){
				/*--------------------------------------------------------------*/
				/*			Check that the next argument exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,"FATAL ERROR: World file header name not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the world file name.						*/
				/*--------------------------------------------------------------*/
				command_line[0].world_header_flag = 1;
				strcpy(command_line->world_header_filename,main_argv[i]);
				i++;
			} /*end if*/
			/*--------------------------------------------------------------*/
			/*		Check if the tec file is next.							*/
			/*--------------------------------------------------------------*/
			else if( strcmp(main_argv[i],"-t") == 0 ){
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1) ){
					fprintf(stderr,"FATAL ERROR: TEC file name not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Read in the tec file name.							*/
				/*--------------------------------------------------------------*/
				command_line[0].tec_flag = 1;
				strcpy(command_line[0].tec_filename,main_argv[i]);
				i++;
			} /*end if*/
			/*--------------------------------------------------------------*/
			/*		Check if the output prefix is next.						*/
			/*--------------------------------------------------------------*/
			else if( strcmp(main_argv[i],"-pre") == 0 ){
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if ((i == main_argc) || (valid_option(main_argv[i])==1)){
					fprintf(stderr,"FATAL ERROR: Output prefix not specified\n");
					exit(EXIT_FAILURE);
				} /*end if*/
				/*--------------------------------------------------------------*/
				/*			Allocate an array for the output prefix and			*/
				/*			Read in the output prefix .							*/
				/*--------------------------------------------------------------*/
				command_line[0].output_prefix =
					(char *) alloc((1+strlen(main_argv[i]))*sizeof(char),
					"output_prefix","construct_command_line");
				strcpy(command_line[0].output_prefix,main_argv[i]);
				command_line[0].prefix_flag = 1;
				i++;
			}/*end if*/
                        /*--------------------------------------------------------------*/
			/*		Check if the stream_routing output flag is next.    				*/
			/*--------------------------------------------------------------*/
			else if( strcmp(main_argv[i],"-stro") == 0 ){
				/*--------------------------------------------------------------*/
				/*			Allocate the stream_routing output specifier.				*/
				/*--------------------------------------------------------------*/
				command_line[0].stro = (struct stro_option *)
					alloc(sizeof(struct stro_option), "stro","construct_command_line" );
				command_line[0].stro[0].reachID 	= -999;
				
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if (i < main_argc){
					/*----------------------------------------------*/
					/*Check that the next arguement is a reachID		*/
					/*----------------------------------------------*/
					if ( valid_option(main_argv[i]) == 0){
						command_line[0].stro[0].reachID = (int)atoi(main_argv[i]);
						i++;
					}/*end if*/
				} /*end if*/
			} /*end if*/
			
                        /*--------------------------------------------------------------*/
			/*		Check if the basin output flag is next.    				*/
			/*--------------------------------------------------------------*/
			else if( strcmp(main_argv[i],"-b") == 0 ){
				/*--------------------------------------------------------------*/
				/*			Allocate the basin output specifier.				*/
				/*--------------------------------------------------------------*/
				command_line[0].b = (struct b_option *)
					alloc(sizeof(struct b_option), "b","construct_command_line" );
				command_line[0].b[0].basinID 	= -999;
                               
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if (i < main_argc){
					/*----------------------------------------------*/
					/*Check that the next arguement is a basinID		*/
					/*----------------------------------------------*/
					if ( valid_option(main_argv[i]) == 0){
						command_line[0].b[0].basinID = (int)atoi(main_argv[i]);
						i++;
					}/*end if*/
				} /*end if*/
			} /*end if*/
			/*----------------------------------------------------------*/
			/*		Check if the hillslope output flag is next.  			*/
			/*----------------------------------------------------------*/
			else if( strcmp(main_argv[i],"-h") == 0 ){
				/*-------------------------------------------------------*/
				/*			Allocate the hillslope output specifier.			*/
				/*-------------------------------------------------------*/
				command_line[0].h = (struct h_option *)
					alloc(sizeof(struct h_option), "h", "construct_command_line" );
				command_line[0].h[0].basinID 	= -999;
				command_line[0].h[0].hillID 	= -999;
				/*-------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*-------------------------------------------------------*/
				i++;
				if ( i < main_argc ){
					/*--------------------------------------------------------------*/
					/*				Check that the next arguement is a basinID		*/
					/*--------------------------------------------------------------*/
					if ( valid_option(main_argv[i]) == 0 ){
						command_line[0].h[0].basinID = (int)atoi(main_argv[i]);
						i++;
						/*-------------------------------------------------------*/
						/*		Check that the next arguement exists.		*/
						/*-------------------------------------------------------*/
						if (  i < main_argc ){
							/*-------------------------------------------*/
							/*	Check that the next arguement is hillID	*/
							/*-------------------------------------------*/
							if ( valid_option(main_argv[i]) == 0 ){
								command_line[0].h[0].hillID = (int)atoi(main_argv[i]);
								i++;
							}/*end if*/
						}/*end if*/
					} /*end if*/
				} /*end if*/
			} /*end if*/
			/*-------------------------------------------------------*/
			/*		Check if the zone output flag is next.  				*/
			/*-------------------------------------------------------*/
			else if( strcmp(main_argv[i],"-z") == 0 ){
				/*----------------------------------------------------*/
				/*			Allocate the zone output specifier.				*/
				/*----------------------------------------------------*/
				command_line[0].z = (struct z_option *)
					alloc(sizeof(struct z_option), "z", "construct_command_line");
				command_line[0].z[0].basinID 	= -999;
				command_line[0].z[0].hillID 	= -999;
				command_line[0].z[0].zoneID 	= -999;
				/*-------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*-------------------------------------------------------*/
				i++;
				if (  i < main_argc ){
					/*----------------------------------------------------------*/
					/*				Check that the next arguement is a basinID		*/
					/*----------------------------------------------------------*/
					if ( valid_option(main_argv[i]) == 0 ){
						command_line[0].z[0].basinID = (int)atoi(main_argv[i]);
						i++;
						/*-------------------------------------------------------*/
						/*					Check that the next arguement exists.		*/
						/*-------------------------------------------------------*/
						if (  i < main_argc ){
							/*----------------------------------------------------*/
							/*  			Check that the next arguement is hillID	*/
							/*----------------------------------------------------*/
							if ( valid_option(main_argv[i]) == 0 ){
								command_line[0].z[0].hillID = (int)atoi(main_argv[i]);
								i++;
								/*----------------------------------------------------*/
								/*			Check that the next arguement exists.	*/
								/*----------------------------------------------------*/
								if (  i < main_argc ){
									/*-------------------------------------------------*/
									/*				Check that the next arg is a zoneID	*/
									/*-------------------------------------------------*/
									if ( valid_option(main_argv[i]) == 0 ){
										command_line[0].z[0].zoneID=(int)atoi(main_argv[i]);
										i++;
									}/*end if*/
								}/*end if*/
							} /*end if*/
						} /*end  if*/
					} /*end  if*/
				} /*end if*/
			} /*end if*/
			/*--------------------------------------------------------------*/
			/*		Check if the patch output flag is next.  			*/
			/*--------------------------------------------------------------*/
			else if( strcmp(main_argv[i],"-p") == 0 ){
				/*--------------------------------------------------------------*/
				/*			Allocate the patch output specifier.				*/
				/*--------------------------------------------------------------*/
				command_line[0].p = (struct p_option *)
					alloc(sizeof(struct p_option),"p","construct_command_line" );
				command_line[0].p[0].basinID 	= -999;
				command_line[0].p[0].hillID 	= -999;
				command_line[0].p[0].zoneID 	= -999;
				command_line[0].p[0].patchID 	= -999;
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if (  i < main_argc ){
					/*--------------------------------------------------------------*/
					/*				Check that the next arguement is a basinID		*/
					/*--------------------------------------------------------------*/
					if ( valid_option(main_argv[i]) == 0 ){
						command_line[0].p[0].basinID = (int)atoi(main_argv[i]);
						i++;
						/*-------------------------------------------------------*/
						/*					Check that the next arguement exists.		*/
						/*-------------------------------------------------------*/
						if (  i < main_argc ){
							/*-------------------------------------------------------*/
							/*		  			Check that the next arguement is hillID	*/
							/*-------------------------------------------------------*/
							if ( valid_option(main_argv[i]) == 0){
								command_line[0].p[0].hillID = (int)atoi(main_argv[i]);
								i++;
								/*----------------------------------------------------*/
								/*	  Check that the next arguement exists.				*/
								/*----------------------------------------------------*/
								if (  i < main_argc ){
									/*----------------------------------------------*/
									/* 				Check that next arg is a zoneID		*/
									/*-----------------------------------------------*/
									if ( valid_option(main_argv[i]) == 0  ){
										command_line[0].p[0].zoneID = (int)atoi(main_argv[i]);
										i++;
										/*--------------------------------------------*/
										/*			Check that next arguement exists.	*/
										/*-------------------------------------------*/
										if (  i < main_argc ){
											/*------------------------------------------*/
											/*			Check next arg is a patchID		*/
											/*----------------------------------------*/
											if ( valid_option(main_argv[i]) == 0 ){
												command_line[0].p[0].patchID =
													(int)atoi(main_argv[i]);
												i++;
											}/*end if*/
										} /*end if*/
									} /*end if*/
								} /*end if*/
							} /*end if*/
						} /*end if*/
					} /*end if*/
				} /*end if*/
			} /*end if*/
			/*--------------------------------------------------------------*/
			/*		Check if the canopy stratum output flag is next.		*/
			/*--------------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-c") == 0 ){
				/*--------------------------------------------------------------*/
				/*			Allocate the patch output specifier.				*/
				/*--------------------------------------------------------------*/
				command_line[0].c = (struct c_option *)
					alloc(sizeof(struct c_option),"c","construct_command_line" );
				command_line[0].c[0].basinID 	= -999;
				command_line[0].c[0].hillID 	= -999;
				command_line[0].c[0].zoneID 	= -999;
				command_line[0].c[0].patchID 	= -999;
				command_line[0].c[0].stratumID 	= -999;
				/*--------------------------------------------------------------*/
				/*			Check that the next arguement exists.				*/
				/*--------------------------------------------------------------*/
				i++;
				if (  i < main_argc ){
					/*--------------------------------------------------------------*/
					/*				Check that the next arguement is a basinID		*/
					/*--------------------------------------------------------------*/
					if ( valid_option(main_argv[i]) == 0 ){
						command_line[0].c[0].basinID = (int)atoi(main_argv[i]);
						i++;
						/*--------------------------------------------------------------*/
						/*					Check that the next arguement exists.		*/
						/*--------------------------------------------------------------*/
						if (  i < main_argc ){
							/*----------------------------------------------------------*/
							/*						Check that the next arguement is hillID	*/
							/*----------------------------------------------------------*/
							if ( valid_option(main_argv[i]) == 0){
								command_line[0].c[0].hillID = (int)atoi(main_argv[i]);
								i++;
								/*-------------------------------------------------------*/
								/*			Check that the next arguement exists.				*/
								/*-------------------------------------------------------*/
								if (  i < main_argc ){
									/*----------------------------------------------------*/
									/*					Check that next arg is a zoneID		*/
									/*----------------------------------------------------*/
									if ( valid_option(main_argv[i]) == 0  ){
										command_line[0].c[0].zoneID =	(int)atoi(main_argv[i]);
										i++;
										/*-------------------------------------------------*/
										/*				Check that next arguement exists.	*/
										/*-------------------------------------------------*/
										if (  i < main_argc ){
											/*---------------------------------------------*/
											/*				Check next arg is a patchID		*/
											/*---------------------------------------------*/
											if ( valid_option(main_argv[i]) == 0 ){
												command_line[0].c[0].patchID =
													(int)atoi(main_argv[i]);
												i++;
												/*------------------------------------------*/
												/*		Check that next arguement exists.	*/
												/*-----------------------------------------*/
												if (  i < main_argc ){
													/*----------------------------------------*/
													/*	Check next arg is a stratumID	*/
													/*----------------------------------------*/
													if ( valid_option(main_argv[i]) == 0 ){
														command_line[0].c[0].stratumID =
															(int)atoi(main_argv[i]);
														i++;
													}/*end if*/
												} /*end if*/
											}/*end if*/
										} /*end if*/
									} /*end if*/
								} /*end if*/
							} /*end if*/
						   } /*end if*/
						} /*end if*/
					} /*end if*/
				} /*end if*/
			/*--------------------------------------------------------------*/
			/*		Check if the gridded climate input flag is next		    */
			/*--------------------------------------------------------------*/
			else if ( strcmp(main_argv[i],"-asciigrid") == 0 ){
				printf("Setting climate mode to gridded ascii\n");
				command_line[0].gridded_ascii_flag = 1;
				i++;
				} /*end if*/
			else if (strcmp(main_argv[i],"-netcdfgrid") == 0 ){
				command_line[0].gridded_netcdf_flag = 1;
				i++;
			}
			else if (strcmp(main_argv[i], "-longwaveevap") == 0) {
				command_line[0].evap_use_longwave_flag = 1;
				i++;
			}
			/*--------------------------------------------------------------*/
			/*	NOTE:  ADD MORE OPTION PARSING HERE.						*/
			/*--------------------------------------------------------------*/
			/*--------------------------------------------------------------*/
			/*		the option must be invalid.								*/
			/*--------------------------------------------------------------*/
			else{
				fprintf(stderr,
					"FATAL ERROR: in construct_command_line option #%d is invalid.\n",i);
				fprintf(stderr,"for argument %s\n", main_argv[i]);
				exit(EXIT_FAILURE);
			} /*end if*/
		} /*end if*/
	} /*end while*/

	return(command_line);
} /*end construct_command_line*/
