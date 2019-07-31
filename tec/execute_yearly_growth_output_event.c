/*--------------------------------------------------------------*/
/* 																*/
/*					execute_yearly_growth_output_event					*/
/*																*/
/*	execute_yearly_growth_output_event - output_growths yearly data			*/
/*																*/
/*	NAME														*/
/*	execute_yearly_growth_output_event - output_growths yearly data 	.		*/
/*																*/
/*	SYNOPSIS													*/
/*	void	execute_yearly_growth_output_event(						*/
/*					struct	world_object	*world,				*/
/*					struct	command_line_object *command_line,	*/
/*					struct	date	date,  						*/
/*					struct	world_output_file_object 	*outfile)*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	output_growths spatial structure according to commandline			*/
/*	specifications to specific files, for yearly info			*/
/*	a spatial structure is output_growth only if its appropriate		*/
/*	option flag has been set and its ID matches a 				*/
/*	specified ID, or -999 which indicates all					*/
/*	units at that spatial scale are to be output_growth				*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	We only permit one fileset per spatial modelling level.     */
/*	Each fileset has one file for each timestep.  				*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

void	execute_yearly_growth_output_event(
										   struct	world_object	*world,
										   struct	command_line_object *command_line,
										   struct	date	date,
										   struct	world_output_file_object	*outfile)
{
    printf("this function has been removed.\n");
		return;
} /*end execute_yearly_growth_output_event*/
