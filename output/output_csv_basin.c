/*--------------------------------------------------------------*/
/* 																*/
/*					output_csv_basin						*/
/*																*/
/*	output_csv_basin - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_csv_basin - outputs current contents of a basin.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_csv_basin( int routing_flag,										*/	
/*					struct	basin_object	*basin,				*/
/*					struct	date	date,  						*/
/*					FILE 	*outfile)							*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	outputs spatial structure according to commandline			*/
/*	specifications to specific files							*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	We only permit one fileset per spatial modelling level.     */
/*	Each fileset has one file for each timestep.  				*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

void	output_csv_basin(			int routing_flag,
					 struct	basin_object	*basin,
					 struct	date	date,
					 FILE *outfile)
{
	printf("this function has been removed./n");

	return;
} /*end output_csv_basin*/
