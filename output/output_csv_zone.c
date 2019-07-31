/*--------------------------------------------------------------*/
/* 																*/
/*					output_csv_zone						*/
/*																*/
/*	output_csv_zone - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_csv_zone - outputs current contents of a zone.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_csv_zone(										*/
/*					struct	zone_object	*zone,				*/
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

void	output_csv_zone(	int basinID, int hillID,
						struct	zone_object	*zone,
						struct	date	current_date,
						FILE *outfile)
{
	printf("this function has been removed./n");
	return;
} /*end output_csv_zone*/
