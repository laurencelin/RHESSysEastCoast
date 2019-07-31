/*--------------------------------------------------------------*/
/* 																*/
/*					output_csv_patch						*/
/*																*/
/*	output_csv_patch - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_csv_patch - outputs current contents of a patch.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_csv_patch(										*/
/*					struct	patch_object	*patch,				*/
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

void	output_csv_patch(
					 int basinID, int hillID, int zoneID,
					 struct	patch_object	*patch,
					 struct	date	current_date,
					 FILE *outfile)
{
	printf("this function has been removed./n");
	return;
} /*end output_csv_patch*/
