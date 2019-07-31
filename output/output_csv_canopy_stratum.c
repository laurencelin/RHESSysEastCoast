/*--------------------------------------------------------------*/
/* 																*/
/*					output_csv_canopy_stratum						*/
/*																*/
/*	output_csv_canopy_stratum - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_csv_canopy_stratum - outputs current contents of a canopy_stratum.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_csv_canopy_stratum(										*/
/*					struct	canopy_stratum_object	*canopy_stratum,				*/
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

void	output_csv_canopy_stratum( int basinID, int hillID, int zoneID, int patchID,
							  struct	canopy_strata_object	*stratum,
							  struct	date	current_date,
							  FILE *outfile)
{
	printf("this function has been removed./n");
	return;
} /*end output_csv_canopy_stratum*/
