/*--------------------------------------------------------------*/
/* 																*/
/*					output_monthly_hillslope						*/
/*																*/
/*	output_monthly_hillslope - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_monthly_hillslope - outputs current contents of a hillslope.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_monthly_hillslope(										*/
/*					struct	hillslope_object	*hillslope,				*/
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

void	output_monthly_hillslope(	int basinID,
							 struct	hillslope_object	*hillslope,
							 struct	date	current_date,
							 FILE *outfile)
{
	printf("this function has been removed./n");
	return;
} /*end output_monthly_hillslope*/
