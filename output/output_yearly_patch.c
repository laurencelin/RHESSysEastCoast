/*--------------------------------------------------------------*/
/* 																*/
/*					output_yearly_patch						*/
/*																*/
/*	output_yearly_patch - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_yearly_patch - outputs current contents of a patch.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_yearly_patch( int basinID, int hillID, int zoneID,								*/
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

void	output_yearly_patch(
				int basinID, int hillID, int zoneID,
				struct	patch_object	*patch,
				struct	date	current_date,
				FILE *outfile)
{
    int check;
	check = fprintf(outfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
    
    current_date.year,
    patch[0].ID,
    patch[0].acc_year.subQnet*1000.0,
    patch[0].acc_year.surfQnet*1000.0,
    patch[0].acc_year.subQvnet*1000.0,//5
    patch[0].acc_year.precip*1000.0,
    patch[0].acc_year.PET*1000.0,
    patch[0].acc_year.ET*1000.0,
    patch[0].acc_year.sat_deficit_z*1000.0 / patch[0].acc_year.days,
    patch[0].acc_year.peakLAI,//10
    patch[0].acc_year.meanLAI/ patch[0].acc_year.days,
    patch[0].acc_year.psn*1000.0,
    patch[0].acc_year.denitrif*1000.0,
    patch[0].acc_year.mineralization*1000.0,
    patch[0].acc_year.uptake*1000.0,//15
    patch[0].acc_year.subNO3net*1000.0,
    patch[0].acc_year.subNO3vnet*1000.0,
    patch[0].acc_year.subDOCnet*1000.0,
    patch[0].acc_year.no3drain2gw*1000.0,
    patch[0].acc_year.satChance/ patch[0].acc_year.days,//20
    patch[0].acc_year.plantlimitN/ patch[0].acc_year.days,
    patch[0].acc_year.plantlimitQ/ patch[0].acc_year.days
    );
    
	return;
    
} /*end output_csv_yearly_patch*/

