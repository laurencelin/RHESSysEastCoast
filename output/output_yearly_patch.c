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
	check = fprintf(outfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
    
    current_date.year,
    patch[0].ID,//3
    patch[0].acc_year.subQnet*1000.0,
    patch[0].acc_year.surfQnet*1000.0,
    patch[0].acc_year.precip*1000.0,
    patch[0].acc_year.recharge*1000.0,
    patch[0].acc_year.PET*1000.0,
    patch[0].acc_year.ET*1000.0,
    patch[0].acc_year.sat_deficit_z*1000.0 / patch[0].acc_year.days,
    patch[0].acc_year.peakLAI,
    patch[0].acc_year.meanLAI/ patch[0].acc_year.days,
    patch[0].acc_year.psn*1000.0, //12
    patch[0].acc_year.denitrif*1000.0,
    patch[0].acc_year.mineralization*1000.0,
    patch[0].acc_year.uptake*1000.0,
    patch[0].acc_year.subNO3net*1000.0,
    patch[0].acc_year.subDOCnet*1000.0
    );
    
    patch[0].acc_year.subQnet = 0.0;
    patch[0].acc_year.surfQnet = 0.0;
    patch[0].acc_year.precip = 0.0;
    patch[0].acc_year.recharge = 0.0;
    patch[0].acc_year.PET = 0.0;
    patch[0].acc_year.ET = 0.0;
    patch[0].acc_year.sat_deficit_z = 0.0;
    patch[0].acc_year.peakLAI = 0.0;
    patch[0].acc_year.meanLAI = 0.0;
    patch[0].acc_year.psn = 0.0;
    patch[0].acc_year.days = 0.0;
    patch[0].acc_year.denitrif = 0.0;
    patch[0].acc_year.mineralization = 0.0;
    patch[0].acc_year.uptake = 0.0;
    patch[0].acc_year.subNO3net = 0.0;
    patch[0].acc_year.subDOCnet = 0.0;
    
	return;
    
} /*end output_csv_yearly_patch*/

