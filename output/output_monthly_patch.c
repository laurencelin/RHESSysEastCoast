/*--------------------------------------------------------------*/
/* 																*/
/*					output_monthly_patch						*/
/*																*/
/*	output_monthly_patch - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_monthly_patch - outputs current contents of a patch.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_monthly_patch(										*/
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

void	output_monthly_patch(
							 int basinID, int hillID, int zoneID,
							 struct	patch_object	*patch,
							 struct	date	current_date,
							 FILE *outfile)
{
    int check;
    check = fprintf(outfile,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
    
    current_date.year,
    current_date.month,
    patch[0].ID,//3
    patch[0].acc_month.subQnet*1000.0,
    patch[0].acc_month.surfQnet*1000.0,
    patch[0].acc_month.subQvnet*1000.0,
    patch[0].acc_month.precip*1000.0,
    patch[0].acc_month.recharge*1000.0,
    patch[0].acc_month.PET*1000.0,
    patch[0].acc_month.ET*1000.0,
    patch[0].acc_month.sat_deficit_z*1000.0 / patch[0].acc_month.days,
    patch[0].acc_month.peakLAI,
    patch[0].acc_month.meanLAI/ patch[0].acc_month.days,
    patch[0].acc_month.psn*1000.0, //12
    patch[0].acc_month.denitrif*1000.0,
    patch[0].acc_month.mineralization*1000.0,
    patch[0].acc_month.uptake*1000.0,
    patch[0].acc_month.subNO3net*1000.0,
    patch[0].acc_month.subNO3vnet*1000.0,
    patch[0].acc_month.subDOCnet*1000.0
    );
    
    patch[0].acc_month.subQnet = 0.0;
    patch[0].acc_month.surfQnet = 0.0;
    patch[0].acc_month.subQvnet = 0.0;
    patch[0].acc_month.precip = 0.0;
    patch[0].acc_month.recharge = 0.0;
    patch[0].acc_month.PET = 0.0;
    patch[0].acc_month.ET = 0.0;
    patch[0].acc_month.sat_deficit_z = 0.0;
    patch[0].acc_month.peakLAI = 0.0;
    patch[0].acc_month.meanLAI = 0.0;
    patch[0].acc_month.psn = 0.0;
    patch[0].acc_month.days = 0.0;
    patch[0].acc_month.denitrif = 0.0;
    patch[0].acc_month.mineralization = 0.0;
    patch[0].acc_month.uptake = 0.0;
    patch[0].acc_month.subNO3net = 0.0;
    patch[0].acc_month.subNO3vnet = 0.0;
    patch[0].acc_month.subDOCnet = 0.0;
    
	return;
} /*end output_monthly_patch*/
