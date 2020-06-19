/*--------------------------------------------------------------*/
/* 																*/
/*					output_yearly_basin						*/
/*																*/
/*	output_yearly_basin - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_yearly_basin - outputs current contents of a basin.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_yearly_basin(										*/
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

void	output_yearly_basin(
							 struct	basin_object	*basin,
							 struct	date	current_date,
							 FILE *outfile)
{
    
    int z,p,hh;
    struct hillslope_object *hillslope;
    struct    patch_object  *patch;
    struct    zone_object    *zone;
    double aarea = 0.0;
    double subQnet = 0.0;
    double surfQnet = 0.0;
    double subQvnet = 0.0;
    double precip = 0.0;
    double PET = 0.0;
    double ET = 0.0;
    double sat_deficit_z = 0.0;
    double peakLAI = 0.0;
    double meanLAI = 0.0;
    double psn = 0.0;
    double denitrif = 0.0;
    double mineralization = 0.0;
    double uptake = 0.0;
    double subNO3net = 0.0;
    double subNO3vnet = 0.0;
    double subDOCnet = 0.0;
    double no3drain2gw = 0.0;
   
    
    for (hh=0; hh < basin[0].num_hillslopes; hh++){
        hillslope = basin[0].hillslopes[hh];
        for (z=0; z<hillslope[0].num_zones; z++){
            zone = hillslope[0].zones[z];
            for (p=0; p< zone[0].num_patches; p++){
                patch = zone[0].patches[p];
                aarea += patch[0].area;
                subQnet += patch[0].acc_year.subQnet*patch[0].area;
                surfQnet += patch[0].acc_year.surfQnet*patch[0].area;
                subQvnet += patch[0].acc_year.subQvnet*patch[0].area;
                precip += patch[0].acc_year.precip*patch[0].area;
                PET += patch[0].acc_year.PET*patch[0].area;
                ET += patch[0].acc_year.ET*patch[0].area;
                sat_deficit_z += patch[0].acc_year.sat_deficit_z / patch[0].acc_year.days *patch[0].area;
                peakLAI += patch[0].acc_year.peakLAI*patch[0].area;
                meanLAI += patch[0].acc_year.meanLAI/ patch[0].acc_year.days*patch[0].area;
                psn += patch[0].acc_year.psn*patch[0].area;
                denitrif += patch[0].acc_year.denitrif*patch[0].area;
                mineralization += patch[0].acc_year.mineralization*patch[0].area;
                uptake += patch[0].acc_year.uptake*patch[0].area;
                subNO3net += patch[0].acc_year.subNO3net*patch[0].area;
                subNO3vnet += patch[0].acc_year.subNO3vnet*patch[0].area;
                subDOCnet += patch[0].acc_year.subDOCnet*patch[0].area;
                no3drain2gw += patch[0].acc_year.no3drain2gw*patch[0].area;
            }//p
        }//z
    }//hh
    
    aarea = 1.0/aarea;
    subQnet *= aarea;
    surfQnet *= aarea;
    subQvnet *= aarea;
    precip *= aarea;
    PET *= aarea;
    ET *= aarea;
    sat_deficit_z *= aarea;
    peakLAI *= aarea;
    meanLAI *= aarea;
    psn *= aarea;
    denitrif *= aarea;
    mineralization *= aarea;
    uptake *= aarea;
    subNO3net *= aarea;
    subNO3vnet *= aarea;
    subDOCnet *= aarea;
    no3drain2gw *= aarea;
   
    
	fprintf(outfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
        current_date.year,
        subQnet*1000.0,
        surfQnet*1000.0,
        subQvnet*1000.0,
        precip*1000.0,
        PET*1000.0,
        ET*1000.0,
        sat_deficit_z*1000.0,
        peakLAI,
        meanLAI,
        psn*1000.0,
        denitrif*1000.0,
        mineralization*1000.0,
        uptake*1000.0,
        subNO3net*1000.0,
        subNO3vnet*1000.0,
        subDOCnet*1000.0,
        no3drain2gw*1000.0
        );
	return;
} /*end output_yearly_basin*/
