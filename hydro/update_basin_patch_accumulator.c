/*--------------------------------------------------------------------------------------*/
/* 											*/
/*			update_basin_patch_accumulator					*/
/*											*/
/*	NAME										*/
/*	update_basin_patch_accumulator.c - update accumulator variables at the end of day	*/
/*					this process is taken from compute_subsurface_routing.c	*/
/*	SYNOPSIS									*/
/*	void update_basin_patch_accumulator( 						*/
/*					struct command_line_object *command_line,	*/
/*					struct basin_object *basin			*/
/*					struct date current_date)			*/
/*											*/
/* 											*/
/*											*/
/*	OPTIONS										*/
/*											*/
/*											*/
/*	DESCRIPTION									*/
/*	this function is called in basin_daily_F at the end of each day, it was in  	*/
/*	the compute_subsurface_routing, 										*/
/*											*/
/*											*/
/*	PROGRAMMER NOTES								*/
/*											*/
/*											*/	
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"
#define NUM_NORMAL  10

void update_basin_patch_accumulator(
			struct command_line_object 	*command_line,
			struct basin_object 		*basin,
			struct date		 	current_date)
{
    /*----------------------------------------------------------------------*/
    /* Local variables definition                                           */
    /*-----------------------------------------------------------------------*/
    double scale;
    double tmp;
    struct patch_object *patch;
    int b,h,p,z,c,s;
    double alai = 0.0;
    double coverf = 0.0;
    int layer, i;
    double aa, bb, cc, dd;
    double water_scalar0, water_scalar1, water_scalari, thetai, theta;
    double water_scalar2, water_scalar3;
    double NORMAL[10]= { 0.0, -1.283,-0.842,-0.524,-0.253, 0.0, 0.253,0.524,0.842,1.283};
    /*----------------------------------------------------------------------*/
    /* initializations                                                   */
    /*----------------------------------------------------------------------*/

    /*---------------------------------------------------------------------*/
    /*update accumulator variables                                            */
    /*-----------------------------------------------------------------------*/
    
    
    
    for (h=0; h < basin->num_hillslopes; ++h) {
        for(z=0; z < basin->hillslopes[h][0].num_zones; ++z) {
            for (p=0; p < basin->hillslopes[h][0].zones[z][0].num_patches; p++) {
                
                alai = 0.0;
                patch=basin->hillslopes[h]->zones[z]->patches[p];
        
                for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
                    for ( c=0 ; c<patch[0].layers[layer].count; c++ ){
                        coverf = patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cover_fraction;
                        alai += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].epv.proj_lai;
                    }// for c
                }//for layer
                
                // ------------ denitrification
                // from update_denitrif
                if (patch[0].soil_defaults[0][0].soil_type.sand > 0.5) {
                    aa = 1.56; bb=12.0; cc=16.0; dd=2.01;
                } else if (patch[0].soil_defaults[0][0].soil_type.clay > 0.5) {
                    aa = 60.0; bb=18.0; cc=22.0; dd=1.06;
                } else {
                    aa=4.82; bb=14.0; cc=16.0; dd=1.39;
                }
                
                if(patch[0].rootzone.potential_sat>ZERO){
                    if (patch[0].sat_deficit > patch[0].rootzone.potential_sat) theta = min(patch[0].rz_storage/patch[0].rootzone.potential_sat, 1.0);//(1.0-patch[0].basementFrac)
                    else theta = min((patch[0].rz_storage + patch[0].rootzone.potential_sat - patch[0].sat_deficit)/patch[0].rootzone.potential_sat,1.0);//(1.0-patch[0].basementFrac)
                }else{ theta = 0.0; }
                
                water_scalar0 = min( aa*exp(-cc*exp(-dd*theta*log(bb))*log(bb)), 1.0);
                water_scalar1 = 0.0;
                for (i = 1; i< NUM_NORMAL; i++) {
                    thetai = theta + patch[0].theta_std * NORMAL[i];
                    thetai = min(1.0, thetai);
                    thetai = max(0.002, thetai);
                    water_scalari = min( aa*exp(-cc*exp(-dd*thetai*log(bb))*log(bb)), 1.0);
                    water_scalar1 += water_scalari;
                }// for i
                water_scalar1 *= 1.0/NUM_NORMAL;
                
                // ------- nitrification
                if (patch[0].soil_defaults[0][0].soil_type.sand > 0.5) {
                   aa = 0.55; bb=1.7; cc=-0.007; dd=3.22;
                } else {
                   aa=0.6; bb=1.27; cc=0.0012; dd=2.84;
                }// if
                
                water_scalar2 = (theta > cc? exp(dd*(bb-aa)/(aa-cc)*log((thetai -bb) / (aa-bb))) * exp(dd*log((thetai-cc)/ (aa-cc))) : 0.000001);
                water_scalar2 = min(water_scalar2, 1.0);
                water_scalar3 = 0.0;
                for (i=0; i<NUM_NORMAL; i++) {
                    thetai = theta + NORMAL[i]*patch[0].theta_std;
                    thetai = min(1.0, thetai);
                    thetai = max(0.002, thetai);
                    water_scalar3  +=  exp(dd*(bb-aa)/(aa-cc)*log((thetai -bb) / (aa-bb))) * exp(dd*log((thetai-cc)/ (aa-cc)));
                }//for
                water_scalar3 *= 1.0/NUM_NORMAL;
                water_scalar3 = min(water_scalar3, 1.0);
                
                

                
                // annual monthly
                patch[0].acc_month.subQnet += (patch[0].Qout_total - patch[0].Qin_total);
                patch[0].acc_month.surfQnet += (patch[0].surface_Qout_total - patch[0].surface_Qin_total);
                patch[0].acc_month.subQvnet += (patch[0].unsat_drainage - patch[0].cap_rise);
                patch[0].acc_month.precip += basin->hillslopes[h][0].zones[z][0].rain_hourly_total+basin->hillslopes[h][0].zones[z][0].rain+basin->hillslopes[h][0].zones[z][0].snow;
                //patch[0].acc_month.recharge += patch[0].recharge; // need to track net recharge
                patch[0].acc_month.PET += patch[0].PET;
                patch[0].acc_month.ET += (patch[0].transpiration_sat_zone + patch[0].transpiration_unsat_zone + patch[0].evaporation + patch[0].evaporation_surf  + patch[0].exfiltration_sat_zone + patch[0].exfiltration_unsat_zone);
                patch[0].acc_month.sat_deficit_z += patch[0].sat_deficit_z;
                patch[0].acc_month.peakLAI = max(patch[0].acc_month.peakLAI,alai);
                patch[0].acc_month.meanLAI += alai;
                patch[0].acc_month.psn += patch[0].net_plant_psn;
                patch[0].acc_month.days += 1.0;
                patch[0].acc_month.satChance += (patch[0].sat_deficit<=0? 1.0:0.0);
                patch[0].acc_month.plantlimitN += patch[0].soil_ns.fract_potential_uptake;
                patch[0].acc_month.plantlimitQ += patch[0].trans_reduc_perc;
                
                patch[0].acc_month.activeS += theta;
                patch[0].acc_month.denitrifaspQs += water_scalar0;
                patch[0].acc_month.denitrifspQs += water_scalar1;
                patch[0].acc_month.denitrif += patch[0].ndf.denitrif;
                patch[0].acc_month.nitrifaspQs += water_scalar2;
                patch[0].acc_month.nitrifspQs += water_scalar3;
                patch[0].acc_month.mineralization += patch[0].ndf.net_mineralized; // immob = (patch[0].ndf.net_mineralized - patch[0].ndf.mineralized)// positive is mineralization
                patch[0].acc_month.uptake += patch[0].ndf.sminn_to_npool;
                patch[0].acc_month.subNO3net += patch[0].soil_ns.NO3_Qout_total - patch[0].soil_ns.NO3_Qin_total;
                patch[0].acc_month.subNO3vnet += patch[0].sat_NO3 - patch[0].acc_month.subNO3vnet;
                patch[0].acc_month.subDOCnet += patch[0].soil_cs.DOC_Qout_total - patch[0].soil_cs.DOC_Qin_total;
                patch[0].acc_month.no3drain2gw += patch[0].gw_drainage_NO3;
                
                
                // annual
                patch[0].acc_year.subQnet += (patch[0].Qout_total - patch[0].Qin_total);
                patch[0].acc_year.surfQnet += (patch[0].surface_Qout_total - patch[0].surface_Qin_total);
                patch[0].acc_year.subQvnet += (patch[0].unsat_drainage - patch[0].cap_rise);
                patch[0].acc_year.precip += basin->hillslopes[h][0].zones[z][0].rain_hourly_total+basin->hillslopes[h][0].zones[z][0].rain+basin->hillslopes[h][0].zones[z][0].snow;
                //patch[0].acc_year.recharge += patch[0].recharge;
                patch[0].acc_year.PET += patch[0].PET;
                patch[0].acc_year.ET += (patch[0].transpiration_sat_zone + patch[0].transpiration_unsat_zone + patch[0].evaporation + patch[0].evaporation_surf  + patch[0].exfiltration_sat_zone + patch[0].exfiltration_unsat_zone);
                patch[0].acc_year.sat_deficit_z += patch[0].sat_deficit_z;
                patch[0].acc_year.peakLAI = max(patch[0].acc_year.peakLAI,alai);
                patch[0].acc_year.meanLAI += alai;
                patch[0].acc_year.psn += patch[0].net_plant_psn;
                patch[0].acc_year.days += 1.0;
                patch[0].acc_year.satChance += (patch[0].sat_deficit<=0? 1.0:0.0);
                patch[0].acc_year.plantlimitN += patch[0].soil_ns.fract_potential_uptake;
                patch[0].acc_year.plantlimitQ += patch[0].trans_reduc_perc;
                
                patch[0].acc_year.activeS += theta;
                patch[0].acc_year.denitrifaspQs += water_scalar0;
                patch[0].acc_year.denitrifspQs += water_scalar1;
                patch[0].acc_year.denitrif += patch[0].ndf.denitrif;
                patch[0].acc_year.nitrifaspQs += water_scalar2;
                patch[0].acc_year.nitrifspQs += water_scalar3;
                patch[0].acc_year.mineralization += patch[0].ndf.net_mineralized;
                patch[0].acc_year.uptake += patch[0].ndf.sminn_to_npool;
                patch[0].acc_year.subNO3net += patch[0].soil_ns.NO3_Qout_total - patch[0].soil_ns.NO3_Qin_total;
                patch[0].acc_year.subNO3vnet += patch[0].sat_NO3 - patch[0].acc_year.subNO3vnet;
                patch[0].acc_year.subDOCnet += patch[0].soil_cs.DOC_Qout_total - patch[0].soil_cs.DOC_Qin_total;
                patch[0].acc_year.no3drain2gw += patch[0].gw_drainage_NO3;
                
        
            } /* end of p*/
        } /* end of z*/
    } /* end of h*/
    
    
	return;
} /* end of update_basin_patch_accumulator.c */
