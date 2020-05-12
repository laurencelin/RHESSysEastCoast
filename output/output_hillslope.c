/*--------------------------------------------------------------*/
/* 																*/
/*					output_hillslope						*/
/*																*/
/*	output_hillslope - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_hillslope - outputs current contents of a hillslope.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_hillslope(										*/
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

void	output_hillslope(int basinID,
						 struct	basin_object	*basin,
                         int hID,
						 struct	date	date,
						 FILE *outfile)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	int z,p,c,hh;
	int layer;
	double asat_deficit_z;
	double asat_deficit;
	double aunsat_storage;
	double aunsat_drainage;
	double acap_rise;
	double areturn_flow;
	double aevaporation;
	double atranspiration;
	double abase_flow;
    double astormdrain;
    double apsn;
    double alaitree;
    double alaigrass;
	double asat_area;
    double adetention_store;
    double alawnirrigation;
	double aarea;
    double arz_storage;
    double asnowmelt;
    double arecharge;
    double zone_area;
    double apcp;
	struct	patch_object  *patch;
	struct	zone_object	*zone;
    double agwstorage;
    double agwOut;
	/*--------------------------------------------------------------*/
	/*	Initialize Accumlating variables.								*/
	/*--------------------------------------------------------------*/
	asat_deficit_z = 0.0 ;
	asat_deficit = 0.0 ;
	aunsat_storage = 0.0 ;
	aunsat_drainage = 0.0 ;
	acap_rise = 0.0 ;
	areturn_flow = 0.0 ;
	aevaporation = 0.0 ;
	atranspiration = 0.0  ;
	abase_flow = 0.0;
    astormdrain = 0.0;
	apsn = 0.0 ;
	alaitree = 0.0;
    alaigrass = 0.0;
    alawnirrigation = 0.0;
	aarea =  0.0 ;
    arz_storage = 0.0;
    asnowmelt = 0.0;
    adetention_store = 0.0;
    asat_area = 0.0;
    arecharge = 0.0;
    zone_area = 0.0;
    apcp = 0.0;
    agwstorage = 0.0;
    agwOut = 0.0;
    
    struct hillslope_object *hillslope;
    for (hh=0; hh < basin[0].num_hillslopes; hh++){
        hillslope = basin[0].hillslopes[hh];
        if(hID%2==0 && (hillslope[0].ID == hID || hillslope[0].ID == hID-1)){
         
            for (z=0; z<hillslope[0].num_zones; z++){
                    zone = hillslope[0].zones[z];
                    apcp += (zone[0].rain_hourly_total+zone[0].rain+zone[0].snow)*zone[0].area;
                    zone_area += zone[0].area;
                    
                    for (p=0; p< zone[0].num_patches; p++){
                        patch = zone[0].patches[p];
                        asat_deficit_z += patch[0].sat_deficit_z * patch[0].area;
                        asat_deficit += patch[0].sat_deficit * patch[0].area;
                        
                        arz_storage += patch[0].rz_storage * patch[0].area;
                        asnowmelt += patch[0].snow_melt*patch[0].area;
                        adetention_store += patch[0].detention_store*patch[0].area;
                        
                        aunsat_drainage += patch[0].unsat_drainage * patch[0].area;
                        acap_rise += patch[0].cap_rise * patch[0].area;
                        abase_flow += (patch[0].base_flow + hillslope[0].base_flow)* patch[0].area;
                        areturn_flow += patch[0].return_flow * patch[0].area;
                        astormdrain += patch[0].stormdrained * patch[0].area;
                        aevaporation += patch[0].evaporation * patch[0].area;
                        aarea += patch[0].area;
                        arecharge += patch[0].recharge * patch[0].area;
                        if (patch[0].sat_deficit <= ZERO) asat_area += patch[0].area;
                        alawnirrigation += patch[0].grassIrrigation_m * patch[0].area;
                        atranspiration += (patch[0].transpiration_sat_zone
                            + patch[0].transpiration_unsat_zone)  *  patch[0].area;
                        
                        agwstorage += hillslope[0].gw.storage * patch[0].area;
                        agwOut += hillslope[0].gw.Qout * patch[0].area;
                        
                        for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
                            for ( c=0 ; c<patch[0].layers[layer].count; c++ ){
                                
                                apsn += patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cover_fraction
                                    * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cs.net_psn
                                    * patch[0].area;
                                
                                if(patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].epc.veg_type == TREE && patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].ID!=802){
                                    alaitree += patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cover_fraction * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].epv.proj_lai * patch[0].area;
                                }//if tree
                                if(patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].epc.veg_type == GRASS || patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].ID==802){
                                    alaigrass += patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cover_fraction * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].epv.proj_lai * patch[0].area;
                                }//if grass
                                
                            }// for c
                        }// for layer
                    }// for patch
                }// for zone
        }// end of if
    }//end of for loop
    
    if(hID%2==0){
        // after aggregate all hillslopes
        aarea = 1.0/aarea;
        asat_deficit_z *= aarea ;
        asat_deficit *= aarea ;
        aunsat_drainage *= aarea ;
        acap_rise *= aarea ;
        areturn_flow *= aarea ;
        astormdrain *= aarea ;
        aevaporation *= aarea ;
        abase_flow *= aarea;
        atranspiration *= aarea  ;
        apsn *= aarea ;
        alaitree *= aarea ;
        alaigrass *= aarea ;
        arz_storage *= aarea ;
        asnowmelt *= aarea ;
        adetention_store *= aarea ;
        asat_area *= aarea ;
        alawnirrigation *= aarea ;
        arecharge *= aarea;
        agwstorage *= aarea;
        agwOut *= aarea;
        apcp /= zone_area;
        //abase_flow += hillslope[0].base_flow;
        
        fprintf(outfile,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
            date.day, //1
            date.month, //2
            date.year, //3
            hID, //4
            (1.0/aarea), //5
            asat_deficit_z * 1000.0, // 6 mm
            asat_deficit * 1000.0, //7 mm
            adetention_store * 1000.0, //8 mm
            asat_area * 100.0, // 9 %
            arz_storage * 1000.0,// 10
            acap_rise * 1000.0, //11 mm
            aunsat_drainage * 1000.0, //12 mm
            abase_flow * 1000.0, // 13 mm
            areturn_flow * 1000.0, // 14 mm
            (areturn_flow + abase_flow + astormdrain)* 1000.0, // 15 mm
            agwOut *1000.0, // 16 mm
            agwstorage *1000.0, //17 mm <---- this is not correct
            asnowmelt * 1000.0, // 18 mm
            apsn, //19 kgC/m2/d
            aevaporation * 1000.0, // 20 mm
            atranspiration * 1000.0, //21 mm
            alaitree, // 22
            alaigrass, //23
            alawnirrigation*1000.0, // 24 mm
            apcp*1000.0, //25 mm
            arecharge*1000.0 // 26 mm
            );
    }
	return;
} /*end output_hillslope*/
