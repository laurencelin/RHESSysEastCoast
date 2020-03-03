/*--------------------------------------------------------------*/
/*                                                              */
/*		update_denitrif				*/
/*                                                              */
/*  NAME                                                        */
/*		update_denitrif				*/
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/*  void update_denitrif(				*/
/*                                                              */
/*			struct  soil_c_object   *               */
/*                      struct  soil_n_object   *               */
/*                      struct  cdayflux_patch_object *         */
/*                      struct  ndayflux_patch_object *         */
/*			struct	soil_class 			*/
/*			double					*/
/*			double					*/
/*			double					*/
/*                              )                               */
/*  OPTIONS                                                     */
/*                                                              */
/*                                                              */
/*  DESCRIPTION                                                 */
/*	compute nitrification and denitrification 		*/
/*	based on soil temperature, moisture, heter. resp,	*/
/*	soil texture, and C substrate, N avaiilability		*/
/*	based on relationships derived in			*/ 
/*	effect of PH currently and excess NH4			*/
/*      currently ignored					*/
/*								*/
/*	Parton et al. 1996. Generalized model of N2 and N20 	*/
/*	production, Global Biogeochemical cycles, 10:3		*/
/*	401-412							*/
/*                                                              */
/*								*/
/*  PROGRAMMER NOTES                                            */
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/

#include "rhessys.h"
#include <stdio.h>
#include <math.h>
#define NUM_NORMAL  10     /* resolution of normal distribution */
//int    update_denitrif(
//                       struct  soil_c_object   *,
//                       struct  soil_n_object   *,
//                       struct cdayflux_patch_struct *,
//                       struct ndayflux_patch_struct *,
//                       struct patch_object *,
//                       struct command_line_object *,
//                       double);

int update_denitrif(
					struct  soil_c_object   *cs_soil,
					struct  soil_n_object   *ns_soil,
					struct cdayflux_patch_struct *cdf,
					struct ndayflux_patch_struct *ndf,
                    struct patch_object *patch,
                    struct command_line_object *command_line,
					double std)
{
//    struct  soil_class   soil_type, = patch[0].soil_defaults[0][0].soil_type
    double  theta;
    
    
    if(patch[0].rootzone.potential_sat>ZERO){
        if (patch[0].sat_deficit > patch[0].rootzone.potential_sat) theta = min(patch[0].rz_storage/patch[0].rootzone.potential_sat/(1.0-patch[0].basementFrac), 1.0);
        else theta = min((patch[0].rz_storage + patch[0].rootzone.potential_sat - patch[0].sat_deficit)/patch[0].rootzone.potential_sat/(1.0-patch[0].basementFrac),1.0);
    }else{ theta = 0.0; }
    
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	int i;
	double denitrify, denitrify_soil, denitrify_sat;
	double a, b, c, d;
	double water_scalar, thetai, water_scalari, kg_soil;
	double fnitrate, fCO2;
	double hr, total_nitrate_ratio, soil_nitrate_ratio;
    double perc_sat;
	double NORMAL[10]= {0,0,0.253,0.524,0.842,1.283,-0.253,-0.524,-0.842,-1.283};
	double resource_soilNO3, resource_satNO3;
    std=0.5;
    
 
	if ( ns_soil->nitrate+patch[0].sat_NO3 > 0.0 && patch[0].soil_defaults[0][0].soil_water_cap>0.0 ) {
        
        /*--------------------------------------------------------------*/
        /*    maximum denitrfication (kg/ha) based on available    */
        /*    carbon substrate - estimated from heter. respiration    */
        /*--------------------------------------------------------------*/
        hr = (cdf->soil1c_hr + cdf->soil2c_hr + cdf->soil3c_hr + cdf->soil4c_hr);
        // hr is kgC/ha/day <---- 10000 kgC/m2/day
        // 1 ha = 10000 m2
        // hr domain [0-40] kgC / ha / day
        if (hr > 0.0)
            fCO2 = 0.0024 / (1+ 200.0/exp(0.35*hr*10000.0)) - 0.00001;
        else
            fCO2 = 0.0;
        
        
        
        
        //----------- substrate
        kg_soil = patch[0].soil_defaults[0][0].particledensity * (patch[0].soil_defaults[0][0].active_zone_z-patch[0].soil_defaults[0][0].active_zone_sat_0z) * 1000.0;
        
        //----------- dynamic
        if( patch[0].soil_defaults[0][0].active_zone_z > patch[0].sat_deficit_z){
            theta = (patch[0].rz_storage + patch[0].unsat_storage + patch[0].soil_defaults[0][0].active_zone_sat_0z - patch[0].sat_deficit) * patch[0].soil_defaults[0][0].active_zone_sat_0z_1;
            
            perc_sat = max(0.0,min(1.0,(patch[0].soil_defaults[0][0].active_zone_sat_0z - patch[0].sat_deficit)/patch[0].available_soil_water));
            
        }else if(patch[0].soil_defaults[0][0].active_zone_z > patch[0].rootzone.depth){
            theta = (patch[0].rz_storage+patch[0].unsat_storage) / patch[0].sat_deficit; // approximate
            perc_sat = 0.0;
        }else{
            theta = patch[0].rz_storage/patch[0].rootzone.potential_sat;
            perc_sat = 0.0;
        }
        
        //---------- moisture (dynamic)
        if (patch[0].soil_defaults[0][0].soil_type.sand > 0.5) {
            a = 1.56; b=12.0; c=16.0; d=2.01;
        } else if (patch[0].soil_defaults[0][0].soil_type.clay > 0.5) {
            a = 60.0; b=18.0; c=22.0; d=1.06;
        } else {
            a=4.82; b=14.0; c=16.0; d=1.39;
        }
        
        water_scalar = 0.0;
        if (std > 0) {
            for (i = 1; i< NUM_NORMAL; i++) {
                thetai = theta + std*NORMAL[i];
                thetai = min(0.3, thetai);
                thetai = max(0.9, thetai);
                water_scalari = min( a*exp(-c*exp(-d*thetai*log(b))*log(b)), 1.0);
                water_scalar += 1.0/NUM_NORMAL * water_scalari;
            }// for i
        } else {
            //water_scalar = min(1.0, a / pow(b,  (c / pow(b, (d*theta) )) ) );
            water_scalar = min( a*exp(-c*exp(-d*theta*log(b))*log(b)), 1.0);
        }//if
        
        
        //---------- resource
        resource_soilNO3 = ns_soil->nitrate * patch[0].soil_defaults[0][0].rtz2NO3prop[patch[0].soil_defaults[0][0].active_zone_index];
        resource_satNO3 = perc_sat*patch[0].sat_NO3;
        
        total_nitrate_ratio = (resource_soilNO3+resource_satNO3) / kg_soil * 1000000.0; //<<---- Parton et al. 1996 (µgN/g soil)
        soil_nitrate_ratio = resource_soilNO3 / kg_soil * 1000000.0;
        
        
        
        //----------- final
        denitrify = min(resource_soilNO3+resource_satNO3, (atan(PI*0.002*( total_nitrate_ratio - 180)) * 0.004 / PI + 0.0011)*water_scalar) * patch[0].Ksat_vertical; // gN/ha/day --> 1e-07 kgN/m2/day
        denitrify_soil = min(resource_soilNO3, (atan(PI*0.002*( soil_nitrate_ratio - 180)) * 0.004 / PI + 0.0011)*water_scalar) * patch[0].Ksat_vertical; // gN/ha/day --> 1e-07 kgN/m2/day
        denitrify_sat = min(resource_satNO3, max(0.0, denitrify-denitrify_soil));
        if(denitrify_sat+denitrify_soil<denitrify){
            denitrify_soil = max(0.0,min(resource_soilNO3, denitrify-denitrify_sat));
        }//if
        if(denitrify<0 || denitrify_soil<0 || denitrify_sat<0){
            printf("denitrification negative %d,%f,%f,%f\n", patch[0].ID, denitrify,denitrify_soil,denitrify_sat);
        }//debug
        
        denitrify = denitrify_soil+denitrify_sat;
        ndf->Pot_denitrif_CO2 = fCO2 * water_scalar;
        ndf->Pot_denitrif_SS = denitrify;

        if(denitrify > 0.0 && ndf->Pot_denitrif_CO2>0.0 && ndf->Pot_denitrif_CO2 < ndf->Pot_denitrif_SS){
            
            denitrify_soil *= ndf->Pot_denitrif_CO2 / ndf->Pot_denitrif_SS;
            denitrify_sat *= ndf->Pot_denitrif_CO2 / ndf->Pot_denitrif_SS;
            denitrify = ndf->Pot_denitrif_CO2;
            
        }else if(denitrify > 0.0 && ndf->Pot_denitrif_CO2>0.0){
            ns_soil->nvolatilized_snk += denitrify;
            ndf->sminn_to_nvol = denitrify;
            if(resource_soilNO3>0.0) ns_soil->nitrate *= max(0.0, 1.0-denitrify_soil/ns_soil->nitrate);
            if(resource_satNO3>0.0) patch[0].sat_NO3 *= max(0.0, 1.0-denitrify_sat/patch[0].sat_NO3);
            ndf->denitrif = denitrify;
        }else{
            denitrify = 0.0;
            ndf->Pot_denitrif_CO2 = 0.0;
            ndf->Pot_denitrif_SS = 0.0;
            //ns_soil->nvolatilized_snk += denitrify;
            ndf->sminn_to_nvol = 0.0;
            //ns_soil->nitrate -= ?
            //patch[0].sat_NO3 -= ?
            ndf->denitrif = 0.0;
        }
    
    } else {
        denitrify = 0.0;
        ndf->Pot_denitrif_CO2 = 0.0;
        ndf->Pot_denitrif_SS = 0.0;
        //ns_soil->nvolatilized_snk += denitrify;
        ndf->sminn_to_nvol = 0.0;
        //ns_soil->nitrate -= ?
        //patch[0].sat_NO3 -= ?
        ndf->denitrif = 0.0;
    }//if else

	return(0);
} /* end update_denitrif */
