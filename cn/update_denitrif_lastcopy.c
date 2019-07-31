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
//    double  theta, = patch[0].rootzone.S
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	int ok,i, water_scalarCount;
	double denitrify, water_scalarSUM;
	double a, b, c, d;
	double water_scalar, theta, thetai, water_scalari, kg_soil, depth;
    double active_zone_z;
    double decay_rate;
	double fnitrate, fCO2;
	double hr, nitrate_ratio;
    #define  PARTICLE_DENSITY    2.65    /* soil particle density g/cm3 (Dingman) */
	#define NUM_NORMAL  10 	/* resolution of normal distribution */
	double NORMAL[10]= {0,0,0.253,0.524,0.842,1.283,-0.253,-0.524,-0.842,-1.283};
	
    double z1 = patch[0].sat_deficit_z>0? patch[0].sat_deficit_z : 0.0; //<<--- count for negative
    double z1_water = patch[0].sat_deficit_z>0? patch[0].sat_deficit : 0.0;
    double active_zone_z_sat = max(z1 + 0.33, patch[0].rootzone.depth);
    
    double p_decayRate = 1.0 / patch[0].soil_defaults[0][0].porosity_decay;
    double resource_soilNO3, resource_satNO3;
    double perc_rtz = 0.0;
    double perc_sat = 0.0;
    double tmp_ratio;
    
    std=0.5;
	ok = 1;
    active_zone_z = (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z));
    decay_rate = (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NO3decayRate : patch[0].soil_defaults[0][0].N_decay_rate);
    
    double constantHold1 = exp(-z1*p_decayRate);
    double constantHold2 = exp(-active_zone_z*decay_rate);
    double constantHold3 = exp(-active_zone_z_sat*p_decayRate);
    
    //depth = patch[0].rootzone.depth * 1.2;
    depth = min(patch[0].rootzone.depth, 0.6); //min(2.0, active_zone_z_sat);
    denitrify = 0.0;
    water_scalarSUM = 0.0;
    water_scalarCount = 0;
    fCO2 = 0.0;
    
    if( ns_soil->nitrate!=ns_soil->nitrate ){printf("ns_soil->nitrate is nan");}
	if ( depth>0 && (ns_soil->nitrate > 0.0) && patch[0].soil_defaults[0][0].soil_water_cap>0.0 ) {
        // there is nirate and positive soil depth
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
        

        
        if (patch[0].soil_defaults[0][0].soil_type.sand > 0.5) {
            a = 1.56; b=12.0; c=16.0; d=2.01;
        }
        else if (patch[0].soil_defaults[0][0].soil_type.clay > 0.5) {
            a = 60.0; b=18.0; c=22.0; d=1.06;
        }
        else {
            a=4.82; b=14.0; c=16.0; d=1.39;
        }
        
    
        // double active_zone_z_sat = max(z1 + 0.33, patch[0].rootzone.depth);
        // depth = min(2.0, active_zone_z_sat);
        if(z1>patch[0].rootzone.depth){
            
            // wtz > rtz --> unsat >0
            // above sat zone [0 - wtz]
            double layerZ = min(depth, patch[0].rootzone.depth);
            double thickZ = 0.1*layerZ;
            double tmpZ;
            theta = max(min(patch[0].rootzone.S, 1.0),0.0); //domin [0.3 - 0.9]
            water_scalar = 0.0;
            if (std > 0) {
                for (i =1; i< NUM_NORMAL; i++) {
                    thetai = theta + std*NORMAL[i];
                    thetai = min(0.9, thetai);
                    thetai = max(0.3, thetai);
                    if (thetai > ZERO)
                    water_scalari = min(1.0,a / pow(b,  (c / pow(b, (d*thetai) )) ));
                    water_scalar += 1.0/NUM_NORMAL * water_scalari;
                }
            }
            else
            water_scalar = min(1.0,a / pow(b,  (c / pow(b, (d*theta) )) ));
            
            fnitrate = 0.0;
            for(i=0; i<10; i++){
                tmpZ = thickZ*i;
                kg_soil = PARTICLE_DENSITY * (thickZ + patch[0].soil_defaults[0][0].porosity_0*patch[0].soil_defaults[0][0].porosity_decay * (exp(-(tmpZ+thickZ)*p_decayRate)- exp(-tmpZ*p_decayRate))) * 1000.0;
                resource_soilNO3 = (exp(-decay_rate*tmpZ)-exp(-decay_rate*(tmpZ+thickZ))) / (1.0 - constantHold2) * ns_soil->nitrate;
                nitrate_ratio = resource_soilNO3 / kg_soil * 1000000.0; //<<---- Parton et al. 1996 (µgN/g soil)
                fnitrate += min(resource_soilNO3, (atan(PI*0.002*( nitrate_ratio - 180)) * 0.004 / PI + 0.0011)*water_scalar); // gN/ha/day --> 1e-07 kgN/m2/day
            }// end of i loop
            denitrify *= (patch[0].Ksat_vertical+(1.0-patch[0].Ksat_vertical)*0.05);
            perc_rtz += denitrify;
            water_scalarSUM += water_scalar; water_scalarCount++;
            
            
            
            if(depth>patch[0].rootzone.depth){
                layerZ = min(depth, z1);
                thickZ = 0.2*(layerZ - patch[0].rootzone.depth);
                theta = patch[0].unsat_storage/(z1_water-patch[0].rootzone.potential_sat);
                water_scalar = 0.0;
                if (std > 0) {
                    for (i =1; i< NUM_NORMAL; i++) {
                        thetai = theta + std*NORMAL[i];
                        thetai = min(0.9, thetai);
                        thetai = max(0.3, thetai);
                        if (thetai > ZERO)
                        water_scalari = min(1.0,a / pow(b,  (c / pow(b, (d*thetai) )) ));
                        water_scalar += 1.0/NUM_NORMAL * water_scalari;
                    }
                }
                else
                water_scalar = min(1.0,a / pow(b,  (c / pow(b, (d*theta) )) ));
                
                fnitrate = 0.0;
                for(i=0; i<5; i++){
                    // (thickZ*(i+1)+patch[0].rootzone.depth)
                    tmpZ = (thickZ*i+patch[0].rootzone.depth);
                    kg_soil = PARTICLE_DENSITY * (thickZ + patch[0].soil_defaults[0][0].porosity_0*patch[0].soil_defaults[0][0].porosity_decay * (exp(-(tmpZ+thickZ)*p_decayRate)- exp(-tmpZ*p_decayRate))) * 1000.0;
                    resource_soilNO3 = (exp(-decay_rate*tmpZ)-exp(-decay_rate*(tmpZ+thickZ))) / (1.0 - constantHold2) * ns_soil->nitrate;
                    nitrate_ratio = resource_soilNO3 / kg_soil * 1000000.0; //<<---- Parton et al. 1996 (µgN/g soil)
                    fnitrate += min(resource_soilNO3, (atan(PI*0.002*( nitrate_ratio - 180)) * 0.004 / PI + 0.0011)*water_scalar); // gN/ha/day --> 1e-07 kgN/m2/day
                }// end of i loop
                denitrify *= (patch[0].Ksat_vertical+(1.0-patch[0].Ksat_vertical)*0.05);
                perc_rtz += denitrify;
                water_scalarSUM += water_scalar; water_scalarCount++;
            }//unsat
            
            if(depth>z1){
                layerZ = depth;
                thickZ = 0.2*(layerZ - z1);
                theta = 1.0;
                water_scalar = 0.0;
                if (std > 0) {
                    for (i =1; i< NUM_NORMAL; i++) {
                        thetai = theta + std*NORMAL[i];
                        thetai = min(0.9, thetai);
                        thetai = max(0.3, thetai);
                        if (thetai > ZERO)
                        water_scalari = min(1.0,a / pow(b,  (c / pow(b, (d*thetai) )) ));
                        water_scalar += 1.0/NUM_NORMAL * water_scalari;
                    }
                }
                else
                water_scalar = min(1.0,a / pow(b,  (c / pow(b, (d*theta) )) ));
                
                resource_soilNO3 = 0.0;
                resource_satNO3 = 0.0;
                double resource_soilNO3_;
                double resource_satNO3_;
                fnitrate = 0.0;
                for(i=0; i<5; i++){
                    // (thickZ*(i+1)+patch[0].rootzone.depth)
                    tmpZ = (thickZ*i+z1);
                    kg_soil = PARTICLE_DENSITY * (thickZ + patch[0].soil_defaults[0][0].porosity_0*patch[0].soil_defaults[0][0].porosity_decay * (exp(-(tmpZ+thickZ)*p_decayRate)- exp(-tmpZ*p_decayRate))) * 1000.0;
                    resource_soilNO3_ = (exp(-decay_rate*tmpZ)-exp(-decay_rate*(tmpZ+thickZ))) / (1.0 - constantHold2) * ns_soil->nitrate;
                    resource_satNO3_ = (exp(-p_decayRate*tmpZ)-exp(-p_decayRate*(tmpZ+thickZ)))/(constantHold1-constantHold3) * patch[0].sat_NO3;
                    resource_soilNO3 += resource_soilNO3_;
                    resource_satNO3 += resource_satNO3_;
                    nitrate_ratio = (resource_soilNO3_ + resource_satNO3_) / kg_soil * 1000000.0; //<<---- Parton et al. 1996 (µgN/g soil)
                    fnitrate += min(resource_soilNO3_ + resource_satNO3_, (atan(PI*0.002*( nitrate_ratio - 180)) * 0.004 / PI + 0.0011)*water_scalar); // gN/ha/day --> 1e-07 kgN/m2/day
                }// end of i loop
                denitrify *= (patch[0].Ksat_vertical+(1.0-patch[0].Ksat_vertical)*0.05);
                
                if(denitrify>0){
                    tmp_ratio = resource_soilNO3/(resource_soilNO3+resource_satNO3);
                    perc_rtz += min(resource_soilNO3, denitrify*tmp_ratio);
                    perc_sat += denitrify - min(resource_soilNO3, denitrify*tmp_ratio);
                }// partitioning
                water_scalarSUM += water_scalar; water_scalarCount++;
            }//satzone
            
        }else{
            // wtz<= rtz --> unsat =0
            
            // above sat zone [0 - wtz]
            double layerZ = min(depth, z1);
            double thickZ = 0.1*layerZ;
            double tmpZ;
            theta = max(min(patch[0].rootzone.S, 1.0),0.0);
            water_scalar = 0.0;
            if (std > 0) {
                for (i =1; i< NUM_NORMAL; i++) {
                    thetai = theta + std*NORMAL[i];
                    thetai = min(0.9, thetai);
                    thetai = max(0.3, thetai);
                    if (thetai > ZERO)
                    water_scalari = min(1.0,a / pow(b,  (c / pow(b, (d*thetai) )) ));
                    water_scalar += 1.0/NUM_NORMAL * water_scalari;
                }
            }
            else
            water_scalar = min(1.0,a / pow(b,  (c / pow(b, (d*theta) )) ));
            
            fnitrate = 0.0;
            for(i=0; i<10; i++){
                tmpZ = thickZ*i;
                kg_soil = PARTICLE_DENSITY * (thickZ + patch[0].soil_defaults[0][0].porosity_0*patch[0].soil_defaults[0][0].porosity_decay * (exp(-(tmpZ+thickZ)*p_decayRate)- exp(-tmpZ*p_decayRate))) * 1000.0;
                resource_soilNO3 = (exp(-decay_rate*tmpZ)-exp(-decay_rate*(tmpZ+thickZ))) / (1.0 - constantHold2) * ns_soil->nitrate;
                nitrate_ratio = resource_soilNO3 / kg_soil * 1000000.0; //<<---- Parton et al. 1996 (µgN/g soil)
                fnitrate += min(resource_soilNO3, (atan(PI*0.002*( nitrate_ratio - 180)) * 0.004 / PI + 0.0011)*water_scalar); // gN/ha/day --> 1e-07 kgN/m2/day
            }// end of i loop
            denitrify *= (patch[0].Ksat_vertical+(1.0-patch[0].Ksat_vertical)*0.05);
            perc_rtz += denitrify;
            water_scalarSUM += water_scalar; water_scalarCount++;
            
            if(depth > z1){
                layerZ = depth;
                thickZ = 0.2*(layerZ - z1);
                theta = 1.0;
                water_scalar = 0.0;
                if (std > 0) {
                    for (i =1; i< NUM_NORMAL; i++) {
                        thetai = theta + std*NORMAL[i];
                        thetai = min(0.9, thetai);
                        thetai = max(0.3, thetai);
                        if (thetai > ZERO)
                        water_scalari = min(1.0,a / pow(b,  (c / pow(b, (d*thetai) )) ));
                        water_scalar += 1.0/NUM_NORMAL * water_scalari;
                    }
                }
                else
                water_scalar = min(1.0,a / pow(b,  (c / pow(b, (d*theta) )) ));
                
                resource_soilNO3 = 0.0;
                resource_satNO3 = 0.0;
                double resource_soilNO3_;
                double resource_satNO3_;
                fnitrate = 0.0;
                for(i=0; i<5; i++){
                    // (thickZ*(i+1)+patch[0].rootzone.depth)
                    tmpZ = (thickZ*i+z1);
                    kg_soil = PARTICLE_DENSITY * (thickZ + patch[0].soil_defaults[0][0].porosity_0*patch[0].soil_defaults[0][0].porosity_decay * (exp(-(tmpZ+thickZ)*p_decayRate)- exp(-tmpZ*p_decayRate))) * 1000.0;
                    resource_soilNO3_ = (exp(-decay_rate*tmpZ)-exp(-decay_rate*(tmpZ+thickZ))) / (1.0 - constantHold2) * ns_soil->nitrate;
                    resource_satNO3_ = (exp(-p_decayRate*tmpZ)-exp(-p_decayRate*(tmpZ+thickZ)))/(constantHold1-constantHold3) * patch[0].sat_NO3;
                    resource_soilNO3 += resource_soilNO3_;
                    resource_satNO3 += resource_satNO3_;
                    nitrate_ratio = (resource_soilNO3_ + resource_satNO3_) / kg_soil * 1000000.0; //<<---- Parton et al. 1996 (µgN/g soil)
                    fnitrate += min(resource_soilNO3_ + resource_satNO3_, (atan(PI*0.002*( nitrate_ratio - 180)) * 0.004 / PI + 0.0011)*water_scalar); // gN/ha/day --> 1e-07 kgN/m2/day
                }// end of i loop
                denitrify *= (patch[0].Ksat_vertical+(1.0-patch[0].Ksat_vertical)*0.05);
                
                if(denitrify>0){
                    tmp_ratio = resource_soilNO3/(resource_soilNO3+resource_satNO3);
                    perc_rtz += min(resource_soilNO3, denitrify*tmp_ratio);
                    perc_sat += denitrify - min(resource_soilNO3, denitrify*tmp_ratio);
                }// partitioning
                water_scalarSUM += water_scalar; water_scalarCount++;
            }//satzone
        }
        
        
        
	} else denitrify = 0.0;
    
    denitrify = perc_sat + perc_rtz;
	/*--------------------------------------------------------------*/
	/*	update state and flux variables				*/
	/*--------------------------------------------------------------*/
    if(denitrify!=denitrify || isinf(denitrify)){
        printf("update_denitrif has infinite or nan problem {%e,%e,%e,%e,%e}\n",
               ns_soil->nitrate,
               depth,
               z1,
               patch[0].rootzone.depth,
               patch[0].rootzone.S);
    }//debug
    
    
    
    ndf->Pot_denitrif_CO2 = fCO2*water_scalarSUM/(1.0*water_scalarCount);
    ndf->Pot_denitrif_SS = denitrify;
    denitrify = min(ndf->Pot_denitrif_SS, ndf->Pot_denitrif_CO2);
    
	//denitrify = min(denitrify, ns_soil->nitrate);
	//denitrify = max(0.0, denitrify);
    
	ns_soil->nvolatilized_snk += denitrify;
	ndf->sminn_to_nvol = denitrify;
	ns_soil->nitrate -= min(ns_soil->nitrate, perc_rtz);
    patch[0].sat_NO3 -= min(patch[0].sat_NO3, perc_sat);
	ndf->denitrif = denitrify;
	ok = 0;
	return(ok);
} /* end update_denitrif */
