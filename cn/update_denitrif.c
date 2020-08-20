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
//#define NUM_NORMAL  10     /* resolution of normal distribution */
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
        if (patch[0].sat_deficit > patch[0].rootzone.potential_sat) theta = min(patch[0].rz_storage/patch[0].rootzone.potential_sat, 1.0);//(1.0-patch[0].basementFrac)
        else theta = min((patch[0].rz_storage + patch[0].rootzone.potential_sat - patch[0].sat_deficit)/patch[0].rootzone.potential_sat,1.0);//(1.0-patch[0].basementFrac)
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
	double NORMAL[10]= { 0.0, -1.283,-0.842,-0.524,-0.253, 0.0, 0.253,0.524,0.842,1.283};
    //double NORMAL_GROWSEASON[50]= { 0.0, -1.9245,-1.8589,-1.6667,-1.3608,-1.263,-1.22,-1.0938,-0.9622,-0.8931,-0.786,-0.7592,-0.6807,-0.6315,-0.5558,-0.4981,-0.393,-0.3795,-0.3666,-0.3287,-0.3269,-0.2683,-0.2034,-0.1897,-0.0982,0.0,0.0982,0.1897,0.2034,0.2683,0.3269,0.3287,0.3666,0.3795,0.393,0.4981,0.5558,0.6315,0.6807,0.7592,0.786,0.8931,0.9622,1.0938,1.22,1.263,1.3608,1.6667,1.8589,1.9245};
    //double NORMAL_DORMSEASON[50]= { 0.0, -1.283,-1.2393,-1.1111,-0.9072,-0.842,-0.8133,-0.7292,-0.6415,-0.5954,-0.524,-0.5061,-0.4538,-0.421,-0.3705,-0.3321,-0.262,-0.253,-0.2444,-0.2191,-0.2179,-0.1789,-0.1356,-0.1265,-0.0655,0,0.0655,0.1265,0.1356,0.1789,0.2179,0.2191,0.2444,0.253,0.262,0.3321,0.3705,0.421,0.4538,0.5061,0.524,0.5954,0.6415,0.7292,0.8133,0.842,0.9072,1.1111,1.2393,1.283};
    //double NORMAL_WEIGHT[50] = {0.0, 2.0,4.0,4.0,4.0,2.0,4.0,4.0,4.0,4.0,2.0,4.0,4.0,4.0,4.0,4.0,4.0,2.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,40.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,2.0,4.0,4.0,4.0,4.0,4.0,4.0,2.0,4.0,4.0,4.0,4.0,2.0,4.0,4.0,4.0,2};
	double resource_soilNO3, resource_satNO3;
    
 
	if ( ns_soil->nitrate+patch[0].sat_NO3 > 0.0 && patch[0].soil_defaults[0][0].soil_water_cap>0.0 ) {
        
        /*--------------------------------------------------------------*/
        /*    maximum denitrfication (kg/ha) based on available    */
        /*    carbon substrate - estimated from heter. respiration    */
        /*--------------------------------------------------------------*/
        hr = (cdf->soil1c_hr + cdf->soil2c_hr+cdf->soil3c_hr); // need to be at this 0.0001 scale to be sense; take out respiration from soil3 and soil4 because these are very slow respiration, which cnanot be capture by the lab incubation.
        hr *= patch[0].soil_defaults[0][0].active_zone_omProp;
        // if(hr > 0.0004670731) printf("patch[%d] hr %e [%e %e %e %e]\n",patch[0].ID, hr, cdf->soil1c_hr , cdf->soil2c_hr , cdf->soil3c_hr , cdf->soil4c_hr);
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
            // change here
            
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
        if (std > 0 ) {
            for (i = 1; i<10; i++) {
                thetai = theta + std*NORMAL[i];
                thetai = min(1.0, thetai);
                thetai = max(0.002, thetai);
                water_scalari = min( a*exp(-c*exp(-d*thetai*log(b))*log(b)), 1.0);
                water_scalar += water_scalari;
            }// for i
            water_scalar *= 0.1111111; // 1/9
        } else {
            //water_scalar = min(1.0, a / pow(b,  (c / pow(b, (d*theta) )) ) );
            water_scalar = min( a*exp(-c*exp(-d*theta*log(b))*log(b)), 1.0);
        }//if
        
//        if (std > 0 && ) {
//            for (i = 1; i<50; i++) {
//                thetai = theta + std*NORMAL_GROWSEASON[i];
//                thetai = min(1.0, thetai);
//                thetai = max(0.002, thetai);
//                water_scalari = min( a*exp(-c*exp(-d*thetai*log(b))*log(b)), 1.0);
//                water_scalar += water_scalari * NORMAL_WEIGHT[i];
//            }// for i
//            water_scalar *= 0.00462963; // 1/216
//        }else if (std > 0) {
//            for (i = 1; i<50; i++) {
//                thetai = theta + std*NORMAL_DORMSEASON[i];
//                thetai = min(1.0, thetai);
//                thetai = max(0.002, thetai);
//                water_scalari = min( a*exp(-c*exp(-d*thetai*log(b))*log(b)), 1.0);
//                water_scalar += water_scalari * NORMAL_WEIGHT[i];
//            }// for i
//            water_scalar *= 0.00462963; // 1/216
//        } else {
//            //water_scalar = min(1.0, a / pow(b,  (c / pow(b, (d*theta) )) ) );
//            water_scalar = min( a*exp(-c*exp(-d*theta*log(b))*log(b)), 1.0);
//        }//if
        
        
        //---------- resource
        resource_soilNO3 = ns_soil->nitrate * patch[0].soil_defaults[0][0].rtz2NO3prop[patch[0].soil_defaults[0][0].active_zone_index];
        resource_satNO3 = perc_sat*patch[0].sat_NO3;
        
        total_nitrate_ratio = (resource_soilNO3+resource_satNO3) / kg_soil * 1000000.0; //<<---- Parton et al. 1996 (µgN/g soil)
        soil_nitrate_ratio = resource_soilNO3 / kg_soil * 1000000.0;
        
        
        
        //----------- final
        // (patch[0].drainage_type == STREAM?, 0.0 : 1.0);
        double MaxDenitrify =
            max(0.0,patch[0].Ksat_vertical-patch[0].aeratedSoilFrac) * 5e-05 + //4e-06 kgN/m2/day Vermes and Myrold et al., 1991 for forest soil in Oregon (max of 1.150685e-05 kgN/m2/day from many forest denitrification studies); //max 5e-05 kgN/m2/day Groffman and Tiedje 1989a forest soil in MI; max 1.150685e-05 kgN/m2/day Groffman and Tiedje 1989b
            patch[0].aeratedSoilFrac * (theta*patch[0].soil_defaults[0][0].sat_def_0zm[patch[0].soil_defaults[0][0].active_zone_index] >0.35? 6.134247e-05 : 8.082192e-07);// 6.134247e-05 kgN/m2/day  Raciti et al., 2011
        
        // MaxDenitrify *= (patch[0].drainage_type>0 && patch[0].drainage_type % actionRIPARIAN==0? 1.7:1.0); // new developing wetland feature
        
        // Parton equation is based on lab incubation at optimized condition. serves max denitrification rate
        denitrify = min(MaxDenitrify, min(resource_soilNO3+resource_satNO3, (atan(PI*0.002*( total_nitrate_ratio - 180)) * 0.004 / PI + 0.0011)*water_scalar)) * patch[0].Ksat_vertical; // gN/ha/day --> 1e-07 kgN/m2/day
        denitrify_soil = min(MaxDenitrify,min(resource_soilNO3, (atan(PI*0.002*( soil_nitrate_ratio - 180)) * 0.004 / PI + 0.0011)*water_scalar)) * patch[0].Ksat_vertical; // gN/ha/day --> 1e-07 kgN/m2/day
        
        
       
        denitrify_sat = min(resource_satNO3, max(0.0, denitrify-denitrify_soil));
        if(denitrify_sat+denitrify_soil<denitrify){
            denitrify_soil = max(0.0,min(resource_soilNO3, denitrify-denitrify_sat));
        }//if
        if(denitrify<0 || denitrify_soil<0 || denitrify_sat<0){
            printf("denitrification negative %d,%f,%f,%f\n", patch[0].ID, denitrify,denitrify_soil,denitrify_sat);
        }//debug
        
        denitrify = denitrify_soil + denitrify_sat;
        ndf->Pot_denitrif_CO2 = fCO2 * water_scalar; // respiration
        ndf->Pot_denitrif_SS = denitrify; // resource

        if(ndf->Pot_denitrif_SS>0.0 && ndf->Pot_denitrif_CO2>0.0 && ndf->Pot_denitrif_CO2 < ndf->Pot_denitrif_SS){
            // limited by respiration factor
            denitrify_soil *= ndf->Pot_denitrif_CO2 / ndf->Pot_denitrif_SS;
            denitrify_sat *= ndf->Pot_denitrif_CO2 / ndf->Pot_denitrif_SS;
            denitrify = ndf->Pot_denitrif_CO2;
            
            ns_soil->nvolatilized_snk += denitrify;
            ndf->sminn_to_nvol = denitrify;
            if(resource_soilNO3>0.0) ns_soil->nitrate *= max(0.0, 1.0-denitrify_soil/ns_soil->nitrate);
            if(resource_satNO3>0.0) patch[0].sat_NO3 *= max(0.0, 1.0-denitrify_sat/patch[0].sat_NO3);
            ndf->denitrif = denitrify;
            
        }else if(ndf->Pot_denitrif_SS>0.0 && ndf->Pot_denitrif_CO2>0.0){
            //limited by resource
            ns_soil->nvolatilized_snk += denitrify;
            ndf->sminn_to_nvol = denitrify;
            if(resource_soilNO3>0.0) ns_soil->nitrate *= max(0.0, 1.0-denitrify_soil/ns_soil->nitrate);
            if(resource_satNO3>0.0) patch[0].sat_NO3 *= max(0.0, 1.0-denitrify_sat/patch[0].sat_NO3);
            ndf->denitrif = denitrify;
        }else{
            denitrify = 0.0;
            ndf->Pot_denitrif_CO2 = 0.0;
            ndf->Pot_denitrif_SS = 0.0;
            ndf->sminn_to_nvol = 0.0;
            ndf->denitrif = 0.0;
        }
    
    } else {
        // no nitrate
        denitrify = 0.0;
        ndf->Pot_denitrif_CO2 = 0.0;
        ndf->Pot_denitrif_SS = 0.0;
        ndf->sminn_to_nvol = 0.0;
        ndf->denitrif = 0.0;
    }//if else

    if(ndf->denitrif<ndf->Pot_denitrif_CO2 && ndf->denitrif<ndf->Pot_denitrif_SS){
        printf("denitrification prob %d,%f,%f, %f\n", patch[0].ID,denitrify_soil,denitrify_sat,
               ndf->denitrif, ndf->Pot_denitrif_CO2, ndf->Pot_denitrif_SS );
    }
    
    
	return(0);
} /* end update_denitrif */
