/*--------------------------------------------------------------*/
/*                                                              */
/*		update_nitrif				*/
/*                                                              */
/*  NAME                                                        */
/*		update_nitrif				*/
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/*  void update_nitrif(				*/
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
/*	compute nitrification and nitrification 		*/
/*	based on soil temperature, moisture, heter. resp,	*/
/*	soil texture, and C substrate, N avaiilability		*/
/*	based on relationships derived in			*/ 
/*								*/
/*	effect of pH from Parton et al 2004			*/
/*      pH equation from tropical acidic soils                  */
/*      effect of excess NH4                                    */
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
#include "phys_constants.h"
#include <stdio.h>
#include <math.h>

//#define  PARTICLE_DENSITY	2.65	/* soil particle density g/cm3 (Dingman) */
//#define	 MAX_PERC		0.1	/* fraction of amonium that goes to nitrate */
//#define  MAX_RATE		120    /* mgN/kg/day twice groffman values for ag soils */
#define NUM_NORMAL  10 	/* resolution of normal distribution */
double NORMAL[10]= {0,0,0.253,0.524,0.842,1.283,-0.253,-0.524,-0.842,-1.283};


int update_nitrif(
				  struct  soil_c_object   *cs_soil,
				  struct  soil_n_object   *ns_soil,
				  struct cdayflux_patch_struct *cdf,
				  struct ndayflux_patch_struct *ndf,
				  struct patch_object *patch,
                  struct command_line_object *command_line,
				  double std)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	int i;
    double MAX_RATE;
	double nitrify_total, nitrify_soil, nitrify_sat;
	double a, b, c, d;
	double nh4_conc, kg_soil;
	double N_scalar_total, N_scalar_soil, water_scalar, T_scalar, pH_scalar;
	double theta, thetai;
	double max_nit_rate; /* kg/m2/day */
    double perc_sat;
    double resource_soilNH4, resource_satNH4;
    
    MAX_RATE = 1.0 * (1.0 + 59.0*patch[0].aeratedSoilFrac);
    // 0.5606 from coweeta forest (mgN/kg soil/day) WS7 with 8.9 mg NH4-N / kg soil (Coweeta forest book)
    // extract from Table 1 on page 75, the guessed max nitrificaiton in forest is 1 mgN/kg soil/day.
    std = 0.5;
   
   
    
	if( patch[0].soil_defaults[0][0].active_zone_z>0 && ns_soil->sminn + patch[0].sat_NH4 > 0.0) {
        
        // ----- Temperature ad pH factors
        if (patch[0].soil_defaults[0][0].soil_type.sand > 0.5) {
            a = 0.55; b=1.7; c=-0.007; d=3.22;
        } else {
            a=0.6; b=1.27; c=0.0012; d=2.84;
        }// if
        T_scalar = min(-0.06 + 0.13 * exp(0.07 * patch[0].Tsoil),1.0);
        pH_scalar = 0.56 + (atan(PI*0.45*(-5+patch[0].PH))/PI); // default 7.0, input by climate series.
        // forest may be lower in pH (should look from SSURGO)
        
        //----------- substrate
        kg_soil = patch[0].soil_defaults[0][0].particledensity * (patch[0].soil_defaults[0][0].active_zone_z-patch[0].soil_defaults[0][0].active_zone_sat_0z) * 1000.0; //within active zone
        max_nit_rate = kg_soil * MAX_RATE * 1e-06; // (mgN/kg soil/day) * (kg soil/m2) ==> mgN/m2/day ==> 1e-06 kgN/m2/day
        
        
        
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
        
        //----------- moisture in active_zone (dynamic)
        if (std > ZERO) {
            for (i=0; i<NUM_NORMAL; i++) {
                thetai = theta + NORMAL[i]*std;
                thetai = min(1.0, thetai);
                thetai = max(0.002, thetai);
                water_scalar  += 1.0/NUM_NORMAL * exp(d*(b-a)/(a-c)*log((thetai -b) / (a-b))) * exp(d*log((thetai-c)/ (a-c)));
            }//for
            
        } else {
            if (theta  > c)
                water_scalar  = exp(d*(b-a)/(a-c)*log((thetai -b) / (a-b))) * exp(d*log((thetai-c)/ (a-c)));
            else
                water_scalar = 0.000001;
        }//if
        water_scalar = min(water_scalar,1.0);
        
        //----------- resource NH4 (soil and sat) in active_zone (dynamic)
        resource_soilNH4 = ns_soil->sminn * patch[0].soil_defaults[0][0].rtz2NH4prop[patch[0].soil_defaults[0][0].active_zone_index];
        resource_satNH4 = perc_sat*patch[0].sat_NH4;
        
        //nh4_conc = 1000000.0 * (resource_soilNH4+resource_satNH4) / kg_soil; // (kgN/m2) / kg soil/m2) -> [conc.] (kgN/Kg soil)
        N_scalar_total = 1.0 - exp(-0.0105 * 1000000.0 * (resource_soilNH4+resource_satNH4) / kg_soil); // domain nh4_conc [0-50] ugN/gSoil
        N_scalar_soil = 1.0 - exp(-0.0105 * 1000000.0 * resource_soilNH4 / kg_soil);
    
        //----------- final
        nitrify_total = min(resource_soilNH4+resource_satNH4, water_scalar * pH_scalar * T_scalar * N_scalar_total * max_nit_rate * patch[0].Ksat_vertical); // boundary by imperivous
        nitrify_soil = min(resource_soilNH4, water_scalar * pH_scalar * T_scalar * N_scalar_soil * max_nit_rate * patch[0].Ksat_vertical);
        nitrify_sat = min(resource_satNH4, max(0.0, nitrify_total-nitrify_soil));
        if(nitrify_soil + nitrify_sat < nitrify_total){
            nitrify_soil = max(0.0,min(resource_soilNH4, nitrify_total-nitrify_sat));
        }//if
        
        if(nitrify_total!=nitrify_total || isinf(nitrify_total) || nitrify_total<0 || nitrify_soil<0 || nitrify_sat<0 || resource_satNH4<0){
            printf("update_nitrif has infinite or nan problem [%d]{%e(%e),%e,%e(%e[%e %e %e] %e),%e,%e, %e %e %e}\n",
                   patch[0].ID,
                   ns_soil->sminn, patch[0].soil_defaults[0][0].rtz2NH4prop[patch[0].soil_defaults[0][0].active_zone_index],
                   nitrify_soil,
                   nitrify_sat, perc_sat,patch[0].soil_defaults[0][0].active_zone_sat_0z, patch[0].sat_deficit,patch[0].available_soil_water, patch[0].sat_NH4, // perc_sat being negative
                   patch[0].rootzone.depth,
                   patch[0].rootzone.SatPct, nitrify_total, nitrify_soil, nitrify_sat);
        }//debug
        
  
        //nitrify = max(min(nitrify, ns_soil->sminn),0.0);
        nitrify_total = nitrify_soil + nitrify_sat;
        ndf->sminn_to_nitrate = nitrify_total;
        
        if(resource_soilNH4>0.0) ns_soil->sminn *= max(0.0, 1.0-nitrify_soil/ns_soil->sminn);
        ns_soil->nitrate += nitrify_soil;
        
        if(resource_satNH4>0.0) patch[0].sat_NH4 *= max(0.0, 1.0-nitrify_sat/patch[0].sat_NH4);
        patch[0].sat_NO3 += nitrify_sat;
        //kg_soil > =0
    }else{
        //kg_soil ==0
        ndf->sminn_to_nitrate = 0.0;
    }
	return(0);
} /* end update_nitrif */

