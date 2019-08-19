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

#define  PARTICLE_DENSITY	2.65	/* soil particle density g/cm3 (Dingman) */
//#define	 MAX_PERC		0.1	/* fraction of amonium that goes to nitrate */
//#define  MAX_RATE		(15*0.1353535)  //times by 6; 30 is max; balaning net mineralization
//old entry:  120	/* mgN/kg soil/day twice groffman values for ag soils */
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
	int ok,i;
    
    double MAX_RATE;
	double nitrify;
	double a, b, c, d;
	double nh4_conc, kg_soil;
	double N_scalar, water_scalar, T_scalar, pH_scalar;
	double theta, thetai;
	double max_nit_rate; /* kg/m2/day */
    double active_zone_z;
    double decay_rate;
    double depth;
    
    double z1 = patch[0].sat_deficit_z>0? patch[0].sat_deficit_z : 0.0; //<<--- count for negative
    double z1_water = patch[0].sat_deficit_z>0? patch[0].sat_deficit : 0.0;
	double active_zone_z_sat = max(z1 + 0.33, patch[0].rootzone.depth);
    
    active_zone_z = (command_line[0].NH4root2active>0.0? patch[0].rootzone.depth * command_line[0].NH4root2active : (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active>0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)));
    decay_rate = (command_line[0].NH4root2active>0.0? patch[0].soil_defaults[0][0].N_decay_rate : (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NH4decayRate : patch[0].soil_defaults[0][0].N_decay_rate));
    
    //depth = max( min( z1, patch[0].rootzone.depth)*1.2, 0.12);
    //depth = max( min( z1, patch[0].rootzone.depth), 1.0);
    depth = min(patch[0].rootzone.depth, 1.5);
    
    double p_decayRate = 1.0 / patch[0].soil_defaults[0][0].porosity_decay;
    double const_satexp_wtz = exp(-z1*p_decayRate);
    double const_satexp_act = exp(-active_zone_z_sat*p_decayRate);
    double const_satexp_dep = exp(-depth*p_decayRate);

    double const_soilexp_dep = exp(-depth*decay_rate);
    double const_soilexp_act = exp(-active_zone_z*decay_rate);
    double const_soilexp_wtz = exp(-z1*decay_rate);
    
    double resource_soilNH4, resource_satNH4;
    double perc_rtz = 0.0;
    double perc_sat = 0.0;
    double tmp_ratio;

    
    // depth = 0.12 (when fully saturated), 1.2rtz, or 1.2sat_def_z
    //**** case A: depth = 0.12, i.e., wtz ≤ 0.1 and rtz ≤ 0.1; no unsat
    //**** case B: depth = 1.2 rtz, i.e., wtz ≤ rtz; no unsat
    //**** case C: depth = 1.2 wtz, i.e., wtz > rtz; yes unsat
    //  --> sat (special case: depth = 1.2rtz = wtz) [wtz - depth] turns to a zero zone
    ok = 1;
    water_scalar = 0.0;
    nitrify = 0.0;
    std = 0.5;
    
    MAX_RATE = (15*0.1353535) * (1.0 + patch[0].aeratedSoilFrac*2.0);
    // Soil aeration should increase nitrification by THREE times?
    
	if( depth>0 && ns_soil->sminn + patch[0].sat_NH4 > 0.0 && patch[0].soil_defaults[0][0].soil_water_cap>0.0) {
        
       
        
        //N0 = N_decay_rate * total_nitrate / (1.0 - exp(-N_decay_rate * activedepthz));
        if (patch[0].soil_defaults[0][0].soil_type.sand > 0.5) {
            a = 0.55; b=1.7; c=-0.007; d=3.22;
        } else {
            a=0.6; b=1.27; c=0.0012; d=2.84;
        }
        
        T_scalar = min(-0.06 + 0.13 * exp(0.07 * patch[0].Tsoil),1.0); //domain [-5,30]
        pH_scalar = 0.56 + (atan(PI*0.45*(-5+patch[0].PH))/PI); // default 7.0, input by climate series.
        // acidic = less;  domain [3-7]
        
    
        if(z1>patch[0].rootzone.depth){
        	// wtz > rtz --> unsat>0
            // 0 - depth or rtz
            double layerZ = min(depth, patch[0].rootzone.depth);
        	theta = max(min(patch[0].rootzone.S, 1.0),0.0);
        	if (std > ZERO) {
                for (i=0; i<NUM_NORMAL; i++) {
                    thetai = theta + NORMAL[i]*std;
                    thetai = min(1.0, thetai);
                    thetai = max(0.002, thetai);
                    water_scalar  += 1.0/NUM_NORMAL * (
                                                       pow( ((thetai -b) / (a-b)), d*(b-a)/(a-c)) * pow( ((thetai-c)/ (a-c)), d) );
                }//for
            } else {
                if (theta  > c)
                    water_scalar  = pow( ((theta -b) / (a-b)), d*(b-a)/(a-c)) * pow( ((theta-c)/ (a-c)), d);
                else
                    water_scalar = 0.000001;
            }//if
            water_scalar = min(water_scalar,1.0);
            
            kg_soil = PARTICLE_DENSITY * (layerZ + patch[0].soil_defaults[0][0].porosity_0*patch[0].soil_defaults[0][0].porosity_decay * (exp(-layerZ*p_decayRate)-1.0)) * 1000.0; // correct, integrated from 0 to layer1
            max_nit_rate = kg_soil * MAX_RATE * 1e-06; // (mgN/kg soil/day) * (kg soil/m2) ==> mgN/m2/day ==> 1e-06 kgN/m2/day
            
            resource_soilNH4 = (1.0-exp(-layerZ*decay_rate)) / (1.0 - const_soilexp_act) * ns_soil->sminn;
            if(kg_soil>0){
                nh4_conc = resource_soilNH4 / kg_soil * 1000000.0; // (kgN/m2) / kg soil/m2) -> [conc.] (kgN/Kg soil)
                N_scalar = 1.0 - exp(-0.0105 * nh4_conc); // domain nh4_conc [0-50] ugN/gSoil
            }else{
                N_scalar = 0.0;
            }
            
            nitrify = min(resource_soilNH4, water_scalar * pH_scalar * T_scalar * N_scalar * max_nit_rate * (patch[0].Ksat_vertical+(1.0-patch[0].Ksat_vertical)*0.05));
            perc_rtz += nitrify;
        	
            //unsat = [rtz-wtz] -> depth = max(1,rtz) -> rtz < 1? <below> : <do nothing>
        	if(depth>patch[0].rootzone.depth){
                layerZ = min(depth, z1);
                theta = patch[0].unsat_storage/(z1_water-patch[0].rootzone.potential_sat);
                water_scalar = 0.0;
                if (std > ZERO) {
                    for (i=0; i<NUM_NORMAL; i++) {
                        thetai = theta + NORMAL[i]*std;
                        thetai = min(1.0, thetai);
                        thetai = max(0.002, thetai);
                        water_scalar  += 1.0/NUM_NORMAL * (
                                                           pow( ((thetai -b) / (a-b)), d*(b-a)/(a-c)) * pow( ((thetai-c)/ (a-c)), d) );
                    }//for
                } else {
                    if (theta  > c)
                    water_scalar  = pow( ((theta -b) / (a-b)), d*(b-a)/(a-c)) * pow( ((theta-c)/ (a-c)), d);
                    else
                    water_scalar = 0.000001;
                }//if
                water_scalar = min(water_scalar,1.0);
                
                kg_soil = PARTICLE_DENSITY * (layerZ-patch[0].rootzone.depth + patch[0].soil_defaults[0][0].porosity_0*patch[0].soil_defaults[0][0].porosity_decay * (exp(-layerZ*p_decayRate) - exp(-patch[0].rootzone.depth*p_decayRate) )) * 1000.0; // check equ *****
                max_nit_rate = kg_soil * MAX_RATE * 1e-06; // (mgN/kg soil/day) * (kg soil/m2) ==> mgN/m2/day ==> 1e-06 kgN/m2/day
                
                resource_soilNH4 = (exp(-decay_rate*patch[0].rootzone.depth)-exp(-decay_rate*layerZ)) / (1.0 - const_soilexp_act) * ns_soil->sminn;
                if(kg_soil>0){
                    nh4_conc = resource_soilNH4 / kg_soil * 1000000.0; // (kgN/m2) / kg soil/m2) -> [conc.] (kgN/Kg soil)
                    N_scalar = 1.0 - exp(-0.0105 * nh4_conc);
                }else{
                    N_scalar = 0.0;
                }
                
                nitrify = min(resource_soilNH4, water_scalar * pH_scalar * T_scalar * N_scalar * max_nit_rate * (patch[0].Ksat_vertical+(1.0-patch[0].Ksat_vertical)*0.05));
                perc_rtz += nitrify; //partitioning
        	}//unsat
        	
            //satzone = [wtz-depth] -> wtz=z1 > 1? <below> : <do nothing>
        	if(depth>z1){
                theta = 1.0;
                if (std > ZERO) {
                    for (i=0; i<NUM_NORMAL; i++) {
                        thetai = theta + NORMAL[i]*std;
                        thetai = min(1.0, thetai);
                        thetai = max(0.002, thetai);
                        water_scalar  += 1.0/NUM_NORMAL * (
                                                           pow( ((thetai -b) / (a-b)), d*(b-a)/(a-c)) * pow( ((thetai-c)/ (a-c)), d) );
                    }//for
                } else {
                    if (theta  > c)
                    water_scalar  = pow( ((theta -b) / (a-b)), d*(b-a)/(a-c)) * pow( ((theta-c)/ (a-c)), d);
                    else
                    water_scalar = 0.000001;
                }//if
                water_scalar = min(water_scalar,1.0);
                
                kg_soil = PARTICLE_DENSITY * (depth-z1 + patch[0].soil_defaults[0][0].porosity_0*patch[0].soil_defaults[0][0].porosity_decay * (const_satexp_dep-const_satexp_wtz)) * 1000.0; // correct, integrated from 0 to layer1
                max_nit_rate = kg_soil * MAX_RATE * 1e-06; // (mgN/kg soil/day) * (kg soil/m2) ==> mgN/m2/day ==> 1e-06 kgN/m2/day
                
                resource_soilNH4 = (const_soilexp_wtz-const_soilexp_dep) / (1.0 - const_soilexp_act) * ns_soil->sminn;
                resource_satNH4 = min(1.0,(const_satexp_wtz-const_satexp_dep)/(const_satexp_wtz-const_satexp_act)) * patch[0].sat_NH4;//can "depth" > activedepthz_sat? Yes
                // can "depth" < wtz? No, under the condition
                // condition: wtz<=rtz --> unsat=0
                // interval: satzone = [wtz-depth]
                // double active_zone_z_sat = max(z1 + 0.33, patch[0].rootzone.depth);
                // depth = max( min( z1, patch[0].rootzone.depth), 1.0);
                if(kg_soil>0){
                    nh4_conc = (resource_soilNH4+resource_satNH4) / kg_soil * 1000000.0; // (kgN/m2) / kg soil/m2) -> [conc.] (kgN/Kg soil)
                    N_scalar = 1.0 - exp(-0.0105 * nh4_conc);
                }else{
                    N_scalar = 0.0;
                }
                
                nitrify = min((resource_soilNH4+resource_satNH4), water_scalar * pH_scalar * T_scalar * N_scalar * max_nit_rate * (patch[0].Ksat_vertical+(1.0-patch[0].Ksat_vertical)*0.05));
                if(nitrify>0){
                    tmp_ratio = resource_satNH4/(resource_soilNH4+resource_satNH4);
                    perc_sat += nitrify * tmp_ratio;
                    perc_rtz += nitrify * (1.0-tmp_ratio);
                }// partitioning
        	}//satzone
        	// can "depth" < wtz? possible
        	// condition: wtz > rtz
        	// double active_zone_z_sat = max(z1 + 0.33, patch[0].rootzone.depth);
            // depth = max( min( z1, patch[0].rootzone.depth), 1.0);
        	
        }else{
        	// wtz<=rtz --> unsat=0 --> rtz > 1 ?
        	//rtz = [0-wtz] -> depth = 1 or wtz
            double layerZ = min(depth, z1);
            theta = min(1.0,(z1>0? patch[0].rz_storage/z1_water : 1.0));//<<----- problem
            if (std > ZERO) {
                for (i=0; i<NUM_NORMAL; i++) {
                    thetai = theta + NORMAL[i]*std;
                    thetai = min(1.0, thetai);
                    thetai = max(0.002, thetai);
                    water_scalar  += 1.0/NUM_NORMAL * (
                                                       pow( ((thetai -b) / (a-b)), d*(b-a)/(a-c)) * pow( ((thetai-c)/ (a-c)), d) );
                }//for
            } else {
                if (theta  > c)
                water_scalar  = pow( ((theta -b) / (a-b)), d*(b-a)/(a-c)) * pow( ((theta-c)/ (a-c)), d);
                else
                water_scalar = 0.000001;
            }//if
            water_scalar = min(water_scalar,1.0);
            //water_scalar = 1.0; //<<--------- debugging
            
            kg_soil = PARTICLE_DENSITY * (layerZ + patch[0].soil_defaults[0][0].porosity_0*patch[0].soil_defaults[0][0].porosity_decay * (exp(-layerZ*p_decayRate)-1.0)) * 1000.0; // correct, integrated from 0 to layer1
            max_nit_rate = kg_soil * MAX_RATE * 1e-06; // (mgN/kg soil/day) * (kg soil/m2) ==> mgN/m2/day ==> 1e-06 kgN/m2/day
            
            resource_soilNH4 = (1.0-exp(-layerZ*decay_rate)) / (1.0 - const_soilexp_act) * ns_soil->sminn;
            if(kg_soil>0){
                nh4_conc = resource_soilNH4 / kg_soil * 1000000.0; // (kgN/m2) / kg soil/m2) -> [conc.] (kgN/Kg soil)
                N_scalar = 1.0 - exp(-0.0105 * nh4_conc);
            }else{
                N_scalar = 0.0;
            }
            
            nitrify = min(resource_soilNH4, water_scalar * pH_scalar * T_scalar * N_scalar * max_nit_rate * (patch[0].Ksat_vertical+(1.0-patch[0].Ksat_vertical)*0.05));
            perc_rtz += nitrify;
            
        	
        	//sat = [wtz-depth]
            if(depth > z1){
                theta = 1.0;
                if (std > ZERO) {
                    for (i=0; i<NUM_NORMAL; i++) {
                        //NORMAL[10]= {0,0,0.253,0.524,0.842,1.283,-0.253,-0.524,-0.842,-1.283};
                        thetai = theta + NORMAL[i]*std;
                        thetai = min(1.0, thetai);
                        thetai = max(0.002, thetai);
                        water_scalar  += 1.0/NUM_NORMAL * (
                                                           pow( ((thetai -b) / (a-b)), d*(b-a)/(a-c)) * pow( ((thetai-c)/ (a-c)), d) );
                    }//for
                } else {
                    if (theta  > c)
                        water_scalar  = pow( ((theta -b) / (a-b)), d*(b-a)/(a-c)) * pow( ((theta-c)/ (a-c)), d);
                    else
                        water_scalar = 0.000001;
                }//if
                water_scalar = min(water_scalar,1.0);
                //water_scalar = 1.0; //<<--------- debugging
                
                kg_soil = PARTICLE_DENSITY * (depth-z1 + patch[0].soil_defaults[0][0].porosity_0*patch[0].soil_defaults[0][0].porosity_decay * (const_satexp_dep-const_satexp_wtz)) * 1000.0; // correct, integrated from 0 to layer1
                max_nit_rate = kg_soil * MAX_RATE * 1e-06; // (mgN/kg soil/day) * (kg soil/m2) ==> mgN/m2/day ==> 1e-06 kgN/m2/day
                
                resource_soilNH4 = (const_soilexp_wtz-const_soilexp_dep) / (1.0 - const_soilexp_act) * ns_soil->sminn;
                resource_satNH4 = min(1.0,(const_satexp_wtz-const_satexp_dep)/(const_satexp_wtz-const_satexp_act)) * patch[0].sat_NH4;//can "depth" > activedepthz_sat? Yes
                    // can "depth" < wtz? No, under the condition
                    // condition: wtz<=rtz --> unsat=0
                    // interval: satzone = [wtz-depth]
                    // double active_zone_z_sat = max(z1 + 0.33, patch[0].rootzone.depth);
                    // depth = max( min( z1, patch[0].rootzone.depth), 1.0);
                if(kg_soil>0){
                    nh4_conc = (resource_soilNH4+resource_satNH4) / kg_soil * 1000000.0; // (kgN/m2) / kg soil/m2) -> [conc.] (kgN/Kg soil)
                    N_scalar = 1.0 - exp(-0.0105 * nh4_conc);
                }else{
                    N_scalar = 0.0;
                }
                //N_scalar = 1.0; //<<--------- debugging
                
                nitrify = min((resource_soilNH4+resource_satNH4), water_scalar * pH_scalar * T_scalar * N_scalar * max_nit_rate * (patch[0].Ksat_vertical+(1.0-patch[0].Ksat_vertical)*0.05));
                if(nitrify>0){
                    tmp_ratio = resource_satNH4/(resource_soilNH4+resource_satNH4);
                    perc_sat += nitrify * tmp_ratio;
                    perc_rtz += nitrify * (1.0-tmp_ratio);
                }// partitioning
            }//satzone
        }//end of ifelse
        nitrify = perc_sat + perc_rtz;
        
        if(nitrify!=nitrify || isinf(nitrify)){
            printf("update_nitrif has infinite or nan problem [%d]{%e,%e,%e,%e,%e}\n",
                   patch[0].ID,
                   ns_soil->sminn,
                   depth,
                   z1,
                   patch[0].rootzone.depth,
                   patch[0].rootzone.S);
        }//debug
        
  
        //nitrify = max(min(nitrify, ns_soil->sminn),0.0);

        ndf->sminn_to_nitrate = nitrify;
        ns_soil->sminn -= min(ns_soil->sminn, perc_rtz);
        patch[0].sat_NH4 -= min(patch[0].sat_NH4, perc_sat);
        ns_soil->nitrate += perc_rtz;
        patch[0].sat_NO3 += perc_sat;
        //kg_soil > =0
    }else{
        //kg_soil ==0
        ndf->sminn_to_nitrate = 0.0;
    }
    ok = 0;
	return(ok);
} /* end update_nitrif */

