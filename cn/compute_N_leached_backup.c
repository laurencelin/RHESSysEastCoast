/*--------------------------------------------------------------*/
/*                                                              */
/*		compute_N_leached				*/
/*                                                              */
/*  NAME                                                        */
/*		compute_N_leached				*/
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/*  void compute_N_leached(int					*/
/*					double	,		*/
/*					double	,		*/
/*					double	,		*/
/*					double	,		*/
/*					double	,		*/
/*					double	,		*/
/*					double	,		*/
/*					double	);		*/
/*                                                              */
/*  OPTIONS                                                     */
/*                                                              */
/*                                                              */
/*  DESCRIPTION                                                 */
/*                                                              */
/*								*/
/*  PROGRAMMER NOTES                                            */
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#define  PARTICLE_DENSITY    2.65    /* soil particle density g/cm3 (Dingman) */

double	compute_N_leached(
            int verbose_flag,
			double total_nitrate, 
			double Qout, 
			double satdef, //satdef --> 0
			double soil_water_capacity, //soil_water_cap --> 0
			double ndistz, // no use
			double patchID, // no use
			double p_0, //porosity_0
			double p, // porosity_decay in depth
			double N_decay_rate, // solute vertical decay
			double activedepthz, //active_z (real depth in m)
			double soildepthz, //soil_depth (real depth in m)
			double N_absorption_rate,
			double *transmissivity) 
			
	{ 
	/*------------------------------------------------------*/ 
	/*	Local Function Declarations.						*/ 
	/*------------------------------------------------------*/
    	double  compute_delta_water(
                int,
                double,
                double,
                double,
                double,
                double);


	double  compute_z_final(
		int,
		double,
		double,
		double,
		double,
		double);
	double  compute_N_absorbed(
		int,
		double,
		double,
		double,
		double,
		double);
		
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	
    double nleached = 0.0;
	//double theta, sat_deficit;
	double z1 = 0.0, z2 = 0.0;
    
    double availableQ = 0.0;
    double N0 = 0.0;
    double absorptionConst = 0.0;
    double critialZ = 0.0;
	//double septic_depth;

	
	

	/*------------------------------------------------------*/
	/* nitrate export only occurs when Qout > 0.0		*/ 
	/*------------------------------------------------------*/
    if (Qout > 0 && total_nitrate>0 && activedepthz>0) {
        if (satdef < 0.0) satdef = 0.0;
        if (soil_water_capacity < satdef) soil_water_capacity = satdef;
        /*------------------------------------------------------*/
        /*    first look at the case of return flow        */
        /*    for return flow we must estimate the sat_deficit */
        /*    that would account for the flow            */
        /*    (assuming all water leaves, so Qout/theta here */
        /*    is 1)                        */
        /*------------------------------------------------------*/
        if ((satdef == 0.0) && (soil_water_capacity == 0.0)) {

            //this is for returnflow
            //z2 = -1.0 * p * log (1 - (Qout) / (p * p_0)); //printf("compute_N_leached %e, %e, %e\n",z2,activedepthz,N_decay_rate);
            z2 = soildepthz;
            z1 = 0.0;


        } else {
            /*------------------------------------------------------*/
            /*    now for regular subsurface flow            */
            /*    integrate through the saturated zone        */
            /*------------------------------------------------------*/
            z2 = soildepthz; //soil depth
            z1 = compute_z_final(
                    verbose_flag,
                    p_0,
                    p,
                    soildepthz, //soil_depth
                    0.0, // inital depth
                    -satdef); //satdef --> sat_def_z

        }//end of else

        // available water below water table
        availableQ = compute_delta_water(
                                         verbose_flag,
                                         p_0, //<<----- porosity P_0
                                         p, // porosity decay (m)
                                         soildepthz, //soil_depth
                                         z2,  //<<--- soil depth (ini)
                                         z1); //<<--- sat_def_z (end) ndistz

        
        if(availableQ<=0){
            printf("N_leached{%e,%e}->[%e,%e|%e]->{%e,%e}\n",
                   satdef,soil_water_capacity,
                   Qout,availableQ,(Qout/availableQ),
                   z1,z2);
        }
        
//        if(N_decay_rate == 0.0 || absorptionConst == 0.0){
//            // what happened if N_decay_rate is zero?
//           // nleached = (Qout/availableQ) * max(total_nitrate - absorptionConst*(z2-z1),0.0);
//            nleached = 0.0000001 * max(total_nitrate - absorptionConst*(z2-z1),0.0);// what is this?
//            //if(nleached!=nleached) nleached = 0.0;
//
//        }else
        
        
        
        if(N_decay_rate>0){
            //---------------------------------------------------------- N decay distribution
            N0 = N_decay_rate * total_nitrate / (1.0 - exp(-N_decay_rate * activedepthz));
            absorptionConst = p_0*(1.0-p_0)*PARTICLE_DENSITY * 1000.0 * N_absorption_rate;
            critialZ =  absorptionConst / N0; //assume p (porosity decay - meter) is large 4000 m (default)
           
            if(critialZ>0 && critialZ < 1.0 && absorptionConst>0 && N0 > 0){
                critialZ = -log(critialZ)/N_decay_rate;
                if(critialZ > z1){ //sat_def_z
                    // integrating from sat_def_z(=z1) to critialZ; approaximation when (porosity decay - meter) is large 4000 m (default)
                    nleached = (Qout/availableQ) * (absorptionConst*(z1-critialZ) + N0/N_decay_rate*(exp(-N_decay_rate*z1)-exp(-N_decay_rate*critialZ)) );
                    nleached = max(min(nleached,total_nitrate), 0.0);
                    //if(nleached!=nleached) nleached = 0.0;
                }else{
                    //nothing comes out; absorptionConst is too great
                    nleached = 0.0;
                }
            }else if(N0 > 0 && !(absorptionConst>0) ){
                // NO3
                // critialZ=inf; absorptionConst=0
                // integrating from sat_def_z(=z1) to critialZ(=inf)
                nleached = (Qout/availableQ) * total_nitrate * (exp(-N_decay_rate*(z1+ndistz))-exp(-N_decay_rate*(activedepthz+ndistz))) / (1.0 - exp(-N_decay_rate * (activedepthz+ndistz)));
                
                nleached = max(min(nleached,total_nitrate), 0.0);
                
            }else{ nleached = 0.0; }// if critialZ
            
        }else{
            //---------------------------------------------------------- N uniform distribution
            //N0 = total_nitrate/activedepthz;
            //absorptionConst = p_0*(1.0-p_0)*PARTICLE_DENSITY * 1000.0 * N_absorption_rate;
            //critialZ = ?
            // solve N0 - #
            // not going to happen in this case; just assign to zero for now.
            nleached = 0.0;
        }// decay
        
    }else{nleached = 0.0;}//end of Qout > ZERO
        
	
	/*------------------------------------------------------*/
	/* there may be enough flow to leach out more than 	*/
	/*	availabe nitrate, so limit export by available	*/
	/*------------------------------------------------------*/

    //if (nleached > total_nitrate){nleached = total_nitrate;}
        
        
    if(nleached<0 || nleached!=nleached) printf("leaching[%f,%f]: ->(%e,%e,%e), N0(%e)=(%e), absorptionConst(%e)=(%e), critialZ_ini(%e), critialZ(%e)>(%e), nleached(%e)<(%e)\n",
           patchID,ndistz,
           Qout,availableQ,Qout/availableQ,
           N0, N_decay_rate,
           N_absorption_rate, absorptionConst,
           absorptionConst/N0,
           critialZ, z1,
           nleached,total_nitrate);
        
//    printf("leaching[%f+%f]: water(%e,%e,%e,%e,%e)(%e), nleached(%e)<(%e)\n",
//           patchID,m,
//           availableQ, Qout,z1,z2,critialZ,N_decay_rate,
//           nleached,total_nitrate);
	
        
	return(nleached);
} /* end compute_N_leached */

