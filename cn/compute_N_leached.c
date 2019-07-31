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
            double N_decay_rate,
            double activedepthz_,
            double N_absorption_rate,
            int signal,
            struct patch_object *patch)
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
	
		
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
        
        double nleached = 0.0;
        double z1 = 0.0;
        double Qout_frac = 0.0;
        double availableN = 0.0;
        double N0 = 0.0;
        double absorptionConst = 0.0;
        double critialZ = 0.0;
        double p_0 = patch[0].soil_defaults[0][0].porosity_0;
        double p_decayRate = 1.0 / patch[0].soil_defaults[0][0].porosity_decay;
        double constantHold1, constantHold2;
			//double N_decay_rate = (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NO3decayRate : patch[0].soil_defaults[0][0].N_decay_rate);
			//double activedepthz = (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z));
        
        z1 = (patch[0].sat_deficit_z>0? patch[0].sat_deficit_z : 0.0);
        z1 = z1>patch[0].constraintWaterTableTopDepth? z1 : patch[0].constraintWaterTableTopDepth; // do this correction for basement
        double activedepthz = max(z1 + 0.33, patch[0].rootzone.depth);
	
        if((patch[0].drainage_type == ROAD) && signal<0){
            z1 = patch[0].road_cut_depth<activedepthz? patch[0].road_cut_depth : (activedepthz*0.9); // activedepthz must be > road_cut_depth
        }
    

	/*------------------------------------------------------*/
	/* nitrate export only occurs when Qout > 0.0		*/ 
	/*------------------------------------------------------*/
    if (Qout > 0 && total_nitrate>0 && activedepthz>0 && z1<activedepthz) {
        constantHold1 = exp(-p_decayRate * activedepthz);
        constantHold2 = exp(-p_decayRate * z1);
        
        Qout_frac = Qout;
        Qout_frac /= patch[0].soil_defaults[0][0].porosity_0 * patch[0].soil_defaults[0][0].porosity_decay * (constantHold2 - constantHold1);
        Qout_frac = min(Qout_frac,1.0);
        //------------------------------------------------------< sat_solute detailed > ----------------------------------------------------------------//
        if(signal%3==2){
            //---------------------------------------------------------- N decay distribution
            // assume "total_nitrate" is "sat_solute" integrated from z1 to activedepthz.
            N0 = p_decayRate * total_nitrate / (constantHold2 - constantHold1); //p_decayRate * total_nitrate / (1.0 - constantHold1);
            availableN = total_nitrate;//*(constantHold2 - constantHold1)/(1.0 - constantHold1);
            absorptionConst = PARTICLE_DENSITY * 1000.0 * N_absorption_rate;
            if(absorptionConst>0){
               critialZ = (absorptionConst*p_0 - N0)/absorptionConst/p_0/p_0;
                if(critialZ<=0){
                    // any z is possible
                    nleached = patch[0].soil_defaults[0][0].porosity_decay * ((N0-absorptionConst*p_0)*(constantHold2-constantHold1) + absorptionConst*p_0*p_0*0.5*(constantHold2*constantHold2-constantHold1*constantHold1));
                    nleached = Qout_frac * max(min(nleached,availableN), 0.0);
                }else if(critialZ<1.0){
                    critialZ = -log(critialZ)*patch[0].soil_defaults[0][0].porosity_decay;
                    critialZ = min(critialZ,activedepthz);
                    if(critialZ > z1){ //sat_deficit_z
                        constantHold1 = exp(-p_decayRate * critialZ);
                        nleached = patch[0].soil_defaults[0][0].porosity_decay * ((N0-absorptionConst*p_0)*(constantHold2-constantHold1) + absorptionConst*p_0*p_0*0.5*(constantHold2*constantHold2-constantHold1*constantHold1));
                        nleached = Qout_frac *  max(min(nleached,availableN), 0.0);
                    }else{
                        //nothing comes out; absorptionConst is too great
                        nleached = 0.0;
                    }
                }else{ nleached = 0.0; }
            }else{
                // absorptionConst = 0
                nleached = Qout_frac * availableN;
            }
        }
        
        //------------------------------------------------------< sat_solute simplified > ----------------------------------------------------------------//
        else if(signal%3==1){
            constantHold1 = exp(-p_decayRate * activedepthz);
            constantHold2 = exp(-p_decayRate * z1);
            //---------------------------------------------------------- N decay distribution
            // assume "total_nitrate" is "sat_solute" integrated from z1 to activedepthz.
            N0 = p_decayRate * total_nitrate / (constantHold2 - constantHold1);
            availableN = total_nitrate; //*(constantHold2 - constantHold1)/(1.0 - constantHold1);
            absorptionConst = p_0*(1.0-p_0)*PARTICLE_DENSITY * 1000.0 * N_absorption_rate;
            availableN -= absorptionConst;
            if( availableN >0){
                nleached = availableN*(activedepthz-z1);
            }else{
                nleached = 0.0;
            }
        }
        
        //-------------------------------------------------------------------------------------------------------------------------------------
        else if(signal%3==0 && N_decay_rate>0){
            constantHold1 = exp(-N_decay_rate * activedepthz);
            constantHold2 = exp(-N_decay_rate * z1);
            //---------------------------------------------------------- N decay distribution
            N0 = N_decay_rate * total_nitrate / (1.0 - constantHold1);
            availableN = total_nitrate*(constantHold2 - constantHold1)/(1.0 - constantHold1);
            absorptionConst = p_0*(1.0-p_0)*PARTICLE_DENSITY * 1000.0 * N_absorption_rate;
            critialZ =  absorptionConst / N0; //assume p (porosity decay - meter) is large 4000 m (default)
           
            if(critialZ>0 && critialZ < 1.0 && absorptionConst>0 && N0 > 0){
                critialZ = -log(critialZ)/N_decay_rate;
                if(critialZ > z1){ //sat_deficit_z
                    // integrating from sat_deficit_z(=z1) to critialZ; approaximation when (porosity decay - meter) is large 4000 m (default)
                    nleached = Qout_frac * (absorptionConst*(z1-critialZ) + N0/N_decay_rate*(constantHold2-exp(-N_decay_rate*critialZ)) );
                    nleached = max(min(nleached,availableN), 0.0);
                    //if(nleached!=nleached) nleached = 0.0;
                }else{
                    //nothing comes out; absorptionConst is too great
                    nleached = 0.0;
                }
            }else if(N0 > 0 && !(absorptionConst>0) ){
                // NO absorption, e.g., NO3
                // critialZ=inf; absorptionConst=0
                // integrating from sat_deficit_z(=z1) to critialZ(=inf)
                constantHold1 = exp(-N_decay_rate * activedepthz);
                constantHold2 = exp(-N_decay_rate * z1);
                nleached = Qout_frac * total_nitrate * (constantHold2-constantHold1) / (1.0 - constantHold1);
                
                nleached = max(min(nleached,availableN), 0.0);
                
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
        //-------------------------------------------------------------------------------------------------------------------------------------
        
        
        
        
    }else{nleached = 0.0;}//end of Qout > ZERO
        
	
	/*------------------------------------------------------*/
	/* there may be enough flow to leach out more than 	*/
	/*	availabe nitrate, so limit export by available	*/
	/*------------------------------------------------------*/

    //if (nleached > total_nitrate){nleached = total_nitrate;}
        
        
    if(nleached<0 || nleached!=nleached) printf("leaching[%d]: ->(%e,%e,%e), N0(%e)=(%e), absorptionConst(%e)=(%e), critialZ_ini(%e), critialZ(%e)>(%e), nleached(%e)<(%e)\n",
           patch[0].ID,
           Qout,Qout_frac,Qout/Qout_frac,
           N0, N_decay_rate,
           N_absorption_rate, absorptionConst,
           absorptionConst/N0,
           critialZ, z1,
           nleached,total_nitrate);
        
        
	return(nleached);
} /* end compute_N_leached */

