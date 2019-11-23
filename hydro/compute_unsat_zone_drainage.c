/*--------------------------------------------------------------*/
/* 								*/
/*		compute_unsat_zone_drainage			*/
/*								*/
/*	NAME							*/
/*	compute_unsat_zone_drainage - estimates vertical 	*/
/*		drainage from the unsat to sat zone.		*/
/*								*/
/*								*/
/*	SYNOPSIS						*/
/*	compute_unsat_zone_drainage(				*/
/*				int	,			*/	
/*				double	,			*/
/*				double	,			*/
/*				double	,			*/
/*				double	)			*/
/*								*/
/*	returns:						*/
/*	unsat_zone_drainage - (m water) drainage from unsat to	*/
/*			sat zone.				*/
/*								*/
/*	OPTIONS							*/
/*	int verbose_flag 					*/
/*	double	m - Ksat decay parameter			*/
/*	double	z - (m) depth to the water table		*/
/*	double ksat_0 - (m/day) sat. hydraulic conductivity	*/
/*				at the surface.			*/
/*	double	potential_drainage - (m water) difference	*/
/*		between current unsat_zone_storage and current	*/
/*		field capacity.					*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	This routine is designed to estimate unsaturated zone	*/
/*	drainage from soils with "TOPMODEL" properties under	*/
/*	the assumptions that:					*/
/*								*/
/*	i) a field capacity for the unsat zone is known.	*/
/*	ii) the drainage will bethe minimum of drainage to	*/ 
/*		field capacity or the amount resulting from	*/
/*		the current Ksat at the water table.		*/
/*								*/
/*	iii) the change in sat_deficit due to drainage does	*/
/*		not significantly affect field capacity.	*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*	To avoid confilcts with modified versions of TOPMODEL	*/
/*	we call a Ksat(z) curve.				*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

double	compute_unsat_zone_drainage(
									int	verbose_flag,
									int	curve,
									double	PSI,//<<--- soil_defaults[0][0].pore_size_index
									double	storage_capacity,
									double	ksat_0,
                                    double  Ksat1,
									double	storage,
									double	resist_drainage,
                                    double  sat_def)
{
	/*--------------------------------------------------------------*/
	/*	Local function declaration									*/
	/*--------------------------------------------------------------*/
//	double	Ksat_z_curve(
//		int,
//		double,
//		double,
//		double);
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	double	coef;
	double	Scurrent,Smin;
	double	Ksat2;
	double	unsat_zone_drainage;
	/*--------------------------------------------------------------*/
	/*	Compute Ksat 		(m water/day)			*/
	/*								*/
	/*	This undersestimates drainage since Ksat will increase	*/
	/*	as drainage progresses.					*/
	/*--------------------------------------------------------------*/
//	Ksat1  = Ksat_z_curve(
//		verbose_flag,
//		m_z, // consistent to other vertical k processes, caprise, infiltrate
//		z,
//		ksat_0);
    
    // problem: ksat_0 is at surface with unit m/day; it does not take long for S to decrease and limit the drainage, right?
    // solution: given this storage, we have S = w/[x]
    //           then S' = w'/w * S, for any w'
    //           for every second, w' = w - [Ksat2 by sec] ==> S'
    //           and it step at Sstop = resist_drainage/w * S;
    // log(S') = log(w'/w * S) = log(w') - log(w) + log(S);
    
    if( storage_capacity<ZERO || resist_drainage>=storage){
        return(0.0);
    }else{
        if (curve == 1){
            coef = (2/PSI+3)+1;//+1 for integration // (2/PSI+3); originally
            //Scurrent = min(1.0, storage / storage_capacity);
            if(storage_capacity>sat_def) Scurrent = min(1.0,(storage + storage_capacity-sat_def) / storage_capacity);
            else Scurrent = min(1.0, storage / storage_capacity);
            Smin = resist_drainage / storage_capacity;
            // Ksat2 = ksat_0 * exp(coef*log(S)); //<<------ Ksat * (S ^ coef)
            // S = ini.S --> S_fc
            // integrate = ksat/(coef+1) * S ^(coef+1) from Smin to S, where Smin = resist_drainage/storage*S
            Ksat2 = ksat_0 * ( exp(log(Scurrent)*coef) - exp(log(Smin)*coef) )/coef;
        }else{
            Scurrent = min(1.0, storage / storage_capacity);
            Ksat2 = ksat_0 * sqrt(Scurrent) * exp(2*log( 1-exp(PSI*log( 1-exp(log(Scurrent)/PSI)))));
            // (1-(1-S^(1/PSI))^PSI) ^ 2; it's integration is not obvious yet. ... working
        }
        unsat_zone_drainage = min(storage-resist_drainage, max(0.0,min(Ksat1,Ksat2)) );
        return(unsat_zone_drainage);
    }//else
	
    

	
    
	/*--------------------------------------------------------------*/
	/*	Compute unsat zone drainage.				*/
	/*--------------------------------------------------------------*/
	//unsat_zone_drainage = max(min( potential_drainage,	Ksat ),0);

	
} /*compute_unsat_zone_drainage*/
