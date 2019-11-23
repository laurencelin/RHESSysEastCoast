/*--------------------------------------------------------------*/
/* 								*/
/*		compute_infiltration				*/
/*								*/
/*	NAME							*/
/*	compute_infiltration - estimates vertical 		*/
/*		drainage into soil		.		*/
/*								*/
/*								*/
/*	SYNOPSIS						*/
/*	compute_infiltration(				*/
/*				int	,			*/	
/*				double	,			*/
/*				double	,			*/
/*				double	,			*/
/*				double	)			*/
/*								*/
/*	returns:						*/
/*	infiltration - (m water) infiltration 			*/
/*								*/
/*	OPTIONS							*/
/*	int verbose_flag 					*/
/*	double	z - (m) depth to the water table		*/
/*	double S - soil moisture storage			*/
/*	double Ksat_0 - (m/day) sat. hydraulic conductivity	*/
/*				at the surface.			*/
/*	double m_z - decay of conductivity with depth		*/
/*	double p_0 - porosity at surface			*/
/*	double p - porosity decay parameter			*/ 
/*	double precip - incoming precip				*/
/*	double duration - duration				*/
/*      double psi_air_entry - air entry pressure		*/
/*								*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	computes infiltration based on Phillip approach  	*/
/*			(Philip, 1957)				*/
/*	sorptivity parameter is estimated			*/
/*	from Ksat_0 and air entry pressure following		*/
/*	Manley (1977)						*/
/*	infiltration is dependent on rain fall rate which	*/
/*	is calculated from rainfall duration; note that		*/
/*	if rainfall duration is not input (see zone_daily_F)	*/
/*	it is estimated as daylength; user should be 		*/
/*	aware of the potential for error in using these		*/
/*	average daily rainfall rates				*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

double	compute_infiltration(int verbose_flag,
							 double z, // water table depth; sat_def_z
							 double no_use, // rtzS or patch.S
							 double Ksat_vertical,
							 double Ksat,
							 double storage,
							 double POR,
							 double storage_capacity, //<<---- decay of p_0
							 double precip,
							 double duration,
							 double psi_air_entry)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	
	//double porosity;
	//double Ksat;
	double Sp;
	double psi_f;
	double theta_1;
	double intensity, tp;
	double infiltration;
    double Satpct;
	/*--------------------------------------------------------------*/
	/* only infiltrate for on unsaturated soil			*/
	/*--------------------------------------------------------------*/
    // 1) intensity = precip/duration; Vs Ksat --> tp
    // 2) tp vs duration --> infiltration
    // if (intensity <= Ksat) then infiltration = precip * (%pervious in patch)
    // if (intensity > Ksat) then lots of calculations -> infiltration * (%pervious in patch)
    // re-program below
    if(storage_capacity>ZERO && storage<storage_capacity){
        //Ksat =  m_z>0? Ksat_0*m_z*(1-exp(-z/m_z))/z : Ksat_0; // averaged vksat
        intensity = precip/duration;
        
        if(intensity <= Ksat) infiltration = precip;
        else{
            Satpct = storage / storage_capacity;
            theta_1 = (1.0-Satpct)*POR; // averaged POR? // integrated S?
            psi_f = 0.76 * psi_air_entry;
            Sp = sqrt(2 * Ksat * psi_f);
            
            tp = Ksat * psi_f * theta_1 / (intensity * (intensity-Ksat)); //<<--- key
            
            infiltration = max(min( Sp*sqrt(duration-tp)+0.5*Ksat*(duration-tp) + tp*intensity,  precip),0.0);
        }//else
    }else{
        infiltration = 0.0;
    }
    
	infiltration = infiltration * Ksat_vertical;
	return(infiltration);
} /*compute_infiltration*/
