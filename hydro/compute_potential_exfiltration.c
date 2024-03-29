/*--------------------------------------------------------------*/
/* 								*/
/*			compute_potential_exfiltration		*/
/*								*/
/*								*/
/*	NAME							*/
/*	compute_potential_exfiltration - estimates maximum	*/
/*		poteential exfiltration  due to soil control	*/
/*								*/
/*	SYNOPSIS						*/
/*	double	compute_potential_exfiltration(			*/
/*			int	,				*/
/*			double	,				*/
/*			double	,				*/
/*			double	,				*/
/*			double	,				*/
/*			double	,				*/
/*			double	,				*/
/*			double	,				*/
/*			double	)				*/
/*								*/
/*	returns:						*/
/*	potential_exfiltration - (m water/ day ) maximum pot 	*/
/*		exfiltration due to soil control.		*/
/*								*/
/*	OPTIONS							*/
/*	int verbose_flag 					*/
/*	double	S (DIM) - relative saturation			*/
/*	double	sat_deficit_z (m) - depth to water table	*/
/*	double	Ksat_0 (m/day) - saturate hydraulic 		*/
/*			conductivity at the surface.		*/
/*	double	m_z (m-1) - decay parameter for Ksat_0		*/	
/*	double	psi_air_entry (Pa) -  air entry pressure.	*/
/*	double	pore_size_index (DIM) - Brooks Corey pore size  */
/*			index parameter.			*/
/*	double 	p  (dim) - porosity  decay rate	        	*/
/*	double 	p_0  (m-1) - porosity   at surface 		*/
/*								*/
/*	DESCRIPTION						*/
/*	Estimates the maximum possible rate of exfiltration of	*/
/*	soil water due to soil control only.  This rate should	*/
/*	be limited to the climatic potential evaporation.	*/
/*								*/
/*	Reference: 						*/
/*	Wigmosta,et.al. 1994 A distributed hydrology-vegetation */
/*	model for complex terrain.  WRR v. 30, no. 6.		*/
/*								*/
/*	Which is a somewhat modified result of the development:	*/
/*								*/
/*	Eagleson, P. S. (1978) "Climate, Soil and Vegetation 	*/
/*	3.  A simplified Model of Soil Moisture Movement in the	*/
/*	Liquid Phase.  WRR, Vol 14, No 5.			*/
/*								*/
/*	We use the geometric mean hydraulic conductivity	*/
/*		of the top and bottom of the unsat zone.	*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*	RAF, Jan 21 1997 					*/
/*	New code for OJP.					*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

double	compute_potential_exfiltration(
									   int 	verbose_flag,
									   double	no_use,
									   double 	wilting_point,
									   double	Ksat_average,
									   double	coef,
									   double	storage,
									   double	storage_capacity,
									   double 	porosity_average, // porosity decay
									   double	S_pow)
{
	/*--------------------------------------------------------------*/
	/*	Local function declaration				*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
	//double	porosity_average;
	//double	Ksat_average, wilting_point;
	
	
	/*--------------------------------------------------------------*/
	/*	Estimate mean porosity.					*/
	/*--------------------------------------------------------------*/
	//porosity_average = p * p_0 * (1-exp(-1*sat_deficit_z/p)) / sat_deficit_z;
	
	/*--------------------------------------------------------------*/
	/*	Estimate mean saturated conductivity.			*/
	/*--------------------------------------------------------------*/
//	if (m_z > ZERO)
//		Ksat_average = m_z* Ksat_0 *(1-exp(-1*sat_deficit_z/m_z)) / sat_deficit_z;
//	else
//		Ksat_average = Ksat_0;
    
    //wilting_point = exp(-pore_size_index * log(2.5/psi_air_entry)); // what's that 2.5?
    // should make this "wilting_point" as exfiltration_wilting_point
    
	
	/*--------------------------------------------------------------*/
	/*	Plug everything into the equation for max infiltration  */
	/*--------------------------------------------------------------*/
    if(storage_capacity<=ZERO){
        return(0.0); // water above
    }else{
        // water table is deep or at surface
        double    potential_exfiltration;
        double    satPct;
        
        //if(storage_capacity>ZERO) satPct = min(1.0, storage / storage_capacity);
        //else satPct = 1.0;
        satPct = min(1.0, storage / storage_capacity);
        
        potential_exfiltration = exp( S_pow * log(satPct)) * coef;
        //sqrt( (8 * porosity_average * Ksat_average * psi_air_entry ) / (3 * (1 + 3 * pore_size_index) * (1 + 4 * pore_size_index)));
        
        potential_exfiltration = min( max(0.0,(satPct-wilting_point)*porosity_average), potential_exfiltration);
        return(potential_exfiltration);
    }
} /*potential_exfiltration*/
