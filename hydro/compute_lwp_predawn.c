/*--------------------------------------------------------------*/
/* 								*/
/*					compute_lwp_predawn	*/
/*								*/
/* 	compute_lwp_predawn - computes pre-dawn leaf water potential*/
/*							*/
/*	NAME						*/
/*							*/
/*	SYNOPSIS					*/
/*	double	compute_lwp_predawn(			*/
/*			int	,			*/
/*			int	,			*/
/*			double	,			*/
/*			double	,			*/
/*			double	,			*/
/*			double	,			*/
/*			double	,			*/
/*			double	,			*/
/*			double	,			*/
/*			double	,			*/
/*			double	,			*/
/*			double	,			*/
/*			double	);			*/
/*							*/
/*	OPTIONS						*/
/*	int	verbose_flag,				*/
/*	int	curve,					*/
/*	double	Tsoil - (deg C)				*/
/*	double	LWP_min_spring - (Mpa)			*/
/*	double	LWP_stom_closure - (Mpa)		*/
/*  	double  psi_air_entry - curve parameter 1                  */
/*      double  pore_size_index - curve parameter 2                  */
/*	double  p3					*/
/*	double  p4					*/
/*      double  p_0 - porosity		                */
/*      double  p porosity decay                  	*/
/*	double	sat_deficit- (m)			*/
/*	double	unsat_storage - (m water)		*/
/*	double	rooting_depth				*/
/*	double	field_capacity - (m water)		*/
/*							*/
/*							*/
/*	DESCRIPTION					*/
/*							*/
/* assume lwp is the same as soil water potential 	*/
/* in the rooting zone					*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

double	compute_lwp_predawn(
							int	verbose_flag,
							int	curve,
							double	Tsoil,
							double	LWP_min_spring, // epc.psi_open
							double	LWP_stom_closure, // epc.psi_close
							double	psi_air_entry,
							double 	pore_size_index,
							double  p3,
							double  p4,
							double  storage,
							double  storage_capacity,
							double	sat_def)
{
	/*--------------------------------------------------------------*/
	/*	Local function declaration				*/
	/*--------------------------------------------------------------*/

	double  compute_soil_water_potential(
		int,
		int,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double);

	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
	double	LWP_predawn;

	LWP_predawn = compute_soil_water_potential(verbose_flag,
		curve,
		Tsoil,
		LWP_min_spring,
		LWP_stom_closure,
		psi_air_entry,
		pore_size_index,
		p3,
		p4,
		storage,
		storage_capacity,
		sat_def);

	return(LWP_predawn);
} /*end compute_lwp_predawn)*/
