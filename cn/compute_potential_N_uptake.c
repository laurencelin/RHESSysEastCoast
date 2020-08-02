/*--------------------------------------------------------------*/
/* 								*/
/*		compute_potential_N_uptake					*/
/*								*/
/*								*/
/*	NAME							*/
/*	compute_potential_N_uptake -				*/
/*		computes potential N uptake from soil		*/
/*		for this strata without mineralization		*/
/*		limitation					*/
/*								*/
/*	SYNOPSIS						*/
/*	int compute_potential_N_uptake(					*/
/*                          struct epconst_struct,              */
/*			    struct epvar_struct *epv,		*/
/*                          struct cstate_struct *,             */
/*                          struct nstate_struct *,             */
/*                          struct cdayflux_struct *)           */
/*								*/
/*								*/
/*	returns int:						*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*								*/
/*              modified from Peter Thornton (1998)             */
/*                      dynamic - 1d-bgc ver4.0                 */
/*--------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"
double compute_potential_N_uptake(
								  struct	epconst_struct epc,
								  struct	epvar_struct *epv,
								  struct cstate_struct *cs,
								  struct nstate_struct *ns,
								  struct cdayflux_struct *cdf,
                                  struct canopy_strata_object     *stratum)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	//double day_gpp;     /* daily gross production */
	//double day_mresp;   /* daily total maintenance respiration */
	double f1;          /* RATIO   new fine root C : new leaf C     */
	double f2;          /* RATIO   new coarse root C : new stem C   */
	double f3;          /* RATIO   new stem C : new leaf C          */
	double f4;          /* RATIO   new live wood C : new wood C     */
	double g1;          /* RATIO   C respired for growth : C grown  */
	double cnl;         /* RATIO   leaf C:N      */
	double cnfr;        /* RATIO   fine root C:N */
	double cnlw;        /* RATIO   live wood C:N */
	double cndw;        /* RATIO   dead wood C:N */
	//double cnmax;       /* RATIO   max of root and leaf C:N      */
	double c_allometry, n_allometry;
	double plant_ndemand, transfer;
    double growthAdjust;
	/*---------------------------------------------------------------
	Assess the carbon availability on the basis of this day's
	gross production and maintenance respiration costs
	----------------------------------------------------------------*/
	cs->availc = cdf->psn_to_cpool - cdf->total_mr;
	/* no allocation when the daily C balance is negative */
	if (cs->availc < 0.0) cs->availc = 0.0;
	/* test for cpool deficit */
	if (cs->cpool < 0.0){
	/*--------------------------------------------------------------
	running a deficit in cpool, so the first priority
	is to let today's available C accumulate in cpool.  The actual
	accumulation in the cpool is resolved in day_carbon_state().
		--------------------------------------------------------------*/
        transfer = min(cs->availc, -cs->cpool);
        cs->availc -= transfer;
        cs->cpool += transfer;
        
	} /* end if negative cpool */
    if(cs->availc <0) printf("compute_potential_N_uptake: (%e,%e,%e)\n",
                             cdf->psn_to_cpool,
                             cdf->total_mr,
                             cs->availc);
    
    growthAdjust = (1.0+(2.0-stratum[0].phen.daily_allocation * stratum[0].defaults[0][0].epc.storage_transfer_prop)*epc.gr_perc);
    
	/* assign local values for the allocation control parameters */
	f1 = epc.alloc_frootc_leafc; //2
	f2 = epc.alloc_crootc_stemc; //0
	f3 = epc.alloc_stemc_leafc; //0
	f4 = epc.alloc_livewoodc_woodc; //0
	g1 = epc.gr_perc;
	cnl = epc.leaf_cn;//e.g., 25.9 (from read in)
	cnfr = epc.froot_cn;//e.g., 63.5 (from read in)
	cnlw = epc.livewood_cn;//e.g., 75.6 (from read in)
    cndw = epc.deadwood_cn;//e.g., = (epc->deadwood_fucel + epc->deadwood_fscel) * CEL_CN + (epc->deadwood_flig) * LIG_CN;
	/*---------------------------------------------------------------
	given the available C, use constant allometric relationships to
	determine how much N is required to meet this potential growth
	demand */
	/* ----------------------------------------------------------------*/
    // problem, need to count for "excess_lai" Jan 11th. 2019
    //double MAX_LAI = stratum[0].local_max_lai; //epc.max_lai;
    
    
	if (epc.veg_type == TREE){
		c_allometry = growthAdjust*(1.0 + f1 + f3 + f3*f2); //<<------------ correct
		n_allometry = (1.0/cnl + f1/cnfr + (f3*f4*(1.0+f2))/cnlw
			+ (f3*(1.0-f4)*(1.0+f2))/cndw); //<<------------ correct
	}
	else{
//        c_allometry = (1.0 + g1 + f1 + f1*g1);
//        n_allometry = (1.0/cnl + f1/cnfr);
        c_allometry = growthAdjust*(1.0 + f1); //<<------------ correct
        n_allometry = (1.0/cnl + f1/cnfr); //<<------------ correct
	}
    plant_ndemand = cs->availc * (n_allometry / c_allometry);// - max(ns->retransn,0.0);

    cdf->fleaf =  1.0/(1.0+f1+f3+f2*f3);
    cdf->froot = cdf->fleaf*f1;
    cdf->fwood = cdf->fleaf*f3*(1.0+f2);
    cdf->fstem = cdf->fleaf*f3; 
    cdf->fcroot = cdf->fstem *f2;
    
    
    //printf("plant uptake (default): %f,%f,%f,%f,%f\n",cdf->psn_to_cpool,cdf->total_mr,cs->cpool,cs->availc,plant_ndemand);
    
	return(plant_ndemand);
} /* 	end compute_potential_N_uptake */
