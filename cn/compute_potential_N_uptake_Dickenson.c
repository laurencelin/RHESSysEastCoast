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
/*	uses Dickenson et al. 1998, J of Climate		*/
/*	where fraction of allocation to leaves 			*/
/*	is a function of LAI					*/
/*								*/
/*	SYNOPSIS						*/
/*	int compute_potential_N_uptake_Dickenson( int		*/
/*                          struct epconst_struct,              */
/*                          struct epv_var_struct *,              */
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
/*--------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"
double compute_potential_N_uptake_Dickenson( 
								  struct	epconst_struct epc,
								  struct epvar_struct	*epv,
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
	double day_gpp;     /* daily gross production */
	double day_mresp;   /* daily total maintenance respiration */
	double cnl;         /* RATIO   leaf C:N      */
	double cnfr;        /* RATIO   fine root C:N */
	double cnlw;        /* RATIO   live wood C:N */
	double cndw;        /* RATIO   dead wood C:N */
	double fwood, fleaf, froot, flive, fdead; 	/* fraction allocate to each component */
	//double total_wood,total_nonwood, B,C;
	double mean_cn, transfer;
	double plant_ndemand;
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
    if(cs->availc <0) printf("compute_potential_N_uptake_Dick: (%e,%e,%e)\n",
                             cdf->psn_to_cpool,
                             cdf->total_mr,
                             cs->availc);
    
	/* assign local values for the allocation control parameters */
    flive = epc.alloc_livewoodc_woodc;
    fdead = 1.0-flive;
	cnl = epc.leaf_cn;
	cnfr = epc.froot_cn;
	cnlw = epc.livewood_cn;
	cndw = epc.deadwood_cn;
	/*---------------------------------------------------------------*/
	/* constant B and C are currently set for forests from Dickenson et al. */	
	/*----------------------------------------------------------------*/
    
    // problem, need to count for "excess_lai" Jan 11th. 2019
    double MAX_LAI = stratum[0].local_max_lai; //epc.max_lai;
    
	fleaf = exp(-1.0*epc.dickenson_pa * epv->proj_lai); // 3.665163 is the turning point
    fleaf = min(fleaf, 1.0) * (epv->proj_lai<MAX_LAI? 1.0 : 0.0);
    
    //total_nonwood = (cs->dead_leafc + cs->leafc + cs->frootc);
	//total_wood = (cs->live_crootc + cs->dead_crootc + cs->live_stemc + cs->dead_stemc);
    // this calculation here only varies the allocation ratio at daily basis within a year. does not actually simulating the aging effect. stem allocation is less as tree gets old.
    

	if (epc.veg_type==TREE) {
		if (2*fleaf < 0.8) {
            //when LAI is greater than the turning point
            froot = (fleaf>0, fleaf, 0.39); //max(fleaf,0.39);
			fwood= 1.0-fleaf-froot; // --> greater wood growth
        }else{
            //when LAI is less than the turning point
			fleaf = min(fleaf, 0.6); //0.6
			froot = 0.5*(1-fleaf); //0.2
			fwood = 0.5*(1-fleaf); //0.2
        }
    }else{
        //never happened
		fwood = 0;
		froot = (1-fleaf);
    }
	
	

	if ((fleaf+froot) > ZERO) {
        if (epc.veg_type == TREE){
            mean_cn = 1.0 / (fleaf/cnl + froot/cnfr + flive*fwood/cnlw + fwood*fdead/cndw); // one is the sum of froot + fwood + fleaf
        }
        else{
           mean_cn = 1.0 / (fleaf/cnl + froot/cnfr);
        }
	}
	else mean_cn = 1.0;



	cdf->fleaf = fleaf; // leaf
	cdf->froot = froot;// fine root
	cdf->fwood = fwood;// stem and croof together

	/* add in nitrogen for plants and for nitrogen deficit in pool */
    plant_ndemand = cs->availc / (1.0+epc.gr_perc) / mean_cn;// - max(ns->retransn,0.0);
	
//    printf("plant uptake (dickenson): %f,%f,%f,%f,%f\n",cdf->psn_to_cpool,cdf->total_mr,cs->cpool,cs->availc,plant_ndemand);
	return(plant_ndemand);
} /* 	end compute_potential_N_uptake */
