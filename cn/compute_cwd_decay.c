/*--------------------------------------------------------------*/
/* 								*/
/*			compute_cwd_decay.c			*/
/*								*/
/*								*/
/*	NAME							*/
/*	compute_cwd_decay - 					*/
/*		computes physical fragmentation of coarse	*/
/*		woody debris into litter			*/
/*								*/
/*	SYNOPSIS						*/
/*	double	compute_cwd_decay( 				*/
/*					);			*/	
/*								*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*								*/
/*	source from Peter Thornton, 1d_bgc, 1997		*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

#define wood_fungi_Nenrich_to_woodCN 0.5

int	compute_cwd_decay(
					  struct epconst_struct *epc,
					  double cover_fraction,
					  struct cstate_struct *cs,
					  struct nstate_struct *ns,
					  struct litter_c_object *cs_litr,
					  struct litter_n_object *ns_litr,
					  struct cdayflux_patch_struct *cdf,
					  struct ndayflux_patch_struct *ndf,
					  struct ndayflux_struct *ndf_stratum,
                      struct canopy_strata_object *stratum,
                      struct patch_object    *patch,
                      double availfrac, // 0.02?
                      int BGC_flag)
{
	/*------------------------------------------------------*/
	/*	Local function declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/* calculate the flux from CWD to litter lignin and cellulose   */
	/*						 compartments, due to physical fragmentation */
	/*--------------------------------------------------------------*/
	/*--------------------------------------------------------------*/
	/*	for now use temperature and water scaleris calculated	*/
	/*	the previous day for soil decomposition limitations	*/
	/*--------------------------------------------------------------*/
        /* "physical" fragmentation of coarse woody debris */ //----- but need biological process first to soften and N enrich the wood first.
    //    kfrag = epc->kfrag_base * rate_scalar;
    //    cwdc_loss = kfrag * cs->cwdc;
    
    double cwdc_loss;
    double rate_scalar;
    double ndemand, nscaler, navail, holdingN;
    double CN1, NC1;
    double CN2, NC2;
    double physical_breakdown_c;
    double physical_breakdown_n;
    
    if(BGC_flag>0 && cs->cwdc>0 && ns->cwdn>0){
        
        rate_scalar = max(cs_litr->w_scalar * cs_litr->t_scalar, 0.0);
        physical_breakdown_c = epc->kfrag_base*cs->cwdc;
        physical_breakdown_n = epc->kfrag_base*(ns->cwdN_stored + ns->cwdn);
        ns->cwdN_stored *= 1.0-epc->kfrag_base;
        ns->cwdn *= 1.0-epc->kfrag_base;
        cs->cwdc *= 1.0-epc->kfrag_base;
        
        // dC/CN1 +exN = dC/CN2 && CN2 = #CN1
        // => dC = exN/(1/#-1) * CN1 <--eq. [1]
        // CN1 = 1/( (epc->deadwood_fucel + epc->deadwood_fscel)/WOODY_CEL_CN + (epc->deadwood_flig)/WOODY_LIG_CN )
        // CN2 = #/( (epc->deadwood_fucel + epc->deadwood_fscel)/WOODY_CEL_CN + (epc->deadwood_flig)/WOODY_LIG_CN )
        nscaler = 1.0/wood_fungi_Nenrich_to_woodCN - 1.0;
        
        navail = ns->cwdN_stored + patch[0].litter.NO3_stored + patch[0].surface_NO3 + patch[0].soil_ns.sminn + patch[0].soil_ns.nitrate; // = exN
        NC1 = (epc->deadwood_fucel + epc->deadwood_fscel)/WOODY_CEL_CN + (epc->deadwood_flig)/WOODY_LIG_CN ;
        CN1 = 1.0 / NC1;
        CN2 = wood_fungi_Nenrich_to_woodCN * CN1;
        NC2 = 1.0 / CN2;
        cwdc_loss = min(navail*CN1/nscaler, epc->kfrag_base * rate_scalar * cs->cwdc);
        ndemand = cwdc_loss * NC1 * nscaler; // cwdc_loss/CN2 - cwdc_loss/CN1
        
        cdf->cwdc_to_litr2c = (cwdc_loss+physical_breakdown_c) * epc->deadwood_fucel;
        cdf->cwdc_to_litr3c = (cwdc_loss+physical_breakdown_c) * epc->deadwood_fscel;
        cdf->cwdc_to_litr4c = (cwdc_loss+physical_breakdown_c) * epc->deadwood_flig;
        ndf->cwdn_to_litr2n = (cwdc_loss*NC2+physical_breakdown_n) * epc->deadwood_fucel;
        ndf->cwdn_to_litr3n = (cwdc_loss*NC2+physical_breakdown_n) * epc->deadwood_fscel;
        ndf->cwdn_to_litr4n = (cwdc_loss*NC2+physical_breakdown_n) * epc->deadwood_flig;
        
        ns->cwdn *= cwdc_loss/cs->cwdc;
        ns_litr->litr2n += (ndf->cwdn_to_litr2n * cover_fraction);
        ns_litr->litr3n += (ndf->cwdn_to_litr3n * cover_fraction);
        ns_litr->litr4n += (ndf->cwdn_to_litr4n * cover_fraction);

        cs->cwdc -= cwdc_loss;
        cs_litr->litr2c += (cdf->cwdc_to_litr2c * cover_fraction);
        cs_litr->litr3c += (cdf->cwdc_to_litr3c * cover_fraction);
        cs_litr->litr4c += (cdf->cwdc_to_litr4c * cover_fraction);
        
        
        holdingN = ns->cwdN_stored; ns->cwdN_stored -= max(0.0,min(holdingN, ndemand)); ndemand -= max(0.0,holdingN);
        if(ndemand>0.0){holdingN = stratum[0].NO3_stored*availfrac; stratum[0].NO3_stored -= max(0.0, min(holdingN, ndemand)); ndemand -= max(0.0,holdingN);}
        if(ndemand>0.0){holdingN = patch[0].litter.NO3_stored; patch[0].litter.NO3_stored -= max(0.0, min(holdingN, ndemand)); ndemand -= max(0.0,holdingN);}
        if(ndemand>0.0){holdingN = patch[0].surface_NO3; patch[0].surface_NO3 -= max(0.0, min(holdingN, ndemand)); ndemand -= max(0.0,holdingN);}
        if(ndemand>0.0){holdingN = patch[0].soil_ns.nitrate; patch[0].soil_ns.nitrate -= max(0.0, min(holdingN, ndemand)); ndemand -= max(0.0,holdingN);}///<<------- future upgrade
        if(ndemand>0.0){holdingN = patch[0].soil_ns.sminn; patch[0].soil_ns.sminn -= max(0.0, min(holdingN, ndemand)); ndemand -= max(0.0,holdingN);}
        if(ndemand>ZERO) printf("compute_cwd_decay: not enough N %e at %d\n",ndemand,patch[0].ID);

        
    }//if
    

    
	return(0);
} /*end compute_cwd_decay*/
