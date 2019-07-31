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
                      double availfrac,
                      int BGC_flag)
{
	/*------------------------------------------------------*/
	/*	Local function declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	
	int ok=1;
	double cwdc_loss, kfrag;
	double rate_scalar;
    double ndemand, nscaler, navail, holdingN;
    //double availNO3_stored, cwdCN;
	/*--------------------------------------------------------------*/
	/* calculate the flux from CWD to litter lignin and cellulose   */
	/*						 compartments, due to physical fragmentation */
	/*--------------------------------------------------------------*/
	/*--------------------------------------------------------------*/
	/*	for now use temperature and water scaleris calculated	*/
	/*	the previous day for soil decomposition limitations	*/
	/*--------------------------------------------------------------*/

	rate_scalar = cs_litr->w_scalar * cs_litr->t_scalar;
	rate_scalar = max(rate_scalar, 0.0);
	/* "physical" fragmentation of coarse woody debris */ //----- but need biological process first to soften and N enrich the wood first.
	kfrag = epc->kfrag_base * rate_scalar;
	cwdc_loss = kfrag * cs->cwdc;
    
    if(BGC_flag==0){ //BGC_flag==0
        cdf->cwdc_to_litr2c = cwdc_loss * epc->deadwood_fucel;
        cdf->cwdc_to_litr3c = cwdc_loss * epc->deadwood_fscel;
        cdf->cwdc_to_litr4c = cwdc_loss * epc->deadwood_flig;
        ndf->cwdn_to_litr2n = cdf->cwdc_to_litr2c/WOODY_CEL_CN;
        ndf->cwdn_to_litr3n = cdf->cwdc_to_litr3c/WOODY_CEL_CN;
        ndf->cwdn_to_litr4n = cdf->cwdc_to_litr4c/WOODY_LIG_CN;
    }else{
        //BGC_flag is ON
        // assume fungus working on this dead log and so the decayed (-> litter) parts are N enriched.
        ndemand = (epc->deadwood_fucel + epc->deadwood_fscel)*1.5/WOODY_CEL_CN + epc->deadwood_flig*1.5/WOODY_LIG_CN; // perc decay carbon
        ndemand *= cwdc_loss;
        navail = ns->cwdN_stored + stratum[0].NO3_stored*availfrac + patch[0].litter.NO3_stored + patch[0].surface_NO3 + patch[0].soil_ns.sminn + patch[0].soil_ns.nitrate;
        nscaler = min(1.0, navail>0.0? navail/ndemand : 0.0);

        cwdc_loss *= nscaler;
        cdf->cwdc_to_litr2c = cwdc_loss * epc->deadwood_fucel;
        cdf->cwdc_to_litr3c = cwdc_loss * epc->deadwood_fscel;
        cdf->cwdc_to_litr4c = cwdc_loss * epc->deadwood_flig;
        ndf->cwdn_to_litr2n = cdf->cwdc_to_litr2c/WOODY_CEL_CN*2.5;
        ndf->cwdn_to_litr3n = cdf->cwdc_to_litr3c/WOODY_CEL_CN*2.5; //CN 100; (new setting 180)
        ndf->cwdn_to_litr4n = cdf->cwdc_to_litr4c/WOODY_LIG_CN*2.5; //CN 142.8571; (new setting 240)

//        if(patch[0].ID==239202 && nscaler<1.0) printf("cwdc decay: %lf %lf %e, %e %e %e %e, %e %e\n",
//                                       epc->kfrag_base, nscaler, ndemand,
//                                       ns->cwdN_stored, stratum[0].NO3_stored, patch[0].litter.NO3_stored, patch[0].surface_NO3,
//                                       patch[0].soil_ns.sminn, patch[0].soil_ns.nitrate);
        
        ndemand *= nscaler;
        holdingN = ns->cwdN_stored; ns->cwdN_stored -= max(0.0,min(holdingN, ndemand)); ndemand -= max(0.0,holdingN);
        if(ndemand>0.0){holdingN = stratum[0].NO3_stored*availfrac; stratum[0].NO3_stored -= max(0.0, min(holdingN, ndemand)); ndemand -= max(0.0,holdingN);}
        if(ndemand>0.0){holdingN = patch[0].litter.NO3_stored; patch[0].litter.NO3_stored -= max(0.0, min(holdingN, ndemand)); ndemand -= max(0.0,holdingN);}
        if(ndemand>0.0){holdingN = patch[0].surface_NO3; patch[0].surface_NO3 -= max(0.0, min(holdingN, ndemand)); ndemand -= max(0.0,holdingN);}
        if(ndemand>0.0){holdingN = patch[0].soil_ns.nitrate; patch[0].soil_ns.nitrate -= max(0.0, min(holdingN, ndemand)); ndemand -= max(0.0,holdingN);}///<<------- future upgrade
        if(ndemand>0.0){holdingN = patch[0].soil_ns.sminn; patch[0].soil_ns.sminn -= max(0.0, min(holdingN, ndemand)); ndemand -= max(0.0,holdingN);}
        if(ndemand>ZERO) printf("compute_cwd_decay: not enough N %e at %d\n",ndemand,patch[0].ID);
        
        //patch[0].soil_ns.sminn += kfrag*nscaler * ns->cwdN_stored * cover_fraction;
        //ns->cwdN_stored *= 1-kfrag*nscaler;
    }
    
    //// see how cwdn was generate
    // m_livestemn_to_cwdn = m_livestemc_to_cwdc / epc.deadwood_cn;
    // epc->deadwood_cn = (epc->deadwood_fucel + epc->deadwood_fscel) * CEL_CN + (epc->deadwood_flig) * LIG_CN;
    
	/*--------------------------------------------------------------*/
	/*	update carbon state variables				*/
	/*--------------------------------------------------------------*/
	cs->cwdc -= (cdf->cwdc_to_litr2c + cdf->cwdc_to_litr3c
		+ cdf->cwdc_to_litr4c);
	cs_litr->litr2c += (cdf->cwdc_to_litr2c * cover_fraction);
	cs_litr->litr3c += (cdf->cwdc_to_litr3c * cover_fraction);
	cs_litr->litr4c += (cdf->cwdc_to_litr4c * cover_fraction);
	/*--------------------------------------------------------------*/
	/*	update nitrogen state variables				*/
	/*--------------------------------------------------------------*/
	ns->cwdn -= (ndf->cwdn_to_litr2n + ndf->cwdn_to_litr3n
		+ ndf->cwdn_to_litr4n);
	ns_litr->litr2n += (ndf->cwdn_to_litr2n * cover_fraction);
	ns_litr->litr3n += (ndf->cwdn_to_litr3n * cover_fraction);
	ns_litr->litr4n += (ndf->cwdn_to_litr4n * cover_fraction);
	return(!ok);
} /*end compute_cwd_decay*/
