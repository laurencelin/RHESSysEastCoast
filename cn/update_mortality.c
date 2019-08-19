
/*--------------------------------------------------------------*/
/*                                                              */ 
/*		update_mortality									*/
/*                                                              */
/*  NAME                                                        */
/*		update_mortality									*/
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/* 	void update_mortality( 
/*                      struct epconst_struct,			*/
/*                      struct phenology_struct *,		*/
/*                      struct cstate_struct *,			*/
/*                      struct cdayflux_struct *,		*/
/*                      struct cdayflux_patch_struct *,		*/
/*                      struct nstate_struct *,			*/
/*                      struct ndayflux_struct *,		*/
/*                      struct ndayflux_patch_struct *,		*/
/*                      struct litter_c_object *,		*/
/*                      struct litter_n_object *,		*/
/*                      double);				*/
/*								*/
/*  OPTIONS                                                     */
/*                                                              */
/*  DESCRIPTION                                                 */
/*                                                              */
/*                                                              */
/*	calculated daily mortality losses and updates 		*/
/*	carbon and nitrogen pools				*/
/*                                                              */
/*  PROGRAMMER NOTES                                            */
/*                                                              */
/*	from P.Thornton (1997) version of 1d_bgc		*/
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"
#include "phys_constants.h"
void update_mortality(
					  struct epconst_struct	epc	,
					  struct cstate_struct *cs,
					  struct cdayflux_struct	*cdf,
					  struct cdayflux_patch_struct *cdf_patch,
					  struct nstate_struct *ns,
					  struct ndayflux_struct	*ndf,
					  struct ndayflux_patch_struct *ndf_patch,
					  struct litter_c_object *cs_litr,
					  struct litter_n_object *ns_litr,
                      struct date current_date,
                      struct patch_object *patch,
                      struct canopy_strata_object *stratum,
					  int thintyp,
					  struct mortality_struct mort,
                      int BGC_flag,
                      int ID)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	//long    yearday( struct date);
    
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	double m_cpool, m_npool;
	double m_leafc_to_litr1c, m_leafc_to_litr2c;
	double m_leafc_to_litr3c, m_leafc_to_litr4c;
	double m_deadleafc_to_litr1c, m_deadleafc_to_litr2c;
	double m_deadleafc_to_litr3c, m_deadleafc_to_litr4c;
	double m_frootc_to_litr1c, m_frootc_to_litr2c;
	double m_frootc_to_litr3c, m_frootc_to_litr4c;
	double m_leafc_store_to_litr1c;
	double m_frootc_store_to_litr1c;
	double m_livestemc_store_to_litr1c;
	double m_deadstemc_store_to_litr1c;
	double m_livecrootc_store_to_litr1c;
	double m_deadcrootc_store_to_litr1c;
	double m_leafc_transfer_to_litr1c;
	double m_frootc_transfer_to_litr1c;
	double m_livestemc_transfer_to_litr1c;
	double m_deadstemc_transfer_to_litr1c;
	double m_livecrootc_transfer_to_litr1c;
	double m_deadcrootc_transfer_to_litr1c;
	double m_gresp_store_to_litr1c;
	double m_gresp_transfer_to_litr1c;
	double m_leafn_to_litr1n, m_leafn_to_litr2n;
	double m_leafn_to_litr3n, m_leafn_to_litr4n;
	double m_deadleafn_to_litr1n, m_deadleafn_to_litr2n;
	double m_deadleafn_to_litr3n, m_deadleafn_to_litr4n;
	double m_frootn_to_litr1n, m_frootn_to_litr2n;
	double m_frootn_to_litr3n, m_frootn_to_litr4n;
	double m_livestemn_to_litr1n, m_livecrootn_to_litr1n;
	double m_leafn_store_to_litr1n;
	double m_frootn_store_to_litr1n;
	double m_livestemn_store_to_litr1n;
	double m_deadstemn_store_to_litr1n;
	double m_livecrootn_store_to_litr1n;
	double m_deadcrootn_store_to_litr1n;
	double m_leafn_transfer_to_litr1n;
	double m_frootn_transfer_to_litr1n;
	double m_livestemn_transfer_to_litr1n;
	double m_deadstemn_transfer_to_litr1n;
	double m_livecrootn_transfer_to_litr1n;
	double m_deadcrootn_transfer_to_litr1n;
	double m_retransn_to_litr1n;
	double m_livestemc_to_cwdc;
	double m_deadstemc_to_cwdc;
	double m_livecrootc_to_cwdc;
	double m_deadcrootc_to_cwdc;
	double m_livestemn_to_cwdn;
	double m_deadstemn_to_cwdn;
	double m_livecrootn_to_cwdn;
	double m_deadcrootn_to_cwdn;
    double avgCN;
    
    double add2L1c, add2L1n;
    double add2DOC, add2DON;
	/******************************************************************/
	/* Non-fire mortality: these fluxes all enter litter or CWD pools */
	/******************************************************************/
    
//    printf("update mortality daily (before plantCN): %lf, %lf, %lf, %lf, %lf, %lf\n",
//           cs->leafc/ns->leafn, cs->frootc/ns->frootn,
//           cs->live_stemc/ns->live_stemn, cs->dead_stemc/ns->dead_stemn,
//           cs->live_crootc/ns->live_crootn, cs->dead_crootc/ns->dead_crootn);
//    printf("update mortality (before plant cs): %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
//           cs->leafc, cs->dead_leafc, cs->frootc,
//           cs->live_stemc, cs->dead_stemc, cs->live_crootc, cs->dead_crootc);
//    printf("update mortality (before plant ns): %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
//           ns->leafn, ns->dead_leafn, ns->frootn,
//           ns->live_stemn, ns->dead_stemn, ns->live_crootn, ns->dead_crootn);
//
//    printf("update branch mortality (before cs): %lf, %lf, %lf, %lf\n",
//           cs_litr->litr1c, cs_litr->litr2c,
//           cs_litr->litr3c, cs_litr->litr4c);
//    printf("update branch mortality (before ns): %lf, %lf, %lf, %lf\n",
//           ns_litr->litr1n, ns_litr->litr2n,
//           ns_litr->litr3n, ns_litr->litr4n);
    
//    if( patch[0].ID == 209782){
//        printf("starting mortality@(%d, %d): %e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n",
//               patch[0].ID, ID,
//               mort.mort_leafc, cs->leafc, cs->leafc_store, cs->leafc_transfer, cs->leafc/ns->leafn, cs->leafc_store/ns->leafn_store, cs->leafc_transfer/ns->leafn_transfer,
//               mort.mort_frootc, cs->frootc, cs->frootc_store, cs->frootc_transfer, cs->frootc/ns->frootn, cs->frootc_store/ns->frootn_store, cs->frootc_transfer/ns->frootn_transfer,
//               mort.mort_livestemc, cs->live_stemc, cs->livestemc_store, cs->livestemc_transfer, cs->live_stemc/ns->live_stemn, cs->livestemc_store/ns->livestemn_store, cs->livestemc_transfer/ns->livestemn_transfer,
//               mort.mort_deadstemc, cs->dead_stemc, cs->deadstemc_store, cs->deadstemc_transfer, cs->dead_stemc/ns->dead_stemn, cs->deadstemc_store/ns->deadstemn_store, cs->deadstemc_transfer/ns->deadstemn_transfer,
//               mort.mort_livecrootc, cs->live_crootc, cs->livecrootc_store, cs->livecrootc_transfer, cs->live_crootc/ns->live_crootn, cs->livecrootc_store/ns->livecrootn_store, cs->livecrootc_transfer/ns->livecrootn_transfer,
//               mort.mort_deadcrootc, cs->dead_crootc, cs->deadcrootc_store, cs->deadcrootc_transfer, cs->dead_crootc/ns->dead_crootn, cs->deadcrootc_store/ns->deadcrootn_store, cs->deadcrootc_transfer/ns->deadcrootn_transfer
//               );
//    }
    
    
	/* daily carbon fluxes due to mortality */
	/* mortality fluxes out of leaf and fine root pools */
	/* carbon depts in cpool dies with the plant - could result in a carbon balance issue */
    // cpool is not a real part of the plant; it's just a numerical variable reflecting the stage of carbon allocation.
	// m_cpool = mort.mort_cpool * cs->cpool;
    m_cpool = 0.0;
    m_gresp_store_to_litr1c = 0.0; // = mort.mort_cpool * cs->gresp_store; //<-------------------------
    m_gresp_transfer_to_litr1c = 0.0; //= mort.mort_cpool * cs->gresp_transfer; //<-------------------------
    // ns->retransn, on the other hands, it's a surplus storage of N, not just a transit/stage variable.
    // npool is related to "ns->retransn", and they work in pair
    m_npool = 0.0; //mort.mort_cpool * ns->npool; // these are pure N output (decay rate is too high)
    m_retransn_to_litr1n = 0.0; //mort.mort_cpool * ns->retransn; // these are pure N output (decay rate is too high)
	
    
	m_leafc_to_litr1c = mort.mort_leafc * cs->leafc * epc.leaflitr_flab;
	m_leafc_to_litr2c = mort.mort_leafc * cs->leafc * epc.leaflitr_fucel;
	m_leafc_to_litr3c = mort.mort_leafc * cs->leafc * epc.leaflitr_fscel;
	m_leafc_to_litr4c = mort.mort_leafc * cs->leafc * epc.leaflitr_flig;
	m_deadleafc_to_litr1c = mort.mort_deadleafc * cs->dead_leafc * epc.leaflitr_flab;
	m_deadleafc_to_litr2c = mort.mort_deadleafc * cs->dead_leafc * epc.leaflitr_fucel;
	m_deadleafc_to_litr3c = mort.mort_deadleafc * cs->dead_leafc * epc.leaflitr_fscel;
	m_deadleafc_to_litr4c = mort.mort_deadleafc * cs->dead_leafc * epc.leaflitr_flig;	
	m_frootc_to_litr1c = mort.mort_frootc * cs->frootc * epc.frootlitr_flab;
	m_frootc_to_litr2c = mort.mort_frootc * cs->frootc * epc.frootlitr_fucel;
	m_frootc_to_litr3c = mort.mort_frootc * cs->frootc * epc.frootlitr_fscel;
	m_frootc_to_litr4c = mort.mort_frootc * cs->frootc * epc.frootlitr_flig;
	
    /* mortality fluxes out of storage and transfer pools */
	/* Assumes cpool mortality fraction applies to all non-structural stores and transfers */
    // I would rather think these store and transfer would go to DOC instead; (Laurence Nov 1, 2018)
	m_leafc_store_to_litr1c  = mort.mort_cpool * cs->leafc_store;
    m_leafc_transfer_to_litr1c = mort.mort_cpool * cs->leafc_transfer;
	m_frootc_store_to_litr1c  = mort.mort_cpool * cs->frootc_store;
	m_frootc_transfer_to_litr1c = mort.mort_cpool * cs->frootc_transfer;
	
    
    
	/* TREE-specific carbon fluxes */
	if (epc.veg_type==TREE){
		m_livestemc_to_cwdc = mort.mort_livestemc * cs->live_stemc;
		m_deadstemc_to_cwdc = mort.mort_deadstemc * cs->dead_stemc;
		m_livecrootc_to_cwdc = mort.mort_livecrootc * cs->live_crootc;
		m_deadcrootc_to_cwdc = mort.mort_deadcrootc * cs->dead_crootc;
		
        /* Assumes cpool mortality fraction applies to all non-structural stores and transfers */
        // I would rather think these store and transfer would go to DOC instead; (Laurence Nov 1, 2018)
		m_livestemc_store_to_litr1c  = mort.mort_cpool * cs->livestemc_store;
		m_livestemc_transfer_to_litr1c = mort.mort_cpool * cs->livestemc_transfer;
        
        m_deadstemc_store_to_litr1c  = mort.mort_cpool * cs->deadstemc_store;
		m_deadstemc_transfer_to_litr1c = mort.mort_cpool * cs->deadstemc_transfer;
        
        m_livecrootc_store_to_litr1c  = mort.mort_cpool * cs->livecrootc_store;
		m_livecrootc_transfer_to_litr1c = mort.mort_cpool * cs->livecrootc_transfer;
        
        m_deadcrootc_store_to_litr1c  = mort.mort_cpool * cs->deadcrootc_store;
		m_deadcrootc_transfer_to_litr1c = mort.mort_cpool * cs->deadcrootc_transfer;
        // --->> all to litter 1 ?
	}

    
	/* daily nitrogen fluxes due to mortality */
	/* mortality fluxes out of leaf and fine root pools */
    if (epc.leaf_cn > ZERO && BGC_flag==1){
        m_leafn_to_litr1n = m_leafc_to_litr1c / epc.leaf_cn;
        m_leafn_to_litr2n = m_leafc_to_litr2c / epc.leaf_cn;
        m_leafn_to_litr3n = m_leafc_to_litr3c / epc.leaf_cn;
        m_leafn_to_litr4n = m_leafc_to_litr4c / epc.leaf_cn;
    
    } else if (epc.leaf_cn > ZERO) {
		m_leafn_to_litr2n = m_leafc_to_litr2c / CEL_CN;
		m_leafn_to_litr3n = m_leafc_to_litr3c / CEL_CN;
		m_leafn_to_litr4n = m_leafc_to_litr4c / LIG_CN;
		m_leafn_to_litr1n = mort.mort_leafc*ns->leafn - (m_leafn_to_litr2n+m_leafn_to_litr3n+m_leafn_to_litr4n);
		m_leafn_to_litr1n = max(m_leafn_to_litr1n, 0.0);
    } else {
		m_leafn_to_litr1n = 0.0;
		m_leafn_to_litr2n = 0.0;
		m_leafn_to_litr3n = 0.0;
		m_leafn_to_litr4n = 0.0;
    }
    
    if (epc.leaflitr_cn > ZERO && BGC_flag==1){
        m_deadleafn_to_litr1n = m_deadleafc_to_litr1c / epc.leaflitr_cn;
        m_deadleafn_to_litr2n = m_deadleafc_to_litr2c / epc.leaflitr_cn;
        m_deadleafn_to_litr3n = m_deadleafc_to_litr3c / epc.leaflitr_cn;
        m_deadleafn_to_litr4n = m_deadleafc_to_litr4c / epc.leaflitr_cn;
        
    } else if (epc.leaflitr_cn > ZERO) {
		m_deadleafn_to_litr2n = m_deadleafc_to_litr2c / CEL_CN;
		m_deadleafn_to_litr3n = m_deadleafc_to_litr3c / CEL_CN;
		m_deadleafn_to_litr4n = m_deadleafc_to_litr4c / LIG_CN;
		m_deadleafn_to_litr1n = mort.mort_deadleafc * ns->dead_leafn - (m_deadleafn_to_litr2n + m_deadleafn_to_litr3n + m_deadleafn_to_litr4n);
		m_deadleafn_to_litr1n = max(m_deadleafn_to_litr1n, 0.0);
	} else {
		m_deadleafn_to_litr1n = 0.0;
		m_deadleafn_to_litr2n = 0.0;
		m_deadleafn_to_litr3n = 0.0;
		m_deadleafn_to_litr4n = 0.0;
	}
    
    if (epc.froot_cn > ZERO && BGC_flag==1){
        m_frootn_to_litr1n = m_frootc_to_litr1c / epc.froot_cn;
        m_frootn_to_litr2n = m_frootc_to_litr2c / epc.froot_cn;
        m_frootn_to_litr3n = m_frootc_to_litr3c / epc.froot_cn;
        m_frootn_to_litr4n = m_frootc_to_litr4c / epc.froot_cn;
        
    }else if (epc.froot_cn > ZERO) {
		m_frootn_to_litr2n = m_frootc_to_litr2c / CEL_CN;
		m_frootn_to_litr3n = m_frootc_to_litr3c / CEL_CN;
		m_frootn_to_litr4n = m_frootc_to_litr4c / LIG_CN;
		m_frootn_to_litr1n = mort.mort_frootc*ns->frootn - (m_frootn_to_litr2n+m_frootn_to_litr3n+m_frootn_to_litr4n);
		m_frootn_to_litr1n = max(m_frootn_to_litr1n, 0.0);
		}
	else {
		m_frootn_to_litr1n = 0.0;
		m_frootn_to_litr2n = 0.0;
		m_frootn_to_litr3n = 0.0;
		m_frootn_to_litr4n = 0.0;		
		}

	/* mortality fluxes out of storage and transfer pools */
	/* Assumes same mortality fractions as for c pools */
	/* Assumes cpool mortality fraction applies to all non-structural stores and transfers */
    // I would rather think these store and transfer would go to DOC instead; (Laurence Nov 1, 2018)
	m_leafn_store_to_litr1n  = mort.mort_cpool * ns->leafn_store;
	m_frootn_store_to_litr1n  = mort.mort_cpool * ns->frootn_store;
	m_leafn_transfer_to_litr1n = mort.mort_cpool * ns->leafn_transfer;
	m_frootn_transfer_to_litr1n = mort.mort_cpool * ns->frootn_transfer;
    

	/* TREE-specific nitrogen fluxes */
	if (epc.veg_type==TREE){

        m_livestemn_to_cwdn = m_livestemc_to_cwdc / epc.deadwood_cn;
        m_livecrootn_to_cwdn = m_livecrootc_to_cwdc / epc.deadwood_cn;
        /* Assumes same mortality fractions as for c pools */
        if(BGC_flag==0){
            m_livestemn_to_litr1n = (mort.mort_livestemc * ns->live_stemn) - m_livestemn_to_cwdn;
            m_livecrootn_to_litr1n = (mort.mort_livecrootc * ns->live_crootn) - m_livecrootn_to_cwdn;
        }else{
            //BGC_flag is ON
            m_livestemn_to_litr1n = (mort.mort_livestemc * ns->live_stemn) - m_livestemn_to_cwdn;
            m_livecrootn_to_litr1n = (mort.mort_livecrootc * ns->live_crootn) - m_livecrootn_to_cwdn;
        }
		
        m_deadstemn_to_cwdn = mort.mort_deadstemc * ns->dead_stemn;
        m_deadcrootn_to_cwdn = mort.mort_deadcrootc * ns->dead_crootn;
        
        // I would rather think these store and transfer would go to DOC instead; (Laurence Nov 1, 2018)
        m_livestemn_store_to_litr1n  = mort.mort_cpool * ns->livestemn_store; //
        m_livestemn_transfer_to_litr1n = mort.mort_cpool * ns->livestemn_transfer; //
        m_deadstemn_store_to_litr1n  = mort.mort_cpool * ns->deadstemn_store;
        m_deadstemn_transfer_to_litr1n = mort.mort_cpool * ns->deadstemn_transfer;
        
        m_livecrootn_store_to_litr1n  = mort.mort_cpool * ns->livecrootn_store; //
        m_livecrootn_transfer_to_litr1n = mort.mort_cpool * ns->livecrootn_transfer; //
        m_deadcrootn_store_to_litr1n  = mort.mort_cpool * ns->deadcrootn_store;
        m_deadcrootn_transfer_to_litr1n = mort.mort_cpool * ns->deadcrootn_transfer;
	}//TREE
	

    
	// debugging
    add2DOC = m_gresp_store_to_litr1c + m_gresp_transfer_to_litr1c + m_leafc_store_to_litr1c + m_leafc_transfer_to_litr1c + m_frootc_store_to_litr1c + m_frootc_transfer_to_litr1c + m_livestemc_store_to_litr1c + m_livestemc_transfer_to_litr1c + m_livecrootc_store_to_litr1c + m_livecrootc_transfer_to_litr1c;
    
    add2DON = m_npool + m_retransn_to_litr1n + m_leafn_store_to_litr1n + m_leafn_transfer_to_litr1n +m_frootn_store_to_litr1n + m_frootn_transfer_to_litr1n + m_livestemn_store_to_litr1n + m_livestemn_transfer_to_litr1n + m_livecrootn_store_to_litr1n + m_livecrootn_transfer_to_litr1n;
    
    add2L1c = m_leafc_to_litr1c + m_deadleafc_to_litr1c + m_frootc_to_litr1c + m_deadstemc_store_to_litr1c + m_deadstemc_transfer_to_litr1c + m_deadcrootc_store_to_litr1c + m_deadcrootc_transfer_to_litr1c;
    
    add2L1n = m_leafn_to_litr1n + m_deadleafn_to_litr1n + m_frootn_to_litr1n + m_deadstemn_store_to_litr1n + m_deadstemn_transfer_to_litr1n + m_deadcrootn_store_to_litr1n + m_deadcrootn_transfer_to_litr1n;
    
    if(add2DOC<0 || add2DON<0){
        printf("update mortality DOC issue [%d,%d: %d,%d,%d]BGC_flag[%d], (%e{%e,%e,%e,%e}, %e{%e,%e,%e,%e})\n",
               patch[0].ID,stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year, BGC_flag,
               add2DOC, cs->leafc_store, cs->leafc_transfer, cs->frootc_store, cs->frootc_transfer,
               add2DON, ns->leafn_store, ns->leafn_transfer, ns->frootn_store, ns->frootn_transfer);
    }
    
	double cn_l1,cn_l3, cn_l2,cn_l4,cn_s1,cn_s2,cn_s3,cn_s4;
    if ((cs_litr->litr1c > 0.0) && (ns_litr->litr1n > 0.0))    cn_l1 = cs_litr->litr1c/ns_litr->litr1n; else cn_l1 = 0.0;//LIVELAB_CN;
    if ((cs_litr->litr2c > 0.0) && (ns_litr->litr2n > 0.0))    cn_l2 = cs_litr->litr2c/ns_litr->litr2n; else cn_l2 = 0.0; //CEL_CN;
    if ((cs_litr->litr3c > 0.0) && (ns_litr->litr3n > 0.0))    cn_l3 = cs_litr->litr3c/ns_litr->litr3n; else cn_l3 = 0.0; //CEL_CN;
    if ((cs_litr->litr4c > 0.0) && (ns_litr->litr4n > 0.0))    cn_l4 = cs_litr->litr4c/ns_litr->litr4n; else cn_l4 = 0.0; // LIG_CN;
    
    
    
    if( (add2L1c/add2L1n < 18 && add2L1c>0) || (add2DOC/add2DON < 18 && add2DOC>0) ){
        printf("update mortality[%d,%d: %d,%d,%d]BGC_flag[%d], Begin_litterpool_CN[%e,%e,%e,%e], add2L1(%e)[%e,%e,%e,%e,%e,%e,%e], add2DOC(%e)[%e,%e,%e,%e,%e,%e,%e,%e](%e,%e,%e,%e)\n",
               patch[0].ID,stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year, BGC_flag,
               cn_l1,cn_l2,cn_l3,cn_l4,//l1CN
               // add2L1
               add2L1c/add2L1n, // why this becomes -9.134378e+17?
               m_leafc_to_litr1c/m_leafn_to_litr1n,//1 * nan
               m_deadleafc_to_litr1c/m_deadleafn_to_litr1n,//2 45
               m_frootc_to_litr1c/m_frootn_to_litr1n,//3 * 50
               m_deadstemc_store_to_litr1c/m_deadstemn_store_to_litr1n,//4 nan
               m_deadstemc_transfer_to_litr1c/m_deadstemn_transfer_to_litr1n,//5 nan
               m_deadcrootc_store_to_litr1c/m_deadcrootn_store_to_litr1n,//6 nan
               m_deadcrootc_transfer_to_litr1c/m_deadcrootn_transfer_to_litr1n,//7 nan
               // add2DOC
               add2DOC/add2DON,
               m_leafc_store_to_litr1c/m_leafn_store_to_litr1n,// turns negative
               m_leafc_transfer_to_litr1c/m_leafn_transfer_to_litr1n,
               m_frootc_store_to_litr1c/m_frootn_store_to_litr1n,
               m_frootc_transfer_to_litr1c/m_frootn_transfer_to_litr1n,
               m_livestemc_store_to_litr1c/m_livestemn_store_to_litr1n,
               m_livestemc_transfer_to_litr1c/m_livestemn_transfer_to_litr1n,
               m_livecrootc_store_to_litr1c/m_livecrootn_store_to_litr1n,
               m_livecrootc_transfer_to_litr1c/m_livecrootn_transfer_to_litr1n,
               //last 4
               m_gresp_store_to_litr1c,
               m_gresp_transfer_to_litr1c,
               m_npool,//3
               m_retransn_to_litr1n//4
               );
        // --- Fix this bug!?
        if(m_leafc_to_litr1c<=0 || m_leafn_to_litr1n<=0){ m_leafc_to_litr1c=0.0; m_leafn_to_litr1n=0.0;}
        if(m_leafc_transfer_to_litr1c<=0 || m_leafn_transfer_to_litr1n<=0){ m_leafc_transfer_to_litr1c=0.0; m_leafn_transfer_to_litr1n=0.0;}
        if(m_frootc_store_to_litr1c<=0 || m_frootn_store_to_litr1n<=0){ m_frootc_store_to_litr1c=0.0; m_frootn_store_to_litr1n=0.0;}
        if(m_frootc_transfer_to_litr1c<=0 || m_frootn_transfer_to_litr1n<=0){ m_frootc_transfer_to_litr1c=0.0; m_frootn_transfer_to_litr1n=0.0;}
        
        if(m_livestemc_store_to_litr1c<=0 || m_livestemn_store_to_litr1n<=0){ m_livestemc_store_to_litr1c=0.0; m_livestemn_store_to_litr1n=0.0;}
        if(m_livestemc_transfer_to_litr1c<=0 || m_livestemn_transfer_to_litr1n<=0){ m_livestemc_transfer_to_litr1c=0.0; m_livestemn_transfer_to_litr1n=0.0;}
        if(m_livecrootc_store_to_litr1c<=0 || m_livecrootn_store_to_litr1n<=0){ m_livecrootc_store_to_litr1c=0.0; m_livecrootn_store_to_litr1n=0.0;}
        if(m_livecrootc_transfer_to_litr1c<=0 || m_livecrootn_transfer_to_litr1n<=0){ m_livecrootn_transfer_to_litr1n=0.0; m_livestemn_store_to_litr1n=0.0;}
        
    }//debug
    
    
    
    
    
    
    
    
    if(thintyp == 3){
        // thinning ID 3 process
        // currently missing
        
    }else if(thintyp != 2 && BGC_flag==1){
        // regular: thintyp = 1
//        cs->gresp_store -= m_gresp_store_to_litr1c; //set as zero
//            patch[0].soil_cs.DOC    += m_gresp_store_to_litr1c;
//
//        cs->gresp_transfer -= m_gresp_transfer_to_litr1c; // set as zero
//            patch[0].soil_cs.DOC    += m_gresp_transfer_to_litr1c;
//
//        ns->npool -= m_npool; // set as zero
//            patch[0].soil_ns.DON += m_npool;
//
//        ns->retransn  -= m_retransn_to_litr1n; // set as zero
//            patch[0].soil_ns.DON += m_retransn_to_litr1n;
        
        // live leaf
        cs->leafc -= m_leafc_to_litr1c;
        cs->leafc -= m_leafc_to_litr2c;
        cs->leafc -= m_leafc_to_litr3c;
        cs->leafc -= m_leafc_to_litr4c; cs->leafc = cs->leafc<0? 0.0 : cs->leafc;
        ns->leafn -= m_leafn_to_litr1n;
        ns->leafn -= m_leafn_to_litr2n;
        ns->leafn -= m_leafn_to_litr3n;
        ns->leafn -= m_leafn_to_litr4n; ns->leafn = ns->leafn<0? 0.0 : ns->leafn;
            cs_litr->litr1c    += m_leafc_to_litr1c;
            cs_litr->litr2c    += m_leafc_to_litr2c;
            cs_litr->litr3c    += m_leafc_to_litr3c;
            cs_litr->litr4c    += m_leafc_to_litr4c;
            ns_litr->litr1n    += m_leafn_to_litr1n;
            ns_litr->litr2n    += m_leafn_to_litr2n;
            ns_litr->litr3n    += m_leafn_to_litr3n;
            ns_litr->litr4n    += m_leafn_to_litr4n;
        
        cs->leafc_store -= m_leafc_store_to_litr1c; cs->leafc_store = cs->leafc_store<0? 0.0 : cs->leafc_store;
        ns->leafn_store -= m_leafn_store_to_litr1n; ns->leafn_store = ns->leafn_store<0? 0.0 : ns->leafn_store;
            patch[0].soil_cs.DOC += m_leafc_store_to_litr1c;
            patch[0].soil_ns.DON += m_leafn_store_to_litr1n;
        
        cs->leafc_transfer    -= m_leafc_transfer_to_litr1c; cs->leafc_transfer = cs->leafc_transfer<0? 0.0 : cs->leafc_transfer;
        ns->leafn_transfer    -= m_leafn_transfer_to_litr1n; ns->leafn_transfer = ns->leafn_transfer<0? 0.0 : ns->leafn_transfer;
            patch[0].soil_cs.DOC    += m_leafc_transfer_to_litr1c;
            patch[0].soil_ns.DON    += m_leafn_transfer_to_litr1n;
        
        
        //dead leaf
        cs->dead_leafc  -= m_deadleafc_to_litr1c;
        cs->dead_leafc  -= m_deadleafc_to_litr2c;
        cs->dead_leafc  -= m_deadleafc_to_litr3c;
        cs->dead_leafc  -= m_deadleafc_to_litr4c; cs->dead_leafc = cs->dead_leafc<0? 0.0 : cs->dead_leafc;
        ns->dead_leafn  -= m_deadleafn_to_litr1n;
        ns->dead_leafn  -= m_deadleafn_to_litr2n;
        ns->dead_leafn  -= m_deadleafn_to_litr3n;
        ns->dead_leafn  -= m_deadleafn_to_litr4n; ns->dead_leafn = ns->dead_leafn<0? 0.0 : ns->dead_leafn;
            cs_litr->litr1c    += m_deadleafc_to_litr1c;
            cs_litr->litr2c    += m_deadleafc_to_litr2c;
            cs_litr->litr3c    += m_deadleafc_to_litr3c;
            cs_litr->litr4c    += m_deadleafc_to_litr4c;
            ns_litr->litr1n    += m_deadleafn_to_litr1n;
            ns_litr->litr2n    += m_deadleafn_to_litr2n;
            ns_litr->litr3n    += m_deadleafn_to_litr3n;
            ns_litr->litr4n    += m_deadleafn_to_litr4n;
        
        
        // froot
        cs->frootc  -= m_frootc_to_litr1c;
        cs->frootc  -= m_frootc_to_litr2c;
        cs->frootc  -= m_frootc_to_litr3c;
        cs->frootc  -= m_frootc_to_litr4c; cs->frootc = cs->frootc<0? 0.0 : cs->frootc;
        ns->frootn  -= m_frootn_to_litr1n;
        ns->frootn  -= m_frootn_to_litr2n;
        ns->frootn  -= m_frootn_to_litr3n;
        ns->frootn  -= m_frootn_to_litr4n; ns->frootn = ns->frootn<0? 0.0 : ns->frootn;
            cs_litr->litr1c    += m_frootc_to_litr1c;
            cs_litr->litr2c    += m_frootc_to_litr2c;
            cs_litr->litr3c    += m_frootc_to_litr3c;
            cs_litr->litr4c    += m_frootc_to_litr4c;
            ns_litr->litr1n    += m_frootn_to_litr1n;
            ns_litr->litr2n    += m_frootn_to_litr2n;
            ns_litr->litr3n    += m_frootn_to_litr3n;
            ns_litr->litr4n    += m_frootn_to_litr4n;
        
        cs->frootc_store     -= m_frootc_store_to_litr1c; cs->frootc_store = cs->frootc_store<0? 0.0 : cs->frootc_store;
        ns->frootn_store     -= m_frootn_store_to_litr1n; ns->frootn_store = ns->frootn_store<0? 0.0 : ns->frootn_store;
            patch[0].soil_cs.DOC    += m_frootc_store_to_litr1c;
            patch[0].soil_ns.DON    += m_frootn_store_to_litr1n;
        
        cs->frootc_transfer  -= m_frootc_transfer_to_litr1c; cs->frootc_transfer = cs->frootc_transfer<0? 0.0 : cs->frootc_transfer;
        ns->frootn_transfer  -= m_frootn_transfer_to_litr1n; ns->frootn_transfer = ns->frootn_transfer<0? 0.0 : ns->frootn_transfer;
            patch[0].soil_cs.DOC    += m_frootc_transfer_to_litr1c;
            patch[0].soil_ns.DON    += m_frootn_transfer_to_litr1n;
        
        
        if(epc.veg_type == TREE){
            
            // live stem
            cs->live_stemc  -= m_livestemc_to_cwdc; cs->live_stemc = cs->live_stemc<0? 0.0 : cs->live_stemc;
            ns->live_stemn  -= m_livestemn_to_cwdn;
            ns->live_stemn  -= m_livestemn_to_litr1n; ns->live_stemn = ns->live_stemn<0? 0.0 : ns->live_stemn;
                cs->cwdc    += m_livestemc_to_cwdc;
                ns->cwdn    += m_livestemn_to_cwdn;
                ns->cwdN_stored     += m_livestemn_to_litr1n; //<<--------------
            
            cs->livestemc_store   -= m_livestemc_store_to_litr1c; cs->livestemc_store = cs->livestemc_store<0? 0.0 : cs->livestemc_store;
            ns->livestemn_store   -= m_livestemn_store_to_litr1n; ns->livestemn_store = ns->livestemn_store<0? 0.0 : ns->livestemn_store;
                patch[0].soil_cs.DOC    += m_livestemc_store_to_litr1c;
                patch[0].soil_ns.DON    += m_livestemn_store_to_litr1n;
            
            cs->livestemc_transfer  -= m_livestemc_transfer_to_litr1c; cs->livestemc_transfer = cs->livestemc_transfer<0? 0.0 : cs->livestemc_transfer;
            ns->livestemn_transfer  -= m_livestemn_transfer_to_litr1n; ns->livestemn_transfer = ns->livestemn_transfer<0? 0.0 : ns->livestemn_transfer;
                patch[0].soil_cs.DOC    += m_livestemc_transfer_to_litr1c;
                patch[0].soil_ns.DON    += m_livestemn_transfer_to_litr1n;
            
            
            // dead stem
            cs->dead_stemc  -= m_deadstemc_to_cwdc; cs->dead_stemc = cs->dead_stemc<0? 0.0 : cs->dead_stemc;
            ns->dead_stemn  -= m_deadstemn_to_cwdn; ns->dead_stemn = ns->dead_stemn<0? 0.0 : ns->dead_stemn;
                cs->cwdc    += m_deadstemc_to_cwdc;
                ns->cwdn    += m_deadstemn_to_cwdn;
            
            cs->deadstemc_store   -= m_deadstemc_store_to_litr1c; cs->deadstemc_store = cs->deadstemc_store<0? 0.0 : cs->deadstemc_store;
            ns->deadstemn_store   -= m_deadstemn_store_to_litr1n; ns->deadstemn_store = ns->deadstemn_store<0? 0.0 : ns->deadstemn_store;
                cs_litr->litr1c    += m_deadstemc_store_to_litr1c;
                ns_litr->litr1n    += m_deadstemn_store_to_litr1n;
            
            cs->deadstemc_transfer  -= m_deadstemc_transfer_to_litr1c; cs->deadstemc_transfer = cs->deadstemc_transfer<0? 0.0 : cs->deadstemc_transfer;
            ns->deadstemn_transfer  -= m_deadstemn_transfer_to_litr1n; ns->deadstemn_transfer = ns->deadstemn_transfer<0? 0.0 : ns->deadstemn_transfer;
                cs_litr->litr1c    += m_deadstemc_transfer_to_litr1c;
                ns_litr->litr1n    += m_deadstemn_transfer_to_litr1n;
            
            
            // live croot
            cs->live_crootc -= m_livecrootc_to_cwdc; cs->live_crootc = cs->live_crootc<0? 0.0 : cs->live_crootc;
            ns->live_crootn -= m_livecrootn_to_cwdn;
            ns->live_crootn -= m_livecrootn_to_litr1n; ns->live_crootn = ns->live_crootn<0? 0.0 : ns->live_crootn;
                cs->cwdc       += m_livecrootc_to_cwdc;
                ns->cwdn       += m_livecrootn_to_cwdn;
                ns->cwdN_stored     += m_livecrootn_to_litr1n;  //<<--------------
            
            cs->livecrootc_store  -= m_livecrootc_store_to_litr1c; cs->livecrootc_store = cs->livecrootc_store<0? 0.0 : cs->livecrootc_store;
            ns->livecrootn_store  -= m_livecrootn_store_to_litr1n; ns->livecrootn_store = ns->livecrootn_store<0? 0.0 : ns->livecrootn_store;
                patch[0].soil_cs.DOC += m_livecrootc_store_to_litr1c;
                patch[0].soil_ns.DON += m_livecrootn_store_to_litr1n;
            
            cs->livecrootc_transfer -= m_livecrootc_transfer_to_litr1c; cs->livecrootc_transfer = cs->livecrootc_transfer<0? 0.0 : cs->livecrootc_transfer;
            ns->livecrootn_transfer -= m_livecrootn_transfer_to_litr1n; ns->livecrootn_transfer = ns->livecrootn_transfer<0? 0.0 : ns->livecrootn_transfer;
                patch[0].soil_cs.DOC += m_livecrootc_transfer_to_litr1c;
                patch[0].soil_ns.DON += m_livecrootn_transfer_to_litr1n;
            
            
            // dead croot
            cs->dead_crootc -= m_deadcrootc_to_cwdc; cs->dead_crootc = cs->dead_crootc<0? 0.0 : cs->dead_crootc;
            ns->dead_crootn -= m_deadcrootn_to_cwdn; ns->dead_crootn = ns->dead_crootn<0? 0.0 : ns->dead_crootn;
                cs->cwdc    += m_deadcrootc_to_cwdc;
                ns->cwdn    += m_deadcrootn_to_cwdn;
            
            cs->deadcrootc_store  -= m_deadcrootc_store_to_litr1c; cs->deadcrootc_store = cs->deadcrootc_store<0? 0.0 : cs->deadcrootc_store;
            ns->deadcrootn_store  -= m_deadcrootn_store_to_litr1n; ns->deadcrootn_store = ns->deadcrootn_store<0? 0.0 : ns->deadcrootn_store;
                cs_litr->litr1c += m_deadcrootc_store_to_litr1c;
                ns_litr->litr1n += m_deadcrootn_store_to_litr1n;
            
            cs->deadcrootc_transfer -= m_deadcrootc_transfer_to_litr1c; cs->deadcrootc_transfer = cs->deadcrootc_transfer<0? 0.0 : cs->deadcrootc_transfer;
            ns->deadcrootn_transfer -= m_deadcrootn_transfer_to_litr1n; ns->deadcrootn_transfer = ns->deadcrootn_transfer<0? 0.0 : ns->deadcrootn_transfer;
                cs_litr->litr1c += m_deadcrootc_transfer_to_litr1c;
                ns_litr->litr1n += m_deadcrootn_transfer_to_litr1n;
            
        }else{
            m_livestemc_to_cwdc = 0.0;
            m_livestemn_to_cwdn = 0.0;
            m_livestemn_to_litr1n = 0.0;
            m_livestemc_store_to_litr1c = 0.0;
            m_livestemn_store_to_litr1n = 0.0;
            m_livestemc_transfer_to_litr1c = 0.0;
            m_livestemn_transfer_to_litr1n = 0.0;
            
            m_deadstemc_to_cwdc = 0.0;
            m_deadstemn_to_cwdn = 0.0;
            m_deadstemc_store_to_litr1c = 0.0;
            m_deadstemn_store_to_litr1n = 0.0;
            m_deadstemc_transfer_to_litr1c = 0.0;
            m_deadstemn_transfer_to_litr1n = 0.0;
            
            m_livecrootc_to_cwdc = 0.0;
            m_livecrootn_to_cwdn = 0.0;
            m_livecrootn_to_litr1n = 0.0;
            m_livecrootc_store_to_litr1c = 0.0;
            m_livecrootn_store_to_litr1n = 0.0;
            m_livecrootc_transfer_to_litr1c = 0.0;
            m_livecrootn_transfer_to_litr1n = 0.0;
            
            m_deadcrootc_to_cwdc = 0.0;
            m_deadcrootn_to_cwdn = 0.0;
            m_deadcrootc_store_to_litr1c = 0.0;
            m_deadcrootn_store_to_litr1n = 0.0;
            m_deadcrootc_transfer_to_litr1c = 0.0;
            m_deadcrootn_transfer_to_litr1n = 0.0;
        }
        
    }else if(thintyp != 2){
        // regular
        cs->gresp_store -= m_gresp_store_to_litr1c;
            cs_litr->litr1c    += m_gresp_store_to_litr1c;
        
        cs->gresp_transfer -= m_gresp_transfer_to_litr1c;
            cs_litr->litr1c    += m_gresp_transfer_to_litr1c;
        
        ns->npool -= m_npool;
            ns_litr->litr1n += m_npool;
        
        ns->retransn  -= m_retransn_to_litr1n;
            ns_litr->litr1n += m_retransn_to_litr1n;
        
        // live leaf
        cs->leafc -= m_leafc_to_litr1c;
        cs->leafc -= m_leafc_to_litr2c;
        cs->leafc -= m_leafc_to_litr3c;
        cs->leafc -= m_leafc_to_litr4c;
        ns->leafn -= m_leafn_to_litr1n;
        ns->leafn -= m_leafn_to_litr2n;
        ns->leafn -= m_leafn_to_litr3n;
        ns->leafn -= m_leafn_to_litr4n;
            cs_litr->litr1c    += m_leafc_to_litr1c;
            cs_litr->litr2c    += m_leafc_to_litr2c;
            cs_litr->litr3c    += m_leafc_to_litr3c;
            cs_litr->litr4c    += m_leafc_to_litr4c;
            ns_litr->litr1n    += m_leafn_to_litr1n;
            ns_litr->litr2n    += m_leafn_to_litr2n;
            ns_litr->litr3n    += m_leafn_to_litr3n;
            ns_litr->litr4n    += m_leafn_to_litr4n;
        
        cs->leafc_store -= m_leafc_store_to_litr1c;
        ns->leafn_store -= m_leafn_store_to_litr1n;
            cs_litr->litr1c += m_leafc_store_to_litr1c;
            ns_litr->litr1n += m_leafn_store_to_litr1n;
        
        cs->leafc_transfer    -= m_leafc_transfer_to_litr1c;
        ns->leafn_transfer    -= m_leafn_transfer_to_litr1n;
            cs_litr->litr1c    += m_leafc_transfer_to_litr1c;
            ns_litr->litr1n    += m_leafn_transfer_to_litr1n;
        
        
        //dead leaf
        cs->dead_leafc  -= m_deadleafc_to_litr1c;
        cs->dead_leafc  -= m_deadleafc_to_litr2c;
        cs->dead_leafc  -= m_deadleafc_to_litr3c;
        cs->dead_leafc  -= m_deadleafc_to_litr4c;
        ns->dead_leafn  -= m_deadleafn_to_litr1n;
        ns->dead_leafn  -= m_deadleafn_to_litr2n;
        ns->dead_leafn  -= m_deadleafn_to_litr3n;
        ns->dead_leafn  -= m_deadleafn_to_litr4n;
            cs_litr->litr1c    += m_deadleafc_to_litr1c;
            cs_litr->litr2c    += m_deadleafc_to_litr2c;
            cs_litr->litr3c    += m_deadleafc_to_litr3c;
            cs_litr->litr4c    += m_deadleafc_to_litr4c;
            ns_litr->litr1n    += m_deadleafn_to_litr1n;
            ns_litr->litr2n    += m_deadleafn_to_litr2n;
            ns_litr->litr3n    += m_deadleafn_to_litr3n;
            ns_litr->litr4n    += m_deadleafn_to_litr4n;
        
        
        // froot
        cs->frootc  -= m_frootc_to_litr1c;
        cs->frootc  -= m_frootc_to_litr2c;
        cs->frootc  -= m_frootc_to_litr3c;
        cs->frootc  -= m_frootc_to_litr4c;
        ns->frootn  -= m_frootn_to_litr1n;
        ns->frootn  -= m_frootn_to_litr2n;
        ns->frootn  -= m_frootn_to_litr3n;
        ns->frootn  -= m_frootn_to_litr4n;
            cs_litr->litr1c    += m_frootc_to_litr1c;
            cs_litr->litr2c    += m_frootc_to_litr2c;
            cs_litr->litr3c    += m_frootc_to_litr3c;
            cs_litr->litr4c    += m_frootc_to_litr4c;
            ns_litr->litr1n    += m_frootn_to_litr1n;
            ns_litr->litr2n    += m_frootn_to_litr2n;
            ns_litr->litr3n    += m_frootn_to_litr3n;
            ns_litr->litr4n    += m_frootn_to_litr4n;
        
        cs->frootc_store     -= m_frootc_store_to_litr1c;
        ns->frootn_store      -= m_frootn_store_to_litr1n;
            cs_litr->litr1c    += m_frootc_store_to_litr1c;
            ns_litr->litr1n    += m_frootn_store_to_litr1n;
        
        cs->frootc_transfer  -= m_frootc_transfer_to_litr1c;
        ns->frootn_transfer  -= m_frootn_transfer_to_litr1n;
            cs_litr->litr1c    += m_frootc_transfer_to_litr1c;
            ns_litr->litr1n    += m_frootn_transfer_to_litr1n;
        
        
        if(epc.veg_type == TREE){
            
            // live stem
            cs->live_stemc  -= m_livestemc_to_cwdc;
            ns->live_stemn  -= m_livestemn_to_cwdn;
            ns->live_stemn  -= m_livestemn_to_litr1n;
                cs->cwdc    += m_livestemc_to_cwdc;
                ns->cwdn    += m_livestemn_to_cwdn;
                ns_litr->litr1n     += m_livestemn_to_litr1n;
            
            cs->livestemc_store   -= m_livestemc_store_to_litr1c;
            ns->livestemn_store   -= m_livestemn_store_to_litr1n;
                cs_litr->litr1c    += m_livestemc_store_to_litr1c;
                ns_litr->litr1n    += m_livestemn_store_to_litr1n;
            
            cs->livestemc_transfer  -= m_livestemc_transfer_to_litr1c;
            ns->livestemn_transfer  -= m_livestemn_transfer_to_litr1n;
                cs_litr->litr1c    += m_livestemc_transfer_to_litr1c;
                ns_litr->litr1n    += m_livestemn_transfer_to_litr1n;
            

            // dead stem
            cs->dead_stemc  -= m_deadstemc_to_cwdc;
            ns->dead_stemn  -= m_deadstemn_to_cwdn;
                cs->cwdc    += m_deadstemc_to_cwdc;
                ns->cwdn    += m_deadstemn_to_cwdn;
            
            cs->deadstemc_store   -= m_deadstemc_store_to_litr1c;
            ns->deadstemn_store   -= m_deadstemn_store_to_litr1n;
                cs_litr->litr1c    += m_deadstemc_store_to_litr1c;
                ns_litr->litr1n    += m_deadstemn_store_to_litr1n;
            
            cs->deadstemc_transfer  -= m_deadstemc_transfer_to_litr1c;
            ns->deadstemn_transfer  -= m_deadstemn_transfer_to_litr1n;
                cs_litr->litr1c    += m_deadstemc_transfer_to_litr1c;
                ns_litr->litr1n    += m_deadstemn_transfer_to_litr1n;
            

            // live croot
            cs->live_crootc -= m_livecrootc_to_cwdc;
            ns->live_crootn -= m_livecrootn_to_cwdn;
            ns->live_crootn -= m_livecrootn_to_litr1n;
                cs->cwdc       += m_livecrootc_to_cwdc;
                ns->cwdn       += m_livecrootn_to_cwdn;
                ns_litr->litr1n     += m_livecrootn_to_litr1n;
            
            cs->livecrootc_store  -= m_livecrootc_store_to_litr1c;
            ns->livecrootn_store  -= m_livecrootn_store_to_litr1n;
                cs_litr->litr1c += m_livecrootc_store_to_litr1c;
                ns_litr->litr1n += m_livecrootn_store_to_litr1n;
            
            cs->livecrootc_transfer -= m_livecrootc_transfer_to_litr1c;
            ns->livecrootn_transfer -= m_livecrootn_transfer_to_litr1n;
                cs_litr->litr1c += m_livecrootc_transfer_to_litr1c;
                ns_litr->litr1n += m_livecrootn_transfer_to_litr1n;
            
        
            // dead croot
            cs->dead_crootc -= m_deadcrootc_to_cwdc;
            ns->dead_crootn -= m_deadcrootn_to_cwdn;
                cs->cwdc    += m_deadcrootc_to_cwdc;
                ns->cwdn    += m_deadcrootn_to_cwdn;
            
            cs->deadcrootc_store  -= m_deadcrootc_store_to_litr1c;
            ns->deadcrootn_store  -= m_deadcrootn_store_to_litr1n;
                cs_litr->litr1c += m_deadcrootc_store_to_litr1c;
                ns_litr->litr1n += m_deadcrootn_store_to_litr1n;
            
            cs->deadcrootc_transfer -= m_deadcrootc_transfer_to_litr1c;
            ns->deadcrootn_transfer -= m_deadcrootn_transfer_to_litr1n;
                cs_litr->litr1c += m_deadcrootc_transfer_to_litr1c;
                ns_litr->litr1n += m_deadcrootn_transfer_to_litr1n;
        
        }else{
            m_livestemc_to_cwdc = 0.0;
            m_livestemn_to_cwdn = 0.0;
            m_livestemn_to_litr1n = 0.0;
            m_livestemc_store_to_litr1c = 0.0;
            m_livestemn_store_to_litr1n = 0.0;
            m_livestemc_transfer_to_litr1c = 0.0;
            m_livestemn_transfer_to_litr1n = 0.0;
            
            m_deadstemc_to_cwdc = 0.0;
            m_deadstemn_to_cwdn = 0.0;
            m_deadstemc_store_to_litr1c = 0.0;
            m_deadstemn_store_to_litr1n = 0.0;
            m_deadstemc_transfer_to_litr1c = 0.0;
            m_deadstemn_transfer_to_litr1n = 0.0;
            
            m_livecrootc_to_cwdc = 0.0;
            m_livecrootn_to_cwdn = 0.0;
            m_livecrootn_to_litr1n = 0.0;
            m_livecrootc_store_to_litr1c = 0.0;
            m_livecrootn_store_to_litr1n = 0.0;
            m_livecrootc_transfer_to_litr1c = 0.0;
            m_livecrootn_transfer_to_litr1n = 0.0;
            
            m_deadcrootc_to_cwdc = 0.0;
            m_deadcrootn_to_cwdn = 0.0;
            m_deadcrootc_store_to_litr1c = 0.0;
            m_deadcrootn_store_to_litr1n = 0.0;
            m_deadcrootc_transfer_to_litr1c = 0.0;
            m_deadcrootn_transfer_to_litr1n = 0.0;
        }
        
    }else{
        // thinning ID thintyp=2 processes (grazing)
        // currently missing
        
        
    }// thinning type
	
    
    
    
    
    
	return;
}/*end update_mortality*/

