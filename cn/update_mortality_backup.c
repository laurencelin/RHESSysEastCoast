
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
                      struct patch_object *patch,
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
    
    double leaffrootLiter1c, leaffrootLiter1n;
    double stemcrootLiter1c, stemcrootLiter1n;
    double add_to_cwdc, add_to_cwdn;
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
	/* carbon depts in cpool have to die with the plant - could result in a carbon balance issue */
	m_cpool = mort.mort_cpool * cs->cpool;
	m_npool = mort.mort_cpool * ns->npool;
	
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
	m_leafc_store_to_litr1c  = mort.mort_cpool * cs->leafc_store;
	m_frootc_store_to_litr1c  = mort.mort_cpool * cs->frootc_store;
	m_leafc_transfer_to_litr1c = mort.mort_cpool * cs->leafc_transfer;
	m_frootc_transfer_to_litr1c = mort.mort_cpool * cs->frootc_transfer;
	m_gresp_store_to_litr1c = mort.mort_cpool * cs->gresp_store; //<-------------------------
	m_gresp_transfer_to_litr1c = mort.mort_cpool * cs->gresp_transfer; //<-------------------------
    // --->> all to litter 1 ?
    
	/* TREE-specific carbon fluxes */
	if (epc.veg_type==TREE){
		m_livestemc_to_cwdc = mort.mort_livestemc * cs->live_stemc;
		m_deadstemc_to_cwdc = mort.mort_deadstemc * cs->dead_stemc;
		m_livecrootc_to_cwdc = mort.mort_livecrootc * cs->live_crootc;
		m_deadcrootc_to_cwdc = mort.mort_deadcrootc * cs->dead_crootc;
		
        /* Assumes cpool mortality fraction applies to all non-structural stores and transfers */ //<<--- why all goes to litter1?
        add_to_cwdc;
		m_livestemc_store_to_litr1c  = mort.mort_cpool * cs->livestemc_store;
		m_deadstemc_store_to_litr1c  = mort.mort_cpool * cs->deadstemc_store;
		m_livecrootc_store_to_litr1c  = mort.mort_cpool * cs->livecrootc_store;
		m_deadcrootc_store_to_litr1c  = mort.mort_cpool * cs->deadcrootc_store;
		m_livestemc_transfer_to_litr1c = mort.mort_cpool * cs->livestemc_transfer;
		m_deadstemc_transfer_to_litr1c = mort.mort_cpool * cs->deadstemc_transfer;
		m_livecrootc_transfer_to_litr1c = mort.mort_cpool * cs->livecrootc_transfer;
		m_deadcrootc_transfer_to_litr1c = mort.mort_cpool * cs->deadcrootc_transfer;
        // --->> all to litter 1 ?
	}
    leaffrootLiter1c = 0.0;
    leaffrootLiter1c += m_leafc_to_litr1c + m_deadleafc_to_litr1c + m_frootc_to_litr1c;
    leaffrootLiter1c += m_leafc_store_to_litr1c + m_frootc_store_to_litr1c + m_leafc_transfer_to_litr1c + m_frootc_transfer_to_litr1c;
    leaffrootLiter1c += m_gresp_store_to_litr1c + m_gresp_transfer_to_litr1c;
    
    stemcrootLiter1c = 0.0;
    stemcrootLiter1c += m_livestemc_store_to_litr1c + m_deadstemc_store_to_litr1c + m_livecrootc_store_to_litr1c + m_deadcrootc_store_to_litr1c;
    stemcrootLiter1c += m_livestemc_transfer_to_litr1c + m_deadstemc_transfer_to_litr1c;
    stemcrootLiter1c += m_livecrootc_transfer_to_litr1c + m_deadcrootc_transfer_to_litr1c;
    
    
    
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
	/* Assumes cpool mortality fraction applies to all non-structural stores and transfers */ //<<--- why all goes to litter1?
	m_leafn_store_to_litr1n  = mort.mort_cpool * ns->leafn_store;
	m_frootn_store_to_litr1n  = mort.mort_cpool * ns->frootn_store;
	m_leafn_transfer_to_litr1n = mort.mort_cpool * ns->leafn_transfer;
	m_frootn_transfer_to_litr1n = mort.mort_cpool * ns->frootn_transfer;
    m_retransn_to_litr1n = 0.0; //mort.mort_cpool * ns->retransn; // ************************************

	/* TREE-specific nitrogen fluxes */
	if (epc.veg_type==TREE){

        m_livestemn_to_cwdn = m_livestemc_to_cwdc / epc.deadwood_cn;
        m_livecrootn_to_cwdn = m_livecrootc_to_cwdc / epc.deadwood_cn;
        /* Assumes same mortality fractions as for c pools */
        if(BGC_flag==0){
            m_livestemn_to_litr1n = (mort.mort_livestemc * ns->live_stemn) - m_livestemn_to_cwdn;
            m_livecrootn_to_litr1n = (mort.mort_livecrootc * ns->live_crootn) - m_livecrootn_to_cwdn;
            
            m_livestemn_store_to_litr1n  = mort.mort_cpool * ns->livestemn_store; //
            m_livestemn_transfer_to_litr1n = mort.mort_cpool * ns->livestemn_transfer; //
            
            m_livecrootn_store_to_litr1n  = mort.mort_cpool * ns->livecrootn_store; //
            m_livecrootn_transfer_to_litr1n = mort.mort_cpool * ns->livecrootn_transfer; //
            
            m_deadstemn_to_cwdn = mort.mort_deadstemc * ns->dead_stemn;
            m_deadstemn_store_to_litr1n  = mort.mort_cpool * ns->deadstemn_store;
            m_deadstemn_transfer_to_litr1n = mort.mort_cpool * ns->deadstemn_transfer;
            
            m_deadcrootn_to_cwdn = mort.mort_deadcrootc * ns->dead_crootn;
            m_deadcrootn_store_to_litr1n  = mort.mort_cpool * ns->deadcrootn_store;
            m_deadcrootn_transfer_to_litr1n = mort.mort_cpool * ns->deadcrootn_transfer;
            
            add_to_cwdn = 0.0;
        }else{
            //BGC_flag is ON
            m_livestemn_to_litr1n = 0.0;
            m_livestemn_store_to_litr1n = 0.0;
            m_livestemn_transfer_to_litr1n = 0.0;
            m_deadstemn_store_to_litr1n = 0.0;
            m_deadstemn_transfer_to_litr1n = 0.0;
            
            m_livecrootn_to_litr1n = 0.0;
            m_livecrootn_store_to_litr1n = 0.0;
            m_livecrootn_transfer_to_litr1n = 0.0;
            m_deadcrootn_store_to_litr1n = 0.0;
            m_deadcrootn_transfer_to_litr1n = 0.0;
            
            m_deadstemn_to_cwdn = mort.mort_deadstemc * ns->dead_stemn;
            m_deadcrootn_to_cwdn = mort.mort_deadcrootc * ns->dead_crootn;
            
            //patch[0].soil_ns.sminn += (mort.mort_livestemc * ns->live_stemn) - m_livestemn_to_cwdn + (mort.mort_livecrootc * ns->live_crootn) - m_livecrootn_to_cwdn; // missing cover fraction!
            ns->cwdN_stored += (mort.mort_livestemc * ns->live_stemn) - m_livestemn_to_cwdn;
            ns->cwdN_stored += mort.mort_cpool * ns->livestemn_store; // cwd does not have store or transfer model compartments
            ns->cwdN_stored += mort.mort_cpool * ns->livestemn_transfer;
            
            ns->cwdN_stored += (mort.mort_livecrootc * ns->live_crootn) - m_livecrootn_to_cwdn;
            ns->cwdN_stored += mort.mort_cpool * ns->livecrootn_store;
            ns->cwdN_stored += mort.mort_cpool * ns->livecrootn_transfer;
            
            ns->cwdN_stored += mort.mort_cpool * ns->deadstemn_store;
            ns->cwdN_stored += mort.mort_cpool * ns->deadcrootn_store;
            ns->cwdN_stored += mort.mort_cpool * ns->deadstemn_transfer;
            ns->cwdN_stored += mort.mort_cpool * ns->deadcrootn_transfer;
            //if(patch[0].ID==239202) printf("update_mortality: %e\n",(mort.mort_livestemc * ns->live_stemn) - m_livestemn_to_cwdn + (mort.mort_livecrootc * ns->live_crootn) - m_livecrootn_to_cwdn);
        }
		
       
	}//TREE
	
    leaffrootLiter1n = 0.0;
    leaffrootLiter1n += m_leafn_to_litr1n + m_deadleafn_to_litr1n + m_frootn_to_litr1n;
    leaffrootLiter1n += m_leafn_store_to_litr1n + m_frootn_store_to_litr1n + m_leafn_transfer_to_litr1n + m_frootn_transfer_to_litr1n;
    
    stemcrootLiter1n = 0.0;
    stemcrootLiter1n += m_livestemn_to_litr1n + m_livecrootn_to_litr1n;
    stemcrootLiter1n += m_livestemn_store_to_litr1n + m_deadstemn_store_to_litr1n;
    stemcrootLiter1n += m_livecrootn_store_to_litr1n + m_deadcrootn_store_to_litr1n;
    stemcrootLiter1n += m_livestemn_transfer_to_litr1n + m_deadstemn_transfer_to_litr1n;
    stemcrootLiter1n += m_livecrootn_transfer_to_litr1n + m_deadcrootn_transfer_to_litr1n;
    
    
	/* update state variables */
	
	/* ---------------------------------------- */
	/* CARBON mortality state variable update   */
	/* ---------------------------------------- */
	double cn_l1,cn_l3, cn_l2,cn_l4,cn_s1,cn_s2,cn_s3,cn_s4;
    if ((cs_litr->litr1c > 0.0) && (ns_litr->litr1n > 0.0))    cn_l1 = cs_litr->litr1c/ns_litr->litr1n; else cn_l1 = 0.0;//LIVELAB_CN;
    if ((cs_litr->litr2c > 0.0) && (ns_litr->litr2n > 0.0))    cn_l2 = cs_litr->litr2c/ns_litr->litr2n; else cn_l2 = 0.0; //CEL_CN;
    if ((cs_litr->litr3c > 0.0) && (ns_litr->litr3n > 0.0))    cn_l3 = cs_litr->litr3c/ns_litr->litr3n; else cn_l3 = 0.0; //CEL_CN;
    if ((cs_litr->litr4c > 0.0) && (ns_litr->litr4n > 0.0))    cn_l4 = cs_litr->litr4c/ns_litr->litr4n; else cn_l4 = 0.0; // LIG_CN;
    

	/* ABOVEGROUND C POOLS */
	
	/* Only add dead leaf and stem c to litter and cwd pools if thintyp   */
	/* is not 2. If thintyp is 2, harvest aboveground c.   */
	if (thintyp != 2) {
//        cs_litr->litr1c    += m_cpool>0? m_cpool : 0.0;
		/*    Leaf mortality */
		cs_litr->litr1c    += m_leafc_to_litr1c;
		cs_litr->litr2c    += m_leafc_to_litr2c;
		cs_litr->litr3c    += m_leafc_to_litr3c;
		cs_litr->litr4c    += m_leafc_to_litr4c;
		cs_litr->litr1c    += m_deadleafc_to_litr1c;
		cs_litr->litr2c    += m_deadleafc_to_litr2c;
		cs_litr->litr3c    += m_deadleafc_to_litr3c;
		cs_litr->litr4c    += m_deadleafc_to_litr4c;		
		cs_litr->litr1c    += m_leafc_store_to_litr1c;
		cs_litr->litr1c    += m_leafc_transfer_to_litr1c;
		if (epc.veg_type == TREE) {
			/*    Stem wood mortality */
			/*	  Transfer to DEADWOOD if standing dead */
			if (thintyp == 3) {
				cs->dead_stemc       += m_livestemc_to_cwdc;
				cs->dead_stemc       += m_deadstemc_to_cwdc;
				}
			/*	  Transfer to CWD otherwise */
			else {
				cs->cwdc       += m_livestemc_to_cwdc;
				cs->cwdc       += m_deadstemc_to_cwdc;
				}
			cs_litr->litr1c    += m_livestemc_store_to_litr1c;
			cs_litr->litr1c    += m_deadstemc_store_to_litr1c;
			cs_litr->litr1c    += m_livestemc_transfer_to_litr1c;
			cs_litr->litr1c    += m_deadstemc_transfer_to_litr1c;			
        }
		/* gresp... group in with aboveground? */
		cs_litr->litr1c         += m_gresp_store_to_litr1c;
		cs_litr->litr1c         += m_gresp_transfer_to_litr1c;
    }
	/* Remove aboveground dead c from carbon stores in all cases. */
	cs->cpool -= m_cpool;
	/*    Leaf mortality */
	cs->leafc          -= m_leafc_to_litr1c;
	cs->leafc          -= m_leafc_to_litr2c;
	cs->leafc          -= m_leafc_to_litr3c;
	cs->leafc          -= m_leafc_to_litr4c;
	cs->dead_leafc          -= m_deadleafc_to_litr1c;
	cs->dead_leafc          -= m_deadleafc_to_litr2c;
	cs->dead_leafc          -= m_deadleafc_to_litr3c;
	cs->dead_leafc          -= m_deadleafc_to_litr4c;	
	cs->leafc_store       -= m_leafc_store_to_litr1c;
	cs->leafc_transfer      -= m_leafc_transfer_to_litr1c;
	if (epc.veg_type == TREE){
		/*    Stem wood mortality */
		cs->live_stemc  -= m_livestemc_to_cwdc;
		cs->dead_stemc  -= m_deadstemc_to_cwdc;
		cs->livestemc_store   -= m_livestemc_store_to_litr1c;
		cs->deadstemc_store   -= m_deadstemc_store_to_litr1c;
		cs->livestemc_transfer  -= m_livestemc_transfer_to_litr1c;
		cs->deadstemc_transfer  -= m_deadstemc_transfer_to_litr1c;	
		}
	/* gresp... group in with aboveground? */
	cs->gresp_store       -= m_gresp_store_to_litr1c;
	cs->gresp_transfer      -= m_gresp_transfer_to_litr1c;
	
	/* BELOWGROUND C POOLS */
	/* Belowground dead c goes to litter and cwd in all cases. */
	/*   Fine root mortality */
	cs_litr->litr1c    += m_frootc_to_litr1c;
	cs->frootc         -= m_frootc_to_litr1c;
	cs_litr->litr2c    += m_frootc_to_litr2c;
	cs->frootc         -= m_frootc_to_litr2c;
	cs_litr->litr3c    += m_frootc_to_litr3c;
	cs->frootc         -= m_frootc_to_litr3c;
	cs_litr->litr4c    += m_frootc_to_litr4c;
	cs->frootc         -= m_frootc_to_litr4c;
	cs_litr->litr1c    += m_frootc_store_to_litr1c;
	cs->frootc_store   -= m_frootc_store_to_litr1c;
	cs_litr->litr1c         += m_frootc_transfer_to_litr1c;
	cs->frootc_transfer     -= m_frootc_transfer_to_litr1c;
	if (epc.veg_type == TREE){
		/* Coarse root wood mortality */
		cs->cwdc       += m_livecrootc_to_cwdc;
		cs->live_crootc -= m_livecrootc_to_cwdc;
		cs->cwdc       += m_deadcrootc_to_cwdc;
		cs->dead_crootc -= m_deadcrootc_to_cwdc;		
		cs_litr->litr1c       += m_livecrootc_store_to_litr1c;
		cs->livecrootc_store  -= m_livecrootc_store_to_litr1c;
		cs_litr->litr1c       += m_deadcrootc_store_to_litr1c;
		cs->deadcrootc_store  -= m_deadcrootc_store_to_litr1c;
		cs_litr->litr1c         += m_livecrootc_transfer_to_litr1c;
		cs->livecrootc_transfer -= m_livecrootc_transfer_to_litr1c;
		cs_litr->litr1c         += m_deadcrootc_transfer_to_litr1c;
		cs->deadcrootc_transfer -= m_deadcrootc_transfer_to_litr1c;
    }
	
	/* ---------------------------------------- */
	/* NITROGEN mortality state variable update */
	/* ---------------------------------------- */
	
    /* ABOVEGROUND N POOLS */
	
	/* Only add dead leaf and stem n to litter and cwd pools if thintyp   */
	/* is not 2 (harvest case). If thintyp is 2, harvest aboveground n.   */
	if (thintyp != 2) {
//        ns_litr->litr1n    += m_cpool>0? max(0.0, m_npool) : 0.0;
		/*    Leaf mortality */
		ns_litr->litr1n    += m_leafn_to_litr1n;
		ns_litr->litr2n    += m_leafn_to_litr2n;
		ns_litr->litr3n    += m_leafn_to_litr3n;
		ns_litr->litr4n    += m_leafn_to_litr4n;
		ns_litr->litr1n    += m_deadleafn_to_litr1n;
		ns_litr->litr2n    += m_deadleafn_to_litr2n;
		ns_litr->litr3n    += m_deadleafn_to_litr3n;
		ns_litr->litr4n    += m_deadleafn_to_litr4n;		
		ns_litr->litr1n    += m_leafn_store_to_litr1n;
		ns_litr->litr1n    += m_leafn_transfer_to_litr1n;
		ns_litr->litr1n    += m_retransn_to_litr1n;
		if (epc.veg_type == TREE){
            /*    Stem wood mortality */
			ns_litr->litr1n     += m_livestemn_to_litr1n;
			if (thintyp != 3) {
                /*      Transfer to CWD if normal thinning */
				ns->cwdn       += m_livestemn_to_cwdn;
				ns->cwdn       += m_deadstemn_to_cwdn;
            }else{
                /*      Transfer to DEADWOOD if standing dead */
				ns->dead_stemn       += m_livestemn_to_cwdn;
				ns->dead_stemn       += m_deadstemn_to_cwdn;
            }
			ns_litr->litr1n    += m_livestemn_store_to_litr1n;
			ns_litr->litr1n    += m_deadstemn_store_to_litr1n;
			ns_litr->litr1n    += m_livestemn_transfer_to_litr1n;
			ns_litr->litr1n    += m_deadstemn_transfer_to_litr1n;
        }// end of TREE if
    }// end of thintyp != 2
    
	/* Remove aboveground dead n from n stores in all cases. */
	ns->npool -= m_npool;
	/*    Leaf mortality */
	ns->leafn          -= m_leafn_to_litr1n;
	ns->leafn          -= m_leafn_to_litr2n;
	ns->leafn          -= m_leafn_to_litr3n;
	ns->leafn          -= m_leafn_to_litr4n;
	ns->dead_leafn          -= m_deadleafn_to_litr1n;
	ns->dead_leafn          -= m_deadleafn_to_litr2n;
	ns->dead_leafn          -= m_deadleafn_to_litr3n;
	ns->dead_leafn          -= m_deadleafn_to_litr4n;	
	ns->leafn_store       -= m_leafn_store_to_litr1n;
	ns->leafn_transfer      -= m_leafn_transfer_to_litr1n;
	ns->retransn            -= m_retransn_to_litr1n;
	if (epc.veg_type == TREE){
		/*    Stem wood mortality */
		ns->live_stemn  -= m_livestemn_to_litr1n;
		ns->live_stemn  -= m_livestemn_to_cwdn;
		ns->dead_stemn  -= m_deadstemn_to_cwdn;
		ns->livestemn_store   -= m_livestemn_store_to_litr1n;
		ns->deadstemn_store   -= m_deadstemn_store_to_litr1n;
		ns->livestemn_transfer  -= m_livestemn_transfer_to_litr1n;
		ns->deadstemn_transfer  -= m_deadstemn_transfer_to_litr1n;
    }

	/* BELOWGROUND N POOLS */
	/* Belowground dead n goes to litter and cwd in all cases. */
	/*   Fine root mortality */
	ns_litr->litr1n    += m_frootn_to_litr1n;
	ns->frootn         -= m_frootn_to_litr1n;
	ns_litr->litr2n    += m_frootn_to_litr2n;
	ns->frootn         -= m_frootn_to_litr2n;
	ns_litr->litr3n    += m_frootn_to_litr3n;
	ns->frootn         -= m_frootn_to_litr3n;
	ns_litr->litr4n    += m_frootn_to_litr4n;
	ns->frootn         -= m_frootn_to_litr4n;
	ns_litr->litr1n         += m_frootn_store_to_litr1n;
	ns->frootn_store      -= m_frootn_store_to_litr1n;
	ns_litr->litr1n         += m_frootn_transfer_to_litr1n;
	ns->frootn_transfer     -= m_frootn_transfer_to_litr1n;
	if (epc.veg_type == TREE){
		/* Coarse root mortality */
		ns_litr->litr1n     += m_livecrootn_to_litr1n;
		ns->live_crootn -= m_livecrootn_to_litr1n;
		ns->cwdn       += m_livecrootn_to_cwdn;
		ns->live_crootn -= m_livecrootn_to_cwdn;
		ns->cwdn       += m_deadcrootn_to_cwdn;
		ns->dead_crootn -= m_deadcrootn_to_cwdn;
		ns_litr->litr1n         += m_livecrootn_store_to_litr1n;
		ns->livecrootn_store  -= m_livecrootn_store_to_litr1n;
		ns_litr->litr1n         += m_deadcrootn_store_to_litr1n;
		ns->deadcrootn_store  -= m_deadcrootn_store_to_litr1n;
		ns_litr->litr1n         += m_livecrootn_transfer_to_litr1n;
		ns->livecrootn_transfer -= m_livecrootn_transfer_to_litr1n;
		ns_litr->litr1n         += m_deadcrootn_transfer_to_litr1n;
		ns->deadcrootn_transfer -= m_deadcrootn_transfer_to_litr1n;		
	}
    
//    printf("update mortality (after plant cs): %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
//           cs->leafc, cs->dead_leafc, cs->frootc,
//           cs->live_stemc, cs->dead_stemc, cs->live_crootc, cs->dead_crootc);
//    printf("update mortality (after plant ns): %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
//           ns->leafn, ns->dead_leafn, ns->frootn,
//           ns->live_stemn, ns->dead_stemn, ns->live_crootn, ns->dead_crootn);
    
//    printf("update mortality daily (before): %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
//           epc.leaflitr_fscel,epc.leaflitr_fscel,epc.frootlitr_fscel,
//           m_leafc_to_litr3c,m_deadleafc_to_litr3c,m_frootc_to_litr3c,
//           m_leafn_to_litr3n,m_deadleafn_to_litr3n,m_frootn_to_litr3n);
//    printf("update branch mortality (after cs): %lf, %lf, %lf, %lf\n",
//           cs_litr->litr1c, cs_litr->litr2c,
//           cs_litr->litr3c, cs_litr->litr4c);
//    printf("update branch mortality (after ns): %lf, %lf, %lf, %lf\n",
//           ns_litr->litr1n, ns_litr->litr2n,
//           ns_litr->litr3n, ns_litr->litr4n);
//
//    printf("update branch mortality (END): %lf, %lf, %lf, %lf\n",
//           cs_litr->litr1c/ns_litr->litr1n, cs_litr->litr2c/ns_litr->litr2n,
//           cs_litr->litr3c/ns_litr->litr3n, cs_litr->litr4c/ns_litr->litr4n);
    

    if( leaffrootLiter1c/leaffrootLiter1n < 18 || (stemcrootLiter1c/stemcrootLiter1n < 18 && stemcrootLiter1c>0)){
        printf("patch %d update mortality[%d]BGC_flag[%d], Begin_litterpool_CN[%e,%e,%e,%e], leaffrootLiter1c[%e],%e, leaf(%e~%e)={%e}, leafstorage(%e)={%e}, leaftransfer(%e)={%e}, froot(%e)={%e}, frootstorage(%e)={%e},froottransfer(%e)={%e},dleaf(%e)={%e}, stemcroot(%e)={%e}, resp=(%e,%e,%e,%e)\n",
               patch[0].ID, ID, BGC_flag,
               cn_l1,cn_l2,cn_l3,cn_l4,//l1CN
               (leaffrootLiter1c/leaffrootLiter1n), leaffrootLiter1c, //leaffroot (being negative or too low)
                    (m_leafc_to_litr1c/m_leafn_to_litr1n), epc.leaf_cn, m_leafc_to_litr1c, //leaf
                    (m_leafc_store_to_litr1c/m_leafn_store_to_litr1n),m_leafc_store_to_litr1c,//leafstorage
                    (m_leafc_transfer_to_litr1c)/(m_leafn_transfer_to_litr1n),m_leafc_transfer_to_litr1c,//leaftransfer
                    (m_frootc_to_litr1c/m_frootn_to_litr1n), m_frootc_to_litr1c,//froot
                    (m_frootc_store_to_litr1c)/(m_frootn_store_to_litr1n),m_frootc_store_to_litr1c,//frootst
                    (m_frootc_transfer_to_litr1c)/(m_frootn_transfer_to_litr1n),m_frootc_transfer_to_litr1c,//frootst
                    (m_deadleafc_to_litr1c/m_deadleafn_to_litr1n),m_deadleafc_to_litr1c, //dleaf
               stemcrootLiter1c/stemcrootLiter1n, stemcrootLiter1c, //stemcroot
               cs->gresp_store,m_gresp_store_to_litr1c, cs->gresp_transfer,m_gresp_transfer_to_litr1c
               
               );
    }//debug
    
//    if( patch[0].ID == 209782 && ID==804){
//        printf("[%ld] mortality@(%d, %d):\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n[%e,%e,%e,%e]\n",
//               day, patch[0].ID, ID,
//               mort.mort_leafc, cs->leafc, cs->leafc_store, cs->leafc_transfer, cs->leafc/ns->leafn, cs->leafc_store/ns->leafn_store, cs->leafc_transfer/ns->leafn_transfer,
//               mort.mort_frootc, cs->frootc, cs->frootc_store, cs->frootc_transfer, cs->frootc/ns->frootn, cs->frootc_store/ns->frootn_store, cs->frootc_transfer/ns->frootn_transfer,
//
//               mort.mort_livestemc, cs->live_stemc, cs->livestemc_store, cs->livestemc_transfer, cs->live_stemc/ns->live_stemn, cs->livestemc_store/ns->livestemn_store, cs->livestemc_transfer/ns->livestemn_transfer,
//               mort.mort_deadstemc, cs->dead_stemc, cs->deadstemc_store, cs->deadstemc_transfer, cs->dead_stemc/ns->dead_stemn, cs->deadstemc_store/ns->deadstemn_store, cs->deadstemc_transfer/ns->deadstemn_transfer,
//
//               mort.mort_livecrootc, cs->live_crootc, cs->livecrootc_store, cs->livecrootc_transfer, cs->live_crootc/ns->live_crootn, cs->livecrootc_store/ns->livecrootn_store, cs->livecrootc_transfer/ns->livecrootn_transfer,
//               mort.mort_deadcrootc, cs->dead_crootc, cs->deadcrootc_store, cs->deadcrootc_transfer, cs->dead_crootc/ns->dead_crootn, cs->deadcrootc_store/ns->deadcrootn_store, cs->deadcrootc_transfer/ns->deadcrootn_transfer,
//               leaffrootLiter1c/leaffrootLiter1n, leaffrootLiter1c,
//               stemcrootLiter1c/stemcrootLiter1n, stemcrootLiter1c
//               );
//    }
//    if( patch[0].ID == 209782 && ID==802){
//        printf("[%ld] mortality@(%d, %d):\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n%e, %e,%e,%e [%e,%e,%e]\n[%e,%e,%e,%e]\n\n",
//               day, patch[0].ID, ID,
//               mort.mort_leafc, cs->leafc, cs->leafc_store, cs->leafc_transfer, cs->leafc/ns->leafn, cs->leafc_store/ns->leafn_store, cs->leafc_transfer/ns->leafn_transfer,
//               mort.mort_frootc, cs->frootc, cs->frootc_store, cs->frootc_transfer, cs->frootc/ns->frootn, cs->frootc_store/ns->frootn_store, cs->frootc_transfer/ns->frootn_transfer,
//
//               mort.mort_livestemc, cs->live_stemc, cs->livestemc_store, cs->livestemc_transfer, cs->live_stemc/ns->live_stemn, cs->livestemc_store/ns->livestemn_store, cs->livestemc_transfer/ns->livestemn_transfer,
//               mort.mort_deadstemc, cs->dead_stemc, cs->deadstemc_store, cs->deadstemc_transfer, cs->dead_stemc/ns->dead_stemn, cs->deadstemc_store/ns->deadstemn_store, cs->deadstemc_transfer/ns->deadstemn_transfer,
//
//               mort.mort_livecrootc, cs->live_crootc, cs->livecrootc_store, cs->livecrootc_transfer, cs->live_crootc/ns->live_crootn, cs->livecrootc_store/ns->livecrootn_store, cs->livecrootc_transfer/ns->livecrootn_transfer,
//               mort.mort_deadcrootc, cs->dead_crootc, cs->deadcrootc_store, cs->deadcrootc_transfer, cs->dead_crootc/ns->dead_crootn, cs->deadcrootc_store/ns->deadcrootn_store, cs->deadcrootc_transfer/ns->deadcrootn_transfer,
//               leaffrootLiter1c/leaffrootLiter1n, leaffrootLiter1c,
//               stemcrootLiter1c/stemcrootLiter1n, stemcrootLiter1c
//               );
//    }
    
//    if ((cs_litr->litr1c > 0.0) && (ns_litr->litr1n > 0.0))
//        if( cn_l1>5 && cs_litr->litr1c/ns_litr->litr1n<cn_l1 && cs_litr->litr1c/ns_litr->litr1n < 5.0) printf("update mortality[%d], l1CN = %e(%e), %e, %e\n",ID,cn_l1, cs_litr->litr1c/ns_litr->litr1n, cs_litr->litr1c, ns_litr->litr1n);
//    
	return;
}/*end update_mortality*/

