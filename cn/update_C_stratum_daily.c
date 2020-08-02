/*--------------------------------------------------------------*/
/* 								*/
/*		update_C_stratum_daily				*/
/*								*/
/*								*/
/*	NAME							*/
/*	update_C_stratum_daily -					*/
/*								*/
/*		updates daily C stores to			*/
/*		account for psn and respiration 		*/
/*	SYNOPSIS						*/
/*	double	update_C_stratum_daily(					*/
/*								*/
/*	returns:						*/
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

int update_C_stratum_daily(struct epconst_struct epc,
						   struct cstate_struct *cs,
                           struct nstate_struct *ns,
                           struct cdayflux_struct *cdf,
                           struct soil_c_object *cs_soil,
                           struct canopy_strata_object *stratum,
                           struct date current_date,
                           struct patch_object *patch)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
    void    update_mortality(
        struct epconst_struct,
        struct cstate_struct *,
        struct cdayflux_struct *,
        struct cdayflux_patch_struct *,
        struct nstate_struct *,
        struct ndayflux_struct *,
        struct ndayflux_patch_struct *,
        struct litter_c_object *,
        struct litter_n_object *,
        struct date,
        struct patch_object *,
        struct canopy_strata_object *stratum,
        int,
        struct mortality_struct,
        int,
        int);
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	int ok=1;
    double checkCpool = cs->cpool;
	
	/* Daily photosynthesis */
	cs->gpsn_src	+= cdf->psn_to_cpool;
	cs->cpool	+= cdf->psn_to_cpool;
	/* Daily maintenance respiration */
	cs->cpool	-= cdf->leaf_day_mr;
	cs->cpool	-= cdf->leaf_night_mr;
	cs->leaf_mr_snk += cdf->leaf_day_mr + cdf->leaf_night_mr;
	cs->cpool	-= cdf->froot_mr;
	cs->froot_mr_snk += cdf->froot_mr;
	if (epc.veg_type == TREE) {
		cs->cpool	-= cdf->livestem_mr;
		cs->cpool	-= cdf->livecroot_mr;
		cs->livestem_mr_snk += cdf->livestem_mr;
		cs->livecroot_mr_snk += cdf->livecroot_mr;
	}
	cs->net_psn = cdf->psn_to_cpool - cdf->total_mr - cdf->total_gr;

      
    
    
    double total_store = 0.0;
    double TOtotal_X = 0.0;
    double TOtotal_store = 0.0;
	/* Daily allocation fluxes */
	/* daily leaf allocation fluxes */
	cs->leafc          += cdf->cpool_to_leafc;
	cs->cpool          -= cdf->cpool_to_leafc;
	cs->leafc_store  += cdf->cpool_to_leafc_store;
	cs->cpool          -= cdf->cpool_to_leafc_store;
    total_store += cs->leafc_store;
    TOtotal_X += cdf->cpool_to_leafc;
    TOtotal_store += cdf->cpool_to_leafc_store;
	/* Daily fine root allocation fluxes */
	cs->frootc         += cdf->cpool_to_frootc;
	cs->cpool          -= cdf->cpool_to_frootc;
	cs->frootc_store += cdf->cpool_to_frootc_store;
	cs->cpool          -= cdf->cpool_to_frootc_store;
    total_store += cs->frootc_store;
    TOtotal_X += cdf->cpool_to_frootc;
    TOtotal_store += cdf->cpool_to_frootc_store;
	if (epc.veg_type == TREE){
		/* Daily live stem wood allocation fluxes */
		cs->live_stemc          += cdf->cpool_to_livestemc;
		cs->cpool              -= cdf->cpool_to_livestemc;
		cs->livestemc_store  += cdf->cpool_to_livestemc_store;
		cs->cpool              -= cdf->cpool_to_livestemc_store;
        total_store += cs->livestemc_store;
        TOtotal_X += cdf->cpool_to_livestemc;
        TOtotal_store += cdf->cpool_to_livestemc_store;
		/* Daily dead stem wood allocation fluxes */
		cs->dead_stemc          += cdf->cpool_to_deadstemc;
		cs->cpool              -= cdf->cpool_to_deadstemc;
		cs->deadstemc_store  += cdf->cpool_to_deadstemc_store;
		cs->cpool              -= cdf->cpool_to_deadstemc_store;
        total_store += cs->deadstemc_store;
        TOtotal_X += cdf->cpool_to_deadstemc;
        TOtotal_store += cdf->cpool_to_deadstemc_store;
		/* Daily live coarse root wood allocation fluxes */
		cs->live_crootc         += cdf->cpool_to_livecrootc;
		cs->cpool              -= cdf->cpool_to_livecrootc;
		cs->livecrootc_store += cdf->cpool_to_livecrootc_store;
		cs->cpool              -= cdf->cpool_to_livecrootc_store;
        total_store += cs->livecrootc_store;
        TOtotal_X += cdf->cpool_to_livecrootc;
        TOtotal_store += cdf->cpool_to_livecrootc_store;
		/* Daily dead coarse root wood allocation fluxes */
		cs->dead_crootc         += cdf->cpool_to_deadcrootc;
		cs->cpool              -= cdf->cpool_to_deadcrootc;
		cs->deadcrootc_store += cdf->cpool_to_deadcrootc_store;
		cs->cpool              -= cdf->cpool_to_deadcrootc_store;
        total_store += cs->deadcrootc_store;
        TOtotal_X += cdf->cpool_to_deadcrootc;
        TOtotal_store += cdf->cpool_to_deadcrootc_store;
	}
	/* Daily allocation for transfer growth respiration */
    cs->gresp_store  += cdf->cpool_to_gresp_store; // allocate and save for the future second growth respiration (from allocate_daily_growth; real time)
	cs->cpool        -= cdf->cpool_to_gresp_store;
	/* Daily growth respiration fluxes */
	/* Leaf growth respiration */
	cs->leaf_gr_snk     += cdf->cpool_leaf_gr; // tracking
	cs->cpool           -= cdf->cpool_leaf_gr; // first growth respiration (real time)
	cs->leaf_gr_snk     += cdf->transfer_leaf_gr; // tracking
	cs->gresp_transfer  -= cdf->transfer_leaf_gr; // second growth respiration (real time)
	/* Fine root growth respiration */
	cs->froot_gr_snk    += cdf->cpool_froot_gr;
	cs->cpool           -= cdf->cpool_froot_gr;
	cs->froot_gr_snk    += cdf->transfer_froot_gr;
	cs->gresp_transfer  -= cdf->transfer_froot_gr;
	if (epc.veg_type == TREE){
		/* Live stem growth respiration */
		cs->livestem_gr_snk  += cdf->cpool_livestem_gr;
		cs->cpool            -= cdf->cpool_livestem_gr;
		cs->livestem_gr_snk  += cdf->transfer_livestem_gr;
		cs->gresp_transfer   -= cdf->transfer_livestem_gr;
		/* Dead stem growth respiration */
		cs->deadstem_gr_snk  += cdf->cpool_deadstem_gr;
		cs->cpool            -= cdf->cpool_deadstem_gr;
		cs->deadstem_gr_snk  += cdf->transfer_deadstem_gr;
		cs->gresp_transfer   -= cdf->transfer_deadstem_gr;
		/* Live coarse root growth respiration */
		cs->livecroot_gr_snk += cdf->cpool_livecroot_gr;
		cs->cpool            -= cdf->cpool_livecroot_gr;
		cs->livecroot_gr_snk += cdf->transfer_livecroot_gr;
		cs->gresp_transfer   -= cdf->transfer_livecroot_gr;
		/* Dead coarse root growth respiration */
		cs->deadcroot_gr_snk += cdf->cpool_deadcroot_gr;
		cs->cpool            -= cdf->cpool_deadcroot_gr;
		cs->deadcroot_gr_snk += cdf->transfer_deadcroot_gr;
		cs->gresp_transfer   -= cdf->transfer_deadcroot_gr;
	}

    if( cs->gresp_transfer < 0.0){
        cs->cpool += cs->gresp_transfer;
        cs->gresp_transfer = 0.0;
    }
    
    
    // self regulating, maint. resp. is too high
    // this mess up things
    if(cs->cpool<0 && cdf->total_mr>0){
        double reduction = -cs->cpool/cdf->total_mr;
        double reduced = 0.0;
        double leafred = min(cs->leafc_store*0.5, reduction*(cdf->leaf_day_mr + cdf->leaf_night_mr));
        reduced += leafred;
        if(leafred>0) leafred /= cs->leafc_store;
        double frootred = min(cs->frootc_store*0.5, reduction*cdf->froot_mr);
        reduced += frootred;
        if(frootred>0) frootred /= cs->frootc_store;
        double stemred = min(cs->livestemc_store*0.5, reduction*cdf->livestem_mr);
        reduced += stemred;
        if(stemred>0) stemred /= cs->livestemc_store;
        double crootred = min(cs->livecrootc_store*0.5, reduction*cdf->livecroot_mr);
        reduced += crootred;
        if(crootred>0) crootred /= cs->livecrootc_store;
        
        cs->leafc_store *= 1.0-leafred;
        cs->frootc_store *= 1.0-frootred;
        cs->livestemc_store *= 1.0-stemred;
        cs->livecrootc_store *= 1.0-crootred;
        ns->npool += ns->leafn_store *leafred;
        ns->npool += ns->frootn_store * frootred;
        ns->npool += ns->livestemn_store * stemred;
        ns->npool += ns->livecrootn_store * crootred;
        ns->leafn_store *= 1.0-leafred;
        ns->frootn_store *= 1.0-frootred;
        ns->livestemn_store *= 1.0-stemred;
        ns->livecrootn_store *= 1.0-crootred;

        // adjusting biomass for resp.
        if(reduction > 0 && cs->leafc+cs->live_stemc+cs->live_crootc+cs->frootc>0){
            // cdf->froot_mr = ns->frootn * mrpern * t1;
            // ns->frootn = (cs->frootc / epc.froot_cn)
            // cdf->froot_mr = (cs->frootc / epc.froot_cn) * mrpern * t1;
            // cdf->froot_mr' = (cs->frootc' / epc.froot_cn) * mrpern * t1;
            // cs->frootc' / cs->frootc = cdf->froot_mr' / cdf->froot_mr
            // mort rate = 1 - cs->frootc' / cs->frootc = 1 - cdf->froot_mr' / cdf->froot_mr
            
            // reduction*cdf->froot_mr = extra resp.
            // cdf->froot_mr' = cdf->froot_mr - extra resp.
            // 1 - cdf->froot_mr' / cdf->froot_mr = 1 - (cdf->froot_mr - extra resp.) / cdf->froot_mr
            // = 1 - (cdf->froot_mr - reduction * cdf->froot_mr) / cdf->froot_mr
            // = 1 - (1 - reduction) = reduction
            reduction = min(0.1, reduction);
            struct mortality_struct mort;
            mort.mort_cpool = 0.0; // some how it does not set to zero!!
            mort.mort_leafc = reduction;
            mort.mort_deadleafc = 0.0;
            mort.mort_livestemc = reduction;
            mort.mort_deadstemc = 0.0;
            mort.mort_livecrootc = reduction;
            mort.mort_deadcrootc = 0.0;
            mort.mort_frootc = reduction;

            //update_mortality() resets all flux related to liter and updates the liter and plantX variables
            update_mortality(
                stratum[0].defaults[0][0].epc,
                &(stratum[0].cs),
                &(stratum[0].cdf),
                &(patch[0].cdf),
                &(stratum[0].ns),
                &(stratum[0].ndf),
                &(patch[0].ndf),
                &(patch[0].litter_cs),
                &(patch[0].litter_ns),
                current_date,
                patch,
                stratum,
                1,
                mort,
                1,
                stratum[0].defaults[0][0].ID
                );
        }
        cs->cpool += reduced;
    }// end of self regulating
    
    // checking
    if(cs->cpool<0 && checkCpool>0 ){ //&& current_date.year>1950
        printf("update C stratum daily %d %d %d %d %d %d %e = %e +%e -%e -%e -%e -%e -%e\n",
               patch[0].ID,
               current_date.day, current_date.month, current_date.year,
               stratum[0].ID,
               stratum->defaults[0][0].ID,
               cs->cpool,
               checkCpool, cdf->psn_to_cpool,
               cdf->total_gr, cdf->total_mr, cdf->cpool_to_gresp_store, TOtotal_store, TOtotal_X);
    }
    
    // drain out excessive cpool, cap at 80% of stores
    if(stratum[0].phen.gwseasonday > epc.ndays_expand && cs->cpool-total_store*0.8>0){
//        if(patch[0].ID==42534) printf("update C stratum daily %d %d %d %d %d %d %e %e\n",
//               patch[0].ID,
//               current_date.day, current_date.month, current_date.year,
//               stratum[0].ID,
//               stratum->defaults[0][0].ID,
//               cs->cpool,
//               total_store);
        total_store *= 0.8;
        total_store -= cs->cpool;
        cs_soil->DOC -= total_store*0.01;
        cs->cpool += total_store*0.01;
    }//end of if
    
	return (!ok);
}		

