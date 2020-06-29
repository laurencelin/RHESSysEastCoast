/*--------------------------------------------------------------*/
/* 								*/
/*								*/
/*	allocate_daily_growth				*/
/*								*/
/*	NAME							*/
/*		allocate_daily_growth			*/
/*								*/
/*	SYNOPSIS						*/
/*	int allocate_daily_growth( int , double, double, double,double,	*/
/*			    struct cdayflux_struct *,  		*/
/*			    struct cstate_struct *,		*/
/*			    struct ndayflux_struct *,	 	*/
/*			    struct nstate_struct *, 		*/
/*			    struct ndayflux_patch_struct *, 	*/
/*			    struct epconst_struct)		*/
/*								*/
/*	returns:						*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*		calculates daily C and N allocations		*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*								*/
/*              modified from Peter Thornton (1998)             */
/*                      dynamic - 1d-bgc ver4.0                 */
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

int allocate_daily_growth(int nlimit,
						  double pnow,
						  double Tsoil,
						  double total_soil_frootc,
						  double cover_fraction,
						  struct cdayflux_struct *cdf,
						  struct cstate_struct *cs,
						  struct ndayflux_struct *ndf,
						  struct nstate_struct *ns,
						  struct ndayflux_patch_struct *ndf_patch,
						  struct epvar_struct *epv,
						  struct epconst_struct epc,
						  struct date current_date,
                          struct patch_object *patch,
                          struct canopy_strata_object *stratum,
                          struct command_line_object *command_line)
{
	/*------------------------------------------------------*/
	/*	Local function declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	
	int ok=1;
	double fleaf;          /* RATIO   new leaf C: new total C     */
	double froot;          /* RATIO   new fine root C : new total C     */
	double flive, fdead;	/* RATIO  live/dead C : new total C */
	double fwood;          /* RATIO   wood          */
    double fstem, fcroot;
	double f3;		
	double g1;          /* RATIO   C respired for growth : C grown  */
	double cnl;         /* RATIO   leaf C:N      */
	double cnfr;        /* RATIO   fine root C:N */
	double cnlw;        /* RATIO   live wood C:N */
	double cndw;        /* RATIO   dead wood C:N */
	//double nlc;         /* actual new leaf C, minimum of C and N limits   */
	double amt_fix, cost_fix, closs;
	double gresp_store, total_wood;
	double mean_cn;
	double sum_plant_nsupply, soil_nsupply;
	double plant_nalloc=0.0;
	double plant_calloc;
	double plant_remaining_ndemand=0.0;
	double excess_allocation_to_leaf, excess_c, excess_lai;
	double sminn_to_npool;
	double totalc_used;
	

	/* assign local values for the allocation control parameters */
	flive = epc.alloc_livewoodc_woodc;
	f3 = epc.alloc_stemc_leafc;
	fdead = (1.0-flive);
	g1 = epc.gr_perc;
	cnl = epc.leaf_cn;
	cnfr = epc.froot_cn;
	cnlw = epc.livewood_cn;
	cndw = epc.deadwood_cn;
	excess_c = 0.0;
	sminn_to_npool = 0.0;


//    if(patch[0].ID == 139835 && stratum->defaults[0][0].ID==102){
//        printf("allocate,%d,%d,%d,%d,%d,%d\n",
//               patch[0].ID, stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year,
//               stratum[0].phen.gwseasonday
//               );
//    }
    if(stratum[0].phen.gwseasonday>0){
        stratum[0].gDayCount ++;
        stratum[0].nFactor += patch[0].soil_ns.fract_potential_uptake;
    }//growth season
	/*--------------------------------------------------------------*/
	/*	allocation partitioning			*/
	/* computed in compute_N_uptake routines for allocation specific models */
	/*--------------------------------------------------------------*/
	flive = epc.alloc_livewoodc_woodc;
	fdead = 1-flive;
	fleaf = cdf->fleaf;
	froot = cdf->froot;
	fwood = cdf->fwood;
	fstem = cdf->fstem;
    fcroot = cdf->fcroot;
    
	if ((fleaf + froot) > ZERO) {
        if (epc.veg_type == TREE){
            mean_cn = (fleaf+froot+fwood)>0? (fleaf+froot+fwood)/(fleaf/cnl + froot/cnfr + flive*fwood/cnlw + fwood*fdead/cndw) : 0.0;
            /*mean_cn = fleaf * cnl + froot * cnfr + flive * fwood * cnlw + fdead * fwood * cndw;*/
        } else{
            mean_cn = (fleaf+froot)>0? (fleaf+froot)/(fleaf/cnl + froot/cnfr) : 0.0;
        }
	} else mean_cn = 0.0;


	// note that cs->availc is calculated in compute_potential_N_uptakeXXX = cdf->psn_to_cpool-cdf->total_mr + ifelse(cpool<0, cpool, 0);
    // ndf_patch->plant_potential_ndemand is cs->availc / (1.0+epc.gr_perc) / mean_cn @ potential_uptake_XXX.c
    if(ndf_patch->plant_potential_ndemand>0 && ndf->potential_N_uptake>0){
        if (nlimit == 1){
            // Yes patch N is limited

            
            // no cover fraction here because it's unit meter. cover fraction is used when communicate to patch level
            soil_nsupply = patch[0].soil_ns.fract_potential_uptake * ndf->potential_N_uptake; //(same as below)
            //soil_nsupply = ndf_patch->plant_avail_uptake / ndf_patch->plant_potential_ndemand * ndf->potential_N_uptake; // no biomass assumption; fairly dividied
            
            if(soil_nsupply!=soil_nsupply || soil_nsupply<0)
                printf("allocate_daily_growth has negative [%d: %d,%d,%d] Nlimit=1 (%e, %e, %e, %e)\n",
                       stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year,
                       soil_nsupply,
                       ndf_patch->plant_avail_uptake, // it negative, it's in resolve_sminn_competition.c
                       ndf_patch->plant_potential_ndemand,
                       ndf->potential_N_uptake); //debug

            sminn_to_npool = soil_nsupply;
            plant_remaining_ndemand = ndf->potential_N_uptake - sminn_to_npool;

            if (plant_remaining_ndemand <= ns->retransn){
                ndf->retransn_to_npool = plant_remaining_ndemand;
                plant_nalloc = ndf->retransn_to_npool + sminn_to_npool;
                plant_calloc = cs->availc/(1.0+epc.gr_perc); // correct for less growth
                ns->nlimit = 0;
            }else{
                ndf->retransn_to_npool = ns->retransn;
                plant_nalloc = ndf->retransn_to_npool + sminn_to_npool;//<---
                plant_calloc = min(cs->availc/(1.0+epc.gr_perc), plant_nalloc * mean_cn);
                if(plant_calloc>0 && (fleaf+froot+fwood)>0) plant_calloc /= (fleaf+froot+fwood);//correction
                excess_c = max(cs->availc - (plant_calloc*(1.0+epc.gr_perc)),0.0); // adjust to cdf->psn_to_cpool

                if  (epc.nfix == 1){
                    // disable for now; fix this later
//                    // N fixing: the code does not make sense!
//                    cost_fix = -0.625*(exp(-3.62 + 0.27 * Tsoil*(1 - 0.5 * Tsoil / 25.15)) - 2);// reference?
//                    if (cost_fix > ZERO)
//                        amt_fix = 0.5*cost_fix * excess_c / mean_cn; // N or C?
//                    else
//                        amt_fix = 0.0;
//
//                    //amt_fix = min(excess_c, amt_fix);
//                    plant_calloc = plant_calloc + excess_c - amt_fix; //non-sense
//                    plant_nalloc = plant_calloc/mean_cn;
//
//                    ndf_patch->nfix_to_sminn = plant_nalloc - ndf->retransn_to_npool - sminn_to_npool;
//                    excess_c = excess_c - amt_fix;//
//                    if (excess_c > ZERO) {
//                        cdf->psn_to_cpool -= excess_c;
//                        ns->nlimit = 1;
//                    }else{
//                        ns->nlimit=0;
//                    }
                }else {
                    // no N-fix
                    cdf->psn_to_cpool -= excess_c; // if not enough N, it will reduce fraq_psn result.
                    ns->nlimit = 1;
                }// N fixer
            }//ns->retransn

        }else{
            // patch N is enough!
            soil_nsupply = ndf->potential_N_uptake;
            if(soil_nsupply!=soil_nsupply || soil_nsupply<0)
                printf("allocate_daily_growth has negative [%d: %d,%d,%d] Nlimit=0 (%e, %e, %e, %e)\n",
                       stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year,
                       soil_nsupply,
                       ndf_patch->plant_avail_uptake,
                       ndf_patch->plant_potential_ndemand,
                       ndf->potential_N_uptake); //debug
            
            
            sum_plant_nsupply = ns->retransn + soil_nsupply;
            if( ns->retransn<0 || sum_plant_nsupply<0)
                printf("allocate_daily_growth has negative [%d: %d,%d,%d] (%e, %e, %e, %e)\n",
                       stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year,
                       soil_nsupply,
                       ndf_patch->plant_avail_uptake,
                       ns->retransn,
                       sum_plant_nsupply); //debug
            
            //ndf->potential_N_uptake in this case is what plant needs, ns->retransn serves as extra
            //ndf->potential_N_uptake * (ns->retransn/sum_plant_nsupply) is fraction of "ns->retransn" used for growth
            ndf->retransn_to_npool = min(ns->retransn, ndf->potential_N_uptake * (ns->retransn/sum_plant_nsupply) );
            sminn_to_npool = ndf->potential_N_uptake - ndf->retransn_to_npool;
            plant_nalloc = ndf->retransn_to_npool + sminn_to_npool; // should == ndf->potential_N_uptake
            plant_calloc = cs->availc/(1.0+epc.gr_perc); // correct for less growth
            ns->nlimit = 0;
        }// nlimit
    }else{
       // ndf_patch->plant_potential_ndemand == 0 (patch level)
        ndf->retransn_to_npool = 0.0;
        plant_nalloc = 0.0;
        plant_calloc = 0.0;
        ns->nlimit = 0;
        sminn_to_npool = 0.0;
        // --> nlc = 0
        // leafc, frootc and all other c parts = 0
    }//


	plant_nalloc = max(plant_nalloc, 0.0); // could be more than needed!
	plant_calloc = max(plant_calloc, 0.0);
	
	

	/* pnow is the proportion of this day's growth that is displayed now,
	the remainder going into storage for display next year through the
	transfer pools */
	//nlc = plant_calloc * fleaf;

	/* daily C fluxes out of cpool and into new growth or storage */
	/* daily N fluxes out of npool and into new growth or storage */
	
	

    if(plant_calloc>0){
        cdf->cpool_to_leafc            = plant_calloc * fleaf * pnow;
        cdf->cpool_to_leafc_store      = plant_calloc * fleaf * (1.0-pnow);
        ndf->npool_to_leafn            = cdf->cpool_to_leafc / cnl;
        ndf->npool_to_leafn_store      = cdf->cpool_to_leafc_store / cnl;
        
        cdf->cpool_to_frootc           = froot * plant_calloc * pnow;
        cdf->cpool_to_frootc_store     = froot * plant_calloc * (1.0-pnow);
        ndf->npool_to_frootn           = cdf->cpool_to_frootc / cnfr;
        ndf->npool_to_frootn_store     = cdf->cpool_to_frootc_store / cnfr;
        
        if (epc.veg_type == TREE){
            cdf->cpool_to_livestemc        = plant_calloc * flive * fstem * pnow;
            cdf->cpool_to_livestemc_store  = plant_calloc * flive * fstem * (1.0-pnow);
            ndf->npool_to_livestemn        = cdf->cpool_to_livestemc / cnlw;
            ndf->npool_to_livestemn_store  = cdf->cpool_to_livestemc_store / cnlw;
            
            cdf->cpool_to_deadstemc        = plant_calloc * fdead * fstem  * pnow;
            cdf->cpool_to_deadstemc_store  = plant_calloc * fdead * fstem * (1.0-pnow);
            ndf->npool_to_deadstemn        = cdf->cpool_to_deadstemc / cndw;
            ndf->npool_to_deadstemn_store  = cdf->cpool_to_deadstemc_store / cndw;
            
            cdf->cpool_to_livecrootc       = plant_calloc * fcroot * flive  * pnow;
            cdf->cpool_to_livecrootc_store = plant_calloc * fcroot * flive  * (1.0-pnow);
            ndf->npool_to_livecrootn       = cdf->cpool_to_livecrootc / cnlw;
            ndf->npool_to_livecrootn_store = cdf->cpool_to_livecrootc_store / cnlw;
            
            cdf->cpool_to_deadcrootc       = plant_calloc * fcroot * fdead *  pnow;
            cdf->cpool_to_deadcrootc_store = plant_calloc * fcroot * fdead *  (1.0-pnow);
            ndf->npool_to_deadcrootn       = cdf->cpool_to_deadcrootc / cndw;
            ndf->npool_to_deadcrootn_store = cdf->cpool_to_deadcrootc_store / cndw;
            
        }else{
            cdf->cpool_to_livestemc        = 0.0;
            cdf->cpool_to_livestemc_store  = 0.0;
            cdf->cpool_to_deadstemc          = 0.0;
            cdf->cpool_to_deadstemc_store  = 0.0;
            ndf->npool_to_livestemn        = 0.0;
            ndf->npool_to_livestemn_store  = 0.0;
            ndf->npool_to_deadstemn        = 0.0;
            ndf->npool_to_deadstemn_store  = 0.0;
            
            cdf->cpool_to_livecrootc         = 0.0;
            cdf->cpool_to_livecrootc_store = 0.0;
            cdf->cpool_to_deadcrootc         = 0.0;
            cdf->cpool_to_deadcrootc_store = 0.0;
            ndf->npool_to_livecrootn        = 0.0;
            ndf->npool_to_livecrootn_store  = 0.0;
            ndf->npool_to_deadcrootn        = 0.0;
            ndf->npool_to_deadcrootn_store  = 0.0;
        }
    }else{
        // not enough carbon to make a positive availc due to dormin season or deficit in cpool
        cdf->cpool_to_leafc            = 0.0;
        cdf->cpool_to_leafc_store      = 0.0;
        cdf->cpool_to_frootc           = 0.0;
        cdf->cpool_to_frootc_store     = 0.0;
        ndf->npool_to_leafn            = 0.0;
        ndf->npool_to_leafn_store      = 0.0;
        ndf->npool_to_frootn           = 0.0;
        ndf->npool_to_frootn_store     = 0.0;
        
        cdf->cpool_to_livestemc        = 0.0;
        cdf->cpool_to_livestemc_store  = 0.0;
        cdf->cpool_to_deadstemc          = 0.0;
        cdf->cpool_to_deadstemc_store  = 0.0;
        ndf->npool_to_livestemn        = 0.0;
        ndf->npool_to_livestemn_store  = 0.0;
        ndf->npool_to_deadstemn        = 0.0;
        ndf->npool_to_deadstemn_store  = 0.0;
        
        cdf->cpool_to_livecrootc         = 0.0;
        cdf->cpool_to_livecrootc_store = 0.0;
        cdf->cpool_to_deadcrootc         = 0.0;
        cdf->cpool_to_deadcrootc_store = 0.0;
        ndf->npool_to_livecrootn        = 0.0;
        ndf->npool_to_livecrootn_store  = 0.0;
        ndf->npool_to_deadcrootn        = 0.0;
        ndf->npool_to_deadcrootn_store  = 0.0;
    }

	 
	// growth respiration for storages
	if (epc.veg_type == TREE){
		gresp_store = (cdf->cpool_to_leafc_store + cdf->cpool_to_frootc_store
			+ cdf->cpool_to_livestemc_store + cdf->cpool_to_deadstemc_store
			+ cdf->cpool_to_livecrootc_store + cdf->cpool_to_deadcrootc_store)
			* g1;
	}
	else{
		gresp_store = (cdf->cpool_to_leafc_store+cdf->cpool_to_frootc_store) * g1;
	}
	cdf->cpool_to_gresp_store = gresp_store;
    // Laurence, Nov 14, 2018
    // I see this resp. as a cost to build store.
    // it should be an instantaneous cost rather than deposit and do it later.
    // so I made this a) "cdf->cpool_to_gresp_store" substract to cs->cpool in "update_C_stratum_daily"
    // b) keep "cs->gresp_store" zero the whole time @ "update_C_stratum_daily"
    // c) take out the "cs->gresp_store" and "cs->gresp_transfer" process in "allocate_annual_growth"
    
	/*---------------------------------------------------------------------------	*/
	/*	create a maximum lai							*/
	/*---------------------------------------------------------------------------	*/
    excess_lai = 0.0;
    double MAX_LAI = stratum[0].local_max_lai; //epc.max_lai;
    double MAX_LAI_c = MAX_LAI/epc.proj_sla;
    
    // nned to fix below; it prevent allocation to leaf_store for next year if current year LAI reaches its max.
    excess_lai = (cs->leafc + cs->leafc_transfer + cdf->cpool_to_leafc) * epc.proj_sla - MAX_LAI; 
    if ( excess_lai > ZERO){

        // cs->leafc_transfer + cdf->cpool_to_leafc need to go away
        excess_allocation_to_leaf = cs->leafc + cs->leafc_transfer + cdf->cpool_to_leafc - MAX_LAI_c; // if exiting > max_LAI
        
        double tmpRatio = min(1.0, excess_allocation_to_leaf / (cs->leafc_transfer + cdf->cpool_to_leafc));
        cs->cpool += cs->leafc_transfer * tmpRatio;
        ns->npool += ns->leafn_transfer * tmpRatio;
        cdf->cpool_to_leafc_store += cdf->cpool_to_leafc * tmpRatio;
        ndf->npool_to_leafn_store += ndf->npool_to_leafn * tmpRatio;
        
        tmpRatio = max(0.0, 1.0 - tmpRatio);
        cs->leafc_transfer *= tmpRatio;
        ns->leafn_transfer *= tmpRatio;
        cdf->cpool_to_leafc *= tmpRatio;
        ndf->npool_to_leafn *= tmpRatio;
        
    }//excess_lai
    

    
    // final mass balance
    ndf->actual_N_uptake  =
        ndf->npool_to_leafn +
        ndf->npool_to_leafn_store+
        ndf->npool_to_frootn +
        ndf->npool_to_frootn_store +
        ndf->npool_to_livestemn  +
        ndf->npool_to_livestemn_store  +
        ndf->npool_to_deadstemn         +
        ndf->npool_to_deadstemn_store  +
        ndf->npool_to_livecrootn         +
        ndf->npool_to_livecrootn_store +
        ndf->npool_to_deadcrootn         +
        ndf->npool_to_deadcrootn_store;
    
    totalc_used =
        cdf->cpool_to_leafc +
        cdf->cpool_to_leafc_store+
        cdf->cpool_to_frootc +
        cdf->cpool_to_frootc_store +
        cdf->cpool_to_livestemc  +
        cdf->cpool_to_livestemc_store  +
        cdf->cpool_to_deadstemc         +
        cdf->cpool_to_deadstemc_store  +
        cdf->cpool_to_livecrootc         +
        cdf->cpool_to_livecrootc_store +
        cdf->cpool_to_deadcrootc         +
        cdf->cpool_to_deadcrootc_store;
    
    if( ndf->actual_N_uptake!=ndf->actual_N_uptake || sminn_to_npool!=sminn_to_npool ){
        printf("allocate_daily_growth [%d: %d,%d,%d]: (%e{%e,%e}[%e,%e],%e[%e],%d,%e,%e>=%e,%e)\n",
               stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year,
               plant_nalloc, ns->retransn, ndf->potential_N_uptake,
               ndf->retransn_to_npool,sminn_to_npool, //nan
               plant_calloc,cs->availc,
               nlimit,
               patch[0].soil_ns.fract_potential_uptake,
               ndf_patch->plant_potential_ndemand,
               ndf_patch->plant_avail_uptake,
               ndf->actual_N_uptake);
    }//debug
    
    
    // fabs(ndf->actual_N_uptake-plant_nalloc) > 1e-8
    // fabs(totalc_used-plant_calloc*(fleaf+froot+fwood)) > 1e-8
    plant_calloc *= (fleaf+froot+fwood);
    if( fabs(ndf->actual_N_uptake-plant_nalloc)>ZERO || fabs(totalc_used-plant_calloc)>ZERO || sminn_to_npool>patch[0].soil_ns.fract_potential_uptake * ndf->potential_N_uptake){
        // actual_N_uptake > plant_nalloc problem!
            
        printf("allocation_daily[%d: %d,%d,%d]::%d, C{%e = %e = %e, %e %e %e} N{(%f)%e = %e = %e(%e), %e, %e -> %e/%e, %e/%e, %e/%e, %e/%e, %e/%e, %e/%e, %e/%e, %e/%e, %e/%e, %e/%e, %e/%e, %e/%e} \n",
               stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year,
               nlimit,
               totalc_used, plant_calloc, cs->availc/(1+epc.gr_perc), fleaf, froot, fwood,
               mean_cn, ndf->actual_N_uptake, plant_nalloc, ndf->potential_N_uptake, patch[0].soil_ns.fract_potential_uptake,
               sminn_to_npool, ndf->retransn_to_npool,
               //-----
                cdf->cpool_to_leafc,ndf->npool_to_leafn,
                cdf->cpool_to_leafc_store,ndf->npool_to_leafn_store,
                cdf->cpool_to_frootc,ndf->npool_to_frootn,
                cdf->cpool_to_frootc_store,ndf->npool_to_frootn_store, //<----
                cdf->cpool_to_livestemc,ndf->npool_to_livestemn,
                cdf->cpool_to_livestemc_store,ndf->npool_to_livestemn_store,
                cdf->cpool_to_deadstemc,ndf->npool_to_deadstemn,
                cdf->cpool_to_deadstemc_store,ndf->npool_to_deadstemn_store,
                cdf->cpool_to_livecrootc,ndf->npool_to_livecrootn,
                cdf->cpool_to_livecrootc_store,ndf->npool_to_livecrootn_store,
                cdf->cpool_to_deadcrootc,ndf->npool_to_deadcrootn,
                cdf->cpool_to_deadcrootc_store,ndf->npool_to_deadcrootn_store
               );
     
        // why "ndf->actual_N_uptake" is so large, and yet make sense value?
    }//debug

    excess_c = max(cs->availc - (totalc_used*(1+epc.gr_perc)),0.0); // adjust to cdf->psn_to_cpool
    cdf->psn_to_cpool -= excess_c; // if not enough N, it will reduce fraq_psn result.
//    if(plant_nalloc>0 && plant_calloc>0){
//        ///<<----------- bad move?
//        ndf->retransn_to_npool *= ndf->actual_N_uptake/plant_nalloc;
//        sminn_to_npool *= ndf->actual_N_uptake/plant_nalloc;
//        plant_nalloc = ndf->retransn_to_npool + sminn_to_npool;
//    }
    
    
    
    ndf->sminn_to_npool = sminn_to_npool;
    ndf_patch->sminn_to_npool += sminn_to_npool * cover_fraction;
    
	return(!ok);
} /* end daily_allocation.c */

