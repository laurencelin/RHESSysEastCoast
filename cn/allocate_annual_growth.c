/*--------------------------------------------------------------*/
/* 								*/
/*								*/
/*	allocate_annual_growth				*/
/*								*/
/*	NAME							*/
/*		allocate_annual_growth			*/
/*								*/
/*	SYNOPSIS						*/
/*	int allocate_annual_growth( int , int, double, double,	*/
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

int allocate_annual_growth(				int id,
							int default_ID,
							int vmort_flag,
							double cover_fraction,
							double cpool_mort_fract,
						   struct epvar_struct *epv,
						   struct cdayflux_struct *cdf,
						   struct cstate_struct *cs,
						   struct ndayflux_struct *ndf,
						   struct nstate_struct *ns,
						   struct cdayflux_patch_struct *cdf_patch,
						   struct ndayflux_patch_struct *ndf_patch,
						   struct litter_c_object *cs_litr,
						   struct litter_n_object *ns_litr,
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

    double totalLeafCarbon, totalFrootCarbon, totalLiveStemCarbon, totalLiveCrootCarbon;
	double rem_excess_carbon, excess_nitrogen;
	double total_store, total_livebiomass;
    double npoolwithdraw, fillN, fillC;
    
    double f1 = epc.alloc_frootc_leafc; //2
    double f2 = epc.alloc_crootc_stemc; //0
    double f3 = epc.alloc_stemc_leafc; //0
    double f4 = epc.alloc_livewoodc_woodc; //0
    double fleaf =  1.0/(1.0+f1+f3+f2*f3); //1/3
    double froot = fleaf*f1; //2/3
    double fwood = fleaf*f3*(1.0+f2);//0
    
    double flive = epc.alloc_livewoodc_woodc; //0
    double fdead = (1-flive);
	double fcroot = f2 / (1.0+f2); // 0
	
	double cnl = epc.leaf_cn;
	double cnfr = epc.froot_cn;
	double cnlw = epc.livewood_cn;
	double cndw = epc.deadwood_cn;

    
    
    
    
    //annual plant report for its growth period performance
    double what = 1.0/(1.0*stratum->gDayCount);
    double evaluate = stratum[0].cs.cpool + (stratum->gwPSN - stratum->gwMResp)*what*(1.0*stratum->gDayCount) - stratum->gwMResp*what*(365.0-1.0*stratum->gDayCount);
    if((epc.veg_type == TREE && evaluate < -0.09) || (epc.veg_type == GRASS && evaluate < -0.04)){
        printf("report,%d,%d,%d,%d,%d,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",
               patch[0].ID, stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year,
               stratum->gDayCount,
               stratum->nFactor *what, // 1 is good
               stratum->wFactor *what, // 1 is good
               stratum->lFactor *what, // large is good
               stratum->gFactor *what, // 1 is good
               patch[0].constraintWaterTableTopDepth,
               stratum[0].cover_fraction,
               stratum[0].cs.cpool,
               stratum->gwPSN *what,   // flux
               stratum->gwMResp *what, // flux
               stratum->gwAPAR *what, // 1 is good
               stratum->gwLWP *what, // 1 is good
               stratum->gwVPD *what, // 1 is good
               patch[0].sat_deficit_z,          //<<---------------- instant
               patch[0].basementSideAdjustWTZ,  //<<---------------- instant
               patch[0].basementSideAdjustH2O,  //<<---------------- instant
               patch[0].Ksat_vertical, //(stratum[0].cs.leafc + stratum[0].cs.leafc_store + stratum[0].cs.leafc_transfer),
               (stratum[0].cs.frootc + stratum[0].cs.frootc_store + stratum[0].cs.frootc_transfer)
               );
    }//print report
    stratum->gDayCount=0;
    stratum->nFactor=0.0; // tracking @ allocate_daily_growth  <<--------- need to check
    stratum->wFactor=0.0; // tracking @ patch_daily_F          <<--------- need to check
    stratum->lFactor=0.0; // tracking @ canopy_stratum_daily_F <<--------- need to check
    stratum->gFactor=0.0; // tracking @ canopy_stratum_daily_F <<--------- actual gl vs potential gl (what's difference?)
    stratum->gwPSN=0.0; // tracking @ canopy_stratum_daily_F
    stratum->gwMResp=0.0; // tracking @ canopy_stratum_daily_F
    stratum->gwAPAR=0.0; // tracking @ canopy_stratum_daily_F
    stratum->gwLWP=0.0; // tracking @ canopy_stratum_daily_F
    stratum->gwVPD=0.0; // tracking @ canopy_stratum_daily_F
    // patch vegid day month year gday nfactor wfactor lfactor gfactor basement cover cpool gwPSN gwMresp gwAPAR gwLWP gwVPD wtz adjustWTZ adjustH2O leafc frootc
    
    
    
 
    /* for deciduous system, force leafc and frootc to exactly 0.0 on the last day */
    if (epc.phenology_type != EVERGREEN){
        if (ns->leafn < 1e-10)  {
            ns->leafn = 0.0;
            cs->leafc = 0.0;
        }
        if (ns->frootn < 1e-10) {
            ns->frootn = 0.0;
            cs->frootc = 0.0;
        }
    }
    
    totalLeafCarbon = cs->leafc + cs->leafc_store + cs->leafc_transfer;
    totalFrootCarbon = cs->frootc + cs->frootc_store + cs->frootc_transfer;
    if (epc.veg_type == TREE) {
        totalLiveStemCarbon = cs->live_stemc + cs->livestemc_store + cs->livestemc_transfer;
        totalLiveCrootCarbon = cs->live_crootc + cs->livecrootc_store + cs->livecrootc_transfer;
    }else{
        totalLiveStemCarbon = 0.0;
        totalLiveCrootCarbon = 0.0;
    }
    total_livebiomass = totalLeafCarbon + totalFrootCarbon + totalLiveStemCarbon + totalLiveCrootCarbon;
    
    double MAX_LAI = stratum->local_max_lai; //epc.max_lai;
    double allometric_leafRatio = fleaf/(fleaf+froot+fwood*flive);
    double allometric_frootRatio = froot/(fleaf+froot+fwood*flive);
    double allometric_woodRatio = fwood*flive/(fleaf+froot+fwood*flive); // what is this for GRASS?
    double minLeafCarbon = total_livebiomass * allometric_leafRatio;
    double minFrootCarbon = total_livebiomass * allometric_frootRatio;
    double minWoodCarbon = total_livebiomass * allometric_woodRatio;
    double excess_carbon = totalLeafCarbon - max(MAX_LAI/epc.proj_sla, minLeafCarbon);//leafc
    double miss_Leafcarbon = totalLeafCarbon - minLeafCarbon;
    double miss_Frootcarbon = totalFrootCarbon - minFrootCarbon;
    double miss_Woodcarbon = totalLiveStemCarbon + totalLiveCrootCarbon - minWoodCarbon;
    
    if( excess_carbon > 0){
        //too much leaf carbon
        rem_excess_carbon = excess_carbon;
        if (epc.veg_type == TREE) {
//            printf("allocate_annual_growth excess_lai tree[%d,%d: %d,%d,%d]: (%e:leaf[%e,%e,%e])\n",
//                   patch[0].ID,stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year,
//                   rem_excess_carbon, cs->leafc,cs->leafc_store,cs->leafc_transfer);
            /* remove excess carbon from storage, transfer and then leaf carbon until gone */
            if (cs->leafc_store > excess_carbon) {
                // reduce store
                cs->leafc_store -= rem_excess_carbon;
                ns->leafn_store -= rem_excess_carbon / cnl;
                
                cs->deadstemc_store += fdead*excess_carbon;
                cs->livestemc_store += flive*excess_carbon;
                ns->deadstemn_store += fdead*excess_carbon / cndw;
                ns->livestemn_store += flive*excess_carbon / cnlw;
                excess_nitrogen = excess_carbon/cnl - fdead*excess_carbon/cndw - flive*excess_carbon/cnlw;
                ns->npool += excess_nitrogen;
            } else {
                // reduce transfer
                rem_excess_carbon -= cs->leafc_store;
                cs->leafc_store = 0.0;
                ns->leafn_store = 0.0;
                if (cs->leafc_transfer > rem_excess_carbon) {
                    cs->leafc_transfer -= rem_excess_carbon;
                    ns->leafn_transfer -= rem_excess_carbon / cnl;
                    
                    cs->deadstemc_store += fdead*excess_carbon;
                    cs->livestemc_store += flive*excess_carbon;
                    ns->deadstemn_store += fdead*excess_carbon / cndw;
                    ns->livestemn_store += flive*excess_carbon / cnlw;
                    excess_nitrogen = excess_carbon/cnl - fdead*excess_carbon/cndw - flive*excess_carbon/cnlw;
                    ns->npool += excess_nitrogen;
                } else {
                    // make extra leaf fall
                    rem_excess_carbon -= cs->leafc_transfer;
                    cs->leafc_transfer = 0.0;
                    ns->leafn_transfer = 0.0;
                    cs->leafc -= rem_excess_carbon;
                    ns->leafn -= rem_excess_carbon / cnl;
                    
                    //litter pool
                    
                }
            }//if
            //<<------------------------------------------------------- pollen
            // pollen/seedling CN ratio is about 40
            // (https://nph.onlinelibrary.wiley.com/doi/full/10.1111/j.1469-8137.2007.02176.x)
            // (https://link.springer.com/article/10.1007/s11738-015-1965-x)
        }else{
            // not tree below
//            printf("allocate_annual_growth excess_lai non-tree[%d,%d: %d,%d,%d]: (%e)\n",
//                   patch[0].ID,stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year,
//                   rem_excess_carbon);
            
            if (cs->leafc_store > excess_carbon) {
                // reduce store
                cs->leafc_store -= rem_excess_carbon;
                ns->leafn_store -= rem_excess_carbon / cnl;
                
                cs->frootc_store += excess_carbon;
                ns->frootn_store += excess_carbon / cnfr;
                excess_nitrogen = excess_carbon/cnl - excess_carbon/cnfr;
                ns->npool += excess_nitrogen;
            }
            else {
                rem_excess_carbon -= cs->leafc_store;
                cs->leafc_store = 0.0;
                ns->leafn_store = 0.0;
                if (cs->leafc_transfer > rem_excess_carbon) {
                    // reduce transfer
                    cs->leafc_transfer -= rem_excess_carbon;
                    ns->leafn_transfer -= rem_excess_carbon / cnl;
                    
                    cs->frootc_store += excess_carbon;
                    ns->frootn_store += excess_carbon / cnfr;
                    excess_nitrogen = excess_carbon/cnl - excess_carbon/cnfr;
                    ns->npool += excess_nitrogen;
                } else {
                    // make extra leaf fall
                    rem_excess_carbon -= cs->leafc_transfer;
                    cs->leafc_transfer = 0.0;
                    ns->leafn_transfer = 0.0;
                    cs->leafc -= rem_excess_carbon;
                    ns->leafn -= rem_excess_carbon / cnl;
                    
                    // litter pool
                }
            }//if
        }// tree
    }else{
        // not enough for some parts
        if(miss_Leafcarbon<0 && miss_Leafcarbon < miss_Frootcarbon && miss_Leafcarbon < miss_Woodcarbon){
            // need to fix leafc
            if(miss_Frootcarbon>0){
                double withdraw = min(miss_Frootcarbon, cs->frootc_store);
                if(withdraw){
                    npoolwithdraw = min(-miss_Leafcarbon/cnl - withdraw/cnfr, ns->npool);
                    fillN = withdraw/cnfr + npoolwithdraw;
                    fillC = fillN*cnl;
                    
                    miss_Leafcarbon += fillC;
                    cs->leafc_store += fillC;
                    ns->leafn_store += fillN;
                    ns->npool -= npoolwithdraw; if(ns->npool<0) ns->npool=0.0;
                    
                    cs->frootc_store -= fillC;
                    ns->frootn_store -= fillC / cnfr;
                    if(cs->frootc_store<0){ cs->frootc_store=0.0; ns->frootn_store=0.0;}
                    //cs->frootc_transfer;
                    //ns->frootn_transfer;
                }
            }//froot
            if(miss_Woodcarbon>0){
                double withdraw = min(miss_Woodcarbon, cs->livecrootc_store+cs->livestemc_store);
                if(withdraw>0){
                    npoolwithdraw = min(-miss_Leafcarbon/cnl - withdraw/cnlw, ns->npool);
                    fillN = withdraw/cnlw + npoolwithdraw;
                    fillC = fillN*cnl;
                    
                    miss_Leafcarbon += fillC;
                    cs->leafc_store += fillC;
                    ns->leafn_store += fillN;
                    ns->npool -= npoolwithdraw; if(ns->npool<0) ns->npool=0.0;
                    
                    double hold = cs->livecrootc_store/(cs->livecrootc_store+cs->livestemc_store);
                    cs->livecrootc_store -= fillC*hold;
                    ns->livecrootn_store -= fillC*hold/cnlw;
                    if(cs->livecrootc_store<0){ cs->livecrootc_store=0.0; ns->livecrootn_store=0.0;}
                    
                    cs->livestemc_store -= fillC*(1.0-hold);
                    ns->livestemn_store -= fillC*(1.0-hold)/cnlw;
                    if(cs->livestemc_store<0){ cs->livestemc_store=0.0; ns->livestemn_store=0.0;}
                }
            }
        }//miss_Frootcarbon
//        }else if( miss_Frootcarbon<0 && miss_Frootcarbon < miss_Leafcarbon && miss_Frootcarbon < miss_Woodcarbon){
//            // need to fix froot
//        }else if(miss_Woodcarbon<0 && miss_Woodcarbon < miss_Leafcarbon && miss_Woodcarbon < miss_Frootcarbon){
//            // need to fix wood
//        }//if
        
    }//if
    
    
    
    epv->prev_leafcalloc = cs->leafc + cs->leafc_transfer;
    cdf->leafc_store_to_leafc_transfer = cs->leafc_store * epc.storage_transfer_prop; // first time
    ndf->leafn_store_to_leafn_transfer = ns->leafn_store * epc.storage_transfer_prop;
    cs->leafc_transfer    += cdf->leafc_store_to_leafc_transfer;
    ns->leafn_transfer    += ndf->leafn_store_to_leafn_transfer;
        cs->cpool += cs->leafc_store * (1.0-epc.storage_transfer_prop);
        ns->npool += ns->leafn_store * (1.0-epc.storage_transfer_prop);
        cs->leafc_store = 0.0;
        ns->leafn_store = 0.0;
    cdf->frootc_store_to_frootc_transfer = cs->frootc_store* epc.storage_transfer_prop;
    ndf->frootn_store_to_frootn_transfer = ns->frootn_store* epc.storage_transfer_prop;
    cs->frootc_transfer   += cdf->frootc_store_to_frootc_transfer;
    ns->frootn_transfer   += ndf->frootn_store_to_frootn_transfer;
        cs->cpool += cs->frootc_store * (1.0-epc.storage_transfer_prop);
        ns->npool += ns->frootn_store * (1.0-epc.storage_transfer_prop);
        cs->frootc_store = 0.0;
        ns->frootn_store = 0.0;
    if (epc.veg_type == TREE){
        cdf->livestemc_store_to_livestemc_transfer = cs->livestemc_store * epc.storage_transfer_prop;
        ndf->livestemn_store_to_livestemn_transfer = ns->livestemn_store * epc.storage_transfer_prop;
        cs->livestemc_transfer  += cdf->livestemc_store_to_livestemc_transfer;
        ns->livestemn_transfer  += ndf->livestemn_store_to_livestemn_transfer;
            cs->cpool += cs->livestemc_store * (1.0-epc.storage_transfer_prop);
            ns->npool += ns->livestemn_store * (1.0-epc.storage_transfer_prop);
            cs->livestemc_store = 0.0;
            ns->livestemn_store = 0.0;
        cdf->deadstemc_store_to_deadstemc_transfer = cs->deadstemc_store * epc.storage_transfer_prop;
        ndf->deadstemn_store_to_deadstemn_transfer = ns->deadstemn_store * epc.storage_transfer_prop;
        cs->deadstemc_transfer  += cdf->deadstemc_store_to_deadstemc_transfer;
        ns->deadstemn_transfer  += ndf->deadstemn_store_to_deadstemn_transfer;
            cs->cpool += cs->deadstemc_store * (1.0-epc.storage_transfer_prop);
            ns->npool += ns->deadstemn_store * (1.0-epc.storage_transfer_prop);
            cs->deadstemc_store = 0.0;
            ns->deadstemn_store = 0.0;
        cdf->livecrootc_store_to_livecrootc_transfer = cs->livecrootc_store * epc.storage_transfer_prop;
        ndf->livecrootn_store_to_livecrootn_transfer = ns->livecrootn_store * epc.storage_transfer_prop;
        cs->livecrootc_transfer += cdf->livecrootc_store_to_livecrootc_transfer;
        ns->livecrootn_transfer += ndf->livecrootn_store_to_livecrootn_transfer;
            cs->cpool += cs->livecrootc_store * (1.0-epc.storage_transfer_prop);
            ns->npool += ns->livecrootn_store * (1.0-epc.storage_transfer_prop);
            cs->livecrootc_store = 0.0;
            ns->livecrootn_store = 0.0;
        cdf->deadcrootc_store_to_deadcrootc_transfer = cs->deadcrootc_store * epc.storage_transfer_prop;
        ndf->deadcrootn_store_to_deadcrootn_transfer = ns->deadcrootn_store * epc.storage_transfer_prop;
        cs->deadcrootc_transfer += cdf->deadcrootc_store_to_deadcrootc_transfer;
        ns->deadcrootn_transfer += ndf->deadcrootn_store_to_deadcrootn_transfer;
            cs->cpool += cs->deadcrootc_store * (1.0-epc.storage_transfer_prop);
            ns->npool += ns->deadcrootn_store * (1.0-epc.storage_transfer_prop);
            cs->deadcrootc_store = 0.0;
            ns->deadcrootn_store = 0.0;
    }
    totalLeafCarbon = cs->leafc + cs->leafc_store + cs->leafc_transfer;
    

    //    cs->gresp_transfer    += cdf->gresp_store_to_gresp_transfer;//<<-------------------- gresp from store; but does not make sense anyway
    //    if(cs->gresp_transfer <= 0){
    //        //printf("allocate_annual_growth [%d,%d],%e\n",id,default_ID,cs->gresp_transfer);
    //        //cs->gresp_transfer = 0.0;
    //    }//debug
    
	/*--------------------------------------------------------------*/
	/* finally if there is really nothing restart with small amount of growth   */
	/*	we allow only a certain amount of resprouting based on 	*/
	/*	a stratum default file parameterization 		*/
	/*--------------------------------------------------------------*/
//  if((epc.veg_type == TREE && evaluate < -0.09) || (epc.veg_type == GRASS && evaluate < -0.04)) --- same conditions for "report"
    if ((epc.veg_type == TREE && evaluate < -0.09) || (epc.veg_type == GRASS && evaluate < -0.04)){  //totalLeafCarbon < epc.min_leaf_carbon //epc.min_leaf_carbon; 0.1*minLeafCarbon
        if (cs->num_resprout < (epc.veg_type == GRASS? 100 : epc.max_years_resprout) ) {
            // SLB whole basin grass 0.307539*1.5=0.4613085 LAI; tree 0.406821*4.5=1.830694 LAI
            // SLB subbain grass 0.257967*1.5=0.3869505 basin LAI; tree 0.368434*4.5=1.657953 basin LAI
//            printf("Resprouting stratum [%d,%d: %d,%d,%d]{%e,%e} (%e,%e[%e,%e,%e]:%e:%e:%e)->(%e,%e,%e);[%e,%e,%e,%e][%e,%e,%e,%e]\n",
//                   patch[0].ID,stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year,
//                   patch[0].Ksat_vertical,stratum[0].cover_fraction,
//                   0.1*minLeafCarbon, totalLeafCarbon, cs->leafc, cs->leafc_store, cs->leafc_transfer,
//                   totalFrootCarbon, totalLiveStemCarbon, totalLiveCrootCarbon,
//                   miss_Leafcarbon, miss_Frootcarbon, miss_Woodcarbon,
//                   patch[0].sat_deficit_z, patch[0].rootzone.S, patch[0].soil_ns.nitrate, patch[0].soil_ns.sminn,
//                   patch[0].soil_cs.soil1c,patch[0].soil_cs.soil2c,patch[0].soil_cs.soil3c,patch[0].soil_cs.soil4c
//                   //stratum->mult_conductance.APAR, stratum->mult_conductance.LWP, stratum->mult_conductance.vpd
//                   );
            cs->num_resprout += 1;
            cs->age = 0;
            cs->cpool = 0.0;
            ns->npool = 0.0;
            cs->leafc_store = epc.resprout_leaf_carbon;
            cs->frootc_store = cs->leafc_store * epc.alloc_frootc_leafc;
            cdf->leafc_store_to_leafc_transfer = cs->leafc_store;
            cdf->frootc_store_to_frootc_transfer = cs->frootc_store;
            ns->leafn_store = cs->leafc_store / epc.leaf_cn;
            ns->frootn_store = cs->frootc_store / epc.froot_cn;
            ndf->leafn_store_to_leafn_transfer = ns->leafn_store;
            ndf->frootn_store_to_frootn_transfer = ns->frootn_store;
            cs->leafc_transfer = 0.0;
            cs->leafc = 0.0;
            cs->frootc_transfer = 0.0;
            ns->leafn_transfer = 0.0;
            ns->frootn_transfer = 0.0;
            ns->leafn = 0.0;
            epv->prev_leafcalloc = epc.resprout_leaf_carbon;

            if (epc.veg_type == TREE) {
                cs->livecrootc_store = epc.resprout_leaf_carbon * epc.alloc_stemc_leafc;
                ns->livecrootn_store = cs->livecrootc_store / epc.livewood_cn;
                cdf->livecrootc_store_to_livecrootc_transfer = cs->livecrootc_store;
                ndf->livecrootn_store_to_livecrootn_transfer = ns->livecrootn_store;

                cs->livestemc_store = epc.resprout_leaf_carbon * epc.alloc_stemc_leafc;
                ns->livestemn_store = cs->livestemc_store / epc.livewood_cn;
                cdf->livestemc_store_to_livestemc_transfer = cs->livestemc_store;
                ndf->livestemn_store_to_livestemn_transfer = ns->livestemn_store;

                cs->live_stemc = 0.0;
                cs->dead_stemc = 0.0;
                cs->live_crootc = 0.0;
                cs->dead_crootc = 0.0;
                ns->live_stemn = 0.0;
                ns->dead_stemn = 0.0;
                ns->live_crootn = 0.0;
                ns->dead_crootn = 0.0;

                cs->deadstemc_store =  0.0;
                cs->deadcrootc_store =  0.0;
                ns->deadstemn_store =  0.0;
                ns->deadcrootn_store =  0.0;

                cs->livestemc_transfer =  0.0;
                cs->deadstemc_transfer =  0.0;
                cs->livecrootc_transfer =  0.0;
                cs->deadcrootc_transfer =  0.0;
                ns->livestemn_transfer =  0.0;
                ns->deadstemn_transfer =  0.0;
                ns->livecrootn_transfer =  0.0;
                ns->deadcrootn_transfer =  0.0;

                cdf->livestemc_store_to_livestemc_transfer = 0.0;
                cdf->deadstemc_store_to_deadstemc_transfer = 0.0;
                cdf->livecrootc_store_to_livecrootc_transfer = 0.0;
                cdf->deadcrootc_store_to_deadcrootc_transfer = 0.0;

                ndf->livestemn_store_to_livestemn_transfer = 0.0;
                ndf->deadstemn_store_to_deadstemn_transfer = 0.0;
                ndf->livecrootn_store_to_livecrootn_transfer = 0.0;
                ndf->deadcrootn_store_to_deadcrootn_transfer = 0.0;
                
            } /* end if TREE */
        } /* end if resprout */
    } else {
        cs->num_resprout = max(cs->num_resprout-1,0);
        cs->age += 1;
    }// min leafc

    

	

	return (!ok);
} /* end allocate_annual_growth */
