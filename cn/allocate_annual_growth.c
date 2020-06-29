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

    double existingLeafCarbon, existingFrootCarbon, existingLiveStemCarbon, existingLiveCrootCarbon;
	double excess_nitrogen;
	double total_store, existingLivebiomass;
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
    // gwPSN = GPP;
    double what = 1.0/(1.0*stratum->gDayCount);
    double evaluate = stratum[0].cs.cpool + (stratum->gwPSN/(1.0+stratum[0].defaults[0][0].epc.gr_perc) - stratum->gwMResp);
    evaluate -= epc.phenology_type != EVERGREEN? (stratum->gwMResp - stratum->gwMRespLeaf)*what*(365.0-1.0*stratum->gDayCount) : stratum->gwMResp*what*(365.0-1.0*stratum->gDayCount);
    
    
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
    }// count for slight negative in leaffall process
    
    
    existingLeafCarbon = cs->leafc + cs->leafc_store + cs->leafc_transfer;
    existingFrootCarbon = cs->frootc + cs->frootc_store + cs->frootc_transfer;
    if (epc.veg_type == TREE) {
        existingLiveStemCarbon = cs->live_stemc + cs->livestemc_store + cs->livestemc_transfer;
        existingLiveCrootCarbon = cs->live_crootc + cs->livecrootc_store + cs->livecrootc_transfer;
    }else{
        existingLiveStemCarbon = 0.0;
        existingLiveCrootCarbon = 0.0;
    }
    existingLivebiomass = existingLeafCarbon + existingFrootCarbon + existingLiveStemCarbon + existingLiveCrootCarbon;
    
    double MAX_LAI = epc.max_lai; // stratum->local_max_lai; //epc.max_lai;
    double allometric_leafRatio = fleaf/(fleaf+froot+fwood*flive);
    double allometric_frootRatio = froot/(fleaf+froot+fwood*flive);
    double allometric_livewoodRatio = fwood*flive/(fleaf+froot+fwood*flive);
    double minLeafCarbon = existingLivebiomass * allometric_leafRatio;
    double minFrootCarbon = existingLivebiomass * allometric_frootRatio;
    double minWoodCarbon = existingLivebiomass * allometric_livewoodRatio;
    
    // this function is call at the last day of leaf fall, right?
    //stratum[0].cdf.leaf_day_mr reset everyday
    double totalResp = epc.veg_type == TREE? 0.002739726/(ns->live_stemn + ns->live_crootn + ns->frootn + (epc.phenology_type != EVERGREEN? 0.0 : ns->leafn)) : 0.002739726/(ns->frootn + ns->leafn);// conver annual to daily
    double potReduceLeafCarbon = evaluate>0? 0.0 : -evaluate*(epc.phenology_type != EVERGREEN? 0.0 : ns->leafn)*totalResp;
    double potReduceCrootCarbon = evaluate>0? 0.0 : -evaluate*ns->live_crootn*totalResp;
    double potReduceStemCarbon = evaluate>0? 0.0 : -evaluate*ns->live_stemn*totalResp;
    double potReduceFrootCarbon = evaluate>0? 0.0 : -evaluate*ns->frootn*totalResp;
    potReduceLeafCarbon = max(existingLeafCarbon - max(MAX_LAI/epc.proj_sla, minLeafCarbon), potReduceLeafCarbon);//issue?
    //for debug
//    if(potReduceLeafCarbon>0 && epc.phenology_type==DECID){
//        printf("annual allocation %d,%d,%f,%f,%f,%f,%f\n",patch[0].ID, stratum->defaults[0][0].ID,
//               existingLeafCarbon, minLeafCarbon, MAX_LAI, epc.proj_sla, potReduceLeafCarbon);
//    }//end of if
    
    double rem_excess_LeafCarbon;
    double rem_excess_StemCarbon=0.0;
    double exceedN=0.0;
    double exceedC=0.0;
    double tmpRatio, tmpMultilper;
    
    // --------- stem
    if(epc.veg_type == TREE && potReduceStemCarbon>0){
        if(cs->livestemc_transfer>0 && cs->livestemc_transfer > potReduceStemCarbon){
            tmpRatio = potReduceStemCarbon/cs->livestemc_transfer;
            exceedC += potReduceStemCarbon;
            exceedN += ns->livestemn_transfer * tmpRatio;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->livestemc_transfer *= tmpRatio;
            ns->livestemn_transfer *= tmpRatio;
            potReduceStemCarbon=0.0;
        }else{
            exceedC += cs->livestemc_transfer;
            exceedN += ns->livestemn_transfer;
            potReduceStemCarbon -= cs->livestemc_transfer;
            cs->livestemc_transfer = 0.0;
            ns->livestemn_transfer = 0.0;
        }// end of if
    }//end of if
    if(epc.veg_type == TREE && potReduceStemCarbon>0){
        if(cs->livestemc_store >0 && cs->livestemc_store > potReduceStemCarbon){
            tmpRatio = potReduceStemCarbon/cs->livestemc_store;
            exceedC += potReduceStemCarbon;
            exceedN += ns->livestemn_transfer * tmpRatio;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->livestemc_store *= tmpRatio;
            ns->livestemn_store *= tmpRatio;
            potReduceStemCarbon=0.0;
        }else{
            exceedC += cs->livestemc_store;
            exceedN += ns->livestemn_store;
            potReduceStemCarbon -= cs->livestemc_store;
            cs->livestemc_store = 0.0;
            ns->livestemn_store = 0.0;
        }// end of if
    }// end of if
    if(epc.veg_type == TREE && potReduceStemCarbon>0){
        if(cs->live_stemc >0 && cs->live_stemc*0.2 > potReduceStemCarbon){
            tmpRatio = potReduceStemCarbon/cs->live_stemc;
            exceedN += ns->live_stemn * tmpRatio - potReduceStemCarbon/epc.deadwood_cn;
            cs->dead_stemc += potReduceStemCarbon;
            ns->dead_stemn += potReduceStemCarbon/epc.deadwood_cn;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->live_stemc *= tmpRatio;
            ns->live_stemn *= tmpRatio;
            potReduceStemCarbon=0.0;
        }else if(cs->live_stemc >0){
            double tmpC = 0.2*cs->live_stemc;
            tmpRatio = 0.2;
            exceedN += ns->live_stemn * tmpRatio - tmpC/epc.deadwood_cn;
            cs->dead_stemc += tmpC;
            ns->dead_stemn += tmpC/epc.deadwood_cn;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->live_stemc *= tmpRatio;
            ns->live_stemn *= tmpRatio;
            potReduceStemCarbon -= tmpC;
        }// end of if
    }// end of if
    
    // ---- croot
    if(epc.veg_type == TREE && potReduceCrootCarbon>0){
        if(cs->livecrootc_transfer>0 && cs->livecrootc_transfer > potReduceCrootCarbon){
            tmpRatio = potReduceCrootCarbon/cs->livecrootc_transfer;
            exceedC += potReduceCrootCarbon;
            exceedN += ns->livecrootn_transfer * tmpRatio;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->livecrootc_transfer *= tmpRatio;
            ns->livecrootn_transfer *= tmpRatio;
            potReduceCrootCarbon=0.0;
        }else{
            exceedC += cs->livecrootc_transfer;
            exceedN += ns->livecrootn_transfer;
            potReduceCrootCarbon -= cs->livecrootc_transfer;
            cs->livecrootc_transfer = 0.0;
            ns->livecrootn_transfer = 0.0;
        }// end of if
    }// end of if
    if(epc.veg_type == TREE && potReduceCrootCarbon>0){
        if(cs->livecrootc_store >0 && cs->livecrootc_store > potReduceCrootCarbon){
            tmpRatio = potReduceCrootCarbon/cs->livecrootc_store;
            exceedC += potReduceCrootCarbon;
            exceedN += ns->livecrootn_store * tmpRatio;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->livecrootc_store *= tmpRatio;
            ns->livecrootn_store *= tmpRatio;
            potReduceCrootCarbon=0.0;
        }else{
            exceedC += cs->livecrootc_store;
            exceedN += ns->livecrootn_store;
            potReduceCrootCarbon -= cs->livecrootc_store;
            cs->livecrootc_store = 0.0;
            ns->livecrootn_store = 0.0;
        }// end of if
    }// end of if
    if(epc.veg_type == TREE && potReduceCrootCarbon>0){
        if(cs->live_crootc >0 && cs->live_crootc*0.2 > potReduceCrootCarbon){
            tmpRatio = potReduceCrootCarbon/cs->live_crootc;
            exceedN += ns->live_crootn * tmpRatio - potReduceCrootCarbon/epc.deadwood_cn;
            cs->dead_crootc += potReduceCrootCarbon;
            ns->dead_crootn += potReduceCrootCarbon/epc.deadwood_cn;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->live_crootc *= tmpRatio;
            ns->live_crootn *= tmpRatio;
            potReduceCrootCarbon=0.0;
        }else if(cs->live_crootc >0){
            double tmpC = 0.2 * cs->live_crootc;
            tmpRatio = 0.2;
            exceedN += ns->live_crootn * tmpRatio - tmpC/epc.deadwood_cn;
            cs->dead_crootc += tmpC;
            ns->dead_crootn += tmpC/epc.deadwood_cn;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->live_crootc *= tmpRatio;
            ns->live_crootn *= tmpRatio;
            potReduceCrootCarbon -= tmpC;
        }// end of if
    }// end of if
    
    // ---- leafc
    double maxReducePropLeaf = epc.phenology_type == EVERGREEN? 0.4 : 0.2;
    if(epc.veg_type != NON_VEG && potReduceLeafCarbon>0){
        if(cs->leafc_transfer>0 && cs->leafc_transfer*maxReducePropLeaf > potReduceLeafCarbon){
            tmpRatio = potReduceLeafCarbon/cs->leafc_transfer;
            exceedC += potReduceLeafCarbon;
            exceedN += ns->leafn_transfer * tmpRatio;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->leafc_transfer *= tmpRatio;
            ns->leafn_transfer *= tmpRatio;
            potReduceLeafCarbon=0.0;
        }else{
            exceedC += cs->leafc_transfer*maxReducePropLeaf;
            exceedN += ns->leafn_transfer*maxReducePropLeaf;
            potReduceLeafCarbon -= cs->leafc_transfer*maxReducePropLeaf;
            cs->leafc_transfer *= 1.0-maxReducePropLeaf;
            ns->leafn_transfer *= 1.0-maxReducePropLeaf;
        }// end of if
    }// end of if
    if(epc.veg_type != NON_VEG && potReduceLeafCarbon>0){
        if(cs->leafc_store >0 && cs->leafc_store*maxReducePropLeaf > potReduceLeafCarbon){
            tmpRatio = potReduceLeafCarbon/cs->leafc_store;
            exceedC += potReduceLeafCarbon;
            exceedN += ns->leafn_store * tmpRatio;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->leafc_store *= tmpRatio;
            ns->leafn_store *= tmpRatio;
            potReduceLeafCarbon=0.0;
        }else{
            exceedC += cs->leafc_store*maxReducePropLeaf;
            exceedN += ns->livecrootn_store*maxReducePropLeaf;
            potReduceLeafCarbon -= cs->leafc_store*maxReducePropLeaf;
            cs->leafc_store *= 1.0-maxReducePropLeaf;
            ns->leafn_store *= 1.0-maxReducePropLeaf;
        }// end of if
    }// end of if
    if(epc.veg_type != NON_VEG && potReduceLeafCarbon>0){
        if(cs->leafc >0 && cs->leafc*maxReducePropLeaf > potReduceLeafCarbon){
            tmpRatio = potReduceLeafCarbon/cs->leafc;
            double N2liter = potReduceLeafCarbon/epc.leaflitr_cn;
            exceedN += ns->leafn * tmpRatio - N2liter;
            
            tmpMultilper = epc.leaflitr_flab * stratum[0].cover_fraction;
            patch[0].litter_cs.litr1c += potReduceLeafCarbon * tmpMultilper;
            patch[0].litter_ns.litr1n += N2liter * tmpMultilper;
            tmpMultilper = epc.leaflitr_fucel * stratum[0].cover_fraction;
            patch[0].litter_cs.litr2c += potReduceLeafCarbon * tmpMultilper;
            patch[0].litter_ns.litr2n += N2liter * tmpMultilper;
            tmpMultilper = epc.leaflitr_fscel * stratum[0].cover_fraction;
            patch[0].litter_cs.litr3c += potReduceLeafCarbon * tmpMultilper;
            patch[0].litter_ns.litr3n += N2liter * tmpMultilper;
            tmpMultilper = epc.leaflitr_flig * stratum[0].cover_fraction;
            patch[0].litter_cs.litr4c += potReduceLeafCarbon * tmpMultilper;
            patch[0].litter_ns.litr4n += N2liter * tmpMultilper;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->leafc *= tmpRatio;
            ns->leafn *= tmpRatio;
            potReduceLeafCarbon=0.0;
        }else if(cs->leafc >0){
            //printf();
            double tmpC = maxReducePropLeaf * cs->leafc;
            tmpRatio = maxReducePropLeaf;
            double N2liter = tmpC/epc.leaflitr_cn;
            exceedN += ns->leafn * tmpRatio - N2liter;
            
            tmpMultilper = epc.leaflitr_flab * stratum[0].cover_fraction;
            patch[0].litter_cs.litr1c += tmpC * tmpMultilper;
            patch[0].litter_ns.litr1n += N2liter * tmpMultilper;
            tmpMultilper = epc.leaflitr_fucel * stratum[0].cover_fraction;
            patch[0].litter_cs.litr2c += tmpC * tmpMultilper;
            patch[0].litter_ns.litr2n += N2liter * tmpMultilper;
            tmpMultilper = epc.leaflitr_fscel * stratum[0].cover_fraction;
            patch[0].litter_cs.litr3c += tmpC * tmpMultilper;
            patch[0].litter_ns.litr3n += N2liter * tmpMultilper;
            tmpMultilper = epc.leaflitr_flig * stratum[0].cover_fraction;
            patch[0].litter_cs.litr4c += tmpC * tmpMultilper;
            patch[0].litter_ns.litr4n += N2liter * tmpMultilper;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->leafc *= tmpRatio;
            ns->leafn *= tmpRatio;
            potReduceLeafCarbon -= tmpC;
        }// end of if
    }// end of if
    
    // ---- froot
    if(epc.veg_type != NON_VEG && potReduceFrootCarbon>0){
        if(cs->frootc_transfer>0 && cs->frootc_transfer > potReduceFrootCarbon){
            tmpRatio = potReduceFrootCarbon/cs->frootc_transfer;
            exceedC += potReduceFrootCarbon;
            exceedN += ns->frootn_transfer * tmpRatio;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->frootc_transfer *= tmpRatio;
            ns->frootn_transfer *= tmpRatio;
            potReduceFrootCarbon=0.0;
        }else{
            exceedC += cs->frootc_transfer;
            exceedN += ns->frootn_transfer;
            potReduceFrootCarbon -= cs->frootc_transfer;
            cs->frootc_transfer = 0.0;
            ns->frootn_transfer = 0.0;
        }// end of if
    }// end of if
    if(epc.veg_type != NON_VEG && potReduceFrootCarbon>0){
        if(cs->frootc_store >0 && cs->frootc_store > potReduceFrootCarbon){
            tmpRatio = potReduceFrootCarbon/cs->frootc_store;
            exceedC += potReduceFrootCarbon;
            exceedN += ns->frootn_store * tmpRatio;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->frootc_store *= tmpRatio;
            ns->frootn_store *= tmpRatio;
            potReduceFrootCarbon=0.0;
        }else{
            exceedC += cs->frootc_store;
            exceedN += ns->frootn_store;
            potReduceFrootCarbon -= cs->frootc_store;
            cs->frootc_store = 0.0;
            ns->frootn_store = 0.0;
        }// end of if
    }
    if(epc.veg_type != NON_VEG && potReduceFrootCarbon>0){
        if(cs->frootc >0 && cs->frootc*0.2 > potReduceFrootCarbon){
            tmpRatio = potReduceFrootCarbon/cs->frootc;
            double N2liter = potReduceFrootCarbon/epc.froot_cn;
            
            tmpMultilper = epc.frootlitr_flab * stratum[0].cover_fraction;
            patch[0].litter_cs.litr1c += potReduceFrootCarbon * tmpMultilper;
            patch[0].litter_ns.litr1n += N2liter * tmpMultilper;
            tmpMultilper = epc.frootlitr_fucel * stratum[0].cover_fraction;
            patch[0].litter_cs.litr2c += potReduceFrootCarbon * tmpMultilper;
            patch[0].litter_ns.litr2n += N2liter * tmpMultilper;
            tmpMultilper = epc.frootlitr_fscel * stratum[0].cover_fraction;
            patch[0].litter_cs.litr3c += potReduceFrootCarbon * tmpMultilper;
            patch[0].litter_ns.litr3n += N2liter * tmpMultilper;
            tmpMultilper = epc.frootlitr_flig * stratum[0].cover_fraction;
            patch[0].litter_cs.litr4c += potReduceFrootCarbon * tmpMultilper;
            patch[0].litter_ns.litr4n += N2liter * tmpMultilper;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->frootc *= tmpRatio;
            ns->frootn *= tmpRatio;
            potReduceFrootCarbon=0.0;
        }else if(cs->frootc >0){
            double tmpC = 0.2 * cs->frootc;
            tmpRatio = 0.2;
            double N2liter = tmpC/epc.froot_cn;
            
            tmpMultilper = epc.frootlitr_flab * stratum[0].cover_fraction;
            patch[0].litter_cs.litr1c += tmpC * tmpMultilper;
            patch[0].litter_ns.litr1n += N2liter * tmpMultilper;
            tmpMultilper = epc.frootlitr_fucel * stratum[0].cover_fraction;
            patch[0].litter_cs.litr2c += tmpC * tmpMultilper;
            patch[0].litter_ns.litr2n += N2liter * tmpMultilper;
            tmpMultilper = epc.frootlitr_fscel * stratum[0].cover_fraction;
            patch[0].litter_cs.litr3c += tmpC * tmpMultilper;
            patch[0].litter_ns.litr3n += N2liter * tmpMultilper;
            tmpMultilper = epc.frootlitr_flig * stratum[0].cover_fraction;
            patch[0].litter_cs.litr4c += tmpC * tmpMultilper;
            patch[0].litter_ns.litr4n += N2liter * tmpMultilper;
            tmpRatio = max(0.0, 1.0 - tmpRatio);
            cs->frootc *= tmpRatio;
            ns->frootn *= tmpRatio;
            potReduceFrootCarbon -= tmpC;
        }
    }// end of if
    
    
    /*--------------------------------------------------------------*/
    /*  maintain allocation ratios  */
    /*--------------------------------------------------------------*/
    existingLeafCarbon = cs->leafc + cs->leafc_store + cs->leafc_transfer;
    existingFrootCarbon = cs->frootc + cs->frootc_store + cs->frootc_transfer;
    if (epc.veg_type == TREE) {
        existingLiveStemCarbon = cs->live_stemc + cs->livestemc_store + cs->livestemc_transfer;
        existingLiveCrootCarbon = cs->live_crootc + cs->livecrootc_store + cs->livecrootc_transfer;
    }else{
        existingLiveStemCarbon = 0.0;
        existingLiveCrootCarbon = 0.0;
    }
    existingLivebiomass = existingLeafCarbon + existingFrootCarbon + existingLiveStemCarbon + existingLiveCrootCarbon;
    
    minLeafCarbon = existingLivebiomass * allometric_leafRatio;
    minFrootCarbon = existingLivebiomass * allometric_frootRatio;
    minWoodCarbon = existingLivebiomass * allometric_livewoodRatio;

    double miss_Leafcarbon = minLeafCarbon - existingLeafCarbon;
    double miss_Frootcarbon = minFrootCarbon - existingFrootCarbon;
    double miss_Woodcarbon = minWoodCarbon - existingLiveStemCarbon - existingLiveCrootCarbon;
    
    if(miss_Leafcarbon>0 && exceedC>0 && exceedN>0){
        double withdrawC = min(min(exceedC, miss_Leafcarbon), (ns->npool+exceedN)*epc.leaf_cn); // constrained by carbon and nitrogen
        double withdrawN = min(ns->npool+exceedN,withdrawC/epc.leaf_cn);
        double withdrawNfromExceedN = min(exceedN, withdrawN);
        
        cs->leafc_store += withdrawC;
        ns->leafn_store += withdrawN;
        ns->npool -= withdrawN - withdrawNfromExceedN;
        exceedN -= withdrawNfromExceedN;
        exceedC -= withdrawC;
    }//end of if
    
    if(miss_Frootcarbon>0 && exceedC>0 && exceedN>0){
        double withdrawC = min(min(exceedC, miss_Frootcarbon), (ns->npool+exceedN)*epc.froot_cn); // constrained by carbon and nitrogen
        double withdrawN = min(ns->npool+exceedN,withdrawC/epc.froot_cn);
        double withdrawNfromExceedN = min(exceedN, withdrawN);
        
        cs->frootc_store += withdrawC;
        ns->frootn_store += withdrawN;
        ns->npool -= withdrawN - withdrawNfromExceedN;
        exceedN -= withdrawNfromExceedN;
        exceedC -= withdrawC;
    }//end of if
    
    if(epc.veg_type == TREE && miss_Woodcarbon>0 && exceedC>0 && exceedN>0){
        double withdrawC = min(min(exceedC, miss_Woodcarbon), (ns->npool+exceedN)*epc.livewood_cn); // constrained by carbon and nitrogen
        double withdrawN = min(ns->npool+exceedN,withdrawC/epc.livewood_cn);
        double withdrawNfromExceedN = min(exceedN, withdrawN);
        double stemRatio = existingLiveStemCarbon/(existingLiveStemCarbon+existingLiveCrootCarbon);
        
        cs->livestemc_store += withdrawC * stemRatio;
        ns->livestemn_store += withdrawN * stemRatio;
        stemRatio = max(0.0,min(1.0,1.0 - stemRatio));
        cs->livecrootc_store += withdrawC * stemRatio;
        ns->livecrootn_store += withdrawN * stemRatio;
    
        ns->npool -= withdrawN - withdrawNfromExceedN;
        exceedN -= withdrawNfromExceedN;
        exceedC -= withdrawC;
    }//end of if
    
    if(exceedC>0 && exceedN>0){
        patch[0].soil_cs.DOC += exceedC;
        patch[0].soil_ns.DON += exceedN;
    }// end of if
    
    
    /*--------------------------------------------------------------*/
    /*  push store to transfer  */
    /*--------------------------------------------------------------*/
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
    existingLeafCarbon = cs->leafc + cs->leafc_store + cs->leafc_transfer;
    

    
	/*--------------------------------------------------------------*/
	/* finally if there is really nothing, restart with small amount of growth   */
	/*	we allow only a certain amount of resprouting based on 	*/
	/*	a stratum default file parameterization 		*/
	/*--------------------------------------------------------------*/
    // annual carbon budget evaluate = cpool + (GPP-respG-repsM) = cpool + PSN - repM --> PSN >= repM + cpool --> need to adjust "repM"
    
    if( (epc.veg_type == TREE && (potReduceLeafCarbon>1e-8 || potReduceCrootCarbon>1e-8 || potReduceStemCarbon>1e-8 || potReduceFrootCarbon>1e-8)) ||
        (epc.veg_type == GRASS &&(potReduceLeafCarbon>1e-8 || potReduceFrootCarbon>1e-8))){
        
        //------- print out report
        printf("report,%d,%d,%d,%d,%d,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e\n",
               patch[0].ID, stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year,
               stratum->gDayCount,
               stratum->nFactor *what, // 1 is good
               stratum->wFactor *what, // 1 is good
               stratum->lFactor *what, // large is good
               stratum->gFactor *what, // 1 is good; actual/potential
               patch[0].sat_def_head,
               stratum[0].cover_fraction,
               stratum[0].cs.cpool,
               stratum->gwPSN *what,   // flux
               stratum->gwMResp *what, // flux
               stratum->gwAPAR *what, // 1 is good. stratum[0].gwAPAR += stratum[0].mult_conductance.APAR;
               stratum->gwLWP *what, // 1 is good. stratum[0].gwLWP += stratum[0].mult_conductance.LWP;
               stratum->gwVPD *what, // 1 is good. stratum[0].gwVPD += stratum[0].mult_conductance.vpd;
               patch[0].sat_deficit_z,          //<<---------------- instant
               patch[0].Ksat_vertical, //
               cs->num_resprout,
               potReduceLeafCarbon, potReduceCrootCarbon, potReduceStemCarbon, potReduceFrootCarbon,
               (stratum[0].cs.leafc + stratum[0].cs.leafc_store + stratum[0].cs.leafc_transfer),
               (stratum[0].cs.live_crootc + stratum[0].cs.livecrootc_store + stratum[0].cs.livecrootc_transfer),
               (stratum[0].cs.live_stemc + stratum[0].cs.livestemc_store + stratum[0].cs.livestemc_transfer),
               (stratum[0].cs.frootc + stratum[0].cs.frootc_store + stratum[0].cs.frootc_transfer)
               );
        
        // report,patch,vegid,day,month,year,gday,nfactor,wfactor,lfactor,gfactor,basement,cover,cpool,gwPSN,gwMresp,gwAPAR,gwLWP,gwVPD,wtz,perivnss,,num_resprout,redLeafC,redCrootc,redStemC,redFrootC,leafc,rootc,stemc,frootc
        cs->num_resprout += 1;
        
        if (cs->num_resprout > 5 ) {
            // SLB whole basin grass 0.307539*1.5=0.4613085 LAI; tree 0.406821*4.5=1.830694 LAI
            // SLB subbain grass 0.257967*1.5=0.3869505 basin LAI; tree 0.368434*4.5=1.657953 basin LAI
//            printf("Resprouting stratum [%d,%d: %d,%d,%d]{%e,%e} (%e,%e[%e,%e,%e]:%e:%e:%e)->(%e,%e,%e);[%e,%e,%e,%e][%e,%e,%e,%e]\n",
//                   patch[0].ID,stratum->defaults[0][0].ID, current_date.day, current_date.month, current_date.year,
//                   patch[0].Ksat_vertical,stratum[0].cover_fraction,
//                   0.1*minLeafCarbon, existingLeafCarbon, cs->leafc, cs->leafc_store, cs->leafc_transfer,
//                   existingFrootCarbon, existingLiveStemCarbon, existingLiveCrootCarbon,
//                   miss_Leafcarbon, miss_Frootcarbon, miss_Woodcarbon,
//                   patch[0].sat_deficit_z, patch[0].rootzone.S, patch[0].soil_ns.nitrate, patch[0].soil_ns.sminn,
//                   patch[0].soil_cs.soil1c,patch[0].soil_cs.soil2c,patch[0].soil_cs.soil3c,patch[0].soil_cs.soil4c
//                   //stratum->mult_conductance.APAR, stratum->mult_conductance.LWP, stratum->mult_conductance.vpd
//                   );
            
            // no litter or OM for decay!
            
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

    
    stratum->gDayCount=0;
    stratum->nFactor=0.0; // tracking @ allocate_daily_growth  <<--------- need to check
    stratum->wFactor=0.0; // tracking @ patch_daily_F          <<--------- need to check
    stratum->lFactor=0.0; // tracking @ canopy_stratum_daily_F <<--------- need to check
    stratum->gFactor=0.0; // tracking @ canopy_stratum_daily_F <<--------- actual gl vs potential gl (what's difference?)
    stratum->gwPSN=0.0; // tracking @ canopy_stratum_daily_F
    stratum->gwMResp=0.0; // tracking @ canopy_stratum_daily_F
    stratum->gwMRespLeaf=0.0; // tracking @ canopy_stratum_daily_F
    stratum->gwAPAR=0.0; // tracking @ canopy_stratum_daily_F
    stratum->gwLWP=0.0; // tracking @ canopy_stratum_daily_F
    stratum->gwVPD=0.0; // tracking @ canopy_stratum_daily_F
    

	

	return (!ok);
} /* end allocate_annual_growth */
