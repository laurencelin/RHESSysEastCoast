/*--------------------------------------------------------------*/
/*                                                              */ 
/*		update_phenology				*/
/*                                                              */
/*  NAME                                                        */
/*		update_phenology				*/
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/* 	void update_phenology(										*/
/*			struct zone *										*/
/*			struct epv_struct *,		                        */
/*			struct epconst_struct,   		*/
/*                      struct phenology_struct *,              */
/*                      struct cstate_struct *,                 */
/*                      struct cdayflux_struct *,               */
/*                      struct cdayflux_patch_struct *,         */
/*                      struct nstate_struct *,                 */
/*                      struct ndayflux_struct *,               */
/*                      struct ndayflux_patch_struct *,         */
/*                      struct litter_c_object *,               */
/*                      struct litter_n_object *,               */
/*                      struct soil_c_object *,                 */
/*                      struct soil_n_object *,                 */
/*                      double,                                 */
/*			struct date current_date,		*/
/*		 	int	grow_flag);			*/
/*                                  				*/
/*                                                              */
/*  OPTIONS                                                     */
/*                                                              */
/*  DESCRIPTION                                                 */
/*                                                              */
/*	performs seasonal leaf drop and budding   		*/
/*	updates annnual allocates during leaf out		*/
/*	computes leaf and fine root litfall			*/
/*                                                              */
/*                                                              */
/*  PROGRAMMER NOTES                                            */
/*                                                              */
/*	modifed from phenology and prephenology in		*/
/*	P.Thornton (1997) version of 1d_bgc			*/
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"
#include "phys_constants.h"

void update_phenology(struct zone_object  *zone,
                      struct patch_object  *patch,
                      struct canopy_strata_object *stratum,
					  struct epvar_struct	*epv ,
					  struct epconst_struct	epc	,
					  struct phenology_struct *phen,	
					  struct cstate_struct *cs,
					  struct cdayflux_struct	*cdf,
					  struct cdayflux_patch_struct *cdf_patch,
					  struct nstate_struct *ns,
					  struct ndayflux_struct	*ndf,
					  struct ndayflux_patch_struct *ndf_patch,
					  struct litter_c_object *cs_litr,
					  struct litter_n_object *ns_litr,
					  struct litter_object *litter,
					  struct soil_c_object *cs_soil,
					  struct soil_n_object *ns_soil,
					  struct rooting_zone_object *rootzone,
					  double	effective_soil_depth,
					  double	cover_fraction,
					  double	gap_fraction,
					  double	theta_noon,
					  struct date current_date,
					  int grow_flag,
                      int dynRtZoff_flag,
                      int BGC_flag)
{
	/*--------------------------------------------------------------*/
	/*  Local function declaration                                  */
	/*--------------------------------------------------------------*/
	long	yearday( struct date);
	int	update_rooting_depth(
		struct rooting_zone_object *,
		double,
		double,
		double,
		double);
	void	update_litter_interception_capacity (double, 
		double,
		struct litter_c_object *,
		struct litter_object *);
	int	compute_annual_litfall(
		struct epconst_struct,
		struct phenology_struct *,
		struct cstate_struct *,
		int);
	int     compute_cwd_decay(	struct epconst_struct *,
		double,
		struct cstate_struct *,
		struct nstate_struct *,
		struct litter_c_object *,
		struct litter_n_object *,
		struct cdayflux_patch_struct *,
		struct ndayflux_patch_struct *,
		struct ndayflux_struct *,
        struct canopy_strata_object *stratum,
        struct patch_object    *patch,
        double,
        int);
	int	compute_deadleaf_turnover(
		struct epconst_struct,
		struct epvar_struct *,
		double,
		struct cstate_struct *,
		struct nstate_struct *,
		struct litter_c_object *,
		struct litter_n_object *,
        struct patch_object *,
		struct cdayflux_patch_struct *,
		struct ndayflux_patch_struct *,
		int, int);
	int	compute_leaf_litfall(
		struct epconst_struct,
		double ,
		double ,
		struct cstate_struct *,
		struct nstate_struct *,
		struct litter_c_object *,
		struct litter_n_object *,
        struct patch_object *,
		struct cdayflux_patch_struct *,
		struct ndayflux_patch_struct *,
		struct cdayflux_struct *,
		struct ndayflux_struct *,
		int, int);
	int	compute_froot_litfall(
		struct epconst_struct,
		double ,
		double ,
		struct cstate_struct *,
		struct nstate_struct *,
		struct litter_c_object *,
		struct litter_n_object *,
        struct patch_object *,
		struct cdayflux_patch_struct *,
		struct ndayflux_patch_struct *,
        int);
	
    //// ------------- need update
	double compute_growingseason_index(
		struct zone_object *,
	 	struct epvar_struct	*epv ,
		struct epconst_struct,
        struct phenology_struct *phen,
        double avg21tmin,
        double avg21dayl,
        double avg21vpd
		);
		
	/*--------------------------------------------------------------*/
	/*  Local variable definition.                                  */
	/*--------------------------------------------------------------*/

	int ok=1;
	long day;
	double perc_sunlit, leaflitfallc, frootlitfallc;
	double	rootc, sai, new_proj_lai_sunlit;
	double excess_n;
    double NO3_stored_transfer, plantcarbon; //these are for surface detention of N deposition in atmosphere. NO3_stored_transfer = from canopy to litter (surface N detention)
	int remdays_transfer;
	int expand_flag, litfall_flag;


	expand_flag = 0; 
	leaflitfallc = 0.0;
	frootlitfallc = 0.0;
	litfall_flag = 0;
	day = yearday(current_date);


 /*--------------------------------------------------------------*/
 /* static phenology                      */  
 /*--------------------------------------------------------------*/

	
  if (epc.phenology_flag == STATIC) {

        if (day == phen->expand_startday){
            expand_flag = 1;
            phen->gwseasonday = 1;
            phen->lfseasonday = -1;
        }

        else if (day == phen->litfall_startday){
            litfall_flag = 1;
            phen->lfseasonday = 1;
            phen->gwseasonday = -1;
        }

        else if (phen->gwseasonday > -1 && phen->gwseasonday <= epc.ndays_expand){
            expand_flag = 1;
        }

        else if (phen->lfseasonday > -1 && phen->lfseasonday <= epc.ndays_litfall){
            litfall_flag = 1;
        }

  } else {
  /*--------------------------------------------------------------*/
  /* dynamic phenology                      */  
  /*--------------------------------------------------------------*/
      int ii;
      for(ii=NUM_fday_Pred-1; ii>=0; ii--){
          // from future to present!
          if(   zone[0].METV_tmin_21ravgNext60[ii] != -999.0 &&
                zone[0].METV_dayl_21ravgNext60[ii] != -999.0 &&
                zone[0].METV_vpd_21ravgNext60[ii] != -999.0 ){
              
              phen->future_gsi[ii] = compute_growingseason_index(zone,
                                                           epv,
                                                           epc,
                                                           phen,
                                                           zone[0].METV_tmin_21ravgNext60[ii],
                                                           zone[0].METV_dayl_21ravgNext60[ii],
                                                           zone[0].METV_vpd_21ravgNext60[ii]); //future
          }else{
              phen->future_gsi[ii] = -999.0;
          }
          
          
          //if(zone[0].ID == 153354){printf("%ld,%ld,%f,%f,%f,%f\n",current_date.month,current_date.day, zone[0].METV_tmin_21ravgNext60[ii],zone[0].METV_dayl_21ravgNext60[ii],zone[0].METV_vpd_21ravgNext60[ii],future_gsi[ii] );}//debug
          
      }//ii
      
      
      
//      phen->gsi = compute_growingseason_index(zone,
//                                              epv,
//                                              epc,
//                                              zone[0].metv.tmin_21ravg[0],
//                                              zone[0].metv.vpd_21ravg[0],
//                                              zone[0].metv.dayl_21ravg[0]); //present GSI

      phen->gsi = phen->future_gsi[0];//present
      
      //------------ how to handle end of series? future_gsi[ii] = -999.0;
      //------------
      int numdayAhead;
      numdayAhead = epc.ndays_litfall-1;
      if( numdayAhead>=60) {numdayAhead = 59;}
      int GSI_fallSignal;
      if( phen->future_gsi[numdayAhead] > 0){
          // have data and firmly yes
          GSI_fallSignal = (phen->future_gsi[numdayAhead] < epc.gs_threshold_off);
      }else{
          // no data and use neighbor information (7 previous days)
          GSI_fallSignal = 0;
          for(ii=numdayAhead-7; ii<=numdayAhead; ii++){
              GSI_fallSignal += (phen->future_gsi[ii] < epc.gs_threshold_off); //T = 1; F = 0
          }//ii
          if(GSI_fallSignal > 3){GSI_fallSignal = 1;}
      }
      
      
      
      
      //------------------------------------------------------------ growth season starts
      //Q1: phen->gwseasonday = ? initially? what is it? it does not look like a flag!!
      //A1: it is not a flag, it's a counting.
      /* Advances seasonday variables */
//      if (phen->gwseasonday > 0)
//          phen->gwseasonday += 1;
//      if (phen->lfseasonday > 0)
//          phen->lfseasonday += 1;

      
      //  phen->gsi = phen->future_gsi[0];//present
      /* first are we before last possible date for leaf onset */
        /* are we already in a leaf onset condition */
        // unkonwn phen->gwseasonday initial value
          if (phen->gwseasonday > -1 ) { 
              if  (phen->gwseasonday <= epc.ndays_expand) // gwseasonday is increasing at the end of this page
                  expand_flag=1;
              }   
          //else if (phen->gsi > epc.gs_threshold_on && labs(phen->expand_startday - day)<epc.gs_window_on  ) { //stratum[0].defaults[0][0].epc.day_leafon
          else if (phen->gsi > epc.gs_threshold_on && labs(stratum[0].defaults[0][0].epc.day_leafon - day)<epc.gs_window_on  ) {
              // somehow spring out way earlier [don't understand at all]
                  phen->gwseasonday = 1;
                  phen->lfseasonday = -1;
                  expand_flag=1;
                  phen->expand_startday = day; //<<----- what's the use of this? it's a local variable (no use)
                  phen->expand_stopday = day + epc.ndays_expand; //<<----- what's the use of this?
              }   

    /* now determine if we are before the last possible date of leaf drop */
         
          /* are we already in a leaf offset */
          if (phen->lfseasonday > -1 ) { 

              phen->gwseasonday = -1; 
         
              if  (phen->lfseasonday <= epc.ndays_litfall)
                  litfall_flag=1;
              }
     
//        else if ((phen->gsi < epc.gs_threshold_off) && (phen->gwseasonday > epc.ndays_expand) && labs(phen->litfall_startday - day)<epc.gs_window_off ){
          //else if ( (GSI_fallSignal > 0) && labs(phen->litfall_startday - day)<epc.gs_window_off ){ //stratum[0].defaults[0][0].epc.day_leafoff
          else if ( (GSI_fallSignal > 0) && labs(stratum[0].defaults[0][0].epc.day_leafoff - day)<epc.gs_window_off ){
                    phen->lfseasonday = 1;
                    phen->gwseasonday = -1;
                    litfall_flag=1;
                    phen->litfall_startday = day;
                    phen->litfall_stopday = day + epc.ndays_litfall;
          }   

      

  } /* end dynamic phenology set up */
    

	
	phen->daily_allocation = epc.alloc_prop_day_growth;
	phen->annual_allocation = 0;
//    if (((cs->frootc + cs->frootc_transfer + cs->frootc_store) < ZERO) && (epc.phenology_flag != STATIC)) {
//        //printf("\n calling annual allocation because fine roots are zero %d %d %lf %lf %lf", zone->ID, epc.veg_type, cs->frootc, cs->frootc_transfer, cs->frootc_store);
//        //phen->annual_allocation=1;
//    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*--------------------------------------------------------------*/
	/*	Leaf expansion - spring leaf out			*/
	/*--------------------------------------------------------------*/

	if (expand_flag == 1) {
		remdays_transfer = max(1.0,(epc.ndays_expand + 1 - phen->gwseasonday)); ///<<-------- this is what Larry was talking about
        cdf->leafc_transfer_to_leafc = 2.0*cs->leafc_transfer / remdays_transfer;
        ndf->leafn_transfer_to_leafn = 2.0*ns->leafn_transfer / remdays_transfer;
        cdf->frootc_transfer_to_frootc=2.0*cs->frootc_transfer / remdays_transfer;
        ndf->frootn_transfer_to_frootn=2.0*ns->frootn_transfer / remdays_transfer;
        
//        double expRate;
//        double coef;
//        coef = 0.08333333 * epc.ndays_expand;
//        expRate = 1.0/(1.0 + exp(-coef*phen->gwseasonday+0.5*coef*epc.ndays_expand));
//        cdf->leafc_transfer_to_leafc = 2.0*cs->leafc_transfer * expRate;
//        ndf->leafn_transfer_to_leafn = 2.0*ns->leafn_transfer * expRate;
//        cdf->frootc_transfer_to_frootc=2.0*cs->frootc_transfer * expRate;
//        ndf->frootn_transfer_to_frootn=2.0*ns->frootn_transfer * expRate;
		
		if (epc.veg_type == TREE) {
			cdf->livestemc_transfer_to_livestemc
				= 2.0*cs->livestemc_transfer / remdays_transfer;
			ndf->livestemn_transfer_to_livestemn
				= 2.0*ns->livestemn_transfer / remdays_transfer;
			cdf->deadstemc_transfer_to_deadstemc
				= 2.0*cs->deadstemc_transfer / remdays_transfer;
			ndf->deadstemn_transfer_to_deadstemn
				= 2.0*ns->deadstemn_transfer / remdays_transfer;
			cdf->livecrootc_transfer_to_livecrootc
				= 2.0*cs->livecrootc_transfer / remdays_transfer;
			ndf->livecrootn_transfer_to_livecrootn
				= 2.0*ns->livecrootn_transfer / remdays_transfer;
			cdf->deadcrootc_transfer_to_deadcrootc
				= 2.0*cs->deadcrootc_transfer / remdays_transfer;
			ndf->deadcrootn_transfer_to_deadcrootn
				= 2.0*ns->deadcrootn_transfer / remdays_transfer;

		}
    }else{
        cdf->leafc_transfer_to_leafc = 0.0;
        ndf->leafn_transfer_to_leafn = 0.0;
        cdf->frootc_transfer_to_frootc = 0.0;
        ndf->frootn_transfer_to_frootn = 0.0;
        
        cdf->livestemc_transfer_to_livestemc = 0.0;
        ndf->livestemn_transfer_to_livestemn = 0.0;
        cdf->deadstemc_transfer_to_deadstemc = 0.0;
        ndf->deadstemn_transfer_to_deadstemn = 0.0;
        cdf->livecrootc_transfer_to_livecrootc = 0.0;
        ndf->livecrootn_transfer_to_livecrootn = 0.0;
        cdf->deadcrootc_transfer_to_deadcrootc = 0.0;
        ndf->deadcrootn_transfer_to_deadcrootn = 0.0;
    }// expend_flag
	/*--------------------------------------------------------------*/
	/*	Leaf drop - fall litter fall				*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/* at beginning of litter fall figure out how much to drop */
	/*--------------------------------------------------------------*/
	if (phen->lfseasonday == 1)  {
		ok = compute_annual_litfall(epc, phen, cs, grow_flag);
	}

	/*--------------------------------------------------------------*/
	/*	compute daily litter fall				*/
	/*--------------------------------------------------------------*/
	if (litfall_flag == 1) {
		remdays_transfer = max(1.0,(epc.ndays_litfall + 1 - phen->lfseasonday));

		leaflitfallc = 2.0*phen->leaflitfallc / remdays_transfer;
		if (leaflitfallc > cs->leafc)
			leaflitfallc = cs->leafc;
		frootlitfallc = 2.0*phen->frootlitfallc / remdays_transfer;
		if (frootlitfallc > cs->frootc)
			frootlitfallc = cs->frootc;
	}
	/*--------------------------------------------------------------*/
	/*	on the last day of litterfall make sure that deciduous no longer */
	/*	have any leaves left					*/
	/*	this is also considered to be the end of the growing season */
	/*	so do annual allcation					*/
	/*--------------------------------------------------------------*/
	if (phen->lfseasonday == epc.ndays_litfall){
		if (epc.phenology_type == DECID) {
			leaflitfallc = cs->leafc;
			phen->daily_allocation = 0;
			}
		phen->annual_allocation = 1;
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*--------------------------------------------------------------*/
	/*	update growth variables					*/
	/*		this is expression of allocation from 		*/
	/*		last seasons stored photosynthesis		*/
	/*	note all cdf and ndf variables are zero'd at the 	*/
	/*	start of the day, so only values set above are used	*/
	/*--------------------------------------------------------------*/
	 
	/* Leaf carbon transfer growth */
	cs->leafc            += cdf->leafc_transfer_to_leafc;
    cs->leafc_transfer   -= cdf->leafc_transfer_to_leafc; if(cs->leafc_transfer<0) cs->leafc_transfer=0.0;
	/* Leaf nitrogen transfer growth */
	ns->leafn           += ndf->leafn_transfer_to_leafn;
	ns->leafn_transfer  -= ndf->leafn_transfer_to_leafn; if(ns->leafn_transfer<0) ns->leafn_transfer=0.0;
	if (grow_flag > 0) {
		cs->frootc           += cdf->frootc_transfer_to_frootc;
        cs->frootc_transfer  -= cdf->frootc_transfer_to_frootc; if(cs->frootc_transfer<0) cs->frootc_transfer=0.0;
		ns->frootn          += ndf->frootn_transfer_to_frootn;
        ns->frootn_transfer -= ndf->frootn_transfer_to_frootn; if(ns->frootn_transfer<0) ns->frootn_transfer=0.0;
		if (epc.veg_type == TREE){
			/* Stem and coarse root transfer growth */
			cs->live_stemc             += cdf->livestemc_transfer_to_livestemc;
            cs->livestemc_transfer    -= cdf->livestemc_transfer_to_livestemc; if(cs->livestemc_transfer<0) cs->livestemc_transfer=0.0;
			cs->dead_stemc             += cdf->deadstemc_transfer_to_deadstemc;
            cs->deadstemc_transfer    -= cdf->deadstemc_transfer_to_deadstemc; if(cs->deadstemc_transfer<0) cs->deadstemc_transfer=0.0;
			cs->live_crootc            += cdf->livecrootc_transfer_to_livecrootc;
            cs->livecrootc_transfer   -= cdf->livecrootc_transfer_to_livecrootc; if(cs->livecrootc_transfer<0) cs->livecrootc_transfer=0.0;
			cs->dead_crootc            += cdf->deadcrootc_transfer_to_deadcrootc;
            cs->deadcrootc_transfer   -= cdf->deadcrootc_transfer_to_deadcrootc; if(cs->deadcrootc_transfer<0) cs->deadcrootc_transfer=0.0;
			/* nitrogen transfer */
			ns->live_stemn           += ndf->livestemn_transfer_to_livestemn;
            ns->livestemn_transfer  -= ndf->livestemn_transfer_to_livestemn; if(ns->livestemn_transfer<0) ns->livestemn_transfer=0.0;
			ns->dead_stemn           += ndf->deadstemn_transfer_to_deadstemn;
            ns->deadstemn_transfer  -= ndf->deadstemn_transfer_to_deadstemn; if(ns->deadstemn_transfer<0) ns->deadstemn_transfer=0.0;
			ns->live_crootn          += ndf->livecrootn_transfer_to_livecrootn;
            ns->livecrootn_transfer -= ndf->livecrootn_transfer_to_livecrootn; if(ns->livecrootn_transfer<0) ns->livecrootn_transfer=0.0;
			ns->dead_crootn          += ndf->deadcrootn_transfer_to_deadcrootn;
            ns->deadcrootn_transfer -= ndf->deadcrootn_transfer_to_deadcrootn; if(ns->deadcrootn_transfer<0) ns->deadcrootn_transfer=0.0;

		}
        
	}	/* end of grow processing */
    
    
	/*--------------------------------------------------------------*/
	/* check for leaf and fine root litfall for this day */
	/*--------------------------------------------------------------*/
    plantcarbon = cs->leafc + cs->live_stemc + cs->livestemc_store + cs->dead_stemc + cs->deadstemc_store;
	if ((leaflitfallc > ZERO ) && (cs->leafc > ZERO) && (ns->leafn > ZERO)){
		/*--------------------------------------------------------------*/
		/* set daily flux variables */
		/*--------------------------------------------------------------*/
		/*	compute leaf litter fall				*/
		/*--------------------------------------------------------------*/
		if (ok && compute_leaf_litfall(epc,
			leaflitfallc,cover_fraction,
			cs,ns, cs_litr,ns_litr,patch,cdf_patch,ndf_patch, cdf,ndf, grow_flag, BGC_flag)){
			fprintf(stderr,
				"FATAL ERROR: in leaf_litfall() from update_phenology()\n");
			exit(EXIT_FAILURE);
		}
		phen->leaflitfallc -= leaflitfallc;
		if (phen->leaflitfallc < 0)
			phen->leaflitfallc = 0;
        //plantcarbon = cs->leafc + cs->live_stemc + cs->livestemc_store + cs->dead_stemc + cs->deadstemc_store + cs->cwdc;
        if(plantcarbon>0){NO3_stored_transfer = stratum[0].NO3_stored * leaflitfallc/(plantcarbon);}else{NO3_stored_transfer = 0.0;}
        litter->NO3_stored += NO3_stored_transfer * stratum[0].cover_fraction;
        stratum[0].NO3_stored -= NO3_stored_transfer;
	}

	/*--------------------------------------------------------------*/
	/*	add additional leaf litterfall if water stress conditions */
	/*	occur							*/ 
	/*	only drop when accumulated litterfall due to water stress */
	/*	is greater than usual litterfall			*/
	/*--------------------------------------------------------------*/
	/*
	 if (epv->psi < epc.psi_open)	{
		perc_leaflitfall = (1.0 / (epc.psi_close - epc.psi_open) *
				epv->psi + (-epc.psi_open) /
				(epc.psi_close - epc.psi_open)) / 100.0;
		leaflitfallc = (perc_leaflitfall * cs->leafc);
		phen->leaflitfallc_wstress += leaflitfallc;
		if ((phen->leaflitfallc_wstress > phen->leaflitfallc) && 
			(phen->leaflitfallc_wstress < 1.5 * phen->leaflitfallc)) {
			if (ok && compute_leaf_litfall(epc,
				leaflitfallc,cover_fraction,
				cs,ns, cs_litr,ns_litr,cdf_patch,ndf_patch, cdf,ndf, grow_flag)){
				fprintf(stderr,
					"FATAL ERROR: in leaf_litfall() from update_phenology()\n");
				exit(EXIT_FAILURE);
				}
		}
	}
	*/

	if ((frootlitfallc > 0.0)  && (grow_flag > 0) && (cs->frootc > ZERO) && (ns->frootn > ZERO)){
		/*--------------------------------------------------------------*/
		/*	compute fine root turnover				*/
		/*--------------------------------------------------------------*/
		if (ok && compute_froot_litfall(epc,frootlitfallc,
			cover_fraction, cs,ns,cs_litr,
			ns_litr,patch,cdf_patch,ndf_patch, BGC_flag)){
			fprintf(stderr,
				"FATAL ERROR: in froot_litfall() from update_phenology()\n");
			exit(EXIT_FAILURE);
		}
		phen->frootlitfallc -= frootlitfallc;
		if (phen->frootlitfallc < 0)
			phen->frootlitfallc = 0;
	}//froot
	/*--------------------------------------------------------------*/
	/*	tree wood turnovers and dead standing grass turnovers	*/
	/*	- note this is not mortality but turnover rates		*/
	/*--------------------------------------------------------------*/

	if (((epc.veg_type == GRASS) || (epc.veg_type == C4GRASS)) && (grow_flag > 0) && (cs->dead_leafc > ZERO)) {
		if (ok && compute_deadleaf_turnover(epc,epv, cover_fraction, cs,ns,
			cs_litr,ns_litr,patch,cdf_patch,ndf_patch,grow_flag,BGC_flag)){
			fprintf(stderr,"FATAL ERROR: in compute_deadleaf_turnover() from update_phenology()\n");
			exit(EXIT_FAILURE);
		}
	}//GRASS

	if ((epc.veg_type == TREE) && (grow_flag > 0)) {
		/*--------------------------------------------------------------*/
		/*	compute coarse woody debris fragmentation		*/
		/*--------------------------------------------------------------*/
       //plantcarbon = cs->leafc + cs->live_stemc + cs->livestemc_store + cs->dead_stemc + cs->deadstemc_store + cs->cwdc;
       if ((cs->cwdc > ZERO) && (cover_fraction > ZERO)) {
			if (ok && compute_cwd_decay(&(epc),cover_fraction, cs,ns,cs_litr,
				ns_litr,cdf_patch,ndf_patch, ndf, stratum, patch, 0.02, BGC_flag)){
				fprintf(stderr,
					"FATAL ERROR: in cwd_decay() from update_phenology()\n");
				exit(EXIT_FAILURE);
			}
		}
		/*--------------------------------------------------------------*/
		/*	compute live steam and coarse root turnover		*/
		/* 	excess n from live stem turnover is added to retranslocated N */
		/*    nambiar et al., (1991) tree physiology			*/
		/*--------------------------------------------------------------*/

		if (cs->live_stemc > ZERO) {
			cdf->livestemc_to_deadstemc = min(epv->day_livestem_turnover, cs->live_stemc);
			
			ndf->livestemn_to_deadstemn= min(cdf->livestemc_to_deadstemc 
						/ epc.livewood_cn, ns->live_stemn);
		
			excess_n = max(0.0, ndf->livestemn_to_deadstemn -
						(cdf->livestemc_to_deadstemc / epc.deadwood_cn ) );
			ns->retransn += excess_n;
			cs->live_stemc -= cdf->livestemc_to_deadstemc;
			cs->dead_stemc += cdf->livestemc_to_deadstemc;
			ns->live_stemn -= ndf->livestemn_to_deadstemn;
			ns->dead_stemn += (ndf->livestemn_to_deadstemn - excess_n);
		}
		if (cs->live_crootc > ZERO) {
			cdf->livecrootc_to_deadcrootc = min(epv->day_livecroot_turnover, cs->live_crootc);
			ndf->livecrootn_to_deadcrootn= min(cdf->livecrootc_to_deadcrootc 
						/ epc.livewood_cn, ns->live_crootn);
		
			excess_n = max(0.0, ndf->livecrootn_to_deadcrootn -
						(cdf->livecrootc_to_deadcrootc / epc.deadwood_cn ) );
			ns->retransn += excess_n;
			cs->live_crootc -= cdf->livecrootc_to_deadcrootc;
			cs->dead_crootc += cdf->livecrootc_to_deadcrootc;
			ns->live_crootn -= ndf->livecrootn_to_deadcrootn;
			ns->dead_crootn += (ndf->livecrootn_to_deadcrootn - excess_n);
		}

	} /* end tree processing */
	if (grow_flag == 0) { /* non-grow processing */
		/*--------------------------------------------------------------*/
		/* update state variables assumes no retranslocation */
		/*	this is done in literfall routine			*/
		/*--------------------------------------------------------------*/
		/*
		cs->leafc        -= leaflitfallc;
		cs->leafc_store  += leaflitfallc;
		ns->leafn        -= leaflitfallc / epc.leaf_cn;
		ns->leafn_store  += leaflitfallc / epc.leaf_cn;
		*/
	}
	/*--------------------------------------------------------------*/
	/* check for rounding errors on end of litfall season */
	/*--------------------------------------------------------------*/
	if (fabs(cs->leafc) <= 1e-13){
		cs->leafc = 0.0;
		ns->leafn = 0.0;
	}
	if (fabs(cs->frootc) <= 1e-13){
		cs->frootc = 0.0;
		ns->frootn = 0.0;
	}



	/*--------------------------------------------------------------*/
	/*	compute new rooting depth based on current root carbon  */
	/*--------------------------------------------------------------*/
	rootc = cs->frootc+cs->live_crootc+cs->dead_crootc;
	if ((grow_flag > 0) && (rootc > ZERO) && (dynRtZoff_flag == 0) ){
		if (ok && update_rooting_depth(
			rootzone, rootc, epc.root_growth_direction, epc.root_distrib_parm,
			effective_soil_depth)){
			fprintf(stderr,
				"FATAL ERROR: in compute_rooting_depth() from update_phenology()\n");
			exit(EXIT_FAILURE);
		}
	}//rootdepth

	/*--------------------------------------------------------------*/
	/* now figure out a sunlit and shaded flux density		*/
	/* use Chen et al. 1999 Ecological Modelling 124: 99-119	*/
	/* to estimate shaded and sunlit fractions			*/
	/* then update lai based on sunlit/shaded sla			*/
	/* we need to do a predictor-corrector type convergence here	*/
	/*	since sunlit/shaded fraction depend on total proj_lai	*/
	/*--------------------------------------------------------------*/

	perc_sunlit = 0.0;
	if ((cs->leafc > ZERO) && (epc.veg_type != NON_VEG)) {
        // epv->proj_sla_sunlit = epc.proj_sla
        // epv->proj_sla_shade = epc.proj_sla * (1)
		epv->proj_lai = max( (cs->leafc * (epv->proj_sla_sunlit * perc_sunlit + epv->proj_sla_shade * (1-perc_sunlit))), 0.0);
        
		new_proj_lai_sunlit = 2.0 * cos(theta_noon) * (1.0 - exp(-0.5*(1-gap_fraction) * epv->proj_lai / cos(theta_noon)));
        
		while (fabs(epv->proj_lai_sunlit - new_proj_lai_sunlit) > 0.00001*new_proj_lai_sunlit )  {
			epv->proj_lai_sunlit = new_proj_lai_sunlit;
			epv->proj_lai_shade = epv->proj_lai - epv->proj_lai_sunlit;
			if ((epv->proj_lai_sunlit + epv->proj_lai_shade) > ZERO)
				perc_sunlit = (epv->proj_lai_sunlit) / (epv->proj_lai_sunlit + epv->proj_lai_shade);
			else
				perc_sunlit = 1.0;
			epv->proj_lai = max((cs->leafc * (epv->proj_sla_sunlit * perc_sunlit + 
				epv->proj_sla_shade * (1-perc_sunlit))), 0.0);
			new_proj_lai_sunlit = 2.0 * cos(theta_noon) *
					(1.0 - exp(-0.5*(1-gap_fraction)*
					epv->proj_lai / cos(theta_noon)));
			}
	}
	else {
		epv->proj_lai = 0.0;
		epv->proj_lai_sunlit = 0.0;
		epv->proj_lai_shade = 0.0;
	}


	/*--------------------------------------------------------------*/
	/* update lai based on sla			*/
	/* use sunlit sla for lai up to 1 and shaded sla for lai above that */
	/*--------------------------------------------------------------*/
	if ((epv->proj_lai_sunlit + epv->proj_lai_shade) > ZERO)
		perc_sunlit = (epv->proj_lai_sunlit) / (epv->proj_lai_sunlit + epv->proj_lai_shade);
	else
		perc_sunlit = 1.0;

	epv->all_lai = epv->proj_lai * epc.lai_ratio;

	if (epc.veg_type == TREE)  {
		sai = epc.proj_swa*(1.0-exp(-0.175*(cs->live_stemc+cs->dead_stemc)));
		epv->proj_pai = max(epv->proj_lai + sai, 0.0);
		epv->all_pai = max(epv->all_lai + sai, 0.0);
	}
	else {
		epv->proj_pai = epv->proj_lai;
		epv->all_pai = epv->all_lai;
	}
	/*--------------------------------------------------------------*/
	/*	update height						*/
	/*--------------------------------------------------------------*/
	if (epc.veg_type == TREE)
		if ( (cs->live_stemc + cs->dead_stemc) > ZERO)
			epv->height = epc.height_to_stem_coef
				* pow ( (cs->live_stemc + cs->dead_stemc), epc.height_to_stem_exp);
		else
			epv->height = 0.0;
	else
		if (epc.veg_type == NON_VEG) {
			epv->height = 0.0;
			}
		else {	if (cs->leafc > ZERO)
				epv->height = epc.height_to_stem_coef
					* pow ( (cs->leafc), epc.height_to_stem_exp);
			else
				epv->height = 0.0;
			}
	/*--------------------------------------------------------------*/
	/*	keep a seasonal max_lai for outputing purposes		*/
	/*--------------------------------------------------------------*/
	if (phen->gwseasonday == 1){
		epv->max_proj_lai = 0.0;
	}
	if (phen->gwseasonday > 1){
		if (epv->proj_lai > epv->max_proj_lai){
			epv->max_proj_lai = epv->proj_lai;
		}
	}
	/*--------------------------------------------------------------*/
	/*	update litter interception capacity			*/
	/*--------------------------------------------------------------*/
	update_litter_interception_capacity(
		litter->moist_coef,
		litter->density,	
		cs_litr,
		litter);


//    if(phen->gwseasonday > epc.ndays_expand && grow_flag>0){
//        if(cs->cpool>0){cs_soil->DOC += cs->cpool*0.05; cs->cpool*=0.95; }//1-0.05
//        if(ns->npool>0){ns_soil->DON += ns->npool*0.05; ns->npool*=0.95; }
//    }// postive CPOOL --> DOC
    

	/* Advances seasonday variables */
	  if (phen->gwseasonday > 0)
	      phen->gwseasonday += 1;
	  if (phen->lfseasonday > 0)
	      phen->lfseasonday += 1;
	

	return;

}/*end update phenology*/
