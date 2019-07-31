/*--------------------------------------------------------------*/
/* 								*/
/*		update_dissolved_organic_losses					*/
/*								*/
/*								*/
/*	NAME							*/
/*	update_dissolved_organic_losses -  					*/
/*		performs decomposition and updates soil/litter	*/
/*		carbon and nitrogen stores			*/
/*								*/
/*	SYNOPSIS						*/
/*	int update_dissolved_organic_losses(					*/
/*			double,					*/
/*			double,					*/
/*			double,					*/
/*			double,					*/
/*			struct	soil_c_object	*		*/
/*			struct	soil_n_object	*		*/
/*			struct	litter_c_object	*		*/
/*			struct	litter_n_object	*		*/
/*			struct	cdayflux_patch_object *		*/
/*			struct	ndayflux_patch_object *		*/
/*				)				*/
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
/*	assume DON at 3% mineralization (Vitousek et al,2000	*/
/*	Within-system element cycles, Input-output budgets      */
/*	and nutrient limitation)				*/
/*	DOC losses then follow DON				*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

int update_dissolved_organic_losses(
				  struct	date	current_date,
				  double 	DON_production_rate,
				  struct  soil_c_object   *cs_soil,
				  struct  soil_n_object   *ns_soil,
				  struct  litter_c_object *cs_litr,
				  struct  litter_n_object *ns_litr,
				  struct cdayflux_patch_struct *cdf,
				  struct ndayflux_patch_struct *ndf,
                  struct patch_object *patch,
                  int soilCNadaptation_falg)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	double cn_l1,cn_l2,cn_l3,cn_l4,cn_s1,cn_s2,cn_s3,cn_s4;
    double decayRate;
	int ok;

	ok=1;
	/* calculate litter and soil compartment C:N ratios */
    if ((cs_litr->litr1c > 0.0) && (ns_litr->litr1n > 0.0))    cn_l1 = cs_litr->litr1c/ns_litr->litr1n; else cn_l1 = 0.0;//LIVELAB_CN;
    if ((cs_litr->litr2c > 0.0) && (ns_litr->litr2n > 0.0))    cn_l2 = cs_litr->litr2c/ns_litr->litr2n; else cn_l2 = 0.0; //CEL_CN;
    if ((cs_litr->litr3c > 0.0) && (ns_litr->litr3n > 0.0))    cn_l3 = cs_litr->litr3c/ns_litr->litr3n; else cn_l3 = 0.0; //CEL_CN;
    if ((cs_litr->litr4c > 0.0) && (ns_litr->litr4n > 0.0))    cn_l4 = cs_litr->litr4c/ns_litr->litr4n; else cn_l4 = 0.0; // LIG_CN;
    
    if(soilCNadaptation_falg == 1 ){
        cn_s1 = patch[0].patch_SOIL1_CN;
        cn_s2 = patch[0].patch_SOIL2_CN;
        cn_s3 = patch[0].patch_SOIL3_CN;
        cn_s4 = patch[0].patch_SOIL4_CN;
    }else{
        cn_s1 = SOIL1_CN;
        cn_s2 = SOIL2_CN;
        cn_s3 = SOIL3_CN;
        cn_s4 = SOIL4_CN;
    }

    //--------------- DOM decay (turn on for upper 7 test)
    if( ns_soil->soil4n>0.0 && ns_soil->DON>0.0 && cs_soil->DOC>0.0){
        // > 0
        //decayRate = fabs(ndf->soil4n_to_sminn);
        //decayRate /= (ns_soil->soil4n + decayRate); // soil4n_to_sminn is negative, right?
//        if(patch[0].ID == 239202) printf("dissolveOM: %e %e %e\n", decayRate, ns_soil->DON * decayRate, cs_soil->DOC);
       
        decayRate = cdf->ks43;
        ns_soil->sminn += ns_soil->DON * decayRate;
        ns_soil->DON *= (1.0 - decayRate);
        cs_soil->DOC *= (1.0 - decayRate);
        
        patch[0].sat_NH4 += patch[0].sat_DON * decayRate;
        patch[0].sat_DON *= (1.0 - decayRate);
        patch[0].sat_DOC *= (1.0 - decayRate);
    }//
    
    // --------------- generate DOC
	/* labile litter fluxes */
	if (ndf->pmnf_l1s1 < 0.0 && cn_l1>0.0) {
		ndf->do_litr1n_loss = -DON_production_rate * ndf->pmnf_l1s1; //mineralization (negative value)
		cdf->do_litr1c_loss = max(0,min(ndf->do_litr1n_loss * cn_l1, cs_litr->litr1c)); // adding boundary mechanism here
        ndf->do_litr1n_loss = cdf->do_litr1c_loss/cn_l1;
    }else{
        ndf->do_litr1n_loss = 0.0;
        cdf->do_litr1c_loss = 0.0;
    }
	/* cellulose litter fluxes */
	if (ndf->pmnf_l2s2 <0.0 && cn_l2>0.0) {
		ndf->do_litr2n_loss = -DON_production_rate * ndf->pmnf_l2s2;
		cdf->do_litr2c_loss = max(0,min(ndf->do_litr2n_loss * cn_l2, cs_litr->litr2c));
        ndf->do_litr2n_loss = cdf->do_litr2c_loss/cn_l2;
    }else{
        ndf->do_litr2n_loss = 0.0;
        cdf->do_litr2c_loss = 0.0;
    }

	/* shielded cellulose litter fluxes */
	/* note these are based on lignan decay */
	if (ndf->pmnf_l4s3 <0.0 && cn_l3>0.0) {
		ndf->do_litr3n_loss = -DON_production_rate * ndf->pmnf_l3l2;
		cdf->do_litr3c_loss = max(0,min(ndf->do_litr3n_loss * cn_l3, cs_litr->litr3c));
        ndf->do_litr3n_loss = cdf->do_litr3c_loss/cn_l3;
    }else{
        ndf->do_litr3n_loss = 0.0;
        cdf->do_litr3c_loss = 0.0;
    }

	/* lignan litter fluxes  */
	if (ndf->pmnf_l4s3 <0.0 && cn_l4>0.0) {
		ndf->do_litr4n_loss = -DON_production_rate * ndf->pmnf_l4s3;
		cdf->do_litr4c_loss = max(0,min(ndf->do_litr4n_loss * cn_l4, cs_litr->litr4c));
        ndf->do_litr4n_loss = cdf->do_litr4c_loss/cn_l4;
    }else{
        ndf->do_litr4n_loss = 0.0;
        cdf->do_litr4c_loss = 0.0;
    }
	
	/* fast microbial recycling pool */
	if (ndf->pmnf_s1s2 <0.0 && cn_s1>0.0) {
		ndf->do_soil1n_loss = -DON_production_rate * ndf->pmnf_s1s2;
		cdf->do_soil1c_loss = max(0,min(ndf->do_soil1n_loss * cn_s1, cs_soil->soil1c));
        ndf->do_soil1n_loss = cdf->do_soil1c_loss/cn_s1;
    }else{
        ndf->do_soil1n_loss = 0.0;
        cdf->do_soil1c_loss = 0.0;
    }

	/* medium microbial recycling pool */
	if (ndf->pmnf_s2s3 <0.0 && cn_s2>0.0) {
		ndf->do_soil2n_loss = -DON_production_rate * ndf->pmnf_s2s3;
		cdf->do_soil2c_loss = max(0,min(ndf->do_soil2n_loss * cn_s2, cs_soil->soil2c));
        ndf->do_soil2n_loss = cdf->do_soil2c_loss/cn_s2;
    }else{
        ndf->do_soil2n_loss = 0.0;
        cdf->do_soil2c_loss = 0.0;
    }

	/* slow microbial recycling pool */
	if (ndf->pmnf_s3s4 <0.0 && cn_s3>0.0) {
		ndf->do_soil3n_loss = -DON_production_rate * ndf->pmnf_s3s4;
		cdf->do_soil3c_loss = max(0,min(ndf->do_soil3n_loss * cn_s3, cs_soil->soil3c));
        ndf->do_soil3n_loss = cdf->do_soil3c_loss/cn_s3;
    }else{
        ndf->do_soil3n_loss = 0.0;
        cdf->do_soil3c_loss = 0.0;
    }

	/* recalcitrant SOM pool (rf = 1.0, always mineralizing) */
	if (ndf->soil4n_to_sminn <0.0 && cn_s4>0.0) {
		ndf->do_soil4n_loss = -DON_production_rate * ndf->soil4n_to_sminn;
		cdf->do_soil4c_loss = max(0,min(ndf->do_soil4n_loss * cn_s4, cs_soil->soil4c));
        ndf->do_soil4n_loss = cdf->do_soil4c_loss/cn_s4;
    }else{
        ndf->do_soil4n_loss = 0.0;
        cdf->do_soil4c_loss = 0.0;
    }

	/* update soild and litter stores */
    cs_litr->litr1c       -= cdf->do_litr1c_loss;
    cs_litr->litr3c       -= cdf->do_litr3c_loss;
    cs_litr->litr2c       -= cdf->do_litr2c_loss;
    cs_litr->litr4c       -= cdf->do_litr4c_loss;
    cs_soil->soil1c       -= cdf->do_soil1c_loss;
    cs_soil->soil2c       -= cdf->do_soil2c_loss;
    cs_soil->soil3c       -= cdf->do_soil3c_loss;
    cs_soil->soil4c       -= cdf->do_soil4c_loss;

    ns_litr->litr1n       -= ndf->do_litr1n_loss;
    ns_litr->litr2n       -= ndf->do_litr2n_loss;
    ns_litr->litr3n       -= ndf->do_litr3n_loss;
    ns_litr->litr4n       -= ndf->do_litr4n_loss;
    ns_soil->soil1n       -= ndf->do_soil1n_loss;
    ns_soil->soil2n       -= ndf->do_soil2n_loss;
    ns_soil->soil3n       -= ndf->do_soil3n_loss;
    ns_soil->soil4n       -= ndf->do_soil4n_loss;


    if(cdf->do_soil1c_loss + cdf->do_soil2c_loss + cdf->do_soil3c_loss + cdf->do_soil4c_loss < 0){
        printf("update_dissolved_organic: %e\n", cdf->do_soil1c_loss + cdf->do_soil2c_loss + cdf->do_soil3c_loss + cdf->do_soil4c_loss );
    }
    
    if(cdf->do_litr1c_loss + cdf->do_litr2c_loss + cdf->do_litr3c_loss + cdf->do_litr4c_loss < 0){
        printf("update_dissolved_organic: %e\n", cdf->do_litr1c_loss + cdf->do_litr2c_loss + cdf->do_litr3c_loss + cdf->do_litr4c_loss );
    }
    
//    if(cs_soil->DOC<0 || ns_soil->DON<0) printf("update_dissolved_organic[%d] (%e,%e) + (%e,%e)\n", patch[0].ID, cs_soil->DOC, ns_soil->DON,
//                                                (cdf->do_soil1c_loss + cdf->do_soil2c_loss + cdf->do_soil3c_loss + cdf->do_soil4c_loss),
//                                                ( ndf->do_soil1n_loss + ndf->do_soil3n_loss + ndf->do_soil2n_loss + ndf->do_soil4n_loss)
//                                                );
    
	cs_soil->DOC +=  (cdf->do_soil1c_loss + cdf->do_soil2c_loss + cdf->do_soil3c_loss + cdf->do_soil4c_loss);
        // patch[0].surface_DOC += cdf->do_litr1c_loss + cdf->do_litr2c_loss + cdf->do_litr3c_loss + cdf->do_litr4c_loss  is done outside of this function
    
	ns_soil->DON  +=  ( ndf->do_soil1n_loss + ndf->do_soil3n_loss + ndf->do_soil2n_loss + ndf->do_soil4n_loss);
       // patch[0].surface_DON += ndf->do_litr1n_loss + ndf->do_litr2n_loss + ndf->do_litr3n_loss + ndf->do_litr4n_loss is done outside of this function


    

	return (!ok);
} /* end update_dissolved_organic_losses.c */
