/*--------------------------------------------------------------*/
/* 								*/
/*			compute_leaf_litfall.c			*/
/*								*/
/*								*/
/*	NAME							*/
/*	compute_leaf_litfall - computes leaf litfall and 	*/
/*				updates carbon stores		*/
/*				for patch.			*/
/*								*/
/*	SYNOPSIS						*/
/*	double	compute_leaf_litfall( 				*/
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
#include "rhessys.h"
#include "phys_constants.h"

int	compute_leaf_litfall(
						 struct epconst_struct epc,
						 double litfallc,
						 double cover_fraction,
						 struct cstate_struct *cs,
						 struct nstate_struct *ns,
						 struct litter_c_object *cs_litr,
						 struct litter_n_object *ns_litr,
                         struct patch_object *patch,
						 struct cdayflux_patch_struct *cdf,
						 struct ndayflux_patch_struct *ndf,
						 struct cdayflux_struct *cdf_stratum,
						 struct ndayflux_struct *ndf_stratum,
						 int grow_flag,
                         int BGC_flag)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
//    printf("compute leaf litfall: %lf, %lf, %lf, %lf, %lf, %lf\n", cs->leafc/ns->leafn, cs->frootc/ns->frootn,
//           cs->live_stemc/ns->live_stemn, cs->dead_stemc/ns->dead_stemn,
//           cs->live_crootc/ns->live_crootn, cs->dead_crootc/ns->dead_crootn);
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	int ok=1;
	double c1,c2,c3,c4;
	double n1,n2,n3,n4;
	double nretrans, nloss;
	double avg_cn;
	
	avg_cn = cs->leafc/ns->leafn;

	/*------------------------------------------------------*/
	/*	Don't allow more leaves to fall than exist	*/
	/*------------------------------------------------------*/
	litfallc = max(0,min(cs->leafc, litfallc));
	/*------------------------------------------------------*/
	/*	determine carbon and nitrgoen to labile, cellulose and lignan pools */
	/*------------------------------------------------------*/
	
    
    if (grow_flag > 0){
        
        c1 = litfallc * epc.leaflitr_flab;
        c2 = litfallc * epc.leaflitr_fucel;
        n2 = BGC_flag==1? c2/epc.leaflitr_cn : c2/CEL_CN;
        c3 = litfallc * epc.leaflitr_fscel;
        n3 = BGC_flag==1? c3/epc.leaflitr_cn : c3/CEL_CN;
        c4 = litfallc * epc.leaflitr_flig;
        n4 = BGC_flag==1? c4/epc.leaflitr_cn : c4/LIG_CN;
        n1 = BGC_flag==1? c1/avg_cn : ((c1+c2+c3+c4)/avg_cn)-n2-n3-n4;
        
		nretrans = (litfallc/avg_cn) - (n1+n2+n3+n4);
        nretrans= max(nretrans, 0.0);


        if ((epc.veg_type == GRASS) || (epc.veg_type == C4GRASS)){
            cdf_stratum->leafc_to_deadleafc = litfallc;
            cs->leafc -= cdf_stratum->leafc_to_deadleafc;
            cs->dead_leafc += cdf_stratum->leafc_to_deadleafc;
            
            nloss = n1+n2+n3+n4+nretrans;
            ns->retransn += nretrans;
            ns->leafn -= nloss;
            ndf_stratum->leafn_to_deadleafn = nloss - nretrans ;
            ns->dead_leafn += ndf_stratum->leafn_to_deadleafn;
        }else{
            cdf->leafc_to_litr1c += c1 * cover_fraction;
            cdf->leafc_to_litr2c += c2 * cover_fraction;
            cdf->leafc_to_litr3c += c3 * cover_fraction;
            cdf->leafc_to_litr4c += c4 * cover_fraction;
            ndf->leafn_to_litr1n += n1 * cover_fraction;
            ndf->leafn_to_litr2n += n2 * cover_fraction;
            ndf->leafn_to_litr3n += n3 * cover_fraction;
            ndf->leafn_to_litr4n += n4 * cover_fraction;
            cs->leafc -= litfallc;
            cs_litr->litr1c += c1 * cover_fraction;// same as flux above; updating stage variable here
            cs_litr->litr2c += c2 * cover_fraction;
            cs_litr->litr3c += c3 * cover_fraction;
            cs_litr->litr4c += c4 * cover_fraction;
            nloss = n1+n2+n3+n4+nretrans;
            ns->retransn += nretrans;
            ns->leafn -= nloss;
            ns_litr->litr1n += n1 * cover_fraction;
            ns_litr->litr2n += n2 * cover_fraction;
            ns_litr->litr3n += n3 * cover_fraction;
            ns_litr->litr4n += n4 * cover_fraction;
            if(BGC_flag==0) patch[0].soil_ns.sminn += ((c2+c3+c4)/epc.leaflitr_cn -n2 -n3 -n4)* cover_fraction; // this should be zero when BGC_flag is on
        }
    }else{
        // not growth mode
		nretrans = 0.0;
        if ((epc.veg_type == GRASS) || (epc.veg_type == C4GRASS)){
            cdf_stratum->leafc_to_deadleafc = litfallc;
            /* update state variables */
            cs->leafc -= cdf_stratum->leafc_to_deadleafc;
            cs->dead_leafc += cdf_stratum->leafc_to_deadleafc;
            cs->leafc_store += litfallc;
            
            nloss = litfallc/avg_cn;
            ns->retransn += nretrans;
            ns->leafn -= nloss;
            ns->leafn_store += nloss;
            ndf_stratum->leafn_to_deadleafn = nloss - nretrans ;
            ns->dead_leafn += ndf_stratum->leafn_to_deadleafn;
        }else{
            cs->leafc -= litfallc;
            cs->leafc_store += litfallc;
            
            nloss = litfallc/avg_cn;
            ns->retransn += nretrans;
            ns->leafn -= nloss;
            ns->leafn_store += nloss;
        }
        
    }
	
	return(0);
} /*compute_leaf_litfall*/ 
