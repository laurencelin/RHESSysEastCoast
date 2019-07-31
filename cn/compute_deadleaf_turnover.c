/*--------------------------------------------------------------*/
/* 								*/
/*			compute_deadleaf_turnover.c			*/
/*								*/
/*								*/
/*	NAME							*/
/*	compute_deadleaf_turnover - computes turnover of standing but dead grass	*/
/*								(so that it enters litter pool		*/
/*								*/
/*	SYNOPSIS						*/
/*	double	compute_deadleaf_turnover( 				*/
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

int	compute_deadleaf_turnover(
							  struct epconst_struct epc,
							  struct epvar_struct *epv,
							  double cover_fraction,
							  struct cstate_struct *cs,
							  struct nstate_struct *ns,
							  struct litter_c_object *cs_litr,
							  struct litter_n_object *ns_litr,
                              struct patch_object *patch,
							  struct cdayflux_patch_struct *cdf,
							  struct ndayflux_patch_struct *ndf,
							  int grow_flag,
                              int BGC_flag)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
//    printf("compute deadleaf turnover: %lf, %lf, %lf, %lf, %lf, %lf\n", cs->leafc/ns->leafn, cs->frootc/ns->frootn,
//           cs->live_stemc/ns->live_stemn, cs->dead_stemc/ns->dead_stemn,
//           cs->live_crootc/ns->live_crootn, cs->dead_crootc/ns->dead_crootn);
//
    // in theory, when dead grass leaf passed into "deadleaf", the CN ratio should be leaflitrCN (when BGC_flag is on)
	//int ok=1;
	double c1,c2,c3,c4;
	double n1,n2,n3,n4;
	double turnover, avg_cn;
	turnover = epv->day_deadleaf_turnover;
    avg_cn = cs->dead_leafc/ns->dead_leafn;
    
	c1 = turnover * epc.leaflitr_flab;
	c2 = turnover * epc.leaflitr_fucel;
	n2 = BGC_flag==1? c2/avg_cn : c2/CEL_CN;
	c3 = turnover * epc.leaflitr_fscel;
	n3 = BGC_flag==1? c3/avg_cn : c3/CEL_CN;
	c4 = turnover * epc.leaflitr_flig;
	n4 = BGC_flag==1? c4/avg_cn : c4/LIG_CN;
    n1 = BGC_flag==1? c1/avg_cn : ((c1+c2+c3+c4)/epc.leaflitr_cn)-n2-n3-n4;
    
	/* set fluxes in daily flux structure */
	cdf->leafc_to_litr1c += c1 * cover_fraction;
	cdf->leafc_to_litr2c += c2 * cover_fraction;
	cdf->leafc_to_litr3c += c3 * cover_fraction;
	cdf->leafc_to_litr4c += c4 * cover_fraction;
	ndf->leafn_to_litr1n += n1 * cover_fraction;
	ndf->leafn_to_litr2n += n2 * cover_fraction;
	ndf->leafn_to_litr3n += n3 * cover_fraction;
	ndf->leafn_to_litr4n += n4 * cover_fraction;
	/* update state variables */
	cs->dead_leafc -= (c1 + c2 + c3 + c4);
	cs_litr->litr1c += c1 * cover_fraction;
	cs_litr->litr2c += c2 * cover_fraction;
	cs_litr->litr3c += c3 * cover_fraction;
	cs_litr->litr4c += c4 * cover_fraction;
	/* nitrogen state variable updates */
	ns->dead_leafn -= (n1 + n2 + n3 + n4);
	ns_litr->litr1n += n1 * cover_fraction;
	ns_litr->litr2n += n2 * cover_fraction;
	ns_litr->litr3n += n3 * cover_fraction;
	ns_litr->litr4n += n4 * cover_fraction;
    
    if(BGC_flag==0) patch[0].soil_ns.sminn += ((c2+c3+c4)/epc.leaflitr_cn -n2 -n3 -n4)* cover_fraction; // this should be zero when BGC_flag is on
    
	return(0);
} /*compute_deadleaf_turnover*/ 
