/*--------------------------------------------------------------*/
/* 								*/
/*		update_decomp					*/
/*								*/
/*								*/
/*	NAME							*/
/*	update_decomp -  					*/
/*		performs decomposition and updates soil/litter	*/
/*		carbon and nitrogen stores			*/
/*								*/
/*	SYNOPSIS						*/
/*	int update_decomp(					*/
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
/*		modified from Peter Thornton (1998)		*/
/*			dynamic - 1d-bgc ver4.0			*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

int update_decomp(
				  struct	date	current_date,
				  struct  soil_c_object   *cs_soil,
				  struct  soil_n_object   *ns_soil,
				  struct  litter_c_object *cs_litr,
				  struct  litter_n_object *ns_litr,
				  struct cdayflux_patch_struct *cdf,
				  struct ndayflux_patch_struct *ndf,
				  struct patch_object	*patch,
                  struct command_line_object *command_line
                  )
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	int ok = 1;
	double rfl1s1, rfl2s2,rfl4s3,rfs1s2,rfs2s3,rfs3s4;
	double cn_l1,cn_l3, cn_l2,cn_l4,cn_s1,cn_s2,cn_s3,cn_s4;
    double daily_net_nmin;
    int nlimit;
    double fpi;
	double total_N, total_preday_N, balance;
	double nitrate_immob, N_uptake, remaining_uptake;
    double microbal_immob, microbal_immoblitter, microbal_mineralize, holdingN;
    double checkl1, checkl2, checkl3, checkl4;
	double checks1, checks2, checks3, checks4;

    
    double totalNO3 = patch[0].rtzNO3 + patch[0].rtzSatNO3; //ns_soil->sminn;
    double totalNH4 = patch[0].rtzNH4 + patch[0].rtzSatNH4; //ns_soil->nitrate;
    
    
    
    
	total_preday_N = ns_litr->litr1n + ns_litr->litr2n +  ns_litr->litr3n
		+ ns_litr->litr4n + ns_soil->soil1n + ns_soil->soil2n + ns_soil->soil3n
		+ ns_soil->soil4n + totalNH4 + totalNO3;
	nlimit = ns_soil->nlimit;
	fpi = ns_soil->fract_potential_immob;
    if(fpi>1.0){fpi=1.0;}
    if(fpi<0.0){fpi=0.0;}
	/* now use the N limitation information fpi to assess the final decomposition
	fluxes. Mineralizing fluxes (pmnf* < 0.0) occur at the potential rate
	regardless of the competing N demands between microbial processes and
	plant uptake, but immobilizing fluxes are reduced when soil mineral
	N is limiting */
	/* calculate litter and soil compartment C:N ratios */
    if ((cs_litr->litr1c > 0.0) && (ns_litr->litr1n > 0.0))    cn_l1 = cs_litr->litr1c/ns_litr->litr1n; else cn_l1 = 0.0;//LIVELAB_CN;
    if ((cs_litr->litr2c > 0.0) && (ns_litr->litr2n > 0.0))    cn_l2 = cs_litr->litr2c/ns_litr->litr2n; else cn_l2 = 0.0; //CEL_CN;
    if ((cs_litr->litr3c > 0.0) && (ns_litr->litr3n > 0.0))    cn_l3 = cs_litr->litr3c/ns_litr->litr3n; else cn_l3 = 0.0; //CEL_CN;
    if ((cs_litr->litr4c > 0.0) && (ns_litr->litr4n > 0.0))    cn_l4 = cs_litr->litr4c/ns_litr->litr4n; else cn_l4 = 0.0; // LIG_CN;
    
    if(command_line[0].soilCNadaptation_falg == 1 ){
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
	/* respiration fractions for fluxes between compartments */
	rfl1s1 = 0.39;
	rfl2s2 = 0.55;
	rfl4s3 = 0.29;
	rfs1s2 = 0.28;
	rfs2s3 = 0.46;
	rfs3s4 = 0.55;
    daily_net_nmin = 0.0;
    microbal_immob = 0.0; microbal_immoblitter = 0.0; microbal_mineralize = 0.0;
    
//    if(fpi<1.0){
//        printf("actual decomp and uptake (before): %d, %lf,(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf), (%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf), (%lf,%lf,%lf,%lf)\n",
//                       patch[0].ID,fpi,
//                       cs_litr->litr1c, cs_litr->litr2c, cs_litr->litr3c, cs_litr->litr4c,
//                       cs_soil->soil1c, cs_soil->soil2c, cs_soil->soil3c, cs_soil->soil4c,
//                       ndf->pmnf_l1s1, ndf->pmnf_l2s2, ndf->pmnf_l3l2, ndf->pmnf_l4s3,
//                       ndf->pmnf_s1s2, ndf->pmnf_s2s3, ndf->pmnf_s3s4, ndf->pmnf_s4,
//                       microbal_immoblitter, microbal_immob, microbal_mineralize, ndf->sminn_to_npool);
//        printf("actual decomp and uptake (before): %d, %lf, (%lf,%lf,%lf), %lf, %lf, %lf\n",
//              patch[0].ID, fpi,
//              patch[0].NO3_throughfall,
//              patch[0].surface_NO3,
//              patch[0].sat_deficit_z,
//              patch[0].soil_ns.sminn,
//              patch[0].soil_ns.nitrate,
//              patch[0].soil_ns.fract_potential_uptake
//              );
//    }
    
//    if(patch[0].ID==239202) printf("soil2: patch[%d, %lf], CN[%lf,%lf,%lf,%lf, %lf,%lf,%lf,%lf], C[%e,%e,%e,%e, %e,%e,%e,%e], nflux[%e,%e,%e,%e, %e,%e,%e,%e]\n",
//           patch[0].ID, fpi,
//           cs_litr->litr1c/ns_litr->litr1n, cs_litr->litr2c/ns_litr->litr2n, cs_litr->litr3c/ns_litr->litr3n, cs_litr->litr4c/ns_litr->litr4n,
//           cs_soil->soil1c/ns_soil->soil1n, cs_soil->soil2c/ns_soil->soil2n, cs_soil->soil3c/ns_soil->soil3n, cs_soil->soil4c/ns_soil->soil4n,
//           cs_litr->litr1c, cs_litr->litr2c, cs_litr->litr3c, cs_litr->litr4c, cs_soil->soil1c, cs_soil->soil2c, cs_soil->soil3c, cs_soil->soil4c,
//           ndf->pmnf_l1s1, ndf->pmnf_l2s2, ndf->pmnf_l3l2, ndf->pmnf_l4s3, ndf->pmnf_s1s2, ndf->pmnf_s2s3, ndf->pmnf_s3s4, ndf->pmnf_s4
//        );

    //////////////////////////////////////////////// litter 1-4
	/* labile litter fluxes */
    checkl1 = ns_litr->litr1n;
	if ((cs_litr->litr1c > 0.0) && (ns_litr->litr1n > 0.0)){
        if(ndf->pmnf_l1s1>0.0){
            cdf->plitr1c_loss *= fpi;
            ndf->pmnf_l1s1 *= fpi;
            microbal_immoblitter += ndf->pmnf_l1s1;
        }else microbal_mineralize -= ndf->pmnf_l1s1;
		cdf->litr1c_hr = rfl1s1 * cdf->plitr1c_loss;
		cdf->litr1c_to_soil1c = (1.0 - rfl1s1) * cdf->plitr1c_loss;
		ndf->litr1n_to_soil1n = cdf->plitr1c_loss / cn_l1; // deltaN from decomposition corresponding to deltaC
        ndf->sminn_to_soil1n_l1 = ndf->pmnf_l1s1; //<<---- what is this?
        daily_net_nmin -= ndf->pmnf_l1s1; // just a variable to track NET_mineralization: negative is immob this time
        
        cs_litr->litr1c -= cdf->plitr1c_loss; if(cs_litr->litr1c<0.0) cs_litr->litr1c=0.0;
        ns_litr->litr1n -= ndf->litr1n_to_soil1n; if(ns_litr->litr1n<0.0) ns_litr->litr1n=0.0;
        if ((cs_litr->litr1c > 0.0) && (ns_litr->litr1n > 0.0))
            if( fabs(cs_litr->litr1c/ns_litr->litr1n - cn_l1) > 0.000001) printf("update demp, l1CN = %e(%e), %e, %e\n",cn_l1, cs_litr->litr1c/ns_litr->litr1n, cs_litr->litr1c, ns_litr->litr1n);
        
    }else{
//        cs_litr->litr1c = 0.0;
//        ns_litr->litr1n = 0.0;
        if ((cs_litr->litr1c < 0.0) || (ns_litr->litr1n < 0.0)) printf("lc1=%e, ln1=%e\n",cs_litr->litr1c, ns_litr->litr1n);
        cdf->litr1c_hr = 0.0;
        cdf->litr1c_to_soil1c = 0.0;
        ndf->litr1n_to_soil1n = 0.0;
        ndf->sminn_to_soil1n_l1 = 0.0;
    }
    checkl1 -= ns_litr->litr1n;
    checkl1 -= ndf->litr1n_to_soil1n;
    if( fabs(checkl1) > ZERO) printf("checkl1 = %e\n",checkl1);
    
    
	/* cellulose litter fluxes */
    checkl2 = ns_litr->litr2n;
	if ((cs_litr->litr2c > 0.0) && (ns_litr->litr2n > 0.0)){
		if(ndf->pmnf_l2s2>0.0){
            cdf->plitr2c_loss *= fpi;
            ndf->pmnf_l2s2 *= fpi;
            microbal_immoblitter += ndf->pmnf_l2s2;
        }else microbal_mineralize -= ndf->pmnf_l2s2;
		cdf->litr2c_hr = rfl2s2 * cdf->plitr2c_loss;
		cdf->litr2c_to_soil2c = (1.0 - rfl2s2) * cdf->plitr2c_loss;
		ndf->litr2n_to_soil2n = cdf->plitr2c_loss / cn_l2;
		ndf->sminn_to_soil2n_l2 = ndf->pmnf_l2s2;
		daily_net_nmin -= ndf->pmnf_l2s2;
        
        cs_litr->litr2c -= cdf->plitr2c_loss;
        ns_litr->litr2n -= ndf->litr2n_to_soil2n;
    }else{
//        cs_litr->litr2c = 0.0;
//        ns_litr->litr2n = 0.0;
        if ((cs_litr->litr2c < 0.0) || (ns_litr->litr2n < 0.0)) printf("lc2=%e, ln2=%e\n",cs_litr->litr2c, ns_litr->litr2n);
        cdf->litr2c_hr = 0.0;
        cdf->litr2c_to_soil2c = 0.0;
        ndf->litr2n_to_soil2n = 0.0;
        ndf->sminn_to_soil2n_l2 = 0.0;
    }
    checkl2 -= ns_litr->litr2n;
    checkl2 -= ndf->litr2n_to_soil2n;
    if( fabs(checkl2) > ZERO) printf("checkl2 = %e when %e\n",checkl2, ns_litr->litr2n);
    checkl2 = ns_litr->litr2n;
    
	/* release of shielded cellulose litter, tied to the decay rate of
	lignin litter */
    checkl3 = ns_litr->litr3n;
	if ((cs_litr->litr3c > 0.0) && (ns_litr->litr3n > 0.0)){
		if(ndf->pmnf_l3l2>0.0){
            cdf->plitr3c_loss *= fpi;
            ndf->pmnf_l3l2 *= fpi;
            microbal_immoblitter += ndf->pmnf_l3l2;
        }else microbal_mineralize -= ndf->pmnf_l3l2;
		cdf->litr3c_hr = rfl4s3 * cdf->plitr3c_loss;
		cdf->litr3c_to_litr2c = (1.0 - rfl4s3) * cdf->plitr3c_loss;
		ndf->litr3n_to_litr2n = cdf->plitr3c_loss / cn_l3;
		ndf->sminn_to_soil2n_l3 = ndf->pmnf_l3l2;
		daily_net_nmin -= ndf->pmnf_l3l2;
        
        cs_litr->litr3c -= cdf->plitr3c_loss;
        ns_litr->litr3n -= ndf->litr3n_to_litr2n;
            //putting here is ok since flux to l2 is zero if l3 is zero
        cs_litr->litr2c += cdf->litr3c_to_litr2c;
        ns_litr->litr2n += ndf->litr3n_to_litr2n;
        ns_litr->litr2n += ndf->sminn_to_soil2n_l3;
    }else{
//        cs_litr->litr3c = 0.0;
//        ns_litr->litr3n = 0.0;
        if ((cs_litr->litr3c < 0.0) || (ns_litr->litr3n < 0.0)) printf("lc3=%e, ln3=%e\n",cs_litr->litr3c, ns_litr->litr3n);
        cdf->litr3c_hr = 0.0;
        cdf->litr3c_to_litr2c = 0.0;
        ndf->litr3n_to_litr2n = 0.0;
        ndf->sminn_to_soil2n_l3 = 0.0;
    }
    checkl3 -= ns_litr->litr3n;
    checkl3 -= ndf->litr3n_to_litr2n;
    if( fabs(checkl3) > ZERO) printf("checkl3 = %e\n",checkl3);
    checkl3 = ns_litr->litr2n;
    checkl3 -= ndf->litr3n_to_litr2n;
    checkl3 -= ndf->sminn_to_soil2n_l3;
	if( fabs(checkl2-checkl3) > ZERO) printf("checkl2 in l3 process = %e\n",(checkl2-checkl3) );
    
    
	/* lignin litter fluxes */
    checkl4 = ns_litr->litr4n;
	if ((cs_litr->litr4c > 0.0) && (ns_litr->litr4n > 0.0)){
		if(ndf->pmnf_l4s3>0.0){
            cdf->plitr4c_loss *= fpi;
            ndf->pmnf_l4s3 *= fpi;
            microbal_immoblitter += ndf->pmnf_l4s3;
        }else microbal_mineralize -= ndf->pmnf_l4s3;
		cdf->litr4c_hr = rfl4s3 * cdf->plitr4c_loss;
		cdf->litr4c_to_soil3c = (1.0 - rfl4s3) * cdf->plitr4c_loss;
		ndf->litr4n_to_soil3n = cdf->plitr4c_loss / cn_l4;
		ndf->sminn_to_soil3n_l4 = ndf->pmnf_l4s3;
		daily_net_nmin -= ndf->pmnf_l4s3;
        
        cs_litr->litr4c -= cdf->plitr4c_loss;
        ns_litr->litr4n -= ndf->litr4n_to_soil3n;
    }else{
//        cs_litr->litr4c = 0.0;
//        ns_litr->litr4n = 0.0;
        if ((cs_litr->litr4c < 0.0) || (ns_litr->litr4n < 0.0)) printf("lc4=%e, ln4=%e\n",cs_litr->litr4c, ns_litr->litr4n);
        cdf->litr4c_hr = 0.0;
        cdf->litr4c_to_soil3c = 0.0;
        ndf->litr4n_to_soil3n = 0.0;
        ndf->sminn_to_soil3n_l4 = 0.0;
    }
    checkl4 -= ns_litr->litr4n;
    checkl4 -= ndf->litr4n_to_soil3n;
    if( fabs(checkl4) > ZERO) printf("checkl4 = %e\n",checkl4);
    
    
	//////////////////////////////////////////////// soil 1-4
	/* fast microbial recycling pool */
    checks1 = ns_soil->soil1n;
	if (cs_soil->soil1c > 0.0 && ns_soil->soil1n > 0.0){
		if(ndf->pmnf_s1s2>0.0){
            cdf->psoil1c_loss *= fpi;
            ndf->pmnf_s1s2 *= fpi;
            microbal_immob += ndf->pmnf_s1s2;
        }else microbal_mineralize -= ndf->pmnf_s1s2;
		cdf->soil1c_hr = rfs1s2 * cdf->psoil1c_loss;
		cdf->soil1c_to_soil2c = (1.0 - rfs1s2) * cdf->psoil1c_loss;
		ndf->soil1n_to_soil2n = cdf->psoil1c_loss / cn_s1;
		ndf->sminn_to_soil2n_s1 = ndf->pmnf_s1s2;
		daily_net_nmin -= ndf->pmnf_s1s2;
        
        cs_soil->soil1c -= cdf->psoil1c_loss;
        ns_soil->soil1n -= ndf->soil1n_to_soil2n;
    }else{
//        cs_soil->soil1c = 0.0;
//        ns_soil->soil1n = 0.0;
        cdf->soil1c_hr = 0.0;
        cdf->soil1c_to_soil2c = 0.0;
        ndf->soil1n_to_soil2n = 0.0;
        ndf->sminn_to_soil2n_s1 = 0.0;
    }
    checks1 -= ns_soil->soil1n;
    checks1 -= ndf->soil1n_to_soil2n;
    if( fabs(checks1) > ZERO) printf("checks1 = %e\n",checks1);
    checkl1 = ns_soil->soil1n;
        cs_soil->soil1c += cdf->litr1c_to_soil1c;
        ns_soil->soil1n += ndf->litr1n_to_soil1n;
        ns_soil->soil1n += ndf->sminn_to_soil1n_l1;
        checks1 = ns_soil->soil1n;
        checks1 -= ndf->litr1n_to_soil1n;
        checks1 -= ndf->sminn_to_soil1n_l1;
    if( fabs(checks1 - checkl1) > ZERO) printf("checks1 in process of adding ltr1 = %e when litr1n_to_soil1n = (%e) and sminn_to_soil1n_l1 = %e and soil1 [%e,%e,%e] and litr1 [%e,%e] and fpi = %lf [%d]\n",(checks1 - checkl1), ndf->litr1n_to_soil1n, ndf->sminn_to_soil1n_l1,
        cs_soil->soil1c, ns_soil->soil1n,checkl1,
        ns_litr->litr1n, cs_litr->litr1c, fpi, nlimit);
    
    
    
	/* medium microbial recycling pool */
    checks2 = ns_soil->soil2n;
	if (cs_soil->soil2c > 0.0 && ns_soil->soil2n > 0.0){
		if(ndf->pmnf_s2s3>0.0){
            cdf->psoil2c_loss *= fpi;
            ndf->pmnf_s2s3 *= fpi;
            microbal_immob += ndf->pmnf_s2s3;
        }else microbal_mineralize -= ndf->pmnf_s2s3;
		cdf->soil2c_hr = rfs2s3 * cdf->psoil2c_loss;
		cdf->soil2c_to_soil3c = (1.0 - rfs2s3) * cdf->psoil2c_loss;
		ndf->soil2n_to_soil3n = cdf->psoil2c_loss / cn_s2;
		ndf->sminn_to_soil3n_s2 = ndf->pmnf_s2s3;
		daily_net_nmin -= ndf->pmnf_s2s3;
        
        cs_soil->soil2c -= cdf->psoil2c_loss;
        ns_soil->soil2n -= ndf->soil2n_to_soil3n;
    }else{
//        cs_soil->soil2c = 0.0;
//        ns_soil->soil2n = 0.0;
        cdf->soil2c_hr = 0.0;
        cdf->soil2c_to_soil3c = 0.0;
        ndf->soil2n_to_soil3n = 0.0;
        ndf->sminn_to_soil3n_s2 = 0.0;
    }
    checks2 -= ns_soil->soil2n;
    checks2 -= ndf->soil2n_to_soil3n;
    if( fabs(checks2) > ZERO) printf("checks2 = %e\n",checks2);
    checkl2 = ns_soil->soil2n;
        cs_soil->soil2c += cdf->litr2c_to_soil2c;
        ns_soil->soil2n += ndf->litr2n_to_soil2n;
        ns_soil->soil2n += ndf->sminn_to_soil2n_l2;

        cs_soil->soil2c += cdf->soil1c_to_soil2c;
        ns_soil->soil2n += ndf->soil1n_to_soil2n;
        ns_soil->soil2n += ndf->sminn_to_soil2n_s1;
    checks2 = ns_soil->soil2n;
    checks2 -= ndf->litr2n_to_soil2n;
    checks2 -= ndf->soil1n_to_soil2n;
    checks2 -= ndf->sminn_to_soil2n_l2 + ndf->sminn_to_soil2n_s1;
    if( fabs(checks2-checkl2) > ZERO) printf("checks2 in process of adding litr2 and soil1 = %e when soil2n = %e, ndf->sminn_to_soil2n_l2 = %e, ndf->sminn_to_soil2n_s1 = %e\n",(checks2-checkl2), ns_soil->soil2n,ndf->sminn_to_soil2n_l2, ndf->sminn_to_soil2n_s1 );
//    if(patch[0].ID==239202) printf("soil2: patch[%d, %lf], decay[%e,%e], from l2 [%e,%e,%e], from s1 [%e, %e, %e] -> soil(%e, %e: %lf), soilCN[%lf,%lf,%lf,%lf, %lf,%lf,%lf,%lf]\n",
//           patch[0].ID, fpi,
//           cdf->psoil2c_loss, ndf->soil2n_to_soil3n,
//           cdf->litr2c_to_soil2c,ndf->litr2n_to_soil2n,ndf->sminn_to_soil2n_l2,
//           cdf->soil1c_to_soil2c,ndf->soil1n_to_soil2n,ndf->sminn_to_soil2n_s1,
//           cs_soil->soil2c,ns_soil->soil2n,cs_soil->soil2c/ns_soil->soil2n,
//           cs_litr->litr1c/ns_litr->litr1n, cs_litr->litr2c/ns_litr->litr2n, cs_litr->litr3c/ns_litr->litr3n, cs_litr->litr4c/ns_litr->litr4n,
//           cs_soil->soil1c/ns_soil->soil1n, cs_soil->soil2c/ns_soil->soil2n, cs_soil->soil3c/ns_soil->soil3n, cs_soil->soil4c/ns_soil->soil4n);
    
	/* slow microbial recycling pool */
	if (cs_soil->soil3c > 0.0 && ns_soil->soil3n > 0.0){
		if(ndf->pmnf_s3s4>0.0){
            cdf->psoil3c_loss *= fpi;
            ndf->pmnf_s3s4 *= fpi;
            microbal_immob += ndf->pmnf_s3s4;
        }else microbal_mineralize -= ndf->pmnf_s3s4;
		cdf->soil3c_hr = rfs3s4 * cdf->psoil3c_loss;
		cdf->soil3c_to_soil4c = (1.0 - rfs3s4) * cdf->psoil3c_loss;
		ndf->soil3n_to_soil4n = cdf->psoil3c_loss / cn_s3;//<<-------
		ndf->sminn_to_soil4n_s3 = ndf->pmnf_s3s4;
		daily_net_nmin -= ndf->pmnf_s3s4;
        
        cs_soil->soil3c -= cdf->psoil3c_loss;
        ns_soil->soil3n -= ndf->soil3n_to_soil4n;
    }else{
//        cs_soil->soil3c = 0.0;
//        ns_soil->soil3n = 0.0;
        cdf->soil3c_hr = 0.0;
        cdf->soil3c_to_soil4c = 0.0;
        ndf->soil3n_to_soil4n = 0.0;//<<-------
        ndf->sminn_to_soil4n_s3 = 0.0;
    }
    cs_soil->soil3c += cdf->litr4c_to_soil3c;
    ns_soil->soil3n += ndf->litr4n_to_soil3n;
    ns_soil->soil3n += ndf->sminn_to_soil3n_l4;
    cs_soil->soil3c += cdf->soil2c_to_soil3c;
    ns_soil->soil3n += ndf->soil2n_to_soil3n;
    ns_soil->soil3n += ndf->sminn_to_soil3n_s2;
    
	/* recalcitrant SOM pool (rf = 1.0, always mineralizing) */
	if (cs_soil->soil4c > 0.0 && ns_soil->soil4n > 0.0){
		cdf->soil4c_hr = cdf->psoil4c_loss;
		ndf->soil4n_to_sminn = cdf->psoil4c_loss / cn_s4;
		daily_net_nmin += ndf->soil4n_to_sminn;
        microbal_mineralize += ndf->soil4n_to_sminn;
        
        cs_soil->soil4c -= cdf->psoil4c_loss;
        ns_soil->soil4n -= ndf->soil4n_to_sminn;
    }else{
//        cs_soil->soil4c = 0.0;
//        ns_soil->soil4n = 0.0;
        cdf->soil4c_hr = 0.0;
        ndf->soil4n_to_sminn = 0.0;
    }
    cs_soil->soil4c += cdf->soil3c_to_soil4c;
    ns_soil->soil4n += ndf->soil3n_to_soil4n;
    ns_soil->soil4n += ndf->sminn_to_soil4n_s3;
    

    
    // need to deal with daily_net_nminlitter, daily_net_nmin, and ndf->sminn_to_npool at the same time, as well as totalNH4. (need to work)
    // daily_net_nmin = 0.0; microbal_immob = 0.0; microbal_immoblitter = 0.0; microbal_mineralize = 0.0;
    double Nnow = patch[0].litter.NO3_stored + patch[0].surface_NO3 + patch[0].surface_NH4 + totalNO3 + totalNH4;
    double wtfNow = microbal_immoblitter + microbal_immob + ndf->sminn_to_npool;
    
    //daily_net_nmin = 0.0; microbal_immob = 0.0; microbal_immoblitter = 0.0; //microbal_mineralize = 0.0; //<-- disable immob
    // fpi = ns_soil->fract_potential_immob;
    nitrate_immob = microbal_immoblitter; // litter microbail N demand is satisfied by surface and stored N
    holdingN = patch[0].litter.NO3_stored; patch[0].litter.NO3_stored -= max(0.0,min(holdingN, nitrate_immob)); nitrate_immob -= max(0.0,holdingN);
    if(nitrate_immob>0){holdingN = patch[0].surface_NO3; patch[0].surface_NO3 -= max(0.0, min(holdingN, nitrate_immob)); nitrate_immob -= max(0.0,holdingN);}else{nitrate_immob=0.0;}
    if(nitrate_immob>0){holdingN = patch[0].surface_NH4; patch[0].surface_NH4 -= max(0.0, min(holdingN, nitrate_immob)); nitrate_immob -= max(0.0,holdingN);}else{nitrate_immob=0.0;}
    if(nitrate_immob<0){nitrate_immob=0.0;}
    nitrate_immob += microbal_immob; // soilOM microbail N demand is satisfied by soil N
    nitrate_immob += ndf->sminn_to_npool; // plant N uptake
    
//    double totalNO3 = patch[0].rtzNO3 + patch[0].rtzSatNO3; //ns_soil->sminn;
//    double totalNH4 = patch[0].rtzNH4 + patch[0].rtzSatNH4; //ns_soil->nitrate;
    double delta_totalNO3=0.0; double delta_totalNH4=0.0;
    if(nitrate_immob>0){holdingN = totalNO3; delta_totalNO3 = max(0.0,min(holdingN, nitrate_immob)); nitrate_immob -= max(0.0,holdingN);}else{nitrate_immob=0.0;}
    if(nitrate_immob<0){nitrate_immob=0.0;}
    nitrate_immob -= microbal_mineralize; //any left is satisfied by mineralized and soil NH4
    if(nitrate_immob>0){holdingN = totalNH4; delta_totalNH4 = max(0.0,min(holdingN, nitrate_immob)); nitrate_immob -= max(0.0,holdingN);}else{ns_soil->sminn -= nitrate_immob;} // put excessive mineralized N to totalNH4
    //change in totalNO3 and totalNH4
    if(totalNO3>0 & delta_totalNO3>0){
        ns_soil->nitrate -= min(ns_soil->nitrate, delta_totalNO3 * patch[0].rtzNO3/totalNO3);
        patch[0].sat_NO3 -= min(patch[0].sat_NO3, delta_totalNO3 * patch[0].rtzSatNO3/totalNO3);
    }
    if(totalNH4>0 & delta_totalNH4>0){
        ns_soil->sminn -= min(ns_soil->sminn, delta_totalNH4 * patch[0].rtzNH4/totalNH4);
        patch[0].sat_NH4 -= min(patch[0].sat_NH4, delta_totalNH4 * patch[0].rtzSatNH4/totalNH4);
    }

    if (nitrate_immob > ZERO) printf("N balance issue: %f, %f, %f(%f), %f(%f) [%f, %f]:: ->%e<-%e (%f)[%f : %f]\n",
                                     microbal_immoblitter, microbal_immob,
                                  ndf->plant_potential_ndemand, ndf->sminn_to_npool,
                                  ndf->mineralized, microbal_mineralize,
                                     totalNO3, totalNH4,
                                     nitrate_immob,-23.3, fpi, Nnow, wtfNow);
    
    //fpi<1.0 ||
//    if((patch[0].ID==239202)&&(fpi<1.0 || ns_litr->litr1n != ns_litr->litr1n || ns_litr->litr2n != ns_litr->litr2n || ns_litr->litr3n != ns_litr->litr3n || ns_litr->litr4n != ns_litr->litr4n)) printf("actual decomp and uptake (after): patch %d, fpi %lf, carbon(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf), nflux(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf), immob(%lf,%lf,%lf,%lf), N(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf)\n",
//           patch[0].ID,fpi,
//           cs_litr->litr1c, cs_litr->litr2c, cs_litr->litr3c, cs_litr->litr4c, cs_soil->soil1c, cs_soil->soil2c, cs_soil->soil3c, cs_soil->soil4c,
//           ndf->pmnf_l1s1, ndf->pmnf_l2s2, ndf->pmnf_l3l2, ndf->pmnf_l4s3, ndf->pmnf_s1s2, ndf->pmnf_s2s3, ndf->pmnf_s3s4, ndf->pmnf_s4,
//           microbal_immoblitter, microbal_immob, microbal_mineralize, ndf->sminn_to_npool,
//           ns_litr->litr1n, ns_litr->litr2n, ns_litr->litr3n, ns_litr->litr4n, ns_soil->soil1n, ns_soil->soil2n, ns_soil->soil3n, ns_soil->soil4n);
//
    
    if(command_line[0].soilCNadaptation_falg == 1 ){
        if ((cs_soil->soil1c > 0.0) && (ns_soil->soil1n > 0.0)) patch[0].patch_SOIL1_CN = cs_soil->soil1c/ns_soil->soil1n;
        if ((cs_soil->soil2c > 0.0) && (ns_soil->soil2n > 0.0)) patch[0].patch_SOIL2_CN = cs_soil->soil2c/ns_soil->soil2n;
        if ((cs_soil->soil3c > 0.0) && (ns_soil->soil3n > 0.0)) patch[0].patch_SOIL3_CN = cs_soil->soil3c/ns_soil->soil3n;
        if ((cs_soil->soil4c > 0.0) && (ns_soil->soil4n > 0.0)) patch[0].patch_SOIL4_CN = cs_soil->soil4c/ns_soil->soil4n;
    }
    
    
//    total_preday_N = ns_litr->litr1n + ns_litr->litr2n +  ns_litr->litr3n
//    + ns_litr->litr4n + ns_soil->soil1n + ns_soil->soil2n + ns_soil->soil3n
//    + ns_soil->soil4n + totalNH4 + totalNO3;
    
    
	ndf->net_mineralized = daily_net_nmin;
	total_N = ns_litr->litr1n + ns_litr->litr2n +  ns_litr->litr3n
		+ ns_litr->litr4n + ns_soil->soil1n + ns_soil->soil2n
		+ ns_soil->soil3n + ns_soil->soil4n + totalNH4 + totalNO3;
    
    // 1) litter C increased
    // 2) mineralization decreased
    // 3) initial immobilization is too big
    // 4) soil C is increased
    
	balance = total_preday_N  - total_N - ndf->sminn_to_npool; // is this true?
//    if ( fabs(balance) > 0.0)
//        printf("\n Decomp N doesn't balance by %lf X1000000 ", balance*1000000.0); /// this balance check is very wrong. should be using fabs
	
	return (!ok);
} /* end update_decomp.c */
