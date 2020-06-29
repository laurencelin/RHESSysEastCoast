/*--------------------------------------------------------------*/
/* 								*/
/*		compute_potential_decomp					*/
/*								*/
/*								*/
/*	NAME							*/
/*	compute_potential_decomp -  					*/
/*		performs decomposition and updates soil/litter	*/
/*		carbon and nitrogen stores			*/
/*								*/
/*	SYNOPSIS						*/
/*	int compute_potential_decomp(					*/
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

int compute_potential_decomp(double tsoil,
                             double maxpsi,
							 double minpsi,
							 double theta, //<<----
							 double std,
							 struct  soil_c_object   *cs_soil,
							 struct  soil_n_object   *ns_soil,
							 struct  litter_c_object *cs_litr,
							 struct  litter_n_object *ns_litr,
							 struct cdayflux_patch_struct *cdf,
							 struct ndayflux_patch_struct *ndf,
                             struct patch_object *patch,
                             double soilDecayScalar,
                             int soilCNadaptation_falg,
                             int vegtype)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	int ok;
	double rate_scalar, t_scalar, w_scalar;
	double a,b,c,d;
	double tk, thetai;
	double rfl1s1, rfl2s2,rfl4s3,rfs1s2,rfs2s3,rfs3s4;
	double kl1_base,kl2_base,kl4_base,ks1_base,ks2_base,ks3_base,ks4_base;
	double kl1,kl2,kl4,ks1,ks2,ks3,ks4;
	double cn_l1,cn_l2,cn_l3,cn_l4,cn_s1,cn_s2,cn_s3,cn_s4;
    double cn_ls1,cn_ls2,cn_ls3,cn_ss4;
	double plitr1c_loss, plitr3c_loss, plitr2c_loss, plitr4c_loss;
	double psoil1c_loss, psoil2c_loss, psoil3c_loss, psoil4c_loss;
	double pmnf_l1s1,pmnf_l2s2,pmnf_l3l2, pmnf_l4s3,pmnf_s1s2,pmnf_s2s3,pmnf_s3s4,pmnf_s4;
	double potential_immob,potential_immoblitter,mineralized;
	double weight1, weight2, theta1, theta2;
    double w_scalar2;
	int nlimit, i;
	#define NUM_NORMAL  10 	/* resolution of normal distribution */
	double NORMAL[10]= {0,0,0.253,0.524,0.842,1.283,-0.253,-0.524,-0.842,-1.283};

	ok = 0;
	/* calculate the rate constant scalar for soil temperature,
	assuming that the base rate constants are assigned for non-moisture
	limiting conditions at 25 C. The function used here is taken from
	Lloyd, J., and J.A. Taylor, 1994. On the temperature dependence of
	soil respiration. Functional Ecology, 8:315-323.
	This equation is a modification of their eqn. 11, changing the base
	temperature from 10 C to 25 C, since most of the microcosm studies
	used to get the base decomp rates were controlled at 25 C. */
	if (tsoil < -10.0){
		/* no decomp processes for tsoil < -10.0 C */
		t_scalar = 0.0;
	}
	else{
		tk = tsoil + 273.15;
		t_scalar = exp(308.56*((1.0/71.02)-(1.0/(tk-227.13))));
	}
	/* calculate the rate constant scalar for soil water content.
	use same empirical function as control on nitrification from NGas model
	but set parameters to reduce sensitivity to water stress
	(Parton et al, 1996 Global Biogeochemical cycles, 10:3, 401-412 ) */

	a=0.68; b=2.5; c=0.0012; d=2.84;
	w_scalar = 0.0;
	if (std > 0.0) {
		for (i=0; i<NUM_NORMAL; i++) {
			thetai = theta + NORMAL[i]*std;
			thetai = min(1.0, thetai);
			thetai = max(0.1, thetai);
			w_scalar  += exp(d*(b-a)/(a-c)*log((thetai -b)/(a-b))) * exp(d*log((thetai-c)/(a-c)));
			}//end of for i
        w_scalar *= 1.0/NUM_NORMAL;
	} else {
		if ((theta <= 0.0) || (theta > 1.0))
			theta = 1.0;
		if (theta  > c) {
			w_scalar  = exp(d*(b-a)/(a-c)*log((theta -b)/(a-b))) * exp(d*log((theta-c)/(a-c)));
        }else{
            w_scalar = 0.000001;
        }//end of ifelse
	}//end of ifelse
    if (w_scalar > 0.0)
        w_scalar = sqrt(w_scalar);
    else
        w_scalar = 0.0;
    
    
    if (patch[0].soil_defaults[0][0].soil_type.sand > 0.5) {
        a = 1.56; b=12.0; c=16.0; d=2.01;
    } else if (patch[0].soil_defaults[0][0].soil_type.clay > 0.5) {
        a = 60.0; b=18.0; c=22.0; d=1.06;
    } else {
        a=4.82; b=14.0; c=16.0; d=1.39;
    }
    w_scalar2 = 0.0;
    if (std > 0) {
        for (i = 1; i< NUM_NORMAL; i++) {
            thetai = theta + std*NORMAL[i];
            thetai = min(1.0, thetai);
            thetai = max(0.002, thetai);
            w_scalar2 += min( a*exp(-c*exp(-d*thetai*log(b))*log(b)), 1.0);;
        }// for i
        w_scalar2 *= 1.0/NUM_NORMAL;
    } else {
        //water_scalar = min(1.0, a / pow(b,  (c / pow(b, (d*theta) )) ) );
        w_scalar2 = min( a*exp(-c*exp(-d*theta*log(b))*log(b)), 1.0);
    }//if

	rate_scalar = w_scalar * t_scalar;
	/* assign output variables */
	cs_litr->t_scalar = t_scalar;
	cs_litr->w_scalar = w_scalar;
	/* calculate compartment C:N ratios */
	if ((cs_litr->litr1c > 0.0) && (ns_litr->litr1n > 0.0))	cn_l1 = cs_litr->litr1c/ns_litr->litr1n; else cn_l1 = 0.0;//LIVELAB_CN;
    if ((cs_litr->litr2c > 0.0) && (ns_litr->litr2n > 0.0))	cn_l2 = cs_litr->litr2c/ns_litr->litr2n; else cn_l2 = 0.0; //CEL_CN;
    if ((cs_litr->litr3c > 0.0) && (ns_litr->litr3n > 0.0)) cn_l3 = cs_litr->litr3c/ns_litr->litr3n; else cn_l3 = 0.0; //CEL_CN;
    if ((cs_litr->litr4c > 0.0) && (ns_litr->litr4n > 0.0))	cn_l4 = cs_litr->litr4c/ns_litr->litr4n; else cn_l4 = 0.0; // LIG_CN;
	
    if( cn_l1 < 5.0 && (cs_litr->litr1c > 1e-10) && (ns_litr->litr1n > 1e-10)) printf("l1CN@%d = %e, %e, %e\n",patch[0].ID,cn_l1, cs_litr->litr1c, ns_litr->litr1n);
    if( cn_l2 < 5.0 && (cs_litr->litr2c > 1e-10) && (ns_litr->litr2n > 1e-10)) printf("l2CN@%d = %e, %e, %e\n",patch[0].ID,cn_l2, cs_litr->litr2c, ns_litr->litr2n);
    if( cn_l3 < 5.0 && (cs_litr->litr3c > 1e-10) && (ns_litr->litr3n > 1e-10)) printf("l3CN@%d = %e, %e, %e\n",patch[0].ID,cn_l3, cs_litr->litr3c, ns_litr->litr3n);
    if( cn_l4 < 5.0 && (cs_litr->litr4c > 1e-10) && (ns_litr->litr4n > 1e-10)) printf("l4CN@%d = %e, %e, %e\n",patch[0].ID,cn_l4, cs_litr->litr4c, ns_litr->litr4n);
    
    if(soilCNadaptation_falg == 1 ){
        cn_ls1 = patch[0].patch_liter1_soil1_ratio * cn_l1;
        cn_ls2 = patch[0].patch_liter2_soil2_ratio * cn_l2;
        cn_ls3 = patch[0].patch_liter4_soil3_ratio * cn_l4;
        cn_ss4 = patch[0].patch_SOIL4_CN; //patch[0].patch_soil3_soil4_ratio * patch[0].patch_SOIL3_CN;//(fix soil4CN for upper 7 testing)
        //double cn_ls1,cn_ls2,cn_ls3,cn_ss4;
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
	
	/* calculate the corrected rate constants from the rate scalar and their
	base values. All rate constants are (1/day) */
    double adjust_rate = 0.45;
    //adjust_rate *= (patch[0].drainage_type>0 && patch[0].drainage_type % actionRIPARIAN==0? 1.7:1.0);
	kl1_base = 0.7 * adjust_rate;      /* labile litter pool */
	kl2_base = 0.07 * adjust_rate;     /* cellulose litter pool */
	kl4_base = 0.014 * adjust_rate;     /* lignin litter pool */
	ks1_base = 0.07 * adjust_rate * soilDecayScalar;    /* fast microbial recycling pool */
	ks2_base = 0.014 * adjust_rate * soilDecayScalar;    /* medium microbial recycling pool */
	ks3_base = 0.0014 * adjust_rate * soilDecayScalar;   /* slow microbial recycling pool */
    ks4_base = 0.0001 * adjust_rate * soilDecayScalar;   /* recalcitrant SOM (humus) pool */ // scale for wetland (new developing feature)
	
    //(vegtype>0? (patch[0].Ksat_vertical+(1.0-patch[0].Ksat_vertical)*0.5) : 0.0);
    kl1 = kl1_base * rate_scalar * patch[0].Ksat_vertical; //(vegtype>0? 1.0 : 0.01);
	kl2 = kl2_base * rate_scalar * patch[0].Ksat_vertical; //(vegtype>0? 1.0 : 0.01);
	kl4 = kl4_base * rate_scalar * patch[0].Ksat_vertical; //(vegtype>0? 1.0 : 0.01);
	ks1 = ks1_base * rate_scalar * patch[0].Ksat_vertical; //(vegtype>0? 1.0 : 0.01);
	ks2 = ks2_base * rate_scalar * patch[0].Ksat_vertical; //(vegtype>0? 1.0 : 0.01);
	ks3 = ks3_base * rate_scalar * patch[0].Ksat_vertical; //(vegtype>0? 1.0 : 0.01);
	ks4 = ks4_base * rate_scalar * patch[0].Ksat_vertical; //(vegtype>0? 1.0 : 0.01);
	
    //printf("%f,%f,%f\n",cn_l1,cn_l2,cn_l4);
//    printf("potential decomp: %lf, (%lf,%lf,%lf,%lf, %lf,%lf,%lf,%lf), (%lf,%lf,%lf,%lf, %lf,%lf,%lf,%lf)\n", rate_scalar,
//           psoil1c_loss, psoil2c_loss, psoil3c_loss, psoil4c_loss,
//           psoil1c_loss, psoil2c_loss, psoil3c_loss, psoil4c_loss,
//           pmnf_l1s1, pmnf_l2s2, pmnf_l3l2, pmnf_l4s3,
//           pmnf_s1s2, pmnf_s2s3, pmnf_s3s4, pmnf_s4
//           );
//
    
	/* initialize the potential loss and mineral N flux variables */
	plitr1c_loss = plitr2c_loss = plitr3c_loss = plitr4c_loss = 0.0;
	psoil1c_loss = psoil2c_loss = psoil3c_loss = psoil4c_loss = 0.0;
	pmnf_l1s1 = pmnf_l2s2 = pmnf_l3l2 = pmnf_l4s3 = 0.0;
	pmnf_s1s2 = pmnf_s2s3 = pmnf_s3s4 = pmnf_s4 = 0.0;
	
	/* calculate the non-nitrogen limited fluxes between litter and
	soil compartments. These will be ammended for N limitation if it turns
	out the potential gross immobilization is greater than potential gross
	mineralization. */
	/* 1. labile litter to fast microbial recycling pool */
	if ((cs_litr->litr1c > 0.0) && (ns_litr->litr1n > 0.0)) {
		plitr1c_loss = kl1 * cs_litr->litr1c;
        pmnf_l1s1 = plitr1c_loss * (soilCNadaptation_falg==1? (1.0 - rfl1s1 - patch[0].patch_liter1_soil1_ratio)/cn_ls1 : (1.0 - rfl1s1 - (cn_s1/cn_l1))/cn_s1); // check forumla (correct)
        //pmnf_l1s1 is positive when immobilizing
	}
	/* 2. cellulose litter to medium microbial recycling pool */
	if ((ns_litr->litr2n > 0.0) && (cs_litr->litr2c > 0.0)) {
		plitr2c_loss = kl2 * cs_litr->litr2c;
        pmnf_l2s2 = plitr2c_loss * (soilCNadaptation_falg==1? (1.0 - rfl2s2 - patch[0].patch_liter2_soil2_ratio)/cn_ls2 : (1.0 - rfl2s2 - (cn_s2/cn_l2))/cn_s2);
	}
	/* 2b. shield cellulose litter to goes to cellulose litter pool */
	/* respiration fractions not available to assume the same as for lignan (the "shield") */
	if ((ns_litr->litr3n > 0.0) && (cs_litr->litr3c > 0.0)) {
		plitr3c_loss = kl4 * cs_litr->litr3c;
		if ((ns_litr->litr2n > 0.0) && (cs_litr->litr2c > 0.0)) pmnf_l3l2 = plitr3c_loss * (1.0 - rfl4s3 - (cn_l2/cn_l3))/cn_l2;
        else pmnf_l3l2 = -plitr3c_loss * rfl4s3/cn_l3; //<<-----
	}
	/* 3. lignin litter to slow microbial recycling pool */
	if ((ns_litr->litr4n > 0.0) && (cs_litr->litr4c > 0.0)) {
		plitr4c_loss = kl4 * cs_litr->litr4c;
        pmnf_l4s3 = plitr4c_loss * (soilCNadaptation_falg==1? (1.0 - rfl4s3 - patch[0].patch_liter4_soil3_ratio)/cn_ls3 : (1.0 - rfl4s3 - (cn_s3/cn_l4))/cn_s3);
	}
	/* 4. fast microbial recycling pool to medium microbial recycling pool */
	if ((ns_soil->soil1n > 0.0) && (cs_soil->soil1c > 0.0)) {
		psoil1c_loss = ks1 * cs_soil->soil1c;
		pmnf_s1s2 = psoil1c_loss * (1.0 - rfs1s2 - (cn_s2/cn_s1))/cn_s2;
	}
	/* 5. medium microbial recycling pool to slow microbial recycling pool */
	if ((ns_soil->soil2n > 0.0) && (cs_soil->soil2c > 0.0)) {
		psoil2c_loss = ks2 * cs_soil->soil2c;
		pmnf_s2s3 = psoil2c_loss * (1.0 - rfs2s3 - (cn_s3/cn_s2))/cn_s3;
	}
	/* 6. slow microbial recycling pool to recalcitrant SOM pool */
	if ((ns_soil->soil3n > 0.0) && (cs_soil->soil3c > 0.0)) {
		psoil3c_loss = ks3 * cs_soil->soil3c;
        pmnf_s3s4 = psoil3c_loss * (soilCNadaptation_falg==1? (1.0 - rfs3s4 - cn_ss4/cn_s3)/cn_ss4 : (1.0 - rfs3s4 - (cn_s4/cn_s3))/cn_s4);// try to fix soil4CN for upper 7 test
	}
	/* 7. mineralization of recalcitrant SOM */
	if ((ns_soil->soil4n > 0.0) && (cs_soil->soil4c > 0.0)) {
		psoil4c_loss = ks4 * cs_soil->soil4c;
		pmnf_s4 = -psoil4c_loss/cn_ss4;
	}
    //printf("%f,%f,%f,%f,%f,%f,%f,%f\n", pmnf_l1s1,pmnf_l2s2, pmnf_l3l2,pmnf_l4s3,  pmnf_s1s2,pmnf_s2s3,  pmnf_s3s4,pmnf_s4);
    
	/* determine if there is sufficient mineral N to support potential
	immobilization. Immobilization fluxes are positive, mineralization fluxes
	are negative */
	nlimit = 0;
	potential_immob = 0.0;
    potential_immoblitter = 0.0;
	mineralized = 0.0;
    if (pmnf_l1s1 > 0.0) potential_immoblitter += pmnf_l1s1; //potential_immob += pmnf_l1s1; //pmnf_l1s1 is positive when immobilizing
	else mineralized += -pmnf_l1s1;
	
	if (pmnf_l2s2 > 0.0) potential_immoblitter += pmnf_l2s2; //potential_immob += pmnf_l2s2;
	else mineralized += -pmnf_l2s2;

	if (pmnf_l3l2 > 0.0) potential_immoblitter += pmnf_l3l2; //potential_immob += pmnf_l3l2;
	else mineralized += -pmnf_l3l2;
	
	if (pmnf_l4s3 > 0.0) potential_immoblitter += pmnf_l4s3; //potential_immob += pmnf_l4s3;
	else mineralized += -pmnf_l4s3;
	
	if (pmnf_s1s2 > 0.0) potential_immob += pmnf_s1s2;
	else mineralized += -pmnf_s1s2;
	
	if (pmnf_s2s3 > 0.0) potential_immob += pmnf_s2s3;
	else mineralized += -pmnf_s2s3;
	
	if (pmnf_s3s4 > 0.0) potential_immob += pmnf_s3s4;
	else mineralized += -pmnf_s3s4;
    
	mineralized += -pmnf_s4;
    
    if(pmnf_s4>0){printf("potential_decomp error here, %f,%f,%f\n", pmnf_s4, ks4, cs_soil->soil4c);}
    
	/* save the potential fluxes until plant demand has been assessed,
	to allow competition between immobilization fluxes and plant growth
	demands */
	ndf->mineralized = mineralized;
    ndf->potential_immob = potential_immob + potential_immoblitter; //<<---- it is a positive number, right?
    ndf->potential_immoblitter = potential_immoblitter; //<<-------- it is a positive number
    //ndf->potential_immob = 0.0; //<<---- disable immob
    //ndf->potential_immoblitter = 0.0; //<<------- disable immob
    cdf->plitr1c_loss = plitr1c_loss;
	ndf->pmnf_l1s1 = pmnf_l1s1;
	cdf->plitr2c_loss = plitr2c_loss;
	ndf->pmnf_l2s2 = pmnf_l2s2;
	cdf->plitr3c_loss = plitr3c_loss;
	ndf->pmnf_l3l2 = pmnf_l3l2;
	cdf->plitr4c_loss = plitr4c_loss;
	ndf->pmnf_l4s3 = pmnf_l4s3;
	cdf->psoil1c_loss = psoil1c_loss;
	ndf->pmnf_s1s2 = pmnf_s1s2;
	cdf->psoil2c_loss = psoil2c_loss;
	ndf->pmnf_s2s3 = pmnf_s2s3;
	cdf->psoil3c_loss = psoil3c_loss;
	ndf->pmnf_s3s4 = pmnf_s3s4;
	cdf->psoil4c_loss = psoil4c_loss;
	ndf->pmnf_s4 = pmnf_s4;
	cdf->kl4 = kl4;
    cdf->ks43 = (ks3 + ks4)*0.5;
//    if(patch[0].ID==239202) printf("potential decomp: %lf, litr(%e,%e,%e,%e)->(%e,%e,%e,%e), soil(%e,%e,%e,%e)->(%e,%e,%e,%e), litrCN[%e,%e,%e,%e]->[%e,%e,%e,%e], nflux(%e,%e,%e,%e, %e,%e,%e,%e)\n",
//           rate_scalar,
//           cs_litr->litr1c,cs_litr->litr2c,cs_litr->litr3c,cs_litr->litr4c, plitr1c_loss, plitr2c_loss, plitr3c_loss, plitr4c_loss,
//           cs_soil->soil1c, cs_soil->soil2c, cs_soil->soil3c, cs_soil->soil4c, psoil1c_loss, psoil2c_loss, psoil3c_loss, psoil4c_loss,
//           cn_l1,cn_l2,cn_l3,cn_l4,  cn_ls1,cn_ls2,cn_ls3,cn_ss4,
//           pmnf_l1s1, pmnf_l2s2, pmnf_l3l2, pmnf_l4s3, pmnf_s1s2, pmnf_s2s3, pmnf_s3s4, pmnf_s4
//           );
//    if(patch[0].ID==239202) printf("potential decomp: %lf, nflux( %e, %e, %e ,%e, %e, %e, %e, %e, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf )\n",
//                                   rate_scalar,
//                                   pmnf_l1s1, pmnf_l2s2, pmnf_l3l2, pmnf_l4s3, pmnf_s1s2, pmnf_s2s3, pmnf_s3s4, pmnf_s4,
//                                   cn_l1, cn_l2, cn_l3, cn_l4, cn_s1, cn_s2, cn_s3, cn_s4
//                                   );

    
    
    
	return(ok);
} /* end compute_potential_decomp.c */

