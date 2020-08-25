/*--------------------------------------------------------------*/
/* 								*/
/*								*/
/*		compute_farq_psn				*/
/*								*/
/*	NAME							*/
/*		compute_farq_psn				*/
/*								*/
/*								*/
/*	SYNOPSIS						*/
/*		int compute_farq_psn(				*/
/*								*/
/*	returns:						*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*		compute farquhar photosynthesis			*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*	from Peter Thornton, 3/8/95				*/ 
/*								*/
/*	Smitch 2001-10-24					*/
/*	added BGC C4 parameterization from BGC code		*/
/*	updated 2 constants and formulae to match BGC 4.1.1	*/
/*								*/
/*--------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

int compute_farq_psn( struct psnin_struct *in,
					 struct psnout_struct* out,
                     struct epconst_struct epc, int verbose)
{
/*-------------------------------------------------------------------
Farquhar photosynthesis routine

	call with verbose=1 for complete output of calculated PSN parameters
	call with verbose=0 for output of A only
	
	 The following variables are held in the indat struct:
	 
	  in->pa         (Pa) atmospheric pressure
	  in->co2        (ppm) atmospheric [CO2]
	  in->t          (deg C) air temperature
	  in->irad       (umol/m2/s) PAR photon flux density
      in->g          (mm/s * 10**3) conductance to CO2  //<<--- (potential) glmax_vapor * lai_stomatal_fraction * 1000 / 1.6;
                                                        //<<--- (actual) gs_sunlit / lai_sunlit * 1000 / 1.6;
                                                                in which gs_sunlit = compute_vascular_stratum_conductane()
                                                                { m_final = m_APAR * m_tavg * m_LWP * m_CO2 * m_tmin * m_vpd;
                                                                  stomatal_conductance = stomatal_conductance_max * m_final;
                                                                  stratum_conductance = max(stomatal_conductance, cuticular_cond) * LAI * stomatal_fraction; <<---- (* species)
	  in->Rd         (umol/m2/s) dark respiration rate = stratum[0].cdf.leaf_day_mr / (stratum[0].epv.proj_lai * zone[0].metv.dayl*12.011e-9);
	  in->lnc        (kg Nleaf/m2) leaf N concentration, area units <<---- (* species)
	  in->flnr       (kg NRub/kg Nleaf) fraction of leaf N in Rubisco <<---- (* species)

	  smitch added:

  	  in->c3	flag for C3 photosynthesis (alternative = C4)
	  
	   Of these, t and Rd are used multiple times, and so are copied into local
	   variables.  The others are referenced by their struct names.
	   
		The outdat structure is used to pass output from this routine to the
		calling program, and the verbose flag sets the level of output reporting.
		The complete list of outdat parameters is:
		
         out->g         (umol/m2/s/Pa) conductance to CO2  = in->g * 1e3 / (R * tk); tk = abs. temp.; R = 8.3143 gas law constant
		 out->O2        (Pa) atmospheric [O2]
		 out->Ca        (Pa) atmospheric [CO2]
		 out->Ci        (Pa) intercellular [CO2]
		 out->gamma     (Pa) CO2 compensation point, no Rd
		 out->Kc        (Pa) MM constant carboxylation
		 out->Ko        (Pa) MM constant oxygenation
		 out->act       (umol/kg/s) Rubisco activity
		 out->Vmax      (umol/m2/s) max rate carboxylation
		 out->Jmax      (umol/m2/s) max rate electron transport
		 out->J         (umol/m2/s) rate of RuBP regeneration
		 out->Av        (umol/m2/s) carboxylation limited assimilation
		 out->Aj        (umol/m2/s) RuBP regen limited assimilation
		 out->A         (umol/m2/s) final assimilation rate <------ assimilation rate
		 out->dC13	delta C13 discrimination
		 
	-----------------------------------------------------------------*/
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	
//    double fleaf_gr_perc =  epc.gr_perc/(1.0+epc.alloc_frootc_leafc+epc.alloc_stemc_leafc+epc.alloc_crootc_stemc*epc.alloc_stemc_leafc);
    // a fraction of epc.gr_perc !? the "A" is total psn_to_cpool; that's why we need to fraction the leafc
    // assuming, leafc grow resp. takes place at the leaf and others are at their corresponding body parts
    // then the measured "A" is just the GPP - leaf g. resp. - leaf m. resp.
    double Rd_grow = 0.0;
    
	int ok=1;
	double t;      /* (deg C) temperature */
	double tk;     /* (K) absolute temperature */
	double g;      /* (umol/m2/s/Pa) conductance to CO2 */
	double O2;     /* (Pa) atmospheric partial pressure O2 */
	double Ca;     /* (Pa) atmospheric partial pressure CO2 */
	double gamma;  /* (Pa) co2 compensation point, no dark respiration */ //<--- this is not light compensation point
	double Kc;     /* (Pa) MM constant for carboxylase reaction */
	double Ko;     /* (Pa) MM constant for oxygenase reaction */
	double act;    /* (umol/kgRubisco/s) Rubisco activity */
	double Rd;     /* (umol/m2/s) dark respiration rate */
	double Vmax;   /* (umol/m2/s) maximum carboxylation velocity */
	double Jmax;   /* (umol/m2/s) maximum rate of electron transport */
	double J;      /* (umol/m2/s) maximum rate of Rubisco regeneration */ //<<-------
	double Av;     /* (umol/m2/s) Rubisco limited assimilation rate */
	double Aj;     /* (umol/m2/s) RuBP regeneration limited assim rate */
	double A;      /* (umol/m2/s) net assimilation rate */
	double aa,bb,cc,det;
	/* New vbl introduced in BGC 4.1.1 for Jmax calc - ppe - smitch 2001 */
	double ppe;	/* (mol/mol) photons absorbed by PSII per e- xported */
	int c3;		/* c3 flag to match bgc */
	/*------------------------------------------------------------------
	the weight proportion of Rubisco to its nitrogen content, fnr, is
	calculated from the relative proportions of the basic amino acids
	that make up the enzyme, as listed in the Handbook of Biochemistry,
	Proteins, Vol III, p. 510, which references:
	Kuehn and McFadden, Biochemistry, 8:2403, 1969
	--------------------------------------------------------------*/
	static double fnr = 7.16;   /* kg Rub/kg NRub */
	/*-----------------------------------------------------------------
	the following constants are from:
	Woodrow, I.E., and J.A. Berry, 1980. Enzymatic regulation of photosynthetic
	CO2 fixation in C3 plants. Ann. Rev. Plant Physiol. Plant Mol. Biol.,
	39:533-594.
	Note that these values are given in the units used in the paper, and that
	they are converted to units appropriate to the rest of this function before
	they are used.
	----------------------------------------------------------------------*/
	/* Changing Kc and Ko to match changes made by Peter Thornton in BGC
		4.1.1, he cites de Pury and Farquharson (1997).  Simple scaling
		of photosynthesis from leaves to canopies. smitch 2001        */
	/* static double Kc25 = 270.0; */  
	static double Kc25 = 404.0;  /* (ubar) MM const carboxylase, 25 deg C */
	static double q10Kc = 2.1;    /* (DIM) Q_10 for kc */
	/* static double Ko25 = 400.0; */  
	static double Ko25 = 248.0;   /* (mbar) MM const oxygenase, 25 deg C */
	static double q10Ko = 1.2;    /* (DIM) Q_10 for ko */
	static double act25 = 3.6;    /* (umol/mgRubisco/min) Rubisco activity */
	static double q10act = 2.4;   /* (DIM) Q_10 for Rubisco activity */
	/* new constant used in calculating Jmax - smitch 2001 */
	static double pabs = 0.85;    /* (DIM) fPAR effectively absorbed by
					PSII */

    //int kk;
    //for(kk=0; kk<10; kk++){
        /* local variables  */
        Rd = in->Rd;// + Rd_grow;
        t = in->t;
        tk = t + 273.15;
        /* convert conductance from m/s * 10**3 --> umol/m2/s/Pa */
        g = in->g * 1e3 / (R * tk);
        /* convert atmospheric CO2 from ppm --> Pa */
        Ca = in->co2 * in->pa / 1e6;
        /* smitch - Add adjustment of  parameters for C3/C4 psn */
        c3 = 0;
        /* tague for Laura: FOR NOW WE WILL HARDCODE C3 */
        in->c3 = 1; //most trees are c3; known c4 is grass
        
        if (in->c3){
            ppe = 2.6;
            c3 = 1;
        } else {
            /* C4 */
            ppe = 3.5;
            Ca *= 10.0;
        }//if else

        /* calculate atmospheric O2 in Pa, assumes 21% O2 by volume */
        O2 = 0.21 * in->pa;
        /* correct kinetic constants for temperature, and do unit conversions */
        Ko = Ko25 * pow(q10Ko, (t-25.0)/10.0);
        Ko = Ko * 100.0;   /* mbar --> Pa */
        if (t > 15.0){
            Kc = Kc25 * pow(q10Kc, (t-25.0)/10.0);
            act = act25 * pow(q10act, (t-25.0)/10.0);
        }
        else{
            Kc = Kc25 * pow(1.8*q10Kc, (t-15.0)/10.0) / q10Kc;
            act = act25 * pow(1.8*q10act, (t-15.0)/10.0) / q10act;
        }
        Kc = Kc * 0.10;   /* ubar --> Pa */
        act = act * 1e6 / 60.0;     /* umol/mg/min --> umol/kg/s */
        /* calculate gamma (Pa), assumes Vomax/Vcmax = 0.21 */
        gamma = 0.5 * 0.21 * Kc * O2 / Ko;
        /* calculate Vmax from leaf nitrogen data and Rubisco activity */
        /*---------------------------------------------------------------
        kg Nleaf   kg NRub    kg Rub      umol            umol
        -------- X -------  X ------- X ---------   =   --------
        m2         kg Nleaf   kg NRub   kg RUB * s       m2 * s
        
         (lnc)  X  (flnr)  X  (fnr)  X   (act)     =    (Vmax)
        // leafN and flnr should be variable; deciduous default flnr is 0.08. then 0.08 of leafn is available for Rub?
        // wiki: RuBisCO is the most abundant protein in leaves, accounting for 50% of soluble leaf protein in C3 plants (20–30% of total leaf nitrogen) and 30% of soluble leaf protein in C4 plants (5–9% of total leaf nitrogen)
        // In cyanobacteria, inorganic phosphate (Pi) participates in the co-ordinated regulation of photosynthesis. Pi binds to the RuBisCO active site and to another site on the large chain where it can influence transitions between activated and less active conformations of the enzyme. Activation of bacterial RuBisCO might be particularly sensitive to Pi levels, which can act in the same way as RuBisCO activase in higher plants.
         
        ---------------------------------------------------------------*/
        Vmax = in->lnc * in->flnr * fnr * act;
        /*-------------------------------------------------------------
        calculate Jmax = f(Vmax), reference:
        Wullschleger, S.D., 1993.  Biochemical limitations to carbon assimilation
        in C3 plants - A retrospective analysis of the A/Ci curves from
        109 species. Journal of Experimental Botany, 44:907:920.
        --------------------------------------------------------------*/
        /* smitch 2001 - PET reports in BGC 4.1.1 that he had to change
        this equation to get quantum yields right.
        Old:
        Jmax = 29.1 + 1.64*Vmax;
        New:    */
        Jmax = 2.1 * Vmax;

        /*-------------------------------------------------------------
        calculate J = f(Jmax, I), reference:
        Farquhar, G.D., and S. von Caemmerer, 1982.  Modelling of photosynthetic
        response to environmental conditions.  In Encyclopedia of Plant
        Physiology, New Series, Vol. 12B, Physiological Plant Ecology II,
        O.L. Lange, P.S. Nobel, C.B. Osmond, and H. Ziegler, eds, Springer-
        Verlag, Berlin, Germany, pp 549-587.
        --------------------------------------------------------------*/
        /* smitch 2001 - PET has replaced this also in BGC 4.1.1.  Commenting
        this version out, then follow with new.
        if (Jmax > 0.0)
            J = Jmax * in->irad / (in->irad + 2.1*Jmax);
        else
            J = 0.0;
        */
        /* New, substituting our vbl names (note I'm assuming that
            his ppfd in BGC is same as irad in this code.
            (also change a->aa, b->bb, c->cc, rename psn
            structure to in->)  smitch 2001  */
            /* calculate J = f(Jmax, ppfd), reference:
            de Pury and Farquhar 1997
            Plant Cell and Env.
            */
            aa = 0.7;
            bb = -Jmax - (in->irad*pabs/ppe);
            cc = Jmax * in->irad*pabs/ppe;
            J = (-bb - sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa);

        /*---------------------------------------------------------------
        solve for Av and Aj using the quadratic equation, substitution for Ci
        from A = g(Ca-Ci) into the equations from Farquhar and von Caemmerer:
         
             Vmax (Ci - gamma)
             Av =  -------------------   -   Rd
             Ci + Kc (1 + O2/Ko)
         
         
               J (Ci - gamma)
               Aj  =  -------------------  -   Rd
               4.5 Ci + 10.5 gamma
        ---------------------------------------------------------------------*/
        /* quadratic solution for Av */
        aa = -1.0/g;
        bb = Ca + (Vmax - Rd)/g + Kc*(1.0 + O2/Ko);
        cc = Vmax*(gamma - Ca) + Rd*(Ca + Kc*(1.0 + O2/Ko));
        if ((det = bb*bb - 4.0*aa*cc) < 0.0){
            printf("negative root error in psn routine\n");
            return (1);
        }
        Av = (-bb + sqrt(det)) / (2.0*aa);
        /* quadratic solution for Aj */
        aa = -4.5/g;
        bb = 4.5*Ca + 10.5*gamma + J/g - 4.5*Rd/g;
        cc = J*(gamma - Ca) + Rd*(4.5*Ca + 10.5*gamma);
        if ((det = bb*bb - 4.0*aa*cc) < 0.0){
            printf("negative root error in psn routine\n");
            return(1);
        }
        Aj = (-bb + sqrt(det)) / (2.0*aa);
        /* calculate A as the minimum of (Av,Aj) */
        if (Av < Aj)
            A = Av;
        else
            A = Aj;
        /* primary output */
        out->A = A;// + Rd_grow; //(important note!) A is net (original); RHESSys cdf.psn_to_cpool is GPP; leaf_mr is added but not leaf_gr. so adding here
        out->g = g;
        out->Ca = Ca;
        out->Ci = Ca - A/g;
        out->dC13 = out->Ci/out->Ca*(28-4.4)+4.4;
        out->gamma = gamma;
        out->O2 = O2;
        out->Kc = Kc;
        out->Ko = Ko;
        out->act = act;
        out->Vmax = Vmax;
        out->Jmax = Jmax;
        out->J = J;
        out->Av = Av;
        out->Aj = Aj;
        
        //Rd_grow = A * fleaf_gr_perc;
        //Rd_grow = A * epc.gr_perc;
   // }//kk for loop .... this is for putting grow respiration back into dark respiration.
    
    
    
	/*
		printf("\n %lf %lf %lf %lf %lf %lf", A, out->Ci, out->Vmax, out->Jmax, out->dC13, out->Ca); 
	*/
/*
                 printf("psnin: %lf %lf %lf %lf %lf %lf %lf %lf %d %lf\n",
                                in->pa, in->co2, t, in->irad,
                                g, in->Rd, in->lnc, in->flnr,
                                c3, ppe);
                                printf("psnout: %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                                g, O2, Ca, Ca - A/g,
                                out->gamma, Kc, Ko, act,
                                Vmax, Jmax, J, Av,
                                Aj, A);   
*/
	return (!ok);
}	 /* end compute_farq_psn.c */

