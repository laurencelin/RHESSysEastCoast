/*--------------------------------------------------------------*/
/* 											*/
/*					update_drainage_stream			*/
/*											*/
/*	update_drainage_stream.c - creates a patch object				*/
/*											*/
/*	NAME										*/
/*	update_drainage_stream.c - creates a patch object				*/
/*											*/
/*	SYNOPSIS									*/
/*	void update_drainage_stream( 							*/
/*					struct patch_object *patch			*/
/*				 			double,			 	*/
/*				 			double,			 	*/
/*				 			double,			 	*/
/*							int,				*/
/*							int)				*/
/*											*/
/* 											*/
/*											*/
/*	OPTIONS										*/
/*											*/
/*											*/
/*	DESCRIPTION									*/
/*											*/
/*											*/
/*											*/
/*											*/
/*	PROGRAMMER NOTES								*/
/*											*/
/*											*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"


void  update_drainage_stream(
								 struct patch_object *patch,
								 struct command_line_object *command_line,
								 double time_int,
								 int verbose_flag)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.				*/
	/*--------------------------------------------------------------*/
	double  compute_return_flow(
		int,
		double  ,
		double  );
	
	double  compute_delta_water(
		int,
		double,
		double,
		double,
		double,
		double);
	
    double compute_N_leached(
                             int verbose_flag,
                             double total_nitrate,
                             double Qout,
                             double N_decay_rate,
                             double activedepthz,
                             double N_absorption_rate,
                             int signal,
                             struct patch_object *patch);
	
	double compute_varbased_returnflow(
		double,
		double,
		double,
		struct litter_object *);


	double compute_varbased_flow(
		int,
		double,
		double,
		double,
		double,
		double *,
		struct patch_object *patch);

	double recompute_gamma(	
		struct patch_object *,
		double);
	/*--------------------------------------------------------------*/ 
	/*	Local variable definition.				*/ 
	/*--------------------------------------------------------------*/ 
	int i, j,k, d; 
	double m, Ksat; 
	double return_flow;  /* m */ 
	double NO3_leached_total, NO3_leached_to_stream; /* kg/m2 */ 
	double NH4_leached_total, NH4_leached_to_stream; /* kg/m2 */ 
	double DON_leached_total, DON_leached_to_stream; /* kg/m2 */ 
	double DOC_leached_total, DOC_leached_to_stream; /* kg/m2 */ 
	double patch_int_depth;  /* m of H2O */
	double route_to_stream; /* m3 */
	double route_to_surface;
	double Qin, Qout,Qstr_total;  /* m */
	double gamma, total_gamma, percent_tobe_routed;
	double Nin, Nout;  /* kg/m2 */
	double t1,t2,t3;
	double extrawater;
    
	d=0;
	route_to_stream = 0.0;
	return_flow=0.0;
    extrawater = 0.0;
	NO3_leached_to_stream = 0.0;
	NH4_leached_to_stream = 0.0;
	DON_leached_to_stream = 0.0;
	DOC_leached_to_stream = 0.0;
    
    
    double constantHold_no3, constantHold_nh4, constantHold_dom;
    double z1 = patch[0].sat_deficit_z>0? patch[0].sat_deficit_z : 0.0;
    //double p_decayRate = 1.0/patch[0].soil_defaults[0][0].porosity_decay;
    //double constantHold_SAT_proj_S0 = (1.0 - exp(-p_decayRate*activedepthz)) / (exp(-p_decayRate*z1)-exp(-p_decayRate*activedepthz));
    
    double N_decay_rate = (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NO3decayRate : patch[0].soil_defaults[0][0].N_decay_rate);//lookup the term that has been using.
    double activedepthz = (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z));
    constantHold_no3 = (exp(-N_decay_rate*z1)-exp(-N_decay_rate*activedepthz)) / (1.0 - exp(-N_decay_rate*activedepthz)) * 0.04166667;
    patch[0].sat_NO3 += patch[0].soil_ns.nitrate * constantHold_no3;
    patch[0].soil_ns.nitrate *= 1.0 - constantHold_no3;
    //double proj_total_sat_NO3 = patch[0].sat_NO3 * constantHold_SAT_proj_S0;
    
    //N_decay_rate = (command_line[0].NH4root2active>0.0? patch[0].soil_defaults[0][0].N_decay_rate : (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NH4decayRate : patch[0].soil_defaults[0][0].N_decay_rate));
    activedepthz = (command_line[0].NH4root2active>0.0? patch[0].rootzone.depth * command_line[0].NH4root2active : (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active>0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)));
    constantHold_nh4 = (exp(-N_decay_rate*z1)-exp(-N_decay_rate*activedepthz)) / (1.0 - exp(-N_decay_rate*activedepthz)) * 0.04166667;
    patch[0].sat_NH4 += patch[0].soil_ns.sminn * constantHold_nh4;
    patch[0].soil_ns.sminn *= 1.0 - constantHold_nh4;
    //double proj_total_sat_NH4 = patch[0].sat_NH4 * constantHold_SAT_proj_S0;
    
    //N_decay_rate = (command_line[0].rootNdecayRate > 0? patch[0].rootzone.DOMdecayRate : patch[0].soil_defaults[0][0].DOM_decay_rate);
    activedepthz = (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z));
    constantHold_dom = (exp(-N_decay_rate*z1)-exp(-N_decay_rate*activedepthz)) / (1.0 - exp(-N_decay_rate*activedepthz)) * 0.04166667;
    patch[0].sat_DOC += patch[0].soil_cs.DOC * constantHold_dom;
    patch[0].sat_DON += patch[0].soil_ns.DON * constantHold_dom;
    patch[0].soil_cs.DOC *= 1.0 - constantHold_dom;
    patch[0].soil_ns.DON *= 1.0 - constantHold_dom;
    //double proj_total_sat_DOC = patch[0].sat_DOC * constantHold_SAT_proj_S0;
    //double proj_total_sat_DON = patch[0].sat_DON * constantHold_SAT_proj_S0;
    
    


	/*--------------------------------------------------------------*/
	/*	for now there should be no recomputing of gamma for 	*/
	/*	streams because they do not route water to downslope	*/
	/*	neighbours						*/
    /*    calculate amuount of water output to stream as baseflow */
	/*--------------------------------------------------------------*/
    d=0; total_gamma =  patch[0].innundation_list[d].gamma; //total_gamma = recompute_gamma(patch, patch[0].innundation_list[d].gamma);
    if (total_gamma < ZERO ) {
		gamma = patch[0].soil_defaults[0][0].Ksat_0 * m * 2.0 * sqrt(patch[0].area) * time_int;
	} else {
		gamma = total_gamma * time_int;
	}

	route_to_stream = compute_varbased_flow(
		patch[0].num_soil_intervals,
		patch[0].std * command_line[0].std_scale,
		patch[0].sat_deficit,
		gamma,
		patch[0].soil_defaults[0][0].interval_size,
		patch[0].transmissivity_profile,
		patch);
	if (route_to_stream < 0.0) route_to_stream = 0.0;
    
    
    
    //Sept 18
    extrawater = patch[0].rz_storage+patch[0].unsat_storage - patch[0].sat_deficit - route_to_stream/patch[0].area;
    if(extrawater>0){
        // bounded by impervious surface
        // trigger returnflow and fully saturation
        return_flow = compute_varbased_returnflow(
                                                  patch[0].std * command_line[0].std_scale,
                                                  (patch[0].rz_storage+patch[0].unsat_storage+patch[0].constraintWaterTableTopDepth_def)*patch[0].Ksat_vertical,
                                                  (patch[0].sat_deficit+route_to_stream/patch[0].area)*patch[0].Ksat_vertical,
                                                  &(patch[0].litter));
        
        route_to_stream += (extrawater>return_flow? (extrawater-return_flow)*patch[0].area : 0.0);
        // route_to_stream should be bounded by impervious too. but this work is going to be later.
        
        patch[0].detention_store += return_flow;
        //patch[0].sat_deficit = 0.0; // wrong; because Qout will modify sat_deficit, outside of this function call, hourly
        //patch[0].sat_deficit += (return_flow - (patch[0].unsat_storage+patch[0].rz_storage));//original; i think it's wrong. it did not consider the ouside Qout-Qin hourly process.
        // patch[0].sat_deficit = -route_to_stream/patch[0].area;
        // patch[0].unsat_storage = 0.0;
        // patch[0].rz_storage = 0.0;
    }// extra water
		
	/*--------------------------------------------------------------*/
	/* compute Nitrogen leaching amount with baseflow		*/
	/*--------------------------------------------------------------*/
	if (command_line[0].grow_flag > 0) {

		NO3_leached_to_stream = compute_N_leached(
			verbose_flag,
			patch[0].sat_NO3, //patch[0].soil_ns.nitrate,
			route_to_stream / patch[0].area,
			(command_line[0].rootNdecayRate > 0? patch[0].rootzone.NO3decayRate : patch[0].soil_defaults[0][0].N_decay_rate),
			(command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)),
			patch[0].soil_defaults[0][0].NO3_adsorption_rate,
			5,patch);
		patch[0].soil_ns.NO3_Qout += NO3_leached_to_stream;


		NH4_leached_to_stream = compute_N_leached(
			verbose_flag,
			patch[0].sat_NH4, //patch[0].soil_ns.sminn,
			route_to_stream / patch[0].area,
			(command_line[0].NH4root2active>0.0? patch[0].soil_defaults[0][0].N_decay_rate : (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NH4decayRate : patch[0].soil_defaults[0][0].N_decay_rate)),
            (command_line[0].NH4root2active>0.0? patch[0].rootzone.depth * command_line[0].NH4root2active : (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active>0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z))),
			patch[0].soil_defaults[0][0].NH4_adsorption_rate,
			8,patch);
		patch[0].soil_ns.NH4_Qout += NH4_leached_to_stream;

		DON_leached_to_stream = compute_N_leached(
			verbose_flag,
			patch[0].sat_DON, //patch[0].soil_ns.DON,
			route_to_stream / patch[0].area,
			(command_line[0].rootNdecayRate > 0? patch[0].rootzone.DOMdecayRate : patch[0].soil_defaults[0][0].DOM_decay_rate),
			(command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)),
			patch[0].soil_defaults[0][0].DON_adsorption_rate,
			11,patch);
		patch[0].soil_ns.DON_Qout += DON_leached_to_stream;

		DOC_leached_to_stream = compute_N_leached(
			verbose_flag,
			patch[0].sat_DOC, //patch[0].soil_cs.DOC,
			route_to_stream / patch[0].area,
			(command_line[0].rootNdecayRate > 0? patch[0].rootzone.DOMdecayRate : patch[0].soil_defaults[0][0].DOM_decay_rate),
			(command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)),
			patch[0].soil_defaults[0][0].DOC_adsorption_rate,
			14,patch);
        
		patch[0].soil_cs.DOC_Qout += DOC_leached_to_stream;
		patch[0].streamflow_NO3 += NO3_leached_to_stream;
		patch[0].streamNO3_from_sub += NO3_leached_to_stream;
		patch[0].hourly[0].streamflow_NO3 += NO3_leached_to_stream;
		patch[0].hourly[0].streamflow_NO3_from_sub += NO3_leached_to_stream;

		patch[0].streamflow_NH4 += NH4_leached_to_stream;
		patch[0].streamflow_DON += DON_leached_to_stream;
		patch[0].streamflow_DOC += DOC_leached_to_stream;


	}//growth flag

	patch[0].Qout += (route_to_stream / patch[0].area);
	patch[0].base_flow += (route_to_stream / patch[0].area);
	patch[0].hourly_subsur2stream_flow += route_to_stream / patch[0].area;

    if(extrawater>0){
        // set this to be  "-route_to_stream/patch[0].area" and later += Qout = "route_to_stream/patch[0].area"
        //patch[0].sat_deficit = -route_to_stream/patch[0].area + patch[0].constraintWaterTableTopDepth_def;
        patch[0].sat_deficit -= patch[0].unsat_storage + patch[0].rz_storage - extrawater*patch[0].Ksat_vertical;
        patch[0].unsat_storage = 0.0; // converted to be part of sat
        patch[0].rz_storage = 0.0; // converted to be part of sat
    }// extra water


	/*--------------------------------------------------------------*/
	/*	calculate any return flow to the stream in this patch   */
	/*	and route any infiltration excess			*/
	/*--------------------------------------------------------------*/
    // move up
//    if ((patch[0].sat_deficit-patch[0].rz_storage-patch[0].unsat_storage) < -1.0*ZERO) {
//        return_flow = compute_varbased_returnflow(patch[0].std * command_line[0].std_scale,
//            patch[0].rz_storage+patch[0].unsat_storage,
//            patch[0].sat_deficit, &(patch[0].litter));
//        patch[0].detention_store += return_flow;
//        patch[0].sat_deficit += (return_flow - (patch[0].unsat_storage+patch[0].rz_storage));;
//        patch[0].unsat_storage = 0.0;
//        patch[0].rz_storage = 0.0;
//    }

	/*--------------------------------------------------------------*/
	/*	calculated any N-transport associated with return flow  */
	/*	- note available N reduced by what has already been lost  */
	/*	due to subsurface routing above				*/
	/* 	note only nitrate is assumed to follow return flow		*/
	/*--------------------------------------------------------------*/
	if (return_flow > ZERO) {
		Nout = compute_N_leached(
			verbose_flag,
			patch[0].sat_NO3, //patch[0].soil_ns.nitrate - NO3_leached_to_stream,
			return_flow,
			(command_line[0].rootNdecayRate > 0? patch[0].rootzone.NO3decayRate : patch[0].soil_defaults[0][0].N_decay_rate),
			(command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)),
			patch[0].soil_defaults[0][0].NO3_adsorption_rate,
			17, patch);
		patch[0].surface_NO3 += Nout;
		patch[0].soil_ns.NO3_Qout += Nout;

		Nout = compute_N_leached(
			verbose_flag,
			patch[0].sat_NH4, //patch[0].soil_ns.sminn - NH4_leached_to_stream,
			return_flow,
			(command_line[0].NH4root2active>0.0? patch[0].soil_defaults[0][0].N_decay_rate : (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NH4decayRate : patch[0].soil_defaults[0][0].N_decay_rate)),
            (command_line[0].NH4root2active>0.0? patch[0].rootzone.depth * command_line[0].NH4root2active : (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active>0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z))),
			patch[0].soil_defaults[0][0].NH4_adsorption_rate,
			20, patch);
		patch[0].surface_NH4 += Nout;
		patch[0].soil_ns.NH4_Qout += Nout;

		Nout = compute_N_leached(
			verbose_flag,
			patch[0].sat_DON, //patch[0].soil_ns.DON - DON_leached_to_stream,
			return_flow,
			(command_line[0].rootNdecayRate > 0? patch[0].rootzone.DOMdecayRate : patch[0].soil_defaults[0][0].DOM_decay_rate),
			(command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)),
			patch[0].soil_defaults[0][0].DON_adsorption_rate,
			23, patch);
		patch[0].surface_DON += Nout;
		patch[0].soil_ns.DON_Qout += Nout;

		Nout = compute_N_leached(
			verbose_flag,
			patch[0].sat_DOC, //patch[0].soil_cs.DOC - DOC_leached_to_stream,
			return_flow,
			(command_line[0].rootNdecayRate > 0? patch[0].rootzone.DOMdecayRate : patch[0].soil_defaults[0][0].DOM_decay_rate),
			(command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)),
			patch[0].soil_defaults[0][0].DOC_adsorption_rate,
			26, patch);
		patch[0].surface_DOC += Nout;
		patch[0].soil_cs.DOC_Qout += Nout;

	}// return flow
    
    patch[0].overland_flow += return_flow; //max(0.0, patch[0].detention_store - patch[0].soil_defaults[0][0].detention_store_size);
    //<<--- reset by compute_subsurface_routing.c
    
	/*--------------------------------------------------------------*/
	/*	route water and nitrogen lossed due to infiltration excess */
	/*	note we assume that this happens before return_flow losses */
	/*--------------------------------------------------------------*/

	if ( (patch[0].detention_store > patch[0].soil_defaults[0][0].detention_store_size) && (patch[0].detention_store > ZERO) ) {
        
        Qout = (patch[0].detention_store - patch[0].soil_defaults[0][0].detention_store_size);
		Nout = (min(1.0, Qout / patch[0].detention_store)) * patch[0].surface_NO3;
		patch[0].surface_NO3  -= Nout;
		patch[0].streamflow_NO3 += Nout;
		patch[0].hourly[0].streamflow_NO3 += Nout;
		patch[0].streamNO3_from_surface +=Nout;
		patch[0].hourly[0].streamflow_NO3_from_surface +=Nout;

		patch[0].surface_ns_leach += Nout;
		Nout = (min(1.0, Qout / patch[0].detention_store)) * patch[0].surface_DOC;
		patch[0].surface_DOC  -= Nout;
		patch[0].streamflow_DOC += Nout;
		Nout = (min(1.0, Qout / patch[0].detention_store)) * patch[0].surface_DON;
		patch[0].surface_DON  -= Nout;
		patch[0].streamflow_DON += Nout;
		Nout = (min(1.0, Qout / patch[0].detention_store)) * patch[0].surface_NH4;
		patch[0].surface_NH4  -= Nout;
		patch[0].streamflow_NH4 += Nout;
        
		patch[0].detention_store -= Qout;
		patch[0].return_flow += Qout; // <<----------------
		patch[0].hourly_sur2stream_flow += Qout;
		}

} /*end update_drainage_stream.c*/

