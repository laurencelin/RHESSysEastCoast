/*--------------------------------------------------------------*/
/* 											*/
/*					update_drainage_land			*/
/*											*/
/*	update_drainage_land.c - creates a patch object				*/
/*											*/
/*	NAME										*/
/*	update_drainage_land.c - creates a patch object				*/
/*											*/
/*	SYNOPSIS									*/
/*	void update_drainage_land( 							*/
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


void  update_drainage_land(
					struct patch_object *patch,
					 struct command_line_object *command_line,
					 double time_int,
					 int verbose_flag)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.				*/
	/*--------------------------------------------------------------*/
	double  compute_delta_water(
		int,
		double,
		double,
		double,
		double,
		double);
	
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
		struct patch_object *);


	double compute_N_leached(
         int verbose_flag,
         double total_nitrate,
         double Qout,
         double N_decay_rate,
         double activedepthz,
         double N_absorption_rate,
         int signal,
         struct patch_object *patch);
		
	
	double recompute_gamma(	
		struct patch_object *,
		double);


	double compute_infiltration( int,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double);
	
	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
	int j, d;
	double m, Ksat, std_scale;
	double NH4_leached_to_patch, NH4_leached_to_stream;
	double NO3_leached_to_patch, NO3_leached_to_stream;
	double DON_leached_to_patch, DON_leached_to_stream;
	double DOC_leached_to_patch, DOC_leached_to_stream;
	double NO3_leached_to_surface; /* kg/m2 */
	double NH4_leached_to_surface; /* kg/m2 */
	double DON_leached_to_surface; /* kg/m2 */
	double DOC_leached_to_surface; /* kg/m2 */
	double N_leached_total; /* kg/m2 */
	double DON_leached_total; /* kg/m2 */
	double DOC_leached_total; /* kg/m2 */
	double route_to_surface;  /* m3 */
	double return_flow,route_to_patch ;  /* m3 */
	double available_sat_water; /* m3 */
	double Qin, Qout;  /* m */
	double innundation_depth, infiltration; /* m */
	double total_gamma;
	double Nin, Nout; /* kg/m2 */ 
    double extrawater;

	struct patch_object *neigh;
	route_to_patch = 0.0;
	route_to_surface = 0.0;
	return_flow=0.0;
    extrawater = 0.0;
    
	DON_leached_to_patch = 0.0;
	DON_leached_to_stream = 0.0;
	DOC_leached_to_patch = 0.0;
	DOC_leached_to_stream = 0.0;
	NH4_leached_to_patch = 0.0;
	NH4_leached_to_stream = 0.0;
	NO3_leached_to_patch = 0.0;
	NO3_leached_to_stream = 0.0;
	NO3_leached_to_surface = 0.0;
	NH4_leached_to_surface = 0.0;
	DOC_leached_to_surface = 0.0;
	DON_leached_to_surface = 0.0;
	
	// by the way, editing on ipad. hope it's saved.
	// in the mother routine, "update_drainage_land.c" is called hourly for every patch (below). 
	// ideally, I wanna track hourly "leached_solute" and add to the neighbor "SAT_solute" (hourly).
	// at the same time, this "SAT_solute" is used for the local leaching calculation. 
	// by the local vertical movement of solute, local "SAT_solute" is supplied from the upper layer.
    // ** important note: sat_deficit_z is updated hourly to count for later IN/OUT from the outside function before calling this function
    double constantHold_no3, constantHold_nh4, constantHold_dom;
    double z1 = patch[0].sat_deficit_z>0? patch[0].sat_deficit_z : 0.0;
    //double p_decayRate = 1.0/patch[0].soil_defaults[0][0].porosity_decay;
    //double constantHold_SAT_proj_S0 = (1.0 - exp(-p_decayRate*activedepthz)) / (exp(-p_decayRate*z1)-exp(-p_decayRate*activedepthz));
    
    double N_decay_rate = (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NO3decayRate : patch[0].soil_defaults[0][0].N_decay_rate);//lookup the term that has been using.
	double activedepthz = (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z));
    if(activedepthz > z1){
        constantHold_no3 = (exp(-N_decay_rate*z1)-exp(-N_decay_rate*activedepthz)) / (1.0 - exp(-N_decay_rate*activedepthz)) * 0.04166667;
        patch[0].sat_NO3 += patch[0].soil_ns.nitrate * constantHold_no3; //adjust for hourly
        patch[0].soil_ns.nitrate *= 1.0 - constantHold_no3;
    }
    //double proj_total_sat_NO3 = patch[0].sat_NO3 * constantHold_SAT_proj_S0;
    
    N_decay_rate = (command_line[0].NH4root2active>0.0? patch[0].soil_defaults[0][0].N_decay_rate : (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NH4decayRate : patch[0].soil_defaults[0][0].N_decay_rate));
    activedepthz = (command_line[0].NH4root2active>0.0? patch[0].rootzone.depth * command_line[0].NH4root2active : (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active>0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)));
    if(activedepthz > z1){
        constantHold_nh4 = (exp(-N_decay_rate*z1)-exp(-N_decay_rate*activedepthz)) / (1.0 - exp(-N_decay_rate*activedepthz)) * 0.04166667;
        patch[0].sat_NH4 += patch[0].soil_ns.sminn * constantHold_nh4;
        patch[0].soil_ns.sminn *= 1.0 - constantHold_nh4;
    }
    //double proj_total_sat_NH4 = patch[0].sat_NH4 * constantHold_SAT_proj_S0;
    
    N_decay_rate = (command_line[0].rootNdecayRate > 0? patch[0].rootzone.DOMdecayRate : patch[0].soil_defaults[0][0].DOM_decay_rate);
    activedepthz = (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z));
    if(activedepthz > z1){
        constantHold_dom = (exp(-N_decay_rate*z1)-exp(-N_decay_rate*activedepthz)) / (1.0 - exp(-N_decay_rate*activedepthz)) * 0.04166667;
        patch[0].sat_DOC += patch[0].soil_cs.DOC * constantHold_dom;
        patch[0].sat_DON += patch[0].soil_ns.DON * constantHold_dom;
        patch[0].soil_cs.DOC *= 1.0 - constantHold_dom;
        patch[0].soil_ns.DON *= 1.0 - constantHold_dom;
    }
    //double proj_total_sat_DOC = patch[0].sat_DOC * constantHold_SAT_proj_S0;
    //double proj_total_sat_DON = patch[0].sat_DON * constantHold_SAT_proj_S0;
	// we need to know the vertical profile of NO3, not just the total amount!
	// "sat_solute" does not have the profile! 
	// solution 1a: we borrow decay from NO3 in whole soil column
	// solution 1b: we project a virtual sat_no3_0 for the profile to work.
	// ----------------- then we need to update the leaching function below.
	// change 1a: passing "sat_solute", N_decay_rate, and projected sat_no3_0 to the function
	// change 1b: passing "sat_deficit_z"; and we still need to determine the critial pt. for integration
	
	// after neighbor "sat_solute" is filled by the hourly step, we need to sort out resources in SAT to UNSAT, which should be done in "where?"
	
	/*--------------------------------------------------------------*/
	/*	recalculate gamma based on current saturation deficits  */
	/*      to account the effect of changes in water table slope 	*/
	/*--------------------------------------------------------------*/
    d=0; total_gamma = recompute_gamma(patch, patch[0].innundation_list[d].gamma);

	available_sat_water = max(((patch[0].soil_defaults[0][0].soil_water_cap
			- max(patch[0].sat_deficit,0.0))
			* patch[0].area),0.0);

	/*------------------------------------------------------------*/
	/*	calculate amuount of water output to OTHER patches			*/
	/*	this only computes subsurface flow, not overland flow	*/
	/*-----------------------------------------------------------*/

	std_scale = command_line[0].std_scale;

	route_to_patch = time_int * compute_varbased_flow(
		patch[0].num_soil_intervals,
		patch[0].std * std_scale, 
		patch[0].sat_deficit,
		total_gamma, 
		patch[0].soil_defaults[0][0].interval_size,
		patch[0].transmissivity_profile,
		patch);



	if (route_to_patch < 0.0) route_to_patch = 0.0;
	if ( route_to_patch > available_sat_water) 
		route_to_patch *= (available_sat_water)/(route_to_patch);
    
    //Sept 18
    extrawater = patch[0].rz_storage+patch[0].unsat_storage - patch[0].sat_deficit - route_to_patch/patch[0].area + patch[0].constraintWaterTableTopDepth_def;
    if(extrawater>0){
        // bounded by impervious surface
        // trigger returnflow and fully saturation
        return_flow = compute_varbased_returnflow(
                              patch[0].std * std_scale,
                              (patch[0].rz_storage+patch[0].unsat_storage+patch[0].constraintWaterTableTopDepth_def)*patch[0].Ksat_vertical,
                              (patch[0].sat_deficit+route_to_patch/patch[0].area)*patch[0].Ksat_vertical,
                              &(patch[0].litter));// mm
        
        //route_to_patch = max(route_to_patch,extrawater*patch[0].area*(1.0-patch[0].Ksat_vertical)); // volumn
        route_to_patch += (extrawater>return_flow? (extrawater-return_flow)*patch[0].area : 0.0);  // volumn
        
        patch[0].detention_store += return_flow;
        //patch[0].sat_deficit = 0.0; // wrong; because Qout will modify sat_deficit, outside of this function call, hourly
        //patch[0].sat_deficit += (return_flow - (patch[0].unsat_storage+patch[0].rz_storage));//original; i think it's wrong. it did not consider the ouside Qout-Qin hourly process.
        // patch[0].Qout is "route_to_patch" (hourly)
        // ---------- cannot do it here.
        // patch[0].sat_deficit = -route_to_patch/patch[0].area + return_flow; //
        // patch[0].unsat_storage = 0.0;
        // patch[0].rz_storage = 0.0;
    }// extra water
    
	/*--------------------------------------------------------------*/
	/* compute Nitrogen leaching amount				*/
	/*--------------------------------------------------------------*/
	if(command_line[0].grow_flag > 0) {
        
//        Nout = compute_N_leached(
//            verbose_flag,
//            patch[0].soil_ns.nitrate, //total solute
//            route_to_patch / patch[0].area,//Qout
//            (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NO3decayRate : patch[0].soil_defaults[0][0].N_decay_rate), // N_decay_rate
//            (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)), //activedepthz
//            patch[0].soil_defaults[0][0].NO3_adsorption_rate, // N_absorption_rate
//            0, //signal
//            patch); // patch
        
        Nout = compute_N_leached(
             verbose_flag,
             patch[0].sat_NO3, //total solute <--- projected total_sat_solute
             route_to_patch / patch[0].area,//Qout (already time_int adjusted)
             (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NO3decayRate : patch[0].soil_defaults[0][0].N_decay_rate), // N_decay_rate
             (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)), //activedepthz
             patch[0].soil_defaults[0][0].NO3_adsorption_rate, // N_absorption_rate
             5,patch); // signal, patch
		NO3_leached_to_patch = Nout * patch[0].area;
		patch[0].soil_ns.NO3_Qout += Nout;//<<-----------*********
        if(Nout<0 || Nout!=Nout || patch[0].soil_ns.nitrate<0 || patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || (Nout<=0&&patch[0].sat_NO3>0&&route_to_patch>0)) printf("update_drainage_land[%d,%e]: soil NO3 (%e) flux (%e)\n", patch[0].ID, route_to_patch, patch[0].soil_ns.nitrate, Nout);


		Nout = compute_N_leached(
			verbose_flag,
			patch[0].sat_NH4,
			route_to_patch / patch[0].area,
			(command_line[0].NH4root2active>0.0? patch[0].soil_defaults[0][0].N_decay_rate : (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NH4decayRate : patch[0].soil_defaults[0][0].N_decay_rate)),
            (command_line[0].NH4root2active>0.0? patch[0].rootzone.depth * command_line[0].NH4root2active : (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active>0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z))),
			patch[0].soil_defaults[0][0].NH4_adsorption_rate,
            8,patch);
		NH4_leached_to_patch = Nout * patch[0].area;
		patch[0].soil_ns.NH4_Qout += Nout;
        if(Nout<0 || Nout!=Nout ) printf("update_drainage_land[%d,%e]: soil NH4 (%e) flux (%e)\n", patch[0].ID,route_to_patch,patch[0].soil_ns.sminn, Nout);

		Nout = compute_N_leached(
			verbose_flag,
            patch[0].sat_DON,
			route_to_patch / patch[0].area,
			(command_line[0].rootNdecayRate > 0? patch[0].rootzone.DOMdecayRate : patch[0].soil_defaults[0][0].DOM_decay_rate),
			(command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)),
			patch[0].soil_defaults[0][0].DON_adsorption_rate,
			11,patch);
		DON_leached_to_patch = Nout * patch[0].area;
		patch[0].soil_ns.DON_Qout += Nout;
        if(Nout<0 || Nout!=Nout) printf("update_drainage_land[%d,%e]: soil DON (%e) flux (%e)\n", patch[0].ID,route_to_patch,patch[0].soil_ns.DON, Nout);

		Nout = compute_N_leached(
			verbose_flag,
			patch[0].sat_DOC,
			route_to_patch / patch[0].area,
			(command_line[0].rootNdecayRate > 0? patch[0].rootzone.DOMdecayRate : patch[0].soil_defaults[0][0].DOM_decay_rate),
			(command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)),
			patch[0].soil_defaults[0][0].DOC_adsorption_rate,
			14,patch);
		DOC_leached_to_patch = Nout * patch[0].area;
		patch[0].soil_cs.DOC_Qout += Nout;
        if(Nout<0 || Nout!=Nout) printf("update_drainage_land[%d,%e]: soil DOC (%e) flux (%e)\n", patch[0].ID,route_to_patch,patch[0].soil_cs.DOC, Nout);


	}//leaching from subsurface growth flag

	
	patch[0].Qout += (route_to_patch / patch[0].area); //<<------- what is time scale patch[0].Qout tracking? daily? this function is call hourly
    if(extrawater>0){
        // set this to be  "-route_to_patch/patch[0].area" and later += Qout = "route_to_patch/patch[0].area"
        //patch[0].sat_deficit = -route_to_patch/patch[0].area + patch[0].constraintWaterTableTopDepth_def;
        patch[0].sat_deficit -= patch[0].unsat_storage + patch[0].rz_storage - extrawater*patch[0].Ksat_vertical;
        patch[0].unsat_storage = 0.0; // converted to be part of sat
        patch[0].rz_storage = 0.0; // converted to be part of sat
    }// extra water

	/*--------------------------------------------------------------*/
	/*	calculate any return flow associated with this patch	*/
	/*	and route any infiltration excess			*/
	/*	return flow is flow leaving patch (i.e surface_Qout)  	*/
	/*	note that return flow that becomes detention storage   */
	/*	is added to surface_Qin					*/
	/*	similarly with associated nitrogen			*/
	/* 	note we move unsat_storage into saturated storage in this case */
	/*	saturated zone will be updated in compute_subsurface_routing	*/
	/*	i.e becomes part of Qout				*/
	/*--------------------------------------------------------------*/
	// calculation is moved up
   
    // below is the returnflow associated solute processes
	/*--------------------------------------------------------------*/
	/*	calculated any N-transport associated with return flow  */
	/*	-note available N reduced by what has already been 	*/
	/*	we assume that only nitrate follows return flow		*/
	/*	lost in subsurface flow routing				*/
	/*--------------------------------------------------------------*/
    // leaching from surface;
    // when return_flow>0, then extrawater = patch[0].rz_storage+patch[0].unsat_storage - patch[0].sat_deficit - route_to_patch/patch[0].area + patch[0].constraintWaterTableTopDepth_def >0
		if (command_line[0].grow_flag > 0 && return_flow>0) {
			Nout = compute_N_leached(
				verbose_flag,
                /// problem re-project!! && "soil_ns -= going2sat_NO3"
				patch[0].sat_NO3, //patch[0].soil_ns.nitrate - (NO3_leached_to_patch/patch[0].area),
				return_flow,
				(command_line[0].rootNdecayRate > 0? patch[0].rootzone.NO3decayRate : patch[0].soil_defaults[0][0].N_decay_rate),
				(command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)),
				patch[0].soil_defaults[0][0].NO3_adsorption_rate,
                17,patch);
			patch[0].surface_NO3 += Nout;
			patch[0].soil_ns.NO3_Qout += Nout;
            if(Nout<0 || Nout!=Nout || patch[0].soil_ns.nitrate!=patch[0].soil_ns.nitrate || patch[0].soil_ns.nitrate<0 || NO3_leached_to_patch<0 || NO3_leached_to_patch!=NO3_leached_to_patch) printf("update_drainage_land[%d,%e]: return NO3 (%e,%e) flux (%e)\n", patch[0].ID,return_flow,patch[0].soil_ns.nitrate - (NO3_leached_to_patch/patch[0].area),NO3_leached_to_patch, Nout);

			Nout = compute_N_leached(
				verbose_flag,
				patch[0].sat_NH4, //patch[0].soil_ns.sminn - (NH4_leached_to_patch/patch[0].area),
				return_flow,
				(command_line[0].NH4root2active>0.0? patch[0].soil_defaults[0][0].N_decay_rate : (command_line[0].rootNdecayRate > 0? patch[0].rootzone.NH4decayRate : patch[0].soil_defaults[0][0].N_decay_rate)),
                (command_line[0].NH4root2active>0.0? patch[0].rootzone.depth * command_line[0].NH4root2active : (command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active>0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z))),
				patch[0].soil_defaults[0][0].NH4_adsorption_rate,
				20,patch);
			patch[0].surface_NH4 += Nout;
			patch[0].soil_ns.NH4_Qout += Nout;
            if(Nout<0 || Nout!=Nout) printf("update_drainage_land[%d,%e]: return NH4 (%e) flux (%e)\n", patch[0].ID,return_flow,patch[0].soil_ns.sminn - (NH4_leached_to_patch/patch[0].area), Nout);

			Nout = compute_N_leached(
				verbose_flag,
				patch[0].sat_DON, //patch[0].soil_ns.DON - (DON_leached_to_patch/patch[0].area),
				return_flow,
				(command_line[0].rootNdecayRate > 0? patch[0].rootzone.DOMdecayRate : patch[0].soil_defaults[0][0].DOM_decay_rate),
				(command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)),
				patch[0].soil_defaults[0][0].DON_adsorption_rate,
				23,patch);
			patch[0].surface_DON += Nout;
			patch[0].soil_ns.DON_Qout += Nout;
            if(Nout<0 || Nout!=Nout) printf("update_drainage_land[%d,%e]: return DON (%e) flux (%e)\n", patch[0].ID,return_flow,patch[0].soil_ns.DON - (DON_leached_to_patch/patch[0].area), Nout);
            
			Nout = compute_N_leached(
				verbose_flag,
				patch[0].sat_DOC, //patch[0].soil_cs.DOC - (DOC_leached_to_patch/patch[0].area),
				return_flow,
				(command_line[0].rootNdecayRate > 0? patch[0].rootzone.DOMdecayRate : patch[0].soil_defaults[0][0].DOM_decay_rate),
				(command_line[0].rootNdecayRate > 0? patch[0].soil_defaults[0][0].soil_depth : (command_line[0].root2active > 0.0? patch[0].rootzone.depth * command_line[0].root2active : patch[0].soil_defaults[0][0].active_zone_z)),
				patch[0].soil_defaults[0][0].DOC_adsorption_rate,
				26,patch);
			patch[0].surface_DOC += Nout;
			patch[0].soil_cs.DOC_Qout += Nout;
            if(Nout<0 || Nout!=Nout) printf("update_drainage_land[%d,%e]: return DOC (%e) flux (%e)\n", patch[0].ID,return_flow,patch[0].soil_cs.DOC - (DOC_leached_to_patch/patch[0].area), Nout);
		}//
    
    patch[0].overland_flow += return_flow; //max(0.0, patch[0].detention_store - patch[0].soil_defaults[0][0].detention_store_size);
    //<<--- reset by compute_subsurface_routing.c
    
	/*--------------------------------------------------------------*/
	/*	route water and nitrogen lossed due to infiltration excess */
	/*--------------------------------------------------------------*/
	if ( (patch[0].detention_store > patch[0].soil_defaults[0][0].detention_store_size) && (patch[0].detention_store > ZERO) ){

		Qout = (patch[0].detention_store - patch[0].soil_defaults[0][0].detention_store_size);
		if (command_line[0].grow_flag > 0) {
                Nout = (min(1.0, (Qout/ patch[0].detention_store))) * patch[0].surface_DOC;
                DOC_leached_to_surface = Nout * patch[0].area;
                patch[0].surface_DOC -= Nout;
            
                Nout = (min(1.0, (Qout/ patch[0].detention_store))) * patch[0].surface_DON;
                DON_leached_to_surface = Nout * patch[0].area;
                patch[0].surface_DON -= Nout;
            
                Nout = (min(1.0, (Qout/ patch[0].detention_store))) * patch[0].surface_NO3;
                NO3_leached_to_surface = Nout * patch[0].area;
                patch[0].surface_NO3 -= Nout;
            
                Nout = (min(1.0, (Qout/ patch[0].detention_store))) * patch[0].surface_NH4;
                NH4_leached_to_surface = Nout * patch[0].area;
                patch[0].surface_NH4 -= Nout;
        }// end of grow_flag if statement
        
        route_to_surface = (Qout *  patch[0].area);
        patch[0].detention_store -= Qout;
        patch[0].surface_Qout += Qout;

    }// end of if statement
			

	if (NO3_leached_to_surface < 0.0)
		printf("WARNING %d %lf",patch[0].ID, NO3_leached_to_surface);

    
    
    
    
    
    
    
    
	/*--------------------------------------------------------------*/
	/*	route flow to neighbours				*/
	/*	route n_leaching if grow flag specfied			*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/* regular downslope routing */
	/*--------------------------------------------------------------*/
	if (command_line[0].noredist_flag == 0) {
        d=0;
        for (j = 0; j < patch[0].innundation_list[d].num_neighbours; j++) {
            neigh = patch[0].innundation_list[d].neighbours[j].patch;
            /*--------------------------------------------------------------*/
            /* first transfer subsurface water and nitrogen */  // --------- subsurface
            /*--------------------------------------------------------------*/
            Qin = (patch[0].innundation_list[d].neighbours[j].gamma * route_to_patch) / neigh[0].area;
            if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){neigh[0].fromLAND_Q+=Qin; }//reset @subsurface_routing
            if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){neigh[0].fromRIPARIAN_Q+=Qin; }
            if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){neigh[0].fromSTREAM_Q+=Qin; }
            
            if (Qin < 0) printf("\n warning negative routing from patch %d with gamma %lf", patch[0].ID, total_gamma);
            if (command_line[0].grow_flag > 0) {
                
                //Note that: "xxx_leached_to_patch" is flux, i.e., not areal
                //Note that: "Nin" is an areal adjustment to the neigh due to influx
                Nin = (patch[0].innundation_list[d].neighbours[j].gamma * DON_leached_to_patch) / neigh[0].area;
                neigh[0].soil_ns.DON_Qin += Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){neigh[0].fromLAND_DON+=Nin; }//reset @subsurface_routing
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){neigh[0].fromRIPARIAN_DON+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){neigh[0].fromSTREAM_DON+=Nin; }
                
    
                Nin = (patch[0].innundation_list[d].neighbours[j].gamma * DOC_leached_to_patch) / neigh[0].area;
                neigh[0].soil_cs.DOC_Qin += Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){neigh[0].fromLAND_DOC+=Nin; }//reset @subsurface_routing
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){neigh[0].fromRIPARIAN_DOC+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){neigh[0].fromSTREAM_DOC+=Nin; }
                
                Nin = (patch[0].innundation_list[d].neighbours[j].gamma * NO3_leached_to_patch) / neigh[0].area;
                neigh[0].soil_ns.NO3_Qin += Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){neigh[0].fromLAND_NO3+=Nin; }//reset @subsurface_routing
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){neigh[0].fromRIPARIAN_NO3+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){neigh[0].fromSTREAM_NO3+=Nin; }
                
                
                Nin = (patch[0].innundation_list[d].neighbours[j].gamma * NH4_leached_to_patch) / neigh[0].area;
                neigh[0].soil_ns.NH4_Qin += Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){neigh[0].fromLAND_NH4+=Nin; }//reset @subsurface_routing
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){neigh[0].fromRIPARIAN_NH4+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){neigh[0].fromSTREAM_NH4+=Nin; }
                
                }
            neigh[0].Qin += Qin; // subsurface --> neighbour sat_def
        }// end of for subsurface routing loop

        
        
        /*--------------------------------------------------------------*/
        /* surface downslope routing */
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        /* determine which innundation depth to consider		*/
        /*--------------------------------------------------------------*/
        if (patch[0].num_innundation_depths > 0) {
            innundation_depth = patch[0].detention_store + route_to_surface/patch[0].area;
            d=0;
            while ((innundation_depth > patch[0].innundation_list[d].critical_depth) && (d < patch[0].num_innundation_depths-1)) {
                d++;}// while
            }// if
        else d=0;
        
        for (j = 0; j < patch[0].surface_innundation_list[d].num_neighbours; j++) {

            neigh = patch[0].surface_innundation_list[d].neighbours[j].patch;

            /*--------------------------------------------------------------*/
            /* now transfer surface water and nitrogen */ // -------------- surface (updated Spet 12)
            /*	- first nitrogen					*/
            /*--------------------------------------------------------------*/
            // route_to_surface = patch[0].area * (patch[0].detention_store - patch[0].soil_defaults[0][0].detention_store_size);
            if (command_line[0].grow_flag > 0) {
                if(neigh[0].drainage_type==STREAM){
                    
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * NO3_leached_to_surface) / neigh[0].area;
                    neigh[0].streamflow_NO3 += Nin;
                    if(neigh[0].ID==command_line[0].outletPatchID && patch[0].drainage_type!=STREAM){neigh[0].stormdrained_NO3 += Nin;}
                    
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * NH4_leached_to_surface) / neigh[0].area;
                    neigh[0].streamflow_NH4 += Nin;
                    if(neigh[0].ID==command_line[0].outletPatchID && patch[0].drainage_type!=STREAM){neigh[0].stormdrained_NH4 += Nin;}
                    
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * DON_leached_to_surface) / neigh[0].area;
                    neigh[0].streamflow_DON += Nin;
                    if(neigh[0].ID==command_line[0].outletPatchID && patch[0].drainage_type!=STREAM){neigh[0].stormdrained_DON += Nin;}
                    
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * DOC_leached_to_surface) / neigh[0].area;
                    neigh[0].streamflow_DOC += Nin;
                    if(neigh[0].ID==command_line[0].outletPatchID && patch[0].drainage_type!=STREAM){neigh[0].stormdrained_DOC += Nin;}
                }else{
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * NO3_leached_to_surface) / neigh[0].area;
                    neigh[0].surface_NO3 += Nin;
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * NH4_leached_to_surface) / neigh[0].area;
                    neigh[0].surface_NH4 += Nin;
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * DON_leached_to_surface) / neigh[0].area;
                    neigh[0].surface_DON += Nin;
                    Nin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * DOC_leached_to_surface) / neigh[0].area;
                    neigh[0].surface_DOC += Nin;
                }
            }//if growth
            
            /*--------------------------------------------------------------*/
            /*	- now surface water 					*/
            /*	surface stores should be updated to facilitate transfer */
            /* added net surface water transfer to detention store		*/
            /*--------------------------------------------------------------*/

            // re-work here!!
            Qin = (patch[0].surface_innundation_list[d].neighbours[j].gamma * route_to_surface) / neigh[0].area;
            if( neigh[0].drainage_type==STREAM && neigh[0].ID==command_line[0].outletPatchID && patch[0].drainage_type>0 && patch[0].drainage_type % actionSTORMDRAIN==0 ){
                
                //patch[0].drainage_type!=STREAM
                patch[0].stormdrainYield += Qin; //Spet 17 tracking how much is storm drain yielded at "local patch"
                neigh[0].stormdrained += Qin; //Spet 15 tracking how much storm drain to the "outlet" (basin)
                
                neigh[0].surface_Qin += Qin;
                neigh[0].streamflow += Qin; // short-cut by passing the detention and infilration at stream patch; but trying to seperate return and stormdrain
                // perviously was if neigh[0].drainage_type==STREAM, then all surface water becomes streamflow
                
            }else{
                neigh[0].detention_store += Qin;// need fix this ****
                neigh[0].surface_Qin += Qin;
            }
            
            /*--------------------------------------------------------------*/
            /* try to infiltrate this water					*/
            /* use time_int as duration */
            /*--------------------------------------------------------------*/
            if (neigh[0].detention_store > ZERO) {
                if (neigh[0].rootzone.depth > ZERO) {
                infiltration = compute_infiltration(
                    verbose_flag,
                    neigh[0].sat_deficit_z,
                    neigh[0].rootzone.S,
                    neigh[0].Ksat_vertical,
                    neigh[0].soil_defaults[0][0].Ksat_0_v,
                    neigh[0].soil_defaults[0][0].mz_v,
                    neigh[0].soil_defaults[0][0].porosity_0,
                    neigh[0].soil_defaults[0][0].porosity_decay,
                    (neigh[0].detention_store),
                    time_int,
                    neigh[0].soil_defaults[0][0].psi_air_entry);
                }
                else {
                infiltration = compute_infiltration(
                    verbose_flag,
                    neigh[0].sat_deficit_z,
                    neigh[0].S,
                    neigh[0].Ksat_vertical,
                    neigh[0].soil_defaults[0][0].Ksat_0_v,
                    neigh[0].soil_defaults[0][0].mz_v,
                    neigh[0].soil_defaults[0][0].porosity_0,
                    neigh[0].soil_defaults[0][0].porosity_decay,
                    (neigh[0].detention_store),
                    time_int,
                    neigh[0].soil_defaults[0][0].psi_air_entry);
                }
            }
            else infiltration = 0.0;
            /*--------------------------------------------------------------*/
            /* added an surface N flux to surface N pool	and		*/
            /* allow infiltration of surface N				*/
            /*--------------------------------------------------------------*/
            if ((command_line[0].grow_flag > 0 ) && (infiltration > ZERO)) {
                Nin = ((infiltration / neigh[0].detention_store) * neigh[0].surface_DOC);
                neigh[0].soil_cs.DOC_Qin += Nin;
                neigh[0].surface_DOC -= Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){patch[0].fromLAND_surfsubDOC+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){patch[0].fromRIPARIAN_surfsubDOC+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){patch[0].fromSTREAM_surfsubDOC+=Nin; }
                
                Nin = ((infiltration / neigh[0].detention_store) * neigh[0].surface_DON);
                neigh[0].soil_ns.DON_Qin += Nin;
                neigh[0].surface_DON -= Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){patch[0].fromLAND_surfsubDON+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){patch[0].fromRIPARIAN_surfsubDON+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){patch[0].fromSTREAM_surfsubDON+=Nin; }
                
                Nin = ((infiltration / neigh[0].detention_store) * neigh[0].surface_NO3);
                neigh[0].soil_ns.NO3_Qin += Nin;
                neigh[0].surface_NO3 -= Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){patch[0].fromLAND_surfsubNO3+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){patch[0].fromRIPARIAN_surfsubNO3+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){patch[0].fromSTREAM_surfsubNO3+=Nin; }
                
                Nin = ((infiltration / neigh[0].detention_store) * neigh[0].surface_NH4);
                neigh[0].soil_ns.NH4_Qin += Nin;
                neigh[0].surface_NH4 -= Nin;
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 11 ==0){patch[0].fromLAND_surfsubNH4+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 7 ==0){patch[0].fromRIPARIAN_surfsubNH4+=Nin; }
                if(patch[0].aggregate_ID != neigh[0].aggregate_ID && patch[0].aggregate_ID>0 && patch[0].aggregate_ID % 5 ==0){patch[0].fromSTREAM_surfsubNH4+=Nin; }
            }

            if (infiltration > neigh[0].sat_deficit - neigh[0].unsat_storage - neigh[0].rz_storage) {
                neigh[0].sat_deficit -= (infiltration + neigh[0].unsat_storage + neigh[0].rz_storage);
                neigh[0].unsat_storage = 0.0;
                neigh[0].rz_storage = 0.0;
                neigh[0].field_capacity = 0.0;
                neigh[0].rootzone.field_capacity = 0.0;
            }

            else if ((neigh[0].sat_deficit > neigh[0].rootzone.potential_sat) &&
                (infiltration > neigh[0].rootzone.potential_sat - neigh[0].rz_storage)) {
            /*------------------------------------------------------------------------------*/
            /*		Just add the infiltration to the rz_storage and unsat_storage	*/
            /*------------------------------------------------------------------------------*/
                neigh[0].unsat_storage += infiltration - (neigh[0].rootzone.potential_sat - neigh[0].rz_storage);
                neigh[0].rz_storage = neigh[0].rootzone.potential_sat;
            }
            /* Only rootzone layer saturated - perched water table case */
            else if ((neigh[0].sat_deficit > neigh[0].rootzone.potential_sat) &&
                (infiltration <= neigh[0].rootzone.potential_sat - neigh[0].rz_storage)) {
                /*--------------------------------------------------------------*/
                /*		Just add the infiltration to the rz_storage	*/
                /*--------------------------------------------------------------*/
                neigh[0].rz_storage += infiltration;
            }
            else if ((neigh[0].sat_deficit <= neigh[0].rootzone.potential_sat) &&
                (infiltration <= neigh[0].sat_deficit - neigh[0].rz_storage - neigh[0].unsat_storage)) {
                neigh[0].rz_storage += neigh[0].unsat_storage;
                /* transfer left water in unsat storage to rootzone layer */
                neigh[0].unsat_storage = 0;
                neigh[0].rz_storage += infiltration;
                neigh[0].field_capacity = 0;
            }

            neigh[0].detention_store -= infiltration;

        }// surface routing loop
        
        

        
        
        

	} /* end if noredistribution flag */

	return;

} /*end update_drainage_land.c*/

