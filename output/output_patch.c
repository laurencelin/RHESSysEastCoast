/*--------------------------------------------------------------*/
/* 																*/
/*					output_patch						*/
/*																*/
/*	output_patch - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_patch - outputs current contents of a patch.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_patch(										*/
/*					struct	patch_object	*patch,				*/
/*					struct	date	date,  						*/
/*					FILE 	*outfile)							*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	outputs spatial structure according to commandline			*/
/*	specifications to specific files							*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	We only permit one fileset per spatial modelling level.     */
/*	Each fileset has one file for each timestep.  				*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

void	output_patch(
					 int basinID, int hillID, int zoneID,
					 struct	patch_object	*patch,
					 struct	zone_object	*zone,
					 struct	date	current_date,
					 FILE *outfile)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	int check, c, layer;
	double alai, asub, apsn, litterS, aheight, agsi;
    double mean_gl, gl_scalar;
    double coverf, totalcoverf;
    struct mult_conduct_struct *mc;
    double mc_scalar;
    
	if (patch[0].litter.rain_capacity > ZERO)
		litterS = patch[0].litter.rain_stored / patch[0].litter.rain_capacity;
	else
		litterS = 1.0;

	apsn = 0.0;
	asub = 0.0;
	alai = 0.0;
	aheight = 0.0;
    agsi = 0.0;
    mean_gl = 0.0;
    gl_scalar= 0.0;
    totalcoverf = 0.0;
    
	for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
		for ( c=0 ; c<patch[0].layers[layer].count; c++ ){
            
            coverf = patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cover_fraction;
            mc = &patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].mult_conductance;
            mc_scalar = mc->APAR * mc->tavg * mc->LWP * mc->CO2 * mc->tmin * mc->vpd;
            
			apsn += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cs.net_psn ;
			asub += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].sublimation;
			alai += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].epv.proj_lai;
			aheight += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].epv.height;
            
            agsi += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].phen.gsi;

            mean_gl += coverf * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].defaults[0][0].epc.gl_smax * mc_scalar;
            gl_scalar += coverf * mc_scalar;
            totalcoverf += coverf;
            
		}
	}
    if(totalcoverf>1.0){ mean_gl/=totalcoverf; gl_scalar/=totalcoverf;  }
	check = fprintf(outfile,"%d-%d-%d %d %lf %lf %lf %lf %lf %lf \
                             %lf %lf %lf %lf %lf\n", // \
                             //%lf %lf %lf %lf %lf \
                             //%lf %lf %lf %lf %lf \
                             //%lf %lf %lf\n",
					current_date.year, current_date.month, current_date.day, //1,2,3,
					patch[0].ID, //4
                    (patch[0].Qout_total - patch[0].Qin_total) * 1000.0, //5
                    (patch[0].surface_Qout_total - patch[0].surface_Qin_total) * 1000.0, //6
                    patch[0].detention_store*1000.0, //7
                    patch[0].stormdrainYield*1000.0, //8
                    patch[0].overland_flow * 1000.0, //9
                    patch[0].rain_throughfall*1000.0, //10
					(patch[0].rain_throughfall - patch[0].recharge)*1000.0,//11
                    (patch[0].cap_rise - patch[0].unsat_drainage)*1000.0,//12
                    patch[0].sat_deficit_z*1000.0,//13
                    (patch[0].sat_deficit>0)? (patch[0].sat_deficit>patch[0].rootzone.potential_sat? (patch[0].rz_storage+patch[0].unsat_storage)/patch[0].sat_deficit : patch[0].rz_storage/patch[0].sat_deficit) : -1, //14
                    (patch[0].sat_deficit>0)? (patch[0].sat_deficit>patch[0].rootzone.potential_sat? patch[0].rz_storage/patch[0].rootzone.potential_sat : patch[0].rz_storage/patch[0].sat_deficit) : -1); //15
                    //patch[0].sat_deficit*1000.0, //patch[0].sat_DOC*1000.0,
                    //patch[0].rootzone.potential_sat*1000.0, //patch[0].soil_ns.nitrate*1000.0,
                    
                    //---extra
//                    patch[0].sat_NO3*1000.0,
//                    patch[0].sat_NH4*1000.0,
//                    patch[0].soil_ns.sminn*1000.0,
//                    patch[0].soil_cs.DOC,
//                    patch[0].z,
//                    patch[0].ndf.denitrif*1000.0,
//                    patch[0].ndf.sminn_to_nitrate*1000.0,
//                    patch[0].ndf.plant_avail_uptake*1000.0, // plant uptake?
//                    patch[0].ndf.net_mineralized*1000.0, // net from decay
//                    (patch[0].ndf.net_mineralized - patch[0].ndf.mineralized)*1000.0, // decay uptake?
//                    patch[0].soil_ns.NO3_Qout_total*1000.0,
//                    patch[0].soil_ns.NO3_Qin_total*1000.0,
//                    patch[0].sat_deficit_z>0? patch[0].rootzone.depth/patch[0].sat_deficit_z : 1000.0
//                    );
	

	if (check <= 0) {
		fprintf(stdout, "\nWARNING: output error has occured in output_patch, file");
	}
	return;
} /*end output_patch*/
