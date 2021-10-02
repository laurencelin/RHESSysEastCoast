# RHESSysEastCoast

This github repository hold a developing branch of RHESSys (https://github.com/RHESSys/RHESSys).
This is branched out from RHESSys 5.20.

Reasons for this branch:
- RHESSys main branch (https://github.com/RHESSys/RHESSys) serves as a general purpose model, while this branch of RHESSys is repeatedly and heavily tested in several catchments (forested and urban) on the U.S. east coast in terms of hydrology, soil moisture pattern, forest ecosystem, and biochemistry cycle & transport.
- model tunings and new features are specially designed for the east coast catchment studies.
- specific RHESSysEC modeling tools are develoepd and available at https://github.com/laurencelin/GIS2RHESSys (GIS process and model setup) and https://github.com/laurencelin/R-coded-scripts-for-RHESSys-calibration (model calibration on Rivanna/SLURM based computational clusters)



new features / modifications for hydrology and soil moisture
- fully real-time dynamic subsurface routing, i.e., daily water table topology/elevation derives the subsurface water movements in terms of direction and magnitude; required RHESSysEC specific flowtable inputs (https://github.com/laurencelin/GIS2RHESSys)
- house basement interreption on subsurface water movement in some urban catchments; required RHESSysEC specific flowtable inputs (disabled in the current version)
- break down "Roads" as paved road (high impervious and parallel to the slope surface) and contour road on mountain/steep slope (contour feature can interrept subsurface water on slopes); required RHESSysEC specific flowtable inputs
- storm drainage point along the paved road network in urban catchment, as well as quick drainage area (e.g., yards, parks, ..etc)
- streamflow = returnflow + baseflow + stormdrainage + (gw.out) if applied
- subpatch LULC modeling, i.e., multiple LULC in one patch. In this case, the "old" LULC def becomes useless because a patch cannot be defined by a single LULC def. So in this new version of model, LULC def becomes a parameter file holding fertilizer, septic, sewer, and pond depth information (e.g., https://github.com/laurencelin/GIS2RHESSys/blob/master/lulcCollectionEC.csv). The id "11" is just an id, unrelated to the action code (defined below) in the flow table.
- patch scale controls on surface-subsurface routing and interaction (use flowtable actionCodes and model commandline flags to replace drinagetype); actionCodes are prime numbers; multiple actions can be implemented on the same patch simply by multiplying the actionCodes; for example, three actions (irrigation, fertilizer, and exess water drainage) could occur in a lawn patch and its actionCode is 741 = 13 X 19 X 3; actionCodes are listed below:
  - 1 = stream drinage
  - 2 = contour road drainage
  - actionSTORMDRAIN (3) = surface storm drinage
  - -gw_flag & -gwtoriparian_flag & actionRIPARIAN (7) = deep groundwater seeping out to subsurface
  - -gw_flag & actionGWDRAIN (5) = drainage to deep groundwater (constrained by surface imperiousness)
  - -grassIrrigation_flag & actionIRRIGRATION (13) & landuse.def [max amount mm/day] (during growth season only; water stress adjustment) 
  - -fertilizer_flag & actionFERTILIZE (13) & landuse.def [amount (kgN/m2/mo) and freq (30 days)] (during growth season only)
  - -sewer_flag & actionSEWER (11) = subsurface sewer drinage
  - actionPIPEDRAIN (17) = subsurface pipe (non-sewer) drinage
  - septic source (19)
When GIS2RHESSys sees the LULC cover fractions (https://github.com/laurencelin/GIS2RHESSys), it generates the combination of action codes for the patch. https://github.com/laurencelin/GIS2RHESSys/blob/master/libraries/g2w_cf_RHESSysEC_soil_fullextraction.R#L1627 

new features / modifications for forest ecosystem
- sun angle and aspect angle fixed 
- corrected cpool / npool dynamics in vegetation
- excessive cpool/npool goes to DOC/DON
- calculations of litter and soil organic matter decomposition are based on BIOME-BGC ecological stoichiometry
- new metrics are developed to measure vegetation health and cause of death in the model
- new rules to evaluate the need of vegetation respout in the model; use of annual carbon budget GPP + cpool - Resp. If annual carbon balance is negative in multiple years, respout process will be triggered.
- resolved PAR competition within patch-scale canopy


new features / modifications for biochemistry cycle & transport
- new flag "-readinWFdoc_flag" when reading in worldfiles that has hillslope.gwDOC, hillslope.gwDON, patch.soil.DOC, and patch.soil.DON; In the latest model, these DOC related state variables are written out to the worldfile when perform "output_current_state" in the tecfile. Please insect this flag (-readinWFdoc_flag) in the command line when reading in the worldfiles that is created by "output_current_state".
- first version of SAT_solutes code implemention; vertically seperating subsurface solute transport and surface solute input/generation
- updated on nitrification and denitrification; make these calculations based on vertical profiles of solute and hydrologic condition, as well as patch scale soil moisture distribution



new flags in commandline to control certain model behaviors
- -dynRtZoff
https://github.com/laurencelin/RHESSysEastCoast/blob/0a79e1efa1b35b4f20695253c8ccaff4d713ea96/init/construct_command_line.c#L622
https://github.com/laurencelin/RHESSysEastCoast/search?q=dynRtZoff_flag

It’s a flag to turn off the runtime calculation of root depth. The root biomass or carbon storage grows unbounded over time in the old model. The unbounded root mass is greatly over the empirical relationship, resulting in unreal root depth (over 80 m deep in Coweeta cases in some old simulation). No one noticed because not many people heavily use the growth mode.  

This flag may not need for 7.2+ EC model. I have implemented “improved” carbon allocation in the model under Dickenson and Waring. Usr can set max height and depth for different plant species in the veg.def.  These max heights and depths redistribute the carbon allocation. Clare has tested this implementation with me and gave some great inputs. You may contact her for additional details. This implementation is not just fixing the numerical calculation of carbon allocation, but also imply the age effects on plant carbon sequestration and N uptake. As mature plant stand reaches to max depth or height, carbon sequestration starts to slow down and become more net zero to the respiration. Note this is different from the tree stand density. 

https://github.com/laurencelin/RHESSysEastCoast/blob/master/cn/compute_potential_N_uptake_Dickenson.c#L130

After all, you can always turn this off by not using the flag and set max height and depth to some large numbers.

- -fracDirectNdep #
https://github.com/laurencelin/RHESSysEastCoast/search?q=fracDirectNdep

In 100% vegetated patch, all N deposition goes to stratum surface, which is unavailable plant uptake or decomposition. Rain water moves the stratum surface N to ground surface pool.  Option of this flag with a fraction [0-1] allows some fixability of this strong assumption, specially in urbanized area. 


- -soilCNadaptation_flag
https://github.com/laurencelin/RHESSysEastCoast/search?q=soilCNadaptation_flag

This is my attempt to improve the old decomposition. The old decomposition assumed fixed C:N ratios on all four soil pools, which is not true. The last soil pool C:N is seemingly fixed around 17.8 from many long term forest decomposition studies, but the first 1-3 soil pools are not fixed at all, which should be related to above-ground vegetation. 

This flag is to allow dynamic C:N ratios in 1-3 soil pools. It’s just my attempt, not official. Maybe someone can work on this part further. This feature has great implication for landuse change simulation, as new vegetation will slowly modify the below-ground soil C:N. 










