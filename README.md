# RHESSysEastCoast

This github repository hold a developing branch of RHESSys (https://github.com/RHESSys/RHESSys).
This is branched out from RHESSys 5.20.

Reasons for this branch:
- RHESSys main branch (https://github.com/RHESSys/RHESSys) serves as a general purpose model, while this branch of RHESSys is heavily and rotatively tested in several catchments (forested and urban) on the U.S. east coast in terms of hydrology, soil moisture pattern, forest ecosystem, and biochemistry cycle & transport.
- model tunings and new features are specially designed for the east coast catchments

new features / modifications for hydrology and soil moisture
- fully real-time dynamic subsurface routing, i.e., daily water table topology/elevation derives the subsurface water movements in terms of direction and magnitude
- house basement interreption on subsurface water movement
- break down "Roads" as paved road (high impervious and parallel to the slope surface) and contour road on mountain/steep slope (contour feature can interrept subsurface water on slopes)
- patch scale controls on surface-subsurface routing and interaction (use flowtable actionCodes to replace drinagetype):
  - 1 = stream drinage
  - 2 = contour road drainage
  - actionSTORMDRAIN (3) = surface storm drinage
  - gw_flag & -gwtoriparian_flag & actionRIPARIAN (7) = deep groundwater seeping out to subsurface
  - gw_flag & actionGWDRAIN (5) = drainage to deep groundwater (constrained by surface imperiousness)
  - grassIrrigation_flag & actionIRRIGRATION (13) & landuse.def [max amount mm/day] (during growth season only; water stress adjustment) 
  - fertilizer_flag & actionFERTILIZE & landuse.def [amount (kgN/m2/mo) and freq (30 days)] (during growth season only)
  - sewer_flag & actionSEWER (11) = subsurface sewer drinage
  - actionPIPEDRAIN (17) = subsurface pipe (non-sewer) drinage
- storm drainage point along the paved road network in urban catchment, as well as quick drainage area (e.g., yards, parks, ..etc)
- streamflow = returnflow + baseflow + stormdrainage + (gw.out) if applied


new features / modifications for forest ecosystem
- sun angle and aspect angle fixed 
- corrected cpool / npool dynamics in vegetation
- excessive cpool/npool goes to DOC/DON
- new metric to evaluate annual vegetation respout: GPP + cpool - Re 


new features / modifications for biochemistry cycle & transport
- new flag "-readinWFdoc_flag" when reading in worldfiles that has hillslope.gwDOC, hillslope.gwDON, patch.soil.DOC, and patch.soil.DON; In the latest model, these DOC related state variables are written out to the worldfile when perform "output_current_state" in the tecfile. Please insect this flag (-readinWFdoc_flag) in the command line when reading in the worldfiles that is created by "output_current_state".
- first version of SAT_solutes code implemention













