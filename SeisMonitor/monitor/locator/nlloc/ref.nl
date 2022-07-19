# =============================================================================
#  NonLinLoc programs control file
#
#  Global location example
#
#  NonLinLoc Version 3.0 - APR2004
#
#  See "Control File" and "Running the Sample Location" pages
#     in the NonLinLoc on-line documentation: http://www.alomax.net/nlloc
# =============================================================================


# = comment

# non-nested include files allowed, use:
# INCLUDE <include_file_name>


# =============================================================================
# =============================================================================
# Generic control file statements
# =============================================================================
#
#

# control (CONTROL message_flag (0:silent,1:few messages,2:verbose,...),
#               RandomNumSeed)
CONTROL 1 54321

# -----------------------------------------------------------------------------
# lat/long to rect grid transformation
# -----------------------------------------------------------------------------

# map projection / transformation
# (TRANS type <params>)
#    (char[])   search_type (SIMPLE, LAMBERT)
#    <params>:
#       GLOBAL
#       SIMPLE LatOrig  LongOrig  RotCW
#       LAMBERT  RefEllipsoid LatOrig  LongOrig
#                   FirstStdParal  SecondStdParal   RotCW
#
#    RefEllipsoid choices:
#                   WGS-84, GRS-80, WGS-72, Australian, Krasovsky,
#                   International, Hayford-1909, Clarke-1880, Clarke-1866,
#                   Airy, Bessel, Hayford-1830, Sphere
#
#
# whole Earth, spherical coordinates
TRANS  GLOBAL

# maplines (MAPLINE id_num, name, red, green, blue,
#       linestyle (SOLID, DASHED, DOTTED, DASHDOT))
#MAPLINE  GMT_LONLAT ./data_geog/map.prov.line  0.0 0.0 0.0  SOLID
#MAPLINE  XY_LONLAT /maps/brasil/PG.line  0.0 0.0 0.0  SOLID

#
#global/iasp91
# =============================================================================
# END of Generic control file statements
# =============================================================================
# =============================================================================


INCLUDE /opt/seiscomp3/etc/nll/sta_list_global.in

#/opt/seiscomp3/etc/nll
# =============================================================================
# =============================================================================
# NLLoc control file statements
# =============================================================================
#
#


# signature
# (LOCSIG signature)
LOCSIG Anthony Lomax (www.alomax.net)


# comment
# (LOCCOM comment)
LOCCOM NEIC events


# input  grid filenames root, output filename
# (LOCFILES <obs file> obs_type  <travel-time grid files path/root> <output file path/root> [iSwapBytes])
#    (char[])  obs_type : (NLLOC_OBS, HYPO71, HYPOELLIPSE, RENASS_DEP, SEISAN)
#
#
#LOCFILES ./obs/20031203.0737_FIJI.neic NEIC  ./taup/ak135/ak135  ./loc/global 1
LOCFILES ./obs/*.neic NEIC  ./taup/ak135/ak135  ./loc/global 1
#


# output hypocenter file types
# (LOCHYPOUT type1, type2, ...)
#    (char[])   typeN (SAVE_NLLOC_ALL, SAVE_NLLOC_SUM,
                SAVE_HYPO71_ALL, SAVE_HYPO71_SUM, SAVE_HYPOELL_ALL, SAVE_HYPOELL_SUM,
                SAVE_HYPOINV_SUM)
LOCHYPOUT SAVE_NLLOC_ALL


# search type
# (LOCSEARCH search_type <params>)
#    (char[])   search_type (GRID, MET (Metropolis), SA (Simulated Annealing))
#    <params>:
#       GRID NumScatterSamples
#       MET  NumSamples NumLearn NumEquil BeginSave NumSkip
#               StepInit StepMin StepFact ProbMin
#       OCT init_num_cells_x, init_num_cells_y, init_num_cells_z,
#               min_node_size, max_num_nodes, num_scatter,
#               use_stations_density (>=1 = Weights oct-tree cell prob values used for subdivide decision
#                               in proportion to distance to nearest station to oct-tree cell.
#                               Gives higher search priority to cells near and containing stations,
#                               stablises convergence to local events when global search used
#                               with dense cluster of local stations.),
#               stop_on_min_node_size (1 = stop search when first min_node_size reached,
#                               0 = stop subdividing a given cell when min_node_size reached.)
#global
LOCSEARCH  OCT 96 48 6 0.05 50000 10000 4 0
#Euro-Med
#LOCSEARCH  OCT 21 21 21  0.05 50000 5000 4 0


# location grids description
# (LOCGRID  num_grid_x  num_grid_y  num_grid_z
#       orig_grid_x  orig_grid_y  orig_grid_z
#       d_grid_x d_grid_y d_grid_z
#       type save_flag
#    (float) num_grid_x/y/z : number of nodes along x/y/z axis
#    (float)    orig_grid_x : x location of grid origin (0,0,0) in km pos east
#    (float)    orig_grid_y : y location of grid origin (0,0,0) in km pos north
#    (float)    orig_grid_z : z location of grid origin (0,0,0) in km pos down
#    (float)   d_grid_x/y/x : grid spacing along  x/y/z axis
#    (char[])  type : (PROB_DENSITY, MISFIT)
#    (char[])  save_flag : (SAVE, NO_SAVE)
# For Grid search, first grid is used for initial search.  Subsequent grids are
# shifted in x/y/z so that they are centered on the minimum misfit hypocenter
# x/y/z of the previous grid if x/y/z < -1.0e20.
#
#global
LOCGRID  361 181 601  -180.0 -90.0 0.0  1.0 1.0 1.0   PROB_DENSITY  SAVE
#Euro-Med
#LOCGRID  61 41 601  -20.0 20.0 0.0  1.0 1.0 1.0   PROB_DENSITY  SAVE
#Euro-Med (Crust only)
#LOCGRID  61 41 31  -20.0 20.0 0.0  1.0 1.0 1.0   PROB_DENSITY  SAVE


# method
# (LOCMETH method)
#    (char[])   method (GAU_ANALYTIC, EDT, EDT_OT_WT)
#          GAU_ANALYTIC - L2 norm following Tarantola and Valette (1982)
#          EDT - Equal Differential Time (see )
#          EDT_OT_WT - Weights EDT sum prob by variance of OT estimated over all pairs of readings.
#                              Downweights locations with inconsistent OT estimates.
#    (float)   maximum_dist_sta_to_grid (use very large value for no max)
#    (int)   minimum_number_phases for location
#    (int)   maximum_number_phases for location (-1 for no max)
#    (int)   minimum_number_S_phases for location (-1 for no min)
#    (float)   Vp/Vs ratio (< 0.0 to use S travel time grids)
#    (int)   maximum_number_3D_grids to attempt to read into memory (-1 for no max)
#    (float)   minimum_dist_sta_to_grid (-1 for no min)
#    (int)   flag indicating if duplicate arrivals used for location (1=reject, 0=use if time diff < sigma / 2)
#            duplicate arrivals have same station label and phase name
#LOCMETH GAU_ANALYTIC 1.0e6 4 100 -1 -1.7  6
LOCMETH EDT_OT_WT 1.0e6 4 50 -1 -1.7 6 -1.0 1


# fixed origin time
# (LOCFIXOTIME year month day hour min sec)
#    (int)   year month day hour min
#    (float)   sec
#LOCFIXOTIME 1995 04 21 08 02 57.09


# gaussian model error parameters
# (LOCGAU Sigma_T (s), CorrLen (km))
#LOCGAU 2.0 0.0
LOCGAU 0.5 0.0


# travel-time dependent gaussian model error parameters
# (LOCGAU2 SigmaTfraction,  SigmaTmin (s),  SigmaTmax (s))
# travel time error is travel_time*SigmaTfraction, with max/min value = SigmaTmin/SigmaTmax
LOCGAU2 0.01 0.05 2.0


# phase identifier mapping
# (LOCPHASEID phase  phase_id0 phase_id1 ...)
#
# examples for P and S
#LOCPHASEID  P   P p
#LOCPHASEID  S   S s
#
# examples for global, first arriving P and S
#LOCPHASEID  P   P p Pn Pdiff PKP PKiKP PKIKP
#LOCPHASEID  S   S s Sn Sdiff SKS SKiKS SKIKS
#ToIgnoreS#LOCPHASEID  S   $


# quality to error mapping (for HYPO71, etc)
# (LOCQUAL2ERR Err0 Err1 Err2 ... )
#
LOCQUAL2ERR 0.2 0.5 1.0 2.0 99999.9


# phase statistics parameters
# (LOCPHSTAT RMS_Max, NRdgs_Min, Gap_Max, P_ResMax, S_ResMax)
#    (float)   RMS_Max : max hypocenter RMS to include in ave res
#    (float)   NRdgs_Min : min hypocenter num readings to include in ave res
#    (float)   Gap_Max : max hypocenter gap (deg) to include in ave res
#    (float)   P_ResMax : max abs(P res) to include in ave res
#    (float)   S_ResMax : max abs(S res) to include in ave res
LOCPHSTAT 9999.0 -1 9999.0 1.0 1.0


# take-off angles mode & minimum quality
# (LOCANGLES angles_mode, min_quality)
#    (char[])   angles_mode (ANGLES_YES, ANGLES_NO)
#    (integer)   min_quality : minimum quality to use take-off angles
LOCANGLES ANGLES_NO 5


# magnitude calculation method
# (LOCMAG magnitude_type <params>)
#    (char[])   magnitude_type (ML_HB (ML, Hutton Boore))
#    <params>:
#       ML_HB  amp_fact n K
#LOCMAG ML_HB 1.0 1.110 0.00189


# station/inst/comp parameters (for specifying component specific parameters, i.e. constants for magnitude calculation)
# (LOCCMP name inst comp amp_fact sta_corr)
#    (char[])   name  : station identifier (after alias evaluation, without trailing underscore "_")
#    (char[])   name  : inst identifier (use '?' for don't care)
#    (char[])   name  : comp identifier (use '?' for don't care)
#    (float)    amp_fact: amplitude factor, will be multiplied by amplitude
#    (float)    sta_corr: mganitude correction
#
# example:
#
#LOCCMP CDR ? ? 1.0 0.0
#


# station name alias (for aliasing sta names, for date validation and for
#    phase time delays)
# (LOCALIAS name alias year mo day year mo day)
#    (char[])   name  : station identifier on input
#    (char[])   alias : station identifier for travel time grid on output
#                    NOTE: a trailing underscore "_" in aliases will only be
#                          used for time grid identification, not for output
#    (ints)    year mo day : start date of validity (0 0 0 = no start date)
#    (ints)    year mo day : end date of validity  (9999 99 99 = no end date)
#
#   Note:
#       Alias evaluation is applied recursively, beware of infinite recursion!
#       P and S delays from last alias only are used!
#
# example:
#
#LOCALIAS ART ART_      1996 05 29      1996 09 18   0.03  0.08
#
