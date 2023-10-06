# -*- mode:Makefile -*-
#
# To generate all surface data sets, run:
# make -f Makefile.data all
#
# To generate a single dataset, run make with the name of the rule you
# want to build. For example, to generate the crop data set for 1x1_numaIA:
#
# make -f Makefile.data crop-numa
#
# NOTE: The default behavior is to parallelize data set creation using
# the batch system by submitting jobs to the batch queue (on cheyenne).
#

# Set up special characters
null  :=

# Set a few things needed for batch handling
PROJECT = $(shell cat $(HOME)/.cesm_proj)
LOGOUT = $@.stdout.txt
PWD = $(shell pwd)

BATCHJOBS_ch = qsub
ifeq ($(PROJECT),$(null))
   $(error Can NOT find PROJECT number from ~/.cesm_proj file create it and try again)
endif
BATCHJOBS = $(BATCHJOBS_ch)

MKSURFDATA = $(PWD)/gen_mksurfdata_jobscript_multi.py --account $(PROJECT)
SUBSETDATA = $(PWD)/../site_and_regional/subset_data
MODIFYSURF = $(PWD)/../modify_input_files/fsurdat_modifier --overwrite

CDATE = $(shell date +%y%m%d)

# subset_data options
#
SUBSETDATA_POINT = $(SUBSETDATA) point --silent --overwrite --uniform-snowpack --cap-saturation --crop --outdir .
SUBSETDATA_POINT_ALLLU = $(SUBSETDATA_POINT) --include-nonveg
SUBSETDATA_POINT_URBAN = $(SUBSETDATA_POINT) --include-nonveg

# Subset data sites...
SUBSETDATA_1X1_BRAZIL := --lat -7 --lon -55 --site 1x1_brazil
SUBSETDATA_1X1_NUMAIA := --lat 40.6878 --lon 267.0228 --site 1x1_numaIA
SUBSETDATA_1X1_SMALL  := --lat 40.6878 --lon 267.0228 --site 1x1_smallvilleIA \
		--dompft 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 \
		--pctpft 6.5 1.5 1.6 1.7 1.8 1.9 1.5 1.6 1.7 1.8 1.9 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5
# NOTE: The 1850 smallvilleIA site is constructed to start with 100% natural vegetation, so we can test transition to crops
SUBSETDATA_1X1_SMALL1850  := --lat 40.6878 --lon 267.0228 --site 1x1_smallvilleIA --dompft 13 --pctpft 100

MEXICOCITY_TMP_FNAME       := surfdata_1x1_mexicocityMEX_hist_2000_78pfts.nc
MEXICOCITY_FNAME           := surfdata_1x1_mexicocityMEX_hist_2000_78pfts_c$(CDATE).nc
VANCOUVER_TMP_FNAME        := surfdata_1x1_vancouverCAN_hist_2000_78pfts.nc
VANCOUVER_FNAME            := surfdata_1x1_vancouverCAN_hist_2000_78pfts_c$(CDATE).nc
URBALPHA_TMP_FNAME         := surfdata_1x1_urbanc_alpha_hist_2000_78pfts.nc
URBALPHA_TMP2_FNAME        := surfdata_1x1_urbanc_alpha_hist_2000_78pfts_tmp.nc
URBALPHA_FNAME             := surfdata_1x1_urbanc_alpha_hist_2000_78pfts_c$(CDATE).nc
SUBSETDATA_1X1_MEXICOCITY  := --lat 19.5 --lon 260.5 --site 1x1_mexicocityMEX --out-surface $(MEXICOCITY_TMP_FNAME)
SUBSETDATA_1X1_VANCOUVER   := --lat 49.5 --lon 236.5 --site 1x1_vancouverCAN --out-surface  $(VANCOUVER_TMP_FNAME)
SUBSETDATA_1X1_URBALPHA    := --lat -37.7308 --lon 0 --site 1x1_urbanc_alpha --out-surface $(URBALPHA_TMP_FNAME)

# ne120np4 and hi-res are for high resolution, ne16np4 is for mid-resolution testing
# low-res is for low resolutions for testing
# nldas is for NWP working with WRF
STANDARD = \
	global-present \
	global-present-low-res \
	global-present-nldas

TROPICS = \
	crop-tropics-present

CROP = \
	crop-global-present \
	crop-global-present-low-res \
	crop-global-present-ne16np4 \
	crop-global-present-ne120np4 \
	crop-global-1850 \
	crop-global-1850-low-res \
	crop-global-1850-ne16np4 \
	crop-global-1850-ne120np4 \
	crop-global-future
all : standard tropics crop urban landuse-timeseries all-subset

all-subset : \
	1x1_brazil-tropics-present \
	crop-tropics-historical \
	crop-tropics-transient \
	crop-numa-present \
	crop-numa-historical \
	crop-smallville \
	crop-smallville-historical \
	urban-present urban-alpha

DEBUG:
	@echo "HOST             := $(HOST)"
	@echo "PROJECT          := $(PROJECT)"
	@echo "MKSURFDATA       := $(MKSURFDATA)"
	@echo "SUBSETDATA_POINT := $(SUBSETDATA_POINT)"
	@echo "MODIFYSURF       := $(MODIFYSURF)"
	@echo "BATCHJOBS        := $(BATCHJOBS)"
	@echo "CDATE            := $(CDATE)"
	@echo "RM               := $(RM)"

#
# standard
#
standard : $(STANDARD)

global-present : FORCE
	$(MKSURFDATA)  --number-of-nodes 24 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

global-present-low-res : FORCE
	$(MKSURFDATA)  --number-of-nodes 10 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

global-present-nldas : FORCE
	$(MKSURFDATA)  --number-of-nodes 24 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

#
# tropics
#
tropics : $(TROPICS)

crop-tropics-present : FORCE
	$(MKSURFDATA)  --number-of-nodes 1 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

1x1_brazil-tropics-present : FORCE
	$(SUBSETDATA_POINT_ALLLU) --create-surface $(SUBSETDATA_1X1_BRAZIL)

crop-tropics-historical : FORCE
	$(SUBSETDATA_POINT_ALLLU) --create-surface $(SUBSETDATA_1X1_BRAZIL) --cfg-file default_data_1850.cfg

crop-tropics-transient : FORCE
	$(SUBSETDATA_POINT_ALLLU) --create-landuse $(SUBSETDATA_1X1_BRAZIL)

#
# crop
#
crop : $(CROP)

crop-global-present : FORCE
	$(MKSURFDATA)  --number-of-nodes 48 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-present-low-res : FORCE
	$(MKSURFDATA)  --number-of-nodes 10 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-present-ne16np4 : FORCE
	$(MKSURFDATA)  --number-of-nodes 24 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-present-ne120np4 : FORCE
	$(MKSURFDATA)  --number-of-nodes 96 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-1850 : FORCE
	$(MKSURFDATA)  --number-of-nodes 48 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-1850-low-res : FORCE
	$(MKSURFDATA)  --number-of-nodes 10 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-1850-ne16np4 : FORCE
	$(MKSURFDATA)  --number-of-nodes 24 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-1850-ne120np4 : FORCE
	$(MKSURFDATA)  --number-of-nodes 96 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-present-0.125 : FORCE
	$(MKSURFDATA)  --number-of-nodes 48 --scenario $@ --jobscript-file $@.sh --walltime 04:00:00
	$(BATCHJOBS) $@.sh

crop-numa-present : FORCE
	$(SUBSETDATA_POINT_ALLLU) --create-surface $(SUBSETDATA_1X1_NUMAIA)

crop-numa-historical : FORCE
	$(SUBSETDATA_POINT_ALLLU) --create-surface $(SUBSETDATA_1X1_NUMAIA) --cfg-file default_data_1850.cfg

crop-smallville : FORCE
	$(SUBSETDATA_POINT) --create-surface $(SUBSETDATA_1X1_SMALL)

# Note that the smallville 1850 dataset is entirely natural vegetation. This
# facilitates testing a transient case that starts with no crop, and then later
# adds crop (to make sure that it works properly to add crop in a grid cell
# where there used to be no crop).
crop-smallville-historical : FORCE
	$(SUBSETDATA_POINT) --create-surface $(SUBSETDATA_1X1_SMALL1850) --cfg-file default_data_1850.cfg

#
# Crop with future scenarios
#
crop-global-future : crop-global-SSP2-4.5 crop-global-SSP2-4.5-low-res crop-global-SSP2-4.5-hi-res

crop-global-SSP2-4.5 : FORCE
	$(MKSURFDATA)  --number-of-nodes 24 --scenario $@ --jobscript-file $@.sh --walltime 04:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-low-res : FORCE
	$(MKSURFDATA)  --number-of-nodes 24 --scenario $@ --jobscript-file $@.sh --walltime 04:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-hi-res : FORCE
	$(MKSURFDATA)  --number-of-nodes 96 --scenario $@ --jobscript-file $@.sh --walltime 04:00:00
	$(BATCHJOBS) $@.sh

#
# urban
#
urban : urban-present urban-alpha

urban-present : mexicocity vancouver

mexicocity : FORCE
	$(SUBSETDATA_POINT_URBAN) --create-surface $(SUBSETDATA_1X1_MEXICOCITY)
	$(MODIFYSURF) modify_1x1_mexicocityMEX.cfg -i $(MEXICOCITY_TMP_FNAME) -o $(MEXICOCITY_FNAME)
	$(RM) $(MEXICOCITY_TMP_FNAME)

vancouver : FORCE
	$(SUBSETDATA_POINT_URBAN) --create-surface $(SUBSETDATA_1X1_VANCOUVER)
	$(MODIFYSURF) modify_1x1_vancouverCAN.cfg -i $(VANCOUVER_TMP_FNAME) -o $(VANCOUVER_FNAME)
	$(RM) $(VANCOUVER_TMP_FNAME)

# NOTE(bja, 2015-01) skip abort on invalid data necessary as of 2015-01. See
# /glade/p/cesm/cseg/inputdata/lnd/clm2/surfdata_map/README_c141219
urban-alpha : FORCE
	$(SUBSETDATA_POINT_URBAN) --create-surface $(SUBSETDATA_1X1_URBALPHA)
	# Set number of nlevurb dimension to 4 for this site
	ncks -O --dmn nlevurb,0,3 -o $(URBALPHA_TMP2_FNAME) $(URBALPHA_TMP_FNAME)
	$(MODIFYSURF) modify_1x1_urbanc_alpha.cfg -i $(URBALPHA_TMP2_FNAME) -o $(URBALPHA_FNAME)
	$(RM) $(URBALPHA_TMP_FNAME) $(URBALPHA_TMP2_FNAME)

#
# landuse timeseries with subset data
#
landuse-timeseries : landuse-timeseries-smallville
# NOTE: TODO: This needs to be chagned to use subset_data when transient configurations are resolved (see Issue #1673
landuse-timeseries-smallville : FORCE
#	$(MKSURFDATA) -no_surfdata -glc_nec 10 -y 1850-1855 -r 1x1_smallvilleIA \
#		-pft_idx 17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78 \
#		-pft_frc 6.5,1.5,1.6,1.7,1.8,1.9,1.5,1.6,1.7,1.8,1.9,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5 \
#		-dynpft single_point_dynpft_files/landuse_timeseries_smallvilleIA_hist_simyr1850-1855.txt \
#		$(BACKGROUND)

#
# clean up the working directory by removing generated files
#
clean : FORCE
	-rm *~

clobber : clean
	-rm surfdata_*.nc surfdata_*.log surfdata_*.namelist

#
# generic rule to force things to happen
#
FORCE :

