# -*- mode:Makefile -*-
#
# Before running "make urban-alpha" or any target that includes it,
# execute "module load nco" first.
#
# To generate all surface data sets, run:
# make all
#
# To generate all surface data sets that use subset_data to be generated:
# make all-subset
#
# To generate a single dataset, run make with the name of the rule you
# want to build. For example, to generate the crop data set for 1x1_numaIA:
#
# make crop-numa
#
# NOTE: The default behavior is to parallelize data set creation using
# the batch system by submitting jobs to the batch queue (on Derecho).
#
# In all cases where "--scenario $@" appears, the code executes the
# recipe for the specific target/scenario that it finds in
# ../../python/ctsm/toolchain/gen_mksurfdata_jobscript_multi.py
# There are a few grids mentioned explicitly here in the Makefile and
# many others that only appear in the gen_mksurfdata_jobscript_multi.py
# Look for resolutions or groups of resolutions in that python code.
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

MKSURFDATA = $(PWD)/gen_mksurfdata_jobscript_multi --account $(PROJECT)
SUBSETDATA = $(PWD)/../site_and_regional/subset_data
MODIFYSURF = $(PWD)/../modify_input_files/fsurdat_modifier --overwrite

CDATE = $(shell date +%y%m%d)

# subset_data options
#
SUBSETDATA_POINT = $(SUBSETDATA) point --silent --overwrite --uniform-snowpack --cap-saturation --crop --outdir .
SUBSETDATA_POINT_ALLLU = $(SUBSETDATA_POINT) --include-nonveg
SUBSETDATA_POINT_URBAN = $(SUBSETDATA_POINT) --include-nonveg

# Subset data sites...
SUBSETDATA_1X1_BRAZIL := --lat -7 --lon 305 --site 1x1_brazil
SUBSETDATA_1X1_NUMAIA := --lat 40.6878 --lon 267.0228 --site 1x1_numaIA
SUBSETDATA_1X1_SMALL_IA  := --lat 40.6878 --lon 267.0228 --site 1x1_smallvilleIA \
		--dompft 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 \
		--pctpft 6.5 1.5 1.6 1.7 1.8 1.9 1.5 1.6 1.7 1.8 1.9 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5
SUBSETDATA_1X1_SMALL_BR  := --lat -12.9952 --lon 305.3233 --site 1x1_cidadinhoBR \
		--dompft 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 \
		--pctpft 6.5 1.5 1.6 1.7 1.8 1.9 1.5 1.6 1.7 1.8 1.9 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5
# NOTE: The 1850 smallvilleIA site is constructed to start with 100% natural vegetation, so we can test transition to crops
SUBSETDATA_1X1_SMALL1850  := --lat 40.6878 --lon 267.0228 --site 1x1_smallvilleIA --dompft 13 --pctpft 100
SUBSETDATA_1X1_SMALLTRANSIENT := --lat 40.6878 --lon 267.0228 --site 1x1_smallvilleIA

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
# STANDARD means no crop, so 16 pfts
# global-hist-1850-f45 is used by FATES and we expect it to be phased out
STANDARD = \
	global-potveg \
	global-present \
	global-present-low-res \
	global-hist-1850-f19 \
	global-hist-1850-f45 \

CROP = \
	crop-global-future \
	crop-global-hist-f09 \
	crop-global-1850-ne120 \
	crop-global-present-ne120 \
	crop-global-present-nldas \
	crop-global-present-ne30 \
	crop-global-present \
	crop-global-present-low-res \
	crop-global-present-ne16 \
	crop-global-present-mpasa480 \
	crop-global-1850 \
	crop-global-1850-low-res \
	crop-global-1850-ne16 \
	crop-global-1850-ne30 \
	crop-global-1850-mpasa480 \

# Start all with all-subset because user is bound to forget to first run
# module load nco
# Usually, include global-present-ultra-hi-res temporarily while
# generating datasets for a new release
# Eventually, include tower, when we get the corresponding scripts to
# use absolute paths for their input files; error message posted here:
# https://github.com/ESCOMP/CTSM/pull/2500 under "PLUMBER2"
all : all-subset global-present-ultra-hi-res crop standard crop-tropics-present

# These are all the surface datasets generated by subset_data
# This runs interactively and does not send jobs to the batch queue
all-subset : \
	1x1_brazil-present \
	1x1_brazil-transient \
	1x1-numa-present \
	1x1-smallville-present \
	1x1-smallville-1850 \
	1x1-smallville-transient \
	1x1-cidadinho-present \
	urban

# These are all the tower site surface datasets
tower : neon plumber2


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

global-potveg : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

global-present : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

global-present-low-res : FORCE
	$(MKSURFDATA)  --number-of-nodes 1 --tasks-per-node 64 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

global-hist-1850-f19 : FORCE
	$(MKSURFDATA)  --number-of-nodes 8 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh
global-hist-1850-f45 : FORCE
	$(MKSURFDATA)  --number-of-nodes 2 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh
#
# Ultra high resolutions (Don't do by default user should select this by hand)
#
global-present-ultra-hi-res : FORCE
	$(MKSURFDATA)  --number-of-nodes 50 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

#
# tropics
#

crop-tropics-present : FORCE
	$(MKSURFDATA)  --number-of-nodes 1 --tasks-per-node 5 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

1x1_brazil-present : FORCE
	$(SUBSETDATA_POINT_ALLLU) --create-surface $(SUBSETDATA_1X1_BRAZIL)

1x1_brazil-transient : FORCE
	$(SUBSETDATA_POINT_ALLLU) --create-surface --create-landuse --surf-year 1850 $(SUBSETDATA_1X1_BRAZIL) --cfg-file ../site_and_regional/default_data_1850.cfg

#
# crop
#
crop : $(CROP)

crop-global-present : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-present-low-res : FORCE
	$(MKSURFDATA)  --number-of-nodes 1 --tasks-per-node 64 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-present-ne16 : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-present-ne30 : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-present-ne120 : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-present-mpasa480 : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-present-nldas : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-1850 : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-1850-low-res : FORCE
	$(MKSURFDATA)  --number-of-nodes 1 --tasks-per-node 64 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-1850-ne16 : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-1850-ne30 : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-1850-ne120 : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-1850-mpasa480 : FORCE
	$(MKSURFDATA)  --number-of-nodes 4 --scenario $@ --jobscript-file $@.sh --walltime 01:00:00
	$(BATCHJOBS) $@.sh

crop-global-hist : FORCE
	$(MKSURFDATA)  --number-of-nodes 72 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-hist-f09 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-hist-low-res : FORCE
	$(MKSURFDATA)  --number-of-nodes 14 --tasks-per-node 32 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-hist-ne16 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-hist-ne30 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

1x1-numa-present : FORCE
	$(SUBSETDATA_POINT_ALLLU) --create-surface $(SUBSETDATA_1X1_NUMAIA)

1x1-smallville-present : FORCE
	$(SUBSETDATA_POINT) --create-surface $(SUBSETDATA_1X1_SMALL_IA)

# Note that the smallville 1850 dataset is entirely natural vegetation. This
# facilitates testing a transient case that starts with no crop, and then later
# adds crop (to make sure that it works properly to add crop in a grid cell
# where there used to be no crop).
1x1-smallville-1850 : FORCE
	$(SUBSETDATA_POINT) --create-surface $(SUBSETDATA_1X1_SMALL1850) --cfg-file ../site_and_regional/default_data_1850.cfg

# Note (slevis): The next line makes the landuse.timeseries from 1850 to 2015, so also run
# modify_smallville.sh to generate three modified landuse.timeseries files needed for testing.
1x1-smallville-transient : FORCE
	$(SUBSETDATA_POINT) --create-surface --create-landuse --surf-year 1850 $(SUBSETDATA_1X1_SMALLTRANSIENT)
	../modify_input_files/modify_smallville.sh

1x1-cidadinho-present : FORCE
	$(SUBSETDATA_POINT) --create-surface $(SUBSETDATA_1X1_SMALL_BR)

#
# Crop with future scenarios
#

crop-global-future : crop-global-SSP1-1.9-f09 \
                     crop-global-SSP1-2.6-f09 \
                     crop-global-SSP2-4.5 \
                     crop-global-SSP2-4.5-low-res \
                     crop-global-SSP3-7.0-f09 \
                     crop-global-SSP4-3.4-f09 \
                     crop-global-SSP4-6.0-f09 \
                     crop-global-SSP5-8.5-f09

crop-global-SSP2-4.5-low-res : crop-global-SSP2-4.5-f10 \
                               crop-global-SSP2-4.5-f45 \
                               crop-global-SSP2-4.5-ne3
crop-global-SSP2-4.5         : crop-global-SSP2-4.5-f09 \
                               crop-global-SSP2-4.5-f19 \
                               crop-global-SSP2-4.5-hcru \
                               crop-global-SSP2-4.5-ne16 \
                               crop-global-SSP2-4.5-ne30 \
                               crop-global-SSP2-4.5-ne0np4 \
                               crop-global-SSP2-4.5-C96 \
                               crop-global-SSP2-4.5-mpasa120

crop-global-SSP1-1.9-f09 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP1-2.6-f09 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-f09 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-f19 : FORCE
	$(MKSURFDATA)  --number-of-nodes 8 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-f10 : FORCE
	$(MKSURFDATA)  --number-of-nodes 8 --tasks-per-node 16 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-f45 : FORCE
	$(MKSURFDATA)  --number-of-nodes 14 --tasks-per-node 32 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-hcru : FORCE
	$(MKSURFDATA)  --number-of-nodes 72 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-ne3 : FORCE
	$(MKSURFDATA)  --number-of-nodes 14 --tasks-per-node 32 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-ne16 : FORCE
	$(MKSURFDATA)  --number-of-nodes 8 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-ne30 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-ne0np4 : FORCE
	$(MKSURFDATA)  --number-of-nodes 2 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-C96 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP2-4.5-mpasa120 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

# Don't do the high-resolution cases by default...
crop-global-SSP2-4.5-hi-res : FORCE
	$(MKSURFDATA)  --number-of-nodes 6 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP3-7.0-f09 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP4-3.4-f09 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP4-6.0-f09 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

crop-global-SSP5-8.5-f09 : FORCE
	$(MKSURFDATA)  --number-of-nodes 9 --scenario $@ --jobscript-file $@.sh --walltime 12:00:00
	$(BATCHJOBS) $@.sh

#
# tower sites: NEON, PLUMBER2
#
neon : FORCE
	../site_and_regional/neon_surf_wrapper

plumber2 : FORCE
	../site_and_regional/plumber2_surf_wrapper

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
# /glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/surfdata_map/README_c141219
urban-alpha : FORCE
	$(SUBSETDATA_POINT_URBAN) --create-surface $(SUBSETDATA_1X1_URBALPHA)
	# Set number of nlevurb dimension to 4 for this site
	ncks -O --dmn nlevurb,0,3 -o $(URBALPHA_TMP2_FNAME) $(URBALPHA_TMP_FNAME)
	$(MODIFYSURF) modify_1x1_urbanc_alpha.cfg -i $(URBALPHA_TMP2_FNAME) -o $(URBALPHA_FNAME)
	$(RM) $(URBALPHA_TMP_FNAME) $(URBALPHA_TMP2_FNAME)


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

