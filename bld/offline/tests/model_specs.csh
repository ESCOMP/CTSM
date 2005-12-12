#!/bin/csh -f

#----------------------------------------------------------------------
# Build the "misc.h" preprocessor token file.
#
# "misc.h" should only have the preprocessor tokens needed for the
# "camclm_share" directory. These are the tokens
# that are shared between the land model and the atmosphere.
#
# SPMD = define if use distributed memory implementation
#-----------------------------------------------------------------------
#
# Based on env variables decide if certain tokens are set or not
#
if ( ($?SPMD) && ($SPMD == "TRUE") ) then
  setenv defSPMD "#define SPMD"
else
  setenv defSPMD "#undef SPMD"
endif

set file = misc.h
set tmpfile = $file.tmp
\cat >! $tmpfile << EOF
#ifndef MISC_SET
#define MISC_SET
$defSPMD
#endif
EOF

if ( ! -e $file )then
  \mv $tmpfile $file
else 
  \cmp $file $tmpfile
  if ( $status ) \mv $tmpfile $file
endif
if ( -e $tmpfile ) \rm $tmpfile

#-----------------------------------------------------------------------
#  Produce preproc.h  preprocessor token file used only by the land model.
#-----------------------------------------------------------------------

set validres = ( T31 T31cn T31casa T31dgvm T31cnall)

if ($RESOLUTION == T31) then
   setenv LSMLON  96
   setenv LSMLAT  48
   setenv MAXPATCH_PFT 4
   set dust     = "#define DUST"
   set voc      = "#define VOC"
   set rtm      = "#define RTM"
   set cn       = "#undef  CN"
   set supln    = "#undef  SUPLN"
   set casa     = "#undef  CASA"
   set dgvm     = "#undef  DGVM"   
else if ($RESOLUTION == T31cn) then
   setenv LSMLON  96
   setenv LSMLAT  48
   setenv MAXPATCH_PFT 4
   set dust     = "#undef  DUST"
   set voc      = "#undef  VOC"
   set rtm      = "#undef  RTM"
   set cn       = "#define CN"
   set supln    = "#define SUPLN"
   set casa     = "#undef  CASA"
   set dgvm     = "#undef  DGVM"   
else if ($RESOLUTION == T31cnall) then
   setenv LSMLON  96
   setenv LSMLAT  48
   setenv MAXPATCH_PFT 17
   set dust     = "#undef  DUST"
   set voc      = "#undef  VOC"
   set rtm      = "#undef  RTM"
   set cn       = "#define CN"
   set supln    = "#undef  SUPLN"
   set casa     = "#undef  CASA"
   set dgvm     = "#undef  DGVM"   
else if ($RESOLUTION == 'T31dgvm') then
   setenv LSMLON  96
   setenv LSMLAT  48
   setenv MAXPATCH_PFT 10
   set dust     = "#undef  DUST"
   set voc      = "#undef  VOC"
   set rtm      = "#undef  RTM"
   set cn       = "#undef  CN"
   set supln    = "#undef  SUPLN"
   set casa     = "#undef  CASA"
   set dgvm     = "#define DGVM"   
else if ($RESOLUTION == 'T31casa') then
   setenv LSMLON  96
   setenv LSMLAT  48
   setenv MAXPATCH_PFT 10
   set dust     = "#undef  DUST"
   set voc      = "#undef  VOC"
   set rtm      = "#undef  RTM"
   set cn       = "#undef  CN"
   set supln    = "#undef  SUPLN"
   set casa     = "#define CASA"
   set dgvm     = "#undef  DGVM"   
else
   echo "Resolution not valid: Valid options are $validres"
   exit 3
endif

set file = preproc.h
set tmpfile = $file.tmp

cat >! $tmpfile << EOF
#ifndef PREPROC_SET
#define PREPROC_SET
#define OFFLINE
#define LSMLON $LSMLON
#define LSMLAT $LSMLAT
#define MAXPATCH_PFT $MAXPATCH_PFT
$dust
$voc
$rtm
$cn
$supln
$casa
$dgvm
#endif
EOF

if ( ! -e $file )then
  \mv $tmpfile $file
else 
  \cmp $file $tmpfile
  if ( $status ) \mv $tmpfile $file
endif
if ( -e $tmpfile ) \rm $tmpfile

#-----------------------------------------------------------------------------
# Generate Filepath
#-----------------------------------------------------------------------------

if ( -f Filepath ) \rm Filepath
touch Filepath

\cat >! Filepath << EOF
$ROOTDIR/src/csm_share/shr
$ROOTDIR/src/main
$ROOTDIR/src/biogeophys
$ROOTDIR/src/biogeochem
$ROOTDIR/src/riverroute
$ROOTDIR/src/utils/timing
$ROOTDIR/src/utils/esmf_wrf_timemgr
EOF
endif






