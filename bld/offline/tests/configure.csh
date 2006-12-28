#!/bin/csh -f

#-----------------------------------------------------------------------
# Configure the CLM model to build an executable on a specific platform
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# 1. Set load and runtime environment variables
#-----------------------------------------------------------------------

if ( ! -d $MODEL_BLDDIR ) mkdir $MODEL_BLDDIR
if ( ! -d $MODEL_EXEDIR ) mkdir $MODEL_EXEDIR
echo "Build directory: $MODEL_BLDDIR"
cd $MODEL_BLDDIR
echo "Get machine specific ENV variables";
source $MODEL_CFGDIR/config_machine_specs.csh
if ( $status != 0 ) exit

#-----------------------------------------------------------------------
# 2. Generate misc.h, misc.h and Filepath files
#-----------------------------------------------------------------------

echo "Make misc.h, misc.h and Filepath";
source $MODEL_CFGDIR/model_specs.csh
if ( $status != 0 ) exit

#-----------------------------------------------------------------------
# 3. Link Makefile to build directory
#-----------------------------------------------------------------------

set makefile = "Makefile"
if ( `pwd` != $MODEL_CFGDIR )then
  \rm $makefile
  \ln -sf $MODEL_CFGDIR/$makefile .
endif

#-----------------------------------------------------------------------
# 4. Write out environment variables that have been set 
#-----------------------------------------------------------------------

if ( $?ROOTDIR ) then
  echo "$ROOTDIR" >! Rootdir
endif

#-----------------------------------------------------------------------
# 5. Write out environment variables that have been set 
#-----------------------------------------------------------------------

set ARCH = `uname -s`
if ( $ARCH == UNICOS/mp ) set ARCH = UNICOSmp
if ( $?SCRIPT_DIR ) then
  set configfile = "$SCRIPT_DIR/config_env_$ARCH.csh"
else
  set configfile = "config_env_$ARCH.csh"
endif
if ( -f $configfile ) \rm $configfile

echo "Save config-file env settings to $configfile"
cat << EOF > $configfile
#!/bin/csh -f
#------------------------------------------------------------------------
# Set env variables for a model simulation set by 
# $0 
# Ran on `date`
#------------------------------------------------------------------------
EOF
if ( $?MODEL_CFGDIR )then
   echo "setenv MODEL_CFGDIR $MODEL_CFGDIR" >> $configfile
endif
if ( $?MODEL_BLDDIR ) then
  echo "setenv MODEL_BLDDIR $MODEL_BLDDIR" >> $configfile
endif
if ( $?MODEL_EXEDIR ) then
  echo "setenv MODEL_EXEDIR $MODEL_EXEDIR" >> $configfile
endif
if ( $?CASE ) then
  echo "setenv CASE $CASE" >> $configfile
endif
if ( $?ROOTDIR ) then
  echo "setenv ROOTDIR $ROOTDIR" >> $configfile
endif
if ( $?RESOLUTION ) then
  echo "setenv RESOLUTION $RESOLUTION" >> $configfile
endif
if ( $?LSMLAT ) then
  echo "setenv LSMLAT $LSMLAT" >> $configfile
endif
if ( $?LSMLON ) then
  echo "setenv LSMLON $LSMLON" >> $configfile
endif
if ( $?SPMD ) then
  echo "setenv SPMD $SPMD" >> $configfile
endif
echo "#CONFIGURE-END" >> $configfile
chmod +x $configfile

#-----------------------------------------------------------------------
# 6. Log
#-----------------------------------------------------------------------

if ( $?SCRIPT_DIR ) then
  set configlog = "$SCRIPT_DIR/config_$ARCH.log"
else
  set configlog = "config_$ARCH.log"
endif
\rm  $configlog
touch $configlog
echo "Filepath: " >>! $configlog
cat $MODEL_BLDDIR/Filepath >>! $configlog
echo "Rootdir: "           >>! $configlog
echo "preprocessor files:" >>! $configlog
cat $MODEL_BLDDIR/*.h      >>! $configlog
echo "Saved config output to $configlog"

# Return a -1 on success
echo "Done with configure";





