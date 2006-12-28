#!/usr/bin/env perl
#=======================================================================
#
#  This is a stand-alone CLM2 model test script
#
# Usage:
#
# perl test-model.pl options
#
#=======================================================================
use 5.004;   # Use at least this version of perl
use Cwd;     # Use current working directory module
#
# Use the shortened syntax for the following ENV variables
#
use Env qw(BUILD_DIR CASE CASE_DIR DEBUG DEPGEN EXEDIR 
           EXENAME GNUMAKE LAB LM_DATDIR LOG_DIR MODEL_BLDDIR MODEL_CFGDIR 
           MODEL_CPRDIR MODEL_DATDIR MODEL_EXEDIR ROOTDIR
           NAMELIST RESOLUTION RUNTYPE 
           SCRIPT_DIR SPMD VPATH);
use strict;

# location of the main script 
if ( defined($SCRIPT_DIR) ) {
  chdir( $SCRIPT_DIR );
} else {
  $SCRIPT_DIR = cwd( );
}

# List of where to look for the Perl modules
use lib ".";   
require "CLM_test.pm";

# Setup tests that will be run and create a new test object
my $clm = CLM_test->new;

# Process input arguments
$clm->process_args;
my $platform = $clm->Platform;# Three letter description of platform
my @cprout;                   # Array of cprnc log files

# Set directory locations
$MODEL_CFGDIR = "$SCRIPT_DIR";              # Configure scripts directory
$ROOTDIR      = "$MODEL_CFGDIR/../../..";   # Model root directory
$MODEL_CPRDIR = "$ROOTDIR/tools/newcprnc";  # Comparision program directory

# Keep track of initial setting of SPMD mode
my $SPMD_prev;
if ( defined($SPMD) ) { 
  $SPMD_prev = $SPMD;
} else {
  $SPMD_prev = undef;
}

#---------------------------------------------------------------------------
# Build support programs  
# (cprnc to compare history files and makdep to build dependencies)
#---------------------------------------------------------------------------

$clm->setup_directories( "support" );     # setup the directories to build the support programs
$clm->machine_specs( "support" );         # ENV variables specific to each machine-type
$EXEDIR = $CASE_DIR;                      # Location to build cprnc and makdep

# Make cprnc
$VPATH  = $MODEL_CPRDIR;                  # Location of source for cprnc
$clm->clean( "$MODEL_CPRDIR" );           # Clean the source directory
$clm->chdir( "$BUILD_DIR/cprncobj" );     # Change to build directory for cprnc
$clm->build_support;                      

# Reset the EXEDIR
$EXEDIR = undef;                          

# Mark how far script has proceeded
$clm->continuation_mark( "0_nothing_done" );  

# Determine resolution
if (defined $main::run_res) {
    $RESOLUTION = $main::run_res;
} else {
    $RESOLUTION = "T31";
}
%main::CLMEXP = {};

# variables needed below
my @files;
my $desc;
my %exedir;
my $newfile;

# Change default namelist settings for all tests

$main::CLMEXP{'start_ymd'} = 19971231;
$main::CLMEXP{'start_tod'} = 77400;
$main::CLMEXP{'hist_dov2xy(1)'} = ".true.";
$main::CLMEXP{'hist_dov2xy(2)'} = ".true.";

#---------------------------------------------------------------------------
# Debug tests
# Configure/build the model with debugging on to catch errors
# Compare the history files with and without SPMD-mode
#---------------------------------------------------------------------------

# Change default namelist settings

$main::CLMEXP{'start_ymd'} = 19971231;
$main::CLMEXP{'start_tod'} = 77400;

if ($RESOLUTION =~ /dgvm/) {
  $main::CLMEXP{'hist_nhtfrq(1)'} = -24;
} else {
  $main::CLMEXP{'hist_nhtfrq(1)'} =  1;
}
$main::CLMEXP{'rtm_nsteps'}     = 2;
$main::CLMEXP{'hist_mfilt(1)'}  = 1;
$main::CLMEXP{'hist_ndens(1)'}  = 1;

foreach $desc ( sort( keys(%main::debug_SPMD) ) ) {

    $main::CLMEXP{'nelapse'} = $main::nelapse{$desc};
    print ("test-model: desc is $desc\n");
    print ("test-model: nelapse is $main::CLMEXP{'nelapse'}\n");
    
    # Set case name
    $SPMD = $main::debug_SPMD{$desc}; # SPMD on or off
    $SPMD =~ /(.)/;                   # Extract first character from SPMD
    my $SPMDC = $1;
    $CASE = "deb${SPMDC}${platform}${RESOLUTION}";
    
    # Model build and executable directories (set later)
    $MODEL_BLDDIR = undef;
    $MODEL_EXEDIR = undef;
    
    # Turn on Debug (compiler error checking)
    $DEBUG = "TRUE";

    # Set runtype to test runtype
    $RUNTYPE = $main::runtypes{$desc};
    print ("desc is $desc and runtype is $RUNTYPE\n");
    
    # Change to the script directory
    $clm->chdir( $SCRIPT_DIR );
    
    # Setup the build/run/source directories
    $clm->setup_directories;
    
    # Change to build directory
    $clm->chdir( "$MODEL_BLDDIR" );
    
    # Run the configure script
    $clm->configure;
    
    # Set build environment variables
    $clm->machine_specs( "build" );
    
    # Build the model executable
    $clm->make;
    
    # Change to executable directory
    $clm->chdir( "$MODEL_EXEDIR" );
    
    # Set exedir descriptor
    $exedir{$desc} = "$MODEL_EXEDIR";
    
    # Remove files before run
    if ( -d "HISTOUT" ) {
      @files = <HISTOUT/*>;
      system ("/bin/rm -f @files");
    } else {
      system ("mkdir HISTOUT");
    }
    
    # Set run-time env vars
    $clm->machine_specs( "run" );
    
    # Build the model namelist
    $clm->build_namelist;
    
    # Run the case
    $clm->run( "$desc: Run $RUNTYPE simulation with DEBUG on for # steps=: ".
	       $main::CLMEXP{'nelapse'}, "$LOG_DIR/$desc.log" );
    
    # Move netcdf files to appropriate directory
    @files = <*clm2.h[012]*.nc>;
    foreach (<@files>) {
      $newfile = $_;
      $newfile =~ s#.*\.clm2#clm2#;
      system ("mv $_ HISTOUT/$newfile");
    }

    # Compare all output netcdf files
    if ( $desc eq "02_debug_run_nonSPMD" ) {
      foreach (<HISTOUT/*>) {
	$_ =~ s#.*/##;
	$clm->compare_files( "$exedir{'01_debug_run_SPMD'}/HISTOUT/$_",  	
			     "$exedir{'02_debug_run_nonSPMD'}/HISTOUT/$_",  
			     "$MODEL_EXEDIR/$_.cprout", 
			     "SPMD and non-SPMD mode are are not bit-for-bit" );
      }
    }
    
    # Clean out BLD and EXE directories
    if ( $clm->do_clean ) { $clm->clean; }   
    
    
    # Mark how far script has proceeded
    $clm->continuation_mark( $desc );  
    
    $SPMD = $SPMD_prev;

}

#---------------------------------------------------------------------------
# Restart tests
# (ensure model starts, restarts and restarts are bit-for-bit)
# (change PE's on restart to make sure restarts not PE dependent)
#---------------------------------------------------------------------------

# Set case name
$CASE  = "test$platform${RESOLUTION}"; 

# Modify namelist
$main::CLMEXP{'start_ymd'} = 19971231;
$main::CLMEXP{'start_tod'} = 0;

$main::CLMEXP{'hist_ndens(1)'}  = 1;
$main::CLMEXP{'hist_ndens(2)'}  = 1;
$main::CLMEXP{'hist_nhtfrq(1)'} = 25;
$main::CLMEXP{'hist_nhtfrq(2)'} = 48;
$main::CLMEXP{'hist_mfilt(1)'}  = 1;
$main::CLMEXP{'hist_mfilt(2)'}  = 2;
$main::CLMEXP{'hist_fincl2(1)'} = "\'TV\'";
$main::CLMEXP{'hist_fincl2(2)'} = "\'TSOI\'";
$main::CLMEXP{'hist_fincl2(3)'} = "\'TSA:I\'";
$main::CLMEXP{'hist_fincl2(4)'} = "\'TREFMNAV\'";
$main::CLMEXP{'hist_fincl2(5)'} = "\'TREFMXAV\'";
$main::CLMEXP{'rtm_nsteps'}     = 6;

# Turn off Debug (compiler error checking) 
$DEBUG = "FALSE";               

# Unset the build  and executable directories (set later) 
$MODEL_BLDDIR = undef;           
$MODEL_EXEDIR = undef;          

# Change to script directory
$clm->chdir( $SCRIPT_DIR ); 

# Setup the build/run/source directories
$clm->setup_directories;         

# Change to build directory
$clm->chdir( "$MODEL_BLDDIR" );  

# Run the configure script
$clm->configure;                 

# Set build environment variables
$clm->machine_specs( "build" );  

# Build the model executable
$clm->make;                     

# Change to executable directory
$clm->chdir( "$MODEL_EXEDIR" );  

# Remove history output files before run
if ( -d "INITIAL" ) {
  @files = <INITIAL/*>;
  system ("/bin/rm -f @files");
} else {
  system ("mkdir INITIAL");
}
if ( -d "RESTART" ) {
  @files = <RESTART/*>;
  system ("/bin/rm -f @files");
} else {
  system ("mkdir RESTART");
}

foreach $desc ( sort( keys(%main::runtypes) ) ) {

  if ($desc =~ /0[345]/) {
    
    # Set exedir descriptor
    $exedir{$desc} = "$MODEL_EXEDIR";

    # Determine run type
    $RUNTYPE = $main::runtypes{$desc};
    $main::CLMEXP{'nelapse'} = $main::nelapse{$desc}; 

    # Lower PE's for restart run 
    if ( $desc eq "04_restart" ) { $clm->lower_PEs; }
    
    # Set run-time env vars
    $clm->machine_specs( "run" );  

    # Build the model namelist
    $clm->build_namelist;          

    # Do run
    $clm->run( "$desc: Run $RUNTYPE simulation for # steps=: ".$main::nelapse{$desc}, 
	       "$LOG_DIR/$desc.log" );

    # Move netcdf files to appropriate directory
    @files = <*clm2.h[012]*.nc>;
    foreach (<@files>) {
      $newfile = $_;
      $newfile =~ s#.*\.clm2#clm2#;
      if ($desc eq "05_norestart_compare_to_restart") {system ("mv $_ INITIAL/$newfile")}  
      if ($desc eq "04_restart"                     ) {system ("mv $_ RESTART/$newfile")}  
    }
    
    # Compare all output netcdf files
    if ( $desc eq "05_norestart_compare_to_restart" ) { 
      foreach (<RESTART/*>) {
	$_ =~ s#.*/##;
	$clm->compare_files( "RESTART/$_", 
			     "INITIAL/$_", 
			     "$MODEL_EXEDIR/$_.cprout", 
			     "Restart or change in PE not exact" );
      }
    }
    
    # Reset PE's for initial run
    if ( $desc eq "04_restart" ) { $clm->reset_PEs; }
    
    # Mark how far script has proceeded
    $clm->continuation_mark( $desc );      

  }
}

# Clean out BLD and EXE directories
if ( $clm->do_clean ) { $clm->clean; }  

#---------------------------------------------------------------------------
# Control tests
# Compare to control source directory.
# always compare if files exist, only rerun if they don't.
#---------------------------------------------------------------------------

# Set Case name
$CASE = "cont${platform}${RESOLUTION}";  

# Model build directory and exec directory (set later)
$MODEL_BLDDIR = undef;                   
$MODEL_EXEDIR = undef;                   

# Change to script directory
$clm->chdir( $SCRIPT_DIR ); 

# Setup the build/run/source directories
$clm->setup_directories( "compare" );    

# Change to build directory
$clm->chdir( "$MODEL_BLDDIR" );          

$desc = $main::control;
if ( $clm->do_compare ) {

  # Run the configure script
  $clm->configure;                    

  # Set build type environment variables
  $clm->machine_specs( "build" );     

  # Build the model executable
  $clm->make;                         

  # Change to executable directory
  $clm->chdir( "$MODEL_EXEDIR" );     

  # Set exedir descriptor
  $exedir{$desc} = "$MODEL_EXEDIR";

  # Remove files before run
  if ( -d "HISTOUT" ) {
    @files = <HISTOUT/*>;
    system ("/bin/rm -f @files");
  } else {
    system ("mkdir HISTOUT");
  }

  # Set run-time environment variables
  $clm->machine_specs( "run" );       

  # Build the model namelist
  $clm->build_namelist;               

  # Do the run
  $clm->run( "$desc: Run $main::runtypes{$desc} simulation with control library for # steps=: ". $main::nelapse{$desc},
	     "$LOG_DIR/$desc.log" );

  # Move netcdf files to appropriate directory
  @files = <*clm2.h[012]*.nc>;
  foreach (<@files>) {
    $newfile = $_;
    $newfile =~ s#.*\.clm2#clm2#;
    system ("mv $_ HISTOUT/$newfile");
  }

  # Compare all output history netcdf files
  foreach (<HISTOUT/*>) {
    $_ =~ s#.*/##;
    $clm->compare_files( "$exedir{'05_norestart_compare_to_restart'}/INITIAL/$_",  	
			 "HISTOUT/$_",  
			 "$MODEL_EXEDIR/$_.cprout", 
			 "comparison with control is not bit-for-bit" );
    push @cprout, "$SCRIPT_DIR/$CASE.$desc.$_.cprout";
  }

  # Clean out BLD and EXE directories
  if ( $clm->do_clean ) { $clm->clean; }   

}

# Mark how far script has proceeded
$clm->continuation_mark( $desc );  

# Clean out BLD and EXE directories
if ( $clm->do_clean ) { $clm->clean; }  

if ( ! $clm->do_nofail ) {
  print "\n\nModel testing was successful!\n\n";
}

# Unset environment variables
$clm->unset_env;

# Report on success of comparisions
$clm->report_compare( @cprout );    

# Mark that script is completed
$clm->continuation_mark( "done" );  


