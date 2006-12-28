#============================================================================
#
#	CLM_run.pm   Mariana Vertenstein (original version by Erik Kluzek)
#
#	Perl5 Object module to handle the building and running
#	of the Community Land Model (CLM2).
#
#	The control is kept by setting of ENV variables as well
#	as CLMEXP associative arrays which control the namelist. 
#       By adding extra keyword value pairs in	the main script 
#       to these arrays new namelist items can be added to the namelists.
#
#	This module extends the basic object "CLM" with important
#	useful functionality. Most likely methods that would need
#	to be edited would be in this module rather than the "CLM"
#	module.
#
#	Description of methods:
#
#	Important "hidden" methods:
#
#	run_time_env ------- Set the Env variables needed for run-time.
#
#	Basic methods to build/run the model.
#
#	setup_directories -- Setup the needed build and run directories and Env variables
#	machine_specs ------ Set the Env variables needed for this specific O/S.
#	configure ---------- Run the configure.csh to configure the model simulation.
#	make --------------- Build the model executable.
#	run ---------------- Run the model simulation.
#
#============================================================================

use 5.004;   # Use at least this version of perl
use strict;
use Cwd;

package CLM_run;
@CLM_run::ISA = "CLM";
#
# Use the shortened syntax for the following ENV variables
#
use Env qw(BUILD_DIR CASE CASE_DIR DEBUG EXENAME GNUMAKE
           INC_MPI INC_NETCDF LAB LIB_MPI LIB_NETCDF LOG_DIR LOGNAME 
           MODEL_BLDDIR MODEL_CFGDIR MODEL_DATDIR 
           MODEL_EXEDIR ROOTDIR NAMELIST  
           RUNTYPE SHMEM_CPUS SPMD SPMD_CPUS_ON_NODE SPMD_NODES SPMD_RUNCMND);
use CLM;
use CLM_namelist;
use CLM_lab;

#============================================================================

sub machine_specs {
#
# Set the machine specific ENV variables
# For both building and running.
#
  my $self = shift;
  my $type = shift;

  # Get run-time ENV variables set for each machine type
  my $OS = $self->{'OS'};

  if ( ! defined($type) || $type eq "build" || $type eq "support" ) {
    print "machine_specs:: Get the build ENV variables needed for this platform\n";
    if ( ! defined($GNUMAKE) ) { $GNUMAKE = $CLM_run::GNUMAKE_default->value($LAB, $OS); }
    if ( ! defined($INC_NETCDF) ) { $INC_NETCDF = $CLM_run::INC_NETCDF_default->value($LAB, $OS); }
    if ( ! defined($LIB_NETCDF) ) { $LIB_NETCDF = $CLM_run::LIB_NETCDF_default->value($LAB, $OS); }
    if ( ! defined($INC_MPI) ) { $INC_MPI = $CLM_run::INC_MPI_default->value($LAB, $OS); }
    if ( ! defined($LIB_MPI) ) { $LIB_MPI = $CLM_run::LIB_MPI_default->value($LAB, $OS); }
    if ( ! defined($MODEL_DATDIR) ) { $MODEL_DATDIR = $CLM_run::MODEL_DATDIR_default->value($LAB, $OS); }
    if ( ! defined($EXENAME) ) { $EXENAME = "clm"; }
  }
  if ( ! defined($type) || $type eq "build" ) {
    $self->append_config_env( "build", ("SPMD",  "GNUMAKE", "INC_NETCDF", "LIB_NETCDF",
					"INC_MPI", "LIB_MPI", "EXENAME") );
  }
  if ( ! defined($type) || $type eq "run" ) {
    print "machine_specs:: Get the run-time ENV variables needed for this platform\n";
    if ( ! defined($RUNTYPE) ) {
      my $run = $self->{'RUNTYPE'};
      my @validrun = @$run;
      $RUNTYPE =  $validrun[0];
      print "Enter simulation type: Options are @validrun [$RUNTYPE]";
      $_= <>; chomp( $_ ); if ( /./ ) { $RUNTYPE = $_; }
    }
    if ( ! defined($NAMELIST) ) {
      $NAMELIST = "lnd.stdin";
      print "Set namelist to $NAMELIST in " . cwd( ) . "\n";
    }
    if ( ! defined($MODEL_DATDIR) ) {
      $MODEL_DATDIR = "/fs/cgd/csm/inputdata";
      print "Enter model data directory: [$MODEL_DATDIR]";
      $_ = <>; chomp( $_ ); if ( /./ ) { $MODEL_DATDIR = $_; }
    }
    if ( ! defined($NAMELIST) ) {
      $NAMELIST = "lnd.stdin";
      print "Set namelist to $NAMELIST in " . cwd( ) . "\n";
    }
    if ( ! defined($SHMEM_CPUS) ) { $SHMEM_CPUS = $CLM_run::SHMEM_CPUS_default->value($LAB, $OS); }
    if ( ! defined($SPMD_NODES) ) { $SPMD_NODES = $CLM_run::SPMD_NODES_default->value($LAB, $OS); }
    # If either set to zero, change to undefined so that defaults will be used
    if ( $SHMEM_CPUS == 0 ) { $SHMEM_CPUS = undef; }
    if ( $SPMD_NODES == 0 ) { $SPMD_NODES = undef; }
    $self->run_time_env;
  }
}

#============================================================================

sub run_time_env {
#
# Run-time Platform dependencies...
#
  my $self = shift;

  my $os = $self->{'OS'};
  my @Env_Vars;
  SWITCH: {
#---------------------------------------------------------------------------------------
# Solaris
#---------------------------------------------------------------------------------------
    ($os =~ /solaris/) && do {
                  #
                  # Prepare hostfilename of valid hosts to use
                  # Must agree with filename used in SPMD_RUNCMND (check CLM_lab.pm)
                  #
                  my $filename = "machine";
                  open( MACH, ">$filename" ) || die "ERROR:: Can not open file: $filename\n";
                  print MACH `hostname`;
                  close( MACH );
                  system( "ulimit -s unlimited" );
               last SWITCH;
             };
    ($os =~ /irix/) && do {
#---------------------------------------------------------------------------------------
# SGI
#---------------------------------------------------------------------------------------
# OMP_DYNAMIC       : False => Don't give up processors even if machine is busy
# _DSM_PLACEMENT    : Ensure efficient memory allocation 
# _DSM_WAIT         : Assume dedicated access to processors (SPIN)
# _DSM_VERBOSE      : Print diagnostic info about system-level software
# MPC_GANG          : Gang scheduling.  Does not work well so turn off
# MP_SLAVE_STACKSIZE: Stack size (bytes) to be used by slave processes
# TRAP_FPE          : Abort on overflow or divide by zero
#---------------------------------------------------------------------------------------
                  if ( ! defined($ENV{'OMP_DYNAMIC'}) ) { $ENV{'OMP_DYNAMIC'} = "FALSE"; }
                  if ( ! defined($ENV{'_DSM_PLACEMENT'}) ) { $ENV{'_DSM_PLACEMENT'} = 
                     "ROUND_ROBIN"; }
                  if ( ! defined($ENV{'_DSM_WAIT'}) ) { $ENV{'_DSM_WAIT'} = "SPIN"; }
                  if ( ! defined($ENV{'MPC_GANG'}) ) { $ENV{'MPC_GANG'} = "OFF"; }
                  if ( ! defined($ENV{'MP_SLAVE_STACKSIZE'}) ) { $ENV{'MP_SLAVE_STACKSIZE'} = 
                       "40000000"; }
                  if ( ! defined($ENV{'TRAP_FPE'}) ) { $ENV{'TRAP_FPE'} = 
                       "UNDERFL=FLUSH_ZERO; OVERFL=ABORT,TRACE; DIVZERO=ABORT,TRACE"; }
                  system( "ulimit -s unlimited" );
                  @Env_Vars = ("TRAP_FPE", "MP_SLAVE_STACKSIZE",
                               "MPC_GANG", "_DSM_WAIT", "_DSM_PLACEMENT", "OMP_DYNAMIC" );
               last SWITCH;
             };
    ($os =~ /linux/) && do {
#---------------------------------------------------------------------------------------
# Linux
#---------------------------------------------------------------------------------------
# MPSTKZ  sets slave stack size
#---------------------------------------------------------------------------------------
                  if ( ! defined($ENV{'MPSTKZ'}) ) { $ENV{'MPSTKZ'} = "128M"; }
                  @Env_Vars = ("MPSTKZ" );
               last SWITCH;
             };
    ($os =~ /aix/) && do {
#---------------------------------------------------------------------------------------
# IBM
#---------------------------------------------------------------------------------------
# XLSMPOPTS             the largest amount of space (bytes) that a thread's 
#                       stack will need 
# MP_EUILIB us          message passing subsystem implementation
# MP_NODES              number of physical nodes to run the parallel tasks
#                       (a task refers specifically to an MPI process)
# MP_TASKS_PER_NODE     number of tasks to run on each of the physical nodes 
# MP_RMPOOL             specifies number of a LoadLeveler pool     
#---------------------------------------------------------------------------------------
                  if ( ! defined($ENV{'XLSMPOPTS'}) ) { $ENV{'XLSMPOPTS'} = "stack=40000000"; }
                  if ( ! defined($ENV{'MP_EUILIB'}) ) { $ENV{'MP_EUILIB'} = "us"; }
                  if ( ! defined($ENV{'MP_TASKS_PER_NODE'}) ) { $ENV{'MP_TASKS_PER_NODE'} = 1; }
                  if ( ! defined($ENV{'MP_RMPOOL'}) ) { $ENV{'MP_RMPOOL'} = 1; }
                  if ( ! defined($ENV{'MP_LABELIO'}) ) { $ENV{'MP_LABELIO'} = "yes"; }
                  if ( ! defined($ENV{'MP_STDOUTMOD'}) ) { $ENV{'MP_STDOUTMODE'} = "ordered"; }
                  $ENV{'MP_NODES'} = $SPMD_NODES;
                  print "MP_NODES = " . $ENV{'MP_NODES'} . "\n";
                  @Env_Vars = ("MP_RMPOOL", "MP_NODES",
                                             "MP_TASKS_PER_NODE", "MP_EUILIB", "XLSMPOPTS" );
               last SWITCH;
             };
    ($os =~ /dec_osf/) && do {
#---------------------------------------------------------------------------------------
# Compaq
#---------------------------------------------------------------------------------------
# MP_STACK_SIZE  slave stack size
#---------------------------------------------------------------------------------------
                  if ( ! defined($ENV{'MP_STACK_SIZE'}) ) { $ENV{'MP_STACK_SIZE'} = "17000000"; }
                  if ( ! defined($SPMD_CPUS_ON_NODE) ) { $SPMD_CPUS_ON_NODE = 1; }
                  system( "ulimit -s unlimited" );
                  @Env_Vars = ("SPMD_CPUS_ON_NODE", "MP_STACK_SIZE");
               last SWITCH;
             };
#---------------------------------------------------------------------------------------
# default
#---------------------------------------------------------------------------------------
  }
  if ( defined($SHMEM_CPUS) ) { 
    $ENV{'OMP_NUM_THREADS'} = $SHMEM_CPUS; 
    print "OMP_NUM_THREADS = $ENV{'OMP_NUM_THREADS'}\n";
    push( @Env_Vars, "OMP_NUM_THREADS" );
  }
  if ( ! defined($SPMD_RUNCMND) ) { $SPMD_RUNCMND = $CLM_run::SPMD_RUNCMND_default->value($LAB, $os); }
  push( @Env_Vars, ("SPMD_RUNCMND", "SPMD", "RUNTYPE", "NAMELIST", "EXENAME") );
  $self->append_config_env( "run", @Env_Vars );
}

#============================================================================

sub configure {
#
# Run the configure script
#
  my $self = shift;
  if ( ! -f "$MODEL_CFGDIR/configure.csh" ) {
    die "Error:: Configure file: $MODEL_CFGDIR/configure.csh not"
                . " available -- MODEL_CFGDIR not correct?";
  }
  print "$MODEL_CFGDIR/configure.csh\n";
  system( "$MODEL_CFGDIR/configure.csh" );    
  if ( $? != 0 ) { die "Error in configure"; }
  $self->config_env();
}

#============================================================================

sub make {
#
# Build the model using GNU make
#
  my $self = shift;

#-----------------------------------------------------------------------
# Build the executable
#-----------------------------------------------------------------------

  print "make:: Make the model executable\n";
  my $log = "$LOG_DIR/compile_log.clm";
  print "Compiling CLM ... see $log for log\n\n";

  open( LOG, ">$log" ) || die "Can not add output to the compile log file: $log\n";
  print LOG << "EOF";
-------------------------------------
`date`
-------------------------------------
EOF
  close( LOG );
  #
  $self->exec( "$GNUMAKE  >> $log 2>&1", "Compile failed look at $log for failure" );
  print "Compile successful!\n\n";
}

#============================================================================

sub setup_directories {
#
# Setup the directories
#
  my $self = shift;

  print "setup_directories:: Setup the directories for building and running the model:\n";
  print "case name and title: $CASE\n";
  my $OS = $self->{'OS'};
  #
  # Case directory (directory where output from different cases are stored)
  #
  if ( ! defined($CASE_DIR) ) {
    $CASE_DIR = $CLM_run::CASE_DIR_default->value( $LAB, $OS ); 
  }
  print "dir: $CASE_DIR lab: $LAB, OS: $OS\n";
  if ( ! -d $CASE_DIR ) { 
     mkdir( $CASE_DIR, 0755 ) || die "Can not mkdir $CASE_DIR"; 
  }
  print "Case directory: $CASE_DIR\n";
  #
  # Config directory
  #
  if ( ! defined($MODEL_CFGDIR) ) {
     $MODEL_CFGDIR = cwd( ); 
  }
  if ( ! -d $MODEL_CFGDIR ) { 
    die "Error:: Configuration file directory: $MODEL_CFGDIR does not exist";
  }
  #
  # Source directory
  #
  if ( ! defined($ROOTDIR) ) {
    $ROOTDIR = cwd( ) . "/../../.."; 
  }
  if ( ! -d $ROOTDIR ) { 
    die "Error:: CLM root directory: $ROOTDIR does not exist";
  }
  print "Source directory: $ROOTDIR/src\n";
  #
  # Execution directory (Need CASE_DIR and CASE set to get it)
  #
  # Case name
  #
  if ( ! defined($CASE) ) {
    die "Error:: Env variable CASE not defined";
  }
  #
  # Now set EXEDIR
  #
  if ( ! defined($MODEL_EXEDIR) ) {
    $MODEL_EXEDIR = "$CASE_DIR/$CASE";
  }
  print "Exec directory: $MODEL_EXEDIR\n";
  #
  # Build directory
  #
  if ( ! defined($BUILD_DIR) ) {
    $BUILD_DIR = $CLM_run::BUILD_DIR_default->value( $LAB, $OS );
  }
  if ( ! defined($MODEL_BLDDIR) ) {
    $MODEL_BLDDIR = $CLM_run::MODEL_BLDDIR_default->value( $LAB, $OS );
  }
  #
  # Directory where log-files go
  #
  $LOG_DIR = $CLM_run::LOG_DIR_default->value( $LAB, $OS ); 
  print "dir: $LOG_DIR lab: $LAB, OS: $OS\n";
  if ( ! -d $LOG_DIR ) { 
     mkdir( $LOG_DIR, 0755 ) || die "Can not mkdir $LOG_DIR"; 
  }
  print "Log-file directory: $LOG_DIR\n";
  #
  # Clean blddir and exedir if requested
  #
  if ( $self->do_clean ) {
    $self->clean;
  }
  #
  # Create exedir and blddir
  #
  if ( ! -d $MODEL_EXEDIR ) { 
    mkdir( $MODEL_EXEDIR, 0755 ) || die "Can not mkdir $MODEL_EXEDIR"; 
  }
  if ( ! -d $MODEL_BLDDIR ) { 
    mkdir( $MODEL_BLDDIR, 0755 ) || die "Can not mkdir $MODEL_BLDDIR"; 
  }
  print "Build directory: $MODEL_BLDDIR\n";
}

#============================================================================

sub run {
#
# Run the model
#
  my $self = shift;
  my $desc = shift;
  my $logfile = shift;

  print "run:: Run the model\n";
  if ( defined($desc) ) { print "$desc\n"; }
  #
  # Check that in $MODEL_EXEDIR
  #
  my $pwd = cwd;
  $self->checkdir( $MODEL_EXEDIR, 
     "Error: Not calling run from within \$MODEL_EXEDIR($MODEL_EXEDIR)" );
  #
  # Check that Namelist exists
  #
  if ( ! -f $NAMELIST ) { die "Error: Namelist $NAMELIST does not exist"; }
  if ( -z $NAMELIST ) { die "Error: Namelist $NAMELIST is empty!"; }

  my $exec = "./$EXENAME";
  if ( (ref($self) ne "CLM_script") && (! -f $exec) ) {
    die "Error: Executable ($exec) does not exist!";
  }
  my $OS = $self->{'OS'};
  my $log; 
  #
  # Prepare log file, put namelist and ChangeLog at top of log
  #
  if ( $self->do_log ) {
    if ( ! defined($logfile) ) {
      my $logid = `/bin/date +%y%m%d-%H%M%S`; chomp( $logid );
      $logfile = "$LOG_DIR/clm.$CASE.log.$logid";
    }
    $self->{'LOGFILE'} = $logfile;
    $log = ">> $logfile 2>&1";
    if ( -f $logfile ) { unlink( $logfile ); }
    open( LOG, ">$logfile" ) || die "Can not open logfile $logfile";
    open( NAMELIST , "<$NAMELIST" ) || die "Can not open namelist $NAMELIST";
    print LOG "Namelist:\n";
    while( defined($_ = <NAMELIST>) ) { print LOG $_; }
    close( NAMELIST );
    my $change = "$ROOTDIR/doc/ChangeLog";
    print "Change: $change\n";
    if ( -f $change ) {
      open( CHANGELOG , "<$change" ) || die "Can not open ChangeLog $change";
      print LOG "ChangeLog:\n";
      <CHANGELOG>;
      while( defined($_ = <CHANGELOG>) && (!/====/) ) { print LOG $_; }
      close( CHANGELOG );
    }
    close( LOG );
  #
  # If no-log file option used
  #
  } else {
    $log = " ";
    $logfile = " ";
  }
  print "Running CLM, $RUNTYPE... directory " . cwd() . "\n\n";

  #
  # Actually run the model
  #
  my $die_msg = "Error running the model";
  if ( $self->do_log ) { $die_msg = "$die_msg - see $logfile for log \n\n"; }
  if ( $SPMD ne "TRUE" ) {
    $self->exec( "time $exec $log", $die_msg, "echo" );
  } elsif ( $SPMD_RUNCMND =~ /sh -c/ ) {
    $self->exec( "time $SPMD_RUNCMND \'$exec $log\'", $die_msg, "echo" );
  } else {
    $self->exec( "time $SPMD_RUNCMND $exec $log", $die_msg, "echo" );
  }
  print "\n\nCLM2 Finished";
  if ( $self->do_log ) { print " - see $logfile for log \n\n"; }
  print "\n\n";

}

1   # To make use or require happy
