#
#	CLM_test.pm			Erik Kluzek (modified by Mariana Vertenstein)
#
#	Perl5 Object module extending the CLM_run.pm module
#	to handle testing of the Community Atmospheric Model (CLM).
#
#	Methods:
#
#	new ---------------- Constructor
#	process_args ------- Process main script command line arguments.
#	list_tests --------- List the tests performed by the test suite.
#	usage -------------- Terminate and return proper syntax for arguments.
#	setup_directories -- Setup directories.
#	exec --------------- Extend CLM_run exec method so that commands aren't
#			     executed if you are continuing from a previous run and
#			     haven't caught up to the point where you started.
#	do_compare --------- Return true if you want to compare to another program library.
#	do_continue -------- Return true if the continue option is set.
#	do_nofail ---------- Return true is the nofail option is set.
#	numerically -------- Use with "sort" to sort values numerically.
#	completed ---------- Return true if this part of the script was already completed.
#	compare_files ------ Compare history files of two simulations (uses cprnc)
#	report_compare ----- Summary report on success of file comparisions.
#	mark --------------- Mark to describe where at in process.
#	start_log ---------- Initialize the test-log file.
#	read_mark ---------- Read the continuation mark.
#	write_mark --------- Write the continuation mark.
#	continuation_mark -- Mark how far the script has gone.
#	lower_PEs ---------- Lower the number of CPU's being used (as a test).
#	reset_PEs ---------- Reset the number of CPU's being used.
#	die ---------------- Terminate script or continue if no-fail option set.
#
#	$Id$
#
use 5.004;   # Use at least this version of perl
use strict;
use Cwd;

package CLM_test;
@CLM_test::ISA = "CLM_run";
use CLM_run;
#
# Use the shortened syntax for the following Env variables
#
use Env qw(BUILD_DIR CASE CASE_DIR CONT_ROOTDIR LAB LOGNAME 
           MODEL_CPRDIR MODEL_BLDDIR MODEL_EXEDIR ROOTDIR
           RESOLUTION SCRIPT_DIR SHMEM_CPUS SPMD SPMD_NODES SPMD_RUNCMND);
#
# Use all of the CLM_run methods and data
# Override the constructor, process_args, and usage methods.
#

sub new {
#
# Constructor
#
  my $class = shift;

  my $self = $class->SUPER::new;
  my $OS = $self->{'OS'};
  $self->{'COMPARE'} = "no";          # Compare to control source
  $self->{'ROOTDIR'}  = undef;        # Original ROOTDIR if it changes
  $self->{'CLM_CONT_DIR'} = undef;    # Directory with CLM model source to compare
  $self->{'CONTINUE'} = "no";         # If you want to continue a previous script run
  $self->{'COMPLETED'} = "yes";       # If this part has already been completed or not
  $self->{'NOFAIL'} = "no";           # Don't die on errors, write log and continue
  $self->{'STATUS'} = "ran";          # Status of test (ran or FAIL)
  $self->{'continue_file'} = undef;   # filename that keeps track of how far script has gone
  $self->{'TESTLOG'} = "$SCRIPT_DIR/test.$OS.log"; # Log of results
  $self->{'LASTTEST'} = undef;        # Last test to run
  $self->{'BEGIN'} = "yes";           # Mark that script has just begun
  $self->{'REMOTE_MACHINE'}  = undef; # Remote machine to use for error-growth
  $self->{'REMOTE_PLATFORM'} = undef; # Remote platform to use for error-growth
  bless( $self, $class );
  $self->{'continue_file'} = "script_restart_" . $self->arch;
  return( $self );
}

sub process_args {
#
# Process the input arguments to the script
#
  my $self = shift;

  print "process_args:: Process the input arguments\n";
  while ( defined($ARGV[0]) ) {
    $_ = shift( @ARGV );
    # Help message
    if ( /^-h/ ) {              
      $self->usage;
    # List tests
    } elsif ( /^-t/ ) {
      $self->list_tests;
    # Clean the old directories out before running
    } elsif ( /^-clean/ ) {
      $self->do_clean( "yes" );
    # If you want to continue a previous run of the script
    } elsif ( /^-res/ ) {
      $_ = shift( @ARGV );
      if ( ! defined($_) ) { $self->usage; }
      $main::run_res = $_;
      print ("run_res is $main::run_res\n"); 
    # If you want to skip to a given test (or just a set of tests) (or pick a resolution)
    } elsif ( /^-s/ ) {
      $_ = shift( @ARGV );
      if ( ! defined($_) ) { $self->usage; }
      if ( ! /([0-9]+)-*([0-9]*)/ ) { $self->usage; }
      my $test = $1; 
      # Check that starting test is not the second or higher test 
      $test = $test - 1;
      $self->{'LASTTEST'} = $2;
      if ( $self->{'LASTTEST'} eq "" ) { 
        $self->{'LASTTEST'} = undef; 
      } else {
        if ( $self->{'LASTTEST'} <= $test ) { $self->usage; }
      }
      $self->do_continue( "yes" );
      $self->write_mark( "${test}_" );
    # If you want to set the no fail option and continue even if something goes wrong
    } elsif ( /^-nofail/ ) {
      $self->do_nofail( "yes" );
    # Compare to this model version
    } elsif ( /^-c/ ) {
      $_ = shift( @ARGV );
      if ( ! defined($_) ) { $self->usage; }
      $self->do_compare( "yes" ); 
      $self->{'CLM_CONT_DIR'} = $_;
      $CONT_ROOTDIR = $_;
    # Set the Lab running at
    } elsif ( /^-l/ ) {              
      $_ = shift( @ARGV );
      if ( ! defined($_) ) { $self->usage; }
      $LAB = $_;     
    } else { 
      print "Invalid option $_\n\n";
      $self->usage;
    }
  }

  # Setup tests
  $self->setup_tests;

  # If Env variable $CONT_ROOTDIR defined use it for comparision
  if ( ! $self->do_compare && defined($CONT_ROOTDIR) ) {
    $self->do_compare( "yes" ); 
    $self->{'CLM_CONT_DIR'} = $CONT_ROOTDIR;
  }
  $self->do_build( "yes" ); 
  $self->do_build_namelist( "yes" ); 
  $self->do_log( "yes" ); 
  $self->unset_env;

  print "Build and run test simulations\n";
}

sub usage {
#
# Usage statement if command arguments are done correctly
#
  my $self = shift;

  my @os_list = $self->OS_list;
  print <<EOF;
Usage: perl $0 [options]

Options are:

	-h      = Help (this message)
	-t      = List the tests that are performed.
	-clean  = Clean the old directories out
	-nofail = Continue even if errors are found
	-res    = Resolution (T31 T31cn T31casa T31dgvm T31cnall)
	-s n    = Skip to given test n (or range of numbers)
                  (example -s 4 start with test no. 4)
	 	  (or      -s 2-4 do tests 2 through 4)
	-c dir  = Compare to given version of the model in this directory
                  (example -c /home/user/clm2/src)
	-l lab  = Set the lab you are running at 

EOF
  die "Terminating";
}

sub setup_tests {

  # Setup the tests that will be run
  # If you change the tests here you also may need to change test-model.pl
  # Tests to run and run-type for each, also give the number of steps to run for each
  # (The tests will run in alphabetical order, so include a letter to order it)
  # Save/restart tests should be timed to maximize the number of restart files
  # that need to be saved. Hence, odd times should be used, so that both abs/ems
  # restart datasets are saved as well as history restart datasets.
  # (If you change the keys you need to change the script later as well)

  my $self = shift;

  %main::debug_SPMD = ( '01_debug_run_SPMD'    => "TRUE",
			'02_debug_run_nonSPMD' => "FALSE" );

  %main::runtypes   = ( '01_debug_run_SPMD'              => "initial",
			'02_debug_run_nonSPMD'           => "initial",
			'03_start'                       => "initial",
			'04_restart'                     => "restart",
			'05_norestart_compare_to_restart'=> "initial",
			'06_control'                     => "initial");

  %main::nelapse    = ( '01_debug_run_SPMD'              => 8,
			'02_debug_run_nonSPMD'           => 8,
			'03_start'                       => 10,
			'04_restart'                     => 86,
			'05_norestart_compare_to_restart'=> 96,
			'06_control'                     => 96);

  if ($main::run_res =~ /DGVM/) {

      %main::runtypes   = ( '01_debug_run_SPMD'              => "branch" ,
			    '02_debug_run_nonSPMD'           => "branch" ,
			    '03_start'                       => "branch" ,
			    '04_restart'                     => "restart",
			    '05_norestart_compare_to_restart'=> "branch" ,
			    '06_control'                     => "branch" );

      %main::nelapse    = ( '01_debug_run_SPMD'              => 346,
			    '02_debug_run_nonSPMD'           => 346,
			    '03_start'                       => 100,
			    '04_restart'                     => 246,
			    '05_norestart_compare_to_restart'=> 346,
                            '06_control'                     => 346);

  } 

  $main::control= "06_control";

}

sub list_tests {
#
# List the tests that are performed.
#
  my $self = shift;

  my @list, sort( keys( %main::runtypes  ) );
  push @list, $main::control;

  print <<EOF;
  List of tests that are performed:

EOF
  foreach my $i ( @list ) {
    print "$i\n";
  }
  die "\n\nTerminating";
}

sub setup_directories {
#
# Extend the set_directories method to handle comparing to previous
# versions.
#
  my $self = shift;
  my $comp = shift;

  if ( defined($comp) && ($comp eq "compare") ) {
    print "Setup directories for comparision: ";
    if ( defined($self->{'CLM_CONT_DIR'}) ) {
      #
      # Change ROOTDIR to comparision directory, change back when continuation_mark called
      #
      $self->{'ROOTDIR'} = $ROOTDIR;
      $ROOTDIR = $self->{'CLM_CONT_DIR'};
      print "compare to $ROOTDIR\n";
    } else {
      print "\n";
    }
    $MODEL_BLDDIR = undef;
    $MODEL_EXEDIR = undef;
    $self->SUPER::setup_directories;
  } elsif ( defined($comp) && ($comp eq "support") ) {
    print "Setup support directories: \n";
    #
    # Build and Case directory (where to build code and run cases)
    #
    my $OS = $self->OS;
    if ( ! defined($CASE_DIR) ) {
      $CASE_DIR = $CLM_run::CASE_DIR_default->value( $LAB, $OS );
    }
    if ( ! -d $CASE_DIR ) {
       mkdir( $CASE_DIR, 0755 ) || die "Can not mkdir $CASE_DIR";
    }
    print "Case directory: $CASE_DIR\n";
    if ( ! defined($BUILD_DIR) ) {
      $BUILD_DIR = $CLM_run::BUILD_DIR_default->value( $LAB, $OS );
    }
    if ( ! -d $BUILD_DIR ) {
      mkdir( $BUILD_DIR, 0755 ) || die "Can not mkdir $BUILD_DIR";
    }
    my $dir = "$BUILD_DIR/cprncobj";
    if ( ! -d $dir ) {
      mkdir( $dir, 0755 ) || die "Can not mkdir $dir";
    }
    system( "/bin/rm $dir/Makefile" );
    symlink( "$MODEL_CPRDIR/Makefile", "$dir/Makefile" );
  } else {
    $self->SUPER::setup_directories;
  }
}

sub do_nofail {
#
#  Test if the no fail option is set
#
  my $self = shift;
  my $value = shift;

  return( $self->do('NOFAIL', $value) );
}

sub do_compare {
#
#  Test if should compare to a control source directory or not
#
  my $self = shift;
  my $value = shift;

  return( $self->do('COMPARE', $value) );
}

sub do_continue {
#
#  Test if should continue a previous run of the script or not
#
  my $self = shift;
  my $value = shift;

  return( $self->do('CONTINUE', $value) );
}

sub completed {
#
#  Test if this part has already been completed or not
#
  my $self = shift;
  my $value = shift;

  #
  # If just begining mark completed as false so that 
  # everything will be done until it gets to the point of
  # running tests.
  #
  if ( $self->do('BEGIN') ) {
    return( 0 );
  # Now return completed based on the value of COMPLETED 
  } else {
    return( $self->do('COMPLETED', $value) );
  }
}

sub numerically {
#
# Use with "sort" to sort values in numerical order
#
  $a <=> $b;
}


sub exec {
#
# Extend the exec method so that if continuation option is set and this
# section already done, it doesn't actually execute.
#
  my $self = shift;
  my $command = shift;
  my $die_msg = shift;
  my $echo = shift;

  if ( $self->do_continue && $self->completed ) {
    print "Continue option set, skip running: $command\n";
    return;
  }
  $self->SUPER::exec( $command, $die_msg, $echo );
}

sub run {
#
# Extend the run method so that if continuation option is set and this
# section already done, it doesn't actually do anything.
#
  my $self = shift;
  my $desc = shift;
  my $logfile = shift;

  print("DEBUG: logfile is $logfile\n");
 
  if ( $self->do_continue && $self->completed ) {
    print "Continue option set, skip running model: $desc\n";
    if ( defined($logfile) ) { $self->{'LOGFILE'} = $logfile; }
    return;
  }
  $self->SUPER::run( $desc, $logfile );
}

sub compare_files {
#
# Compare the history files using "cprnc"
#
  my $self = shift;
  my $file1 = shift;
  my $file2 = shift;
  my $log = shift;
  my $die_msg = shift;

  print "Compare the two history tapes $file1 and $file2 output to $log:\n";
  if ( $self->do_continue && $self->completed ) {
    print "Continue option set, skip this file comparision\n";
    return;
  }
  $self->exec( "$CASE_DIR/cprnc $file1 $file2 > $log", "Error:: running cprnc" );
  open( LOG, "<$log" ) || $self->die( "Error, opening cprout file: $log" );
  my $different = 0;
  while( defined($_ = <LOG>) ) {
    if ( /RMS ([a-zA-Z0-9_-]+) (.+)/ ) {
      $CLM_test::cprout{$log}{$1} = $2;
      if ( $2 !~ /0.[0]+[eE]\+[0]+/ ) { 
        $different = 1; 
      }
    }
  }
  close( LOG );
  if ( defined($die_msg) ) {
    if ( $different ) {
      $self->die( "Error: $die_msg, see $log" );
    } else {
      print "RMS difference between all fields on the history files are identical!" .
            " Comparision successful, continuing...\n";
    }
  }
  return( $different );
}

sub report_compare {
#
# Report a summary of the comparisions
#
  my $self = shift;
  my @cprlogs = @_;

  my $log; my $var; my %diff;
  #
  # Loop through comparision logs and report which fields bit-for-bit or different
  #
  foreach $log ( @cprlogs ) {
    my $ref = $CLM_test::cprout{$log};
    if ( defined($ref) ) {
      print "For $log:\n";
      my %var_diffs = %$ref;
      my @var_list = sort( keys(%var_diffs) );
      my @bit4bit_list;
      my @diff_list;
      my $max_diff = 0.0;
      foreach $var ( sort( keys(%var_diffs) ) ) {
        if ( $var_diffs{$var} > $max_diff ) { $max_diff = $var_diffs{$var}; }
        if ( $var_diffs{$var} =~ /0.[0]+[eE]\+[0]+/ ) {
          push @bit4bit_list, $var;
        } else {
          push @diff_list, $var;
        }
      }
      print "Fields on files: @var_list\n";
      $diff{$log} = "NOT";
      if ( $#bit4bit_list == $#var_list ) { 
        print "All fields are bit for bit\n"; 
        $diff{$log} = "b4b";
      } elsif ( $#diff_list == $#var_list ) { 
        print "No fields are bit for bit\n"; 
        print "Max difference: $max_diff\n";
      } else {
        print $#bit4bit_list+1 . " fields are bit for bit\n";
        print $#diff_list+1 . " fields are different: @diff_list\n";
        print "Max difference: $max_diff\n";
      }
    }
  }
  #
  # Now Loop through logs and report if bit-for-bit with control library
  #
  my $b4b = "yes";
  foreach $log ( @cprlogs ) {
    my $ref = $CLM_test::cprout{$log};
    if ( defined($ref) && ($log =~ /con.+control.+/) ) {
      if ( $diff{$log} eq "NOT" ) { $b4b = "NO"; }
    }
  }
  if ( $b4b eq "yes" ) {
    print "Code is bit-for-bit with control library\n";
  } else {
    print "\n\nWARNING!!!:: Code is NOT bit-for-bit with control library\n";
    print "             Verify that this is ok or validate that at least within roundoff.\n\n";
  }
}

sub lower_PEs {
#
# Change the number of PE's for testing
#
  my $self = shift;

  print "Lower the number of PE\'s used to ensure that different configurations match results:\n";
  if ( defined($SHMEM_CPUS) ) {
    $CLM_run::SAVE_SHMEM_CPUS = "$SHMEM_CPUS";   # Save number it was set to
    print "Current number of shared memory CPUS: " . $SHMEM_CPUS. "\n";
    $SHMEM_CPUS = $SHMEM_CPUS / 2;
    if ( $SHMEM_CPUS < 1 ) { 
      $SHMEM_CPUS = 1; 
      print "Leave number at:                      " . $SHMEM_CPUS. "\n";
    } else {
      print "Change number to:                     " . $SHMEM_CPUS. "\n";
    }
  } else {
    $CLM_run::SAVE_SHMEM_CPUS = 0;
    $SHMEM_CPUS = 1;
    print "Change number of shared memory CPUS:  " . $SHMEM_CPUS. "\n";
  }
  if ( ($SPMD eq "TRUE") && defined($SPMD_NODES) ) {
    $CLM_run::SAVE_SPMD_NODES = "$SPMD_NODES";   # Save number it was set to
    print "Current number of SPMD nodes: " . $SPMD_NODES . "\n";
    $SPMD_NODES = $SPMD_NODES / 2;
    if ( $SPMD_NODES < 2 ) {
      $SPMD_NODES = 2;
      if ( ! defined($SHMEM_CPUS) ) {
        $CLM_run::SAVE_SHMEM_CPUS = 0;
      }
      $SHMEM_CPUS = 1;
      print "Keep SPMD nodes at $SPMD_NODES change shared-memory CPU\'s to $SHMEM_CPUS\n";
    } else {
      print "Change number to            : " . $SPMD_NODES . "\n";
    }
    $SPMD_RUNCMND = undef; # un-define the SPMD run command so that node change will take effect
                           # this happens when machine_specs( "run" ) method is invoked
  }
}

sub start_log {
#
# Start the log-file
#
  my $self = shift;

  my $testlog = $self->{'TESTLOG'};
  `/bin/mv $testlog $testlog.last`;
  open( TESTLOG, ">$testlog" ) || die "Can not open $testlog";
  my $hostname = `/bin/hostname`; chomp( $hostname );
  my $date = `date`;
  my $OS = $self->{'OS'};
  print TESTLOG "Log of test results for $hostname a $OS on $date\n\n";
  close( TESTLOG );
}

sub read_mark {
#
# Read the restart filemark
#
  my $self = shift;

  my $filename = $SCRIPT_DIR . "/" . $self->{continue_file};
  open( MARK, "<$filename" ) || die "Could not open $filename";
  $_ = <MARK>;
  my $read_mark;
  if ( ! /([0-9]+_)/ ) {
    die "ERROR:: Continuation mark on $filename is not in correct format\n";
  }
  my $desc = $1;
  $read_mark = $self->mark( $desc );
  close( MARK );
  return( $read_mark );
}

sub write_mark {
#
# Write the restart file-mark out
#
  my $self = shift;
  my $desc = shift;

  my $filename = $SCRIPT_DIR . "/" . $self->{continue_file};
  my $mark = $self->mark( $desc );
  open( MARK, ">$filename" ) || die "Could not open $filename";
  print MARK "$mark\n";
  close( MARK );
}

sub mark {
#
# Return mark to designate how far the testing has proceeded
# description assumed to be a integer number followed by an optional description.
#
  my $self = shift;
  my $desc = shift;

  # Check that desc set
  if ( ! defined($desc) ) { die "Error: desc not sent to mark method"; }
  # Make sure begining number printed out as two characters
  my $number = $desc;
  my $rest;
  if ( $desc =~ /(^[0-9]+)(.*)$/ ) {
    $number = sprintf "%2.2d", $1;
    $rest   = $2;
  }
  return( "$number$rest" );
}

sub continuation_mark {
#
# Mark and/or check how far the script has completed.
#
  my $self = shift;
  my $desc = shift;

  if ( defined($self->{'ROOTDIR'}) ) {
    $ROOTDIR = $self->{'ROOTDIR'};
    $self->{'ROOTDIR'} = undef;
  }
  my $filename = $SCRIPT_DIR . "/" . $self->{continue_file};
  my $mark = $self->mark( $desc );
  $self->{'STATUS'} = "ran";
  #
  # If script is completed
  #
  if ( $desc eq "done" ) {
    `/bin/rm $filename`;
    my $testlog = $self->{'TESTLOG'};
    print "Results of the tests are in $testlog\n";
    open( TESTLOG, "<$testlog" ) || die "Can not open $testlog";
    while( defined($_ = <TESTLOG>) ) {
      print $_;
    }
    close( TESTLOG );
  #
  # If script is just starting
  #
  } elsif ($desc =~ /^0_/) {
    $self->do('BEGIN',"no");
    $self->start_log;
  #
  # If you are continuing and haven't found the part that isn't done yet
  #
  } elsif ( $self->do_continue && $self->completed ) {
    my $read = $self->read_mark;
    #
    # Set completed to no if reached the section where it stopped
    #
    if ( $mark =~ /^$read/ ) {
      $self->completed( "no" );
    }
    $self->{'STATUS'} = "skipped";
  } else {
    $self->write_mark( $desc );
    #
    # Print test results to log file
    #
    if ( $desc !~ /^0_/) {
      my $testlog = $self->{'TESTLOG'};
      open( TESTLOG, ">>$testlog" ) || die "Could not access testlog file: " . 
                                            $testlog;
      print TESTLOG  "$mark: " . "${RESOLUTION} " . $self->{'STATUS'} . "\n";
      close( TESTLOG );
    }
    #
    # If just doing part of the tests stop if reached the last test
    #
    my $lasttest = $self->{'LASTTEST'};
    if ( defined($lasttest) ) {
      my $last = $self->mark( $lasttest );
      $self->{'STATUS'} = "ran";
      # If at last test mark completed to yes and skip the rest of the tests
      if ( $mark =~ /^$last/ ) {
        $self->completed( "yes" );
      }
    }
  }
}

sub reset_PEs {
#
# Reset the number of PE's changed by "lower_PEs" to the original values
#
  my $self = shift;

  print "Reset the PE\'s\n";
  if ( defined($SPMD_RUNCMND) ) {
    $SPMD_RUNCMND  = undef;
  }
  if ( defined($SHMEM_CPUS) ) {
    $SHMEM_CPUS = $CLM_run::SAVE_SHMEM_CPUS;
    print "Shared memory CPU\'s: $SHMEM_CPUS\n";
  }
  if ( defined($SPMD_NODES) ) {
    $SPMD_NODES = $CLM_run::SAVE_SPMD_NODES;
    print "SPMD nodes: $SPMD_NODES\n";
  }

}

sub die {
#
# Terminate gracifully if found an error
#
  my $self = shift;
  my $desc = shift;

  if ( $self->do_nofail ) {
    my $testlog = $self->{'TESTLOG'};
    open( TESTLOG, ">>$testlog" ) || die "Could not access testlog file: " . 
                                          "$testlog";
    print TESTLOG "TEST FAILED -- $desc\n";
    close( TESTLOG );
    $self->{'STATUS'} = "FAIL";
  } else {
    die "$desc";
  }
}

1   # To make use or require happy






