#============================================================================
#
#	CLM.pm   Mariana Vertenstein (original version by Erik Kluzed)
#
#	Perl5 Base Object module to handle the building and running
#	of the Community Land Model (CLM).
#
#	Other objects are used to extend the very basic functionality
#	of this object. The most important features will be carried
#	in the objects that extend this one.
#
#	Description of methods:
#
#	Important basic methods:
#
#	new ---------------- Constructor of CLM_run object.
#
#	Basic methods to build/run the model. (callable from a script)
#
#	chdir -------------- Change the directory.
#	checkdir ----------- Check that directory is as expected.
#	build_support ------ Build a support program  (e.g. makdep, the dependency generator).
#	build_namelist ----- Build the model namelist.
#
#	Hidden and utility methods
#
#	OS ----------------- Return the Operating system being used.
#	Namelist ----------- Return the namelist object.
#	OS_list ------------ Return the list of OS's build system ported to.
#	duplicate ---------- Create a duplicate object with same data as current object.
#	append_config_env -- Append build or run-time env varible settings to 
#                            the config_env.csh script.
#	print -------------- Print out contents of object.
#	clean -------------- Delete object and dependency files in BLDDIR and output files 
#    			     in EXEDIR.
#	arch --------------- Return the architecture type (same as preprocessor token used)
#	Platform ----------- Return a 3-letter description of Operating System 
#			     (O/S) running on.
#	config_env --------- Apply the env settings from the "config_env.csh" script 
#			     (created by configure.csh).
#	exec --------------- Execute the given system command.
#	unset_env ---------- Unset all the Env variables needed for a CLM model run.
#	do ----------------- Test if a given option is set (returns 0 or 1 (false or true).
#	do_clean ----------- Return true if option to clean old directories is set.
#	do_build_namelist -- Return true if option to build a namelist is set.
#	die ---------------- Do a clean termination of the script.
#                            (The primary purpose of this method is to allow
#                            it to be over-ridden or extended by other objects)
#
#============================================================================

use 5.004;   # Use at least this version of perl
use strict;
#use diagnostics;
use Cwd;

package CLM;
#
# Use the shortened syntax for the following ENV variables
#
use Env qw(CASE EXENAME GNUMAKE LAB MODEL_BLDDIR MODEL_EXEDIR RESOLUTION SCRIPT_DIR);
use CLM_namelist;
use CLM_lab;

#============================================================================

sub new {
#
# Constructor
#
  my $class = shift;
  my $self = {};

  $self->{'OS'} = $^O;                             # Operating system running on
  $self->{'RUNTYPE'} = ["initial", "restart", "branch"]; # Valid run types
  $self->{'NAMELIST'} = undef;                     # Namelist object
  $self->{'BUILD'} = "yes";                        # If should configure and build or not
  $self->{'CLEAN'} = "no";                         # If should clean old dirs
  $self->{'BUILD_NAMELIST'} = "yes";               # Build model namelist or not
  $self->{'LOG'} = "yes";                          # Create logfile or not
  $self->{'LOGFILE'} = undef;                      # Name of logfile
  if ( ! defined($SCRIPT_DIR) ) { $SCRIPT_DIR = "."; }   
  # List of Env variables, to unset when running with prompts
  # Don't include variables that may need to be configured differently at another site.
  # Only include variables that define the model configuration.
  $self->{'ENV_LIST'} = ["RESOLUTION", "MODEL_EXEDIR", "MODEL_BLDDIR", 
			 "CASE", "RUNTYPE", "NAMELIST", "DEBUG" ];
  bless( $self, $class );
  $self->{'NAMELIST'} = CLM_namelist->new( $self );
  # File with ENV vars from configure.csh
  $self->{'CONFIG_ENV_FILE'} = "$SCRIPT_DIR/config_env_". $self->arch . ".csh";
  #
  # Environment variables to set
  #
  if ( ! defined($LAB) ) { $LAB = "ncar"; }        # Lab
  return( $self );
}

#============================================================================

sub arch {
#
# Return appropropriate architecture name
#
  my $self = shift;

  $_ = `uname -s`;
  my $ARCH;
  if ( /UNICOS/ ) {
     $ARCH = "UNICOSmp";
  } else {
     $ARCH = $_;
  }
  chomp( $ARCH );
  return( $ARCH );
}

#============================================================================

sub OS {
#
# Return the OS
#
  my $self = shift;

  my $OS = $self->{'OS'};
  return( $OS );
}

#============================================================================

sub Namelist {
#
# Set the namelist object
#
  my $self = shift;

  my $NL = $self->{'NAMELIST'};
  return( $NL );
}

#============================================================================

sub OS_list {
#
# Return the list of OS's build system is ported to.
#
  my $self = shift;

  my @list = ( "dec_osf", "linux", "solaris", "irix", "aix", "unicosmp" );
  return( @list );
}

#============================================================================

sub duplicate {
#
# Create a duplicate of the current object with the same data
#
  my $self = shift;

  my $new = ref($self)->new;
  foreach my $key ( keys(%$self) ) {
    $new->{$key} = $self->{$key};
  }
  return( $new );
}

#============================================================================

sub print {
#
# Print out contents of object
#
  my $self = shift;

  my $OS = $self->{'OS'};
  print "Lab running at: " . $LAB . " O/S: $OS\n";
  if ( $self->do_build ) { print "DO Build the model\n"; }
  else                   { print "Do NOT Build the model\n"; }
  if ( $self->do_build_namelist ) { print "DO Build the model namelist\n"; }
  else                            { print "Do NOT Build the model namelist\n"; }
  print "Logfile: " . $self->{'LOGFILE'} . "\n";
}

#============================================================================

sub clean {
#
# Delete the object and dependency files in the MODEL_BLDDIR and the output files
# in the MODEL_EXEDIR.
#
  my $self = shift;
  my $input = shift;

  if ( defined($input) ) {
    print "Clean out the directory: $input\n";
    $self->chdir( $input );
    if ( -f "Makefile" && defined($GNUMAKE) ) { 
      $self->exec( "$GNUMAKE clean" );
    } else {
      $self->exec( "/bin/rm -f *.[do]" );
    }
  } else {
    if ( -d $MODEL_BLDDIR ) {
      print "Clean out the build directory: $MODEL_BLDDIR\n";
      $self->chdir( $MODEL_BLDDIR );
      if ( -f "Makefile" && defined($GNUMAKE) ) { 
        $self->exec( "$GNUMAKE realclean" );
      } else {
        $self->exec( "/bin/rm -f *.[do] *.mod" );
      }
    }
    if ( -d $MODEL_EXEDIR ) {
      print "Clean out the executable directory: $MODEL_EXEDIR\n";
      $self->chdir( $MODEL_EXEDIR );
      $self->exec( "/bin/rm -f clm* timing* $EXENAME" );
    }
  }
}

#============================================================================

sub Platform {
#
# Return a three letter string describing the platform
#
  my $self = shift;

  $_ = $self->{'OS'};
  my $platform;
  if ( /dec_osf/ ) {
    $platform = "osf";
  } elsif ( /linux/ ) {
    $platform = "PC";
  } elsif ( /solaris/ ) {
    $platform = "sun";
  } elsif ( /irix/ ) {
    $platform = "sgi";
  } elsif ( /aix/ ) {
    $platform = "aix";
  } elsif ( /unicos/ ) {
    $platform = "unicosmp";
  } else {
    $platform = "xxx";
  }
  return( $platform );
}

#============================================================================

sub unset_env {
#
# Unset all the associated ENV variables
#
  my $self = shift;

  print "Unset the ENV variables for next pass:: \n";
  my $ref = $self->{'ENV_LIST'};
  my @list = @$ref;
  foreach my $i ( @list ) {
    print "Unset env variable: $i\n";
    if ( defined($ENV{$i}) ) { $ENV{$i} = undef; }
  }
}

#============================================================================

sub append_config_env {
#
# Append the build or run time ENV vars onto the end of the config_env.csh file.
#
  my $self = shift;
  my $type = shift;
  my @vars = @_;

  my $file = $self->{'CONFIG_ENV_FILE'};
  my $tmp = "$file.tmp";
  open( FILE, "<$file" ) || die "Could not open file: $file";
  open( OUT, ">$tmp" ) || die "Could not open file: $tmp";
  while( defined($_ = <FILE>) ) {
    print OUT $_;
    if ( /^#CONFIGURE-END/ ) {
      last;
    }
  }
  #
  # If appending run-time information skip past build info if exists
  #
  if ( $type eq "run" ) {
    $_ = <FILE>;
    if ( defined($_) ) {
      if ( /^#build/ ) {
        while( defined($_ = <FILE>) ) {
          print OUT $_;
          if ( /^#build-END/ ) {
            last;
          }
        }
      }
    }
  }
  # Make sure stacksize is set to unlimited
  print OUT "limit stacksize unlimited\n";
  close( FILE );
  print OUT "#$type\n";
  print OUT "# Add $type-time env vars: from $0\n#  " . `date`;
  foreach my $var ( @vars ) {
    if ( ! defined($ENV{$var}) )  { 
      print "Warning: The ENV variable: $var not set in the main script\n"; 
    } else {
      # If variable has quotes around it
      if ( $var =~ /\"(.+)\"/ ) {
        print OUT "setenv $var $ENV{$var}\n";
      } else {
        print OUT "setenv $var \"$ENV{$var}\"\n";
      }
    }
  }
  print OUT "#$type-END\n";
  close( OUT );
  `/bin/mv $tmp $file`;
}

#============================================================================

sub do {
#
# Test if given option is turned on or not
# (Or set if parameter given)
#
  my $self = shift;
  my $option = shift;
  my $value = shift;

  if ( defined($value) ) {
    $self->{$option} = $value;
  } elsif ( $self->{$option} =~ /[Yy][eE]*[Ss]*/ ) {
    return( 1 );
  } else {
    return( 0 );
  }
}

#============================================================================

sub do_build {
#
#  Test if should build/make model or not
#
  my $self = shift;
  my $value = shift;

  return( $self->do('BUILD', $value) );
}

#============================================================================

sub do_clean {
#
#  Test if should clean out the old directories or not.
#
  my $self = shift;
  my $value = shift;

  return( $self->do('CLEAN', $value) );
}

#============================================================================

sub do_build_namelist {
#
#  Test if should build model-namelist or not
#
  my $self = shift;
  my $value = shift;

  return( $self->do('BUILD_NAMELIST', $value) );
}

#============================================================================

sub do_log {
#
#  Test if should create log file for run or send output to screen
#
  my $self = shift;
  my $value = shift;

  return( $self->do('LOG', $value) );
}

#============================================================================

sub chdir {
#
# Change to given directory, return error if can't
#
  my $self = shift;
  my $dir = shift;

  print "Change to $dir\n";
  if ( ! defined( $dir ) ) {
    die "Can not change to an empty directory";
  }
  if ( ! chdir( $dir ) ) {
    die "Could not change to directory: $dir\n";
  }
}

#============================================================================

sub checkdir {
#
# Check that current directory is same as expected directory
#
  my $self = shift;
  my $dir = shift;
  my $die = shift;

  my $pwd = cwd( );
  if ( ! chdir( $dir ) ) {
    die "checkdir: Could not change to directory: $dir\n";
  }
  my $newpwd = cwd( );
  if ( $newpwd ne $pwd ) {
    die "$die\n\n";
  }
  if ( ! chdir( $pwd ) ) {
    die "checkdir: Could not change back to directory: $pwd\n";
  }
}

#============================================================================

sub exec {
#
# Method to execute given command line to Bourne shell
#
  my $self = shift;
  my $command = shift;
  my $die_msg = shift;
  my $echo = shift;

  if ( ! defined($command)  ) {
    die "Command not given to exec method";
  }
  print "$command\n";
  # Use system, so that results will be echoed out, if echo option set
  if ( defined($echo) && ($echo == "echo") ) {
    system( $command );
  # Otherwise, use backtics so that results of command will not be seen
  } else {
    `$command`;
  }
  #
  # Die on an error
  #
  if ( $? != 0 ) {
    if ( defined($die_msg) ) { print "$die_msg\n"; }
    $self->die( "ERROR:: Trouble executing: $command" );
    return;
  }
}

#============================================================================

sub build_support {
#
# Build a support program needed for the process
#
  my $self = shift;
  my $flag = shift;

  print "build_support:: Build a support program in directory:" . cwd() . "\n";
  if ( $flag eq "clean" ) {
    $self->exec( "$GNUMAKE clean" );
  }
  $self->exec( $GNUMAKE );
}

#============================================================================

sub build_namelist {
#
# Create the namelist based on the OS ENV variables and the values set
# in the CLMEXP associative arrays
#
  my $self = shift;

  if ( $self->do_build_namelist ) {
    print "build_namelist:: Build the model namelist\n";
    $self->{'NAMELIST'}->build;
  }
}

#============================================================================

sub config_env {
#
# Apply the settings in the config_env.csh script created by configure.csh
#
  my $self = shift;

  print "Apply the settings from the last configure: \n";
  my $file = $self->{'CONFIG_ENV_FILE'};
  open( FILE, "<$file" ) || die "Could not open file: $file";
  while( defined($_ = <FILE>) ) {
    if ( /setenv\s+(\S+)\s+\"*([^\s"]+)\"*/ ) {
      print "Set env $1 = $2 \n";
      $ENV{$1} = $2;
    }
  }
  close( FILE );
}

#============================================================================

sub die {
#
# Terminate gracifully if found an error
#
  my $self = shift;
  my $desc = shift;

  die "$desc";
}

1   # To make use or require happy
