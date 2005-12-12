#
#	CLM_namelist.pm			Erik Kluzek (modified by Mariana Vertenstein)
#
#	Perl module to create a namelist for the offline CLM2 model.
#
#------------------------------------------------------------------------
#
# Note: Implied inputs are: %main::CLMEXP 
#
#	This object interacts with the %main::CLMEXP associative arrays 
#       set in the main program. The main program sets these associative 
#       arrays as input to this object. This object then converts the keys 
#       to lowercase, and adds the needed default values that aren't set 
#       from the values set in the main program.
#
#------------------------------------------------------------------------
#
#	Description of methods:
#
#	new --------------------------- Constructor
#	build ------------------------- Build namelist.
#	checkinputfile ---------------- Check if the given input file exists.
#	convert_case ------------------ Convert the case of the associative array keys to lowercase.
#	set_default_values ------------ Set defaults for all configurations.
#	write_screen -------------------Print namelist to screen.
#	write_file ---------------------Write the namelist out to a file
#
#	$Id$
#
use strict;
#use diagnostics;
use Cwd;

package CLM_namelist;
#
# Use the shortened expression to address the following ENV variables
#
use Env qw(CASE LAB OS LM_DATDIR MODEL_DATDIR NAMELIST RESOLUTION RUNTYPE MODEL_EXEDIR SCRIPT_DIR);
#
# This script takes the %main::CLMEXP associative array and stores the keys in
# Lower case to the following lowercase copies. It uses the values passed in
# and sets needed default values based on ENV variables. Then it writes out
# a namelist according to the corresponding resultant associative array.
#

sub new {
#
# Constructor
#
  my $class = shift;
  my $run = shift;
  my $self  = {};
  $self->{'FILENAME'} = undef;                    # Filename of output namelist
  $self->{'CLM_RUN'} = $run;                      # CLM_run object
  bless( $self, $class );
  $self->setdir;
  return( $self );
}

#============================================================================

sub setdir {
#
# Set the default directories if they aren't already set
#
  my $self = shift;

  if ( ! defined($MODEL_DATDIR) ) {
    $MODEL_DATDIR = $CLM_run::MODEL_DATDIR_default->value($LAB, $OS);
  }
  if ( ! defined $LM_DATDIR ) {
    $LM_DATDIR    = "$MODEL_DATDIR/lnd/clm2";     # clm2 land-model data directory
  }
}

#============================================================================

sub build {
#
# Build the namelist
#
  my $self = shift;

  # Reset namelist when rebuilding
  %CLM_namelist::CLMEXP = {};
  $self->{'FILENAME'} = $NAMELIST;
  $self->convert_case;                     # Convert the main associative arrays to lower case
  $self->set_default_values();             # Set the default values that aren't already set.
  $self->write_file;                       # Write the namelist file out
  my $pwd = Cwd::cwd( );                   # Current directory
  if ( $pwd ne $SCRIPT_DIR ) {
    $self->write_file( $SCRIPT_DIR );      # Write the namelist to the location of the main script
  }
}

#============================================================================

sub convert_case {
#
# Convert the main associative arrays to the local version and change the case
# of the keys to lowercase. Also terminate if there are two keys with the same
# name but different case.
#
  my $self = shift;

  my $key;
  foreach $key ( keys(%main::CLMEXP) ) {
    my $lckey = $key;
    $lckey =~ tr/[A-Z]/[a-z]/;
    if ( defined($CLM_namelist::CLMEXP{$lckey}) ) {
      print "$lckey already defined\n";
      die "Fix your CLMEXP namelist so that two definitions of $lckey do not exist";
    }
    $CLM_namelist::CLMEXP{$lckey} = $main::CLMEXP{$key};
  }
}

#============================================================================

sub write_screen {
#
# Print the namelist out
#
  my $self = shift;
  my $ref = shift;

  my $key;
  my %namelist = %$ref;
  foreach $key ( sort( keys(%namelist) ) ) {
    if ( defined($namelist{$key}) ) {
      print " $key = $namelist{$key}\n";
    }
  }
}

#============================================================================

sub checkinputfile {
#
# Check that input file exists
#
  my $self = shift;
  my $EXPNLref = shift; 
  my $item = shift;

  my %EXPNL = %$EXPNLref;
  my $name = $EXPNL{$item};
  if ( $name !~ /\'(.+)\'/ ) {
    die "$item needs \' around filename";
  }
  my $infile = $1;
  # If file does not exist on disk
  if ( ! -f $infile ) {
    # Check if may be a Mass Store type path
    if ( $infile =~ /^\/[^a-z\/]+/ ) {
      print "Warning: $infile does not exist on disk but may exist on mass store\n";
    } else {
      # Check if exists under exec directory
      $infile =~ /([^\/]+$)/; # Get filename without path info.
      my $file = $1;
      if ( ! defined($MODEL_EXEDIR) ) {
        print "Warning: $infile does not exist on disk but may be under exec directory\n";
      } elsif ( -f "$MODEL_EXEDIR/$file" ) {
        print "Warning: $infile does not exist on disk but is under exec directory\n";
      } else {
        die "$item ($file) does not exist on local disk and does not appear to be a valid mss pathname\n";
      }
    }
  }
}

#============================================================================

sub set_default_values {
#
# Set the default values that haven't already been set.
# The ENV variables: RUNTYPE, LM_DATDIR, and CASE will over-ride
# any parameters already set in the associative array.
#
  my $self = shift;
  my $res = $RESOLUTION;
  my $fsurdat = " "; 
  my $finidat = " ";
  my $nrevsn  = " ";

  # Type of run 
  if (      $RUNTYPE =~ /initial/ ) {
    $CLM_namelist::CLMEXP{'nsrest'} = 0;
  } elsif ( $RUNTYPE =~ /restart/ ) {
    $CLM_namelist::CLMEXP{'nsrest'} = 1;
  } elsif ( $RUNTYPE =~ /branch/ ) {
    $CLM_namelist::CLMEXP{'nsrest'} = 3;
  } else {
    die "Not a valid run type: $RUNTYPE"
  }

  # Case name (truncate if casename is too long)
  my $case = $CASE;
  if ( $case !~ /([^\/]{1,16})([^\/]*)/ ) {
    die "Bad casename: $case\n";
  }
  if ( $2 ne "" ) {
    $case = $1;
    print "\n\nWARNING:: Truncating caseid in namelist from $CASE to $case\n\n";
  }
  $CLM_namelist::CLMEXP{'caseid'} = "\'$case\'";

  # Title
  if ( ! defined($CLM_namelist::CLMEXP{'ctitle'}) ) {
    $CLM_namelist::CLMEXP{'ctitle'} = "\'$CASE :: res: $RESOLUTION\'";
  }

  # Data directory existence check
  if (! defined $LM_DATDIR) {
    die "LM_DATDIR needs to be defined";
  }

  # Surface dataset
  SWITCH: {
      ($res =~ /T31cnall/) && do {
        $fsurdat = "$LM_DATDIR/srfdata/csm/surface-data.096x048_atm.gx3v5_ocn.050203.nc";
	$finidat = "$LM_DATDIR/inidata_2.1/offline/ccsm3_bgc19_F_2.clm2.i.1900-01-01-00000.c050905.nc";
	last SWITCH;
      };
      ($res =~ /T31cn/) && do {
        $fsurdat = "$LM_DATDIR/srfdata/csm/surface-data.096x048_atm.gx3v5_ocn.050203.nc";
	$finidat = "";
	last SWITCH;
      };
      ($res =~ /T31/) && do {
        $fsurdat = "$LM_DATDIR/srfdata/csm/surface-data.096x048_atm.gx3v5_ocn.050203.nc";
	$finidat = "";
	last SWITCH;
      };
  }
  $CLM_namelist::CLMEXP{'fsurdat'}        = "\'$fsurdat\'";
  $CLM_namelist::CLMEXP{'finidat'}        = "\'$finidat\'";
  $CLM_namelist::CLMEXP{'offline_atmdir'} = "\'$LM_DATDIR/NCEPDATA\'";
  $CLM_namelist::CLMEXP{'fpftcon'}        = "\'$LM_DATDIR/pftdata/pft-physiology-cn16.c040719\'";
  $CLM_namelist::CLMEXP{'fndepdat'}       = "\'$LM_DATDIR/ndepdata/1890/regrid_ndep_clm.nc\'";
  $CLM_namelist::CLMEXP{'frivinp_rtm'}    = "\'$LM_DATDIR/rtmdata/rdirc.05\'";
  $CLM_namelist::CLMEXP{'nrevsn'}         = "\'$nrevsn\'";
  $CLM_namelist::CLMEXP{'dtime'}          =  1800;
  $CLM_namelist::CLMEXP{'wrtdia'}         = ".true.";
}

#============================================================================

sub write_file {
#
# Write out the namelist based on values set in the associative arrays
#
  my $self = shift;
  my $dir  = shift;

  my $file = $self->{'FILENAME'};
  if ( defined($dir) ) {
    $file = "$dir/$file";
  }
  print "writing the namelist file $file \n";
  if ( -f $file ) { unlink( $file ); }
  my $key;
  open( OUT, ">$file" ) || die "Can not open namelist file: $file";
  print OUT "&clmexp\n";
  foreach my $key ( sort( keys(%CLM_namelist::CLMEXP) ) ) {
    if ( defined($CLM_namelist::CLMEXP{$key}) ) {
      print OUT " $key\t\t= $CLM_namelist::CLMEXP{$key}\n";
    }
  }
  print OUT "/\n";
  close( OUT );
}

1   # to make use or require happy
