#
#	compnl.pm			Erik Kluzek
#
#	Perl module to deal with namelists for a component namelist.
#
#------------------------------------------------------------------------
#
#	Description of methods:
#
#       new ----------------------- Constructor.
#       Read_Default_XML_Values --- Read the XML file for the default values.
#       default_vals -------------- Return a particular default value.
#       do_interactive ------------ Return true if in interactive mode.
#	checkinputfile ------------ Check if the given input file exists.
#
#	$Id$
#
use strict;
#use diagnostics;
use Cwd;

package compnl;

use queryDefaultXML;

use namelist;
@compnl::ISA = "namelist";
#
# Extend the namelist class to have a method to check that input files
# exist either on local disk, full-path given, or on Mass Store.
#
sub new {
#
# Constructor
#
  my $class         = shift;
  my $name          = shift;
  my $refNL         = shift;
  my $interactive   = shift;
  my $file          = shift;
  my $defaults_file = shift;
  my $config        = shift;
  my $printlev      = shift;

  if ( ! defined($file) ) {
    die "ERROR($class): compnl constructor was not sent the namelist filename\n";
  }

  my $self = $class->SUPER::new( $name, "$file", $refNL, $printlev );
     
  $self->{'DEFAULTS_FILE'} = $defaults_file;       # XML File with default values
 
  if ( ($interactive != 0) && ($interactive != 1) ) {
    die "ERROR($class): interactive option passed in to new was not valid: $interactive\n";
  }
  $self->{'INTERACTIVE'}      = $interactive;    # Interactive mode (0 or 1)
  $self->{'config'}           = $config;
  $self->{'defaults_ref'}     = undef;
  $self->{'lnd_defaults_ref'} = undef;

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub Read_Default_XML_Values {
#
# Read in the default values from the XML file
#
  my ($self, $settings_ref) = @_;

  my %inputopts;
  my $opts_ref = $self->{'optsref'};
  my $res      = $$opts_ref{'RESOLUTION'};
  my $lnd_res  = $$opts_ref{'lnd_RESOLUTION'};
  my $file     = $self->{'DEFAULTS_FILE'};
  my $nm       = $$opts_ref{'ProgName'};
  $inputopts{'file'}     = $file;
  $inputopts{'namelist'} = $self->{'NAME'};
  $inputopts{'printing'} = 0;
  $inputopts{'ProgName'} = $$opts_ref{'ProgName'};
  $inputopts{'cmdline'}  = $$opts_ref{'cmdline'};
  $inputopts{'csmdata'}  = $$opts_ref{'csmdata'};
  $inputopts{'config'}   = $self->{'config'};
  $inputopts{'res'}      = $res;
  my $printing = $$opts_ref{'printlev'};
  print "($nm) Read: $file for resolution $res\n" if $printing;
  $self->{'defaults_ref'} = &queryDefaultXML::ReadDefaultXMLFile( \%inputopts, 
                                            $settings_ref );
  if ( $res ne $lnd_res ) {
     $inputopts{res} = $lnd_res;
     print "($nm) Read: $file for resolution $lnd_res\n" if $printing;
     $self->{'lnd_defaults_ref'} = &queryDefaultXML::ReadDefaultXMLFile( \%inputopts,
                                            $settings_ref );
  } else {
     $self->{'lnd_defaults_ref'} = undef;
  }

}

#============================================================================

sub default_vals {
#
# Return default value from XML file
#
  my $self     = shift;
  my $var      = shift;
  my $noquotes = shift;
  my $lnd_res  = shift;

  my $class = ref($self);
  my $nm    = "$class\:\:default_vals";

  my $use_land = 0;
  if ( defined($lnd_res) && defined($self->{'lnd_defaults_ref'}) ) {
     $use_land = $lnd_res;
  }
  my $quotes = 1;
  if ( defined($noquotes) ) {
     $quotes = ! $noquotes;
  }
  my $defaults_ref;
  if ( ! $use_land ) {
     $defaults_ref = $self->{'defaults_ref'};
  } else {
     $defaults_ref = $self->{'lnd_defaults_ref'};
  }
  if ( ! defined($defaults_ref) ) {
    die "$nm: XML file has NOT been read yet!\n";
  }
  my %defaults = %$defaults_ref;
  my $value    = $defaults{$var}{'value'};
  my $isadir   = $defaults{$var}{'isdir'};
  my $isafile  = $defaults{$var}{'isfile'};
  if ( $quotes &&  ($isafile || $isadir) ) {
     $value = "\'$value\'";
  }
  return( $value );
}

#============================================================================

sub requiredVar {
#
# Make sure given variable is set in the namelist
#
  my $self  = shift;
  my $var   = shift;
  my $quote = shift;

  my $opts_ref = $self->{'optsref'};
  my $name       = $self->{'NAME'};
  my $NLref      = $self->{'NLREF'};
  my $ProgName   = $$opts_ref{'ProgName'};
  my $RESOLUTION = $$opts_ref{'RESOLUTION'};
  my $mask       = $$opts_ref{'mask'};

  my $value = $NLref->{$var};
  if ( !defined($value) || ($value =~ /^['"][ ]+['"]$/) ) {
      if ( $self->do_interactive ) {
          print "Enter a value for the $var variable for the $name namelist\n";
          my $opt = <>; chomp $opt;
          if ( defined($quote) && $quote ) {
             $opt = namelist::quote_string($opt);
          }
          $NLref->{$var} = $opt;
       
      } else {
          die "ERROR ($ProgName): The variable $var MUST be set in namelist $name.\n" .
          "       Set the namelist variable \`$var\` in the $name namelist.\n".
          "       This can be done on the command-line using the -namelist\n".
          "       option or in an input namelist file that is specified\n".
          "       using the -infile option. If the default was being used\n".
          "       perhaps there is not a value for this configuration and resolution\n".
          "       (res=$RESOLUTION and mask=$mask).\n";
      } 
  }
}

#============================================================================

sub do_interactive {
#
# Return true if interactive option set
#
  my $self = shift;

  my $value = $self->{INTERACTIVE};
  return( $value );
}

#============================================================================

sub checkinputfile {

# Check that the namelist value for an initial or boundary datasets is
# properly quoted.  Then check that the file exists on local filesystem.
# If the file is not found by looking at the full filepath, check for it in
# the directory where the sequential CCSM executable was created.

  my $self = shift;
  my $item = shift;

  my $class = ref($self);
  my $nm = "$class\:\:checkinputfile";

  my $EXPNLref = $self->{'NLREF'};
  my %EXPNL = %$EXPNLref;
  my $name = $EXPNL{$item};

  # check for quoting
  if ( $name !~ /["'](.*)['"]/ ) {
    die "$nm: $item needs quotes around filename: value = $name";
  }
  my $infile = $1;

  # Check for empty filename
  if ( $infile =~ /^[ ]+$/ ) {
     return;
  }

  my $found_message = "Found $item dataset on local disk.";

  # check full pathname
  if ( -f $infile ) { 
      print "$found_message\n" if ($self->{'printlev'}>1);
      return;
  }

  # check for file in directory containing sequential CCSM executable
  $infile =~ /([^\/]+$)/;     # strip filename from the path
  my $file = $1;
  my $MODEL_EXEDIR = $self->{'MODEL_EXEDIR'};
  if ( defined($MODEL_EXEDIR) ) {
      if ( -f "$MODEL_EXEDIR/$file" ) { 
	  print "$found_message\n" if ($self->{'printlev'}>1);
	  return;
      }
  }

  print "Warning($nm): $item dataset $infile not found on local disk\n".
        "This dataset must be copied or linked to the run directory.\n";
}

#============================================================================

1   # to make use or require happy
