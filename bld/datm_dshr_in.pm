#-----------------------------------------------------------------------------------------------
#
#	datm_dshr_in.pm			Erik Kluzek
#
#	Perl module to create a namelist for datm7 -- datm_dshr_in.
#
#	Description of methods:
#
#	new ----------------------  Constructor
#	set_output_values --------  Set output values based on precedence of the various input
#                                   values and ensure that a valid namelist is produced.
#
#-----------------------------------------------------------------------------------------------

use strict;
#use diagnostics;
use Cwd;

package datm_dshr_in;

use compnl;
@datm_dshr_in::ISA = qw(compnl  namelist);
%datm_dshr_in::NL = {};

#============================================================================

sub new {
#
# Constructor
#
  my $class = shift;
  my $optsref = shift;
  my $config = shift;

  my $interactive = $$optsref{'interactive'};
  my $file = $$optsref{'dir'} . "/datm_dshr_in";
  my $printlev = $$optsref{'printlev'};
  my $cfgdir      = $$optsref{'cfgdir'};

  my $self = $class->SUPER::new( "dshr_nml", \%datm_dshr_in::NL, $interactive, $file,
                                 "$cfgdir/DefaultDATM_DSHR_NML_Namelist.xml", $config, 
                                 $printlev );

  $self->{'printlev'} = $printlev;
  $self->{'optsref'}  = $optsref;

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub set_output_values {

# Set the datm_dshr_in namelist variables.

  my ($self, $settings_ref) = @_;

  my $class = ref($self);
  my $nm = "$class\:\:set_output_values";

  my $printlev = $self->{'printlev'};

  my $NLref = $self->{'NLREF'};
  my $optsref = $self->{'optsref'};
  my $opt;

  my $dir = $$optsref{'dir'};
  # Read in XML defaults file
  $self->Read_Default_XML_Values( $settings_ref );

  unless (defined($NLref->{'restpfile'})) {
      $NLref->{'restpfile'} = $self->default_vals( "restpfile" );
  }
  $self->requiredVar( "restpfile" );
  unless (defined($NLref->{'domainfile'})) {
      $NLref->{'domainfile'} = $self->default_vals( "domainfile" );
  }
  $self->requiredVar( "domainfile" );
  unless (defined($NLref->{'streams'})) {
      my $streams = $self->default_vals( "streamsxmlfile", 1 );
      my $init_year = $$optsref{'cycle_init_year'};
      my $beg_year  = $$optsref{'cycle_beg_year'};
      my $end_year  = $$optsref{'cycle_end_year'};
      if ( $beg_year > $end_year ) {
         print "\n\ncycle_beg_year=$beg_year cycle_end_year=$end_year\n";
         die "ERROR($nm):: cycle_beg_year greater than cycle_end_year\n"
      }
      $NLref->{'streams'} = "'$streams $init_year $beg_year $end_year'";
  }
  $self->requiredVar( "streams" );

  unless (defined($NLref->{'infodbug'})) {
      $NLref->{'infodbug'} = $self->default_vals( "infodbug" );
  }

  unless (defined($NLref->{'datamode'})) {
      $NLref->{'datamode'} = $self->default_vals( "datamode" );
  }

  # Check that "restSfile" and "restBfile" is set if this is a branch simulation

  if ($$optsref{'start_type'} eq 'branch' ) {
     $self->requiredVar( "restsfile" );
     $self->requiredVar( "restbfile" );
  }

}

#============================================================================

sub print_hash {
    my %h = @_;
    my ($k, $v);
    while ( ($k,$v) = each %h ) { print "$k => $v\n"; }
}


1   # to make use or require happy
