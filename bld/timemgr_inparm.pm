#-----------------------------------------------------------------------------------------------
#
#	timemgr_inparm.pm			Erik Kluzek
#
#	Perl module to create a namelist for sequential CCSM.
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

package timemgr_inparm;

use compnl;
@timemgr_inparm::ISA = qw(compnl  namelist);
%timemgr_inparm::NL;

#============================================================================

sub new {
#
# Constructor
#
  my $class = shift;
  my $optsref = shift;
  my $config = shift;

  my $interactive = $$optsref{'interactive'};
  my $file = $$optsref{'dir'} . "/drv_in";
  my $printlev = $$optsref{'printlev'};
  my $cfgdir      = $$optsref{'cfgdir'};

  my $self = $class->SUPER::new( "seq_timemgr_inparm", \%timemgr_inparm::NL, $interactive, $file,
                                 "$cfgdir/DefaultTIMEMGR_INPARM_Namelist.xml", $config, 
                                 $printlev );

  $self->{'printlev'} = $printlev;
  $self->{'optsref'}  = $optsref;

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub set_output_values {

# Set the CCSM namelist variables.

  my ($self, $settings_ref) = @_;

  my $class = ref($self);
  my $nm = "$class\:\:set_output_values";

  my $NLref = $self->{'NLREF'};
  my $optsref = $self->{'optsref'};
  my $opt;

  # Read in XML defaults file
  $self->Read_Default_XML_Values( $settings_ref );
  # Length of simulation
  if ( ! defined($NLref->{'stop_option'}) ) {
    $opt = $self->default_vals( 'stop_option' );
  } else {
    $opt = $NLref->{'stop_option'};
  }
  $NLref->{'stop_option'} = namelist::quote_string($opt);
  $self->requiredVar( "stop_option" );
  unless ( defined($NLref->{'stop_n'}) ) {
    $NLref->{'stop_n'} = $self->default_vals(  'stop_n' );
  }
  $self->requiredVar( "stop_n" );
 
  # Restart frequency
  unless ( defined($NLref->{'restart_option'}) ) {
    $opt = $self->default_vals( 'restart_option' );
  } else {
    $opt = $NLref->{'restart_option'};
  }
  $NLref->{'restart_option'} = namelist::quote_string($opt);
  unless ( defined($NLref->{'restart_n'}) ) {
    $NLref->{'restart_n'} = $self->default_vals( 'restart_n' );
  }

  # Starting date
  unless ( defined($NLref->{'start_ymd'}) ) {
    $NLref->{'start_ymd'} = $self->default_vals( 'start_ymd' );
  }
  $self->requiredVar( "start_ymd" );

  # Starting time of day
  unless ( defined($NLref->{'start_tod'}) ) {
    $opt = $self->default_vals( 'start_tod' );
    if ( defined($opt) ) {
       $NLref->{'start_tod'} = $opt;
    }
  }

  # Coupling frequency
  unless ( defined($NLref->{'atm_cpl_dt'}) ) {
      $NLref->{'atm_cpl_dt'} = $self->default_vals( 'atm_cpl_dt' );
  }
  $self->requiredVar( "atm_cpl_dt" );

  # end restart
  unless ( defined($NLref->{'end_restart'}) ) {
      $NLref->{'end_restart'} = $self->default_vals( 'end_restart' );
  }

}

#============================================================================

sub print_hash {
    my %h = @_;
    my ($k, $v);
    while ( ($k,$v) = each %h ) { print "$k => $v\n"; }
}


1   # to make use or require happy
