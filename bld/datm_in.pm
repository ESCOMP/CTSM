#-----------------------------------------------------------------------------------------------
#
#	datm_in.pm			Erik Kluzek
#
#	Perl module to create a namelist for datm7 -- datm_in.
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

package datm_in;

use compnl;
@datm_in::ISA = qw(compnl  namelist);
%datm_in::NL = {};

#============================================================================

sub new {
#
# Constructor
#
  my $class = shift;
  my $optsref = shift;
  my $config = shift;

  my $interactive = $$optsref{'interactive'};
  my $file = $$optsref{'dir'} . "/datm_in";
  my $printlev = $$optsref{'printlev'};
  my $cfgdir      = $$optsref{'cfgdir'};

  my $self = $class->SUPER::new( "datm_in", \%datm_in::NL, $interactive, $file,
                                 "$cfgdir/DefaultDATM_IN_Namelist.xml", $config, 
                                 $printlev );

  $self->{'printlev'} = $printlev;
  $self->{'optsref'}  = $optsref;

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub set_output_values {

# Set the datm_in namelist variables.

  my ($self, $settings_ref) = @_;

  my $class = ref($self);
  my $nm = "$class\:\:set_output_values";

  my $NLref = $self->{'NLREF'};
  my $optsref = $self->{'optsref'};
  my $opt;

  # Read in XML defaults file
  $self->Read_Default_XML_Values( $settings_ref );

  unless (defined($NLref->{'tn460_factorFn'})) {
      $NLref->{'tn460_factorFn'} = $self->default_vals( 'tn460_factorfn' );
  }

}

#============================================================================

sub print_hash {
    my %h = @_;
    my ($k, $v);
    while ( ($k,$v) = each %h ) { print "$k => $v\n"; }
}


1   # to make use or require happy
