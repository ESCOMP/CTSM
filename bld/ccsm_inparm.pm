#-----------------------------------------------------------------------------------------------
#
#	ccsm_inparm.pm			Erik Kluzek
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

package ccsm_inparm;

use compnl;
@ccsm_inparm::ISA = qw(compnl  namelist);
%ccsm_inparm::NL;

#============================================================================

sub new {
#
# Constructor
#
  my $class = shift;
  my $optsref = shift;
  my $config = shift;

  my $interactive = $$optsref{'interactive'};
  my $file        = $$optsref{'dir'} . "/drv_in";
  my $printlev    = $$optsref{'printlev'};
  my $cfgdir      = $$optsref{'cfgdir'};

  my $self = $class->SUPER::new( "seq_infodata_inparm", \%ccsm_inparm::NL, $interactive, $file,
                                 "$cfgdir/DefaultCCSM_INPARM_Namelist.xml", $config, 
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

  my $name = $self->{'NAME'};
  my $addquotes = 1;

  $self->Read_Default_XML_Values( $settings_ref );

  # Start type
  if (defined($optsref->{'start_type'})) {
      $opt = $optsref->{'start_type'};
  } elsif (defined($NLref->{'start_type'})) {
      $opt = $NLref->{'start_type'};
  } else {
      $opt = $self->default_vals( 'start_type' );
  }
  if ( $opt eq "arb_ic" ) {
     $opt = "startup";
  }
  $NLref->{'start_type'} = namelist::quote_string($opt);
  $self->requiredVar( "start_type", $addquotes );

  # Case name
  if (defined($optsref->{'case'})) {
      $opt = $optsref->{'case'};
  } elsif (defined($NLref->{'case_name'})) {
      $opt = $NLref->{'case_name'};
  } else {
      $opt = $self->default_vals( 'case_name' );
  }
  my $case = $opt;
  if ( $case !~ /([^\/]{1,79})([^\/]*)/ ) {
    die "$nm Bad casename: $case\n";
  }
  # If casename too long
  if ( $2 ne "" ) {
    $case = $1;
    print "\n\nWARNING:: Truncating case_name in namelist from $opt to $case\n\n" if ($self->{'printlev'});
  }
  $NLref->{'case_name'} = namelist::quote_string($case);
  $self->requiredVar( "case_name", $addquotes );

  if ( ($NLref->{'atm_adiabatic'}   =~ /\.[tT][Rr]*[Uu]*[Ee]*\./) &&
       (($NLref->{'atm_ideal_phys'} =~ /\.[tT][Rr]*[Uu]*[Ee]*\./) ||
        ($NLref->{'aqua_planet'}    =~ /\.[tT][Rr]*[Uu]*[Ee]*\./) )  ) {
    die "$nm only one of atm_adiabatic, atm_ideal_phys and aqua_planet can be defined.\n";
  }
  # Check that "restart_file" is set if this is a branch simulation
  if ( $NLref->{'start_type'} eq "\'branch\'" ) {
     $self->requiredVar( "restart_file", $addquotes );
  }
  # Output path root
  if (defined($NLref->{'outpathroot'})) {
      $opt = $NLref->{'outpathroot'};
  } else {
      $opt = $optsref->{'dir'} . "/";
  }
  $NLref->{'outpathroot'} = namelist::quote_string($opt);
  if ( $opt !~ /^['"]+\// ) {
     $self->requiredVar( "outpathroot", $addquotes );
  }

  # Orbit (if not coupled)
  if ( $self->{'MODE'} ne "ccsm_seq" ) {
     unless ( defined($NLref->{'orb_obliq'}) and defined($NLref->{'orb_eccen'}) and
              defined($NLref->{'orb_mvelp'}) ) {
         unless ( defined($NLref->{'orb_iyear_ad'}) ) {
             $NLref->{'orb_iyear_ad'} = $self->default_vals( 'orb_iyear_ad' );
         }
     }
  }

}

#============================================================================

sub print_hash {
    my %h = @_;
    my ($k, $v);
    while ( ($k,$v) = each %h ) { print "$k => $v\n"; }
}


1   # to make use or require happy
