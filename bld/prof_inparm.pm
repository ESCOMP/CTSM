#-----------------------------------------------------------------------------------------------
#
#	prof_inparm.pm
#
#	Perl module to create prof_inparm namelist
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

package prof_inparm;

use compnl;
@prof_inparm::ISA = qw(compnl  namelist);
%prof_inparm::NL;

#============================================================================

sub new {
#
# Constructor
#
  my $class = shift;
  my $optsref = shift;
  my $config = shift;

  my $interactive = $$optsref{'interactive'};
  my %bld = &queryDefaultXML::read_cfg_file( $$optsref{'config'} );
  my $outfile;
  if ( $bld{'MODE'} eq "ccsm_seq" ) {
     $outfile = "drv_in";
  } else  {
     $outfile = "lnd.stdin";
  }

  my $file = $$optsref{'dir'} . "/$outfile";
  my $printlev = $$optsref{'printlev'};
  my $cfgdir      = $$optsref{'cfgdir'};

  my $self = $class->SUPER::new( "prof_inparm", \%prof_inparm::NL, $interactive, $file,
                                 "$cfgdir/DefaultPROF_INPARM_Namelist.xml", $config, 
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

  # Read in XML defaults file
  $self->Read_Default_XML_Values( $settings_ref );

  my $NLref = $self->{'NLREF'};
  my $optsref = $self->{'optsref'};
  my $opt;

  # If timing is enabled or not
  unless (defined($NLref->{'profile_disable'})) {
      $NLref->{'profile_disable'} = $self->default_vals( 'profile_disable' );
  }

  # If barriers are enabled before each timing call
  unless (defined($NLref->{'profile_barrier'})) {
      $NLref->{'profile_barrier'} = $self->default_vals( 'profile_barrier' );
  }

  # If timing output goes to single file or multitple files per MPI task
  unless (defined($NLref->{'profile_single_file'})) {
      $NLref->{'profile_single_file'} = $self->default_vals( 'profile_single_file' );
  }


}

#============================================================================

sub print_hash {
    my %h = @_;
    my ($k, $v);
    while ( ($k,$v) = each %h ) { print "$k => $v\n"; }
}


1   # to make use or require happy
