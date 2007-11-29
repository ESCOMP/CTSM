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
use Streams::Template;

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

  my $dir     = $$optsref{'dir'};
  my $cfgdir  = $$optsref{'cfgdir'};
  my $csmdata = $$optsref{'csmdata'};

  # Read in XML defaults file
  $self->Read_Default_XML_Values( $settings_ref );

  if (defined($NLref->{'restpfile'})) {
    die "ERROR($nm):: can NOT set restpfile on namelist -- as archiving script depends on filename\n";
  }
  $NLref->{'restpfile'} = $self->default_vals( "restpfile" );
  $self->requiredVar( "restpfile" );
  unless (defined($NLref->{'domainfile'})) {
      $NLref->{'domainfile'} = $self->default_vals( "domainfile" );
  }
  $self->requiredVar( "domainfile" );
  $self->checkinputfile('domainfile') if $$optsref{'test'};
  unless (defined($NLref->{'streams'})) {
      my $init_year = $$optsref{'cycle_init_year'};
      my $beg_year  = $$optsref{'cycle_beg_year'};
      my $end_year  = $$optsref{'cycle_end_year'};
      if ( $beg_year > $end_year ) {
         print "\n\ncycle_beg_year=$beg_year cycle_end_year=$end_year\n";
         die "ERROR($nm):: cycle_beg_year greater than cycle_end_year\n"
      }
      my $outstreams = "$dir/datm7.streams.txt";
      my $template  = $self->default_vals( "streamstemplate" );

      my %inputopts;
      if ( $$optsref{'silent'} ) {
         $inputopts{'printing'} = 0;
      } else {
         $inputopts{'printing'} = 1;
      }
      if ( $$optsref{'source'} eq "default" ) {
         $inputopts{'datasource'}   = $self->default_vals( "source" );
      } else {
         $inputopts{'datasource'}   = $$optsref{'source'};
      }
      my $source = $inputopts{'datasource'};
      if ( $source ne "CLMNCEP" && $source ne "CAMHIST.eul64x128_datm6.01" &&
           $source ne "CAMHIST.GENERIC" ) {
         die "ERROR($nm):: Bad source option -- need to specify source as one of: " .
             "CLMNCEP, CAMHIST.eul64x128_datm6.01, or CAMHIST.GENERIC";
      }
      $inputopts{'ProgName'}   = $$optsref{'ProgName'};
      $inputopts{'ProgDir'}    = "$cfgdir";
      $inputopts{'res'}        = "";
      $inputopts{'yearfirst'}  = $beg_year;
      $inputopts{'yearlast'}   = $end_year;
      #my @keys = ( "datm_data_dir", "datm_dom_dir", "cam_hist_case", "datm_domain" );
      #foreach my $key ( @keys ) {
      #   if ( defined($$optsref{$key}) ) {
      #      $inputopts{'filepath'} = $$optsref{$key};
      #   } else {
      #      if ( $source ne "CAMHIST.eul64x128_datm6.01" ) {
      #         die "ERROR($nm):: Need to specify datm_data_dir directory of where " .
      #             " datm data is for datasource=$source.\n";
      #      }
      #      $inputopts{'filepath'}   = "";
      #   }
      #}
      if ( defined($$optsref{'datm_data_dir'}) ) {
         $inputopts{'filepath'}   = $$optsref{'datm_data_dir'};
      } else {
         if ( $source ne "CAMHIST.eul64x128_datm6.01" ) {
            die "ERROR($nm):: Need to specify datm_data_dir directory of where " .
                " datm data is for datasource=$source.\n";
         }
         $inputopts{'filepath'}   = "";
      }
      if ( ! defined($$optsref{'datm_dom_dir'} ) ) {
         $inputopts{'domainpath'} = "";
         if ( $source eq "CAMHIST.GENERIC" ) {
            die "ERROR($nm):: Need to specify datm_dom_dir directory of where " .
                " datm domain data is for datasource=$source.\n";
         }
      } else {
         $inputopts{'domainpath'} = $$optsref{'datm_dom_dir'};
      }
      if ( ! defined($$optsref{'cam_hist_case'}) ) {
         $inputopts{'case'}    = "";
         if ( $source eq "CAMHIST.GENERIC" ) {
            die "ERROR($nm):: Need to specify datm_dom_dir directory of where " .
                " datm domain data is for datasource=$source.\n";
         }
      } else {
         $inputopts{'case'}    = $$optsref{'cam_hist_case'};
      }
      if ( ! defined($$optsref{'datm_domain'}) ) {
         $inputopts{'domain'} = "";
         if ( $source eq "CAMHIST.GENERIC" ) {
            die "ERROR($nm):: Need to specify datm_domain name of domain file " .
                " for datasource=$source.\n";
         }
      } else {
          $inputopts{'domain'} = $$optsref{'datm_domain'};
      }
      $inputopts{'cmdline'}    = $$optsref{'cmdline'};
      $inputopts{'csmdata'}    = $csmdata;
      $inputopts{'filenames'}  = "";

      my $streams = Streams::Template->new( \%inputopts );
      $streams->Read( "$cfgdir/$template" );
      if ( $$optsref{'test'} ) {
         $streams->TestFilesExist( "data" );
         $streams->TestFilesExist( "domain" );
      }
      $streams->Write( $outstreams );

      $NLref->{'streams'} = "'$outstreams $init_year $beg_year $end_year'";
  }
  $self->requiredVar( "streams" );

  unless (defined($NLref->{'infodbug'})) {
      $NLref->{'infodbug'} = $self->default_vals( "infodbug" );
  }

  unless (defined($NLref->{'datamode'})) {
      $NLref->{'datamode'} =namelist::quote_string( $$settings_ref{'datamode'} );
  }
  if ($NLref->{'datamode'} !~ /$$settings_ref{'datamode'}/ ) {
    die "ERROR($nm):: datamode inconsistent between file and input command line options\n";
  }
  $self->requiredVar( "datamode" );

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
