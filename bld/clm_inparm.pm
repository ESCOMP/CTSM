#
#	clm_inparm.pm			Erik Kluzek
#
#	Perl module to create a namelist for CLM2.
#
#	Description of methods:
#
#	new ----------------------- Constructor
#	set_output_values --------  Set output values based on precedence of the various input
#                                   values and ensure that a valid namelist is produced.
#
#-----------------------------------------------------------------------------------------------

use strict;
#use diagnostics;
use Cwd;

package clm_inparm;

use compnl;
@clm_inparm::ISA = qw(compnl  namelist);
%clm_inparm::NL;

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
     $outfile = "lnd_in";
  } else  {
     $outfile = "lnd.stdin";
  }

  my $file     = $$optsref{'dir'} . "/$outfile";
  my $printlev = $$optsref{'printlev'};
  my $cfgdir   = $$optsref{'cfgdir'};

  my $self = $class->SUPER::new( "clm_inparm", \%clm_inparm::NL, $interactive, $file,
                                 "$cfgdir/DefaultCLM_INPARM_Namelist.xml", $config, 
                                 $printlev );

  $self->{'printlev'} = $printlev;
  $self->{'optsref'}  = $optsref;
  $self->{'MODE'}     = $bld{'MODE'};
  $self->{'RTM'}      = $bld{'RTM'};
  $self->{'BGC'}      = $bld{'BGC'};

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub set_output_values {

# Set the CLM2 namelist variables.

  my ($self, $settings_ref) = @_;

  my $start_type = $$settings_ref{start_type};
  my $class      = ref($self);
  my $nm         = "$class\:\:set_output_values";

  my $NLref      = $self->{'NLREF'};
  my $optsref    = $self->{'optsref'};
  my $opt;
  my $raw_files_spec;

  my $noquote   = 0;
  my $use_land  = 1;
  my $name      = $self->{'NAME'};
  my $addquotes = 1;

  # If not seq-ccsm`

  if ( $self->{'MODE'} ne "ccsm_seq" ) {

     # start_ymd

     if (defined($NLref->{'start_ymd'})) {
         $opt = $NLref->{'start_ymd'};
     } else {
         $opt = 19990101;
     }
     $NLref->{'start_ymd'} = $opt;
     $self->requiredVar( "start_ymd" );
     $$settings_ref{'ic_ymd'} = $opt;

     # start_tod

     if (defined($NLref->{'start_tod'})) {
         $opt = $NLref->{'start_tod'};
     } else {
         $opt = 0;
     }
     $NLref->{'start_tod'}    = $opt;
     $$settings_ref{'ic_tod'} = $opt;
  }

  # Read in XML defaults file
  $self->Read_Default_XML_Values( $settings_ref );

  my $default_vals = $self->{'default_vals'};

  # Check that "nrevsn" is set if this is a branch simulation
  if ($start_type eq 'branch') {
     $self->requiredVar( "nrevsn", $addquotes );
  }

  # Plant function types.
  if ( ! defined($NLref->{'fpftcon'})) {
      $opt = $self->default_vals( 'fpftcon' );
  } else {
      $opt = $NLref->{'fpftcon'};
  }
  $NLref->{'fpftcon'} = namelist::quote_string($opt);
  $self->requiredVar( "fpftcon", $addquotes );
  $self->checkinputfile('fpftcon') if $optsref->{'test'};

  # Surface dataset
  if ( ! defined($NLref->{'fsurdat'})) {
      my $val = $self->default_vals( 'fsurdat' );
      $opt = $self->default_vals( 'fsurdat', $noquote, $use_land );
  } else {
      $opt = $NLref->{'fsurdat'};
  }
  $NLref->{'fsurdat'} = namelist::quote_string($opt);
  $self->requiredVar( "fsurdat", $addquotes );
  $self->checkinputfile('fsurdat') if $optsref->{'test'};

  # Time-step
  unless (defined($NLref->{'dtime'})) {
      if ( defined($$settings_ref{'dtime'}) ) {
	  $NLref->{'dtime'} = $$settings_ref{dtime};
      } else {
          $NLref->{'dtime'} = $self->default_vals( 'dtime' );
      }
  }
  $self->requiredVar( "dtime", $addquotes );

  # Grid
  if ( ! defined($NLref->{'fatmgrid'}) ) {
      $opt = $self->default_vals( 'fatmgrid' );
  } else {
      $opt = $NLref->{'fatmgrid'};
  }
  $NLref->{'fatmgrid'} = namelist::quote_string($opt);
  $self->requiredVar( "fatmgrid", $addquotes );
  $self->checkinputfile('fatmgrid') if $optsref->{'test'};

  # Nitrogen deposition
  if ( $self->{'BGC'} eq "cn" ) {
     unless ( defined($NLref->{'fndepdat'}) ) {
         $opt  = $self->default_vals( 'fndepdat', $noquote, $use_land );
         if ( defined($opt) ) {
            $NLref->{'fndepdat'} = namelist::quote_string($opt);
         }
     }
     if ( defined $NLref->{'fndepdat'} and $optsref->{'test'} ) {
         $self->checkinputfile('fndepdat');
     }
  }

  # RTM
  if ( $self->{'RTM'} eq "on" ) {
     unless ( defined($NLref->{'frivinp_rtm'}) ) {
         $opt                    = $self->default_vals( 'frivinp_rtm' );
         $NLref->{'frivinp_rtm'} = namelist::quote_string($opt);
     }
     $self->requiredVar( "frivinp_rtm", $addquotes );
     if ( defined $NLref->{'frivinp_rtm'} and $optsref->{'test'} ) {
         $self->checkinputfile('frivinp_rtm');
     }
  }

  # CO2
  unless ( defined($NLref->{'co2_ppmv'}) ) {
      $NLref->{'co2_ppmv'} = $self->default_vals( 'co2_ppmv' );
  }
  $self->requiredVar( "co2_ppmv", $addquotes );

  # Land fraction
  unless ( defined($NLref->{'fatmlndfrc'}) ) {
      $opt                   = $self->default_vals( 'fatmlndfrc' );
      $NLref->{'fatmlndfrc'} = namelist::quote_string($opt);
  }
  $self->requiredVar( "fatmlndfrc", $addquotes );
  $self->checkinputfile('fatmlndfrc')  if $optsref->{'test'};

  #
  # Fine mesh grid
  #
  if ( $$optsref{'lnd_RESOLUTION'} ne $$optsref{'RESOLUTION'} ) {
     # Land topography
     unless ( defined($NLref->{'flndtopo'}) ) {
         $opt                 = $self->default_vals( 'flndtopo', $noquote, $use_land );
         $NLref->{'flndtopo'} = namelist::quote_string($opt);
     }
     $self->requiredVar( "flndtopo", $addquotes );
     $self->checkinputfile('flndtopo')  if $optsref->{'test'};
     # Atmosphere topography
     unless ( defined($NLref->{'fatmtopo'}) ) {
         $opt                 = $self->default_vals( 'fatmtopo' );
         $NLref->{'fatmtopo'} = namelist::quote_string($opt);
     }
     $self->requiredVar( "fatmtopo", $addquotes );
     $self->checkinputfile('fatmtopo')  if $optsref->{'test'};
  }

  # If not seq-ccsm`

  if ( $self->{'MODE'} ne "ccsm_seq" ) {

     # nsrest

     my %start = ( arb_ic=>0, startup=>0, continue=>1, branch=>3 );
     if (defined($optsref->{'start_type'})) {
         my $key = $optsref->{'start_type'};
         $opt = $start{$key};
     } elsif (defined($NLref->{'nsrest'})) {
         $opt = $NLref->{'nsrest'};
     } else {
         $opt = 0;
     }
     $NLref->{'nsrest'} = $opt;
     $self->requiredVar( "nsrest", $addquotes );

     # Case name
     if (defined($optsref->{'case'})) {
         $opt = $optsref->{'case'};
     } elsif (defined($NLref->{'caseid'})) {
         $opt = $NLref->{'caseid'};
     } else {
         $opt = "\'clmrun\'";
     }
     $NLref->{'caseid'} = namelist::quote_string($opt);
     $self->requiredVar( "caseid", $addquotes );

     # nelapse
     if ( ! defined($NLref->{'nestep'})) {
        if (defined($NLref->{'nelapse'})) {
            $opt = $NLref->{'nelapse'};
        } else {
            if ( $self->{'MODE'} eq "ext_ccsm_con" ) {
               $opt = -9999;
            } else {
               $opt = -1;
            }
        }
        $NLref->{'nelapse'} = $opt;
        $self->requiredVar( "nelapse" );
     }


     #
     # Concurrent CCSM only
     #
     if (      $self->{'MODE'} eq "ext_ccsm_con" ) {
        # do flux averaging
        if (defined($NLref->{'csm_doflxave'})) {
            $opt = $NLref->{'csm_doflxave'};
        } else {
            $opt = ".true.";
        }
        $NLref->{'csm_doflxave'} = $opt;

     #
     # Offline only
     #
     } elsif ( $self->{'MODE'} eq "offline" ) {
        # offline_atmdir
        if (defined($NLref->{'offline_atmdir'})) {
            $opt = $NLref->{'offline_atmdir'};
        } else {
            $opt = $self->default_vals( 'offline_atmdir' );
        }
        $NLref->{'offline_atmdir'} = namelist::quote_string($opt);
        $self->requiredVar( "offline_atmdir", $addquotes );
        # cycle_begyr
        if (defined($optsref->{'cycle_beg_year'})) {
            $opt = $optsref->{'cycle_beg_year'};
            if (defined($NLref->{'cycle_begyr'})) {
	      die "ERROR: define cycle_begyr using -cycle_beg_year command line option rather than in namelist.\n";
            }
        }
        $NLref->{'cycle_begyr'} = $opt;
        # cycle_nyr
        if (defined($optsref->{'cycle_end_year'})) {
            $opt = $optsref->{'cycle_end_year'} - $NLref->{'cycle_begyr'} + 1;
            if (defined($NLref->{'cycle_nyr'})) {
	      die "ERROR: define cycle_nyr using -cycle_beg_year and cycle_end_year command line options rather than in namelist.\n";
            }
        }
        $NLref->{'cycle_nyr'} = $opt;
     }
  }

  # Initial conditions
  if ( $optsref->{'start_type'} eq "arb_ic" ) {
     if ( ! defined($NLref->{'finidat'}) ) {
        $NLref->{'finidat'} = "' '";
     } else {
        if ( $NLref->{'finidat'} !~ /^[' "]+$/ ) {
	   die "ERROR: start_type set to arb_ic and yet finidat is also set\n";
        }
     }
  } else {
     unless ( defined($NLref->{'finidat'}) ) {
        $opt                = $self->default_vals( 'finidat', $noquote, $use_land );
        if ( ! defined($opt) ) { $opt = "' '"; }
        $NLref->{'finidat'} = namelist::quote_string($opt);
     }
     $self->requiredVar( "finidat", $addquotes ) if $optsref->{'start_type'} eq "startup";
     $self->checkinputfile('finidat') if $optsref->{'test'};
  }
  # Make sure clm_demand list is filled
  if ( defined($optsref->{'clm_demand'}) ) {
     my @demandlist = split( ",", $optsref->{'clm_demand'} );
     foreach my $item ( @demandlist ) {
        unless ( defined($NLref->{$item}) ) {
           $NLref->{$item} = $self->default_vals( $item );
        }
        $self->requiredVar( $item );
        my $isafile = 0;
        if ( $item eq "fpftdyn"  ) { $isafile = 1; }
        if ( $item eq "fndepdyn" ) { $isafile = 1; }
        if ( $isafile ) {
           $NLref->{$item} = namelist::quote_string( $NLref->{$item} );
           $self->checkinputfile($item) if $optsref->{'test'};
        }
     }
  }

}

#============================================================================

1   # to make use or require happy
