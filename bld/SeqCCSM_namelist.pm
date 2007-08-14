#
#	SeqCCSM_namelist.pm			Erik Kluzek
#
#	Perl module to create, parse and check the namelists needed
#	for Sequential CCSM.
#
#------------------------------------------------------------------------
#
#	Description of methods:
#
#	new ----------------------- Constructor
#	build --------------------- Create the namelist
#	interactive --------------- Interactively change the values in the namelist
#	parse --------------------- Parse a previous namelist.
#	print --------------------- Print contents of namelist to terminal.
#
#	$Id$
#
use strict;
#use diagnostics;
use Cwd;

package SeqCCSM_namelist;
use drv_in;
use ccsm_inparm;
use timemgr_inparm;
use prof_inparm;
use clm_inparm;
use datm_dshr_in;
use datm_nml;
use queryDefaultXML;

sub new {
#
# Constructor
#
  my $class = shift;
  my $optsref = shift;

  my $self = {};
  my %opts = %$optsref;
  $self->{'DRVNL'}     = undef;                # Drv_in "virtual" namelist object
  $self->{'CSMNL'}     = undef;                # CCSM namelist object
  $self->{'TIMNL'}     = undef;                # Time-manager namelist object
  $self->{'PROFNL'}    = undef;                # Performance profiling namelist object
  $self->{'ATMNL'}     = undef;                # Atmospheric model namelist object
  $self->{'ATMSHRNL'}  = undef;                # Atmospheric share model namelist object
  $self->{'LNDNL'}     = undef;                # Land-model namelist object

  $self->{'PARSE_FILE'} = $opts{'infile'};     # User supplied input namelist file
  if ( $opts{namelist} ) {                     # namelist settings from the command-line
    my $filename = ".tmpfile.";
    open( TMPNL, ">$filename" ) || die "ERROR:: Can not open temp file : $filename\n";
    my $namelist = $opts{namelist};
    print TMPNL "$namelist\n";
    close( TMPNL );
    $self->{'PARSE_NLLINE'} = $filename;
  } else {
    $self->{'PARSE_NLLINE'} = undef;              # Temporary-Filename of command-line 
  }

  $self->{'INTERACTIVE'} = $opts{interactive};    # Whether to use interactive prompting or not

  $self->{'OUTDIR'} = $opts{dir};                 # Output directory for namelists

  my $config_file      = $opts{config};     # config_cache.xml filename

  $self->{'optsref'} = $optsref;

  $self->{'printlev'} = $opts{'printlev'};   # Print level
  my %inputopts;
  my $cfgdir             = $$optsref{'cfgdir'};
  $inputopts{'file'}     = "$cfgdir/DefaultSettings.xml";
  $inputopts{'namelist'} = "default_settings";
  $inputopts{'printing'} = $self->{'printlev'};
  $inputopts{'ProgName'} = $$optsref{'ProgName'};
  $inputopts{'cmdline'}  = $$optsref{'cmdline'};    
  $inputopts{'csmdata'}  = $$optsref{'csmdata'};
  $inputopts{'config'}   = $config_file;
  $inputopts{'res'}      = $$optsref{'RESOLUTION'};
  my %settings;
  my $default_ref = &queryDefaultXML::ReadDefaultXMLFile( \%inputopts, \%settings );
  my %defaults = %$default_ref;
  foreach my $i ( keys(%defaults) ) {
     if ( $$optsref{$i} eq "default" ) {
        $$optsref{$i} = $defaults{$i}{'value'};
     }
  }
  if ( $$optsref{'lnd_RESOLUTION'} eq "default" ) {
     $$optsref{'lnd_RESOLUTION'} = $$optsref{'RESOLUTION'};
  }

  my %bld = &queryDefaultXML::read_cfg_file( $$optsref{'config'} );
  $self->{'MODE'} = $bld{'MODE'};


  if ( $self->{'MODE'} eq "ccsm_seq" ) {

     $self->{'CSMNL'} = ccsm_inparm->new( $self->{'optsref'}, $config_file );

     # Timemanager namelist
     $self->{'TIMNL'} = timemgr_inparm->new( $self->{'optsref'}, $config_file );

     # Timing Profile namelist
     $self->{'PROFNL'} = prof_inparm->new( $self->{'optsref'}, $config_file );

     # atm namelist
     $self->{'ATMNL'}    = datm_nml->new( $self->{'optsref'}, $config_file );
     $self->{'ATMSHRNL'} = datm_dshr_in->new( $self->{'optsref'}, $config_file );
  } else {
     # Timing Profile namelist
     $self->{'PROFNL'} = prof_inparm->new( $self->{'optsref'}, $config_file );
  }

  # Land-model namelist
  $self->{'LNDNL'} = clm_inparm->new( $self->{'optsref'}, $config_file );

  # Driver "virtual" namelist that will put info in the previous namelists
  $self->{'DRVNL'} = drv_in->new( $self->{'optsref'}, $self->{'LNDNL'},
                                  $self->{'PROFNL'},  $self->{'CSMNL'},
                                  $self->{'TIMNL'},   $self->{'ATMSHRNL'} );

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub interactive {
#
# Ask for additional settings interactively
#
  my $self = shift;

  if ( $self->{'MODE'} eq "ccsm_seq" ) {
     $self->{'CSMNL'}->change;
     $self->{'TIMNL'}->change;
     $self->{'ATMNL'}->change;
     $self->{'ATMSHRNL'}->change;
  }
  $self->{'PROFNL'}->change;
  $self->{'LNDNL'}->change;
}

#============================================================================

sub print {
#
# Print out the resulting namelist
#
  my $self = shift;

  if ( $self->{'MODE'} eq "ccsm_seq" ) {
     $self->{CSMNL}->print;
     $self->{TIMNL}->print;
     $self->{ATMNL}->print;
     $self->{ATMSHRNL}->print;
  }
  $self->{PROFNL}->print;
  $self->{LNDNL}->print;
}

#============================================================================

sub parse {
#
# Parse a namelist
#
  my $self = shift;
  my $file = shift;

  $self->{DRVNL}->parse(   $file );
  if ( $self->{'MODE'} eq "ccsm_seq" ) {
     $self->{CSMNL}->parse(    $file );
     $self->{TIMNL}->parse(    $file );
     $self->{ATMNL}->parse(    $file );
     $self->{ATMSHRNL}->parse( $file );
  }
  $self->{PROFNL}->parse(   $file );
  $self->{LNDNL}->parse(    $file );
}

#============================================================================

sub build {
#
# Build the namelist
#
  my $self = shift;

  my $class = ref($self);
  my $nm = "$class\:\:build";

  my $optsref = $self->{'optsref'};

  # Get values from user specified namelist file
  if ( defined( $self->{'PARSE_FILE'} ) ) {
    $self->parse( $self->{'PARSE_FILE'});
  }

  # Get values from namelist specified on the command-line (these values
  # will overwrite the ones read from the namelist file).
  if ( defined( $self->{'PARSE_NLLINE'} ) ) {
    $self->parse( $self->{'PARSE_NLLINE'});
    my $cmd = "/bin/rm -rf " . $self->{'PARSE_NLLINE'};
    `$cmd`;
  }

  # set all keys to lower case (this should be done where the keys are set)
  $self->{'DRVNL'}->convert_case;
  if ( $self->{'MODE'} eq "ccsm_seq" ) {
     $self->{'CSMNL'}->convert_case;
     $self->{'TIMNL'}->convert_case;
     $self->{'ATMNL'}->convert_case;
     $self->{'ATMSHRNL'}->convert_case;
  }
  $self->{'PROFNL'}->convert_case;
  $self->{'LNDNL'}->convert_case;

  # Run type is set in SeqCCSM namelist.  Land model namelist only needs to know
  # the run type to check that "nrevsn" is set for a branch simulation.
  my %settings;
  $settings{start_type} = $$optsref{start_type};
  $settings{sim_year}   = $$optsref{sim_year};
  $settings{mask}       = $$optsref{mask};
  $settings{ic_ymd}     = undef;
  $settings{dtime}      = undef;

  $settings{ignore_ic_date} = $$optsref{ignore_ic_date};
  $settings{ignore_ic_year} = $$optsref{ignore_ic_year};

  $self->{DRVNL}->set_output_values( );
  if ( $self->{'MODE'} eq "ccsm_seq" ) {
     # Set output values according to precedence of various input values.  Ensure
     # that valid output namelists are produced.
     $self->{'CSMNL'}->set_output_values( \%settings );

     if ( $$optsref{'runlength'} ne "default" ) {
        if ( $$optsref{'runlength'} =~ /^([1-9][0-9]*)([sdy])$/ ) {
           my %options = ( s=>"nsteps", d=>"ndays", y=>"nyears" );
           $self->{'TIMNL'}->SetValue( "stop_n",      $1 );
           $self->{'TIMNL'}->SetValue( "stop_option", $options{$2} );
        } else {
           die "ERROR($nm):: bad input to runlength option\n";
        }
     }
     $self->{'TIMNL'}->set_output_values( \%settings );
     # Get the date and time of day that the initial condition file should have
     $settings{ic_ymd} = $self->{'TIMNL'}->Value( "start_ymd" );
     $settings{ic_tod} = $self->{'TIMNL'}->Value( "start_tod" );
     if ( ! defined($settings{ic_tod}) ) {
        $settings{ic_tod} = 0;
     }
     if ( $$optsref{'cycle_init_year'} eq "default" ) {
        my $ymd = $self->{'TIMNL'}->Value( "start_ymd" );
        $$optsref{'cycle_init_year'} = int( $ymd / 10000 );
     }

     $self->{'ATMNL'}->set_output_values(    \%settings );
     $self->{'ATMSHRNL'}->set_output_values( \%settings );
  } else {
     if ( $$optsref{'runlength'} ne "default" ) {
        if ( $$optsref{'runlength'} =~ /^([1-9][0-9]*)([sdy])$/ ) {
           my $length = $1;
           if (      $2 eq "d" ) {
              $length *= -1;
           } elsif ( $2 eq "y" ) {
              $length *= -365;
           }
           $self->{'LNDNL'}->SetValue( "nelapse", $length );
        } else {
           die "ERROR($nm):: bad input to runlength option\n";
        }
     }
  }
  $self->{'LNDNL'}->set_output_values(  \%settings );

  # Performance profiling namelist
  $self->{'PROFNL'}->set_output_values( \%settings );

  # Allow for additional settings interactively
  if ( defined($self->{'INTERACTIVE'}) && $self->{'INTERACTIVE'} ) {
    $self->interactive;
  }
  
  # Write the final namelists to the output namelist files.
  my $OUTDIR = $self->{OUTDIR};
  my $eol = $$optsref{'eol'};
  if ( $self->{'printlev'} ) { print "Write out namelists to directory: $OUTDIR $eol"; }


  my %bld = &queryDefaultXML::read_cfg_file( $$optsref{'config'} );
  if ( ($bld{'MODE'} eq "ext_ccsm_con")  && defined($self->{'csmdata'}) ) { 
      $self->{'LNDNL'}->Write_prestage(  $self->{'csmdata'} );
      $self->{'PROFNL'}->Write_prestage( $self->{'csmdata'}, 'Append' );
  } else {
      if ( $bld{'MODE'} eq "ccsm_seq") {
         # Driver namelists
         $self->{'CSMNL'}->Write;
         $self->{'TIMNL'}->Write(  'Append' );
         # Atmosphere namelists
         $self->{'ATMNL'}->Write;
         $self->{'ATMSHRNL'}->Write;
      }
      # Land namelist
      $self->{'LNDNL'}->Write;
      $self->{'PROFNL'}->Write( 'Append' );
  }
}

1   # to make use or require happy
