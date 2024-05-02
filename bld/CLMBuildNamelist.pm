# build-namelist
#
# This script builds the namelists for CLM
#
# The simplest use of build-namelist is to execute it from the build directory where configure
# was run.  By default it will use the config_cache.xml file that was written by configure to
# determine the build time properties of the executable, and will write the files that contain
# the output namelists in that same directory.  But if multiple runs are to made using the
# same executable, successive invocations of build-namelist will overwrite previously generated
# namelist files.  So generally the best strategy is to invoke build-namelist from the run
# directory and use the -config option to provide the filepath of the config_cache.xml file.
#
#
# Date        Contributor      Modification
# -------------------------------------------------------------------------------------------
# 2009-01-20  Vertenstein      Original version
# 2010-04-27  Kluzek           Add ndep streams capability
# 2011-07-25  Kluzek           Add multiple ensemble's of namelists
# 2012-03-23  Kluzek           Add megan namelist and do checking on it
# 2012-07-01  Kluzek           Add some common CESM namelist options
# 2013-12     Andre            Refactor everything into subroutines
# 2013-12     Muszala          Add Ecosystem Demography functionality
#--------------------------------------------------------------------------------------------

package CLMBuildNamelist;

require 5;

use strict;
#use warnings;
#use diagnostics;

use Cwd qw(getcwd abs_path);
use File::Basename qw(dirname);
use English;
use Getopt::Long;
use IO::File;
use File::Glob ':bsd_glob';

#-------------------------------------------------------------------------------
#
# Define a small number of global variables
#
#-------------------------------------------------------------------------------

(my $ProgName = $0) =~ s!(.*)/!!;      # name of this script
my $ProgDir  = $1;
$ProgName = "CLM " . "$ProgName";

my $cwd = abs_path(getcwd());  # absolute path of the current working directory
my $log;                       # Log messages object -- will be set in main, declaring it global here means it can be used everywhere

#-------------------------------------------------------------------------------

sub usage {
    die <<EOF;
SYNOPSIS
     build-namelist [options]

     Create the namelist for CLM
REQUIRED OPTIONS
     -cimeroot "directory"    Path to cime directory
     -config "filepath"       Read the given CLM configuration cache file.
                              Default: "config_cache.xml".
     -configuration "cfg"     The overall configuration being used [ clm | nwp ]
                                clm = Climate configuration
                                nwp = Numerical Weather Prediction configuration
     -d "directory"           Directory where output namelist file will be written
                              Default: current working directory.
     -envxml_dir "directory"  Directory name of env_*.xml case files to read in.
                              (if read they allow user_nl_clm and CLM_BLDNML_OPTS to expand
                               variables [for example to use \$DIN_LOC_ROOT])
                              (default current directory)
     -lnd_frac "domainfile"   Land fraction file (the input domain file) (needed for MCT driver and LILAC)
     -res "resolution"        Specify horizontal grid.  Use nlatxnlon for spectral grids;
                              dlatxdlon for fv grids (dlat and dlon are the grid cell size
                              in degrees for latitude and longitude respectively)
                              "-res list" to list valid resolutions.
                              (default: 0.9x1.25)
     -sim_year "year"         Year to simulate for input datasets
                              (i.e. PtVg, 1850, 2000, 2010, 1850-2000, 1850-2100)
                              "-sim_year list" to list valid simulation years
                              (default 2000)
     -structure "structure"   The overall structure being used [ standard | fast ]
OPTIONS
     -driver "value"          CESM driver type you will run with [ mct | nuopc ]
     -bgc "value"             Build CLM with BGC package [ sp | bgc | fates ]
                              (default is sp).
                                CLM Biogeochemistry mode
                                sp    = Satellite Phenology (SP)
                                    This toggles off the namelist variable: use_cn
                                bgc   = Carbon Nitrogen with methane, nitrification, vertical soil C,
                                        CENTURY or MIMICS decomposition
				        This toggles on the namelist variables:
                                          use_cn, use_lch4, use_nitrif_denitrif
				fates = FATES/Ecosystem Demography with below ground BGC
				        CENTURY or MIMICS decomposition
                                        This toggles on the namelist variables:
				        use_fates. use_lch4 and use_nitrif_denitrif are optional
				 
                              (Only for CLM4.5/CLM5.0)
     -[no-]chk_res            Also check [do NOT check] to make sure the resolution and
                              land-mask is valid.
     -clm_accelerated_spinup "on|off" Setup in a configuration to run as fast as possible for doing a throw-away
                              simulation in order to get the model to a spun-up state. So do things like
                              turn off expensive options and setup for a low level of history output.

                              If CLM4.5/CLM5.0 and bgc it also includes a prognostic Carbon model (cn or bgc)
                              , also by default turn on Accelerated Decomposition mode which
                              is controlled by the namelist variable spinup_state.

                              Turn on given spinup mode for BGC setting of CN
                                  on : Turn on Accelerated Decomposition   (spinup_state = 1 or 2)
                                  off : run in normal mode                 (spinup_state = 0)

                              Default is set by clm_accelerated_spinup mode.

                              Spinup is now a two step procedure. First, run the model
                              with clm_accelerated_spinup = "on". Then run the model for a while with
                              spinup_state = 0. The exit spinup step happens automatically
                              on the first timestep when using a restart file from spinup
                              mode.

                              The spinup state is saved to the restart file.
                              If the values match between the model and the restart
                              file it proceeds as directed.

                              If the restart file is in spinup mode and the model is in
                              normal mode, then it performs the exit spinup step
                              and proceeds in normal mode after that.

                              If the restart file has normal mode and the model is in
                              spinup, then it enters spinup. This is useful if you change
                              a parameter and want to rapidly re-equilibrate without doing
                              a cold start.

     -clm_demand "list"       List of variables to require on clm namelist besides the usuals.
                              "-clm_demand list" to list valid options.
                              (can include a list member "null" which does nothing)
     -clm_start_type "type"   Start type of simulation
                              (default, cold, arb_ic, startup, continue, or branch)
                              (default=do the default type for this configuration)
                              (cold=always start with arbitrary initial conditions)
                              (arb_ic=start with arbitrary initial conditions if
                               initial conditions do not exist)
                              (startup=ensure that initial conditions are being used)
     -clm_usr_name     "name" Dataset resolution/descriptor for personal datasets.
                              Default: not used
                              Example: 1x1pt_boulderCO_c090722 to describe location,
                                       number of pts, and date files created
     -co2_type "value"        Set CO2 the type of CO2 variation to use.
     -co2_ppmv "value"        Set CO2 concentration to use when co2_type is constant (ppmv).
     -crop                    Toggle for prognostic crop model. (default is off)
                              (can ONLY be turned on when BGC type is CN or BGC)
                              This turns on the namelist variable: use_crop
     -csmdata "dir"           Root directory of CESM input data.
                              Can also be set by using the CSMDATA environment variable.
     -drydep                  Produce a drydep_inparm namelist that will go into the
                              "drv_flds_in" file for the driver to pass dry-deposition to the atm.
                              Default: -no-drydep
                              (Note: buildnml copies the file for use by the driver)
     -dynamic_vegetation      Toggle for dynamic vegetation model. (default is off)
                              (can ONLY be turned on when BGC type is 'bgc')
                              This turns on the namelist variable: use_cndv
                              (Deprecated, this will be removed)
     -fire_emis               Produce a fire_emis_nl namelist that will go into the
                              "drv_flds_in" file for the driver to pass fire emissions to the atm.
                              (Note: buildnml copies the file for use by the driver)
     -glc_nec <name>          Glacier number of elevation classes [0 | 3 | 5 | 10 | 36]
                              (default is 0) (standard option with land-ice model is 10)
     -glc_use_antarctica      Set defaults appropriate for runs that include Antarctica
     -help [or -h]            Print usage to STDOUT.
     -light_res <value>       Resolution of lightning dataset to use for CN or FATES fire (360x720, 106x174, or 94x192)
                              106x174 can only be used for NEON sites
     -lilac                   If CTSM is being run through LILAC (normally not used)
                              (LILAC is the Lightweight Infrastructure for Land-Atmosphere Coupling)
     -ignore_ic_date          Ignore the date on the initial condition files
                              when determining what input initial condition file to use.
     -ignore_ic_year          Ignore just the year part of the date on the initial condition files
                              when determining what input initial condition file to use.
     -ignore_warnings         Allow build-namelist to continue, rather than stopping on
                              warnings
     -infile "filepath"       Specify a file (or list of files) containing namelists to
                              read values from.

                              If used with a CLM build with multiple ensembles (ninst_lnd>1)
                              and the filename entered is a directory to files of the
                              form filepath/filepath and filepath/filepath_\$n where \$n
                              is the ensemble member number. the "filepath/filepath"
                              input namelist file is the master input namelist file
                              that is applied to ALL ensemble members.

                              (by default for CESM this is setup for files of the
                               form \$CASEDIR/user_nl_clm/user_nl_clm_????)
     -inputdata "filepath"    Writes out a list containing pathnames for required input datasets in
                              file specified.
     -lnd_tuning_mode "value" Use the parameters tuned for the given configuration (CLM version and atmospheric forcing)
     -mask "landmask"         Type of land-mask (default, navy, gx3v5, gx1v5 etc.)
                              "-mask list" to list valid land masks.
     -namelist "namelist"     Specify namelist settings directly on the commandline by supplying
                              a string containing FORTRAN namelist syntax, e.g.,
                                 -namelist "&clm_inparm dt=1800 /"
     -no-megan                DO NOT PRODUCE a megan_emis_nl namelist that will go into the
                              "drv_flds_in" file for the driver to pass VOCs to the atm.
                              MEGAN (Model of Emissions of Gases and Aerosols from Nature)
                              (Note: buildnml copies the file for use by the driver)
     -[no-]note               Add note to output namelist  [do NOT add note] about the
                              arguments to build-namelist.
     -output_reals <file>     Output real parameters to the given output file.
     -ssp_rcp "value"         Shared Socioeconomic Pathway (SSP) and
                              Representative Concentration Pathway (RCP) combination to use for
                              future scenarios.
                              "-ssp_rcp list" to list valid ssp_rcp settings.
     -s                       Turns on silent mode - only fatal messages issued.
     -test                    Enable checking that input datasets exist on local filesystem.
     -use_case "case"         Specify a use case which will provide default values.
                              "-use_case list" to list valid use-cases.
     -verbose [or -v]         Turn on verbose echoing of informational messages.
     -version                 Echo the SVN tag name used to check out this CLM distribution.
     -vichydro                Toggle to turn on VIC hydrologic parameterizations (default is off)
                              This turns on the namelist variable: use_vichydro


Note: The precedence for setting the values of namelist variables is (highest to lowest):
      0. namelist values set by specific command-line options, like, -d, -sim_year
             (i.e.  compset choice and CLM_BLDNML_OPTS, CLM_ACCELERATED_SPINUP, LND_TUNING_MODE env_run variables)
     (NOTE: If you try to contradict these settings by methods below, an error will be triggered)
      1. values set on the command-line using the -namelist option,
             (i.e. CLM_NAMELIST_OPTS env_run variable)
      2. values read from the file(s) specified by -infile,
             (i.e.  user_nl_clm files)
      3. datasets from the -clm_usr_name option,
             (i.e.  CLM_USRDAT_NAME env_run variable)
      4. values set from a use-case scenario, e.g., -use_case
             (i.e.  CLM_NML_USE_CASE env_run variable)
      5. values from the namelist defaults file.
EOF
}

#-------------------------------------------------------------------------------

sub process_commandline {
  # Process command-line options and return the hash
  my ($nl_flags) = @_;

  # Save the command line arguments to the script. NOTE: this must be
  # before GetOptions() is called because items are removed from from
  # the array!
  $nl_flags->{'cmdline'} = "@ARGV\n";

  my %opts = ( cimeroot              => undef,
               config                => "config_cache.xml",
               configuration         => undef,
               csmdata               => undef,
               clm_usr_name          => undef,
               co2_type              => undef,
               co2_ppmv              => undef,
               clm_demand            => "null",
               driver                => "nuopc",
               help                  => 0,
               glc_nec               => "default",
               glc_use_antarctica    => 0,
               light_res             => "default",
               lnd_tuning_mode       => "default",
               lnd_frac              => undef,
               dir                   => "$cwd",
               ssp_rcp               => "default",
               sim_year              => "default",
               structure             => undef,
               clm_accelerated_spinup=> "default",
               chk_res               => undef,
               note                  => undef,
               drydep                => 0,
               lilac                 => 0,
               output_reals_filename => undef,
               fire_emis             => 0,
               megan                 => "default",
               res                   => "default",
               silent                => 0,
               ignore_warnings       => 0,
               mask                  => "default",
               test                  => 0,
               bgc                   => "default",
               crop                  => 0,
               dynamic_vegetation    => 0,
               envxml_dir            => ".",
               vichydro              => 0,
               maxpft                => "default",
             );

  GetOptions(
             "cimeroot=s"                => \$opts{'cimeroot'},
             "driver=s"                  => \$opts{'driver'},
             "clm_demand=s"              => \$opts{'clm_demand'},
             "co2_ppmv=f"                => \$opts{'co2_ppmv'},
             "co2_type=s"                => \$opts{'co2_type'},
             "config=s"                  => \$opts{'config'},
             "configuration=s"           => \$opts{'configuration'},
             "csmdata=s"                 => \$opts{'csmdata'},
             "clm_usr_name=s"            => \$opts{'clm_usr_name'},
             "envxml_dir=s"              => \$opts{'envxml_dir'},
             "drydep!"                   => \$opts{'drydep'},
             "lilac!"                    => \$opts{'lilac'},
             "fire_emis!"                => \$opts{'fire_emis'},
             "ignore_warnings!"          => \$opts{'ignore_warnings'},
             "chk_res!"                  => \$opts{'chk_res'},
             "note!"                     => \$opts{'note'},
             "megan!"                    => \$opts{'megan'},
             "glc_nec=i"                 => \$opts{'glc_nec'},
             "glc_use_antarctica!"       => \$opts{'glc_use_antarctica'},
             "light_res=s"               => \$opts{'light_res'},
             "d:s"                       => \$opts{'dir'},
             "h|help"                    => \$opts{'help'},
             "ignore_ic_date"            => \$opts{'ignore_ic_date'},
             "ignore_ic_year"            => \$opts{'ignore_ic_year'},
             "infile=s"                  => \$opts{'infile'},
             "lnd_frac=s"                => \$opts{'lnd_frac'},
             "lnd_tuning_mode=s"         => \$opts{'lnd_tuning_mode'},
             "inputdata=s"               => \$opts{'inputdata'},
             "mask=s"                    => \$opts{'mask'},
             "namelist=s"                => \$opts{'namelist'},
             "res=s"                     => \$opts{'res'},
             "ssp_rcp=s"                 => \$opts{'ssp_rcp'},
             "s|silent"                  => \$opts{'silent'},
             "sim_year=s"                => \$opts{'sim_year'},
             "structure=s"               => \$opts{'structure'},
             "output_reals=s"            => \$opts{'output_reals_filename'},
             "clm_accelerated_spinup=s"  => \$opts{'clm_accelerated_spinup'},
             "clm_start_type=s"          => \$opts{'clm_start_type'},
             "test"                      => \$opts{'test'},
             "use_case=s"                => \$opts{'use_case'},
             "bgc=s"                     => \$opts{'bgc'},
             "crop!"                     => \$opts{'crop'},
             "dynamic_vegetation"        => \$opts{'dynamic_vegetation'},
             "vichydro"                  => \$opts{'vichydro'},
             "maxpft=i"                  => \$opts{'maxpft'},
             "v|verbose"                 => \$opts{'verbose'},
             "version"                   => \$opts{'version'},
            )  or usage();

  # Give usage message.
  usage() if $opts{'help'};

  # Check for unparsed arguments
  if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
  }
  return %opts;
}

#-------------------------------------------------------------------------------

sub check_for_perl_utils {

  my $cfgdir = shift;
  my $opts_ref = shift;

  # Determine CIME root directory and perl5lib root directory
  my $cimeroot = $opts_ref->{'cimeroot'};
  if ( ! defined($cimeroot) ) {
    $cimeroot = "$cfgdir/../cime";
    if (      -d $cimeroot ) {
    } elsif ( -d "$cfgdir/../../../cime" ) {
      $cimeroot = "$cfgdir/../../../cime";
    } else {
      die <<"EOF";
** Cannot find the root of the cime directory  enter it using the -cimeroot option
   Did you run the checkout_externals scripts?
EOF
    }
  }

  my $perl5lib_dir = "$cimeroot/utils/perl5lib";

  #-----------------------------------------------------------------------------
  # Add $perl5lib_dir to the list of paths that Perl searches for modules
  my @dirs = ( $ProgDir, $cfgdir, "$perl5lib_dir");
  unshift @INC, @dirs;

  require config_files::clm_phys_vers;
  require namelist_files::LogMessages;

  my $locallog = namelist_files::LogMessages->new( $ProgName, $opts_ref );
  # The XML::Lite module is required to parse the XML files.
  (-f "$perl5lib_dir/XML/Lite.pm")  or
      $locallog->fatal_error("Cannot find perl module \"XML/Lite.pm\" in directory\n" .
                "\"$perl5lib_dir\"");

  # The Build::Config module provides utilities to access the configuration information
  # in the config_cache.xml file
  (-f "$perl5lib_dir/Build/Config.pm")  or
      $locallog->fatal_error("Cannot find perl module \"Build/Config.pm\" in directory\n" .
                "\"$perl5lib_dir\"");

  # The Build::NamelistDefinition module provides utilities to validate that the output
  # namelists are consistent with the namelist definition file
  (-f "$perl5lib_dir/Build/NamelistDefinition.pm")  or
      $locallog->fatal_error("Cannot find perl module \"Build/NamelistDefinition.pm\" in directory\n" .
		  "\"$perl5lib_dir\"");

  # The Build::NamelistDefaults module provides a utility to obtain default values of namelist
  # variables based on finding a best fit with the attributes specified in the defaults file.
  (-f "$perl5lib_dir/Build/NamelistDefaults.pm")  or
      $locallog->fatal_error("Cannot find perl module \"Build/NamelistDefaults.pm\" in directory\n" .
		  "\"$perl5lib_dir\"");

  # The Build::Namelist module provides utilities to parse input namelists, to query and modify
  # namelists, and to write output namelists.
  (-f "$perl5lib_dir/Build/Namelist.pm")  or
      $locallog->fatal_error("Cannot find perl module \"Build/Namelist.pm\" in directory\n" .
		  "\"$perl5lib_dir\"");


  # required cesm perl modules
  require XML::Lite;
  require Build::Config;
  require Build::NamelistDefinition;
  require Build::NamelistDefaults;
  require Build::Namelist;
  require Config::SetupTools;
}

#-------------------------------------------------------------------------------

sub read_configure_definition {
  # Read the configure definition and specific config_cache file for this case
  # configure are the build-time settings for CLM
  my ($cfgdir, $opts) = @_;

  $log->verbose_message("Setting CLM configuration script directory to $cfgdir");

  # Create a configuration object from the default config_definition file
  my $configfile;
  if ( -f $opts->{'config'} ) {
    $configfile = $opts->{'config'};
  } else {
    $configfile = "$cfgdir/config_files/config_definition_ctsm.xml";
  }

  # Check that configuration cache file exists.
  $log->verbose_message("Using CLM configuration cache file $opts->{'config'}");
  if ( $configfile ne $opts->{'config'} ) {
    $log->fatal_error("Cannot find configuration cache file: \"$opts->{'config'}\"");
  }

  my $cfg = Build::Config->new("$configfile");

  return $cfg;
}

#-----------------------------------------------------------------------------------------------

sub read_namelist_definition {
  my ($cfgdir, $opts, $nl_flags) = @_;

  # The namelist definition file contains entries for all namelist
  # variables that can be output by build-namelist.
  my @nl_definition_files = ( "$cfgdir/namelist_files/namelist_definition_drv.xml",
                              "$cfgdir/namelist_files/namelist_definition_drv_flds.xml",
                              "$cfgdir/namelist_files/namelist_definition_ctsm.xml" );
  foreach my $nl_defin_file  ( @nl_definition_files ) {
    (-f "$nl_defin_file")  or  $log->fatal_error("Cannot find namelist definition file \"$nl_defin_file\"");

    $log->verbose_message("Using namelist definition file $nl_defin_file");
  }

  # Create a namelist definition object.  This object provides a
  # method for verifying that the output namelist variables are in the
  # definition file, and are output in the correct namelist groups.
  my $definition = Build::NamelistDefinition->new( shift(@nl_definition_files) );
  foreach my $nl_defin_file ( @nl_definition_files ) {
    $definition->add( "$nl_defin_file" );
  }

  return $definition;
}

#-----------------------------------------------------------------------------------------------

sub read_envxml_case_files {
  # read the contents of the env*.xml files in the case directory
  my ($opts) = @_;

  my %envxml = ();
  if ( defined($opts->{'envxml_dir'}) ) {
      (-d $opts->{'envxml_dir'})  or  $log->fatal_error( "envxml_dir is not a directory" );
      my @files = bsd_glob( $opts->{'envxml_dir'}."/env_*xml" );
      ($#files >= 0)              or  $log->fatal_error( "there are no env_*xml files in the envxml_dir" );
      foreach my $file (@files) {
          $log->verbose_message( "Open env.xml file: $file" );
          my $xml = XML::Lite->new( "$file" );
          my @e   = $xml->elements_by_name('entry');
          while ( my $e = shift @e ) {
              my %a = $e->get_attributes();
              $envxml{$a{'id'}} = $a{'value'};
          }
      }
      foreach my $attr (keys %envxml) {
          if ( $envxml{$attr} =~ m/\$/ ) {
             $envxml{$attr} = SetupTools::expand_xml_var( $envxml{$attr}, \%envxml );
          }
      }
  } else {
      $log->fatal_error( "The -envxml_dir option was NOT given and it is a REQUIRED option" );
  }
  return( %envxml );
}

#-----------------------------------------------------------------------------------------------

sub read_namelist_defaults {
  my ($cfgdir, $opts, $nl_flags, $cfg) = @_;

  # The namelist defaults file contains default values for all required namelist variables.
  my @nl_defaults_files = ( "$cfgdir/namelist_files/namelist_defaults_overall.xml",
                            "$cfgdir/namelist_files/namelist_defaults_ctsm.xml",
                            "$cfgdir/namelist_files/namelist_defaults_drv.xml",
                            "$cfgdir/namelist_files/namelist_defaults_fire_emis.xml",
                            "$cfgdir/namelist_files/namelist_defaults_drydep.xml" );

  # Add the location of the use case defaults files to the options hash
  $opts->{'use_case_dir'} = "$cfgdir/namelist_files/use_cases";

  if (defined $opts->{'use_case'}) {
    if ( $opts->{'use_case'} ne "list" ) {
      unshift( @nl_defaults_files, "$opts->{'use_case_dir'}/$opts->{'use_case'}.xml" );
    }
  }

  foreach my $nl_defaults_file ( @nl_defaults_files ) {
    (-f "$nl_defaults_file")  or  $log->fatal_error("Cannot find namelist defaults file \"$nl_defaults_file\"");

    $log->verbose_message("Using namelist defaults file $nl_defaults_file");
  }

  # Create a namelist defaults object.  This object provides default
  # values for variables contained in the input defaults file.  The
  # configuration object provides attribute values that are relevent
  # for the CLM executable for which the namelist is being produced.
  my $defaults = Build::NamelistDefaults->new( shift( @nl_defaults_files ), $cfg);
  foreach my $nl_defaults_file ( @nl_defaults_files ) {
    $defaults->add( "$nl_defaults_file" );
  }
  return $defaults;
}

#-------------------------------------------------------------------------------

sub check_cesm_inputdata {
  # Check that the CESM inputdata root directory has been specified.  This must be
  # a local or nfs mounted directory.

  my ($opts, $nl_flags) = @_;

  $nl_flags->{'inputdata_rootdir'} = undef;
  if (defined($opts->{'csmdata'})) {
    $nl_flags->{'inputdata_rootdir'} = $opts->{'csmdata'};
  }
  elsif (defined $ENV{'CSMDATA'}) {
    $nl_flags->{'inputdata_rootdir'} = $ENV{'CSMDATA'};
  }
  else {
    $log->fatal_error("CESM inputdata root directory must be specified by either -csmdata\n" .
                "argument or by the CSMDATA environment variable.");
  }
  if ( ! defined($ENV{'DIN_LOC_ROOT'}) ) {
    $ENV{'DIN_LOC_ROOT'} = $nl_flags->{'inputdata_rootdir'};
  }

  if ($opts->{'test'}) {
    (-d $nl_flags->{'inputdata_rootdir'})  or  $log->fatal_error("CESM inputdata root is not a directory: \"$nl_flags->{'inputdata_rootdir'}\"");
  }

  $log->verbose_message("CESM inputdata root directory: $nl_flags->{'inputdata_rootdir'}");
}

#-------------------------------------------------------------------------------

sub process_namelist_user_input {
  # Process the user input in general by order of precedence.  At each point
  # we'll only add new values to the namelist and not overwrite
  # previously specified specified values which have higher
  # precedence. The one exception to this rule are the specifc command-line
  # options which are done last as if the user contradicts these settings
  # CLM build-namelist will abort with an error.
  #
  # 1. values set on the command-line using the -namelist option,
  #         (i.e. CLM_NAMELIST_OPTS env_run variable)
  # 2. values read from the file(s) specified by -infile,
  #         (i.e.  user_nl_clm files)
  # After the above are done the command line options are processed and they
  # are made sure the user hasn't contradicted any of their settings with
  # anything above. Because of this they are condsidered to have the highest
  # precedence.
  # 0. namelist values set by specific command-line options, like, -d, -sim_year
  #         (i.e.  CLM_BLDNML_OPTS env_run variable)
  # The results of these are needed for the final two user input
  # 3. datasets from the -clm_usr_name option,
  #         (i.e.  CLM_USRDAT_NAME env_run variable)
  # 4. values set from a use-case scenario, e.g., -use_case
  #         (i.e.  CLM_NML_USE_CASE env_run variable)
  #
  # Finally after all the above is done, the defaults are found from the
  # namelist defaults file (outside of this routine).
  #


  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg, $envxml_ref, $physv) = @_;

  # Get the inputs that will be coming from the user...
  process_namelist_commandline_namelist($opts, $definition, $nl, $envxml_ref);
  process_namelist_commandline_infile($opts, $definition, $nl, $envxml_ref);

  # Apply the commandline options and make sure the user didn't change it above
  process_namelist_commandline_options($opts, $nl_flags, $definition, $defaults, $nl, $envxml_ref, $physv);

  # The last two process command line arguments for usr_name and use_case
  # They require that process_namelist_commandline_options was called before this
  process_namelist_commandline_clm_usr_name($opts, $nl_flags, $definition, $defaults, $nl, $cfg, $envxml_ref);
  process_namelist_commandline_use_case($opts, $nl_flags, $definition, $defaults, $nl, $cfg, $envxml_ref);

  # Set the start_type by the command line setting for clm_start_type
  process_namelist_commandline_clm_start_type($opts, $nl_flags, $definition, $defaults, $nl);

}

#-------------------------------------------------------------------------------

sub process_namelist_commandline_options {
  # First process the commandline args that provide specific namelist values.
  #
  # First get the command-line specified overall values or their defaults
  # Obtain default values for the following build-namelist input arguments
  # : res, mask, ssp_rcp, sim_year, sim_year_range, and clm_accelerated_spinup.

  my ($opts, $nl_flags, $definition, $defaults, $nl, $envxml_ref, $physv) = @_;

  setup_cmdl_chk_res($opts, $defaults);
  setup_cmdl_resolution($opts, $nl_flags, $definition, $defaults, $envxml_ref);
  setup_cmdl_mask($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_configuration_and_structure($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_bgc($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_fire_light_res($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_spinup($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_crop($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_maxpft($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_glc_nec($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_ssp_rcp($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_simulation_year($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_dynamic_vegetation($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_fates_mode($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_vichydro($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_lnd_tuning($opts, $nl_flags, $definition, $defaults, $nl, $physv);
  setup_cmdl_run_type($opts, $nl_flags, $definition, $defaults, $nl);
  setup_cmdl_output_reals($opts, $nl_flags, $definition, $defaults, $nl);
}

#-------------------------------------------------------------------------------

sub setup_cmdl_chk_res {
  my ($opts, $defaults) = @_;

  my $var = "chk_res";
  if ( ! defined($opts->{$var}) ) {
    $opts->{$var} = $defaults->get_value($var);
  }
}

sub setup_cmdl_resolution {
  my ($opts, $nl_flags, $definition, $defaults, $envxml_ref) = @_;

  my $var = "res";
  my $val;

  if ( $opts->{$var} ne "default" ) {
    $val = $opts->{$var};
  } else {
    $val= $defaults->get_value($var);
  }

  $nl_flags->{'res'} = $val;
  $log->verbose_message("CLM atm resolution is $nl_flags->{'res'}");
  $opts->{$var} = $val;
  if ( $opts->{'chk_res'} ) {
    $val = &quote_string( $nl_flags->{'res'} );
    if (  ! $definition->is_valid_value( $var, $val ) ) {
      my @valid_values   = $definition->get_valid_values( $var );
      if ( $nl_flags->{'res'} ne "CLM_USRDAT" ) {
        $log->fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values");
      }
    }
  }
  if ( $nl_flags->{'res'} eq "CLM_USRDAT" ) {
    if ( ! defined($opts->{'clm_usr_name'}) ) {
        $log->fatal_error("Resolution is CLM_USRDAT, but --clm_usr_name option is NOT set, and it is required for CLM_USRDAT resolutions");
    }
  }
  #
  # For NEON sites
  #
  $nl_flags->{'neon'} = ".false.";
  $nl_flags->{'neonsite'} = "";
  if ( $nl_flags->{'res'} eq "CLM_USRDAT" ) {
    if ( $opts->{'clm_usr_name'} eq "NEON" ) {
       $nl_flags->{'neon'} = ".true.";
       $nl_flags->{'neonsite'} = $envxml_ref->{'NEONSITE'};
       $log->verbose_message( "This is a NEON site with NEONSITE = " . $nl_flags->{'neonsite'} );
    }
  }
  if ( ! &value_is_true( $nl_flags->{'neon'} ) ) {
    $log->verbose_message( "This is NOT a NEON site" );
  }

}

#-------------------------------------------------------------------------------

sub setup_cmdl_mask {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $var = "mask";
  my $val;

  if ( $opts->{$var} ne "default" ) {
    $val = $opts->{$var};
  } else {
    my %tmp = ( 'hgrid'=>$nl_flags->{'res'} );
    $val = $defaults->get_value($var, \%tmp );
  }

  $nl_flags->{'mask'} = $val;
  $opts->{'mask'} = $nl_flags->{'mask'};
  if ( $opts->{'chk_res'} ) {
    $val = &quote_string( $val );
    my $group = $definition->get_group_name($var);
    $nl->set_variable_value($group, $var, $val);
    if (  ! $definition->is_valid_value( $var, $val ) ) {
      my @valid_values   = $definition->get_valid_values( $var );
      $log->fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values");
    }
  }
  $log->verbose_message("CLM land mask is $nl_flags->{'mask'}");
}

#-------------------------------------------------------------------------------
sub setup_cmdl_fates_mode {
  #
  # call this at least after crop check is called
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "bgc_mode";

  if ( $nl_flags->{'crop'} eq "on" ) {
    if ( $nl_flags->{$var} eq "fates" ) {
       # FATES should not be used with crop
       $log->fatal_error("** Cannot turn fates mode on with crop." );
    }
  } elsif ($nl_flags->{"bgc_mode"} eq "fates" && ! &value_is_true($nl_flags->{"use_fates"}) ) {
    $log->fatal_error("DEV_ERROR: internal logic error: bgc_mode = fates and use_fates = false.");

  } else {

    $var = "use_fates";
    if ( &value_is_true($nl_flags->{$var}) ) {
      # This section is a place-holder to test for modules that are not allowed with FATES
      # the defaults which are set in the logic section of the namelist builder will
      # automatically set these correctly (well that is the assumption), but here we
      # want to set a catch to fail and warn users if they explicitly set incompatible user namelist
      # options

       my $var = "use_crop";
       $val = $nl_flags->{$var};
       if ( defined($nl->get_value($var))  ) {
          if ( &value_is_true($nl->get_value($var)) ) {
             $log->fatal_error("$var was set to .true., which is incompatible when -bgc fates option is used.");
          }
       }

    } else {
       # dis-allow fates specific namelist items with non-fates runs
       my @list  = (  "fates_spitfire_mode", "use_fates_planthydro", "use_fates_ed_st3", "use_fates_ed_prescribed_phys",
                      "use_fates_cohort_age_tracking","use_fates_inventory_init","use_fates_fixed_biogeog",
                      "use_fates_nocomp","use_fates_sp","fates_inventory_ctrl_filename","use_fates_logging",
                      "fates_parteh_mode","use_fates_tree_damage","fates_history_dimlevel","fates_seeddisp_cadence",
		      "use_fates_luh","fluh_timeseries" );
       # dis-allow fates specific namelist items with non-fates runs
       foreach my $var ( @list ) {
          if ( defined($nl->get_value($var)) ) {
              $log->fatal_error("$var is being set, but can ONLY be set when -bgc fates option is used.\n");
          }
       }
    }
  }
}

#-------------------------------------------------------------------------------
sub setup_cmdl_configuration_and_structure {
   # Error-check and set the 'configuration' and 'structure' namelist flags

   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   my $val;

   my $var = "configuration";
   $val = $opts->{$var};
   if (defined($val) && ($val eq "clm" || $val eq "nwp")) {
      $nl_flags->{$var} = $val;
   } else {
      $log->fatal_error("$var has a value (".$val.") that is NOT valid. Valid values are: clm, nwp.");
   }

   $var = "structure";
   $val = $opts->{$var};
   if (defined($val) && ($val eq "standard" || $val eq "fast")) {
      $nl_flags->{$var} = $val;
   } else {
      $log->fatal_error("$var has a value (".$val.") that is NOT valid. Valid values are: standard, fast.");
   }
}

#-------------------------------------------------------------------------------
sub setup_cmdl_bgc {
  # BGC - alias for group of biogeochemistry related use_XXX namelists

  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "bgc";

  $val = $opts->{$var};
  $nl_flags->{'bgc_mode'} = $val;

  my $var = "bgc_mode";
  if ( $nl_flags->{$var} eq "default" ) {
     $nl_flags->{$var} = $defaults->get_value($var);
  }
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, quote_string( $nl_flags->{$var} ) );
  if (  ! $definition->is_valid_value( $var, quote_string( $nl_flags->{$var}) ) ) {
     my @valid_values   = $definition->get_valid_values( $var );
     $log->fatal_error("$var has a value (".$nl_flags->{$var}.") that is NOT valid. Valid values are: @valid_values");
  }
  $log->verbose_message("Using $nl_flags->{$var} for bgc.");

  # now set the actual name list variables based on the bgc alias
  if ($nl_flags->{$var} eq "bgc" ) {
     $nl_flags->{'use_cn'} = ".true.";
     $nl_flags->{'use_fates'} = ".false.";
  } elsif ($nl_flags->{$var} eq "fates" ) {
     $nl_flags->{'use_cn'} = ".false.";
     $nl_flags->{'use_fates'} = ".true.";
  } else {
     $nl_flags->{'use_cn'} = ".false.";
     $nl_flags->{'use_fates'} = ".false.";
  }
  if ( defined($nl->get_value("use_cn")) && ($nl_flags->{'use_cn'} ne $nl->get_value("use_cn")) ) {
     $log->fatal_error("The namelist variable use_cn is inconsistent with the -bgc option");
  }
  if ( defined($nl->get_value("use_fates")) && ($nl_flags->{'use_fates'} ne $nl->get_value("use_fates")) ) {
     $log->fatal_error("The namelist variable use_fates is inconsistent with the -bgc option");
  }

  # Now set use_cn and use_fates
  foreach $var ( "use_cn", "use_fates" ) {
     $val = $nl_flags->{$var};
     $group = $definition->get_group_name($var);
     $nl->set_variable_value($group, $var, $val);
     if (  ! $definition->is_valid_value( $var, $val ) ) {
        my @valid_values   = $definition->get_valid_values( $var );
        $log->fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values");
     }
  }
  #
  # Set FATES-SP mode
  #
  if ( &value_is_true( $nl_flags->{'use_fates'} ) ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_fates_sp', 'use_fates'=>$nl_flags->{'use_fates'} );
     if ( &value_is_true($nl->get_value('use_fates_sp')) ) {
        $nl_flags->{'use_fates_sp'} = ".true.";
     } else {
        $nl_flags->{'use_fates_sp'} = ".false.";
     }
  } else {
     $nl_flags->{'use_fates_sp'} = ".false.";
  }
  #
  # Determine Soil decomposition method
  #
  my $var = "soil_decomp_method";
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
              'phys'=>$nl_flags->{'phys'}, 'use_cn'=>$nl_flags->{'use_cn'}, 'use_fates'=>$nl_flags->{'use_fates'},
              'use_fates_sp'=>$nl_flags->{'use_fates_sp'} );
  my $soil_decomp_method = remove_leading_and_trailing_quotes( $nl->get_value( $var ) );

  if (      &value_is_true($nl_flags->{'use_cn'}) ) {
     if ( $soil_decomp_method eq "None" ) {
        $log->fatal_error("$var must NOT be None if use_cn is on");
     }
  } elsif ( &value_is_true($nl_flags->{'use_fates'}) && (not &value_is_true($nl_flags->{'use_fates_sp'}))  )  {
     if ( $soil_decomp_method eq "None" ) {
        $log->fatal_error("$var must NOT be None if use_fates is on and use_fates_sp is not TRUE");
     }
  } elsif ( $soil_decomp_method ne "None" ) {
     $log->fatal_error("$var must be None if use_cn and use_fates are off");
  }
  #
  # Soil decomposition control variables, methane and Nitrification-Denitrification
  #
  my @list  = (  "use_lch4", "use_nitrif_denitrif" );
  my %settings = ( 'bgc_mode'=>$nl_flags->{'bgc_mode'} );
  foreach my $var ( @list ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                 'phys'=>$nl_flags->{'phys'}, 'soil_decomp_method'=>$soil_decomp_method );
     $nl_flags->{$var} = $nl->get_value($var);
  }
  if ( $soil_decomp_method eq "None" ) {
     foreach my $var ( @list ) {
        if ( &value_is_true($nl_flags->{$var}) ) {
           $log->fatal_error("When soil_decomp_method is None $var can NOT be TRUE");
        }
     }
  } else {
     # nitrif_denitrif can only be .false. if fates is on
     if ( (! &value_is_true($nl_flags->{'use_fates'})) && &value_is_true($nl_flags->{'use_cn'}) ) {
        $var = "use_nitrif_denitrif";
        if ( ! &value_is_true($nl_flags->{$var}) ) {
           $log->warning("$var normally use_nitrif_denitrif should only be FALSE if FATES is on, it has NOT been validated for being off for BGC mode" );
        }
     }
     # if MIMICS is on and use_fates = .true. then use_lch4 must = .true.
     if ( (! &value_is_true($nl_flags->{'use_lch4'})) && &value_is_true($nl_flags->{'use_fates'}) ) {
        if ( $soil_decomp_method eq "MIMICSWieder2015" ) {
           $log->warning("If MIMICS is on and use_fates = .true. then use_lch4 must be .true. and currently it's not" );
        }
     }
  }
  #
  # Set FUN for BGC
  #
  my $var = "use_fun";
  if ( ! defined($nl->get_value($var)) ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                 'phys'=>$nl_flags->{'phys'}, 'use_cn'=>$nl_flags->{'use_cn'},
                 'use_nitrif_denitrif'=>$nl_flags->{'use_nitrif_denitrif'} );
  }
  if ( (! &value_is_true($nl_flags->{'use_cn'}) ) && &value_is_true($nl->get_value('use_fun')) ) {
     $log->fatal_error("When FUN is on, use_cn MUST also be on!");
  }
  if ( (! &value_is_true($nl_flags->{'use_nitrif_denitrif'}) ) && &value_is_true($nl->get_value('use_fun')) ) {
     $log->fatal_error("When FUN is on, use_nitrif_denitrif MUST also be on!");
  }
} # end bgc


#-------------------------------------------------------------------------------
sub setup_cmdl_fire_light_res {
  # light_res - alias for lightning resolution

  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $var = "light_res";
  my $val = $opts->{$var};
  if ( &value_is_true($nl->get_value('use_cn')) ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fire_method');
  }
  my $fire_method = remove_leading_and_trailing_quotes( $nl->get_value('fire_method') );
  if ( $val eq "default" ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                 'phys'=>$nl_flags->{'phys'}, 'use_cn'=>$nl_flags->{'use_cn'},
                 'neon'=>$nl_flags->{'neon'},
                 'fates_spitfire_mode'=>$nl->get_value('fates_spitfire_mode'),
                 'use_fates'=>$nl_flags->{'use_fates'}, fire_method=>$fire_method );
     $val              = remove_leading_and_trailing_quotes( $nl->get_value($var) );
     $nl_flags->{$var} = $val;
  } else {
     if ( defined($fire_method) && $val ne "none" ) {
        if ( $fire_method eq "nofire" ) {
           $log->fatal_error("-$var option used with fire_method='nofire'. -$var can ONLY be used without the nofire option");
        }
     }
     my $stream_fldfilename_lightng = remove_leading_and_trailing_quotes( $nl->get_value('stream_fldfilename_lightng') );
     if ( defined($stream_fldfilename_lightng) && $val ne "none" ) {
        $log->fatal_error("-$var option used while also explicitly setting stream_fldfilename_lightng filename which is a contradiction. Use one or the other not both.");
     }
     if ( ! &value_is_true($nl->get_value('use_cn')) ) {
        if ( &value_is_true($nl_flags->{'use_fates'}) ) {
           if ( $nl->get_value('fates_spitfire_mode') < 2) {
              if ( $val ne "none" ) {
                  $log->fatal_error("-$var option used when FATES is on, but fates_spitfire_mode does NOT use lightning data");
              }
           } else {
              if ( $val eq "none" ) {
                 $log->fatal_error("-$var option is set to none, but FATES is on and fates_spitfire_mode requires lightning data");
              }
           }
        } else {
           $log->fatal_error("-$var option used when FATES off and CN is NOT on. -$var can only be used when BGC is set to bgc or fates");
        }
     } else {
        if ( $val eq "none" and $fire_method ne "nofire" ) {
           $log->fatal_error("-$var option is set to none, but CN is on (with bgc: cn or bgc) which is a contradiction");
        }
     }
     $nl_flags->{$var} = $val;
  }
  # Check that NEON data is only used for NEON sites
  if ( $val eq "106x174" ) {
     if ( ! &value_is_true($nl_flags->{'neon'}) ) {
         if ( defined($opts->{'clm_usr_name'}) ) {
            $log->warning("The NEON lightning dataset does NOT cover the entire globe, make sure it covers the region for your grid");
         } else { 
            $log->fatal_error("The NEON lightning dataset can NOT be used for global grids or regions or points outside of its area as it does NOT cover the entire globe.");
         }
     }
  }
  # check for valid values...
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, quote_string($nl_flags->{$var}) );
  if (  ! $definition->is_valid_value( $var, $nl_flags->{$var}, 'noquotes'=>1 ) ) {
    my @valid_values   = $definition->get_valid_values( $var );
    $log->fatal_error("$var has a value (".$nl_flags->{$var}.") that is NOT valid. Valid values are: @valid_values");
  }
  $log->verbose_message("Using $nl_flags->{$var} for $var.");
  #
  # Set flag if cn-fires are on or not, only for BGC (not FATES)
  #
  $var = "cnfireson";
  my $fire_method = remove_leading_and_trailing_quotes( $nl->get_value('fire_method') );
  if ( defined($fire_method) && ! &value_is_true($nl_flags->{'use_cn'}) && ! &value_is_true($nl_flags->{'use_fates'}) ) {
     $log->fatal_error("fire_method is being set while use_cn and use_fates are both false.");
  }
  if ( defined($fire_method) && $fire_method eq "nofire" ) {
     $nl_flags->{$var} = ".false.";
# } elsif ( &value_is_true($nl->get_value('use_cn')) || $nl_flags->{'fates_spitfire_mode'} > 1 ) {
  } elsif ( &value_is_true($nl->get_value('use_cn')) || &value_is_true($nl->get_value('use_fates')) ) {
     $nl_flags->{$var} = ".true.";
  } else {
     $nl_flags->{$var} = ".false.";
  }
}

#-------------------------------------------------------------------------------

sub setup_cmdl_crop {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  $nl_flags->{'use_crop'} = ".false.";
  my $val;
  my $var = "crop";
  $val = $opts->{$var};
  $nl_flags->{'crop'} = $val;
  if ( $nl_flags->{'crop'} eq 1 ) {
     $nl_flags->{'use_crop'} = ".true.";
  }
  if ( defined($nl->get_value("use_crop")) && ($nl_flags->{'use_crop'} ne $nl->get_value("use_crop")) ) {
     $log->fatal_error("Namelist item use_crop contradicts the command-line option -crop, use the command line option");
  }
  if ( ($nl_flags->{'crop'} eq 1 ) && ($nl_flags->{'bgc_mode'} eq "sp") ) {
     $log->fatal_error("** Cannot turn crop mode on mode bgc=sp\n" .
                       "**\n" .
                       "** Set the bgc mode to 'bgc' by the following means from highest to lowest precedence:\n" .
                       "** * by the command-line options -bgc bgc\n" .
                       "** * by a default configuration file, specified by -defaults");
  }

  $var = "use_crop";
  $val = ".false.";
  if ($nl_flags->{'crop'} eq 1) {
     $val = ".true.";
  }
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, $val);
  if (  ! $definition->is_valid_value( $var, $val ) ) {
     my @valid_values   = $definition->get_valid_values( $var );
     $log->fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values");
  }
}

#-------------------------------------------------------------------------------

sub setup_cmdl_maxpft {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "maxpft";
  my %maxpatchpft;
  $maxpatchpft{'.true.'}   = 79;
  $maxpatchpft{'.false.'} = 17;
  if ( $opts->{$var} ne "default") {
     $val = $opts->{$var};
  } else {
     $val = $maxpatchpft{$nl_flags->{'use_crop'}};
  }
  $nl_flags->{'maxpft'} = $val;

}

#-------------------------------------------------------------------------------

sub setup_cmdl_glc_nec {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "glc_nec";

  if ( $opts->{$var} ne "default" ) {
    $val = $opts->{$var};
  } else {
    $val = $defaults->get_value($var);
  }

  $nl_flags->{'glc_nec'} = $val;
  $opts->{'glc_nec'} = $val;
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, $val);
  if (  ! $definition->is_valid_value( $var, $val ) ) {
    my @valid_values   = $definition->get_valid_values( $var );
    $log->fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values");
  }
  $log->verbose_message("Glacier number of elevation classes is $val");
}

#-------------------------------------------------------------------------------

sub setup_cmdl_ssp_rcp {
  # shared socioeconmic pathway and representative concentration pathway combination
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "ssp_rcp";
  if ( $opts->{$var} ne "default" ) {
    $val = $opts->{$var};
  } else {
    $val = remove_leading_and_trailing_quotes( $defaults->get_value($var) );
  }
  $nl_flags->{'ssp_rcp'} = $val;
  $opts->{'ssp_rcp'} = $nl_flags->{'ssp_rcp'};
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, quote_string($val) );
  if (  ! $definition->is_valid_value( $var, $val, 'noquotes'=>1 ) ) {
    my @valid_values   = $definition->get_valid_values( $var );
    $log->fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values");
  }
  $log->verbose_message("CLM future scenario SSP-RCP combination is $nl_flags->{'ssp_rcp'}");
}

#-------------------------------------------------------------------------------

sub setup_cmdl_spinup {
  # BGC spinup mode controlled from "clm_accelerated_spinup" in build-namelist
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var;
  $nl_flags->{'spinup'} = undef;
  $var = "clm_accelerated_spinup";
  if ( $opts->{$var} ne "default" ) {
    $val = $opts->{$var};
  } else {
    $val = $defaults->get_value($var);
  }
  $nl_flags->{$var} = $val;
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, quote_string($val) );
  if (  ! $definition->is_valid_value( $var, $val , 'noquotes' => 1) ) {
    my @valid_values   = $definition->get_valid_values( $var );
    $log->fatal_error("$var has an invalid value ($val). Valid values are: @valid_values");
  }
  $log->verbose_message("CLM accelerated spinup mode is $val");
  if ( &value_is_true($nl_flags->{'use_cn'}) ) {
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition,
                $defaults, $nl, "spinup_state", clm_accelerated_spinup=>$nl_flags->{$var},
                use_cn=>$nl_flags->{'use_cn'}, use_fates=>$nl_flags->{'use_fates'} );
    if ( $nl->get_value("spinup_state") ne 0 ) {
       $nl_flags->{'bgc_spinup'} = "on";
       if ( $nl_flags->{'clm_accelerated_spinup'} eq "off" ) {
          $log->fatal_error("spinup_state is accelerated, but clm_accelerated_spinup is off, change one or the other");
       }
    } else {
       $nl_flags->{'bgc_spinup'} = "off";
       $val = $defaults->get_value($var);
    }
    # For AD spinup mode by default reseed dead plants
    if ( $nl_flags->{$var} ne "off" ) {
        add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition,
                    $defaults, $nl, "reseed_dead_plants", clm_accelerated_spinup=>$nl_flags->{$var},
                    use_cn=>$nl_flags->{'use_cn'} );
    }
  } else {
    if ( defined($nl->get_value("spinup_state")) ) {
       $log->fatal_error("spinup_state is accelerated (=1 or 2) which is for a BGC mode of CN or BGC," .
                         " but the BGC mode is Satellite Phenology, change one or the other");
    }
  }
  $nl_flags->{$var} = $val;
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, quote_string($val) );
  if (  ! $definition->is_valid_value( $var, $val , 'noquotes' => 1) ) {
     my @valid_values   = $definition->get_valid_values( $var );
     $log->fatal_error("$var has an invalid value ($val). Valid values are: @valid_values");
  }
  if ( $nl_flags->{'bgc_spinup'} eq "on" && (not &value_is_true( $nl_flags->{'use_cn'} ))  && (not &value_is_true($nl_flags->{'use_fates'})) ) {
     $log->fatal_error("$var can not be '$nl_flags->{'bgc_spinup'}' if neither CN nor FATES is turned on (use_cn=$nl_flags->{'use_cn'}, use_fates=$nl_flags->{'use_fates'}).");
  }
  if ( $nl->get_value("spinup_state") eq 0 && $nl_flags->{'bgc_spinup'} eq "on" ) {
     $log->fatal_error("Namelist spinup_state contradicts the command line option bgc_spinup" );
  }
  if ( $nl->get_value("spinup_state") eq 1 && $nl_flags->{'bgc_spinup'} eq "off" ) {
     $log->fatal_error("Namelist spinup_state contradicts the command line option bgc_spinup" );
  }

  $val = $nl_flags->{'bgc_spinup'};
  $log->verbose_message("CLM CN bgc_spinup mode is $val");
}

#-------------------------------------------------------------------------------

sub setup_cmdl_simulation_year {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "sim_year";
  if ( $opts->{$var} ne "default" ) {
    $val = $opts->{$var};
  } else {
    $val = $defaults->get_value($var);
  }

  $nl_flags->{'sim_year_range'} = $defaults->get_value("sim_year_range");
  $nl_flags->{'sim_year'}       = &remove_leading_and_trailing_quotes($val);
  if ( $val =~ /([0-9]+)-([0-9]+)/ ) {
    $nl_flags->{'sim_year'}       = $1;
    $nl_flags->{'sim_year_range'} = $val;
  }
  $val = $nl_flags->{'sim_year'};
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, "'$val'" );
  if (  ! $definition->is_valid_value( $var, $val, 'noquotes'=>1 ) ) {
    my @valid_values   = $definition->get_valid_values( $var );
    $log->fatal_error("$var of $val is NOT valid. Valid values are: @valid_values");
  }
  $nl->set_variable_value($group, $var, "'$val'" );
  $log->verbose_message("CLM sim_year is $nl_flags->{'sim_year'}");

  $var = "sim_year_range";
  $val = $nl_flags->{'sim_year_range'};
  if ( $val ne "constant" ) {
    $opts->{$var}   = $val;
    $group = $definition->get_group_name($var);
    $nl->set_variable_value($group, $var, $val );
    if (  ! $definition->is_valid_value( $var, $val, 'noquotes'=>1 ) ) {
      my @valid_values   = $definition->get_valid_values( $var );
      $log->fatal_error("$var of $val is NOT valid. Valid values are: @valid_values");
    }
    $val = "'".$defaults->get_value($var)."'";
    $nl->set_variable_value($group, $var, $val );
    $log->verbose_message("CLM sim_year_range is $nl_flags->{'sim_year_range'}");
  }
}

#-------------------------------------------------------------------------------

sub setup_cmdl_run_type {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;
  # Set the clm_start_type and the st_year, start year
  # This MUST be done after lnd_tuning_mode is set

  my $val;
  my $var = "clm_start_type";
  my $ic_date = $nl->get_value('start_ymd');
  my $st_year;
  if ( defined($ic_date) ) {
     $st_year = int( $ic_date / 10000);
  } else {
    $st_year = $nl_flags->{'sim_year'};
    $ic_date = $st_year *10000 + 101;
    my $date = 'start_ymd';
    my $group = $definition->get_group_name($date);
    $nl->set_variable_value($group, $date, $ic_date );
  }
  my $set = undef;
  if (defined $opts->{$var}) {
    if ($opts->{$var} ne "default" ) {
      $set = 1;
      my $group = $definition->get_group_name($var);
      $nl->set_variable_value($group, $var, quote_string( $opts->{$var} ) );
    }
  }
  if ( ! defined $set ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                 'use_cndv'=>$nl_flags->{'use_cndv'}, 'use_fates'=>$nl_flags->{'use_fates'},
                 'sim_year'=>$st_year, 'sim_year_range'=>$nl_flags->{'sim_year_range'},
                 'bgc_spinup'=>$nl_flags->{'bgc_spinup'}, 'lnd_tuning_mode'=>$nl_flags->{'lnd_tuning_mode'} );
  }
  $nl_flags->{'clm_start_type'} = $nl->get_value($var);
  $nl_flags->{'st_year'}        = $st_year;
}

#-------------------------------------------------------------------------------

sub setup_cmdl_dynamic_vegetation {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "dynamic_vegetation";
  $val = $opts->{$var};
  $nl_flags->{'dynamic_vegetation'} = $val;
  if ( ($nl_flags->{'dynamic_vegetation'} eq 1 ) && ($nl_flags->{'bgc_mode'} eq "sp") ) {
     $log->fatal_error("** Cannot turn dynamic_vegetation mode on with bgc=sp.\n" .
                       "**\n" .
                       "** Set the bgc mode to 'bgc' by the following means from highest to lowest precedence:" .
                       "** * by the command-line options -bgc bgc\n");
  }

  $var = "use_cndv";
  $nl_flags->{$var} = ".false.";
  if ($nl_flags->{'dynamic_vegetation'} eq 1) {
     $val = ".true.";
     $nl_flags->{$var} = $val;
  }
  if ( defined($nl->get_value($var)) && $nl->get_value($var) ne $val ) {
     $log->fatal_error("$var is inconsistent with the commandline setting of -dynamic_vegetation");
  }
  if ( &value_is_true($nl_flags->{$var}) ) {
     my $group = $definition->get_group_name($var);
     $nl->set_variable_value($group, $var, $val);
     if (  ! $definition->is_valid_value( $var, $val ) ) {
        my @valid_values   = $definition->get_valid_values( $var );
        $log->fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values");
     }
     $log->warning("The use_cndv=T option is deprecated. We do NOT recommend using it." .
                   " It's known to have issues and it's not calibrated.");
  }
}
#-------------------------------------------------------------------------------

sub setup_cmdl_output_reals {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $var = "output_reals_filename";
  my $file = $opts->{$var};
  if ( defined($file) ) {
     # Make sure can open file and if not die with an error
     my $fh = IO::File->new($file, '>') or $log->fatal_error("can't create real parameter filename: $file");
     $fh->close();
  }
}

#-------------------------------------------------------------------------------

sub setup_cmdl_vichydro {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $val;
  my $var = "vichydro";
  $val = $opts->{$var};
  $nl_flags->{'vichydro'} = $val;
  if ($nl_flags->{'vichydro'} eq 1) {
     $log->verbose_message("Using VIC hydrology for runoff calculations.");
  }

  $var = "use_vichydro";
  $val = $nl->get_value($var);
  my $set = undef;
  if ($nl_flags->{'vichydro'} eq 1) {
     my $group = $definition->get_group_name($var);
     $set = ".true.";
     if ( defined($val) && $set ne $val ) {
        $log->fatal_error("$var contradicts the command-line -vichydro option" );
     }
     $nl->set_variable_value($group, $var, $set);
     if ( ! $definition->is_valid_value($var, $val) ) {
        my @valid_values   = $definition->get_valid_values( $var );
        $log->fatal_error("$var has a value ($val) that is NOT valid. Valid values are: @valid_values");
     }
  } else {
     $set = ".false.";
  }
  $nl_flags->{$var} = $set;
}


#-------------------------------------------------------------------------------

sub process_namelist_commandline_namelist {
  # Process the commandline '-namelist' arg.
  my ($opts, $definition, $nl, $envxml_ref) = @_;

  if (defined $opts->{'namelist'}) {
    # Parse commandline namelist
    my $nl_arg = Build::Namelist->new($opts->{'namelist'});

    # Validate input namelist -- trap exceptions
    my $nl_arg_valid;
    eval { $nl_arg_valid = $definition->validate($nl_arg); };
    if ($@) {
      $log->fatal_error("Invalid namelist variable in commandline arg '-namelist'.\n $@");
    }
    # Go through all variables and expand any XML env settings in them
    expand_xml_variables_in_namelist( $nl_arg_valid, $envxml_ref );

    # Merge input values into namelist.  Previously specified values have higher precedence
    # and are not overwritten.
    $nl->merge_nl($nl_arg_valid);
  }
}

#-------------------------------------------------------------------------------

sub process_namelist_commandline_infile {
  # Process the commandline '-infile' arg.
  my ($opts, $definition, $nl, $envxml_ref) = @_;

  if (defined $opts->{'infile'}) {
    my @infiles = split( /,/, $opts->{'infile'} );
    foreach my $infile ( @infiles ) {
      # Make sure a valid file was found
      if (    -f "$infile" ) {
        # Otherwise abort as a valid file doesn't exist
      } else {
        $log->fatal_error("input namelist file does NOT exist $infile.\n $@");
      }
      # Parse namelist input from the next file
      my $nl_infile = Build::Namelist->new($infile);

      # Validate input namelist -- trap exceptions
      my $nl_infile_valid;
      eval { $nl_infile_valid = $definition->validate($nl_infile); };
      if ($@) {
        $log->fatal_error("Invalid namelist variable in '-infile' $infile.\n $@");
      }
      # Go through all variables and expand any XML env settings in them
      expand_xml_variables_in_namelist( $nl_infile_valid, $envxml_ref );

      # Merge input values into namelist.  Previously specified values have higher precedence
      # and are not overwritten.
      $nl->merge_nl($nl_infile_valid);
    }
  }
}

#-------------------------------------------------------------------------------

sub process_namelist_commandline_clm_usr_name {
  # Process the -clm_usr_name argument
  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg, $envxml_ref) = @_;

  if (defined $opts->{'clm_usr_name'}) {
    # The user files definition is contained in an xml file with the same format as the defaults file.

    # The one difference is that variables are expanded.
    # Create a new NamelistDefaults object.
    my $nl_defaults_file = "$nl_flags->{'cfgdir'}/namelist_files/namelist_defaults_usr_files.xml";
    my $uf_defaults = Build::NamelistDefaults->new("$nl_defaults_file", $cfg );
    # Loop over the variables specified in the user files
    # Add each one to the namelist.
    my @vars = $uf_defaults->get_variable_names();
    my %settings;
    $settings{'mask'}           = $nl_flags->{'mask'};
    $settings{'sim_year'}       = $nl_flags->{'sim_year'};
    $settings{'ssp_rcp'}        = $nl_flags->{'ssp_rcp'};
    $settings{'sim_year_range'} = $nl_flags->{'sim_year_range'};
    $settings{'clm_accelerated_spinup'} = $nl_flags->{'clm_accelerated_spinup'};
    $settings{'clm_usr_name'}   = $opts->{'clm_usr_name'};

    if ( $nl_flags->{'inputdata_rootdir'} eq "\$DIN_LOC_ROOT" ) {
      $settings{'csmdata'}     = $ENV{'DIN_LOC_ROOT'};
    } else {
      $settings{'csmdata'}     = $nl_flags->{'inputdata_rootdir'};
    }

    my $nvars = 0;
    my $nl_usrfile = Build::Namelist->new();
    foreach my $var (@vars) {
      my $val = $uf_defaults->get_usr_file($var, $definition, \%settings);

      if ($val) {
        $log->message("adding clm user file defaults for var $var with val $val");
        add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl_usrfile, $var, 'val'=>$val);
        $nvars++;
      }
    }
    # Go through all variables and expand any XML env settings in them
    expand_xml_variables_in_namelist( $nl_usrfile, $envxml_ref );
    # Merge input values into namelist.  Previously specified values have higher precedence
    # and are not overwritten.
    $nl->merge_nl($nl_usrfile);
  }
}

#-------------------------------------------------------------------------------

sub process_namelist_commandline_use_case {
  # Now process the -use_case arg.
  my ($opts, $nl_flags, $definition, $defaults, $nl, $cfg, $envxml_ref) = @_;

  if (defined $opts->{'use_case'}) {

    # The use case definition is contained in an xml file with the same format as the defaults file.
    # Create a new NamelistDefaults object.
    my $uc_defaults = Build::NamelistDefaults->new("$opts->{'use_case_dir'}/$opts->{'use_case'}.xml", $cfg);

    my %settings;
    $settings{'res'}            = $nl_flags->{'res'};
    $settings{'ssp_rcp'}        = $nl_flags->{'ssp_rcp'};
    $settings{'mask'}           = $nl_flags->{'mask'};
    $settings{'sim_year'}       = $nl_flags->{'sim_year'};
    $settings{'sim_year_range'} = $nl_flags->{'sim_year_range'};
    $settings{'phys'}           = $nl_flags->{'phys'};
    $settings{'lnd_tuning_mode'}= $nl_flags->{'lnd_tuning_mode'};
    $settings{'use_cn'}      = $nl_flags->{'use_cn'};
    $settings{'use_cndv'}    = $nl_flags->{'use_cndv'};
    $settings{'use_crop'}    = $nl_flags->{'use_crop'};
    $settings{'cnfireson'}   = $nl_flags->{'cnfireson'};

    # Loop over the variables specified in the use case.
    # Add each one to the namelist.
    my @vars = $uc_defaults->get_variable_names();
    my $nl_usecase = Build::Namelist->new();
    foreach my $var (@vars) {
      my $val = $uc_defaults->get_value($var, \%settings );

      if ( defined($val) ) {
        $log->verbose_message("CLM adding use_case $opts->{'use_case'} defaults for var '$var' with val '$val'");

        add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl_usecase, $var, 'val'=>$val);
      }
    }
    # Go through all variables and expand any XML env settings in them
    expand_xml_variables_in_namelist( $nl_usecase, $envxml_ref );

    # Merge input values into namelist.  Previously specified values have higher precedence
    # and are not overwritten.
    $nl->merge_nl($nl_usecase);
  }
}

#-------------------------------------------------------------------------------

sub process_namelist_commandline_clm_start_type {
  # Set the start_type according to the command line clm_start_type option

  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  # Run type for driver namelist - note that arb_ic implies that the run is startup
  my $var = "start_type";
  if ($nl_flags->{'clm_start_type'} eq "'cold'" || $nl_flags->{'clm_start_type'} eq "'arb_ic'") {
    # Add default is used here, but the value is explicitly set
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'val'=>'startup'   );
  } else {
    # Add default is used here, but the value is explicitly set
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'val'=>$nl_flags->{'clm_start_type'} );
  }
}

#-------------------------------------------------------------------------------

sub process_namelist_inline_logic {
  # Use the namelist default object to add default values for required
  # namelist variables that have not been previously set.
  my ($opts, $nl_flags, $definition, $defaults, $nl, $envxml_ref, $physv) = @_;


  ##############################
  # namelist group: clm_inparm #
  ##############################
  setup_logic_site_specific($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_lnd_frac($opts, $nl_flags, $definition, $defaults, $nl, $envxml_ref);
  setup_logic_co2_type($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_irrigate($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_start_type($opts, $nl_flags, $nl);
  setup_logic_decomp_performance($opts,  $nl_flags, $definition, $defaults, $nl);
  setup_logic_roughness_methods($opts, $nl_flags, $definition, $defaults, $nl, $physv);
  setup_logic_snicar_methods($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_snow($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_glacier($opts, $nl_flags, $definition, $defaults, $nl,  $envxml_ref);
  setup_logic_dynamic_plant_nitrogen_alloc($opts, $nl_flags, $definition, $defaults, $nl, $physv);
  setup_logic_luna($opts, $nl_flags, $definition, $defaults, $nl, $physv);
  setup_logic_hillslope($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_o3_veg_stress_method($opts, $nl_flags, $definition, $defaults, $nl,$physv);
  setup_logic_hydrstress($opts,  $nl_flags, $definition, $defaults, $nl);
  setup_logic_dynamic_roots($opts,  $nl_flags, $definition, $defaults, $nl, $physv);
  setup_logic_params_file($opts,  $nl_flags, $definition, $defaults, $nl);
  setup_logic_create_crop_landunit($opts,  $nl_flags, $definition, $defaults, $nl);
  setup_logic_subgrid($opts,  $nl_flags, $definition, $defaults, $nl);
  setup_logic_fertilizer($opts,  $nl_flags, $definition, $defaults, $nl);
  setup_logic_grainproduct($opts,  $nl_flags, $definition, $defaults, $nl, $physv);
  setup_logic_soilstate($opts,  $nl_flags, $definition, $defaults, $nl);
  setup_logic_demand($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_surface_dataset($opts,  $nl_flags, $definition, $defaults, $nl, $envxml_ref);
  setup_logic_dynamic_subgrid($opts,  $nl_flags, $definition, $defaults, $nl);
  if ( remove_leading_and_trailing_quotes($nl_flags->{'clm_start_type'}) ne "branch" ) {
    setup_logic_initial_conditions($opts, $nl_flags, $definition, $defaults, $nl, $physv);
  }
  setup_logic_spinup($opts,  $nl_flags, $definition, $defaults, $nl);
  setup_logic_supplemental_nitrogen($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_snowpack($opts,  $nl_flags, $definition, $defaults, $nl);
  setup_logic_fates($opts,  $nl_flags, $definition, $defaults, $nl);
  setup_logic_z0param($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_misc($opts, $nl_flags, $definition, $defaults, $nl);

  #########################################
  # namelist group: atm2lnd_inparm
  #########################################
  setup_logic_atm_forcing($opts,  $nl_flags, $definition, $defaults, $nl);

  #########################################
  # namelist group: lnd2atm_inparm
  #########################################
  setup_logic_lnd2atm($opts,  $nl_flags, $definition, $defaults, $nl);

  #########################################
  # namelist group: clm_humanindex_inparm #
  #########################################
  setup_logic_humanindex($opts,  $nl_flags, $definition, $defaults, $nl);

  #################################
  # namelist group: cnfire_inparm #
  #################################
  setup_logic_cnfire($opts,  $nl_flags, $definition, $defaults, $nl);

  ######################################
  # namelist group: cnprecision_inparm #
  ######################################
  setup_logic_cnprec($opts,  $nl_flags, $definition, $defaults, $nl);

  ###############################
  # namelist group: clmu_inparm #
  ###############################
  setup_logic_urban($opts,  $nl_flags, $definition, $defaults, $nl);

  ###############################
  # namelist group: crop_inparm #
  ###############################
  setup_logic_crop_inparm($opts,  $nl_flags, $definition, $defaults, $nl);

  ###############################
  # namelist group: tillage     #
  ###############################
  setup_logic_tillage($opts, $nl_flags, $definition, $defaults, $nl, $physv);

  ###############################
  # namelist group: ch4par_in   #
  ###############################
  setup_logic_methane($opts, $nl_flags, $definition, $defaults, $nl);
  setup_logic_c_isotope($opts, $nl_flags, $definition, $defaults, $nl);

  ###############################
  # namelist group: ndepdyn_nml #
  ###############################
  setup_logic_nitrogen_deposition($opts,  $nl_flags, $definition, $defaults, $nl);

  ##################################
  # namelist group: cnmresp_inparm #
  ##################################
  setup_logic_cnmresp($opts,  $nl_flags, $definition, $defaults, $nl);

  #################################
  # namelist group: nitrif_inparm #
  #################################
  setup_logic_nitrif_params( $nl_flags, $definition, $defaults, $nl );

  #############################################
  # namelist group: mineral_nitrogen_dynamics #
  #############################################
  setup_logic_mineral_nitrogen_dynamics( $opts, $nl_flags, $definition, $defaults, $nl );

  ####################################
  # namelist group: photosyns_inparm #
  ####################################
  setup_logic_photosyns($opts,  $nl_flags, $definition, $defaults, $nl);

  #################################
  # namelist group: popd_streams  #
  #################################
  setup_logic_popd_streams($opts,  $nl_flags, $definition, $defaults, $nl);

  ####################################
  # namelist group: urbantv_streams  #
  ####################################
  setup_logic_urbantv_streams($opts,  $nl_flags, $definition, $defaults, $nl);

  ##################################
  # namelist group: light_streams  #
  ##################################
  setup_logic_lightning_streams($opts,  $nl_flags, $definition, $defaults, $nl);

  #################################
  # namelist group: drydep_inparm #
  #################################
  setup_logic_dry_deposition($opts, $nl_flags, $definition, $defaults, $nl);

  #################################
  # namelist group: fire_emis_nl  #
  #################################
  setup_logic_fire_emis($opts, $nl_flags, $definition, $defaults, $nl);

  ######################################
  # namelist options for dust emissions
  ######################################
  setup_logic_dust_emis($opts, $nl_flags, $definition, $defaults, $nl);

  #################################
  # namelist group: megan_emis_nl #
  #################################
  setup_logic_megan($opts, $nl_flags, $definition, $defaults, $nl);

  ##################################
  # namelist group: lai_streams  #
  ##################################
  setup_logic_lai_streams($opts,  $nl_flags, $definition, $defaults, $nl);

  ##################################
  # namelist group: cropcal_streams  #
  ##################################
  setup_logic_cropcal_streams($opts,  $nl_flags, $definition, $defaults, $nl);

  ##########################################
  # namelist group: soil_moisture_streams  #
  ##########################################
  setup_logic_soilm_streams($opts,  $nl_flags, $definition, $defaults, $nl, $physv);

  ##################################
  # namelist group: bgc_shared
  ##################################
  setup_logic_bgc_shared($opts,  $nl_flags, $definition, $defaults, $nl, $physv);

  ##################################
  # namelist group: cnphenology
  ##################################
  setup_logic_cnphenology($opts,  $nl_flags, $definition, $defaults, $nl, $physv);

  #############################################
  # namelist group: soilwater_movement_inparm #
  #############################################
  setup_logic_soilwater_movement($opts,  $nl_flags, $definition, $defaults, $nl);

  #############################################
  # namelist group: rooting_profile_inparm    #
  #############################################
  setup_logic_rooting_profile($opts,  $nl_flags, $definition, $defaults, $nl);

  #############################
  # namelist group: cngeneral #
  #############################
  setup_logic_cngeneral($opts,  $nl_flags, $definition, $defaults, $nl);
  ####################################
  # namelist group: cnvegcarbonstate #
  ####################################
  setup_logic_cnvegcarbonstate($opts,  $nl_flags, $definition, $defaults, $nl);

  #############################################
  # namelist group: soil_resis_inparm #
  #############################################
  setup_logic_soil_resis($opts,  $nl_flags, $definition, $defaults, $nl);

  #############################################
  # namelist group: canopyfluxes_inparm #
  #############################################
  setup_logic_canopyfluxes($opts,  $nl_flags, $definition, $defaults, $nl);

  ##########################################################
  # namelist group: friction_velocity (after canopyfluxes) #
  ##########################################################
  setup_logic_friction_vel($opts,  $nl_flags, $definition, $defaults, $nl);

  #############################################
  # namelist group: canopyhydrology_inparm #
  #############################################
  setup_logic_canopyhydrology($opts,  $nl_flags, $definition, $defaults, $nl);

  #####################################
  # namelist group: clm_canopy_inparm #
  #####################################
  setup_logic_canopy($opts,  $nl_flags, $definition, $defaults, $nl);

  ########################################
  # namelist group: soilhydrology_inparm #
  ########################################
  setup_logic_hydrology_params($opts,  $nl_flags, $definition, $defaults, $nl);

  #####################################
  # namelist group: irrigation_inparm #
  #####################################
  setup_logic_irrigation_parameters($opts,  $nl_flags, $definition, $defaults, $nl);

  ########################################
  # namelist group: surfacealbedo_inparm #
  ########################################
  setup_logic_surfacealbedo($opts, $nl_flags, $definition, $defaults, $nl);

  ########################################
  # namelist group: water_tracers_inparm #
  ########################################
  setup_logic_water_tracers($opts, $nl_flags, $definition, $defaults, $nl);

  #######################################################################
  # namelist groups: clm_hydrology1_inparm and clm_soilhydrology_inparm #
  #######################################################################
  setup_logic_hydrology_switches($opts, $nl_flags, $definition, $defaults, $nl);

  #######################################################################
  # namelist group: scf_swenson_lawrence_2012_inparm                    #
  #######################################################################
  setup_logic_scf_SwensonLawrence2012($opts, $nl_flags, $definition, $defaults, $nl);

  #########################################
  # namelist group: clm_initinterp_inparm #
  #########################################
  setup_logic_initinterp($opts, $nl_flags, $definition, $defaults, $nl);

  ###############################
  # namelist group: exice_streams   #
  ###############################
  setup_logic_exice($opts, $nl_flags, $definition, $defaults, $nl);
}

#-------------------------------------------------------------------------------

sub setup_logic_site_specific {
  # site specific requirements
  my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

  # res check prevents polluting the namelist with an unnecessary
  # false variable for every run
  if ($nl_flags->{'res'} eq "1x1_vancouverCAN") {
     my $var = "use_vancouver";
     my $val = ".true.";
     my $group = $definition->get_group_name($var);
     $nl->set_variable_value($group, $var, $val);
  }

  # res check prevents polluting the namelist with an unnecessary
  # false variable for every run
  if ($nl_flags->{'res'} eq "1x1_mexicocityMEX") {
     my $var = "use_mexicocity";
     my $val = ".true.";
     my $group = $definition->get_group_name($var);
     $nl->set_variable_value($group, $var, $val);
  }

  if ($nl_flags->{'res'} eq "1x1_smallvilleIA") {
    if (! &value_is_true($nl_flags->{'use_cn'}) || ! &value_is_true($nl_flags->{'use_crop'})) {
      $log->fatal_error("1x1_smallvilleIA grids must use a compset with CN and CROP turned on.");
    }
  }

  if ($nl_flags->{'res'} eq "1x1_numaIA") {
    if (! &value_is_true($nl_flags->{'use_cn'}) || ! &value_is_true($nl_flags->{'use_crop'})) {
      $log->fatal_error("1x1_numaIA grids must use a compset with CN and CROP turned on.");
    }
  }
  #  Set compname
  add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'compname',
              'phys'=>$nl_flags->{'phys'} );
}

#-------------------------------------------------------------------------------

sub setup_logic_lnd_tuning {

  my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

  my $var    = "lnd_tuning_mode";
  if ( $opts->{$var} eq "default" ) {
     my %settings;
     $settings{'phys'} = $nl_flags->{'phys'};
     $nl_flags->{$var} = $defaults->get_value($var, \%settings );
  } else {
     $nl_flags->{$var} = $opts->{$var};
  }
  my $group = $definition->get_group_name($var);
  $nl->set_variable_value($group, $var, quote_string( $nl_flags->{$var} ) );
  if (  ! $definition->is_valid_value( $var, quote_string( $nl_flags->{$var}) ) ) {
    my @valid_values   = $definition->get_valid_values( $var );
    $log->fatal_error("$var has a value (".$nl_flags->{$var}.") that is NOT valid. Valid values are: @valid_values");
  }
  $log->verbose_message("Using $nl_flags->{$var} for lnd_tuning_mode");
  my $phys = $physv->as_string();
  if ( $nl_flags->{$var} !~ /^${phys}_/ ) {
     $log->fatal_error("First part of lnd_tuning_mode MUST match the CLM version you are using.");
  }
}


#-------------------------------------------------------------------------------

sub setup_logic_lnd_frac {

  my ($opts, $nl_flags, $definition, $defaults, $nl, $envxml_ref) = @_;

  #
  # fatmlndfrc is required for the MCT driver (or LILAC), but uneeded for NUOPC
  #
  my $var = "lnd_frac";
  if ( ($opts->{'driver'} eq "mct") || $opts->{'lilac'} ) {
     if ( defined($opts->{$var}) ) {
       if ( defined($nl->get_value('fatmlndfrc')) ) {
         $log->fatal_error("Can NOT set both -lnd_frac option (set via LND_DOMAIN_PATH/LND_DOMAIN_FILE " .
                     "env variables) AND fatmlndfrac on namelist");
       }
       if ( $opts->{$var} =~ /UNSET/ ) {
          $log->fatal_error("-lnd_frac was set as UNSET in the CTSM build-namelist set it with the env variables: LND_DOMAIN_PATH/LND_DOMAIN_FILE.");
       }
       my $lnd_frac = SetupTools::expand_xml_var( $opts->{$var}, $envxml_ref);
       add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fatmlndfrc','val'=>$lnd_frac );
     }

     # Get the fraction file
     if (defined $nl->get_value('fatmlndfrc')) {
       # do nothing - use value provided by config_grid.xml and clm.cpl7.template
     } else {
       $log->fatal_error("fatmlndfrc was NOT sent into CLM build-namelist.");
     }
  #
  # For the NUOPC driver neither lnd_frac nor fatmlndfrc need to be set
  #
  } else {
     if ( defined($opts->{$var}) ) {
       if ( $opts->{$var} !~ /UNSET/ ) {
          $log->fatal_error("$var should NOT be set for the NUOPC driver as it is unused" );
       }
     }
     if ( defined($nl->get_value('fatmlndfrc')) ) {
       $log->fatal_error("fatmlndfrac should NOT be set in the namelist for the NUOPC driver as it is unused" );
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_co2_type {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $var = "co2_type";
  if ( defined($opts->{$var}) ) {
    if ( ! defined($nl->get_value($var)) ) {
      add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'co2_type','val'=>"$opts->{'co2_type'}");
    } else {
      $log->fatal_error("co2_type set on namelist as well as -co2_type option.");
    }
  }
  add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'co2_type');
  if ( $nl->get_value('co2_type') =~ /constant/ ) {
    my $var = 'co2_ppmv';
    if ( defined($opts->{$var}) ) {
      if ( $opts->{$var} <= 0.0 ) {
        $log->fatal_error("co2_ppmv can NOT be less than or equal to zero.");
      }
      my $group = $definition->get_group_name($var);
      $nl->set_variable_value($group, $var, $opts->{$var});
    } else {
      add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'sim_year'=>$nl_flags->{'sim_year'},
                  'ssp_rcp'=>$nl_flags->{'ssp_rcp'} );
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_irrigate {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'irrigate',
                'use_crop'=>$nl_flags->{'use_crop'}, 'use_cndv'=>$nl_flags->{'use_cndv'},
                'sim_year'=>$nl_flags->{'sim_year'}, 'sim_year_range'=>$nl_flags->{'sim_year_range'}, );
  if ( &value_is_true($nl->get_value('irrigate') ) ) {
     $nl_flags->{'irrigate'} = ".true.";
     if ( $nl_flags->{'sim_year'} eq "PtVg" ) {
        $log->fatal_error("irrigate=TRUE does NOT make sense with the Potential Vegetation dataset, leave irrigate=FALSE");
     }
  } else {
     $nl_flags->{'irrigate'} = ".false.";
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_start_type {
  my ($opts, $nl_flags, $nl) = @_;

  my $var = "start_type";
  my $drv_start_type = $nl->get_value($var);
  my $my_start_type  = $nl_flags->{'clm_start_type'};

  if ( $my_start_type =~ /branch/ ) {
    if (not defined $nl->get_value('nrevsn')) {
      $log->fatal_error("nrevsn is required for a branch type.");
    }
    if (defined $nl->get_value('use_init_interp')) {
       if ( &value_is_true($nl->get_value('use_init_interp') ) ) {
         # Always print this warning, but don't stop if it happens
         print "\nWARNING: use_init_interp will NOT happen for a branch case.\n\n";
       }
    }
  } else {
    if (defined $nl->get_value('nrevsn')) {
      $log->fatal_error("nrevsn should ONLY be set for a branch type.");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_decomp_performance {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  # Set the number of segments per clump
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'nsegspc', 'hgrid'=>$nl_flags->{'res'});
}

#-------------------------------------------------------------------------------

sub setup_logic_roughness_methods {
  my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'z0param_method',
              'phys'=>$nl_flags->{'phys'} );

  my $var = remove_leading_and_trailing_quotes( $nl->get_value("z0param_method") );
  if ( $var ne "Meier2022" && $var ne "ZengWang2007" ) {
    $log->fatal_error("$var is incorrect entry for the namelist variable z0param_method; expected Meier2022 or ZengWang2007");
  }
  my $phys = $physv->as_string();
  if ( $phys eq "clm4_5" || $phys eq "clm5_0" ) {
    if ( $var eq "Meier2022" ) {
      $log->fatal_error("z0param_method = $var and phys = $phys, but this method has been tested only with clm5_1 and later versions; to use with earlier versions, disable this error, and add Meier2022 parameters to the corresponding params file");
    }
  }
}
#-------------------------------------------------------------------------------

sub setup_logic_snicar_methods {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snicar_snw_shape' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snicar_solarspec' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snicar_dust_optics' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snicar_numrad_snw' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snicar_snobc_intmix' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snicar_snodst_intmix' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snicar_use_aerosol' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_snicar_frc' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'do_sno_oc' );

  # Error checking in loop
  my %supportedSettings = ( 'snicar_solarspec' => "'mid_latitude_winter'", 'snicar_dust_optics' => "'sahara'", 'snicar_numrad_snw' => '5', 'snicar_snodst_intmix' => '.false.', 'snicar_use_aerosol' => '.true.', 'do_sno_oc' => '.false.' );
  keys %supportedSettings;
  while ( my ($key, $val) = each %supportedSettings ) {
    my $var = $nl->get_value($key);
    if ( $var ne $val ) {
      $log->warning("$key=$val is the supported option; $var is EXPERIMENTAL, UNSUPPORTED, and UNTESTED!");
    }
  }

  # Error checking not in loop
  my $key1 = 'snicar_snw_shape';
  my $var1 = $nl->get_value($key1);
  my $val1a = "'sphere'";  # supported value for this option
  my $val1b = "'hexagonal_plate'";  # supported value for this option
  if (($var1 ne $val1a) && ($var1 ne $val1b)) {
    $log->warning("$key1=$val1a and $val1b are supported; $var1 is EXPERIMENTAL, UNSUPPORTED, and UNTESTED!");
  }

  # snicar_snobc_intmix and snicar_snodst_intmix cannot both be true, however, they can both be false
  my $key1 = 'snicar_snobc_intmix';
  my $key2 = 'snicar_snodst_intmix';
  my $var1 = $nl->get_value($key1);
  my $var2 = $nl->get_value($key2);
  my $val2 = $supportedSettings{$key2};  # supported value for this option
  if (($var1 eq $var2) && ($var2 ne $val2)) {
    $log->warning("$key1 = $var1 and $key2 = $var2 do not work together!");
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_snow {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snow_thermal_cond_method' );

  my $var = $nl->get_value('snow_thermal_cond_method');
  if ( $var ne "'Jordan1991'" && $var ne "'Sturm1997'" ) {
    $log->fatal_error("$var is incorrect entry for the namelist variable snow_thermal_cond_method; expected Jordan1991 or Sturm1997");
  }

  my $numrad_snw = $nl->get_value('snicar_numrad_snw');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fsnowoptics',
                'snicar_numrad_snw' => $numrad_snw);
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fsnowaging' );
}

#-------------------------------------------------------------------------------

sub setup_logic_glacier {
  #
  # Glacier multiple elevation class options
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl, $envxml_ref) = @_;

  my $clm_upvar = "GLC_TWO_WAY_COUPLING";
  # glc_do_dynglacier is set via GLC_TWO_WAY_COUPLING; it cannot be set via
  # user_nl_clm (this is because we might eventually want the coupler and glc
  # to also respond to GLC_TWO_WAY_COUPLING, by not bothering to send / map
  # these fields - so we want to ensure that CLM is truly listening to this
  # shared xml variable and not overriding it)
  my $var = "glc_do_dynglacier";
  my $val = logical_to_fortran($envxml_ref->{$clm_upvar});
  add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'val'=>$val);
  if (lc($nl->get_value($var)) ne lc($val)) {
     $log->fatal_error("glc_do_dynglacier can only be set via the env variable $clm_upvar: it can NOT be set in user_nl_clm");
  }

  my $var = "maxpatch_glc";
  add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'val'=>$nl_flags->{'glc_nec'} );

  my $val = $nl->get_value($var);
  if ( $val != $nl_flags->{'glc_nec'} ) {
    $log->fatal_error("$var set to $val does NOT agree with -glc_nec argument of $nl_flags->{'glc_nec'} (set with GLC_NEC env variable)");
  }

  if ( $nl_flags->{'glc_nec'} < 1 ) {
     $log->fatal_error("GLC_NEC must be at least 1.");
  }

  add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'glc_snow_persistence_max_days');

  add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'albice');
  add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'glacier_region_behavior',
              'glc_use_antarctica'=>$opts->{'glc_use_antarctica'});
  add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'glacier_region_melt_behavior');
  add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'glacier_region_ice_runoff_behavior');
}

#-------------------------------------------------------------------------------

sub setup_logic_params_file {
  # get param data
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'paramfile',
              'phys'=>$nl_flags->{'phys'},
              'use_flexibleCN'=>$nl_flags->{'use_flexibleCN'} );
}

#-------------------------------------------------------------------------------

sub setup_logic_create_crop_landunit {
  # Create crop land unit
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $var = 'create_crop_landunit';
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
              'use_fates'=>$nl_flags->{'use_fates'} );
  if ( &value_is_true($nl_flags->{'use_fates'}) && &value_is_true($nl->get_value($var)) ) {
     $log->fatal_error( "$var is true and yet use_fates is being set, which contradicts that (use_fates requires $var to be .false." );
  }
  if ( (! &value_is_true($nl_flags->{'use_fates'})) && (! &value_is_true($nl->get_value($var))) ) {
     $log->fatal_error( "$var is false which is ONLY allowed when FATES is being used" );
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_subgrid {
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   my $var = 'run_zero_weight_urban';
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var);
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'convert_ocean_to_land');
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'collapse_urban',
               'structure'=>$nl_flags->{'structure'});
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'n_dom_landunits',
               'structure'=>$nl_flags->{'structure'});
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'n_dom_pfts',
               'structure'=>$nl_flags->{'structure'});
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'toosmall_soil');
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'toosmall_crop');
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'toosmall_glacier');
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'toosmall_lake');
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'toosmall_wetland');
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'toosmall_urban');

   if ( &value_is_true($nl_flags->{'use_fates'}) && $nl->get_value('n_dom_pfts') != 0 ) {
      $log->fatal_error( "FATES and n_dom_pfts can NOT be set at the same time" );
   }
}

#-------------------------------------------------------------------------------

sub setup_logic_cnfire {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my @fire_consts = ( "rh_low", "rh_hgh", "bt_min", "bt_max", "cli_scale", "boreal_peatfire_c", "non_boreal_peatfire_c",
                      "pot_hmn_ign_counts_alpha", "cropfire_a1", "occur_hi_gdp_tree", "lfuel", "ufuel",
                      "cmb_cmplt_fact_litter", "cmb_cmplt_fact_cwd" );
  if ( &value_is_true($nl->get_value('use_cn')) ) {
     foreach my $item ( @fire_consts ) {
        if ( ! &value_is_true($nl_flags->{'cnfireson'} ) ) {
           if ( defined($nl->get_value($item)) ) {
              $log->fatal_error( "fire_method is no_fire and yet $item is being set, which contradicts that" );
           }
        } else {
           my $fire_method = remove_leading_and_trailing_quotes( $nl->get_value('fire_method') );
           add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $item,
                       'lnd_tuning_mode'=>$nl_flags->{'lnd_tuning_mode'},
                       'fire_method'=>$fire_method );
        }
     }
  } else {
     foreach my $item ( @fire_consts ) {
        if ( defined($nl->get_value($item)) ) {
           $log->fatal_error( "CN is off which implies that cnfire is off and yet a fire constant ($item) is being set, which contradicts that" );
        }
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_cnprec {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( &value_is_true($nl_flags->{'use_cn'}) ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults,
                 $nl, 'ncrit', 'use_cn'=>$nl_flags->{'use_cn'});
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults,
                 $nl, 'cnegcrit', 'use_cn'=>$nl_flags->{'use_cn'});
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults,
                 $nl, 'nnegcrit', 'use_cn'=>$nl_flags->{'use_cn'});
  }
}
#-------------------------------------------------------------------------------

sub setup_logic_humanindex {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'calc_human_stress_indices');
}

#-------------------------------------------------------------------------------

sub setup_logic_urban {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'building_temp_method');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'urban_hac');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'urban_traffic');
}

#-------------------------------------------------------------------------------

sub setup_logic_crop_inparm {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( &value_is_true($nl->get_value('use_crop')) ) {
     my $maptype = 'baset_mapping';
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $maptype,
                 'use_crop'=>$nl->get_value('use_crop') );
     my $baset_mapping = remove_leading_and_trailing_quotes( $nl->get_value($maptype) );
     if ( $baset_mapping eq "varytropicsbylat" ) {
        add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, "baset_latvary_slope",
                    'use_crop'=>$nl->get_value('use_crop'), $maptype=>$baset_mapping );
        add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, "baset_latvary_intercept",
                    'use_crop'=>$nl->get_value('use_crop'), $maptype=>$baset_mapping );
     } else {
        error_if_set( $nl, "Can only be set if $maptype == varytropicsbylat", "baset_latvary_slope", "baset_latvary_intercept" );
     }
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, "initial_seed_at_planting",
                 'use_crop'=>$nl->get_value('use_crop') );
     
     my $crop_residue_removal_frac = $nl->get_value('crop_residue_removal_frac');
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'crop_residue_removal_frac' );
     if ( $crop_residue_removal_frac < 0.0 or $crop_residue_removal_frac > 1.0 ) {
        $log->fatal_error("crop_residue_removal_frac must be in range [0, 1]");
     }
  } else {
     error_if_set( $nl, "Can NOT be set without crop on", "baset_mapping", "baset_latvary_slope", "baset_latvary_intercept", "crop_residue_removal_frac" );
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'crop_fsat_equals_zero' );
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_tillage {
  my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'tillage_mode',
              'use_crop'=>$nl_flags->{'use_crop'}, 'phys'=>$physv->as_string() );

  my $tillage_mode = remove_leading_and_trailing_quotes( $nl->get_value( "tillage_mode" ) );
  if ( $tillage_mode ne "off" && $tillage_mode ne "" && not &value_is_true($nl_flags->{'use_crop'}) ) {
      $log->fatal_error( "Tillage only works on crop columns, so use_crop must be true if tillage is enabled." );
  }
}

#-------------------------------------------------------------------------------
sub error_if_set {
   # do a fatal_error and exit if any of the input variable names are set
   my ($nl, $error, @list) = @_;
   foreach my $var ( @list ) {
      if ( defined($nl->get_value($var)) ) {
        $log->fatal_error( "$var $error" );
      }
   }
}

#-------------------------------------------------------------------------------

sub setup_logic_soilstate {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'organic_frac_squared' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_bedrock',
              'use_fates'=>$nl_flags->{'use_fates'}, 'vichydro'=>$nl_flags->{'vichydro'} );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_excess_ice'); # excess ice flag should be read before stream vars

  my $var1 = "soil_layerstruct_predefined";
  my $var2 = "soil_layerstruct_userdefined";
  my $var3 = "soil_layerstruct_userdefined_nlevsoi";
  my $soil_layerstruct_predefined = $nl->get_value($var1);
  my $soil_layerstruct_userdefined = $nl->get_value($var2);
  my $soil_layerstruct_userdefined_nlevsoi = $nl->get_value($var3);

  if (defined($soil_layerstruct_userdefined)) {
    if (defined($soil_layerstruct_predefined)) {
      $log->fatal_error("You have set both soil_layerstruct_userdefined and soil_layerstruct_predefined in your namelist; model cannot determine which to use");
    }
    if (not defined($soil_layerstruct_userdefined_nlevsoi)) {
      $log->fatal_error("You have set soil_layerstruct_userdefined and NOT set soil_layerstruct_userdefined_nlevsoi in your namelist; both MUST be set");
    }
  } else {
    if (defined($soil_layerstruct_userdefined_nlevsoi)) {
      $log->fatal_error("You have set soil_layerstruct_userdefined_nlevsoi and NOT set soil_layerstruct_userdefined in your namelist; EITHER set both OR neither; in the latter case soil_layerstruct_predefined will be assigned a default value");
    } else {
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'soil_layerstruct_predefined',
                  'structure'=>$nl_flags->{'structure'});
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_demand {
  #
  # Deal with options that the user has said are required...
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my %settings;
  $settings{'hgrid'}          = $nl_flags->{'res'};
  $settings{'sim_year'}       = $nl_flags->{'sim_year'};
  $settings{'sim_year_range'} = $nl_flags->{'sim_year_range'};
  $settings{'use_vichydro'}   = $nl_flags->{'use_vichydro'};
  $settings{'mask'}           = $nl_flags->{'mask'};
  $settings{'crop'}           = $nl_flags->{'crop'};
  $settings{'ssp_rcp'}        = $nl_flags->{'ssp_rcp'};
  $settings{'glc_nec'}        = $nl_flags->{'glc_nec'};
  # necessary for demand to be set correctly (flanduse_timeseries requires
  # use_crop, maybe other options require other flags?)!
  $settings{'irrigate'}            = $nl_flags->{'irrigate'};
  $settings{'use_cn'}              = $nl_flags->{'use_cn'};
  $settings{'use_cndv'}            = $nl_flags->{'use_cndv'};
  $settings{'use_lch4'}            = $nl_flags->{'use_lch4'};
  $settings{'use_nitrif_denitrif'} = $nl_flags->{'use_nitrif_denitrif'};
  $settings{'use_crop'}            = $nl_flags->{'use_crop'};
  $settings{'neon'}                = $nl_flags->{'neon'};

  my $demand = $nl->get_value('clm_demand');
  if (defined($demand)) {
    $demand =~ s/\'//g;   # Remove quotes
    if ( $demand =~ /.+/ ) {
      $opts->{'clm_demand'} .= ",$demand";
    }
  }

  $demand = $defaults->get_value('clm_demand', \%settings);
  if (defined($demand)) {
    $demand =~ s/\'//g;   # Remove quotes
    if ( $demand =~ /.+/ ) {
      $opts->{'clm_demand'} .= ",$demand";
    }
  }

  my @demandlist = split( ",", $opts->{'clm_demand'} );
  foreach my $item ( @demandlist ) {
    if ( $item eq "null" ) {
      next;
    }
    if ( $item eq "finidat" ) {
        $log->fatal_error( "Do NOT put findat in the clm_demand list, set the clm_start_type=startup so initial conditions are required");
    }
    # For landuse.timeseries try with crop on first eise try with exact settings
    # Logic for this is identical for fsurdat
    if ( $item eq "flanduse_timeseries" ) {
       $settings{'use_crop'} = ".true.";
       $settings{'nofail'}   = 1;
    }
    add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $item, %settings );
    if ( $item eq "flanduse_timeseries" ) {
       $settings{'nofail'}   = 0;
       $settings{'use_crop'} = $nl_flags->{'use_crop'};
       if ( ! defined($nl->get_value( $item )) ) {
          add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $item, %settings );
       }
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_surface_dataset {
  #
  # Get surface dataset after flanduse_timeseries so that we can get surface data
  # consistent with it
  # MUST BE AFTER: setup_logic_demand which is where flanduse_timeseries is set
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl, $xmlvar_ref) = @_;

  $nl_flags->{'flanduse_timeseries'} = "null";
  my $flanduse_timeseries = $nl->get_value('flanduse_timeseries');
  if (defined($flanduse_timeseries)) {
    $flanduse_timeseries =~ s!(.*)/!!;
    $flanduse_timeseries =~ s/\'//;
    $flanduse_timeseries =~ s/\"//;
    if ( $flanduse_timeseries ne "" ) {
      $nl_flags->{'flanduse_timeseries'} = $flanduse_timeseries;
    }
  }
  $flanduse_timeseries = $nl_flags->{'flanduse_timeseries'};

  if ($flanduse_timeseries ne "null" && &value_is_true($nl_flags->{'use_cndv'}) ) {
     $log->fatal_error( "dynamic PFT's (setting flanduse_timeseries) are incompatible with dynamic vegetation (use_cndv=.true)." );
  }
  #
  # Always get the crop version of the datasets now and let the code turn it into the form desired
  # Provided this isn't with FATES on
  #
  my $var = "fsurdat";
  if ( !  &value_is_true($nl_flags->{'use_fates'}) ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                 'hgrid'=>$nl_flags->{'res'}, 'ssp_rcp'=>$nl_flags->{'ssp_rcp'},
                 'neon'=>$nl_flags->{'neon'}, 'neonsite'=>$nl_flags->{'neonsite'},
                 'sim_year'=>$nl_flags->{'sim_year'}, 'use_vichydro'=>$nl_flags->{'use_vichydro'},
                 'use_crop'=>".true.", 'use_fates'=>$nl_flags->{'use_fates'}, 'nofail'=>1);
  }
  # If didn't find the crop version check for the exact match
  my $fsurdat = $nl->get_value($var);
  if ( ! defined($fsurdat) ) {
     if ( !  &value_is_true($nl_flags->{'use_fates'}) ) {
        $log->verbose_message( "Crop version of $var NOT found, searching for an exact match" );
     }
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                 'hgrid'=>$nl_flags->{'res'}, 'ssp_rcp'=>$nl_flags->{'ssp_rcp'}, 'use_vichydro'=>$nl_flags->{'use_vichydro'},
                 'sim_year'=>$nl_flags->{'sim_year'}, 'use_fates'=>$nl_flags->{'use_fates'},
                 'neon'=>$nl_flags->{'neon'}, 'neonsite'=>$nl_flags->{'neonsite'},
                 'use_crop'=>$nl_flags->{'use_crop'} );
  }
  #
  # Expand the XML variables for NEON cases so that NEONSITE will be used
  #
  if ( &value_is_true($nl_flags->{'neon'}) ) {
     my $fsurdat = $nl->get_value($var);
     my $newval = SetupTools::expand_xml_var( $fsurdat, $xmlvar_ref );
     if ( $newval ne $fsurdat ) {
        my $group = $definition->get_group_name($var);
        $nl->set_variable_value($group, $var, $newval);
        $log->verbose_message( "This is a NEON site and the fsurdat file selected is: $newval" );
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_initial_conditions {
  # Initial conditions
  # The initial date is an attribute in the defaults file which should be matched unless
  # the user explicitly requests to ignore the initial date via the -ignore_ic_date option,
  # or just ignore the year of the initial date via the -ignore_ic_year option.
  #
  # MUST BE AFTER: setup_logic_demand   which is where flanduse_timeseries is set
  #         AFTER: setup_logic_irrigate which is where irrigate is set
  my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

  my $var = "finidat";
  my $finidat = $nl->get_value($var);
  if ( $nl_flags->{'clm_start_type'} =~ /cold/ ) {
    if (defined $finidat ) {
      $log->warning("setting $var (either explicitly in your user_nl_clm or by doing a hybrid or branch RUN_TYPE)\n is incomptable with using a cold start" .
              " (by setting CLM_FORCE_COLDSTART=on)." );
      $log->warning("Overridding input $var file with one specifying that this is a cold start from arbitrary initial conditions." );
      my $group = $definition->get_group_name($var);
      $nl->set_variable_value($group, $var, "' '" );
    }
    add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
                $var, 'val'=>"' '", 'no_abspath'=>1);
    $finidat = $nl->get_value($var);
  } elsif ( defined $finidat ) {
    if ( string_is_undef_or_empty($finidat) ) {
      print "You are setting $var to blank which signals arbitrary initial conditions.\n";
      print "But, CLM_FORCE_COLDSTART is off which is a contradiction. For arbitrary initial conditions just use the CLM_FORCE_COLDSTART option\n";
      $log->fatal_error("To do a cold-start set ./xmlchange CLM_FORCE_COLDSTART=on, and remove the setting of $var in the user_nl_clm file");
    }
  }
  my $useinitvar = "use_init_interp";

  my %settings;
  my $use_init_interp_default = $nl->get_value($useinitvar);
  $settings{$useinitvar} = $use_init_interp_default;
  if ( string_is_undef_or_empty( $use_init_interp_default ) ) {
    $use_init_interp_default = $defaults->get_value($useinitvar, \%settings);
    $settings{$useinitvar} = ".false.";
  }
  if (not defined $finidat ) {
    my $ic_date = $nl->get_value('start_ymd');
    my $st_year = $nl_flags->{'st_year'};
    my $nofail = 1;
    $settings{'hgrid'}   = $nl_flags->{'res'};
    $settings{'phys'}    = $physv->as_string();
    $settings{'nofail'}  = $nofail;
    my $fsurdat          = $nl->get_value('fsurdat');
    $fsurdat             =~ s!(.*)/!!;
    $settings{'fsurdat'} = $fsurdat;
    $settings{'do_transient_pfts'} = $nl->get_value('do_transient_pfts');
    #
    # If not transient use sim_year, otherwise use date
    #
    if (string_is_undef_or_empty($nl->get_value('flanduse_timeseries'))) {
       $settings{'sim_year'}     = $nl_flags->{'sim_year'};
       $opts->{'ignore_ic_year'} = 1;
    } else {
       $settings{'sim_year'}     = $st_year;
    }
    foreach my $item ( "mask", "maxpft", "irrigate", "glc_nec", "use_crop", "use_cn", "use_cndv",
                       "use_fates",
                       "lnd_tuning_mode",
                     ) {
       $settings{$item}    = $nl_flags->{$item};
    }
    if ($opts->{'ignore_ic_date'}) {
      if ( &value_is_true($nl_flags->{'use_crop'}) ) {
        $log->warning("using ignore_ic_date is incompatable with crop! If you choose to ignore this error, " .
                      "the counters since planting for crops will be messed up. \nSo you should ignore at " .
                      "least the first season for crops. And since it will impact the 20 year means, ideally the " .
                      "first 20 years should be ignored.");
      }
    } elsif ($opts->{'ignore_ic_year'}) {
       $settings{'ic_md'} = $ic_date;
    } else {
       $settings{'ic_ymd'} = $ic_date;
    }
    my $try = 0;
    my $done = 2;
    do {
       $try++;
       add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, %settings );
       # If couldn't find a matching finidat file, check if can turn on interpolation and try to find one again
       $finidat = $nl->get_value($var);
       if (not defined $finidat) {
          # Delete any date settings, except for crop
          delete( $settings{'ic_ymd'} );
          delete( $settings{'ic_md'}  );
          #if ( &value_is_true($nl_flags->{'use_crop'}) ) {
             #$settings{'ic_md'} = $ic_date;
          #}
          add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, "init_interp_sim_years" );
          add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, "init_interp_how_close" );
          #
          # Figure out which sim_year has a usable finidat file that is closest to the desired one
          #
          my $close = $nl->get_value("init_interp_how_close");
          my $closest_sim_year = undef;
          my @sim_years = split( /,/, $nl->get_value("init_interp_sim_years") );
SIMYR:    foreach my $sim_yr ( @sim_years ) {
             my $how_close = undef;
             if ( $nl_flags->{'sim_year'} eq "PtVg" ) {
                $how_close = abs(1850 - $sim_yr);
             } elsif ( $nl_flags->{'flanduse_timeseries'} eq "null" ) {
                $how_close = abs($nl_flags->{'sim_year'} - $sim_yr);
             } else {
                $how_close = abs($st_year - $sim_yr);
             }
             if ( ($sim_yr == $sim_years[-1]) || (($how_close < $nl->get_value("init_interp_how_close")) && ($how_close < $close)) ) {
                my $group = $definition->get_group_name($useinitvar);
                $settings{'sim_year'} = $sim_yr;
                $settings{$useinitvar} = $defaults->get_value($useinitvar, \%settings);
                if ( ! defined($settings{$useinitvar}) ) {
                   $settings{$useinitvar} = $use_init_interp_default;
                }
                if ( &value_is_true($settings{$useinitvar}) ) {

                   if ( ($how_close < $nl->get_value("init_interp_how_close")) && ($how_close < $close) ) {
                      $close = $how_close;
                      $closest_sim_year = $sim_yr;
                   }
                }
             }
          }    # SIMYR:
          $settings{'sim_year'} = $closest_sim_year;
          # Add options set here to the "$set" variable below...
          add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $useinitvar,
                      'use_cndv'=>$nl_flags->{'use_cndv'}, 'phys'=>$physv->as_string(), 'hgrid'=>$nl_flags->{'res'},
                      'sim_year'=>$settings{'sim_year'}, 'nofail'=>1, 'lnd_tuning_mode'=>$nl_flags->{'lnd_tuning_mode'},
                      'use_fates'=>$nl_flags->{'use_fates'} );
          $settings{$useinitvar} = $nl->get_value($useinitvar);
          if ( ! &value_is_true($nl->get_value($useinitvar) ) ) {
             if ( $nl_flags->{'clm_start_type'} =~ /startup/ ) {
                my $err_msg = "clm_start_type is startup so an initial conditions ($var) file is required,";
                if ( defined($use_init_interp_default) ) {
                   $log->fatal_error($err_msg." but can't find one without $useinitvar being set to true, change it to true in your user_nl_clm file in your case");
                } else {
                   my $set = "Relevent settings: use_cndv = ". $nl_flags->{'use_cndv'} . " phys = " . 
                              $physv->as_string() . " hgrid = " . $nl_flags->{'res'} . " sim_year = " . 
                              $settings{'sim_year'} . " lnd_tuning_mode = " . $nl_flags->{'lnd_tuning_mode'} .
                              "use_fates = " . $nl_flags->{'use_fates'};
                   $log->fatal_error($err_msg." but the default setting of $useinitvar is false, so set both $var to a startup file and $useinitvar==TRUE, or developers should modify the namelist_defaults file".$set);
                }
             }
          } else {
             my $stat = add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, "init_interp_attributes",
                                 'sim_year'=>$settings{'sim_year'}, 'use_cndv'=>$nl_flags->{'use_cndv'},
                                 'glc_nec'=>$nl_flags->{'glc_nec'}, 'use_fates'=>$nl_flags->{'use_fates'},
                                 'hgrid'=>$nl_flags->{'res'},
                                 'use_cn'=>$nl_flags->{'use_cn'}, 'lnd_tuning_mode'=>$nl_flags->{'lnd_tuning_mode'}, 'nofail'=>1 );
             if ( $stat ) {
                $log->fatal_error("$useinitvar is NOT synchronized with init_interp_attributes in the namelist_defaults file, this should be corrected there");
             }
             my $attributes = $nl->get_value("init_interp_attributes");
             my $attributes_string = remove_leading_and_trailing_quotes($attributes);
             foreach my $pair ( split( /\s/, $attributes_string) ) {
                if ( $pair =~ /^([a-z_]+)=([a-zA-Z._0-9]+)$/ ) {
                   $settings{$1} = $2;
                } else {
                   $log->fatal_error("Problem interpreting init_interp_attributes from the namelist_defaults file: $pair");
                }
             }
          }
       } else {
         $try = $done
       }
    } while( ($try < $done) && (not defined $finidat ) );
    if ( not defined $finidat  ) {
      my $group = $definition->get_group_name($var);
      $nl->set_variable_value($group, $var, "' '" );
    }
  }
  $finidat = $nl->get_value($var);
  if ( &value_is_true($nl->get_value($useinitvar) ) && string_is_undef_or_empty($finidat) ) {
     if ( ! defined($use_init_interp_default) ) {
        $log->fatal_error("You set $useinitvar but a $var file could not be found for this case, try setting $var explicitly, and/or removing the setting for $useinitvar" );
     } else {
        $log->fatal_error("$useinitvar is being set for you but a $var was not found, so $useinitvar, init_interp_attributes, and finidat must not be set correctly for this configuration in the namelist_default file" );
     }
  }
} # end initial conditions

#-------------------------------------------------------------------------------

sub setup_logic_dynamic_subgrid {
   #
   # Options controlling which parts of flanduse_timeseries to use
   #
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   setup_logic_do_transient_pfts($opts, $nl_flags, $definition, $defaults, $nl);
   setup_logic_do_transient_crops($opts, $nl_flags, $definition, $defaults, $nl);
   setup_logic_do_transient_lakes($opts, $nl_flags, $definition, $defaults, $nl);
   setup_logic_do_transient_urban($opts, $nl_flags, $definition, $defaults, $nl);
   setup_logic_do_harvest($opts, $nl_flags, $definition, $defaults, $nl);
   setup_logic_do_grossunrep($opts, $nl_flags, $definition, $defaults, $nl);

   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'reset_dynbal_baselines');
   if ( &value_is_true($nl->get_value('reset_dynbal_baselines')) &&
        &remove_leading_and_trailing_quotes($nl_flags->{'clm_start_type'}) eq "branch") {
      $log->fatal_error("reset_dynbal_baselines has no effect in a branch run");
   }
}

sub setup_logic_do_transient_pfts {
   #
   # Set do_transient_pfts default value, and perform error checking on do_transient_pfts
   #
   # Assumes the following are already set in the namelist (although it's okay
   # for them to be unset if that will be their final state):
   # - flanduse_timeseries
   # - use_cndv
   # - use_fates
   #
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   my $var = 'do_transient_pfts';

   # Start by assuming a default value of '.true.'. Then check a number of
   # conditions under which do_transient_pfts cannot be true. Under these
   # conditions: (1) set default value to '.false.'; (2) make sure that the
   # value is indeed false (e.g., that the user didn't try to set it to true).

   my $default_val = ".true.";

   # cannot_be_true will be set to a non-empty string in any case where
   # do_transient_pfts should not be true; if it turns out that
   # do_transient_pfts IS true in any of these cases, a fatal error will be
   # generated
   my $cannot_be_true = "";

   my $n_dom_pfts = $nl->get_value( 'n_dom_pfts' );
   my $n_dom_landunits = $nl->get_value( 'n_dom_landunits' );
   my $toosmall_soil = $nl->get_value( 'toosmall_soil' );
   my $toosmall_crop = $nl->get_value( 'toosmall_crop' );
   my $toosmall_glacier = $nl->get_value( 'toosmall_glacier' );
   my $toosmall_lake = $nl->get_value( 'toosmall_lake' );
   my $toosmall_wetland = $nl->get_value( 'toosmall_wetland' );
   my $toosmall_urban = $nl->get_value( 'toosmall_urban' );

   if (string_is_undef_or_empty($nl->get_value('flanduse_timeseries'))) {
      $cannot_be_true = "$var can only be set to true when running a transient case (flanduse_timeseries non-blank)";
   } elsif (&value_is_true($nl->get_value('use_cndv'))) {
      $cannot_be_true = "$var cannot be combined with use_cndv";
   } elsif (&value_is_true($nl->get_value('use_fates'))) {
      $cannot_be_true = "$var cannot be combined with use_fates";
   } elsif (&value_is_true($nl->get_value('use_hillslope'))) {
      $cannot_be_true = "$var cannot be combined with use_hillslope";
   }

   if ($cannot_be_true) {
      $default_val = ".false.";
   }

   if (!$cannot_be_true) {
      # Note that, if the variable cannot be true, we don't call add_default
      # - so that we don't clutter up the namelist with variables that don't
      # matter for this case
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, val=>$default_val);
   }

   # Make sure the value is false when it needs to be false - i.e., that the
   # user hasn't tried to set a true value at an inappropriate time.

   if (&value_is_true($nl->get_value($var)) && $cannot_be_true) {
      $log->fatal_error($cannot_be_true);
   }

   # if do_transient_pfts is .true. and any of these (n_dom_* or toosmall_*)
   # are > 0 or collapse_urban = .true., then give fatal error
   if (&value_is_true($nl->get_value($var))) {
      if (&value_is_true($nl->get_value('collapse_urban'))) {
         $log->fatal_error("$var cannot be combined with collapse_urban");
      }
      if ($n_dom_pfts > 0 || $n_dom_landunits > 0 || $toosmall_soil > 0 || $toosmall_crop > 0 || $toosmall_glacier > 0 || $toosmall_lake > 0 || $toosmall_wetland > 0 || $toosmall_urban > 0) {
         $log->fatal_error("$var cannot be combined with any of the of the following > 0: n_dom_pfts > 0, n_dom_landunit > 0, toosmall_soi > 0._r8, toosmall_crop > 0._r8, toosmall_glacier > 0._r8, toosmall_lake > 0._r8, toosmall_wetland > 0._r8, toosmall_urban > 0._r8");
      }
   }
}

sub setup_logic_do_transient_crops {
   #
   # Set do_transient_crops default value, and perform error checking on do_transient_crops
   #
   # Assumes the following are already set in the namelist (although it's okay
   # for them to be unset if that will be their final state):
   # - flanduse_timeseries
   # - use_crop
   # - use_fates
   #
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   my $var = 'do_transient_crops';

   # Start by assuming a default value of '.true.'. Then check a number of
   # conditions under which do_transient_crops cannot be true. Under these
   # conditions: (1) set default value to '.false.'; (2) make sure that the
   # value is indeed false (e.g., that the user didn't try to set it to true).

   my $default_val = ".true.";

   # cannot_be_true will be set to a non-empty string in any case where
   # do_transient_crops should not be true; if it turns out that
   # do_transient_crops IS true in any of these cases, a fatal error will be
   # generated
   my $cannot_be_true = "";

   my $n_dom_pfts = $nl->get_value( 'n_dom_pfts' );
   my $n_dom_landunits = $nl->get_value( 'n_dom_landunits' );
   my $toosmall_soil = $nl->get_value( 'toosmall_soil' );
   my $toosmall_crop = $nl->get_value( 'toosmall_crop' );
   my $toosmall_glacier = $nl->get_value( 'toosmall_glacier' );
   my $toosmall_lake = $nl->get_value( 'toosmall_lake' );
   my $toosmall_wetland = $nl->get_value( 'toosmall_wetland' );
   my $toosmall_urban = $nl->get_value( 'toosmall_urban' );

   if (string_is_undef_or_empty($nl->get_value('flanduse_timeseries'))) {
      $cannot_be_true = "$var can only be set to true when running a transient case (flanduse_timeseries non-blank)";
   } elsif (&value_is_true($nl->get_value('use_fates'))) {
      # In principle, use_fates should be compatible with
      # do_transient_crops. However, this hasn't been tested, so to be safe,
      # we are not allowing this combination for now.
      $cannot_be_true = "$var has not been tested with FATES, so for now these two options cannot be combined";
   } elsif (&value_is_true($nl->get_value('use_hillslope'))) {
      $cannot_be_true = "$var cannot be combined with use_hillslope";
   }

   if ($cannot_be_true) {
      $default_val = ".false.";
   }

   if (!$cannot_be_true) {
      # Note that, if the variable cannot be true, we don't call add_default
      # - so that we don't clutter up the namelist with variables that don't
      # matter for this case
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, val=>$default_val);
   }

   # Make sure the value is false when it needs to be false - i.e., that the
   # user hasn't tried to set a true value at an inappropriate time.

   if (&value_is_true($nl->get_value($var)) && $cannot_be_true) {
      $log->fatal_error($cannot_be_true);
   }

   # if do_transient_crops is .true. and any of these (n_dom_* or toosmall_*)
   # are > 0 or collapse_urban = .true., then give fatal error
   if (&value_is_true($nl->get_value($var))) {
      if (&value_is_true($nl->get_value('collapse_urban'))) {
         $log->fatal_error("$var cannot be combined with collapse_urban");
      }
      if ($n_dom_pfts > 0 || $n_dom_landunits > 0 || $toosmall_soil > 0 || $toosmall_crop > 0 || $toosmall_glacier > 0 || $toosmall_lake > 0 || $toosmall_wetland > 0 || $toosmall_urban > 0) {
         $log->fatal_error("$var cannot be combined with any of the of the following > 0: n_dom_pfts > 0, n_dom_landunit > 0, toosmall_soil > 0._r8, toosmall_crop > 0._r8, toosmall_glacier > 0._r8, toosmall_lake > 0._r8, toosmall_wetland > 0._r8, toosmall_urban > 0._r8");
      }
   }

   my $dopft = "do_transient_pfts";
   # Make sure the value agrees with the do_transient_pft flag
   if ( (  &value_is_true($nl->get_value($var))) && (! &value_is_true($nl->get_value($dopft))) ||
        (! &value_is_true($nl->get_value($var))) && (  &value_is_true($nl->get_value($dopft))) ) {
      $log->fatal_error("$var and $dopft do NOT agree and need to");
   }
}

sub setup_logic_do_transient_lakes {
   #
   # Set do_transient_lakes default value, and perform error checking on do_transient_lakes
   #
   # Assumes the following are already set in the namelist (although it's okay
   # for them to be unset if that will be their final state):
   # - flanduse_timeseries
   #
   # NOTE(wjs, 2020-08-23) I based this function on setup_logic_do_transient_crops. I'm
   # not sure if all of the checks here are truly important for transient lakes (in
   # particular, my guess is that collapse_urban could probably be done with transient
   # lakes - as well as transient pfts and transient crops for that matter), but some of
   # the checks probably are needed, and it seems best to keep transient lakes consistent
   # with other transient areas in this respect.
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   my $var = 'do_transient_lakes';

   # Start by assuming a default value of '.true.'. Then check a number of
   # conditions under which do_transient_lakes cannot be true. Under these
   # conditions: (1) set default value to '.false.'; (2) make sure that the
   # value is indeed false (e.g., that the user didn't try to set it to true).

   my $default_val = ".true.";

   # cannot_be_true will be set to a non-empty string in any case where
   # do_transient_lakes should not be true; if it turns out that
   # do_transient_lakes IS true in any of these cases, a fatal error will be
   # generated
   my $cannot_be_true = "";

   my $n_dom_pfts = $nl->get_value( 'n_dom_pfts' );
   my $n_dom_landunits = $nl->get_value( 'n_dom_landunits' );
   my $toosmall_soil = $nl->get_value( 'toosmall_soil' );
   my $toosmall_crop = $nl->get_value( 'toosmall_crop' );
   my $toosmall_glacier = $nl->get_value( 'toosmall_glacier' );
   my $toosmall_lake = $nl->get_value( 'toosmall_lake' );
   my $toosmall_wetland = $nl->get_value( 'toosmall_wetland' );
   my $toosmall_urban = $nl->get_value( 'toosmall_urban' );

   if (string_is_undef_or_empty($nl->get_value('flanduse_timeseries'))) {
      $cannot_be_true = "$var can only be set to true when running a transient case (flanduse_timeseries non-blank)";
   }

   if (!$cannot_be_true) {
      # Note that, if the variable cannot be true, we don't call add_default
      # - so that we don't clutter up the namelist with variables that don't
      # matter for this case
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, val=>$default_val);
   }

   # Make sure the value is false when it needs to be false - i.e., that the
   # user hasn't tried to set a true value at an inappropriate time.

   if (&value_is_true($nl->get_value($var)) && $cannot_be_true) {
      $log->fatal_error($cannot_be_true);
   }

   # if do_transient_lakes is .true. and any of these (n_dom_* or toosmall_*)
   # are > 0 or collapse_urban = .true., then give fatal error
   if (&value_is_true($nl->get_value($var))) {
      if (&value_is_true($nl->get_value('collapse_urban'))) {
         $log->fatal_error("$var cannot be combined with collapse_urban");
      } elsif (&value_is_true($nl->get_value('use_hillslope'))) {
         $log->fatal_error("$var cannot be combined with use_hillslope");
      }
      if ($n_dom_pfts > 0 || $n_dom_landunits > 0 || $toosmall_soil > 0 || $toosmall_crop > 0 || $toosmall_glacier > 0 || $toosmall_lake > 0 || $toosmall_wetland > 0 || $toosmall_urban > 0) {
         $log->fatal_error("$var cannot be combined with any of the of the following > 0: n_dom_pfts > 0, n_dom_landunit > 0, toosmall_soil > 0._r8, toosmall_crop > 0._r8, toosmall_glacier > 0._r8, toosmall_lake > 0._r8, toosmall_wetland > 0._r8, toosmall_urban > 0._r8");
      }
   }
}

sub setup_logic_do_transient_urban {
   #
   # Set do_transient_urban default value, and perform error checking on do_transient_urban
   #
   # Assumes the following are already set in the namelist (although it's okay
   # for them to be unset if that will be their final state):
   # - flanduse_timeseries
   #
   # NOTE(kwo, 2021-08-11) I based this function on setup_logic_do_transient_lakes.
   # As in NOTE(wjs, 2020-08-23) I'm not sure if all of the checks here are truly important
   # for transient urban (in particular, my guess is that collapse_urban could probably be done with transient
   # urban - as well as transient pfts and transient crops for that matter), but some of
   # the checks probably are needed, and it seems best to keep transient urban consistent
   # with other transient areas in this respect.
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   my $var = 'do_transient_urban';

   # Start by assuming a default value of '.true.'. Then check a number of
   # conditions under which do_transient_urban cannot be true. Under these
   # conditions: (1) set default value to '.false.'; (2) make sure that the
   # value is indeed false (e.g., that the user didn't try to set it to true).

   my $default_val = ".true.";

   # cannot_be_true will be set to a non-empty string in any case where
   # do_transient_urban should not be true; if it turns out that
   # do_transient_urban IS true in any of these cases, a fatal error will be
   # generated
   my $cannot_be_true = "";

   my $n_dom_pfts = $nl->get_value( 'n_dom_pfts' );
   my $n_dom_landunits = $nl->get_value( 'n_dom_landunits' );
   my $toosmall_soil = $nl->get_value( 'toosmall_soil' );
   my $toosmall_crop = $nl->get_value( 'toosmall_crop' );
   my $toosmall_glacier = $nl->get_value( 'toosmall_glacier' );
   my $toosmall_lake = $nl->get_value( 'toosmall_lake' );
   my $toosmall_wetland = $nl->get_value( 'toosmall_wetland' );
   my $toosmall_urban = $nl->get_value( 'toosmall_urban' );

   if (string_is_undef_or_empty($nl->get_value('flanduse_timeseries'))) {
      $cannot_be_true = "$var can only be set to true when running a transient case (flanduse_timeseries non-blank)";
   }

   if (!$cannot_be_true) {
      # Note that, if the variable cannot be true, we don't call add_default
      # - so that we don't clutter up the namelist with variables that don't
      # matter for this case
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, val=>$default_val);
   }

   # Make sure the value is false when it needs to be false - i.e., that the
   # user hasn't tried to set a true value at an inappropriate time.

   if (&value_is_true($nl->get_value($var)) && $cannot_be_true) {
      $log->fatal_error($cannot_be_true);
   }

   # if do_transient_urban is .true. and any of these (n_dom_* or toosmall_*)
   # are > 0 or collapse_urban = .true., then give fatal error
   if (&value_is_true($nl->get_value($var))) {
      if (&value_is_true($nl->get_value('collapse_urban'))) {
         $log->fatal_error("$var cannot be combined with collapse_urban");
      } elsif (&value_is_true($nl->get_value('use_hillslope'))) {
         $log->fatal_error("$var cannot be combined with use_hillslope");
      }
      if ($n_dom_pfts > 0 || $n_dom_landunits > 0 || $toosmall_soil > 0 || $toosmall_crop > 0 || $toosmall_glacier > 0 || $toosmall_lake > 0 || $toosmall_wetland > 0 || $toosmall_urban > 0) {
         $log->fatal_error("$var cannot be combined with any of the of the following > 0: n_dom_pfts > 0, n_dom_landunit > 0, toosmall_soil > 0._r8, toosmall_crop > 0._r8, toosmall_glacier > 0._r8, toosmall_lake > 0._r8, toosmall_wetland > 0._r8, toosmall_urban > 0._r8");
      }
   }
}

sub setup_logic_do_harvest {
   #
   # Set do_harvest default value, and perform error checking on do_harvest
   #
   # Assumes the following are already set in the namelist (although it's okay
   # for them to be unset if that will be their final state):
   # - flanduse_timeseries
   # - use_cn
   # - use_fates
   #
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   my $var = 'do_harvest';

   # Start by assuming a default value of '.true.'. Then check a number of
   # conditions under which do_harvest cannot be true. Under these
   # conditions: (1) set default value to '.false.'; (2) make sure that the
   # value is indeed false (e.g., that the user didn't try to set it to true).

   my $default_val = ".true.";

   # cannot_be_true will be set to a non-empty string in any case where
   # do_harvest should not be true; if it turns out that do_harvest IS true
   # in any of these cases, a fatal error will be generated
   my $cannot_be_true = "";

      if (string_is_undef_or_empty($nl->get_value('flanduse_timeseries'))) {
         $cannot_be_true = "$var can only be set to true when running a transient case (flanduse_timeseries non-blank)";
      }

      elsif (!&value_is_true($nl->get_value('use_cn')) && !&value_is_true($nl->get_value('use_fates'))) {
         $cannot_be_true = "$var can only be set to true when running with either CN or FATES";
      }

      if ($cannot_be_true) {
         $default_val = ".false.";
      }

   if (!$cannot_be_true) {
      # Note that, if the variable cannot be true, we don't call add_default
      # - so that we don't clutter up the namelist with variables that don't
      # matter for this case
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, val=>$default_val);
   }

   # Make sure the value is false when it needs to be false - i.e., that the
   # user hasn't tried to set a true value at an inappropriate time.

   if (&value_is_true($nl->get_value($var)) && $cannot_be_true) {
      $log->fatal_error($cannot_be_true);
   }
}

#-------------------------------------------------------------------------------

sub setup_logic_do_grossunrep {
   #
   # Set do_grossunrep default value, and perform error checking on do_grossunrep
   #
   # Assumes the following are already set in the namelist (although it's okay
   # for them to be unset if that will be their final state):
   # - flanduse_timeseries
   # - use_cn
   # - use_fates
   #
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   my $var = 'do_grossunrep';

   # Start by assuming a default value of '.true.'. Then check a number of
   # conditions under which do_grossunrep cannot be true. Under these
   # conditions: (1) set default value to '.false.'; (2) make sure that the
   # value is indeed false (e.g., that the user didn't try to set it to true).

   my $default_val = ".false.";

   # cannot_be_true will be set to a non-empty string in any case where
   # do_grossunrep should not be true; if it turns out that do_grossunrep IS true
   # in any of these cases, a fatal error will be generated
   my $cannot_be_true = "";

   if (string_is_undef_or_empty($nl->get_value('flanduse_timeseries'))) {
      $cannot_be_true = "$var can only be set to true when running a transient case (flanduse_timeseries non-blank)";
   }
   elsif (&value_is_true($nl->get_value('use_fates'))) {
      $cannot_be_true = "$var currently doesn't work with FATES";
   }
   elsif (!&value_is_true($nl->get_value('use_cn'))) {
      $cannot_be_true = "$var can only be set to true when running with CN (use_cn = true)";
   }

   if ($cannot_be_true) {
      $default_val = ".false.";
   }

   if (!$cannot_be_true) {
      # Note that, if the variable cannot be true, we don't call add_default
      # - so that we don't clutter up the namelist with variables that don't
      # matter for this case
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, val=>$default_val);
   }

   # Make sure the value is false when it needs to be false - i.e., that the
   # user hasn't tried to set a true value at an inappropriate time.

   if (&value_is_true($nl->get_value($var)) && $cannot_be_true) {
      $log->fatal_error($cannot_be_true);
   }

}

#-------------------------------------------------------------------------------
sub setup_logic_spinup {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( $nl_flags->{'bgc_mode'} eq "sp" && defined($nl->get_value('override_bgc_restart_mismatch_dump'))) {
     $log->fatal_error("CN must be on if override_bgc_restart_mismatch_dump is set.");
  }
  if ( $nl_flags->{'clm_accelerated_spinup'} eq "on" ) {
     foreach my $var ( "hist_nhtfrq", "hist_fincl1", "hist_empty_htapes", "hist_mfilt" ) {
         add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
                     $var, use_cn=>$nl_flags->{'use_cn'}, use_fates=>$nl_flags->{'use_fates'},
                     use_cndv=>$nl_flags->{'use_cndv'} );
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_bgc_shared {
  my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

  if ( $nl_flags->{'bgc_mode'} ne "sp" ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'constrain_stress_deciduous_onset', 'phys'=>$physv->as_string() );
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_cnphenology {
  my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

  my @list  = (  "onset_thresh_depends_on_veg", "min_critical_dayl_method" );
  foreach my $var ( @list ) {
    if (  &value_is_true($nl_flags->{'use_cn'}) ) {
       add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                   'phys'=>$physv->as_string(), 'use_cn'=>$nl_flags->{'use_cn'} );
    } else {
       if ( defined($nl->get_value($var)) ) {
          $log->fatal_error("$var should only be set if use_cn is on");
       }
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_supplemental_nitrogen {
  #
  # Supplemental Nitrogen for prognostic crop cases
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( $nl_flags->{'bgc_mode'} ne "sp" && $nl_flags->{'bgc_mode'} ne "fates" && &value_is_true($nl_flags->{'use_crop'}) ) {
      # If this is non-fates, non-sp and crop is active
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
		  'suplnitro', 'use_cn'=>$nl_flags->{'use_cn'}, 'use_crop'=>$nl_flags->{'use_crop'});

  } elsif ( $nl_flags->{'bgc_mode'} eq "fates" && not &value_is_true( $nl_flags->{'use_fates_sp'})  ) {
      # Or... if its fates but not fates-sp
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
		  'suplnitro', 'use_fates'=>$nl_flags->{'use_fates'});
  }
  
  #
  # Error checking for suplnitro
  #
  my $suplnitro = $nl->get_value('suplnitro');
  if ( defined($suplnitro) ) {
    if ( $nl_flags->{'bgc_mode'} eq "sp" ) {
      $log->fatal_error("supplemental Nitrogen (suplnitro) is set, but neither CN nor CNDV nor FATES is active!");
    }
    if ( ! &value_is_true($nl_flags->{'use_crop'}) && $suplnitro =~ /PROG_CROP_ONLY/i ) {
      $log->fatal_error("supplemental Nitrogen is set to run over prognostic crops, but prognostic crop is NOT active!");
    }

    if ( $suplnitro =~ /ALL/i ) {
      if ( $nl_flags->{'bgc_spinup'} eq "on" && $nl_flags->{'bgc_mode'} ne "fates" ) {
        $log->warning("There is no need to use a bgc_spinup mode when supplemental Nitrogen is on for all PFT's, as these modes spinup Nitrogen" );
      }
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_hydrology_params {
  #
  # Logic for hydrology parameters
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $lower = $nl->get_value( 'lower_boundary_condition'  );
  my $var   = "baseflow_scalar";
  if ( $lower == 1 || $lower == 2 ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
                 $var, 'lower_boundary_condition' => $lower );
  }
  my $val   = $nl->get_value( $var );
  if ( defined($val) ) {
     if ( $lower != 1 && $lower != 2 ) {
        $log->fatal_error("baseflow_scalar is only used for lower_boundary_condition of flux or zero-flux");
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_irrigation_parameters {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $var;
  foreach $var ("irrig_min_lai", "irrig_start_time", "irrig_length",
                "irrig_target_smp", "irrig_depth", "irrig_threshold_fraction",
                "limit_irrigation_if_rof_enabled","use_groundwater_irrigation",
                "irrig_method_default") {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var);
  }

  if ( &value_is_true($nl->get_value('use_groundwater_irrigation')) &&
       ! &value_is_true($nl->get_value('limit_irrigation_if_rof_enabled'))) {
     $log->fatal_error("use_groundwater_irrigation only makes sense if limit_irrigation_if_rof_enabled is set. (If limit_irrigation_if_rof_enabled is .false., then groundwater extraction will never be invoked.)")
  }

  my $lower = $nl->get_value( 'lower_boundary_condition' );
  if ( ($lower == 3 || $lower == 4) && (&value_is_true($nl->get_value( 'use_groundwater_irrigation' ))) ) {
     $log->fatal_error("use_groundwater_irrigation can only be used when lower_boundary_condition is NOT 3 or 4");
  }

  $var = "irrig_river_volume_threshold";
  if ( &value_is_true($nl->get_value("limit_irrigation_if_rof_enabled")) ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var);
  } else {
     if (defined($nl->get_value($var))) {
        $log->fatal_error("$var can only be set if limit_irrigation_if_rof_enabled is true");
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_surfacealbedo {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snowveg_affects_radiation' );
}

#-------------------------------------------------------------------------------

sub setup_logic_water_tracers {
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   my $var;
   foreach $var ("enable_water_tracer_consistency_checks",
                 "enable_water_isotopes") {
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var);
   }
}

#-------------------------------------------------------------------------------

sub setup_logic_nitrif_params {
  #
  # Logic for nitrification parameters
  #
  my ($nl_flags, $definition, $defaults, $nl) = @_;

  if ( !  &value_is_true($nl_flags->{'use_nitrif_denitrif'}) ) {
    my @vars = ( "k_nitr_max", "denitrif_respiration_coefficient", "denitrif_respiration_exponent");
    foreach my $var ( @vars ) {
       if ( defined($nl->get_value( $var ) ) ) {
         $log->fatal_error("$var is only used when use_nitrif_denitrif is turned on");
       }
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_mineral_nitrogen_dynamics {
  #
  # Logic for mineral_nitrogen_dynamics
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my @vars = ( "freelivfix_slope_wet", "freelivfix_intercept" );
  if (  &value_is_true($nl_flags->{'use_cn'}) && &value_is_true($nl->get_value('use_fun')) ) {
    foreach my $var ( @vars ) {
       add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                'use_cn'=>$nl_flags->{'use_cn'}, 'use_fun'=>$nl->get_value('use_fun') );
    }
  } else {
    foreach my $var ( @vars ) {
       if ( defined($nl->get_value( $var ) ) ) {
         $log->fatal_error("$var is only used when use_cn and use_fun are both turned on");
       }
    }
  }
}


#-------------------------------------------------------------------------------

sub setup_logic_hydrology_switches {
  #
  # Check on Switches for hydrology
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_subgrid_fluxes');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snow_cover_fraction_method');
  my $subgrid    = $nl->get_value('use_subgrid_fluxes' );
  my $h2osfcflag = $nl->get_value('h2osfcflag'  );
  my $scf_method = $nl->get_value('snow_cover_fraction_method');
  if ( $h2osfcflag == 1 && ! &value_is_true($subgrid) ) {
    $log->fatal_error("if h2osfcflag is ON, use_subgrid_fluxes can NOT be off!");
  }
  if ( remove_leading_and_trailing_quotes($scf_method) eq 'NiuYang2007' && &value_is_true($subgrid) ) {
     $log->fatal_error("snow_cover_fraction_method NiuYang2007 is incompatible with use_subgrid_fluxes");
  }
  # Test bad configurations
  my $lower   = $nl->get_value( 'lower_boundary_condition'  );
  my $use_vic = $nl_flags->{'use_vichydro'};
  my $use_bed = $nl->get_value( 'use_bedrock'               );
  my $soilmtd = $nl->get_value( 'soilwater_movement_method' );
  if ( defined($soilmtd) && defined($lower) && $soilmtd == 0 && $lower != 4 ) {
     $log->fatal_error( "If soil water movement method is zeng-decker -- lower_boundary_condition can only be aquifer" );
  }
  if ( defined($soilmtd) && defined($lower) && $soilmtd == 1 && $lower == 4 ) {
     $log->fatal_error( "If soil water movement method is adaptive -- lower_boundary_condition can NOT be aquifer" );
  }
  if ( defined($use_bed) && defined($lower) && (&value_is_true($use_bed)) && $lower != 2 ) {
     $log->fatal_error( "If use_bedrock is on -- lower_boundary_condition can only be flux" );
  }
  if ( defined($use_vic) && defined($lower) && (&value_is_true($use_vic)) && $lower != 3 && $lower != 4) {
     $log->fatal_error( "If use_vichydro is on -- lower_boundary_condition can only be table or aquifer" );
  }
  if ( defined($h2osfcflag) && defined($lower) && $h2osfcflag == 0 && $lower != 4 ) {
     $log->fatal_error( "If h2osfcflag is 0 lower_boundary_condition can only be aquifer" );
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_methane {
  #
  # CH4 model if bgc=CN or CNDV
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( &value_is_true($nl_flags->{'use_lch4'}) ) {
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'finundation_method',
                'use_cn'=>$nl_flags->{'use_cn'}, 'use_fates'=>$nl_flags->{'use_fates'} );
    my $finundation_method = remove_leading_and_trailing_quotes($nl->get_value('finundation_method' ));
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_ch4finundated',
             'finundation_method'=>$finundation_method);
    if ($opts->{'driver'} eq "nuopc" ) {
        add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_meshfile_ch4finundated',
                    'finundation_method'=>$finundation_method);
    }
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_aereoxid_prog',
                'use_cn'=>$nl_flags->{'use_cn'}, 'use_fates'=>$nl_flags->{'use_fates'} );
    #
    # Check if use_aereoxid_prog is set.  If no, then read value of aereoxid from
    # parameters file
    #
    my $use_aereoxid_prog = $nl->get_value('use_aereoxid_prog');
    if ( defined($use_aereoxid_prog) && ! &value_is_true($use_aereoxid_prog) ) {
      $log->warning("Using aereoxid value from parameters file." );
    }
  } else {
    my @vars = $nl->get_variable_names('ch4par_in');
    if ( $#vars >= 0 ) {
      $log->fatal_error("ch4par_in namelist variables were set, but Methane model NOT defined in the configuration (use_lch4)");
    }
  }

  #
  # Ch4 namelist checking
  #
  if ( &value_is_true($nl_flags->{'use_lch4'}) ) {
    my $allowlakeprod = $nl->get_value('allowlakeprod');
    if ( ! defined($allowlakeprod) ||
         (defined($allowlakeprod) && ! &value_is_true($allowlakeprod)) ) {
      if ( defined($nl->get_value('lake_decomp_fact')) ) {
        $log->fatal_error("lake_decomp_fact set without allowlakeprod=.true.");
      }
    }
    my $pftspec_rootprof = $nl->get_value('pftspecific_rootingprofile');
    if ( ! defined($pftspec_rootprof) ||
         (defined($pftspec_rootprof) && &value_is_true($pftspec_rootprof) ) ) {
      if ( defined($nl->get_value('rootprof_exp')) ) {
        $log->fatal_error("rootprof_exp set without pftspecific_rootingprofile=.false.");
      }
    }
  } else {
    my @vars = ( "allowlakeprod", "anoxia", "pftspecific_rootingprofile" );
    foreach my $var ( @vars ) {
      if ( defined($nl->get_value($var)) ) {
        $log->fatal_error("$var set without methane model configuration on (use_lch4)");
      }
    }
    my $var = "use_nitrif_denitrif";
    if ( (! &value_is_true( $nl_flags->{'use_fates'} ) ) && &value_is_true($nl->get_value($var)) ) {
       $log->warning("methane is off (use_lch4=FALSE), but $var is TRUE, both need to be on, unless FATES is also on" );
    }
  }
} # end methane

#-------------------------------------------------------------------------------

sub setup_logic_dynamic_plant_nitrogen_alloc {
  #
  # dynamic plant nitrogen allocation model, bgc=bgc
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

  if ( &value_is_true($nl_flags->{'use_cn'}) ) {
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_flexibleCN',
                'phys'=>$physv->as_string(), 'use_cn'=>$nl_flags->{'use_cn'} );
    $nl_flags->{'use_flexibleCN'} = $nl->get_value('use_flexibleCN');
    if ( &value_is_true($nl->get_value('use_fun') ) && not &value_is_true( $nl_flags->{'use_flexibleCN'}) ) {
       $log->warning("FUN has NOT been extensively tested without use_flexibleCN on, so could result in failures or unexpected results" );
    }

    if ( &value_is_true($nl_flags->{'use_flexibleCN'}) ) {
      # TODO(bja, 2015-04) make this depend on > clm 5.0 and bgc mode at some point.
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'MM_Nuptake_opt',
                  'use_flexibleCN'=>$nl_flags->{'use_flexibleCN'} );
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'CNratio_floating',
                  'use_flexibleCN'=>$nl_flags->{'use_flexibleCN'} );
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'reduce_dayl_factor',
                  'use_flexibleCN'=>$nl_flags->{'use_flexibleCN'} );
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'vcmax_opt',
                  'use_flexibleCN'=>$nl_flags->{'use_flexibleCN'} );
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'CN_evergreen_phenology_opt',
                  'use_flexibleCN'=>$nl_flags->{'use_flexibleCN'} );
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'carbon_resp_opt',
                  'use_flexibleCN'=>$nl_flags->{'use_flexibleCN'}, 'use_fun'=>$nl->get_value('use_fun') );
      if ( $nl->get_value('carbon_resp_opt') == 1 && &value_is_true($nl->get_value('use_fun')) ) {
        $log->fatal_error("carbon_resp_opt should NOT be set to 1 when FUN is also on");
      }
    }
  } else {
     if ( &value_is_true($nl->get_value('use_flexibleCN')) ) {
        $log->fatal_error("use_flexibleCN can ONLY be set if CN is on");
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_o3_veg_stress_method {
  #
  # Ozone vegetation stress method
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

  my $var = 'o3_veg_stress_method';

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var );

  my $val = $nl->get_value($var);

  if (remove_leading_and_trailing_quotes($val) eq "stress_falk" && not (&value_is_true($nl_flags->{'use_luna'})) ) {
    $log->fatal_error(" use_luna=.true. is required for $var='stress_falk'.");
  }

}

#-------------------------------------------------------------------------------

sub setup_logic_luna {
  #
  # LUNA model to calculate photosynthetic capacities based on environmental conditions
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_luna',
              'phys'=>$physv->as_string(), 'use_cn'=>$nl_flags->{'use_cn'}, 'use_fates'=>$nl_flags->{'use_fates'},
              'use_nitrif_denitrif'=>$nl_flags->{'use_nitrif_denitrif'} );

  if ( &value_is_true( $nl_flags->{'use_cn'} ) ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_nguardrail',
                 'use_cn'=>$nl_flags->{'use_cn'} );
  }
  $nl_flags->{'use_luna'} = $nl->get_value('use_luna');

  # LUNA can NOT be on with FATES
  if ( &value_is_true( $nl_flags->{'use_luna'} ) && &value_is_true( $nl_flags->{'use_fates'} )) {
     $log->fatal_error("Cannot turn use_luna to true when bgc=fates" );
  }

  my $vcmax_opt= $nl->get_value('vcmax_opt');
  # lnc_opt only applies if luna is on or for vcmax_opt=3/4
  if ( &value_is_true( $nl_flags->{'use_luna'} ) || $vcmax_opt == 3 || $vcmax_opt == 4 ) {
     # lnc_opt can be set for both CN on and off
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'lnc_opt',
                 'use_cn'=>$nl_flags->{'use_cn'} );
  }
  if ( &value_is_true($nl->get_value('lnc_opt') ) && not &value_is_true( $nl_flags->{'use_cn'}) ) {
     $log->fatal_error("Cannot turn lnc_opt to true when bgc=sp" );
  }
  my $var = "jmaxb1";
  if ( &value_is_true( $nl_flags->{'use_luna'} ) ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var,
                 'use_luna'=>$nl_flags->{'use_luna'} );
  }
  my $val = $nl->get_value($var);
  if ( ! &value_is_true( $nl_flags->{'use_luna'} ) ) {
     if ( defined($val) ) {
        $log->fatal_error("Cannot set $var when use_luna is NOT on" );
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_hillslope {
  #
  # Hillslope model
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_hillslope' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'downscale_hillslope_meteorology' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'hillslope_head_gradient_method' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'hillslope_transmissivity_method' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'hillslope_pft_distribution_method' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'hillslope_soil_profile_method' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_hillslope_routing', 'use_hillslope'=>$nl_flags->{'use_hillslope'} );
  my $use_hillslope = $nl->get_value('use_hillslope');
  my $use_hillslope_routing = $nl->get_value('use_hillslope_routing');
  if ( (! &value_is_true($use_hillslope)) && &value_is_true($use_hillslope_routing) ) {
      $log->fatal_error("Cannot turn on use_hillslope_routing when use_hillslope is off\n" );
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_hydrstress {
  #
  # Plant hydraulic stress model
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_hydrstress',
              'configuration'=>$nl_flags->{'configuration'}, 'use_fates'=>$nl_flags->{'use_fates'} );
  $nl_flags->{'use_hydrstress'} = $nl->get_value('use_hydrstress');
  if ( &value_is_true( $nl_flags->{'use_fates'} ) && &value_is_true( $nl_flags->{'use_hydrstress'} ) ) {
     $log->fatal_error("Cannot turn use_hydrstress on when use_fates is on" );
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_fertilizer {
  #
  # Flags to control fertilizer application
  #
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_fertilizer',
               'use_crop'=>$nl_flags->{'use_crop'} );
   my $use_fert = $nl->get_value('use_fertilizer');
   if ( (! &value_is_true($nl_flags->{'use_crop'})) && &value_is_true($use_fert) ) {
      $log->fatal_error("use_ferilizer can NOT be on without prognostic crop\n" );
   }
}

#-------------------------------------------------------------------------------

sub setup_logic_grainproduct {
  #
  # Flags to control 1-year grain product pool
  #
   my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_grainproduct',
               'use_crop'=>$nl_flags->{'use_crop'}, 'phys'=>$physv->as_string() );
   if ( (! &value_is_true($nl_flags->{'use_crop'})) && &value_is_true($nl->get_value('use_grainproduct') ) ) {
      $log->fatal_error("use_grainproduct can NOT be on without prognostic crop\n" );
   }
}

#-------------------------------------------------------------------------------

sub setup_logic_dynamic_roots {
  #
  # dynamic root model
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_dynroot', 'phys'=>$physv->as_string(), 'bgc_mode'=>$nl_flags->{'bgc_mode'});
  my $use_dynroot = $nl->get_value('use_dynroot');
  if ( &value_is_true($use_dynroot) && ($nl_flags->{'bgc_mode'} eq "sp") ) {
     $log->fatal_error("Cannot turn dynroot mode on mode bgc=sp\n" .
                       "Set the bgc mode to 'bgc'.");
  }
  if ( &value_is_true( $use_dynroot ) && &value_is_true( $nl_flags->{'use_hydrstress'} ) ) {
     $log->fatal_error("Cannot turn use_dynroot on when use_hydrstress is on" );
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_c_isotope {
  #
  # Error checking for C-isotope options
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $use_c13 = $nl->get_value('use_c13');
  my $use_c14 = $nl->get_value('use_c14');
  if ( $nl_flags->{'bgc_mode'} ne "sp" && $nl_flags->{'bgc_mode'} ne "fates" ) {
    if ( $nl_flags->{'bgc_mode'} ne "bgc" ) {
      if ( defined($use_c13) && &value_is_true($use_c13) ) {
        $log->warning("use_c13 is ONLY scientifically validated with the bgc=BGC configuration" );
      }
      if ( defined($use_c14) && &value_is_true($use_c14) ) {
        $log->warning("use_c14 is ONLY scientifically validated with the bgc=BGC configuration" );
      }
    }
    if ( defined($use_c14) ) {
      if ( &value_is_true($use_c14) ) {
        my $use_c14_bombspike = $nl->get_value('use_c14_bombspike');
        if ( defined($use_c14_bombspike) && &value_is_true($use_c14_bombspike) ) {
           add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'atm_c14_filename',
                   'use_c14'=>$use_c14, 'use_cn'=>$nl_flags->{'use_cn'}, 'use_c14_bombspike'=>$nl->get_value('use_c14_bombspike'),
                   'ssp_rcp'=>$nl_flags->{'ssp_rcp'} );
        }
      } else {
        if ( defined($nl->get_value('use_c14_bombspike')) ||
             defined($nl->get_value('atm_c14_filename')) ) {
          $log->fatal_error("use_c14 is FALSE and use_c14_bombspike or atm_c14_filename set");
        }
      }
    } else {
      if ( defined($nl->get_value('use_c14_bombspike')) ||
           defined($nl->get_value('atm_c14_filename')) ) {
        $log->fatal_error("use_c14 NOT set to .true., but use_c14_bompspike/atm_c14_filename defined.");
      }
    }
    if ( defined($use_c13) ) {
      if ( &value_is_true($use_c13) ) {
        my $use_c13_timeseries = $nl->get_value('use_c13_timeseries');
        if ( defined($use_c13_timeseries) && &value_is_true($use_c13_timeseries) ) {
           add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'atm_c13_filename',
                   'use_c13'=>$use_c13, 'use_cn'=>$nl_flags->{'use_cn'}, 'use_c13_timeseries'=>$nl->get_value('use_c13_timeseries'),
                   'ssp_rcp'=>$nl_flags->{'ssp_rcp'} );
        }
      } else {
        if ( defined($nl->get_value('use_c13_timeseries')) ||
             defined($nl->get_value('atm_c13_filename')) ) {
          $log->fatal_error("use_c13 is FALSE and use_c13_timeseries or atm_c13_filename set");
        }
      }
    } else {
      if ( defined($nl->get_value('use_c13_timeseries')) ||
           defined($nl->get_value('atm_c13_filename')) ) {
        $log->fatal_error("use_c13 NOT set to .true., but use_c13_bompspike/atm_c13_filename defined.");
      }
    }
  } else {
    if ( defined($use_c13) ||
         defined($use_c14) ||
         defined($nl->get_value('use_c14_bombspike')) ||
         defined($nl->get_value('atm_c14_filename'))  ||
         defined($nl->get_value('use_c13_timeseries')) ||
         defined($nl->get_value('atm_c13_filename')) ) {
           $log->fatal_error("bgc=sp and C isotope  namelist variables were set, both can't be used at the same time");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_nitrogen_deposition {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  #
  # Nitrogen deposition for bgc=CN or fates
  #
  if ( ($nl_flags->{'bgc_mode'} =~/bgc/) ) {   # or  ($nl_flags->{'bgc_mode'} =~/fates/) ) {
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'ndepmapalgo', 'phys'=>$nl_flags->{'phys'},
                'use_cn'=>$nl_flags->{'use_cn'}, 'hgrid'=>$nl_flags->{'res'},
                'clm_accelerated_spinup'=>$nl_flags->{'clm_accelerated_spinup'} );
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'ndep_taxmode', 'phys'=>$nl_flags->{'phys'},
		'use_cn'=>$nl_flags->{'use_cn'}, 
		'lnd_tuning_mode'=>$nl_flags->{'lnd_tuning_mode'} );
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'ndep_varlist', 'phys'=>$nl_flags->{'phys'},
		'use_cn'=>$nl_flags->{'use_cn'}, 
		'lnd_tuning_mode'=>$nl_flags->{'lnd_tuning_mode'} );
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_first_ndep', 'phys'=>$nl_flags->{'phys'},
                'use_cn'=>$nl_flags->{'use_cn'}, 'sim_year'=>$nl_flags->{'sim_year'},
                'sim_year_range'=>$nl_flags->{'sim_year_range'});
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_last_ndep', 'phys'=>$nl_flags->{'phys'},
                'use_cn'=>$nl_flags->{'use_cn'}, 'sim_year'=>$nl_flags->{'sim_year'},
                'sim_year_range'=>$nl_flags->{'sim_year_range'});
    # Set align year, if first and last years are different
    if ( $nl->get_value('stream_year_first_ndep') != $nl->get_value('stream_year_last_ndep') ) {
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'model_year_align_ndep', 'sim_year'=>$nl_flags->{'sim_year'},
                  'sim_year_range'=>$nl_flags->{'sim_year_range'});
    }
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_ndep', 'phys'=>$nl_flags->{'phys'},
                'use_cn'=>$nl_flags->{'use_cn'}, 'lnd_tuning_mode'=>$nl_flags->{'lnd_tuning_mode'},
                'hgrid'=>"0.9x1.25", 'ssp_rcp'=>$nl_flags->{'ssp_rcp'}, 'nofail'=>1 );
    if ( ! defined($nl->get_value('stream_fldfilename_ndep') ) ) {
        # Also check at f19 resolution
        add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_ndep', 'phys'=>$nl_flags->{'phys'},
                    'use_cn'=>$nl_flags->{'use_cn'}, 'lnd_tuning_mode'=>$nl_flags->{'lnd_tuning_mode'},
                    'hgrid'=>"1.9x2.5", 'ssp_rcp'=>$nl_flags->{'ssp_rcp'}, 'nofail'=>1 );
        # If not found report an error
        if ( ! defined($nl->get_value('stream_fldfilename_ndep') ) ) {
            $log->warning("Did NOT find the Nitrogen-deposition forcing file (stream_fldfilename_ndep) for this ssp_rcp\n" .
                          "One way to get around this is to point to a file for another existing ssp_rcp in your user_nl_clm file.\n" .
                          "If you are running with CAM and WACCM chemistry Nitrogen deposition will come through the coupler.\n" .
                          "This file won't be used, so it doesn't matter what it points to -- but it's required to point to something.\n" )
        }
    }
    if ($opts->{'driver'} eq "nuopc" ) {
        add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_meshfile_ndep', 'phys'=>$nl_flags->{'phys'},
                    'use_cn'=>$nl_flags->{'use_cn'}, 'lnd_tuning_mode'=>$nl_flags->{'lnd_tuning_mode'},
                    'hgrid'=>"0.9x1.25", 'ssp_rcp'=>$nl_flags->{'ssp_rcp'}, 'nofail'=>1 );
        if ( ! defined($nl->get_value('stream_fldfilename_ndep') ) ) {
            # Also check at f19 resolution
            add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_meshfile_ndep', 'phys'=>$nl_flags->{'phys'},
                        'use_cn'=>$nl_flags->{'use_cn'}, 'lnd_tuning_mode'=>$nl_flags->{'lnd_tuning_mode'},
                        'hgrid'=>"1.9x2.5", 'ssp_rcp'=>$nl_flags->{'ssp_rcp'}, 'nofail'=>1 );
            # If not found report an error
            if ( ! defined($nl->get_value('stream_meshfile_ndep') ) ) {
                $log->warning("Did NOT find the Nitrogen-deposition meshfile file (stream_meshfilee_ndep) for this ssp_rcp. \n")
            }
        }
    }
  } else {
    # If bgc is NOT CN/CNDV then make sure none of the ndep settings are set!
    if ( defined($nl->get_value('stream_year_first_ndep')) ||
         defined($nl->get_value('stream_year_last_ndep'))  ||
         defined($nl->get_value('model_year_align_ndep'))  ||
         defined($nl->get_value('ndep_taxmode'         ))  ||
         defined($nl->get_value('ndep_varlist'         ))  ||
         defined($nl->get_value('stream_fldfilename_ndep'))
       ) {
      $log->fatal_error("When bgc is NOT CN or CNDV none of: stream_year_first_ndep," .
                  "stream_year_last_ndep, model_year_align_ndep, ndep_taxmod, " .
                  "ndep_varlist, nor stream_fldfilename_ndep" .
                  " can be set!");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_cnmresp {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  #
  # CN Maintence respiration for bgc=CN
  #
  if ( $nl_flags->{'bgc_mode'} ne "sp" ) {
    # When FUN is on get a default value
    if ( &value_is_true( $nl->get_value('use_fun') ) ) {
       add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults,
                   $nl, 'br_root', 'phys'=>$nl_flags->{'phys'},
                   'use_fun'=>$nl->get_value('use_fun'),
                   'use_cn'=>$nl_flags->{'use_cn'} );
    }
  } else {
    # If bgc is NOT CN/CNDV then make sure not set
    if ( defined($nl->get_value('br_root'))) {
      $log->fatal_error("br_root can NOT be set when bgc_mode==sp!");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_photosyns {
  # MUST be after use_hydrstress is set
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  #
  # Photosynthesis
  #
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults,
              $nl, 'rootstem_acc', 'phys'=>$nl_flags->{'phys'});
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults,
              $nl, 'light_inhibit', 'phys'=>$nl_flags->{'phys'});
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults,
              $nl, 'leafresp_method', 'phys'=>$nl_flags->{'phys'},
              'use_cn'=>$nl_flags->{'use_cn'});
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults,
              $nl, 'modifyphoto_and_lmr_forcrop', 'phys'=>$nl_flags->{'phys'} );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults,
              $nl, 'stomatalcond_method', 'phys'=>$nl_flags->{'phys'},
              'use_hydrstress'=>$nl_flags->{'use_hydrstress'} );
  # When CN on, must NOT be scaled by vcmax25top
  if ( &value_is_true( $nl_flags->{'use_cn'} ) ) {
     if ( $nl->get_value('leafresp_method') == 0 ) {
        $log->fatal_error("leafresp_method can NOT be set to scaled to vcmax (0) when CN is on!");
     }
     # And when CN off, must NOT be anything besides scaled by vxmac25top
  } else {
     if ( $nl->get_value('leafresp_method') != 0 ) {
        $log->fatal_error("leafresp_method can NOT be set to anything besides scaled to vcmax (0) when bgc_mode==sp!");
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_canopy {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;
  #
  # Canopy state
  #
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults,
              $nl, 'leaf_mr_vcm', 'phys'=>$nl_flags->{'phys'} )
}

#-------------------------------------------------------------------------------

sub setup_logic_popd_streams {
  # population density streams require CN/BGC
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( &value_is_true($nl_flags->{'cnfireson'}) ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'popdensmapalgo', 'hgrid'=>$nl_flags->{'res'},
                 'clm_accelerated_spinup'=>$nl_flags->{'clm_accelerated_spinup'}, 'cnfireson'=>$nl_flags->{'cnfireson'}  );
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_first_popdens', 'phys'=>$nl_flags->{'phys'},
                 'cnfireson'=>$nl_flags->{'cnfireson'}, 'sim_year'=>$nl_flags->{'sim_year'},
                 'sim_year_range'=>$nl_flags->{'sim_year_range'});
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_last_popdens', 'phys'=>$nl_flags->{'phys'},
                 'cnfireson'=>$nl_flags->{'cnfireson'}, 'sim_year'=>$nl_flags->{'sim_year'},
                 'sim_year_range'=>$nl_flags->{'sim_year_range'});
     # Set align year, if first and last years are different
     if ( $nl->get_value('stream_year_first_popdens') !=
          $nl->get_value('stream_year_last_popdens') ) {
        add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'model_year_align_popdens', 'sim_year'=>$nl_flags->{'sim_year'},
                    'sim_year_range'=>$nl_flags->{'sim_year_range'}, 'cnfireson'=>$nl_flags->{'cnfireson'});
     }
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_popdens', 'phys'=>$nl_flags->{'phys'},
                 'cnfireson'=>$nl_flags->{'cnfireson'}, 'hgrid'=>"0.5x0.5", 'ssp_rcp'=>$nl_flags->{'ssp_rcp'} );
     #
     # TODO (mvertens, 2021-06-22) the following is needed for MCT since a use case enforces this  - so for now stream_meshfile_popdens will be added to the mct
     # stream namelist but simply not used
    if ($opts->{'driver'} eq "nuopc" ) {
        add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_meshfile_popdens', 'hgrid'=>"0.5x0.5");
        my $inputdata_rootdir = $nl_flags->{'inputdata_rootdir'};
        my $default_value = $nl->get_value('stream_meshfile_popdens');
        my $none_filename = $inputdata_rootdir . '/none';
        my $none_filename = &quote_string($none_filename);
        if ($default_value eq $none_filename) {
            my $var = 'stream_meshfile_popdens';
            my $group = $definition->get_group_name($var);
            my $val = "none";
            $val = &quote_string( $val );
            $nl->set_variable_value($group, $var, $val);
        }
    } else {
        my $var = 'stream_meshfile_popdens';
        my $group = $definition->get_group_name($var);
        my $val = "none";
        $val = &quote_string( $val );
        $nl->set_variable_value($group, $var, $val);
    }
  } else {
     # If bgc is NOT CN/CNDV or fire_method==nofire then make sure none of the popdens settings are set
     if ( defined($nl->get_value('stream_year_first_popdens')) ||
          defined($nl->get_value('stream_year_last_popdens'))  ||
          defined($nl->get_value('model_year_align_popdens'))  ||
          defined($nl->get_value('popdens_tintalgo'        ))  ||
          defined($nl->get_value('stream_fldfilename_popdens'))   ) {
        $log->fatal_error("When bgc is SP (NOT CN or BGC) or fire_method==nofire none of: stream_year_first_popdens,\n" .
                          "stream_year_last_popdens, model_year_align_popdens, popdens_tintalgo nor\n" .
                          "stream_fldfilename_popdens can be set!");
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_urbantv_streams {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'urbantvmapalgo',
              'hgrid'=>$nl_flags->{'res'} );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_first_urbantv', 'phys'=>$nl_flags->{'phys'},
              'sim_year'=>$nl_flags->{'sim_year'},
              'sim_year_range'=>$nl_flags->{'sim_year_range'});
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_last_urbantv', 'phys'=>$nl_flags->{'phys'},
              'sim_year'=>$nl_flags->{'sim_year'},
              'sim_year_range'=>$nl_flags->{'sim_year_range'});
  # Set align year, if first and last years are different
  if ( $nl->get_value('stream_year_first_urbantv') !=
       $nl->get_value('stream_year_last_urbantv') ) {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
                 'model_year_align_urbantv', 'sim_year'=>$nl_flags->{'sim_year'},
                 'sim_year_range'=>$nl_flags->{'sim_year_range'});
  }
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_urbantv', 'phys'=>$nl_flags->{'phys'},
              'hgrid'=>"0.9x1.25" );
  if ($opts->{'driver'} eq "nuopc" ) {
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_meshfile_urbantv', 'phys'=>$nl_flags->{'phys'},
                  'hgrid'=>"0.9x1.25" );
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_lightning_streams {
  # lightning streams require CN/BGC
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

    if ( $nl_flags->{'light_res'} ne "none" ) {
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'lightngmapalgo',
                  'hgrid'=>$nl_flags->{'res'},
                  'clm_accelerated_spinup'=>$nl_flags->{'clm_accelerated_spinup'}  );
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_first_lightng',
                  'sim_year'=>$nl_flags->{'sim_year'},
                  'sim_year_range'=>$nl_flags->{'sim_year_range'});
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_last_lightng',
                  'sim_year'=>$nl_flags->{'sim_year'},
                  'sim_year_range'=>$nl_flags->{'sim_year_range'});
      # Set align year, if first and last years are different
      if ( $nl->get_value('stream_year_first_lightng') !=
           $nl->get_value('stream_year_last_lightng') ) {
        add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'model_year_align_lightng', 'sim_year'=>$nl_flags->{'sim_year'},
                    'sim_year_range'=>$nl_flags->{'sim_year_range'});
     }
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_lightng',
                 'hgrid'=>$nl_flags->{'light_res'} );
      if ($opts->{'driver'} eq "nuopc" ) {
          add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_meshfile_lightng',
                      'hgrid'=>$nl_flags->{'light_res'} );
      }
  } else {
     # If bgc is NOT CN/CNDV then make sure none of the Lightng settings are set
     if ( defined($nl->get_value('stream_year_first_lightng')) ||
          defined($nl->get_value('stream_year_last_lightng'))  ||
          defined($nl->get_value('model_year_align_lightng'))  ||
          defined($nl->get_value('lightng_tintalgo'        ))  ||
          defined($nl->get_value('stream_fldfilename_lightng'))   ) {
        $log->fatal_error("When bgc is SP (NOT CN or BGC or FATES) or fire is turned off none of: stream_year_first_lightng,\n" .
                          "stream_year_last_lightng, model_year_align_lightng, lightng_tintalgo nor\n" .
                          "stream_fldfilename_lightng can be set!");
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_dry_deposition {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ($opts->{'drydep'} ) {
    if ( &value_is_true( $nl_flags->{'use_fates'}) && not &value_is_true( $nl_flags->{'use_fates_sp'}) ) {
       $log->warning("DryDeposition can NOT be on when FATES is also on unless FATES-SP mode is on.\n" .
                     "   Use the '--no-drydep' option when '-bgc fates' is activated");
    }
    add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'drydep_list');
    add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'dep_data_file');
  } else {
    if ( defined($nl->get_value('drydep_list')) ) {
      $log->fatal_error("drydep_list defined, but drydep option NOT set");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_fire_emis {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ($opts->{'fire_emis'} ) {
     if ( &value_is_true( $nl_flags->{'use_fates'} ) ) {
       $log->warning("Fire emission can NOT be on when FATES is also on.\n" .
                   "  DON'T use the '-fire_emis' option when '-bgc fates' is activated");
    }
    add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fire_emis_factors_file');
    add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fire_emis_specifier');
  } else {
    if ( defined($nl->get_value('fire_emis_elevated'))     ||
         defined($nl->get_value('fire_emis_factors_file')) ||
         defined($nl->get_value('fire_emis_specifier')) ) {
      $log->fatal_error("fire_emission setting defined: fire_emis_elevated, fire_emis_factors_file, or fire_emis_specifier, but fire_emis option NOT set");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_dust_emis {
  # Logic to handle the dust emissions
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  # First get the dust emission method
  my $var = "dust_emis_method";
  add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var );

  my $dust_emis_method = remove_leading_and_trailing_quotes( $nl->get_value($var) );

  my @zender_files_in_lnd_opts = ( "stream_fldfilename_zendersoilerod", "stream_meshfile_zendersoilerod",
                                   "zendersoilerod_mapalgo" );
  if ( $dust_emis_method eq "Zender_2003" ) {
     # get the zender_soil_erod_source
     add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
                 "zender_soil_erod_source", 'dust_emis_method'=>$dust_emis_method );

     my $zender_source = remove_leading_and_trailing_quotes( $nl->get_value('zender_soil_erod_source') );
     if ( $zender_source eq "lnd" ) {
        foreach my $option ( @zender_files_in_lnd_opts ) {
           add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $option,
                       'dust_emis_method'=>$dust_emis_method, 'zender_soil_erod_source'=>$zender_source,
                       'hgrid'=>$nl_flags->{'res'}, 'lnd_tuning_mod'=>$nl_flags->{'lnd_tuning_mode'} );
        }
     } else {
        foreach my $option ( @zender_files_in_lnd_opts ) {
           if ( defined($nl->get_value($option)) ) {
             $log->fatal_error("zender_soil_erod_source is NOT lnd, but the file option $option is being set" .
                               " and should NOT be unless you want it handled here in the LAND model, " .
                               "otherwise the equivalent option is set in CAM" );
           }
        }
     }
  } else {
     # Verify that NONE of the Zender options are being set if Zender is NOT being used
     push @zender_files_in_lnd_opts, "zender_soil_erod_source";
     foreach my $option ( @zender_files_in_lnd_opts ) {
        if ( defined($nl->get_value($option)) ) {
          $log->fatal_error("dust_emis_method is NOT set to Zender_2003, but one of it's options " .
                            "$option is being set, need to change one or the other" );
        }
     }
     if ( $dust_emis_method eq "Leung_2023" ) {
        $log->warning("dust_emis_method is Leung_2023 and that option has NOT been brought into CTSM yet");
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_megan {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  my $var   = "megan";

  if ( $opts->{$var} eq "default" ) {
     add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'megan',
                 'clm_accelerated_spinup'=>$nl_flags->{'clm_accelerated_spinup'},
                 'configuration'=>$nl_flags->{'configuration'} );
    $nl_flags->{$var} = $nl->get_value($var);
  } else {
    $nl_flags->{$var} = $opts->{$var};
  }

  if ($nl_flags->{'megan'} ) {
    if ( &value_is_true( $nl_flags->{'use_fates'} ) ) {
       $log->warning("MEGAN can NOT be on when FATES is also on.\n" .
                   "   Use the '-no-megan' option when '-bgc fates' is activated");
    }
    add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'megan_specifier');
    check_megan_spec( $opts, $nl, $definition );
    add_default($opts,  $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'megan_factors_file');
  } else {
    if ( defined($nl->get_value('megan_specifier')) ||
         defined($nl->get_value('megan_factors_file')) ) {
      $log->fatal_error("megan_specifier or megan_factors_file defined, but megan option NOT set");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_soilm_streams {
  my ($opts, $nl_flags, $definition, $defaults, $nl, $physv) = @_;

      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_soil_moisture_streams');
      if ( &value_is_true( $nl->get_value('use_soil_moisture_streams') ) ) {
         add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'soilm_tintalgo',
                     'hgrid'=>$nl_flags->{'res'} );
         add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'soilm_offset',
                     'hgrid'=>$nl_flags->{'res'} );
         add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_first_soilm', 'phys'=>$nl_flags->{'phys'},
                     'sim_year'=>$nl_flags->{'sim_year'},
                     'sim_year_range'=>$nl_flags->{'sim_year_range'});
         add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_last_soilm', 'phys'=>$nl_flags->{'phys'},
                     'sim_year'=>$nl_flags->{'sim_year'},
                     'sim_year_range'=>$nl_flags->{'sim_year_range'});
         # Set align year, if first and last years are different
         if ( $nl->get_value('stream_year_first_soilm') !=
              $nl->get_value('stream_year_last_soilm') ) {
              add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
                          'model_year_align_soilm', 'sim_year'=>$nl_flags->{'sim_year'},
                          'sim_year_range'=>$nl_flags->{'sim_year_range'});
         }
         add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_soilm', 'phys'=>$nl_flags->{'phys'},
                     'hgrid'=>$nl_flags->{'res'} );
         if ( ($opts->{'use_case'} =~ /_transient$/) &&
              (remove_leading_and_trailing_quotes($nl->get_value("soilm_tintalgo")) eq "linear") ) {
             $log->warning("For a transient case, soil moisture streams, should NOT use soilm_tintalgo='linear'" .
                           " since vegetated areas could go from missing to not missing or vice versa" );
         }
      } else {
         if ( defined($nl->get_value('stream_year_first_soilm')) ||
              defined($nl->get_value('model_year_align_soilm')) ||
              defined($nl->get_value('stream_fldfilename_soilm')) ||
              defined($nl->get_value('soilm_tintalgo')) ||
              defined($nl->get_value('soilm_offset')) ||
              defined($nl->get_value('stream_year_last_soilm')) ) {
             $log->fatal_error("One of the soilm streams namelist items (stream_year_first_soilm, " .
                                " model_year_align_soilm, stream_fldfilename_soilm, stream_fldfilename_soilm)" .
                                " soilm_tintalgo soilm_offset" .
                                " is defined, but use_soil_moisture_streams option NOT set to true");
         }
      }
}

#-------------------------------------------------------------------------------

sub setup_logic_lai_streams {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_lai_streams');

  if ( &value_is_true($nl_flags->{'use_crop'}) && &value_is_true($nl->get_value('use_lai_streams'))  ) {
    $log->fatal_error("turning use_lai_streams on is incompatable with use_crop set to true.");
  }
  if ( $nl_flags->{'bgc_mode'} eq "sp" || ($nl_flags->{'bgc_mode'} eq "fates" && &value_is_true($nl->get_value('use_fates_sp')) )) {
     if ( &value_is_true($nl->get_value('use_lai_streams')) ) {
       add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_lai_streams');
       add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'lai_mapalgo',
                   'hgrid'=>$nl_flags->{'res'} );
       add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_first_lai',
                   'sim_year'=>$nl_flags->{'sim_year'},
                   'sim_year_range'=>$nl_flags->{'sim_year_range'});
       add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_last_lai',
                   'sim_year'=>$nl_flags->{'sim_year'},
                   'sim_year_range'=>$nl_flags->{'sim_year_range'});
       # Set align year, if first and last years are different
       if ( $nl->get_value('stream_year_first_lai') !=
            $nl->get_value('stream_year_last_lai') ) {
          add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
                      'model_year_align_lai', 'sim_year'=>$nl_flags->{'sim_year'},
                      'sim_year_range'=>$nl_flags->{'sim_year_range'});
       }
       add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_lai',
                   'hgrid'=>"360x720cru" );
       if ($opts->{'driver'} eq "nuopc" ) {
           add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_meshfile_lai',
                       'hgrid'=>"360x720cru" );
       }
     }
  } else {
     # If bgc is BGC/BGCDV then make sure none of the LAI settings are set
     if ( &value_is_true($nl->get_value('use_lai_streams'))) {
        $log->fatal_error("When not in SP mode use_lai_streams cannot be .true.\n" .
                          "(eg. don't use this option with BGC or non-SP FATES), \n" .
                          "Update compset to use SP)");
     }
     if ( defined($nl->get_value('stream_year_first_lai'))  ||
          defined($nl->get_value('stream_year_last_lai'))   ||
          defined($nl->get_value('model_year_align_lai'))   ||
          defined($nl->get_value('lai_tintalgo'        ))   ||
          defined($nl->get_value('stream_fldfilename_lai'))   ) {
        $log->fatal_error("When not in SP mode none of the following can be set: stream_year_first_lai,\n" .
                          "stream_year_last_lai, model_year_align_lai, lai_tintalgo nor\n" .
                          "stream_fldfilename_lai (eg. don't use this option with BGC or FATES-SP).");
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_cropcal_streams {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  # Set first and last stream years
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_first_cropcal',
              'sim_year'=>$nl_flags->{'sim_year'},
              'sim_year_range'=>$nl_flags->{'sim_year_range'});
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_year_last_cropcal',
              'sim_year'=>$nl_flags->{'sim_year'},
              'sim_year_range'=>$nl_flags->{'sim_year_range'});

  # Set align year, if first and last years are different
  if ( $nl->get_value('stream_year_first_cropcal') !=
      $nl->get_value('stream_year_last_cropcal') ) {
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
                'model_year_align_cropcal', 'sim_year'=>$nl_flags->{'sim_year'},
                'sim_year_range'=>$nl_flags->{'sim_year_range'});
  }

  # Set up other crop calendar parameters
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'generate_crop_gdds');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_mxmat');

  # Option checks
  my $generate_crop_gdds = $nl->get_value('generate_crop_gdds') ;
  my $use_mxmat = $nl->get_value('use_mxmat') ;
  my $swindow_start_file = $nl->get_value('stream_fldFileName_swindow_start') ;
  my $swindow_end_file = $nl->get_value('stream_fldFileName_swindow_end') ;
  my $gdd_file = $nl->get_value('stream_fldFileName_cultivar_gdds') ;
  my $mesh_file = $nl->get_value('stream_meshfile_cropcal') ;
  if ( ($swindow_start_file eq '' and $swindow_start_file ne '') or ($swindow_start_file ne '' and $swindow_start_file eq '') ) {
    $log->fatal_error("When specifying sowing window dates, you must provide both swindow_start_file and swindow_end_file. To specify exact sowing dates, use the same file." );
  }
  if ( $generate_crop_gdds eq '.true.' ) {
      if ( $use_mxmat eq '.true.' ) {
          $log->fatal_error("If generate_crop_gdds is true, you must also set use_mxmat to false" );
      }
      if ( $swindow_start_file eq '' or $swindow_end_file eq '' ) {
          $log->fatal_error("If generate_crop_gdds is true, you must specify stream_fldFileName_swindow_start and stream_fldFileName_swindow_end")
      }
      if ( $swindow_start_file ne $swindow_end_file ) {
          $log->fatal_error("If generate_crop_gdds is true, you must specify exact sowing dates by setting stream_fldFileName_swindow_start and stream_fldFileName_swindow_end to the same file")
      }
      if ( $gdd_file ne '' ) {
          $log->fatal_error("If generate_crop_gdds is true, do not specify stream_fldFileName_cultivar_gdds")
      }
  }
  if ( $mesh_file eq '' and ( $swindow_start_file ne '' or $gdd_file ne '' ) ) {
      $log->fatal_error("If prescribing crop sowing dates and/or maturity requirements, you must specify stream_meshfile_cropcal")
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_soilwater_movement {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'soilwater_movement_method' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'upper_boundary_condition' );

  my $soilmtd = $nl->get_value("soilwater_movement_method");
  my $use_bed = $nl->get_value('use_bedrock'              );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl,
              'lower_boundary_condition', 'vichydro'=>$nl_flags->{'vichydro'},
              'soilwater_movement_method'=>$soilmtd, 'use_bedrock'=>$use_bed
             );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'dtmin' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'verySmall' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'xTolerUpper' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'xTolerLower' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'expensive' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'inexpensive' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'flux_calculation' );
}
#-------------------------------------------------------------------------------

sub setup_logic_cnvegcarbonstate {
  #  MUST be AFTER: setup_logic_dynamic_plant_nitrogen_alloc as depends on mm_nuptake_opt which is set there
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( &value_is_true($nl->get_value('use_cn')) ) {
    my $mmnuptake = $nl->get_value('mm_nuptake_opt');
    if ( ! defined($mmnuptake) ) { $mmnuptake = ".false."; }
    add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'initial_vegC',
                'use_cn' => $nl->get_value('use_cn'), 'mm_nuptake_opt' => $mmnuptake );
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_cngeneral {
  # Must be set after setup_logic_co2_type
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if ( &value_is_true($nl->get_value('use_cn')) ) {
    if ( &value_is_true($nl->get_value('use_crop')) ) {
       add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'dribble_crophrv_xsmrpool_2atm',
                   'co2_type' => remove_leading_and_trailing_quotes($nl->get_value('co2_type')),
                   'use_crop' => $nl->get_value('use_crop')  );
    } else {
      if ( defined($nl->get_value('dribble_crophrv_xsmrpool_2atm')) ) {
        $log->fatal_error("When CROP is NOT on dribble_crophrv_xsmrpool_2atm can NOT be set\n" );
      }
    }
  } else {
    if ( defined($nl->get_value('reseed_dead_plants')) ||
         defined($nl->get_value('dribble_crophrv_xsmrpool_2atm'))   ) {
             $log->fatal_error("When CN is not on none of the following can be set: ,\n" .
                  "dribble_crophrv_xsmrpool_2atm nor reseed_dead_plantsr\n" .
                  "(eg. don't use these options with SP mode).");
    }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_rooting_profile {
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'rooting_profile_method_water' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'rooting_profile_method_carbon' );
}

#-------------------------------------------------------------------------------

sub setup_logic_friction_vel {
  # Must be after canopyfluxes so that use_biomass_heat_storage will be set
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'zetamaxstable',
     'use_biomass_heat_storage'=>$nl_flags->{'use_biomass_heat_storage'}, 'phys'=>$nl_flags->{'phys'} );
}

#-------------------------------------------------------------------------------

sub setup_logic_soil_resis {
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'soil_resis_method' );
}

sub setup_logic_canopyfluxes {
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_undercanopy_stability' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'itmax_canopy_fluxes',
              'structure'=>$nl_flags->{'structure'});
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_biomass_heat_storage',
              'use_fates'=>$nl_flags->{'use_fates'}, 'phys'=>$nl_flags->{'phys'} );
  if ( &value_is_true($nl->get_value('use_biomass_heat_storage') ) && &value_is_true( $nl_flags->{'use_fates'}) ) {
     $log->fatal_error('use_biomass_heat_storage can NOT be set to true when fates is on');
  }
  if ( &value_is_true($nl->get_value('use_biomass_heat_storage')) ) {
     $nl_flags->{'use_biomass_heat_storage'} = ".true.";
  } else {
     $nl_flags->{'use_biomass_heat_storage'} = ".false.";
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_canopyhydrology {
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'interception_fraction' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'maximum_leaf_wetted_fraction' );
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_clm5_fpi' );
}

#-------------------------------------------------------------------------------

sub setup_logic_snowpack {
  #
  # Snowpack related options
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'nlevsno',
              'structure'=>$nl_flags->{'structure'});
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'h2osno_max',
              'structure'=>$nl_flags->{'structure'});
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'wind_dependent_snow_density');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snow_overburden_compaction_method');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'lotmp_snowdensity_method');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'upplim_destruct_metamorph');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'reset_snow');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'reset_snow_glc');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'reset_snow_glc_ela');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snow_dzmin_1');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snow_dzmin_2');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snow_dzmax_l_1');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snow_dzmax_l_2');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snow_dzmax_u_1');
  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'snow_dzmax_u_2');

  my $dzmin1 = $nl->get_value('snow_dzmin_1');
  my $dzmin2 = $nl->get_value('snow_dzmin_2');
  my $dzmax_l1 = $nl->get_value('snow_dzmax_l_1');
  my $dzmax_l2 = $nl->get_value('snow_dzmax_l_2');
  my $dzmax_u1 = $nl->get_value('snow_dzmax_u_1');
  my $dzmax_u2 = $nl->get_value('snow_dzmax_u_2');

  if ($dzmin1 != 0.01 || $dzmin2 != 0.015 || $dzmax_u1 != 0.02 || $dzmax_u2 != 0.05 || $dzmax_l1 != 0.03 || $dzmax_l2 != 0.07) {
     $log->warning("Setting any of the following namelist variables to NON DEFAULT values remains untested as of Sep 6, 2019: snow_dzmin_1 & 2, snow_dzmax_u_1 & 2, snow_dzmax_l_1 & 2." );
     $log->warning("Leave these variables unspecified in user_nl_clm in order to use the default values." );
  }
  if ($dzmin1 <= 0.0 || $dzmin2 <= 0.0 || $dzmax_u1 <= 0.0 || $dzmax_u2 <= 0.0 || $dzmax_l1 <= 0.0 || $dzmax_l2 <= 0.0) {
     $log->fatal_error('One or more of the snow_dzmin_* and/or snow_dzmax_* were set incorrectly to be <= 0');
  }
  if ($dzmin2 <= $dzmin1) {
     $log->fatal_error('snow_dzmin_2 was set incorrectly to be <= snow_dzmin_1');
  }
  if ($dzmax_l2 <= $dzmax_l1) {
     $log->fatal_error('snow_dzmax_l_2 was set incorrectly to be <= snow_dzmax_l_1');
  }
  if ($dzmax_u2 <= $dzmax_u1) {
     $log->fatal_error('snow_dzmax_u_2 was set incorrectly to be <= snow_dzmax_u_1');
  }
  if ($dzmin1 >= $dzmax_u1) {
     $log->fatal_error('snow_dzmin_1 was set incorrectly to be >= snow_dzmax_u_1');
  }
  if ($dzmin2 >= $dzmax_u2) {
     $log->fatal_error('snow_dzmin_2 was set incorrectly to be >= snow_dzmax_u_2');
  }
  if ($dzmax_u1 >= $dzmax_l1) {
     $log->fatal_error('snow_dzmax_u_1 was set incorrectly to be >= snow_dzmax_l_1');
  }
  if ($dzmax_u2 >= $dzmax_l2) {
     $log->fatal_error('snow_dzmax_u_2 was set incorrectly to be >= snow_dzmax_l_2');
  }

  if (remove_leading_and_trailing_quotes($nl->get_value('snow_overburden_compaction_method')) eq 'Vionnet2012') {
     # overburden_compress_tfactor isn't used if we're using the Vionnet2012
     # snow overburden compaction method, so make sure the user hasn't tried
     # to set it
     if (defined($nl->get_value('overburden_compress_tfactor'))) {
        $log->fatal_error('overburden_compress_tfactor is set, but does not apply when using snow_overburden_compaction_method=Vionnet2012');
     }
  } else {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'overburden_compress_tfactor');
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_scf_SwensonLawrence2012 {
   # Options related to the SwensonLawrence2012 snow cover fraction method
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  if (remove_leading_and_trailing_quotes($nl->get_value('snow_cover_fraction_method')) eq 'SwensonLawrence2012') {
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'int_snow_max');
     add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'n_melt_glcmec');
  }
  else {
     if (defined($nl->get_value('int_snow_max'))) {
        $log->fatal_error('int_snow_max is set, but only applies for snow_cover_fraction_method=SwensonLawrence2012');
     }
     if (defined($nl->get_value('n_melt_glcmec'))) {
        $log->fatal_error('n_melt_glcmec is set, but only applies for snow_cover_fraction_method=SwensonLawrence2012');
     }
  }
}

#-------------------------------------------------------------------------------

sub setup_logic_atm_forcing {
   #
   # Options related to atmospheric forcings
   #
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'glcmec_downscale_longwave');
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'repartition_rain_snow');
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'lapse_rate');

   my $var;

   foreach $var ("lapse_rate_longwave",
                 "longwave_downscaling_limit") {
      if ( &value_is_true($nl->get_value("glcmec_downscale_longwave")) ) {
         add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var);
      } else {
         if (defined($nl->get_value($var))) {
            $log->fatal_error("$var can only be set if glcmec_downscale_longwave is true");
         }
      }
   }

   foreach $var ("precip_repartition_glc_all_snow_t",
                 "precip_repartition_glc_all_rain_t",
                 "precip_repartition_nonglc_all_snow_t",
                 "precip_repartition_nonglc_all_rain_t") {
      if ( &value_is_true($nl->get_value("repartition_rain_snow")) ) {
         add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var);
      } else {
         if (defined($nl->get_value($var))) {
            $log->fatal_error("$var can only be set if repartition_rain_snow is true");
         }
      }
   }
}

#-------------------------------------------------------------------------------

sub setup_logic_lnd2atm {
   #
   # Options related to fields sent to atmosphere
   #
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'melt_non_icesheet_ice_runoff');
}

#-------------------------------------------------------------------------------

sub setup_logic_initinterp {
   #
   # Options related to init_interp
   #
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   my $var = 'init_interp_method';
   if ( &value_is_true($nl->get_value("use_init_interp"))) {
      add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var);
   } else {
      if (defined($nl->get_value($var))) {
         $log->fatal_error("$var can only be set if use_init_interp is true");
      }
   }
}

#-------------------------------------------------------------------------------

sub setup_logic_fates {
    #
    # Set some default options related to Ecosystem Demography
    #
    my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

    if (&value_is_true( $nl_flags->{'use_fates'})  ) {
        add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'fates_paramfile', 'phys'=>$nl_flags->{'phys'});
        my @list  = (  "fates_spitfire_mode", "use_fates_planthydro", "use_fates_ed_st3", "use_fates_ed_prescribed_phys",
                       "use_fates_inventory_init","use_fates_fixed_biogeog","use_fates_nocomp","fates_seeddisp_cadence",
                       "use_fates_logging","fates_parteh_mode", "use_fates_cohort_age_tracking","use_fates_tree_damage",
		       "use_fates_luh","fates_history_dimlevel" );
        foreach my $var ( @list ) {
 	  add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'use_fates'=>$nl_flags->{'use_fates'},
                      'use_fates_sp'=>$nl_flags->{'use_fates_sp'} );
        }
        my $suplnitro = $nl->get_value('suplnitro');
        my $parteh_mode = $nl->get_value('fates_parteh_mode');
        if ( ($parteh_mode == 1) &&  ($suplnitro !~ /ALL/) && not &value_is_true( $nl_flags->{'use_fates_sp'}) ) {
          $log->fatal_error("supplemental Nitrogen (suplnitro) is NOT set to ALL, FATES is on, " . 
                            "but and FATES-SP is not active, but fates_parteh_mode is 1, so Nitrogen is not active" . 
                            "Change suplnitro back to ALL");
        }
        #
        # For FATES SP mode make sure no-competetiion, and fixed-biogeography are also set
        # And also check for other settings that can't be trigged on as well
        #
        my $var = "use_fates_sp";
        if ( defined($nl->get_value($var))  ) {
           if ( &value_is_true($nl->get_value($var)) ) {
              my @list = ( "use_fates_nocomp", "use_fates_fixed_biogeog" );
              foreach my $var ( @list ) {
                 if ( ! &value_is_true($nl->get_value($var)) ) {
                    $log->fatal_error("$var is required when FATES SP is on (use_fates_sp)" );
                 }
              }
              # spit-fire can't be on with FATES SP mode is active
              if ( $nl->get_value('fates_spitfire_mode') > 0 ) {
                    $log->fatal_error('fates_spitfire_mode can NOT be set to greater than 0 when use_fates_sp is true');
              }
              # hydro isn't currently supported to work when FATES SP mode is active
              if (&value_is_true( $nl->get_value('use_fates_planthydro') )) {
                    $log->fatal_error('fates sp mode is currently not supported to work with fates hydro');
              }
           }
        }
        my $var = "use_fates_inventory_init";
        if ( defined($nl->get_value($var))  ) {
           if ( &value_is_true($nl->get_value($var)) ) {
              $var = "fates_inventory_ctrl_filename";
              my $fname = remove_leading_and_trailing_quotes( $nl->get_value($var) );
              if ( ! defined($nl->get_value($var))  ) {
                 $log->fatal_error("$var is required when use_fates_inventory_init is set" );
              } elsif ( ! -f "$fname" ) {
                 $log->fatal_error("$fname does NOT point to a valid filename" );
              }
           }
        }
        my $var = "use_fates_luh";
        if ( defined($nl->get_value($var))  ) {
           if ( &value_is_true($nl->get_value($var)) ) {
              $var = "fluh_timeseries";
              add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, $var, 'phys'=>$nl_flags->{'phys'}, 'hgrid'=>$nl_flags->{'res'}, 'sim_year_range'=>$nl_flags->{'sim_year_range'}, nofail=>1 );
              my $fname = remove_leading_and_trailing_quotes( $nl->get_value($var) );
              if ( ! defined($nl->get_value($var))  ) {
                 $log->fatal_error("$var is required when use_fates_luh is set" );
              } elsif ( ! -f "$fname" ) {
                 $log->fatal_error("$fname does NOT point to a valid filename" );
              }
           }
        }
    }
}

#-------------------------------------------------------------------------------
sub setup_logic_exice {
  #
  # excess ice streams
  #
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;
  my $use_exice = $nl->get_value( 'use_excess_ice' );
  my $use_exice_streams = $nl->get_value( 'use_excess_ice_streams' );
  # IF excess ice streams is on
  if (defined($use_exice_streams) && value_is_true($use_exice_streams)) {
     # Can only be true if excess ice is also on, otherwise fail
     if (defined($use_exice) && not value_is_true($use_exice)) {
        $log->fatal_error("use_excess_ice_streams can NOT be TRUE when use_excess_ice is FALSE" );
     }
  # Otherwise if ice streams are off
  } else {
     my @list = ( "stream_meshfile_exice", "stream_fldfilename_exice" );
     # fail is excess ice streams files are set
     foreach my $var ( @list ) {
        if ( defined($nl->get_value($var)) ) {
           $log->fatal_error("$var should NOT be set when use_excess_ice_streams=FALSE" );
        }
     }
     # mapalgo can only be none, if excess ice streams are off
     my $map_algo = $nl->get_value("stream_mapalgo_exice");
     if ( defined($map_algo) && ($map_algo ne "none") ) {
        $log->fatal_error("stream_mapalgo_exice can ONLY be none when use_excess_ice_streams=FALSE" );
     }
  }
  # If excess ice is on
  if (defined($use_exice) && value_is_true($use_exice)) {
     # IF nuopc driver and excess ice streams are on get the stream defaults
     if (defined($use_exice_streams) && value_is_true($use_exice_streams)) {
       add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_fldfilename_exice');
       add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_mapalgo_exice');
       # If excess ice streams on, but NOT the NUOPC driver fail
       if ( not $opts->{'driver'} eq "nuopc" ) {
          $log->fatal_error("nuopc driver is required when use_excess_ice_streams is set to true" );
       # NUOPC driver needs a mesh file
       } else {
          add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'stream_meshfile_exice');
       }
     }
  }


} # end exice streams

#-------------------------------------------------------------------------------

sub setup_logic_z0param {
   #
   # Set default z0 paramterization
   #
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'z0param_method');

   my $z0param_method = remove_leading_and_trailing_quotes($nl->get_value('z0param_method' ));
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'use_z0m_snowmelt',
           'z0param_method'=>$z0param_method );

   my $use_z0m_snowmelt = $nl->get_value( 'use_z0m_snowmelt' );

   if ( $z0param_method eq "ZengWang2007" && defined($use_z0m_snowmelt) && value_is_true($use_z0m_snowmelt)) {
      $log->fatal_error("use_z0m_snowmelt must be .false. when z0param_method = $z0param_method.\n $@");
   }

}

#-------------------------------------------------------------------------------

sub setup_logic_misc {
   #
   # Set some misc options
   #
   my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

   if ( $opts->{'driver'} ne "nuopc" ) {
      my $var = "force_send_to_atm";
      my $val = $nl->get_value($var);
      if ( defined($val) ) {
         $log->fatal_error( "$var can only be set for the nuopc driver" );
      }
   }
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'for_testing_run_ncdiopio_tests');
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'for_testing_use_second_grain_pool');
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'for_testing_use_repr_structure_pool');
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'for_testing_no_crop_seed_replenishment');
   add_default($opts, $nl_flags->{'inputdata_rootdir'}, $definition, $defaults, $nl, 'hist_fields_list_file');
}

#-------------------------------------------------------------------------------

sub write_output_files {
  my ($opts, $nl_flags, $defaults, $nl) = @_;

  my $note = "";
  my $var = "note";
  if ( ! defined($opts->{$var}) ) {
    $opts->{$var} = $defaults->get_value($var);
  }
  if ( $opts->{$var} ) {
    $note = "Comment:\n" .
      "This namelist was created using the following command-line:\n" .
        "    $nl_flags->{'cfgdir'}/$ProgName $nl_flags->{'cmdline'}\n" .
          "For help on options use: $nl_flags->{'cfgdir'}/$ProgName -help";
  }

  # CLM component
  my @groups;

  @groups = qw(clm_inparm ndepdyn_nml popd_streams urbantv_streams light_streams
               soil_moisture_streams lai_streams atm2lnd_inparm lnd2atm_inparm clm_canopyhydrology_inparm cnphenology
               cropcal_streams
               clm_soilhydrology_inparm dynamic_subgrid cnvegcarbonstate
               finidat_consistency_checks dynpft_consistency_checks
               clm_initinterp_inparm century_soilbgcdecompcascade
               soilhydrology_inparm luna friction_velocity mineral_nitrogen_dynamics
               soilwater_movement_inparm rooting_profile_inparm
               soil_resis_inparm  bgc_shared canopyfluxes_inparm aerosol
               clmu_inparm clm_soilstate_inparm clm_nitrogen clm_snowhydrology_inparm hillslope_hydrology_inparm hillslope_properties_inparm
               cnprecision_inparm clm_glacier_behavior crop_inparm irrigation_inparm
               surfacealbedo_inparm water_tracers_inparm tillage_inparm);

  #@groups = qw(clm_inparm clm_canopyhydrology_inparm clm_soilhydrology_inparm
  #             finidat_consistency_checks dynpft_consistency_checks);
  # Eventually only list namelists that are actually used when CN on
  if ( &value_is_true($nl_flags->{'use_lch4'}) ) {
     push @groups, "ch4par_in";
  }
  if ( $opts->{'driver'} eq "nuopc" ) {
     push @groups, "ctsm_nuopc_cap";
  }
  push @groups, "clm_humanindex_inparm";
  push @groups, "cnmresp_inparm";
  push @groups, "photosyns_inparm";
  push @groups, "cnfire_inparm";
  push @groups, "cn_general";
  push @groups, "nitrif_inparm";
  push @groups, "lifire_inparm";
  push @groups, "ch4finundated";
  push @groups, "exice_streams";
  push @groups, "soilbgc_decomp";
  push @groups, "clm_canopy_inparm";
  push @groups, "zendersoilerod";
  if (remove_leading_and_trailing_quotes($nl->get_value('snow_cover_fraction_method')) eq 'SwensonLawrence2012') {
     push @groups, "scf_swenson_lawrence_2012_inparm";
  }

  my $outfile;
  $outfile = "$opts->{'dir'}/lnd_in";
  $nl->write($outfile, 'groups'=>\@groups, 'note'=>"$note" );
  $log->verbose_message("Writing clm namelist to $outfile");

  # Drydep, fire-emission or MEGAN namelist for driver
  @groups = qw(drydep_inparm megan_emis_nl fire_emis_nl carma_inparm);
  $outfile = "$opts->{'dir'}/drv_flds_in";
  $nl->write($outfile, 'groups'=>\@groups, 'note'=>"$note" );
  $log->verbose_message("Writing @groups namelists to $outfile");
}

sub write_output_real_parameter_file {
  my ($opts, $nl_flags, $definition, $defaults, $nl) = @_;

  # Output real parameters
  if ( defined($opts->{'output_reals_filename'}) ) {
     my $file = $opts->{'output_reals_filename'};
     my $fh = IO::File->new($file, '>>') or $log->fatal_error("can't create real parameter filename: $file");
     foreach my $var ( $definition->get_var_names() ) {
        my $type = $definition->get_var_type($var);
        my $doc  = $definition->get_var_doc($var);
        $doc =~ s/\n//g;
        if ( $type =~ /real/ ) {
           my $val = $nl->get_value($var);
           if ( ! defined($val) ) { $val = "?.??"; }
           print $fh "\! $doc\n$var = $val\n";
        }
     }
     $fh->close();
  }
}

#-------------------------------------------------------------------------------

sub add_default {

# Add a value for the specified variable to the specified namelist object.  The variables
# already in the object have the higher precedence, so if the specified variable is already
# defined in the object then don't overwrite it, just return.
#
# This method checks the definition file and adds the variable to the correct
# namelist group.
#
# The value can be provided by using the optional argument key 'val' in the
# calling list.  Otherwise a default value is obtained from the namelist
# defaults object.  If no default value is found this method throws an exception
# unless the 'nofail' option is set true.
#
# Example 1: Specify the default value $val for the namelist variable $var in namelist
#            object $nl:
#
#  add_default($inputdata_rootdir, $definition, $defaults, $nl, $var, 'val'=>$val)
#
# Example 2: Add a default for variable $var if an appropriate value is found.  Otherwise
#            don't add the variable
#
#  add_default($inputdata_rootdir, $definition, $defaults, $nl, $var, 'nofail'=>1)
#
#
# ***** N.B. ***** This routine assumes the following variables are in package main::
#  $definition        -- the namelist definition object
#  $defaults          -- the namelist defaults object
#  $inputdata_rootdir -- CESM inputdata root directory

  my $opts = shift;
  my $inputdata_rootdir = shift;
  my $definition = shift;
  my $defaults = shift;
  my $nl = shift;
  my $var = shift;
  my %settings = @_;

  my $test_files = $opts->{'test'};
  #my $nl = shift;     # namelist object
  #my $var = shift;    # name of namelist variable
  #my %settings = @_;      # options

  # If variable has quotes around it
  if ( $var =~ /'(.+)'/ ) {
    $var = $1;
  }
  # Query the definition to find which group the variable belongs to.  Exit if not found.

  my $group = $definition->get_group_name($var);
  unless ($group) {
    my $fname = $definition->get_file_name();
    $log->fatal_error("variable \"$var\" not found in namelist definition file $fname.");
  }

  # check whether the variable has a value in the namelist object -- if so then skip to end
  my $val = $nl->get_variable_value($group, $var);
  if (! defined $val) {

    # Look for a specified value in the options hash

    if (defined $settings{'val'}) {
      $val = $settings{'val'};
    }
    # or else get a value from namelist defaults object.
    # Note that if the 'val' key isn't in the hash, then just pass anything else
    # in %settings to the get_value method to be used as attributes that are matched
    # when looking for default values.
    else {

      $val = $defaults->get_value($var, \%settings);

      # Truncate model_version appropriately

      if ( $var eq "model_version" ) {
        $val =~ /(URL: https:\/\/[a-zA-Z0-9._-]+\/)([a-zA-Z0-9\/._-]+)(\/bld\/.+)/;
        $val = $2;
      }
    }

    # if no value is found then exit w/ error (unless 'nofail' option set)
    unless ( defined($val) ) {
      unless ($settings{'nofail'}) {
        if ($var eq 'finidat') {
          $log->message("No default value found for $var.\n" .
                  "            Are defaults provided for this resolution and land mask?" );
        } else {
          $log->fatal_error("No default value found for $var.\n" .
                      "            Are defaults provided for this resolution and land mask?");
        }
      }
      else {
        return( 1 );
      }
    }

    # query the definition to find out if the variable is an input pathname
    my $is_input_pathname = $definition->is_input_pathname($var);

    # The default values for input pathnames are relative.  If the namelist
    # variable is defined to be an absolute pathname, then prepend
    # the CESM inputdata root directory.
    if (not defined $settings{'no_abspath'}) {
      if (defined $settings{'set_abspath'}) {
        $val = set_abs_filepath($val, $settings{'set_abspath'});
      } else {
        if ($is_input_pathname eq 'abs') {
          $val = set_abs_filepath($val, $inputdata_rootdir);
          if ( $test_files and ($val !~ /null|none/) and (! -f "$val") ) {
            $log->fatal_error("file not found: $var = $val");
          }
        }
      }
    }

    # query the definition to find out if the variable takes a string value.
    # The returned string length will be >0 if $var is a string, and 0 if not.
    my $str_len = $definition->get_str_len($var);

    # If the variable is a string, then add quotes if they're missing
    if ($str_len > 0) {
      $val = quote_string($val);
    }

    # set the value in the namelist
    $nl->set_variable_value($group, $var, $val);
  }
  return( 0 );

}

#-------------------------------------------------------------------------------

sub expand_xml_variables_in_namelist {
   # Go through all variables in the namelist and expand any XML env settings in them
   my ($nl, $xmlvar_ref) = @_;

   foreach my $group ( $nl->get_group_names() ) {
       foreach my $var ( $nl->get_variable_names($group) ) {
          my $val    = $nl->get_variable_value($group, $var);
          my $newval = SetupTools::expand_xml_var( $val, $xmlvar_ref );
          if ( $newval ne $val ) {
             $nl->set_variable_value($group, $var, $newval);
          }
       }
   }
}

#-------------------------------------------------------------------------------

sub check_input_files {

# For each variable in the namelist which is an input dataset, check to see if it
# exists locally.
#
# ***** N.B. ***** This routine assumes the following variables are in package main::
#  $definition        -- the namelist definition object
#  $nl                -- namelist object
#  $inputdata_rootdir -- if false prints test, else creates inputdata file

    my ($nl, $inputdata_rootdir, $outfile, $definition) = @_;

    open(OUTFILE, ">>$outfile") if defined $inputdata_rootdir;

    # Look through all namelist groups
    my @groups = $nl->get_group_names();
    foreach my $group (@groups) {

        # Look through all variables in each group
        my @vars = $nl->get_variable_names($group);
        foreach my $var (@vars) {

            # Is the variable an input dataset?
            my $input_pathname_type = $definition->is_input_pathname($var);

            # If it is, check whether it exists locally and print status
            if ($input_pathname_type) {

                # Get pathname of input dataset
                my $pathname = $nl->get_variable_value($group, $var);
                # Need to strip the quotes
                $pathname =~ s/['"]//g;
                next if ($pathname =~ /UNSET$/);
                if ($input_pathname_type eq 'abs') {
                    if ($inputdata_rootdir) {
                        if ( $pathname !~ /^\s*$/ ) {   # If pathname isn't blank or null
                           print OUTFILE "$var = $pathname\n";
                        }
                    }
                    else {
                        if (-e $pathname) {  # use -e rather than -f since the absolute pathname
                                             # might be a directory
                            print "OK -- found $var = $pathname\n";
                        }
                        else {
                            print "NOT FOUND:  $var = $pathname\n";
                        }
                    }
                }
                elsif ($input_pathname_type =~ m/rel:(.+)/o) {
                    # The match provides the namelist variable that contains the
                    # root directory for a relative filename
                    my $rootdir_var = $1;
                    my $rootdir = $nl->get_variable_value($group, $rootdir_var);
                    $rootdir =~ s/['"]//g;
                    if ($inputdata_rootdir) {
                        $pathname = "$rootdir/$pathname";
                        #MV $pathname =~ s:$inputdata_rootdir::;
                        if ( $pathname !~ /^\s*$/ ) {   # If pathname isn't blank or null
                           print OUTFILE "$var = $pathname\n";
                        }
                    }
                    else {
                        if (-f "$rootdir/$pathname") {
                            print "OK -- found $var = $rootdir/$pathname\n";
                        }
                        else {
                            print "NOT FOUND:  $var = $rootdir/$pathname\n";
                        }
                    }
                }
            }
        }
    }
    close OUTFILE if defined $inputdata_rootdir;
    return 0 if defined $inputdata_rootdir;
}

#-------------------------------------------------------------------------------

sub set_abs_filepath {

# check whether the input filepath is an absolute path, and if it isn't then
# prepend a root directory

    my ($filepath, $rootdir) = @_;

    # strip any leading/trailing whitespace and quotes
    $filepath = trim($filepath);
    $filepath = remove_leading_and_trailing_quotes($filepath);
    $rootdir  = trim($rootdir);
    $rootdir = remove_leading_and_trailing_quotes($rootdir);

    my $out = $filepath;
    unless ( $filepath =~ /^\// ) {  # unless $filepath starts with a /
        $out = "$rootdir/$filepath"; # prepend the root directory
    }
    return $out;
}

#-------------------------------------------------------------------------------

sub valid_option {

    my ($val, @expect) = @_;

    my $expect;

    $val = trim($val);

    foreach $expect (@expect) {
        if ($val =~ /^$expect$/i) { return $expect; }
    }
    return undef;
}

#-------------------------------------------------------------------------------

sub check_use_case_name {
#
# Check the use-case name and ensure it follows the naming convention.
#
  my ($use_case) = @_;

  my $diestring = "bad use_case name $use_case, follow the conventions " .
                  "in namelist_files/use_cases/README\n";
  my $desc = "[a-zA-Z0-9]*";
  my $ssp_rcp  = "SSP[0-9]-[0-9\.]+";
  if (      $use_case =~ /^[0-9]+-[0-9]+([a-zA-Z0-9_\.-]*)_transient$/ ) {
    my $string = $1;
    if (      $string =~ /^_($ssp_rcp)_*($desc)$/ ) {
       # valid name
    } elsif ( $string =~ /^_*($desc)$/ ) {
       # valid name
    } else {
      $log->fatal_error($diestring);
    }
  } elsif ( $use_case =~ /^20thC([a-zA-Z0-9_\.]*)_transient$/ ) {
    my $string = $1;
    if (      $string =~ /^_($ssp_rcp)_*($desc)$/ ) {
       # valid name
    } elsif ( $string =~ /^_*($desc)$/ ) {
       # valid name
    } else {
      $log->fatal_error($diestring);
    }
  } elsif ( $use_case =~ /^([0-9]+|PI)-PD_*($desc)_transient$/   ) {
     # valid name
  } elsif ( $use_case =~ /^([0-9]+)_*($desc)_control$/   ) {
     # valid name
  } elsif ( $use_case =~ /^($desc)_pd$/   ) {
     # valid name
  } else {
      $log->fatal_error($diestring);
  }
}

#-------------------------------------------------------------------------------

sub validate_options {

# $source -- text string declaring the source of the options being validated
# $cfg    -- configure object
# $opts   -- reference to hash that contains the options

    my ($source, $cfg, $opts) = @_;

    my ($opt, $old, @expect);

    # use_case
    $opt = 'use_case';
    if (defined $opts->{$opt}) {

        if ( $opts->{$opt} ne "list" ) {
           # create the @expect array by listing the files in $use_case_dir
           # and strip off the ".xml" part of the filename
           @expect = ();
           my @files = bsd_glob("$opts->{'use_case_dir'}/*.xml");
           foreach my $file (@files) {
               $file =~ m{.*/(.*)\.xml};
               &check_use_case_name( $1 );
               push @expect, $1;
           }

           $old = $opts->{$opt};
           $opts->{$opt} = valid_option($old, @expect)
               or $log->fatal_error("invalid value of $opt ($old) specified in $source\n" .
                              "expected one of: @expect");
        } else {
           print "Use cases are:...\n\n";
           my @ucases;
           foreach my $file( sort( bsd_glob($opts->{'use_case_dir'}."/*.xml") ) ) {
              my $use_case;
              if ( $file =~ /\/([^\/]+)\.xml$/ ) {
                 &check_use_case_name( $1 );
                 $use_case = $1;
              } else {
                 $log->fatal_error("Bad name for use case file = $file");
              }
              my $uc_defaults = Build::NamelistDefaults->new("$file", $cfg);
              printf "%15s = %s\n", $use_case, $uc_defaults->get_value("use_case_desc");
              push @ucases, $use_case;
           }
           $log->exit_message("use cases : @ucases");
        }
    }
}

#-------------------------------------------------------------------------------

sub list_options {
#
# List the options for different command line values if asked for
#
    my ($opts_cmdl, $definition, $defaults) = @_;

    # options to list values that are in the defaults files
    my @opts_list = ( "res", "mask", "sim_year", "ssp_rcp" );
    my %opts_local;
    foreach my $var ( "res", "mask", "sim_year", "ssp_rcp" ) {
       my $val;
       if (      $opts_cmdl->{$var} eq "list" ) {
         $val = "default";
       } elsif ( $opts_cmdl->{$var} eq "default" ) {
         $val = $defaults->get_value($var, \%opts_local );
       } else {
         $val = $opts_cmdl->{$var};
       }
       my $vname = $var;
       if ( $vname eq "res" ) { $vname = "hgrid"; }
       $opts_local{$vname} = $val;
    }
    foreach my $var ( @opts_list ) {
       if (defined $opts_cmdl->{$var}) {

           if ( $opts_cmdl->{$var} eq "list" ) {
               my @valid_values   = $definition->get_valid_values( $var );
               if ( $var eq "sim_year" ) {
                   unshift( @valid_values,
                            $definition->get_valid_values( "sim_year_range" ) );
               }
               unshift( @valid_values, "default" );
               # Strip out quotes and the constant value
               for( my $i = 0; $i <= $#valid_values; $i++ ) {
                  $valid_values[$i] =~ s/('|')//g;
                  if ( $valid_values[$i] eq "constant" ) { $valid_values[$i] = undef; }
               }
               my $val= $defaults->get_value($var, \%opts_local);
               my $doc = $definition->get_var_doc( $var );
               $doc =~ s/\n//;
               chomp( $doc );
               $log->exit_message("valid values for $var ($doc) :\n" .
                            "    Values: @valid_values\n" .
                            "    Default = $val\n" .
                            "    (NOTE: resolution and mask and other settings may influence what the default is)");
           }
       }
    }
    # clm_demand
    my $var = 'clm_demand';
    if (defined $opts_cmdl->{$var}) {

        if ( $opts_cmdl->{$var} eq "list" ) {
           my @vars = $definition->get_var_names( );
           my @demands = ( "null" );
           foreach my $var ( @vars ) {
              if ( $definition->get_group_name( $var ) ne "clm_inparm" ) { next; }
              if ( defined($defaults->get_value($var, $opts_cmdl ) ) ) {
                 push( @demands, $var );
              }
           }
           my $doc = $definition->get_var_doc( 'clm_demand' );
           $doc =~ s/\n//;
           chomp( $doc );
           $log->exit_message("valid values for $var ($doc) :\n" .
                        "Namelist options to require: @demands\n" .
                        "any valid namelist item for clm_inparm can be set. However, not all are\n" .
                        "available in the clm defaults file. The defaults are also dependent on\n" .
                        "resolution and landmask, as well as other settings. Hence, the list above\n" .
                        "will vary depending on what you set for resolution and landmask.");
        }
    }
}

#-------------------------------------------------------------------------------

sub check_megan_spec {
#
# Check the megan specifier setting
#
    my ($opts, $nl, $definition) = @_;

    my $megan_spec      = $nl->get_value('megan_specifier');
    my @megan_spec_list = split( /\s*,\s*/, $megan_spec );
    foreach $megan_spec ( @megan_spec_list ) {
       if ( $megan_spec =~ /^['"]+[A-Za-z0-9]+\s*\=\s*([\sA-Za-z0-9+_-]+)["']+$/ ) {
          my $megan_list = $1;
          my @megan_cmpds = split( /\s*\+\s*/, $megan_list );
          my $var = "megan_cmpds";
          my $warn = 0;
          foreach my $megan_cmpd ( @megan_cmpds ) {
             if (  ! $definition->is_valid_value( $var, $megan_cmpd, 'noquotes'=>1 ) ) {
                $log->warning("megan_compound $megan_cmpd NOT found in list" );
                $warn++;
             }
          }
          if ( $warn > 0 ) {
             my @valid_values   = $definition->get_valid_values( $var, 'noquotes'=>1 );
             $log->warning("list of megan compounds includes:\n" .
                     "@valid_values\n" .
                     "Does your megan_factors_file include more compounds?\n" .
                     "If NOT your simulation will fail." );
          }
       } else {
          $log->fatal_error("Bad format for megan_specifier = $megan_spec");
       }
    }
}

#-------------------------------------------------------------------------------

sub trim {
   # remove leading and trailing whitespace from a string.
   my ($str) = @_;
   $str =~ s/^\s+//;
   $str =~ s/\s+$//;
   return $str;
}

#-------------------------------------------------------------------------------

sub quote_string {
   # Add quotes around a string, unless they are already there
   my ($str) = @_;
   $str = trim($str);
   unless ($str =~ /^['"]/) {        #"'
      $str = "\'$str\'";
   }
   return $str;
 }

#-------------------------------------------------------------------------------

sub remove_leading_and_trailing_quotes {
   # Remove leading and trailing single and double quotes from a string. Also
   # removes leading spaces before the leading quotes, and trailing spaces after
   # the trailing quotes.

   my ($str) = @_;

   $str = trim($str);

   # strip any leading/trailing quotes
   $str =~ s/^['"]+//;
   $str =~ s/["']+$//;

   return $str;
}

#-------------------------------------------------------------------------------

sub logical_to_fortran {
   # Given a logical variable ('true' / 'false'), convert it to a fortran-style logical ('.true.' / '.false.')
   # The result will be lowercase, regardless of the case of the input.
   my ($var) = @_;
   my $result;

   if (lc($var) eq 'true') {
      $result = ".true.";
   }
   elsif (lc($var) eq 'false') {
      $result = ".false.";
   }
   else {
      $log->fatal_error("Unexpected value in logical_to_fortran: $var");
   }

   return $result;
}

#-------------------------------------------------------------------------------

sub string_is_undef_or_empty {
   # Return true if the given string is undefined or only spaces, false otherwise.
   # A quoted empty string (' ' or " ") is treated as being empty.
   my ($str) = @_;
   if (!defined($str)) {
      return 1;
   }
   else {
      $str = remove_leading_and_trailing_quotes($str);
      if ($str =~ /^\s*$/) {
         return 1;
      }
      else {
         return 0;
      }
   }
}

#-------------------------------------------------------------------------------

sub value_is_true {
   # Return true if the given namelist value is .true.
   # An undefined value is treated as false (with the assumption that false is the default in the code)
   my ($val) = @_;

   # Some regular expressions...
   ###my $TRUE  = qr/\.true\./i;
   ###my $FALSE = qr/\.false\./i;
   # **N.B.** the use of qr// for precompiling regexps isn't supported until perl 5.005.
   my $TRUE  = '\.?true\.?|[t]';
   my $FALSE = '\.?false\.?|[f]';
   my $is_true = 0;
   if (defined($val)) {
      if ($val =~ /$TRUE/i) {
         $is_true = 1;
      }
   }

   return $is_true;
}

#-------------------------------------------------------------------------------

sub version {
# The version is found in CLM ChangeLog file.
# $cfgdir is set by the configure script to the name of its directory.

    my ($cfgdir) = @_;

    my $logfile = "$cfgdir/../doc/ChangeLog";

    my $fh = IO::File->new($logfile, '<') or $log->fatal_error("can't open ChangeLog file: $logfile");

    while (my $line = <$fh>) {

        if ($line =~ /^Tag name:\s*([a-zA-Z0-9_. -]*[clmcesm0-9_.-]+)$/ ) {
            $log->exit_message("$1");
        }
    }
}

#-------------------------------------------------------------------------------

sub main {
  my %nl_flags;
  $nl_flags{'cfgdir'} = dirname(abs_path($0));

  my %opts = process_commandline(\%nl_flags);
  my $cfgdir = $nl_flags{'cfgdir'};
  check_for_perl_utils($cfgdir, \%opts);

  $log     = namelist_files::LogMessages->new( $ProgName, \%opts );   # global
  version($cfgdir) if $opts{'version'};
  my $cfg = read_configure_definition($cfgdir, \%opts);

  my $physv      = config_files::clm_phys_vers->new( $cfg->get('phys') );
  my $definition = read_namelist_definition($cfgdir, \%opts, \%nl_flags);
  my $defaults   = read_namelist_defaults($cfgdir, \%opts, \%nl_flags, $cfg);

  # List valid values if asked for
  list_options(\%opts, $definition, $defaults);

  # Validate some of the commandline option values.
  validate_options("commandline", $cfg, \%opts);

  # Create an empty namelist object.
  my $nl = Build::Namelist->new();

  check_cesm_inputdata(\%opts, \%nl_flags);

  # Read in the env_*.xml files
  my %env_xml    = read_envxml_case_files( \%opts );

  # Process the user inputs
  process_namelist_user_input(\%opts, \%nl_flags, $definition, $defaults, $nl, $cfg, \%env_xml, $physv );
  # Get any other defaults needed from the namelist defaults file
  process_namelist_inline_logic(\%opts, \%nl_flags, $definition, $defaults, $nl, \%env_xml, $physv);

  # Validate that the entire resultant namelist is valid
  $definition->validate($nl);
  write_output_files(\%opts, \%nl_flags, $defaults, $nl);
  write_output_real_parameter_file(\%opts, \%nl_flags, $definition, $defaults, $nl);

  if ($opts{'inputdata'}) {
    check_input_files($nl, $nl_flags{'inputdata_rootdir'}, $opts{'inputdata'}, $definition);
  }
  $log->final_exit("Successfully made CLM namelist file");
}

#-------------------------------------------------------------------------------

1;
