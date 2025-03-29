#!/usr/bin/env perl

# Test command line options of the build-namelist script.
# Try to test that all the different options at least work.
# Test that inconsistentcies are appropriately caught.

#########################

use lib '.';
use Test::More;
use xFail::expectedFail;
use IO::File;

#########################

use strict;
use Getopt::Long;
use NMLTest::CompFiles;
use English;

sub usage {
    die <<EOF;
SYNOPSIS
     build-namelist_test.pl [options]

     Test the the CLM build-namelist
OPTIONS
     -help [or -h]                 Print usage to STDOUT.
     -compare <directory>          Compare namelists for this version to namelists
                                   created by another version.
     -generate                     Leave the namelists in place to do a later compare.
     -no-test                      Do NOT Use the -test option to make sure datasets exist.
     -csmdata "dir"                Root directory of CESM input data.

EOF
}

sub make_env_run {
#
# Create a env_run.xml file to read in
#
    my %settings = @_;

    # Set default settings
    my %env_vars = ( DIN_LOC_ROOT=>"MYDINLOCROOT", GLC_TWO_WAY_COUPLING=>"FALSE",  LND_SETS_DUST_EMIS_DRV_FLDS=>"TRUE", NEONSITE=>"", PLUMBER2SITE=>"" );
    # Set any settings that came in from function call
    foreach my $item ( keys(%settings) ) {
       $env_vars{$item} = $settings{$item};
    }

    # Now write the file out
    my $envfile = "env_run.xml";
    my $fh = IO::File->new($envfile, '>') or die "can't open file: $envfile";
    print $fh <<EOF;
<?xml version="1.0"?>

<config_definition>

EOF
    foreach my $item ( keys(%env_vars) ) {
      print $fh <<EOF;
<entry id="$item"         value="$env_vars{$item}"  />
EOF
    }
    print $fh <<EOF;

</config_definition>
EOF
    $fh->close();
}

sub make_config_cache {
   # Write a config_cache.xml file to read in
   my ($phys) = @_;
   my $config_cachefile = "config_cache.xml";
   my $fh = IO::File->new($config_cachefile, '>') or die "can't open file: $config_cachefile";
   print $fh <<EOF;
<?xml version="1.0"?>
<config_definition>
<commandline></commandline>
<entry id="phys" value="$phys" list="" valid_values="clm4_5,clm5_0,clm6_0">Specifies clm physics</entry>
</config_definition>
EOF
   $fh->close();
}

sub cat_and_create_namelistinfile {
#
# Concatenate the user_nl_clm files together and turn it into a namelist input file
# that can be read in by build-namelist
#
   my ($file1, $file2, $outfile) = @_;

   my $fh    = IO::File->new($file1,   '<') or die "can't open file: $file1";
   my $outfh = IO::File->new($outfile, '>') or die "can't open file: $outfile";
   print $outfh "&clm_settings\n\n";
   while ( my $line = <$fh> ) {
     print $outfh " $line";
   }
   $fh->close();
   if ( defined($file2) ) {
      my $fh    = IO::File->new($file2,   '<') or die "can't open file: $file2";
      while ( my $line = <$fh> ) {
        print $outfh " $line";
      }
   }
   print $outfh "\n/\n";
   $fh->close();
   $outfh->close();
}

#
# Process command-line options.
#
my %opts = ( help     => 0,
             generate => 0,
             test     => 1,
             compare  => undef,
             csmdata  => undef,
            );

GetOptions(
    "h|help"     => \$opts{'help'},
    "compare=s"  => \$opts{'compare'},
    "generate"   => \$opts{'generate'},
    "test!"      => \$opts{'test'},
    "csmdata=s"  => \$opts{'csmdata'},
)  or usage();

# Give usage message.
usage() if $opts{'help'};

# Check that the CESM inputdata root directory has been specified.  This must be
# a local or nfs mounted directory.
my $inputdata_rootdir = undef;
if (defined($opts{'csmdata'})) {
    $inputdata_rootdir = $opts{'csmdata'};
} elsif (defined $ENV{'CSMDATA'} ) {
    $inputdata_rootdir = $ENV{'CSMDATA'};
} else {
   # use yellowstone location as default
   $inputdata_rootdir="/glade/campaign/cesm/cesmdata/cseg/inputdata";
   print("WARNING:  -csmdata nor CSMDATA are set, using default yellowstone location: $inputdata_rootdir\n");
}

###################################
#_# read in expected fail test list
###################################
my $compGen;
if ( $opts{'generate'} eq 1 && !(defined($opts{'compare'}) )) {
   $compGen='generate';
} elsif ( defined($opts{'compare'}) ) {
   $compGen='compare';
} elsif ( defined($opts{'compare'} && ($opts{'generate'} eq 1 ))) {
   #_# if compare and generate are both given, use compare
   $compGen='compare';
}

my $ProgName;
($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!;
my $testType="namelistTest";

#
# Figure out number of tests that will run
#
my $ntests = 3264;

if ( defined($opts{'compare'}) ) {
   $ntests += 1980;
}
plan( tests=>$ntests );

#_# ============================================================
#_# setup for xFail module
#_# ============================================================
my $xFail = xFail::expectedFail->new($ProgName,$compGen,$ntests);
my $captOut="";  #_# variable to capture Test::More output
Test::More->builder->output(\$captOut);
#_# ============================================================
#_#
#_# ============================================================

# Check for unparsed arguments
if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
}
my $phys = "clm5_0";
my $mode = "-phys $phys";
&make_config_cache($phys);

my $DOMFILE = "$inputdata_rootdir/atm/datm7/domain.lnd.fv0.9x1.25_gx1v6.090309.nc";
my $bldnml = "../build-namelist -verbose -csmdata $inputdata_rootdir -configuration clm -structure standard -glc_nec 10 -no-note";
if ( $opts{'test'} ) {
   $bldnml .= " -test";
}

my $tempfile = "temp_file.txt";
if ( -f $tempfile ) {
  system( "/bin/rm $tempfile" );
}

my @files = ( "lnd_in", $tempfile );
my $cwd = `pwd`;
chomp( $cwd );
my $cfiles = NMLTest::CompFiles->new( $cwd, @files );

print "\n==================================================\n";
print "Run simple tests \n";
print "==================================================\n";

# Simple test -- just run build-namelist with -help option
eval{ system( "$bldnml -help > $tempfile 2>&1 " ); };
   is( $@, '', "help" );
   &cleanup();
# Simple test -- just run build-namelist with -version option
eval{ system( "$bldnml -version > $tempfile 2>&1 " ); };
   is( $@, '', "version" );
   system( "/bin/cat $tempfile" );
   &cleanup();
# Simple test -- just run build-namelist
&make_env_run();
eval{ system( "$bldnml > $tempfile 2>&1 " ); };
   is( $@, '', "plain build-namelist" );
   $cfiles->checkfilesexist( "default", $mode );
   # Compare to baseline
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "default", $mode );
      $cfiles->comparefiles( "default", $mode, $opts{'compare'} );
   }

print "\n==================================================\n";
print "Run simple tests with all list options \n";
print "==================================================\n";

$cfiles->copyfiles( "default", $mode );
&cleanup();
# Simple test -- run all the list options
foreach my $options ( "clm_demand", "ssp_rcp",      "res",
                      "sim_year",   "use_case" ) {
   &make_env_run();
   eval{ system( "$bldnml -${options} list > $tempfile 2>&1 " ); };
   my $result = `cat $tempfile`;
   my $expect;
   if ( $options =~ /use_case/ ) {
      $expect = "use cases :";
   } else {
      $expect = "valid values for $options";
   }
   $expect    = "/CLM build-namelist : $expect/";
   like( $result, $expect, "$options list" );
   is( (-f "lnd_in"), undef, "Check that lnd_in file does NOT exist" );
   &cleanup();
}

print "\n==================================================\n";
print "Run simple tests with additional options \n";
print "==================================================\n";

# Exercise a bunch of options
my $options = "-co2_ppmv 250 ";
   $options .= " -res 10x15 -ssp_rcp SSP2-4.5 -envxml_dir .";

   &make_env_run();
   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "options: $options" );
      $cfiles->checkfilesexist( "default", $mode );
      $cfiles->copyfiles( "most_options", $mode );
   # Compare to default
      $cfiles->doNOTdodiffonfile( "lnd_in",    "default", $mode );
      $cfiles->doNOTdodiffonfile( "$tempfile", "default", $mode );
      $cfiles->comparefiles( "default", $mode );
   # Compare to baseline
   if ( defined($opts{'compare'}) ) {
      $cfiles->dodiffonfile(      "lnd_in",    "most_options", $mode );
      $cfiles->doNOTdodiffonfile( "$tempfile", "most_options", $mode );
      $cfiles->comparefiles( "most_options", $mode, $opts{'compare'} );
   }
   &cleanup();

print "\n==================================================\n";
print "Test drydep, fire_emis and megan namelists  \n";
print "==================================================\n";

# drydep and megan namelists
$phys = "clm5_0";
$mode = "-phys $phys";
&make_config_cache($phys);
my @mfiles = ( "lnd_in", "drv_flds_in", $tempfile );
my $mfiles = NMLTest::CompFiles->new( $cwd, @mfiles );
foreach my $options ( "-drydep", "-megan", "-drydep -megan", "-fire_emis", "-drydep -megan -fire_emis" ) {
   &make_env_run();
   eval{ system( "$bldnml -envxml_dir . $options > $tempfile 2>&1 " ); };
   is( $@, '', "options: $options" );
   $mfiles->checkfilesexist( "$options", $mode);
   if ( $options ne "-drydep" ) {
     $mfiles->shownmldiff( "-drydep", $mode );
   }
   if ( defined($opts{'compare'}) ) {
      $mfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $mfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $mfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}
$phys = "clm5_0";
$mode = "-phys $phys";
&make_config_cache($phys);

print "\n===============================================================================\n";
print "Test configuration, structure, irrigate, verbose, clm_demand, ssp_rcp, test, sim_year, use_case\n";
print "=================================================================================\n";

my $startfile = "clmrun.clm2.r.1964-05-27-00000.nc";
foreach my $driver ( "nuopc" ) {
   print "   For $driver driver\n\n";
   # configuration, structure, irrigate, verbose, clm_demand, ssp_rcp, test, sim_year, use_case
   foreach my $options ( "-res 0.9x1.25 -configuration nwp",
                         "-res 0.9x1.25 -structure fast",
                         "-res 0.9x1.25 -namelist '&a irrigate=.true./'", "-res 0.9x1.25 -verbose", "-res 0.9x1.25 -ssp_rcp SSP2-4.5", "-res 0.9x1.25 -test", "-res 0.9x1.25 -sim_year 1850",
                         "-res 0.9x1.25 -namelist '&a use_lai_streams=.true.,use_soil_moisture_streams=.true./'",
                         "-res 0.9x1.25 -namelist '&a use_excess_ice=.true. use_excess_ice_streams=.true./'",
                         "-res 0.9x1.25 --clm_start_type cold -namelist '&a use_excess_ice=.true. use_excess_ice_streams=.true./'",
                         "-res 0.9x1.25 -use_case 1850_control",
                         "-res 1x1pt_US-UMB -clm_usr_name 1x1pt_US-UMB -namelist '&a fsurdat=\"/dev/null\"/'",
                         "-res 1x1_brazil",
                         "-namelist '&a use_matrixcn=F,use_soil_matrixcn=F,hist_wrt_matrixcn_diag=F,spinup_matrixcn=F/' -bgc sp",
                         "-namelist '&a use_matrixcn=T,use_soil_matrixcn=T,hist_wrt_matrixcn_diag=T,spinup_matrixcn=T/' -bgc bgc -crop -clm_accelerated_spinup on",
                         "-namelist \"&a soil_decomp_method='MIMICSWieder2015',use_matrixcn=F/\" -bgc bgc -crop",
                         "-namelist \"&a soil_decomp_method='MIMICSWieder2015',use_matrixcn=T/\" -bgc bgc -crop",
                         "-bgc bgc -crop -clm_accelerated_spinup sasu",
                         "-res 0.9x1.25 -clm_start_type startup", "-namelist '&a irrigate=.false./' -crop -bgc bgc",
                         "-res 0.9x1.25 -infile myuser_nl_clm",
                         "-res 0.9x1.25 -ignore_ic_date -clm_start_type branch -namelist '&a nrevsn=\"thing.nc\"/' -bgc bgc -crop",
                         "-res 0.9x1.25 -clm_start_type branch -namelist '&a nrevsn=\"thing.nc\",use_init_interp=T/'",
                         "-res 0.9x1.25 -ignore_ic_date -clm_start_type startup -namelist '&a finidat=\"thing.nc\"/' -bgc bgc -crop",
                        ) {
      my $file = $startfile;
      &make_env_run();
      my $base_options = "-envxml_dir . -driver $driver";
      if ( $driver eq "nuopc" ) {
         $base_options = "$base_options -namelist '&a force_send_to_atm = .false./'";
      }
      eval{ system( "$bldnml $base_options $options > $tempfile 2>&1 " ); };
      is( $@, '', "options: $base_options $options" );
      $cfiles->checkfilesexist( "$base_options $options", $mode );
      $cfiles->shownmldiff( "default", $mode );
      my $finidat = `grep finidat lnd_in`;
      if ( $options =~ /myuser_nl_clm/ ) {
         my $fsurdat =  `grep fsurdat lnd_in`;
         like( $fsurdat, "/MYDINLOCROOT/lnd/clm2/PTCLMmydatafiles/1x1pt_US-UMB/surfdata_1x1pt_US-UMB_simyr2000_clm4_5_c131122.nc/", "$options" );
      }
      if ( defined($opts{'compare'}) ) {
         $cfiles->doNOTdodiffonfile( "$tempfile", "$base_options $options", $mode );
         $cfiles->comparefiles( "$base_options $options", $mode, $opts{'compare'} );
      }
      if ( defined($opts{'generate'}) ) {
         $cfiles->copyfiles( "$base_options $options", $mode );
      }
      &cleanup();
   }
}
print "\n===============================================================================\n";
print "Test the NEON sites\n";
print "=================================================================================\n";
my $phys = "clm6_0";
$mode = "-phys $phys";
&make_config_cache($phys);
my $neondir      = "../../cime_config/usermods_dirs/clm/NEON";
foreach my $site ( "ABBY", "BLAN", "CPER", "DEJU", "GRSM", "HEAL", "KONA", "LENO", "NIWO",
                   "ONAQ", "PUUM", "SERC", "SRER", "TALL", "TREE", "WOOD", "BARR", "BONA",
                   "DCFS", "DELA", "GUAN", "JERC", "KONZ", "MLBS", "NOGP", "ORNL", "RMNP",
                   "SJER", "STEI", "TEAK", "UKFS", "WREF", "BART", "CLBJ", "DSNY", "HARV",
                   "JORN", "LAJA", "MOAB", "OAES", "OSBS", "SCBI", "SOAP", "STER", "TOOL",
                   "UNDE", "YELL"
 ) {
   &make_env_run( NEONSITE=>"$site" );
   #
   # Concatonate  default usermods and specific sitetogether expanding env variables while doing that
   #
   if ( ! -d "$neondir/$site" ) {
      print "NEON directory is not there: $neondir/$site\n";
      die "ERROR:: NEON site does not exist: $site\n";
   }
   my $neondefaultfile = "$neondir/defaults/user_nl_clm";
   my $neonsitefile = "$neondir/$site/user_nl_clm";
   if ( ! -f $neonsitefile )  {
      $neonsitefile = undef;
   }
   $ENV{'NEONSITE'} = $site;
   my $namelistfile = "temp.namelistinfile_$site";
   &cat_and_create_namelistinfile( $neondefaultfile, $neonsitefile, $namelistfile );
   #
   # Now run  the site for both bgc and non-FATES
   #
   foreach my $bgc ( "bgc", "fates") {
      if ( ($bgc eq "bgc") or ($site ne "STER" and $site ne "KONA")) {
         my $options = "--res CLM_USRDAT --clm_usr_name NEON --no-megan --bgc $bgc --use_case 2018_control --infile $namelistfile";
         eval{ system( "$bldnml -envxml_dir . $options > $tempfile 2>&1 " ); };
         is( $@, '', "options: $options" );
         $cfiles->checkfilesexist( "$options", $mode );
         $cfiles->shownmldiff( "default", $mode );
         if ( defined($opts{'compare'}) ) {
            $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
            $cfiles->dodiffonfile(      "lnd_in",    "$options", $mode );
            $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
         }
         if ( defined($opts{'generate'}) ) {
            $cfiles->copyfiles( "$options", $mode );
         }
      }
   }
   system( "/bin/rm $namelistfile" );
   &cleanup();
}
print "\n===============================================================================\n";
print "Test the PLUMBER2 sites\n";
print "=================================================================================\n";
my $phys = "clm6_0";
$mode = "-phys $phys";
&make_config_cache($phys);
my $plumdir      = "../../cime_config/usermods_dirs/clm/PLUMBER2";
foreach my $site ( 
    "AR-SLu",  "AU-Emr",  "AU-TTE",  "CA-NS1",  "CA-SF3",  "CN-HaM",  "DE-Obe",  "ES-ES1",  "FR-Gri",  "IE-Dri",  "IT-LMa",  "IT-SRo",  "RU-Fyo",  "US-Aud",  "US-Ho1",  "US-Ne2",  "US-Syv",  "ZM-Mon",
    "AT-Neu",  "AU-Gin",  "AU-Tum",  "CA-NS2",  "CH-Cha",  "CN-Qia",  "DE-Seh",  "ES-ES2",  "FR-Hes",  "IT-Amp",  "IT-Mal",  "JP-SMF",  "RU-Zot",  "US-Bar",  "US-KS2",  "US-Ne3",  "US-Ton",
    "AU-ASM",  "AU-GWW",  "AU-Whr",  "CA-NS4",  "CH-Dav",  "CZ-wet",  "DE-SfN",  "ES-LgS",  "FR-LBr",  "IT-BCi",  "IT-MBo",  "NL-Ca1",  "SD-Dem",  "US-Bkg",  "US-Los",  "US-NR1",  "US-Tw4",
    "AU-Cow",  "AU-How",  "AU-Wrr",  "CA-NS5",  "CH-Fru",  "DE-Bay",  "DE-Tha",  "ES-LMa",  "FR-Lq1",  "IT-CA1",  "IT-Noe",  "NL-Hor",  "SE-Deg",  "US-Blo",  "US-Me2",  "US-PFa",  "US-Twt",
    "AU-Cpr",  "AU-Lit",  "AU-Ync",  "CA-NS6",  "CH-Oe1",             "DE-Wet",  "ES-VDA",  "FR-Lq2",  "IT-CA2",  "IT-Non",  "NL-Loo",  "UK-Gri",  "US-Bo1",  "US-Me4",  "US-Prr",  "US-UMB",
    "AU-Ctr",  "AU-Otw",  "BE-Bra",  "CA-NS7",  "CN-Cha",  "DE-Geb",  "DK-Fou",  "FI-Hyy",  "FR-Pue",  "IT-CA3",  "IT-PT1",  "PL-wet",  "UK-Ham",  "US-Cop",  "US-Me6",  "US-SP1",  "US-Var",
    "AU-Cum",  "AU-Rig",  "BE-Lon",  "CA-Qcu",  "CN-Cng",  "DE-Gri",  "DK-Lva",  "FI-Kaa",  "GF-Guy",  "IT-Col",  "IT-Ren",  "PT-Esp",  "UK-PL3",  "US-FPe",  "US-MMS",  "US-SP2",  "US-WCr",
    "AU-DaP",  "AU-Rob",  "BE-Vie",  "CA-Qfo",  "CN-Dan",  "DE-Hai",  "DK-Ris",  "FI-Lom",  "HU-Bug",  "IT-Cpz",  "IT-Ro1",  "PT-Mi1",  "US-AR1",  "US-GLE",  "US-MOz",  "US-SP3",  "US-Whs",
    "AU-DaS",  "AU-Sam",  "BR-Sa3",  "CA-SF1",  "CN-Din",  "DE-Kli",  "DK-Sor",  "FI-Sod",  "ID-Pag",  "IT-Isp",  "IT-Ro2",  "PT-Mi2",  "US-AR2",  "US-Goo",  "US-Myb",  "US-SRG",  "US-Wkg",
    "AU-Dry",  "AU-Stp",  "BW-Ma1",  "CA-SF2",  "CN-Du2",  "DE-Meh",  "DK-ZaH",  "FR-Fon",  "IE-Ca1",  "IT-Lav",  "IT-SR2",  "RU-Che",  "US-ARM",  "US-Ha1",  "US-Ne1",  "US-SRM",  "ZA-Kru"
 ) {
   &make_env_run( PLUMBER2SITE=>"$site" );
   #
   # Concatonate  default usermods and specific sitetogether expanding env variables while doing that
   #
   if ( ! -d "$plumdir/$site" ) {
      print "PLUMBER2 directory is not there: $plumdir/$site\n";
      die "ERROR:: PLUMBER2 site does not exist: $site\n";
   }
   my $plumdefaultfile = "$plumdir/defaults/user_nl_clm";
   my $plumsitefile = "$plumdir/$site/user_nl_clm";
   if ( ! -f $plumsitefile )  {
      $plumsitefile = undef;
   }
   $ENV{'PLUMBER2'} = $site;
   my $namelistfile = "temp.namelistinfile_$site";
   &cat_and_create_namelistinfile( $plumdefaultfile, $plumsitefile, $namelistfile );
   #
   # Now run  the site
   #
   my $options = "--res CLM_USRDAT --clm_usr_name PLUMBER2 --no-megan --bgc sp --infile $namelistfile";
   eval{ system( "$bldnml -envxml_dir . $options > $tempfile 2>&1 " ); };
   is( $@, '', "options: $options" );
   $cfiles->checkfilesexist( "$options", $mode );
   $cfiles->shownmldiff( "default", $mode );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->dodiffonfile(      "lnd_in",    "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   system( "/bin/rm $namelistfile" );
   &cleanup();
}

print "\n===============================================================================\n";
print "Test some CAM specific setups for special grids \n";
print "=================================================================================\n";
foreach my $phys ( "clm4_5", "clm5_0", "clm6_0" ) {
   $mode = "-phys $phys";
   &make_config_cache($phys);
   foreach my $options (
                      "-res ne0np4.ARCTIC.ne30x4 -bgc sp -use_case 2000_control -namelist '&a start_ymd=19790101/' -lnd_tuning_mode ${phys}_cam7.0",
                      "-res ne0np4.ARCTICGRIS.ne30x8 -bgc sp -use_case 1850_control -namelist '&a start_ymd=19790101/' -lnd_tuning_mode ${phys}_cam7.0",
                      "-res 1.9x2.5 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=19790101/' -lnd_tuning_mode ${phys}_cam7.0",
                      "-res 0.9x1.25 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=19790101/' -lnd_tuning_mode ${phys}_cam7.0",
                      "-res 0.9x1.25 -bgc bgc -crop -use_case 20thC_transient -namelist '&a start_ymd=19500101/' -lnd_tuning_mode ${phys}_cam7.0",
                      "-res ne0np4CONUS.ne30x8 -bgc sp -use_case 2000_control  -namelist '&a start_ymd=20130101/' -lnd_tuning_mode ${phys}_cam7.0",
                      "-res 1.9x2.5 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=20030101/' -lnd_tuning_mode ${phys}_cam7.0",
                      "-res 1.9x2.5 -bgc sp -use_case 2010_control -namelist '&a start_ymd=20100101/' -lnd_tuning_mode ${phys}_cam7.0",
                      "-res 1x1_brazil -no-megan -use_case 2000_control -lnd_tuning_mode ${phys}_GSWP3v1",
                      "-res C96 -bgc sp -use_case 2010_control -namelist '&a start_ymd=20100101/' -lnd_tuning_mode ${phys}_cam7.0",
                      "-res ne0np4.ARCTIC.ne30x4 -bgc sp -use_case 2000_control -namelist '&a start_ymd=20130101/' -lnd_tuning_mode ${phys}_cam7.0",
                     ) {
      &make_env_run();
      eval{ system( "$bldnml -envxml_dir . $options > $tempfile 2>&1 " ); };
      is( $@, '', "options: $options" );
      $cfiles->checkfilesexist( "$options", $mode );
      $cfiles->shownmldiff( "default", $mode );
      if ( defined($opts{'compare'}) ) {
         $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
         $cfiles->dodiffonfile(      "lnd_in",    "$options", $mode );
         $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
      }
      if ( defined($opts{'generate'}) ) {
         $cfiles->copyfiles( "$options", $mode );
      }
      &cleanup();
   }
}
print "\n===============================================================================\n";
print "Test setting drv_flds_in fields in CAM, clm60 only";
print "=================================================================================\n";
foreach my $phys ( "clm6_0" ) {
   $mode = "-phys $phys CAM_SETS_DRV_FLDS";
   &make_config_cache($phys);
   foreach my $options (
                      "--res ne0np4.POLARCAP.ne30x4 --mask tx0.1v2 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=19790101/' -lnd_tuning_mode ${phys}_cam7.0 --infile empty_user_nl_clm",
                     ) {
      &make_env_run( 'LND_SETS_DUST_EMIS_DRV_FLDS'=>"FALSE" );
      eval{ system( "$bldnml --envxml_dir . $options > $tempfile 2>&1 " ); };
      is( $@, '', "options: $options" );
      $cfiles->checkfilesexist( "$options", $mode );
      $cfiles->shownmldiff( "default", $mode );
      if ( defined($opts{'compare'}) ) {
         $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
         $cfiles->dodiffonfile(      "lnd_in",    "$options", $mode );
         $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
      }
      if ( defined($opts{'generate'}) ) {
         $cfiles->copyfiles( "$options", $mode );
      }
      &cleanup();
   }
}
print "\n===============================================================================\n";
print "Test setting drv_flds_in fields in CAM";
print "=================================================================================\n";
foreach my $phys ( "clm5_0", "clm6_0" ) {
   $mode = "-phys $phys CAM_SETS_DRV_FLDS";
   &make_config_cache($phys);
   foreach my $options (
                      "--res 1.9x2.5 --mask gx1v7 --bgc sp --use_case 20thC_transient --namelist '&a start_ymd=19790101/' --lnd_tuning_mode ${phys}_cam6.0 --infile empty_user_nl_clm",
                      "--res 1.9x2.5 --mask gx1v7 --bgc sp --use_case 20thC_transient --namelist '&a start_ymd=19790101/' --lnd_tuning_mode ${phys}_cam7.0 --infile empty_user_nl_clm",
                      "--res 1.9x2.5 --mask gx1v7 --bgc sp -no-crop --use_case 20thC_transient --namelist '&a start_ymd=19790101/' --lnd_tuning_mode ${phys}_cam7.0 --infile empty_user_nl_clm",
                      "--res ne0np4.ARCTIC.ne30x4 --mask tx0.1v2 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=19790101/' -lnd_tuning_mode ${phys}_cam7.0 --infile empty_user_nl_clm",
                      "--res ne0np4.ARCTICGRIS.ne30x8 --mask tx0.1v2 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=19790101/' -lnd_tuning_mode ${phys}_cam7.0 --infile empty_user_nl_clm",
                      "--res ne0np4CONUS.ne30x8 --mask tx0.1v2 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=20130101/' -lnd_tuning_mode ${phys}_cam7.0 --infile empty_user_nl_clm",
                     ) {
      &make_env_run( 'LND_SETS_DUST_EMIS_DRV_FLDS'=>"FALSE" );
      eval{ system( "$bldnml --envxml_dir . $options > $tempfile 2>&1 " ); };
      is( $@, '', "options: $options" );
      $cfiles->checkfilesexist( "$options", $mode );
      $cfiles->shownmldiff( "default", $mode );
      if ( defined($opts{'compare'}) ) {
         $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
         $cfiles->dodiffonfile(      "lnd_in",    "$options", $mode );
         $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
      }
      if ( defined($opts{'generate'}) ) {
         $cfiles->copyfiles( "$options", $mode );
      }
      &cleanup();
   }
}
print "\n==============================================================\n";
print "Test several use_cases and specific configurations for clm5_0\n";
print "==============================================================\n";
$phys = "clm5_0";
$mode = "-phys $phys";
&make_config_cache($phys);
foreach my $options (
                      "--res 0.9x1.25 --bgc sp  --use_case 1850-2100_SSP2-4.5_transient --namelist '&a start_ymd=18501223/'",
                      "-bgc fates  -use_case 2000_control -no-megan",
                      "-bgc fates  -use_case 20thC_transient -no-megan",
                      "-bgc fates  -use_case 20thC_transient -no-megan -no-crop --res 4x5",
                      "-bgc fates  -use_case 1850_control -no-megan -namelist \"&a use_fates_sp=T, soil_decomp_method='None'/\"",
                      "-bgc sp  -use_case 2000_control -res 0.9x1.25 -namelist '&a use_soil_moisture_streams = T/'",
                      "--res 1.9x2.5 --bgc bgc --use_case 1850-2100_SSP2-4.5_transient --namelist '&a start_ymd=19101023/'",
                      "-namelist \"&a dust_emis_method='Zender_2003', zender_soil_erod_source='lnd' /'\"",
                      "-bgc bgc -use_case 2000_control -namelist \"&a fire_method='nofire'/\" -crop",
                      "-res 0.9x1.25 -bgc sp -use_case 1850_noanthro_control -drydep -fire_emis",
                      "-res 0.9x1.25 -bgc bgc -use_case 1850_noanthro_control -drydep -fire_emis -light_res 360x720",
                      "--bgc bgc --light_res none --namelist \"&a fire_method='nofire'/\"",
                      "--bgc fates --light_res 360x720 --no-megan --namelist \"&a fates_spitfire_mode=2/\"",
                      "--bgc fates --light_res none --no-megan --namelist \"&a fates_spitfire_mode=1/\"",
                     ) {
   my $file = $startfile;
   &make_env_run();
   eval{ system( "$bldnml -envxml_dir . $options > $tempfile 2>&1 " ); };
   is( $@, '', "options: $options" );
   $cfiles->checkfilesexist( "$options", $mode );
   $cfiles->shownmldiff( "default", $mode );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->dodiffonfile(      "lnd_in",    "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}



print "\n==================================================\n";
print "Start Failure testing.  These should fail \n";
print "==================================================\n";

# Failure testing, do things that SHOULD fail
my $finidat  = "thing.nc";
system( "touch $finidat" );

my %failtest = (
     "coldstart but with IC file"=>{ options=>"-clm_start_type cold -envxml_dir .",
                                     namelst=>"finidat='$finidat'",
                                     phys=>"clm5_0",
                                   },
     "clm_demand on finidat"     =>{ options=>"-clm_demand finidat -envxml_dir .",
                                     namelst=>"",
                                     phys=>"clm5_0",
                                   },
     "blank IC file, not cold"   =>{ options=>"-clm_start_type startup -envxml_dir .",
                                     namelst=>"finidat=' '",
                                     phys=>"clm5_0",
                                   },
     "startup without interp"    =>{ options=>"-clm_start_type startup -envxml_dir . -bgc sp -sim_year 1850",
                                     namelst=>"use_init_interp=.false., start_ymd=19200901",
                                     phys=>"clm5_0",
                                   },
     "use_crop without -crop"    =>{ options=>" -envxml_dir .",
                                     namelst=>"use_crop=.true.",
                                     phys=>"clm4_5",
                                   },
     "LeungDust_WO_Prigent"      =>{ options=>" -envxml_dir . -bgc sp",
                                     namelst=>"use_prigent_roughness=.false.,dust_emis_method='Leung_2023'",
                                     phys=>"clm6_0",
                                   },
     "soilm_stream off w file"      =>{ options=>"-res 0.9x1.25 -envxml_dir .",
                                     namelst=>"use_soil_moisture_streams = .false.,stream_fldfilename_soilm='file_provided_when_off'",
                                     phys=>"clm5_0",
                                   },
     "exice_stream off w file"  =>{ options=>"-res 0.9x1.25 -envxml_dir .",
                                     namelst=>"use_excess_ice=.true., use_excess_ice_streams = .false.,stream_fldfilename_exice='file_provided_when_off'",
                                     phys=>"clm5_0",
                                   },
     "exice_stream off w mesh"  =>{ options=>"-res 0.9x1.25 -envxml_dir .",
                                     namelst=>"use_excess_ice=.true., use_excess_ice_streams = .false.,stream_meshfile_exice='file_provided_when_off'",
                                     phys=>"clm5_0",
                                   },
     "exice off, but stream on"  =>{ options=>"-res 0.9x1.25 -envxml_dir .",
                                     namelst=>"use_excess_ice=.false., use_excess_ice_streams = .true.,stream_fldfilename_exice='file_provided', stream_meshfile_exice='file_provided'",
                                     phys=>"clm5_0",
                                   },
     "exice stream off, but setmap"=>{ options=>"-res 0.9x1.25 -envxml_dir .",
                                     namelst=>"use_excess_ice=.true., use_excess_ice_streams = .false.,stream_mapalgo_exice='bilinear'",
                                     phys=>"clm5_0",
                                   },
     "coldstart exice on wo stream"=>{ options=>"-res 0.9x1.25 -envxml_dir . --clm_start_type cold",
                                     namelst=>"use_excess_ice=.true., use_excess_ice_streams = .false.",
                                     phys=>"clm6_0",
                                   },
     "coldstart exice on bad temp" =>{ options=>"-res 0.9x1.25 -envxml_dir . --clm_start_type cold",
                                     namelst=>"use_excess_ice=.true., use_excess_ice_streams = .true., excess_ice_coldstart_temp=0.0",
                                     phys=>"clm6_0",
                                   },
     "coldstart exice on bad depth" =>{ options=>"-res 0.9x1.25 -envxml_dir . --clm_start_type cold",
                                     namelst=>"use_excess_ice=.true., use_excess_ice_streams = .true., excess_ice_coldstart_depth=0.0",
                                     phys=>"clm6_0",
                                   },
     "clm50CNDVwtransient"       =>{ options=>" -envxml_dir . -use_case 20thC_transient -dynamic_vegetation -res 10x15 -ignore_warnings",
                                     namelst=>"",
                                     phys=>"clm5_0",
                                   },
     "decomp_without_cn"         =>{ options=>" -envxml_dir . -bgc sp",
                                     namelst=>"soil_decomp_method='CENTURYKoven2013'",
                                     phys=>"clm5_0",
                                   },
     "bgc_with_no_decomp"        =>{ options=>" -envxml_dir . -bgc bgc",
                                     namelst=>"soil_decomp_method='None'",
                                     phys=>"clm5_0",
                                   },
     "reseed without CN"         =>{ options=>" -envxml_dir . -bgc sp",
                                     namelst=>"reseed_dead_plants=.true.",
                                     phys=>"clm5_0",
                                   },
     "onset_threh w SP"          =>{ options=>" -envxml_dir . -bgc sp",
                                     namelst=>"onset_thresh_depends_on_veg=.true.",
                                     phys=>"clm6_0",
                                   },
     "dribble_crphrv w/o CN"     =>{ options=>" -envxml_dir . -bgc sp",
                                     namelst=>"dribble_crophrv_xsmrpool_2atm=.true.",
                                     phys=>"clm5_0",
                                   },
     "dribble_crphrv w/o crop"   =>{ options=>" -envxml_dir . -bgc bgc -no-crop",
                                     namelst=>"dribble_crophrv_xsmrpool_2atm=.true.",
                                     phys=>"clm5_0",
                                   },
     "CNDV with flanduse_timeseries - clm4_5"=>{ options=>"-bgc bgc -dynamic_vegetation -envxml_dir . -ignore_warnings",
                                     namelst=>"flanduse_timeseries='my_flanduse_timeseries_file.nc'",
                                     phys=>"clm4_5",
                                   },
     "use_cndv=T without bldnml op"=>{ options=>"-bgc bgc -envxml_dir . -ignore_warnings",
                                     namelst=>"use_cndv=.true.",
                                     phys=>"clm4_5",
                                   },
     "use_cndv=F with dyn_veg op"=>{ options=>"-bgc bgc -dynamic_vegetation -envxml_dir . -ignore_warnings",
                                     namelst=>"use_cndv=.false.",
                                     phys=>"clm4_5",
                                   },
     "crop with use_crop false"  =>{ options=>"-crop -bgc bgc -envxml_dir .",
                                     namelst=>"use_crop=.false.",
                                     phys=>"clm4_5",
                                   },
     "crop without CN"           =>{ options=>"-crop -bgc sp -envxml_dir .",
                                     namelst=>"",
                                     phys=>"clm4_5",
                                   },
     "toosmall soil w trans"     =>{ options=>"-envxml_dir .",
                                     namelst=>"toosmall_soil=10, dyn_transient_pfts=T",
                                     phys=>"clm5_0",
                                   },
     "toosmall lake w trans"     =>{ options=>"-envxml_dir .",
                                     namelst=>"toosmall_lake=10, dyn_transient_pfts=T",
                                     phys=>"clm5_0",
                                   },
     "toosmall crop w trans"     =>{ options=>"-bgc bgc -crop -envxml_dir .",
                                     namelst=>"toosmall_crop=10, dyn_transient_pfts=T",
                                     phys=>"clm5_0",
                                   },
     "toosmall wetl w trans"     =>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"toosmall_wetland=10, dyn_transient_pfts=T",
                                     phys=>"clm5_0",
                                   },
     "toosmall glc  w trans"     =>{ options=>"-bgc sp  -envxml_dir .",
                                     namelst=>"toosmall_glacier=10, dyn_transient_pfts=T",
                                     phys=>"clm5_0",
                                   },
     "toosmall urban w trans"    =>{ options=>"-bgc sp  -envxml_dir .",
                                     namelst=>"toosmall_urban=10, dyn_transient_pfts=T",
                                     phys=>"clm5_0",
                                   },
     "collapse_urban w trans"    =>{ options=>"-bgc sp  -envxml_dir .",
                                     namelst=>"collapse_urban=T, dyn_transient_crops=T",
                                     phys=>"clm5_0",
                                   },
     "n_dom_landunits w trans"    =>{ options=>"-bgc sp  -envxml_dir .",
                                     namelst=>"n_dom_landunits=2, dyn_transient_crops=T",
                                     phys=>"clm5_0",
                                   },
     "n_dom_pfts w trans"         =>{ options=>"-bgc sp  -envxml_dir .",
                                     namelst=>"n_dom_pfts=2, dyn_transient_crops=T",
                                     phys=>"clm5_0",
                                   },
     "baset_map without crop"    =>{ options=>"-bgc bgc -envxml_dir . -no-crop",
                                     namelst=>"baset_mapping='constant'",
                                     phys=>"clm5_0",
                                   },
     "mapvary var w/o varymap"   =>{ options=>"-crop -bgc bgc -envxml_dir . -crop",
                                     namelst=>"baset_mapping='constant', baset_latvary_slope=1.0, baset_latvary_intercept=10.0",
                                     phys=>"clm5_0",
                                   },
     "grainproductWOcrop"       =>{ options=>"-bgc bgc -no-crop -envxml_dir .",
                                    namelst=>"use_grainproduct=.true.",
                                    phys=>"clm4_5",
                                   },
     "interp without finidat"    =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"use_init_interp=.true. finidat=' '",
                                     phys=>"clm5_0",
                                   },
     "sp and c13"                =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"use_c13=.true.",
                                     phys=>"clm4_5",
                                   },
     "sp and c14"                =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"use_c14=.true.",
                                     phys=>"clm4_5",
                                   },
     "bombspike no c14"          =>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"use_c14=.false. use_c14_bombspike=.true.",
                                     phys=>"clm5_0",
                                   },
     "use c13 timeseries no cn"  =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"use_c13_timeseries=.true.",
                                     phys=>"clm4_5",
                                   },
     "use c13 timeseries no c13"=>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"use_c13=.false. use_c13_timeseries=.true.",
                                     phys=>"clm4_5",
                                   },
     "bombspike no cn"           =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"use_c14_bombspike=.true.",
                                     phys=>"clm5_0",
                                   },
     "lightres no cn"            =>{ options=>"-bgc sp -envxml_dir . -light_res 360x720",
                                     namelst=>"",
                                     phys=>"clm5_0",
                                   },
     "NEONlightresButGlobal"     =>{ options=>"--res 4x5 --bgc bgc --envxml_dir . --light_res 106x740",
                                     namelst=>"",
                                     phys=>"clm6_0",
                                   },
     "spno-fire"                 =>{ options=>"-bgc sp -envxml_dir . -use_case 2000_control",
                                     namelst=>"fire_method='nofire'",
                                     phys=>"clm5_0",
                                   },
     "lightres no fire"          =>{ options=>"-bgc bgc -envxml_dir . -light_res 360x720",
                                     namelst=>"fire_method='nofire'",
                                     phys=>"clm5_0",
                                   },
     "lightres none bgc"         =>{ options=>"-bgc bgc -envxml_dir . -light_res none",
                                     namelst=>"",
                                     phys=>"clm5_0",
                                   },
     "lightresnotnone-nofire"    =>{ options=>"-bgc bgc -envxml_dir . -light_res 94x192",
                                     namelst=>"fire_method='nofire'",
                                     phys=>"clm5_0",
                                   },
     "lightresnonenofirelightfil"=>{ options=>"-bgc bgc -envxml_dir . -light_res none",
                                     namelst=>"fire_method='nofire',stream_fldfilename_lightng='build-namelist_test.pl'",
                                     phys=>"clm5_0",
                                   },
     "lightrescontradictlightfil"=>{ options=>"-bgc bgc -envxml_dir . -light_res 360x720",
                                     namelst=>"stream_fldfilename_lightng='build-namelist_test.pl'",
                                     phys=>"clm5_0",
                                   },
     "finundated and not methane"=>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"use_lch4=.false.,finundation_method='h2osfc'",
                                     phys=>"clm5_0",
                                   },
     "use_cn=true bgc=sp"        =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"use_cn=.true.",
                                     phys=>"clm4_5",
                                   },
     "freeliv wo fun"            =>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"freelivfix_intercept=9.",
                                     phys=>"clm4_5",
                                   },
     "use_cn=false bgc=bgc"      =>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"use_cn=.false.",
                                     phys=>"clm4_5",
                                   },
     "lower=aqu-45 with/o Zeng"  =>{ options=>"-envxml_dir .",
                                     namelst=>"lower_boundary_condition=4,soilwater_movement_method=1,use_bedrock=.false.",
                                     phys=>"clm5_0",
                                   },
     "Zeng w lower=flux"         =>{ options=>"-envxml_dir .",
                                     namelst=>"lower_boundary_condition=1,soilwater_movement_method=0,use_bedrock=.false.",
                                     phys=>"clm4_5",
                                   },
     "Zeng w lower=zeroflux"     =>{ options=>"-envxml_dir .",
                                     namelst=>"lower_boundary_condition=2,soilwater_movement_method=0",
                                     phys=>"clm4_5",
                                   },
     "Zeng w lower=table"        =>{ options=>"-envxml_dir .",
                                     namelst=>"lower_boundary_condition=3,soilwater_movement_method=0,use_bedrock=.false.",
                                     phys=>"clm4_5",
                                   },
     "use_vic=F with -vic op"    =>{ options=>"-vichydro -envxml_dir .",
                                     namelst=>"use_vichydro=.false.",
                                     phys=>"clm4_5",
                                   },
     "-vic with l_bnd=flux"      =>{ options=>"-vichydro -envxml_dir .",
                                     namelst=>"lower_boundary_condition=1",
                                     phys=>"clm4_5",
                                   },
     "-vic with l_bnd=zeroflux"  =>{ options=>"-vichydro -envxml_dir .",
                                     namelst=>"lower_boundary_condition=2",
                                     phys=>"clm4_5",
                                   },
     "bedrock with l_bnc=flux"   =>{ options=>"-envxml_dir .",
                                     namelst=>"use_bedrock=.true., lower_boundary_condition=1",
                                     phys=>"clm5_0",
                                   },
     "bedrock with l_bnc=tabl"   =>{ options=>"-envxml_dir .",
                                     namelst=>"use_bedrock=.true., lower_boundary_condition=3",
                                     phys=>"clm5_0",
                                   },
     "bedrock with l_bnc=aqui"   =>{ options=>"-envxml_dir .",
                                     namelst=>"use_bedrock=.true., lower_boundary_condition=4",
                                     phys=>"clm5_0",
                                   },
     "zengdeck with l_bnc=flux"  =>{ options=>"-envxml_dir .",
                                     namelst=>"soilwater_movement_method=0, lower_boundary_condition=1",
                                     phys=>"clm4_5",
                                   },
     "zengdeck with l_bnc=z-flux"=>{ options=>"-envxml_dir .",
                                     namelst=>"soilwater_movement_method=0, lower_boundary_condition=2",
                                     phys=>"clm4_5",
                                   },
     "zengdeck with l_bnc=tabl"  =>{ options=>"-envxml_dir .",
                                     namelst=>"soilwater_movement_method=0, lower_boundary_condition=3",
                                     phys=>"clm4_5",
                                   },
     "l_bnd=tabl with h2osfcfl=0"=>{ options=>"-envxml_dir .",
                                     namelst=>"h2osfcflag=0, lower_boundary_condition=3",
                                     phys=>"clm4_5",
                                   },
     "l_bnd=flux with h2osfcfl=0"=>{ options=>"-envxml_dir .",
                                     namelst=>"h2osfcflag=0, lower_boundary_condition=1",
                                     phys=>"clm4_5",
                                   },
     "l_bnd=zflux with h2osfcfl=0"=>{ options=>"-envxml_dir .",
                                     namelst=>"h2osfcflag=0, lower_boundary_condition=2",
                                     phys=>"clm4_5",
                                   },
     "h2osfcfl=0 with clm5.0"    =>{ options=>"-envxml_dir .",
                                     namelst=>"h2osfcflag=0",
                                     phys=>"clm5_0",
                                   },
     "45bad lnd_tuning_mode value" =>{ options=>"-lnd_tuning_mode clm5_0_GSWP3v1  -envxml_dir .",
                                     namelst=>"",
                                     phys=>"clm4_5",
                                   },
     "50bad lnd_tuning_mode value" =>{ options=>"-lnd_tuning_mode clm4_5_CRUNCEP  -envxml_dir .",
                                     namelst=>"",
                                     phys=>"clm5_0",
                                   },
     "bgc_spinup without cn"     =>{ options=>"-clm_accelerated_spinup on -bgc sp -envxml_dir .",
                                     namelst=>"spinup_state=1",
                                     phys=>"clm4_5",
                                   },
     "spinup=1 without bldnml op"=>{ options=>"-clm_accelerated_spinup off -bgc bgc -envxml_dir .",
                                     namelst=>"spinup_state=1",,
                                     phys=>"clm5_0",
                                   },
     "bgc_spinup without cn"     =>{ options=>"-clm_accelerated_spinup on -bgc sp -envxml_dir .",
                                     namelst=>"spinup_state=1",
                                     phys=>"clm4_5",
                                   },
     "baseflow w aquifer"        =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"baseflow_scalar=1.0, lower_boundary_condition=4,use_bedrock=.false.",
                                     phys=>"clm5_0",
                                   },
     "baseflow w table"          =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"baseflow_scalar=1.0, lower_boundary_condition=3,use_bedrock=.false.",
                                     phys=>"clm5_0",
                                   },
     "br_root and bgc=sp"        =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"br_root=1.0",
                                     phys=>"clm5_0",
                                   },
     "both co2_type and on nml"  =>{ options=>"-co2_type constant -envxml_dir .",
                                     namelst=>"co2_type='prognostic'",
                                     phys=>"clm5_0",
                                   },
     "lnd_frac set but nuopc"    =>{ options=>"-driver nuopc -lnd_frac $DOMFILE -envxml_dir .",
                                     namelst=>"",
                                     phys=>"clm6_0",
                                   },
     "driver is invalid"         =>{ options=>"-driver invalid_name -envxml_dir .",
                                     namelst=>"",
                                     phys=>"clm6_0",
                                   },
     "lnd_frac not set but lilac"=>{ options=>"-driver nuopc -lilac -envxml_dir . -lnd_frac UNSET",
                                     namelst=>"fsurdat='surfdata.nc'",
                                     phys=>"clm6_0",
                                   },
     "fatmlndfrc set but nuopc"  =>{ options=>"-driver nuopc -envxml_dir .",
                                     namelst=>"fatmlndfrc='frac.nc'",
                                     phys=>"clm6_0",
                                   },
     "branch but NO nrevsn"      =>{ options=>"-clm_start_type branch -envxml_dir .",
                                     namelst=>"",
                                     phys=>"clm5_0",
                                   },
     "glc_nec inconsistent"      =>{ options=>"-envxml_dir .",
                                     namelst=>"maxpatch_glc=5",
                                     phys=>"clm5_0",
                                   },
     "NoGLCMec"                  =>{ options=>"-envxml_dir . -glc_nec 0",
                                     namelst=>"",
                                     phys=>"clm4_5",
                                   },
     "UpdateGlcContradict"       =>{ options=>"-envxml_dir .",
                                     namelst=>"glc_do_dynglacier=.false.",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm4_5",
                                   },
     "matrixWOBGC"               =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"use_matrixcn=.true.",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm5_0",
                                   },
     "soilmatrixWOBGC"           =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"use_soil_matrixcn=T",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm5_0",
                                   },
     "soilmatrixWmimics"         =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_soil_matrixcn=T,soil_decomp_method='MIMICSWieder2015'",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm5_0",
                                   },
     "matrixcn_diagWOmatrix"     =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_soil_matrixcn=.false.,use_matrixcn=F,hist_wrt_matrixcn_diag=T",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm5_0",
                                   },
     "spinupWOsoilmatrix"        =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_soil_matrixcn=F,use_matrixcn=T,spinup_matrixcn=T",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm5_0",
                                   },
     "sasuspinupWOsoilmatx"      =>{ options=>"-envxml_dir . -bgc bgc -clm_accelerated_spinup sasu",
                                     namelst=>"use_soil_matrixcn=.false.,use_matrixcn=.false.",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm6_0",
                                   },
     "sasuspinupWOCN"            =>{ options=>"-envxml_dir . -bgc sp  -clm_accelerated_spinup sasu",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm6_0",
                                   },
     "nyrforceWOspinup"          =>{ options=>"-envxml_dir . -bgc bgc -clm_accelerated_spinup sasu",
                                     namelst=>"use_matrixcn=.false.,spinup_matrixcn=F,nyr_forcing=20",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm5_0",
                                   },
     "nyrsasuGTnyrforce"         =>{ options=>"-envxml_dir . -bgc bgc -clm_accelerated_spinup sasu",
                                     namelst=>"use_matrixcn=.false.,spinup_matrixcn=T,nyr_forcing=20,nyr_sasu=21",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm5_0",
                                   },
     "iloopZero"                 =>{ options=>"-envxml_dir . -bgc bgc -clm_accelerated_spinup sasu",
                                     namelst=>"use_matrixcn=.false.,spinup_matrixcn=T,iloop_avg=0",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm5_0",
                                   },
     "matrixspinupWADmode"        =>{ options=>"-envxml_dir . -bgc bgc -clm_accelerated_spinup sasu",
                                     namelst=>"spinup_matrixcn=T,spinup_state=2",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm5_0",
                                   },
     "matrixspinupWclmaccell"     =>{ options=>"-envxml_dir . -bgc bgc -clm_accelerated_spinup off",
                                     namelst=>"use_soil_matrixcn=T,spinup_matrixcn=T",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm5_0",
                                   },
     "fatesWuse_cnmatrix"        =>{ options=>"-envxml_dir . -bgc fates",
                                     namelst=>"use_matrixcn=.true.",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm5_0",
                                   },
     "fatesWuse_soilcnmatrix"    =>{ options=>"-envxml_dir . -bgc fates",
                                     namelst=>"use_soil_matrixcn=.true.",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm5_0",
                                   },
     "useFATESContradict"        =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates=.false.",
                                     phys=>"clm4_5",
                                   },
     "useFATESContradict2"       =>{ options=>"-envxml_dir . -no-megan",
                                     namelst=>"use_fates=.true.",
                                     phys=>"clm4_5",
                                   },
     "useFATESWCN"               =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_cn=.true.",
                                     phys=>"clm5_0",
                                   },
     "useFATESWcrop"             =>{ options=>"-bgc fates -envxml_dir . -no-megan -crop",
                                     namelst=>"",
                                     phys=>"clm6_0",
                                   },
     "useFATESWcreatecrop"       =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"create_crop_landunit=.true.",
                                     phys=>"clm5_0",
                                   },
     "useFATESWn_dom_pft"        =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"n_dom_pfts = 1",
                                     phys=>"clm5_0",
                                   },
     "useFATESWbMH"              =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_biomass_heat_storage=.true.",
                                     phys=>"clm6_0",
                                   },
     "FireNoneButFATESfireon"    =>{ options=>"-bgc fates -envxml_dir . -no-megan -light_res none",
                                     namelst=>"fates_spitfire_mode=4",
                                     phys=>"clm6_0",
                                   },
     "FATESwspitfireOffLigtOn"    =>{ options=>"-bgc fates -envxml_dir . -no-megan -light_res 360x720",
                                     namelst=>"fates_spitfire_mode=0",
                                     phys=>"clm6_0",
                                   },
     "useFATESWluna"             =>{ options=>"--bgc fates --envxml_dir . --no-megan",
                                     namelst=>"use_luna=TRUE",
                                     phys=>"clm6_0",
                                   },
     "useFATESWfun"              =>{ options=>"--bgc fates --envxml_dir . --no-megan",
                                     namelst=>"use_fun=TRUE",
                                     phys=>"clm6_0",
                                   },
     "useFATESWOsuplnitro"       =>{ options=>"--bgc fates --envxml_dir . --no-megan",
                                     namelst=>"suplnitro='NONE'",
                                     phys=>"clm6_0",
                                   },
     "FireNoneButBGCfireon"    =>{ options=>"-bgc bgc -envxml_dir . -light_res none",
                                     namelst=>"fire_method='li2021gswpfrc'",
                                     phys=>"clm6_0",
                                   },
     "createcropFalse"           =>{ options=>"-bgc bgc -envxml_dir . -no-megan",
                                     namelst=>"create_crop_landunit=.false.",
                                     phys=>"clm5_0",
                                   },
     "usespitfireButNOTFATES"    =>{ options=>"-envxml_dir . -no-megan",
                                     namelst=>"fates_spitfire_mode=1",
                                     phys=>"clm4_5",
                                   },
     "usespitfireusefatessp"    =>{ options=>"-envxml_dir . --bgc fates",
                                     namelst=>"fates_spitfire_mode=1,use_fates_sp=.true.",
                                     phys=>"clm5_0",
                                   },
     "usefatesspusefateshydro"   =>{ options=>"-envxml_dir . --bgc fates",
                                     namelst=>"use_fates_sp=.true.,use_fates_planthydro=.true.",
                                     phys=>"clm5_0",
                                   },
     "useloggingButNOTFATES"     =>{ options=>"-envxml_dir . -no-megan",
                                     namelst=>"fates_harvest_mode='event_code'",
                                     phys=>"clm4_5",
                                   },
     "useinventorybutnotfile"    =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_inventory_init=.true.",
                                     phys=>"clm4_5",
                                   },
     "inventoryfileDNE"          =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_inventory_init=.true., fates_inventory_ctrl_filename='zztop'",
                                     phys=>"clm4_5",
                                   },
     "useFATESLUH2butnotfile"    =>{ options=>"--res 0.9x1.25 --bgc fates --envxml_dir . --no-megan",
                                     namelst=>"use_fates_luh=.true.",
                                     phys=>"clm4_5",
                                   },
     "useFATESLUPFTbutnotfile"   =>{ options=>"--res 0.9x1.25 --bgc fates --envxml_dir . --no-megan",
                                     namelst=>"use_fates_lupft=.true.",
                                     phys=>"clm4_5",
                                   },
     "inventoryfileDNE"          =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_luh=.true., fluh_timeseries='zztop'",
                                     phys=>"clm4_5",
                                   },
     "useMEGANwithFATES"         =>{ options=>"-bgc fates -envxml_dir . -megan",
                                     namelst=>"",
                                     phys=>"clm4_5",
                                   },
     "useFIREEMISwithFATES"      =>{ options=>"-bgc fates -envxml_dir . -fire_emis --no-megan",
                                    namelst=>"",
                                    phys=>"clm4_5",
                                 },
     "useDRYDEPwithFATES"        =>{ options=>"--bgc fates --envxml_dir . --no-megan --drydep",
                                     namelst=>"",
                                     phys=>"clm4_5",
                                   },
     "useFATESSPWONOCOMP"        =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_sp=T,use_fates_nocomp=F",
                                     phys=>"clm5_0",
                                   },
     "useFATESSPwithLUH"        =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_sp=T,use_fates_luh=T",
                                     phys=>"clm5_0",
                                   },
     "useFATESPOTVEGwithHARVEST" =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_potentialveg=T,fates_harvest_mode='event_code',use_fates_luh=T",
                                     phys=>"clm5_0",
                                   },
     "useFATESHARVEST3WOLUH"     =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_luh=F,fates_harvest_mode='luhdata_area'",
                                     phys=>"clm5_0",
                                   },
     "useFATESLUPFTWOLUH"        =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_lupft=T,use_fates_luh=F",
                                     phys=>"clm5_0",
                                   },
     "useFATESLUPFTWONOCOMP"     =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_lupft=T,use_fates_nocomp=F",
                                     phys=>"clm5_0",
                                   },
     "useFATESLUPFTWOFBG"        =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_lupft=T,use_fates_fixedbiogeog=F",
                                     phys=>"clm5_0",
                                   },
     "useFATESTRANSWdynPFT"      =>{ options=>"-bgc fates -envxml_dir . -use_case 20thC_transient -no-megan",
                                     namelst=>"do_transient_pfts=T",
                                     phys=>"clm5_0",
                                   },
     "useHYDSTwithFATES"         =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_hydrstress=.true.",
                                     phys=>"clm5_0",
                                   },
     "useMeierwithFATES"         =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"z0param_method=Meier2022",
                                     phys=>"clm5_0",
                                   },
     "noanthro_w_crop"            =>{ options=>"-envxml_dir . -res 0.9x1.25 -bgc bgc -crop -use_case 1850_noanthro_control",
                                     namelst=>"",
                                     phys=>"clm5_0",
                                   },
     "noanthro_w_irrig"           =>{ options=>"-envxml_dir . -res 0.9x1.25 -bgc bgc -use_case 1850_noanthro_control",
                                     namelst=>"irrigate=T",
                                     phys=>"clm5_0",
                                   },
     "spdotransconflict"          =>{ options=>"-envxml_dir . -bgc sp -use_case 20thC_transient",
                                     namelst=>"do_transient_pfts=T,do_transient_crops=.false.",
                                     phys=>"clm5_0",
                                   },
     "dogrossandsp"               =>{ options=>"--envxml_dir . --bgc sp --use_case 20thC_transient",
                                     namelst=>"do_grossunrep=.true.",
                                     phys=>"clm5_0",
                                   },
     "dogrossandfates"            =>{ options=>"--envxml_dir . --bgc fates --use_case 20thC_transient --no-megan",
                                     namelst=>"do_grossunrep=.true.",
                                     phys=>"clm5_0",
                                   },
     "dogrossandnottrans"         =>{ options=>"--envxml_dir . --bgc bgc --use_case 2000_control",
                                     namelst=>"do_grossunrep=.true.",
                                     phys=>"clm5_0",
                                   },
     "nocropwfert"                =>{ options=>"-envxml_dir . -bgc sp -no-crop",
                                     namelst=>"use_fertilizer=T",
                                     phys=>"clm5_0",
                                   },
     "lmr1WOcn"                   =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"leafresp_method=1",
                                     phys=>"clm5_0",
                                   },
     "lmr2WOcn"                   =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"leafresp_method=2",
                                     phys=>"clm5_0",
                                   },
     "lmr0Wcn"                    =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"leafresp_method=0",
                                     phys=>"clm5_0",
                                   },
     "nofireButSetcli_scale"     =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"fire_method='nofire', cli_scale=5.",
                                     phys=>"clm5_0",
                                   },
     "nocnButSetrh_low"          =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"rh_low=5.",
                                     phys=>"clm5_0",
                                   },
     "funWOcn"                   =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"use_fun=.true.",
                                     phys=>"clm5_0",
                                   },
     "flexCNWOcn"                =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"use_flexibleCN=.true.",
                                     phys=>"clm5_0",
                                   },
     "flexCNFUNwcarbonresp"      =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_flexibleCN=.true.,use_FUN=.true.,carbon_resp_opt=1",
                                     phys=>"clm5_0",
                                   },
     "funWOnitrif"               =>{ options=>"-envxml_dir .",
                                     namelst=>"use_fun=.true., use_nitrif_denitrif=.false.",
                                     phys=>"clm5_0",
                                   },
     "SPModeWNitrifNMethane"     =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"use_lch4=.true., use_nitrif_denitrif=.true.",
                                     phys=>"clm5_0",
                                   },
     "knitrmaxWOnitrif"          =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_nitrif_denitrif=.false., k_nitr_max=1.0",
                                     phys=>"clm5_0",
                                   },
     "respcoefWOnitrif"          =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_nitrif_denitrif=.false., denitrif_respiration_coefficient=1.0",
                                     phys=>"clm5_0",
                                   },
     "respexpWOnitrif"           =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_nitrif_denitrif=.false., denitrif_respiration_exponent=1.0",
                                     phys=>"clm5_0",
                                   },
     "lunaWSPandlnctrue"         =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"use_luna=.true., lnc_opt=.true.",
                                     phys=>"clm5_0",
                                   },
     "envxml_not_dir"            =>{ options=>"-envxml_dir myuser_nl_clm",
                                     namelst=>"",
                                     phys=>"clm5_0",
                                   },
     "envxml_emptydir"           =>{ options=>"-envxml_dir xFail",
                                     namelst=>"",
                                     phys=>"clm5_0",
                                   },
     "fates_non_sp_laistreams"   =>{ options=>"--envxml_dir . --bgc fates",
                                     namelst=>"use_lai_streams=.true., use_fates_sp=.false.",
                                     phys=>"clm5_0",
                                     },
     "bgc_non_sp_laistreams"     =>{ options=>"--envxml_dir . -bgc bgc",
                                     namelst=>"use_lai_streams=.true.",
                                     phys=>"clm5_0",
                                     },
     "bgc_laistreams_input"     =>{ options=>"--envxml_dir . --bgc bgc",
                                     namelst=>"stream_year_first_lai=1999",
                                     phys=>"clm5_0",
                                     },
     "crop_laistreams_input"     =>{ options=>"--envxml_dir . --bgc sp --crop",
                                     namelst=>"use_lai_streams=.true.",
                                     phys=>"clm5_0",
                                     },
     "soil_erod_wo_Zender"      =>{ options=>"--envxml_dir . --ignore_warnings",
                                     namelst=>"dust_emis_method='Leung_2023', stream_meshfile_zendersoilerod = '/dev/null'",
                                     phys=>"clm6_0",
                                     },
     "soil_erod_wo_lnd_source"  =>{ options=>"--envxml_dir .",
                                     namelst=>"dust_emis_method='Zender_2003', stream_fldfilename_zendersoilerod = '/dev/null', zender_soil_erod_source='atm'",
                                     phys=>"clm6_0",
                                     },
     "soil_erod_none_w_Zender"  =>{ options=>"--envxml_dir .",
                                     namelst=>"dust_emis_method='Zender_2003', zender_soil_erod_source='none'",
                                     phys=>"clm6_0",
                                     },
     "soil_erod_bad_w_Zender"   =>{ options=>"--envxml_dir .",
                                     namelst=>"dust_emis_method='Zender_2003', zender_soil_erod_source='zztop'",
                                     phys=>"clm6_0",
                                     },
     "Set_Dust_When_CAM_Sets"   =>{ options=>"--envxml_dir .",
                                     namelst=>"dust_emis_method='Zender_2003'",
                                     LND_SETS_DUST_EMIS_DRV_FLDS=>"FALSE",
                                     phys=>"clm6_0",
                                     },
               );
foreach my $key ( keys(%failtest) ) {
   print( "$key\n" );
   my $var;
   foreach $var ( "phys" , "options", "namelst" ) {
      if ( not exists $failtest{$key}{$var} ) {
         die  "ERROR: Subkey $var does not exist for failtest $key\nERROR:Check if you spelled $var correctly\n"
      }
   }

   &make_config_cache($failtest{$key}{"phys"});
   my $options  = $failtest{$key}{"options"};
   my $namelist = $failtest{$key}{"namelst"};
   my %settings;
   foreach my $xmlvar ( "GLC_TWO_WAY_COUPLING", "LND_SETS_DUST_EMIS_DRV_FLDS") {
      if ( defined($failtest{$key}{$xmlvar}) ) {
         $settings{$xmlvar} = $failtest{$key}{$xmlvar};
      }
   }
   &make_env_run( %settings );
   eval{ system( "$bldnml $options -namelist \"&clmexp $namelist /\" > $tempfile 2>&1 " ); };
   isnt( $?, 0, $key );
   system( "cat $tempfile" );
}


print "\n===============================================================================\n";
print "Start Warning testing.  These should fail unless -ignore_warnings option is used \n";
print "=================================================================================\n";

# Warning testing, do things that give warnings, unless -ignore_warnings option is used

my %warntest = (
     # Warnings without the -ignore_warnings option given
     "bgcspin_w_suplnitro"       =>{ options=>"-envxml_dir . -bgc bgc -clm_accelerated_spinup on",
                                     namelst=>"suplnitro='ALL'",
                                     phys=>"clm5_0",
                                   },
     "bgc=bgc WO nitrif_denit"   =>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"use_nitrif_denitrif=.false.",
                                     phys=>"clm4_5",
                                   },
     "methane off W nitrif_denit"=>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"use_nitrif_denitrif=.true.,use_lch4=.false.",
                                     phys=>"clm6_0",
                                   },
     "soilm_stream w transient"  =>{ options=>"-res 0.9x1.25 -envxml_dir . -use_case 20thC_transient",
                                     namelst=>"use_soil_moisture_streams=T,soilm_tintalgo='linear'",
                                     phys=>"clm5_0",
                                   },
     "missing_ndep_file"         =>{ options=>"-envxml_dir . -bgc bgc -ssp_rcp SSP5-3.4",
                                     namelst=>"",
                                     phys=>"clm5_0",
                                   },
     "bad_megan_spec"            =>{ options=>"-envxml_dir . -bgc bgc -megan",
                                     namelst=>"megan_specifier='ZZTOP=zztop%'",
                                     phys=>"clm4_5",
                                   },
     "FUN_wo_flexCN"             =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_fun=.true.,use_flexiblecn=.false.",
                                     phys=>"clm6_0",
                                   },
     "Set coldtemp wo coldstart" =>{ options=>"-envxml_dir . --clm_start_type startup",
                                     namelst=>"use_excess_ice=.true.,excess_ice_coldstart_temp=-10.",
                                     phys=>"clm6_0",
                                   },
     "Set colddepth wo coldstart" =>{ options=>"-envxml_dir . --clm_start_type startup",
                                     namelst=>"use_excess_ice=.true.,excess_ice_coldstart_depth=0.5",
                                     phys=>"clm6_0",
                                   },
     "PrigentOnWOLeung"          =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"use_prigent_roughness=.true.,dust_emis_method='Zender_2003'",
                                     phys=>"clm6_0",
                                   },
     "NotNEONbutNEONlightres"    =>{ options=>"--res CLM_USRDAT --clm_usr_name regional --envxml_dir . --bgc bgc --light_res 106x174",
                                     namelst=>"fsurdat='build-namelist_test.pl'",
                                     phys=>"clm6_0",
                                   },
     "hillslope with init_interp"=>{ options=>"--res 10x15 --bgc bgc --envxml_dir .",
                                     namelst=>"use_init_interp=.true.,use_hillslope=.true.,hillslope_file='/dev/null'",
                                     phys=>"clm6_0",
                                   },
               );
foreach my $key ( keys(%warntest) ) {
   print( "$key\n" );

   my $var;
   foreach $var ( "phys" , "options", "namelst" ) {
      if ( not exists $warntest{$key}{$var} ) {
         die  "ERROR: Subkey $var does not exist for warntest $key\nERROR:Check if you spelled $var correctly\n"
      }
   }

   &make_config_cache($warntest{$key}{"phys"});
   my $options  = $warntest{$key}{"options"};
   my $namelist = $warntest{$key}{"namelst"};
   my %settings;
   foreach my $xmlvar ( "GLC_TWO_WAY_COUPLING" ) {
      if ( defined($failtest{$key}{$xmlvar}) ) {
         $settings{$xmlvar} = $failtest{$key}{$xmlvar};
      }
   }
   &make_env_run( %settings );
   eval{ system( "$bldnml $options -namelist \"&clmexp $namelist /\" > $tempfile 2>&1 " ); };
   isnt( $?, 0, $key );
   system( "cat $tempfile" );
   # Now run with -ignore_warnings and make sure it works
   $options .= " -ignore_warnings";
   eval{ system( "$bldnml $options -namelist \"&clmexp $namelist /\" > $tempfile 2>&1 " ); };
   is( $?, 0, $key );
   is( $@, '', "$options" );
   system( "cat $tempfile" );
}

print "\n===============================================================================\n";
print "Ensure cold starts with finidat are handled properly \n";
print "=================================================================================\n";

my %coldwfinidat = (
     "bgc"   => { options=>"-envxml_dir . -clm_start_type cold",
                  namelst=>"finidat = 'testfile.nc'",
                  phys=>"clm5_0",
                  expected_fail=>1,
                },
     "fates" => { options=>"-envxml_dir . -clm_start_type cold -bgc fates -no-megan",
                  namelst=>"finidat = 'testfile.nc', use_fates = .true.",
                  phys=>"clm5_0",
                  expected_fail=>0,
                },
);
my $finidat;
foreach my $key ( keys(%coldwfinidat) ) {
   print( "$key\n" );

   my $var;
   foreach $var ( "phys" , "options", "namelst", "expected_fail" ) {
      if ( not exists $coldwfinidat{$key}{$var} ) {
         die  "ERROR: Subkey $var does not exist for coldwfinidat $key\nERROR:Check if you spelled $var correctly\n"
      }
   }

   &make_config_cache($coldwfinidat{$key}{"phys"});
   my $options  = $coldwfinidat{$key}{"options"};
   my $namelist = $coldwfinidat{$key}{"namelst"};
   my $expected_fail = $coldwfinidat{$key}{"expected_fail"};
   my %settings;
   &make_env_run( %settings );

   # Should fail if expected to, pass otherwise
   eval{ system( "$bldnml $options -namelist \"&clmexp $namelist /\" > $tempfile 2>&1 " ); };
   is( $? eq 0, $expected_fail eq 0, "coldwfinidat $key run");

   if ( $expected_fail ) {
      # Now run with -ignore_warnings and make sure it still doesn't work
      $options .= " -ignore_warnings";
      eval{ system( "$bldnml $options -namelist \"&clmexp $namelist /\" > $tempfile 2>&1 " ); };
      isnt( $?, 0, "coldwfinidat $key run -ignore_warnings" );
   } else {
      # Check that finidat was correctly set
      $finidat = `grep finidat lnd_in`;
      ok ( $finidat =~ "testfile.nc", "coldwfinidat $key finidat? $finidat" );
   }
}

#
# Loop over all physics versions
#
foreach my $phys ( "clm4_5", "clm5_0", "clm6_0" ) {
$mode = "-phys $phys";
&make_config_cache($phys);

print "\n========================================================================\n";
print "Test ALL resolutions that have surface datasets with SP for 1850 and 2000\n";
print "========================================================================\n";

# Check for ALL resolutions with CLM50SP
my @resolutions = ( "360x720cru", "10x15", "4x5", "0.9x1.25", "1.9x2.5", "ne3np4.pg3", "ne16np4.pg3", "ne30np4", "ne30np4.pg2", "ne30np4.pg3", "ne120np4.pg3", "ne0np4CONUS.ne30x8", "ne0np4.ARCTIC.ne30x4", "ne0np4.ARCTICGRIS.ne30x8", "C96", "mpasa480", "mpasa120" );
my @only2000_resolutions = ( "1x1_numaIA", "1x1_brazil", "1x1_mexicocityMEX", "1x1_vancouverCAN", "1x1_urbanc_alpha", "5x5_amazon", "0.125nldas2", "mpasa60", "mpasa15", "mpasa3p75" );
my @regional;
foreach my $res ( @resolutions ) {
   chomp($res);
   print "=== Test $res === \n";
   foreach my $use_case ( "1850_control", "2000_control" ) {
      # Skip resolutions that only have 2000 versions
      if ( ($use_case eq "1850_control") && ($res ~~ @only2000_resolutions) ) {
         next;
      }
      print "=== Test $use_case === \n";
      my $options  = "-res $res -bgc sp -envxml_dir . --use_case $use_case";

      &make_env_run();
      eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
      is( $@, '', "$options" );

      $cfiles->checkfilesexist( "$options", $mode );

      $cfiles->shownmldiff( "default", "standard" );
      if ( defined($opts{'compare'}) ) {
         $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
         $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
      }

      if ( defined($opts{'generate'}) ) {
         $cfiles->copyfiles( "$options", $mode );
      }
      &cleanup(); print "\n";
   }
}

print "\n==================================================\n";
print " Test important resolutions for BGC and historical\n";
print "==================================================\n";

my @resolutions = ( "4x5", "10x15", "360x720cru", "ne30np4.pg3", "ne3np4.pg3", "1.9x2.5", "0.9x1.25", "C96", "mpasa120" );
my @regional;
my $nlbgcmode = "bgc";
my $mode = "$phys-$nlbgcmode";
foreach my $res ( @resolutions ) {
   chomp($res);
   print "=== Test $res === \n";
   my $options  = "-res $res -envxml_dir . -bgc $nlbgcmode --use_case 20thC_transient";

   &make_env_run();
   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );

   $cfiles->checkfilesexist( "$options", $mode );

   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }

   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup(); print "\n";
}

print "\n==================================================\n";
print " Test all use-cases over all physics options\n";
print "==================================================\n";

# Run over all use-cases for f09 and all physics...
my $list = `$bldnml -use_case list 2>&1 | grep "use case"`;
my @usecases;
if ( $list =~ /build-namelist : use cases : (.+)$/ ) {
  @usecases  = split( / /, $1 );
} else {
  die "ERROR:: Trouble getting list of use-cases\n";
}
if ( $#usecases != 15 ) {
  print "use-cases = @usecases\n";
  die "ERROR:: Number of use-cases isn't what's expected\n";
}
my @expect_fails = ( "1850-2100_SSP5-3.4_transient", "1850-2100_SSP4-3.4_transient", "2018-PD_transient", "1850-2100_SSP1-1.9_transient",
                      "1850-2100_SSP4-6.0_transient", "2018_control" );
foreach my $phys ( "clm4_5", "clm5_0", "clm6_0" ) {
   print "physics = $phys\n";
   &make_config_cache($phys);
   foreach my $usecase ( @usecases ) {
      print "usecase = $usecase\n";
      $options = "-res 0.9x1.25 -use_case $usecase  -envxml_dir .";
      &make_env_run();
      my $expect_fail = undef;
      foreach my $failusecase ( @expect_fails ) {
         if ( $failusecase eq $usecase ) {
            $expect_fail = 1;
            last;
         }
      }
      eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
      if ( ! defined($expect_fail) ) {
         is( $@, '', "options: $options" );
         $cfiles->checkfilesexist( "$options", $mode );
         $cfiles->shownmldiff( "default", "standard" );
         if ( defined($opts{'compare'}) ) {
            $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
            $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
         }
         if ( defined($opts{'generate'}) ) {
            $cfiles->copyfiles( "$options", $mode );
         }
      } else {
         isnt( $@, 0, "options: $options" );
      }
      &cleanup();
   }
}

print "\n=======================================================================================\n";
print " Test the seperate initial condition files, for ones not tested elsewhere\n";
print "=========================================================================================\n";

my %finidat_files = (
     "f091850Clm45BgcGSW"        =>{ phys =>"clm4_5",
                                     atm_forc=>"GSWP3v1",
                                     res => "0.9x1.25",
                                     bgc => "bgc",
                                     crop => "--no-crop",
                                     use_case => "1850_control",
                                     start_ymd => "18500101",
                                     namelist => "irrigate=T",
                                   },
     "f091850Clm45BgcCRU"        =>{ phys =>"clm4_5",
                                     atm_forc=>"CRUv7",
                                     res => "0.9x1.25",
                                     bgc => "bgc",
                                     crop => "--no-crop",
                                     use_case => "1850_control",
                                     start_ymd => "18500101",
                                     namelist => "irrigate=T",
                                   },
     "f091850Clm45BgcCAM6"       =>{ phys =>"clm4_5",
                                     atm_forc=>"cam6.0",
                                     res => "0.9x1.25",
                                     bgc => "bgc",
                                     crop => "--crop",
                                     use_case => "1850_control",
                                     start_ymd => "18500101",
                                     namelist => "irrigate=F",
                                   },
     "f091850Clm50BgcGSW"        =>{ phys =>"clm5_0",
                                     atm_forc=>"GSWP3v1",
                                     res => "0.9x1.25",
                                     bgc => "bgc",
                                     crop => "--crop",
                                     use_case => "1850_control",
                                     start_ymd => "18500101",
                                     namelist => "irrigate=F",
                                   },
     "f091850Clm50SpGSW"         =>{ phys =>"clm5_0",
                                     atm_forc=>"GSWP3v1",
                                     res => "0.9x1.25",
                                     bgc => "sp",
                                     crop => "--no-crop",
                                     use_case => "1850_control",
                                     start_ymd => "18500101",
                                     namelist => "irrigate=T",
                                   },
     "f091850Clm50BgcCRU"        =>{ phys =>"clm5_0",
                                     atm_forc=>"CRUv7",
                                     res => "0.9x1.25",
                                     bgc => "bgc",
                                     crop => "--crop",
                                     use_case => "1850_control",
                                     start_ymd => "18500101",
                                     namelist => "irrigate=F",
                                   },
     "f091850Clm50SpCRU"         =>{ phys =>"clm5_0",
                                     atm_forc=>"CRUv7",
                                     res => "0.9x1.25",
                                     bgc => "sp",
                                     crop => "--no-crop",
                                     use_case => "1850_control",
                                     start_ymd => "18500101",
                                     namelist => "irrigate=T",
                                   },
     "f091850Clm50BgcCAM6"       =>{ phys =>"clm5_0",
                                     atm_forc=>"cam6.0",
                                     res => "0.9x1.25",
                                     bgc => "bgc",
                                     crop => "--crop",
                                     use_case => "1850_control",
                                     start_ymd => "18500101",
                                     namelist => "irrigate=F",
                                   },
   );

foreach my $key ( keys(%finidat_files) ) {
   print( "$key\n" );

   my $var;
   foreach $var ( "phys" , "atm_forc", "res", "bgc", "crop", "use_case", "start_ymd", "namelist" ) {
      if ( not exists $finidat_files{$key}{$var} ) {
         die  "ERROR: Subkey $var does not exist for finidat_file $key\nERROR:Check if you spelled $var correctly\n"
      }
   }

   my $phys = $finidat_files{$key}{'phys'};
   print "physics = $phys\n";
   &make_config_cache($phys);
   my $usecase = $finidat_files{$key}{'use_case'};
   my $bgc = $finidat_files{$key}{'bgc'};
   my $res = $finidat_files{$key}{'res'};
   my $crop = $finidat_files{$key}{'crop'};
   my $namelist = $finidat_files{$key}{'namelist'};
   my $start_ymd = $finidat_files{$key}{'start_ymd'};
   my $lnd_tuning_mode = "${phys}_" . $finidat_files{$key}{'atm_forc'};
   $options = "-bgc $bgc -res $res -use_case $usecase -envxml_dir . $crop --lnd_tuning_mode $lnd_tuning_mode " . 
              "-namelist '&a start_ymd=$start_ymd, $namelist/'";
   &make_env_run();
   eval{ system( "$bldnml $options  > $tempfile 2>&1 " ); };
   is( $@, '', "options: $options" );
   my $finidat = `grep finidat lnd_in`;
   if ( $finidat =~ /initdata_map/ ) {
      my $result;
      eval( $result = `grep use_init_interp lnd_in` );
      is ( $result =~ /.true./, 1, "use_init_interp needs to be true here: $result");
   }
   $cfiles->checkfilesexist( "$options", $mode );
   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
         $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
         $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
         $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}

print "\n==================================================\n";
print "Test crop resolutions \n";
print "==================================================\n";

# Check for crop resolutions
my @crop1850_res = ( "1x1_smallvilleIA", "1x1_cidadinhoBR" );
foreach my $res ( @crop1850_res ) {
   my $use_case = "1850_control";
   if ( $res =~ /1x1_cidadinhoBR/ ) {
      $use_case = "2000_control";
   }
   $options = "-bgc bgc -crop -res $res -use_case $use_case -envxml_dir .";
   &make_env_run();
   eval{ system( "$bldnml $options  > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );
   $cfiles->checkfilesexist( "$options", $mode );
   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}

my @crop_res = ( "1x1_numaIA", "4x5", "10x15", "0.9x1.25", "1.9x2.5", "ne3np4.pg3", "ne30np4", "ne30np4.pg3", "C96", "mpasa120" );
foreach my $res ( @crop_res ) {
   $options = "-bgc bgc -crop -res $res -envxml_dir .";
   &make_env_run();
   eval{ system( "$bldnml $options  > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );
   $cfiles->checkfilesexist( "$options", $mode );
   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}
print "\n==================================================\n";
print " Test glc_mec resolutions \n";
print "==================================================\n";

# Check for glc_mec resolutions
#
# NOTE(wjs, 2017-12-17) I'm not sure if these glc_mec-specific tests are
# still needed: are they covered with other tests now that we always run
# with glc_mec? Some historical notes: (1) The three resolutions listed
# here used to be the only three with which you could run glc_mec; now
# you can run glc_mec with all resolutions. (2) This used to point to
# all of the glacierMEC use cases; now we don't have glacierMEC-specific
# use cases, but I've kept these pointing to the equivalent normal use
# cases; I'm not sure if it's actually important to test this with all
# of the different use cases.
my @glc_res = ( "0.9x1.25", "1.9x2.5" );
my @use_cases = (
                  "1850-2100_SSP2-4.5_transient",
                  "1850_control",
                  "2000_control",
                  "2010_control",
                  "20thC_transient",
                 );
foreach my $res ( @glc_res ) {
   foreach my $usecase ( @use_cases ) {
      my $startymd = undef;
      if ( ($usecase eq "1850_control") || ($usecase eq "20thC_transient") ) {
         $startymd = 18500101;
      } elsif ( $usecase eq "2000_control") {
         $startymd = 20000101;
      } elsif ( $usecase eq "2010_control") {
         $startymd = 20100101;
      } else {
         $startymd = 20150101;
      }
      $options = "-bgc bgc -res $res -use_case $usecase -envxml_dir . -namelist '&a start_ymd=$startymd/'";
      &make_env_run();
      eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
      is( $@, '', "$options" );
      $cfiles->checkfilesexist( "$options", $mode );
      $cfiles->shownmldiff( "default", "standard" );
      if ( defined($opts{'compare'}) ) {
         $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
         $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
      }
      if ( defined($opts{'generate'}) ) {
         $cfiles->copyfiles( "$options", $mode );
      }
      &cleanup();
   }
}
# Transient 20th Century simulations
my @tran_res = ( "0.9x1.25", "1.9x2.5", "ne30np4.pg3", "10x15" );
my $usecase  = "20thC_transient";
my $GLC_NEC         = 10;
foreach my $res ( @tran_res ) {
   $options = "-res $res -use_case $usecase -envxml_dir . -namelist '&a start_ymd=18500101/' -bgc bgc -crop -namelist '&a do_grossunrep=T/'";
   &make_env_run();
   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );
   $cfiles->checkfilesexist( "$options", $mode );
   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}
# Transient ssp_rcp scenarios that work
my @tran_res = ( "4x5", "0.9x1.25", "1.9x2.5", "10x15", "360x720cru", "ne3np4.pg3", "ne16np4.pg3", "ne30np4.pg3", "C96", "mpasa120" );
foreach my $usecase ( "1850-2100_SSP2-4.5_transient" ) {
   my $startymd = 20150101;
   foreach my $res ( @tran_res ) {
      $options = "-res $res -bgc bgc -crop -use_case $usecase -envxml_dir . -namelist '&a start_ymd=$startymd/'";
      &make_env_run();
      eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
      is( $@, '', "$options" );
      $cfiles->checkfilesexist( "$options", $mode );
      $cfiles->shownmldiff( "default", "standard" );
      if ( defined($opts{'compare'}) ) {
         $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
         $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
      }
      if ( defined($opts{'generate'}) ) {
         $cfiles->copyfiles( "$options", $mode );
      }
      &cleanup();
   }
}
}  # End loop over all physics versions
#
# End loop over versions
#

print "\n==================================================\n";
print "Test clm4.5/clm5.0/clm6_0 resolutions \n";
print "==================================================\n";

foreach my $phys ( "clm4_5", "clm5_0", "clm6_0" ) {
  my $mode = "-phys $phys";
  &make_config_cache($phys);
  my @clmoptions = ( "-bgc bgc -envxml_dir .", "-bgc bgc -envxml_dir . -clm_accelerated_spinup=on", "-bgc bgc -envxml_dir . -light_res 360x720",
                     "-bgc sp -envxml_dir . -vichydro", "-bgc bgc -dynamic_vegetation -ignore_warnings",
                     "-bgc bgc -clm_demand flanduse_timeseries -sim_year 1850-2000 -namelist '&a start_ymd=18500101/'",
                     "-bgc bgc -envxml_dir . -namelist '&a use_c13=.true.,use_c14=.true.,use_c14_bombspike=.true./'" );
  foreach my $clmopts ( @clmoptions ) {
     my @clmres = ( "10x15", "4x5", "360x720cru", "0.9x1.25", "1.9x2.5", "ne3np4.pg3", "ne16np4.pg3", "ne30np4.pg3", "C96", "mpasa120" );
     foreach my $res ( @clmres ) {
        $options = "-res $res -envxml_dir . ";
        &make_env_run( );
        eval{ system( "$bldnml $options $clmopts > $tempfile 2>&1 " ); };
        is( $@, '', "$options $clmopts" );
        $cfiles->checkfilesexist( "$options $clmopts", $mode );
        $cfiles->shownmldiff( "default", "standard" );
        if ( defined($opts{'compare'}) ) {
           $cfiles->doNOTdodiffonfile( "$tempfile", "$options $clmopts", $mode );
           $cfiles->comparefiles( "$options $clmopts", $mode, $opts{'compare'} );
        }
        if ( defined($opts{'generate'}) ) {
           $cfiles->copyfiles( "$options $clmopts", $mode );
        }
        &cleanup();
     }
  }
  my @clmoptions = ( "-bgc bgc -envxml_dir .",
                     "-bgc sp -envxml_dir .", );
  foreach my $clmopts ( @clmoptions ) {
     my @clmres = ( "ne16np4.pg3" );
     foreach my $res ( @clmres ) {
        $options = "-res $res -envxml_dir . ";
        &make_env_run( );
        eval{ system( "$bldnml $options $clmopts > $tempfile 2>&1 " ); };
        is( $@, '', "$options $clmopts" );
        $cfiles->checkfilesexist( "$options $clmopts", $mode );
        $cfiles->shownmldiff( "default", "standard" );
        if ( defined($opts{'compare'}) ) {
           $cfiles->doNOTdodiffonfile( "$tempfile", "$options $clmopts", $mode );
           $cfiles->comparefiles( "$options $clmopts", $mode, $opts{'compare'} );
        }
        if ( defined($opts{'generate'}) ) {
           $cfiles->copyfiles( "$options $clmopts", $mode );
        }
        &cleanup();
     }
  }
  my $clmopts = "-bgc bgc -crop";
  my $res = "1.9x2.5";
  $options = "-res $res -namelist '&a irrigate=.true./' -crop -envxml_dir .";
  &make_env_run();
  eval{ system( "$bldnml $options $clmopts  > $tempfile 2>&1 " ); };
  is( $@, '', "$options $clmopts" );
  $cfiles->checkfilesexist( "$options $clmopts", $mode );
  $cfiles->shownmldiff( "default", "standard" );
  if ( defined($opts{'compare'}) ) {
     $cfiles->doNOTdodiffonfile( "$tempfile", "$options $clmopts", $mode );
     $cfiles->comparefiles( "$options $clmopts", "$mode", $opts{'compare'} );
  }
  if ( defined($opts{'generate'}) ) {
     $cfiles->copyfiles( "$options $clmopts", $mode );
  }
  &cleanup();
  # Run FATES mode for several resolutions and configurations
  my $clmoptions = "-bgc fates -envxml_dir . -no-megan";
  my @clmres = ( "4x5", "1.9x2.5" );
  foreach my $res ( @clmres ) {
     $options = "-res $res -clm_start_type cold";
     my @edoptions = ( "-use_case 2000_control",
                       "-use_case 1850_control",
                       "",
                       "-namelist \"&a use_lch4=.true.,use_nitrif_denitrif=.true./\"",
                       "-clm_accelerated_spinup on"
                     );
     foreach my $edop (@edoptions ) {
        if ( $res eq "5x5_amazon" && ($edop =~ /1850_control/) ) {
           next;
        }
        &make_env_run( );
        eval{ system( "$bldnml $options $clmoptions $edop  > $tempfile 2>&1 " ); };
        is( $@, '', "$options $edop" );
        $cfiles->checkfilesexist( "$options $clmoptions $edop", $mode );
        $cfiles->shownmldiff( "default", "standard" );
        if ( defined($opts{'compare'}) ) {
           $cfiles->doNOTdodiffonfile( "$tempfile", "$options $clmoptions $edop", $mode );
           $cfiles->comparefiles( "$options $clmoptions $edop", $mode, $opts{'compare'} );
        }
        if ( defined($opts{'generate'}) ) {
           $cfiles->copyfiles( "$options $clmoptions $edop", $mode );
        }
        &cleanup();
     }
  }
}
#
# Run over the different lnd_tuning modes
#
my $res = "0.9x1.25";
my $mask = "gx1v7";
my $simyr = "1850";
foreach my $phys ( "clm4_5", "clm5_0", "clm6_0" ) {
  my $mode = "-phys $phys";
  &make_config_cache($phys);
  my @forclist = ();
  @forclist = ( "CRUJRA2024", "CRUv7", "GSWP3v1", "cam7.0", "cam6.0", "cam5.0", "cam4.0" );
  foreach my $forc ( @forclist ) {
     foreach my $bgc ( "sp", "bgc" ) {
        my $lndtuningmode = "${phys}_${forc}";
        if ( $lndtuningmode eq "clm6_0_CRUv7" or
             $lndtuningmode eq "clm4_5_CRUJRA2024") {
           next;
        }
        my $clmoptions = "-res $res -mask $mask -sim_year $simyr -envxml_dir . -lnd_tuning_mod $lndtuningmode -bgc $bgc";
        &make_env_run( );
        eval{ system( "$bldnml $clmoptions > $tempfile 2>&1 " ); };
        is( $@, '', "$clmoptions" );
        $cfiles->checkfilesexist( "$clmoptions", $mode );
        $cfiles->shownmldiff( "default", "standard" );
        if ( defined($opts{'compare'}) ) {
           $cfiles->doNOTdodiffonfile( "$tempfile", "$clmoptions", $mode );
           $cfiles->comparefiles( "$clmoptions", $mode, $opts{'compare'} );
        }
        if ( defined($opts{'generate'}) ) {
           $cfiles->copyfiles( "$clmoptions", $mode );
        }
        &cleanup();
     }
  }
}
&cleanup();

system( "/bin/rm $finidat" );

print "\n==================================================\n";
print " Dumping output  \n";
print "==================================================\n";

$xFail->parseOutput($captOut);

print "Successfully ran all testing for build-namelist\n\n";

&cleanup( "config" );
system( "/bin/rm $tempfile" );

sub cleanup {
#
# Cleanup files created
#
  my $type = shift;

  print "Cleanup files created\n";
  if ( defined($type) ) {
     if ( $type eq "config" ) {
        system( "/bin/rm config_cache.xml" );
     }
  } else {
     system( "/bin/rm $tempfile *_in" );
  }
}

