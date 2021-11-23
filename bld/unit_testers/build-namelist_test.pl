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
    my %env_vars = ( DIN_LOC_ROOT=>"MYDINLOCROOT", GLC_TWO_WAY_COUPLING=>"FALSE" );
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
<entry id="phys" value="$phys" list="" valid_values="clm4_5,clm5_0,clm5_1">Specifies clm physics</entry>
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
   $inputdata_rootdir="/glade/p/cesm/cseg/inputdata";
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
my $ntests = 1844;
if ( defined($opts{'compare'}) ) {
   $ntests += 1254;
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

my $DOMFILE = "$inputdata_rootdir/atm/datm7/domain.lnd.T31_gx3v7.090928.nc";
my $real_par_file = "user_nl_ctsm_real_parameters";
my $bldnml = "../build-namelist -verbose -csmdata $inputdata_rootdir -configuration clm -structure standard -glc_nec 10 -no-note -output_reals $real_par_file";
if ( $opts{'test'} ) {
   $bldnml .= " -test";
}

my $tempfile = "temp_file.txt";
if ( -f $tempfile ) {
  system( "/bin/rm $tempfile" );
}

my @files = ( "lnd_in", $tempfile, $real_par_file );
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
   $options .= " -res 0.9x1.25 -ssp_rcp SSP1-2.6 -envxml_dir .";

   &make_env_run();
   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "options: $options" );
      $cfiles->checkfilesexist( "default", $mode );
      $cfiles->copyfiles( "most_options", $mode );
   # Compare to default
      $cfiles->doNOTdodiffonfile( "lnd_in",    "default", $mode );
      $cfiles->doNOTdodiffonfile( "$real_par_file", "default", $mode );
      $cfiles->doNOTdodiffonfile( "$tempfile", "default", $mode );
      $cfiles->comparefiles( "default", $mode );
   # Compare to baseline
   if ( defined($opts{'compare'}) ) {
      $cfiles->dodiffonfile(      "lnd_in",    "most_options", $mode );
      $cfiles->dodiffonfile( "$real_par_file", "most_options", $mode );
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
foreach my $driver ( "mct", "nuopc" ) {
   print "   For $driver driver\n\n";
   # configuration, structure, irrigate, verbose, clm_demand, ssp_rcp, test, sim_year, use_case
   foreach my $options ( "-configuration nwp",
                         "-structure fast",
                         "-namelist '&a irrigate=.true./'", "-verbose", "-ssp_rcp SSP1-2.6", "-test", "-sim_year 1850",
                         "-namelist '&a use_lai_streams=.true.,use_soil_moisture_streams=.true./'",
                         "-use_case 1850_control",
                         "-res 1x1pt_US-UMB -clm_usr_name 1x1pt_US-UMB -namelist '&a fsurdat=\"/dev/null\"/'",
                         "-res 1x1_brazil",
                         "-clm_start_type startup", "-namelist '&a irrigate=.false./' -crop -bgc bgc",
                         "-envxml_dir . -infile myuser_nl_clm",
                         "-ignore_ic_date -clm_start_type branch -namelist '&a nrevsn=\"thing.nc\"/' -bgc bgc -crop",
                         "-clm_start_type branch -namelist '&a nrevsn=\"thing.nc\",use_init_interp=T/'",
                         "-ignore_ic_date -clm_start_type startup -namelist '&a finidat=\"thing.nc\"/' -bgc bgc -crop",
                        ) {
      my $file = $startfile;
      &make_env_run();
      my $base_options = "-res 0.9x1.25 -envxml_dir . -driver $driver";
      if ( $driver eq "mct" ) {
         $base_options = "$base_options -lnd_frac $DOMFILE";
      } else {
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
         $cfiles->dodiffonfile( "$real_par_file", "$base_options $options", $mode );
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
my $phys = "clm5_1";
$mode = "-phys $phys";
&make_config_cache($phys);
my $neondir      = "../../cime_config/usermods_dirs/NEON";
foreach my $site ( "ABBY", "BLAN", "CPER", "DEJU", "GRSM", "HEAL", "KONA", "LENO", "NIWO",
                   "ONAQ", "PUUM", "SERC", "SRER", "TALL", "TREE", "WOOD", "BARR", "BONA",
                   "DCFS", "DELA", "GUAN", "JERC", "KONZ", "MLBS", "NOGP", "ORNL", "RMNP",
                   "SJER", "STEI", "TEAK", "UKFS", "WREF", "BART", "CLBJ", "DSNY", "HARV",
                   "JORN", "LAJA", "MOAB", "OAES", "OSBS", "SCBI", "SOAP", "STER", "TOOL",
                   "UNDE", "YELL"
 ) {
   &make_env_run();
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
   # Now run  the site
   #
   my $options = "-res CLM_USRDAT -clm_usr_name NEON -no-megan -bgc bgc -sim_year 2000 -infile $namelistfile";
   eval{ system( "$bldnml -envxml_dir . $options > $tempfile 2>&1 " ); };
   is( $@, '', "options: $options" );
   $cfiles->checkfilesexist( "$options", $mode );
   $cfiles->shownmldiff( "default", $mode );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->dodiffonfile(      "lnd_in",    "$options", $mode );
      $cfiles->dodiffonfile( "$real_par_file", "$options", $mode );
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
foreach my $phys ( "clm4_5", "clm5_0" ) {
   $mode = "-phys $phys";
   &make_config_cache($phys);
   foreach my $options (
                      "-res ne0np4.ARCTIC.ne30x4 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=19790101/' -lnd_tuning_mode ${phys}_cam6.0",
                      "-res ne0np4.ARCTICGRIS.ne30x8 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=19790101/' -lnd_tuning_mode ${phys}_cam6.0",
                      "-res 1.9x2.5 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=19790101/' -lnd_tuning_mode ${phys}_cam6.0",
                      "-res 0.9x1.25 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=19790101/' -lnd_tuning_mode ${phys}_cam6.0",
                      "-res 0.9x1.25 -bgc bgc -crop -use_case 20thC_transient -namelist '&a start_ymd=19500101/' -lnd_tuning_mode ${phys}_cam6.0",
                      "-res ne0np4CONUS.ne30x8 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=20130101/' -lnd_tuning_mode ${phys}_cam6.0",
                      "-res 1.9x2.5 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=20030101/' -lnd_tuning_mode ${phys}_cam6.0",
                      "-res 1.9x2.5 -bgc sp -use_case 2010_control -namelist '&a start_ymd=20100101/' -lnd_tuning_mode ${phys}_cam6.0",
                      "-res 1x1_brazil -bgc fates -no-megan -use_case 2000_control -lnd_tuning_mode ${phys}_CRUv7",
                      "-res C192 -bgc sp -use_case 2010_control -namelist '&a start_ymd=20100101/' -lnd_tuning_mode ${phys}_cam6.0",
                      "-res ne0np4.ARCTIC.ne30x4 -bgc sp -use_case 20thC_transient -namelist '&a start_ymd=20130101/' -lnd_tuning_mode ${phys}_cam6.0",
                     ) {
      &make_env_run();
      eval{ system( "$bldnml -envxml_dir . $options > $tempfile 2>&1 " ); };
      is( $@, '', "options: $options" );
      $cfiles->checkfilesexist( "$options", $mode );
      $cfiles->shownmldiff( "default", $mode );
      if ( defined($opts{'compare'}) ) {
         $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
         $cfiles->dodiffonfile(      "lnd_in",    "$options", $mode );
         $cfiles->dodiffonfile( "$real_par_file", "$options", $mode );
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
                      "-bgc bgc -use_case 1850-2100_SSP1-2.6_transient -namelist '&a start_ymd=20100101/'",
                      "-bgc sp  -use_case 1850-2100_SSP2-4.5_transient -namelist '&a start_ymd=18501223/'",
                      "-bgc bgc -use_case 1850-2100_SSP3-7.0_transient -namelist '&a start_ymd=20701029/'",
                      "-bgc fates  -use_case 2000_control -no-megan",
                      "-bgc fates  -use_case 20thC_transient -no-megan",
                      "-bgc fates  -use_case 1850_control -no-megan -namelist '&a use_fates_sp=T/'",
                      "-bgc sp  -use_case 2000_control -res 0.9x1.25 -namelist '&a use_soil_moisture_streams = T/'",
                      "-bgc bgc -use_case 1850-2100_SSP5-8.5_transient -namelist '&a start_ymd=19101023/'",
                      "-bgc bgc -use_case 2000_control -namelist \"&a fire_method='nofire'/\" -crop",
                      "-res 0.9x1.25 -bgc sp -use_case 1850_noanthro_control -drydep -fire_emis",
                      "-res 0.9x1.25 -bgc bgc -use_case 1850_noanthro_control -drydep -fire_emis -light_res 360x720",
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
      $cfiles->dodiffonfile( "$real_par_file", "$options", $mode );
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
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "clm_demand on finidat"     =>{ options=>"-clm_demand finidat -envxml_dir .",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "blank IC file, not cold"   =>{ options=>"-clm_start_type startup -envxml_dir .",
                                     namelst=>"finidat=' '",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "startup without interp"    =>{ options=>"-clm_start_type startup -envxml_dir . -bgc sp -sim_year 1850",
                                     namelst=>"use_init_interp=.false., start_ymd=19200901",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "use_crop without -crop"    =>{ options=>" -envxml_dir .",
                                     namelst=>"use_crop=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "soilm_stream wo use"       =>{ options=>"-res 0.9x1.25 -envxml_dir .",
                                     namelst=>"use_soil_moisture_streams = .false.,stream_fldfilename_soilm='missing_file'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "clm50CNDVwtransient"       =>{ options=>" -envxml_dir . -use_case 20thC_transient -dynamic_vegetation -res 10x15 -ignore_warnings",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "decomp_without_cn"         =>{ options=>" -envxml_dir . -bgc sp",
                                     namelst=>"soil_decomp_method='CENTURYKoven2013'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "bgc_with_no_decomp"        =>{ options=>" -envxml_dir . -bgc bgc",
                                     namelst=>"soil_decomp_method='None'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "reseed without CN"         =>{ options=>" -envxml_dir . -bgc sp",
                                     namelst=>"reseed_dead_plants=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "onset_threh w SP"          =>{ options=>" -envxml_dir . -bgc sp",
                                     namelst=>"onset_thresh_depends_on_veg=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_1",
                                   },
     "dribble_crphrv w/o CN"     =>{ options=>" -envxml_dir . -bgc sp",
                                     namelst=>"dribble_crophrv_xsmrpool_2atm=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "dribble_crphrv w/o crop"   =>{ options=>" -envxml_dir . -bgc bgc -no-crop",
                                     namelst=>"dribble_crophrv_xsmrpool_2atm=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "CNDV with flanduse_timeseries - clm4_5"=>{ options=>"-bgc bgc -dynamic_vegetation -envxml_dir . -ignore_warnings",
                                     namelst=>"flanduse_timeseries='my_flanduse_timeseries_file.nc'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "use_cndv=T without bldnml op"=>{ options=>"-bgc bgc -envxml_dir . -ignore_warnings",
                                     namelst=>"use_cndv=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "use_cndv=F with dyn_veg op"=>{ options=>"-bgc bgc -dynamic_vegetation -envxml_dir . -ignore_warnings",
                                     namelst=>"use_cndv=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "crop with use_crop false"  =>{ options=>"-crop -bgc bgc -envxml_dir .",
                                     namelst=>"use_crop=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "crop without CN"           =>{ options=>"-crop -bgc sp -envxml_dir .",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "toosmall soil w trans"     =>{ options=>"-envxml_dir .",
                                     namelst=>"toosmall_soil=10, dyn_transient_pfts=T",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "toosmall lake w trans"     =>{ options=>"-envxml_dir .",
                                     namelst=>"toosmall_lake=10, dyn_transient_pfts=T",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "toosmall crop w trans"     =>{ options=>"-bgc bgc -crop -envxml_dir .",
                                     namelst=>"toosmall_crop=10, dyn_transient_pfts=T",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "toosmall wetl w trans"     =>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"toosmall_wetland=10, dyn_transient_pfts=T",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "toosmall glc  w trans"     =>{ options=>"-bgc sp  -envxml_dir .",
                                     namelst=>"toosmall_glacier=10, dyn_transient_pfts=T",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "toosmall urban w trans"    =>{ options=>"-bgc sp  -envxml_dir .",
                                     namelst=>"toosmall_urban=10, dyn_transient_pfts=T",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "collapse_urban w trans"    =>{ options=>"-bgc sp  -envxml_dir .",
                                     namelst=>"collapse_urban=T, dyn_transient_crops=T",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "n_dom_landunits w trans"    =>{ options=>"-bgc sp  -envxml_dir .",
                                     namelst=>"n_dom_landunits=2, dyn_transient_crops=T",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "n_dom_pfts w trans"         =>{ options=>"-bgc sp  -envxml_dir .",
                                     namelst=>"n_dom_pfts=2, dyn_transient_crops=T",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "baset_map without crop"    =>{ options=>"-bgc bgc -envxml_dir . -no-crop",
                                     namelst=>"baset_mapping='constant'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "mapvary var w/o varymap"   =>{ options=>"-crop -bgc bgc -envxml_dir . -crop",
                                     namelst=>"baset_mapping='constant', baset_latvary_slope=1.0, baset_latvary_intercept=10.0",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "grainproductWOcrop"       =>{ options=>"-bgc bgc -no-crop -envxml_dir .",
                                    namelst=>"use_grainproduct=.true.",
                                    GLC_TWO_WAY_COUPLING=>"FALSE",
                                    phys=>"clm4_5",
                                   },
     "interp without finidat"    =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"use_init_interp=.true. finidat=' '",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "sp and c13"                =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"use_c13=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "sp and c14"                =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"use_c14=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "bombspike no c14"          =>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"use_c14=.false. use_c14_bombspike=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "use c13 timeseries no cn"  =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"use_c13_timeseries=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "use c13 timeseries no c13"=>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"use_c13=.false. use_c13_timeseries=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "bombspike no cn"           =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"use_c14_bombspike=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "lightres no cn"            =>{ options=>"-bgc sp -envxml_dir . -light_res 360x720",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "spno-fire"                 =>{ options=>"-bgc sp -envxml_dir . -use_case 2000_control",
                                     namelst=>"fire_method='nofire'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "lightres no fire"          =>{ options=>"-bgc bgc -envxml_dir . -light_res 360x720",
                                     namelst=>"fire_method='nofire'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "lightres none bgc"         =>{ options=>"-bgc bgc -envxml_dir . -light_res none",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "lightresnotnone-nofire"    =>{ options=>"-bgc bgc -envxml_dir . -light_res 94x192",
                                     namelst=>"fire_method='nofire'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "lightresnonenofirelightfil"=>{ options=>"-bgc bgc -envxml_dir . -light_res none",
                                     namelst=>"fire_method='nofire',stream_fldfilename_lightng='build-namelist_test.pl'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "lightrescontradictlightfil"=>{ options=>"-bgc bgc -envxml_dir . -light_res 360x720",
                                     namelst=>"stream_fldfilename_lightng='build-namelist_test.pl'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "finundated and not methane"=>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"use_lch4=.false.,finundation_method='h2osfc'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "use_cn=true bgc=sp"        =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"use_cn=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "freeliv wo fun"            =>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"freelivfix_intercept=9.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "use_cn=false bgc=bgc"      =>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"use_cn=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "lower=aqu-45 with/o Zeng"  =>{ options=>"-envxml_dir .",
                                     namelst=>"lower_boundary_condition=4,soilwater_movement_method=1,use_bedrock=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "Zeng w lower=flux"         =>{ options=>"-envxml_dir .",
                                     namelst=>"lower_boundary_condition=1,soilwater_movement_method=0,use_bedrock=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "Zeng w lower=zeroflux"     =>{ options=>"-envxml_dir .",
                                     namelst=>"lower_boundary_condition=2,soilwater_movement_method=0",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "Zeng w lower=table"        =>{ options=>"-envxml_dir .",
                                     namelst=>"lower_boundary_condition=3,soilwater_movement_method=0,use_bedrock=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "use_vic=F with -vic op"    =>{ options=>"-vichydro -envxml_dir .",
                                     namelst=>"use_vichydro=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "-vic with l_bnd=flux"      =>{ options=>"-vichydro -envxml_dir .",
                                     namelst=>"lower_boundary_condition=1",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "-vic with l_bnd=zeroflux"  =>{ options=>"-vichydro -envxml_dir .",
                                     namelst=>"lower_boundary_condition=2",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "-vic with origflag=1"      =>{ options=>"-vichydro -envxml_dir .",
                                     namelst=>"origflag=1",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "l_bnd=flux with origflag=0"=>{ options=>"-envxml_dir .",
                                     namelst=>"origflag=0, lower_boundary_condition=1",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "l_bnd=zflux with origflag=0"=>{ options=>"-envxml_dir .",
                                     namelst=>"origflag=0, lower_boundary_condition=2",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "bedrock with l_bnc=flux"   =>{ options=>"-envxml_dir .",
                                     namelst=>"use_bedrock=.true., lower_boundary_condition=1",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "bedrock with l_bnc=tabl"   =>{ options=>"-envxml_dir .",
                                     namelst=>"use_bedrock=.true., lower_boundary_condition=3",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "bedrock with l_bnc=aqui"   =>{ options=>"-envxml_dir .",
                                     namelst=>"use_bedrock=.true., lower_boundary_condition=4",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "zengdeck with l_bnc=flux"  =>{ options=>"-envxml_dir .",
                                     namelst=>"soilwater_movement_method=0, lower_boundary_condition=1",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "zengdeck with l_bnc=z-flux"=>{ options=>"-envxml_dir .",
                                     namelst=>"soilwater_movement_method=0, lower_boundary_condition=2",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "zengdeck with l_bnc=tabl"  =>{ options=>"-envxml_dir .",
                                     namelst=>"soilwater_movement_method=0, lower_boundary_condition=3",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "l_bnd=tabl with h2osfcfl=0"=>{ options=>"-envxml_dir .",
                                     namelst=>"h2osfcflag=0, lower_boundary_condition=3",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "l_bnd=flux with h2osfcfl=0"=>{ options=>"-envxml_dir .",
                                     namelst=>"h2osfcflag=0, lower_boundary_condition=1",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "l_bnd=zflux with h2osfcfl=0"=>{ options=>"-envxml_dir .",
                                     namelst=>"h2osfcflag=0, lower_boundary_condition=2",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "h2osfcfl=0 with clm5.0"    =>{ options=>"-envxml_dir .",
                                     namelst=>"h2osfcflag=0",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "45bad lnd_tuning_mode value" =>{ options=>"-lnd_tuning_mode clm5_0_GSWP3  -envxml_dir .",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "50bad lnd_tuning_mode value" =>{ options=>"-lnd_tuning_mode clm4_5_CRUNCEP  -envxml_dir .",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "bgc_spinup without cn"     =>{ options=>"-clm_accelerated_spinup on -bgc sp -envxml_dir .",
                                     namelst=>"spinup_state=1",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "spinup=1 without bldnml op"=>{ options=>"-clm_accelerated_spinup off -bgc bgc -envxml_dir .",
                                     namelst=>"spinup_state=1",,
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "bgc_spinup without cn"     =>{ options=>"-clm_accelerated_spinup on -bgc sp -envxml_dir .",
                                     namelst=>"spinup_state=1",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "baseflow w aquifer"        =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"baseflow_scalar=1.0, lower_boundary_condition=4,use_bedrock=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "baseflow w table"          =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"baseflow_scalar=1.0, lower_boundary_condition=3,use_bedrock=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "br_root and bgc=sp"        =>{ options=>"-bgc sp -envxml_dir .",
                                     namelst=>"br_root=1.0",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "both co2_type and on nml"  =>{ options=>"-co2_type constant -envxml_dir .",
                                     namelst=>"co2_type='prognostic'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "both lnd_frac and on nml"  =>{ options=>"-driver mct -lnd_frac $DOMFILE -envxml_dir .",
                                     namelst=>"fatmlndfrc='frac.nc'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "lnd_frac set to UNSET"     =>{ options=>"-driver mct -lnd_frac UNSET -envxml_dir .",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_1",
                                   },
     "lnd_frac set but nuopc"    =>{ options=>"-driver nuopc -lnd_frac $DOMFILE -envxml_dir .",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_1",
                                   },
     "lnd_frac not set but lilac"=>{ options=>"-driver nuopc -lilac -envxml_dir . -lnd_frac UNSET",
                                     namelst=>"fsurdat='surfdata.nc'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_1",
                                   },
     "fatmlndfrc set but nuopc"  =>{ options=>"-driver nuopc -envxml_dir .",
                                     namelst=>"fatmlndfrc='frac.nc'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_1",
                                   },
     "force_send but not nuopc"  =>{ options=>"-driver mct -lnd_frac $DOMFILE -envxml_dir .",
                                     namelst=>"force_send_to_atm = .false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_1",
                                   },
     "branch but NO nrevsn"      =>{ options=>"-clm_start_type branch -envxml_dir .",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "glc_nec inconsistent"      =>{ options=>"-envxml_dir .",
                                     namelst=>"maxpatch_glc=5",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "NoGLCMec"                  =>{ options=>"-envxml_dir . -glc_nec 0",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "UpdateGlcContradict"       =>{ options=>"-envxml_dir .",
                                     namelst=>"glc_do_dynglacier=.false.",
                                     GLC_TWO_WAY_COUPLING=>"TRUE",
                                     phys=>"clm4_5",
                                   },
     "useFATESContradict"        =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "useFATESContradict2"       =>{ options=>"-envxml_dir . -no-megan",
                                     namelst=>"use_fates=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "useFATESWCN"               =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_cn=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "useFATESWcrop"             =>{ options=>"-bgc fates -envxml_dir . -no-megan -crop",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_1",
                                   },
     "useFATESWcreatecrop"       =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"create_crop_landunit=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "useFATESWbMH"              =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_biomass_heat_storage=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_1",
                                   },
     "FireNoneButFATESfireon"    =>{ options=>"-bgc fates -envxml_dir . -no-megan -light_res none",
                                     namelst=>"fates_spitfire_mode=4",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_1",
                                   },
     "FireNoneButBGCfireon"    =>{ options=>"-bgc bgc -envxml_dir . -light_res none",
                                     namelst=>"fire_method='li2021gswpfrc'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_1",
                                   },
     "createcropFalse"           =>{ options=>"-bgc bgc -envxml_dir . -no-megan",
                                     namelst=>"create_crop_landunit=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "usespitfireButNOTFATES"    =>{ options=>"-envxml_dir . -no-megan",
                                     namelst=>"fates_spitfire_mode=1",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "useloggingButNOTFATES"     =>{ options=>"-envxml_dir . -no-megan",
                                     namelst=>"use_fates_logging=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "useinventorybutnotfile"    =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_inventory_init=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "inventoryfileDNE"          =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_inventory_init=.true., fates_inventory_ctrl_filename='zztop'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "useMEGANwithFATES"         =>{ options=>"-bgc fates -envxml_dir . -megan",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "useFATESSPWONOCOMP"        =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_fates_sp=T,use_fates_nocomp=F",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "useFATESTRANSWdynPFT"      =>{ options=>"-bgc fates -envxml_dir . -use_case 20thC_transient -no-megan",
                                     namelst=>"do_transient_pfts=T",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "useHYDSTwithFATES"         =>{ options=>"-bgc fates -envxml_dir . -no-megan",
                                     namelst=>"use_hydrstress=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "useHYDSTwithdynroot"       =>{ options=>"-bgc bgc -envxml_dir . -megan",
                                     namelst=>"use_hydrstress=.true., use_dynroot=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "specWOfireemis"            =>{ options=>"-envxml_dir . -no-fire_emis",
                                     namelst=>"fire_emis_specifier='bc_a1 = BC'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "elevWOfireemis"            =>{ options=>"-envxml_dir . -no-fire_emis",
                                     namelst=>"fire_emis_elevated=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "noanthro_w_crop"            =>{ options=>"-envxml_dir . -res 0.9x1.25 -bgc bgc -crop -use_case 1850_noanthro_control",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "noanthro_w_irrig"           =>{ options=>"-envxml_dir . -res 0.9x1.25 -bgc bgc -use_case 1850_noanthro_control",
                                     namelst=>"irrigate=T",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "spdotransconflict"          =>{ options=>"-envxml_dir . -bgc sp -use_case 20thC_transient",
                                     namelst=>"do_transient_pfts=T,do_transient_crops=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "nocropwfert"                =>{ options=>"-envxml_dir . -bgc sp -no-crop",
                                     namelst=>"use_fertilizer=T",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "lmr1WOcn"                   =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"leafresp_method=1",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "lmr2WOcn"                   =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"leafresp_method=2",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "lmr0Wcn"                    =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"leafresp_method=0",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "nofireButSetcli_scale"     =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"fire_method='nofire', cli_scale=5.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "nocnButSetrh_low"          =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"rh_low=5.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "funWOcn"                   =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"use_fun=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "flexCNWOcn"                =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"use_flexibleCN=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "flexCNFUNwcarbonresp"      =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_flexibleCN=.true.,use_FUN=.true.,carbon_resp_opt=1",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "funWOnitrif"               =>{ options=>"-envxml_dir .",
                                     namelst=>"use_fun=.true., use_nitrif_denitrif=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "SPModeWNitrifNMethane"     =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"use_lch4=.true., use_nitrif_denitrif=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "knitrmaxWOnitrif"          =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_nitrif_denitrif=.false., k_nitr_max=1.0",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "respcoefWOnitrif"          =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_nitrif_denitrif=.false., denitrif_respiration_coefficient=1.0",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "respexpWOnitrif"           =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_nitrif_denitrif=.false., denitrif_respiration_exponent=1.0",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "lunaWSPandlnctrue"         =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"use_luna=.true., lnc_opt=.true.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "NOlunabutsetJmaxb1"        =>{ options=>"-envxml_dir . -bgc sp",
                                     namelst=>"use_luna=.false., jmaxb1=1.0",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "envxml_not_dir"            =>{ options=>"-envxml_dir myuser_nl_clm",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "envxml_emptydir"           =>{ options=>"-envxml_dir xFail",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
               );
foreach my $key ( keys(%failtest) ) {
   print( "$key\n" );
   &make_config_cache($failtest{$key}{"phys"});
   my $options  = $failtest{$key}{"options"};
   my $namelist = $failtest{$key}{"namelst"};
   &make_env_run( GLC_TWO_WAY_COUPLING=>$failtest{$key}{"GLC_TWO_WAY_COUPLING"} );
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
     "coldwfinidat"              =>{ options=>"-envxml_dir . -clm_start_type cold",
                                     namelst=>"finidat = 'testfile.nc'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "bgcspin_w_suplnitro"       =>{ options=>"-envxml_dir . -bgc bgc -clm_accelerated_spinup on",
                                     namelst=>"suplnitro='ALL'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "bgc=bgc WO nitrif_denit"   =>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"use_nitrif_denitrif=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "methane off W nitrif_denit"=>{ options=>"-bgc bgc -envxml_dir .",
                                     namelst=>"use_nitrif_denitrif=.true.,use_lch4=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_1",
                                   },
     "soilm_stream w transient"  =>{ options=>"-res 0.9x1.25 -envxml_dir . -use_case 20thC_transient",
                                     namelst=>"use_soil_moisture_streams=T,soilm_tintalgo='linear'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "missing_ndep_file"         =>{ options=>"-envxml_dir . -bgc bgc -ssp_rcp SSP5-3.4",
                                     namelst=>"",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_0",
                                   },
     "bad_megan_spec"            =>{ options=>"-envxml_dir . -bgc bgc -megan",
                                     namelst=>"megan_specifier='ZZTOP=zztop'",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm4_5",
                                   },
     "FUN_wo_flexCN"             =>{ options=>"-envxml_dir . -bgc bgc",
                                     namelst=>"use_fun=.true.,use_flexiblecn=.false.",
                                     GLC_TWO_WAY_COUPLING=>"FALSE",
                                     phys=>"clm5_1",
                                   },
               );
foreach my $key ( keys(%warntest) ) {
   print( "$key\n" );
   &make_config_cache($warntest{$key}{"phys"});
   my $options  = $warntest{$key}{"options"};
   my $namelist = $warntest{$key}{"namelst"};
   &make_env_run( GLC_TWO_WAY_COUPLING=>$warntest{$key}{"GLC_TWO_WAY_COUPLING"} );
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

#
# Loop over all physics versions
#
foreach my $phys ( "clm4_5", "clm5_0", "clm5_1" ) {
$mode = "-phys $phys";
&make_config_cache($phys);

print "\n==================================================\n";
print "Test ALL resolutions with SP\n";
print "==================================================\n";

# Check for ALL resolutions with CLM50SP
my $reslist = `../queryDefaultNamelist.pl -res list -s`;
my @resolutions = split( / /, $reslist );
my @regional;
foreach my $res ( @resolutions ) {
   chomp($res);
   print "=== Test $res === \n";
   my $options  = "-res $res -bgc sp -envxml_dir .";

   # Regional single point resolutions
   if ( $res =~ /^([0-9]+x[0-9]+_[a-zA-Z]+)$/ ) {
      push( @regional, $res );
      next;
   # Resolutions for mksurfdata mapping
   } elsif ( $res eq "0.5x0.5"     ||
             $res eq "0.25x0.25"   ||
             $res eq "3x3min"      ||
             $res eq "5x5min"      ||
             $res eq "10x10min"    ||
             $res eq "0.125x0.125" ||
             $res eq "0.33x0.33"   ||
             $res eq "1km-merge-10min" ) {
      next;
   # Resolutions that were supported in clm40 but NOT clm45/clm50
   } elsif ( $res eq "ne240np4"    ||
             $res eq "ne60np4"     ||
             $res eq "ne4np4"      ||
             $res eq "2.5x3.33"    ||
             $res eq "0.23x0.31"   ||
             $res eq "0.47x0.63"   ||
             $res eq "94x192"      ||
             $res eq "8x16"        ||
             $res eq "32x64"       ||
             $res eq "128x256"     ||
             $res eq "360x720cru"  ||
             $res eq "512x1024" ) {
      next;
   # Resolutions not supported on release branch
   } elsif ( $res eq "ne120np4"    ||
             $res eq "conus_30_x8" ) {
      next;
   }

   &make_env_run();
   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );

   $cfiles->checkfilesexist( "$options", $mode );

   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->dodiffonfile( "$real_par_file", "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }

   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup(); print "\n";
}

print "\n==================================================\n";
print " Test important resolutions for BGC\n";
print "==================================================\n";

my @resolutions = ( "4x5", "10x15", "ne30np4", "ne16np4", "1.9x2.5", "0.9x1.25" );
my @regional;
my $nlbgcmode = "bgc";
my $mode = "$phys-$nlbgcmode";
foreach my $res ( @resolutions ) {
   chomp($res);
   print "=== Test $res === \n";
   my $options  = "-res $res -envxml_dir . -bgc $nlbgcmode";

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
print " Test all use-cases \n";
print "==================================================\n";

# Run over all use-cases...
my $list = `$bldnml -use_case list 2>&1 | grep "use case"`;
my @usecases;
if ( $list =~ /build-namelist : use cases : (.+)$/ ) {
  my @usecases  = split( / /, $list );
} else {
  die "ERROR:: Trouble getting list of use-cases\n";
}
foreach my $usecase ( @usecases ) {
   $options = "-use_case $usecase  -envxml_dir .";
   &make_env_run();
   eval{ system( "$bldnml $options  > $tempfile 2>&1 " ); };
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
   &cleanup();
}

print "\n==================================================\n";
print "Test crop resolutions \n";
print "==================================================\n";

# Check for crop resolutions
my @crop_res = ( "1x1_numaIA", "1x1_smallvilleIA", "4x5", "10x15", "0.9x1.25", "1.9x2.5", "ne30np4" );
foreach my $res ( @crop_res ) {
   $options = "-bgc bgc -crop -res $res -envxml_dir .";
   &make_env_run();
   eval{ system( "$bldnml $options  > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );
   $cfiles->checkfilesexist( "$options", $mode );
   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->dodiffonfile( "$real_par_file", "$options", $mode );
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
my @use_cases = ( "1850-2100_SSP1-2.6_transient",
                  "1850-2100_SSP2-4.5_transient",
                  "1850-2100_SSP3-7.0_transient",
                  "1850-2100_SSP5-8.5_transient",
                  "1850_control",
                  "2000_control",
                  "2010_control",
                  "20thC_transient",
                 );
foreach my $res ( @glc_res ) {
   foreach my $usecase ( @usecases ) {
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
my @tran_res = ( "0.9x1.25", "1.9x2.5", "ne30np4", "10x15" );
my $usecase  = "20thC_transient";
my $GLC_NEC         = 10;
foreach my $res ( @tran_res ) {
   $options = "-res $res -use_case $usecase -envxml_dir . -namelist '&a start_ymd=18500101/'";
   &make_env_run();
   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );
   $cfiles->checkfilesexist( "$options", $mode );
   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->dodiffonfile( "$real_par_file", "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}
# Transient ssp_rcp scenarios that work
my @tran_res = ( "0.9x1.25", "1.9x2.5", "10x15" );
foreach my $usecase ( "1850_control", "1850-2100_SSP5-8.5_transient", "1850-2100_SSP1-2.6_transient", "1850-2100_SSP3-7.0_transient",
                      "1850-2100_SSP2-4.5_transient" ) {
   my $startymd = undef;
   if ( $usecase eq "1850_control") {
      $startymd = 18500101;
   } else {
      $startymd = 20150101;
   }
   foreach my $res ( @tran_res ) {
      $options = "-res $res -bgc bgc -crop -use_case $usecase -envxml_dir . -namelist '&a start_ymd=$startymd/'";
      &make_env_run();
      eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
      is( $@, '', "$options" );
      $cfiles->checkfilesexist( "$options", $mode );
      $cfiles->shownmldiff( "default", "standard" );
      if ( defined($opts{'compare'}) ) {
         $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
         $cfiles->dodiffonfile( "$real_par_file", "$options", $mode );
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

# The SSP's that fail...
$phys = "clm5_0";
$mode = "-phys $phys";
&make_config_cache($phys);
my $res = "0.9x1.25";
foreach my $usecase ( "1850-2100_SSP4-3.4_transient", "1850-2100_SSP5-3.4_transient", "1850-2100_SSP1-1.9_transient",
                      "1850-2100_SSP4-6.0_transient" ) {
      $options = "-res $res -bgc bgc -crop -use_case $usecase -envxml_dir . -namelist '&a start_ymd=20150101/'";
      &make_env_run();
      eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
      isnt( $?, 0, $usecase );
      system( "cat $tempfile" );
}

print "\n==================================================\n";
print "Test clm4.5/clm5.0/clm5_1 resolutions \n";
print "==================================================\n";

foreach my $phys ( "clm4_5", 'clm5_0', 'clm5_1' ) {
  my $mode = "-phys $phys";
  &make_config_cache($phys);
  my @clmoptions = ( "-bgc bgc -envxml_dir .", "-bgc bgc -envxml_dir . -clm_accelerated_spinup=on", "-bgc bgc -envxml_dir . -light_res 360x720",
                     "-bgc sp -envxml_dir . -vichydro", "-bgc bgc -dynamic_vegetation -ignore_warnings",
                     "-bgc bgc -clm_demand flanduse_timeseries -sim_year 1850-2000 -namelist '&a start_ymd=18500101/'",
                     "-bgc bgc -envxml_dir . -namelist '&a use_c13=.true.,use_c14=.true.,use_c14_bombspike=.true./'" );
  foreach my $clmopts ( @clmoptions ) {
     my @clmres = ( "10x15", "0.9x1.25", "1.9x2.5" );
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
     my @clmres = ( "ne16np4" );
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
  my @clmres = ( "1x1_brazil", "5x5_amazon", "4x5", "1.9x2.5" );
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
# Run over the differen lnd_tuning modes
#
my $res = "0.9x1.25";
my $mask = "gx1v6";
my $simyr = "1850";
foreach my $phys ( "clm4_5", 'clm5_0', 'clm5_1' ) {
  my $mode = "-phys $phys";
  &make_config_cache($phys);
  my @forclist = ();
  if ( $phys == "clm5_1" ) {
    @forclist = ( "GSWP3v1" );
  } else {
    @forclist = ( "CRUv7", "GSWP3v1", "cam6.0" );
  }
  foreach my $forc ( @forclist ) {
     foreach my $bgc ( "sp", "bgc" ) {
        my $lndtuningmode = "${phys}_${forc}";
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
  system( "/bin/rm env_run.xml $real_par_file" );
  if ( defined($type) ) {
     if ( $type eq "config" ) {
        system( "/bin/rm config_cache.xml" );
     }
  } else {
     system( "/bin/rm $tempfile *_in" );
  }
}

