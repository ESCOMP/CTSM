#!/usr/bin/env perl

# Test methods of the build-namelist script.
# Just tests for configure mode using CESM.
# Try to test that all the different options at least work.
# Test that inconsistentcies are appropriately caught.

#########################

use Test::More;

#########################

use strict;
use Getopt::Long;
use NMLTest::CompFiles;

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
     -test                         Use the -test option to make sure datasets exist.

EOF
}


#
# Process command-line options.
#
my %opts = ( help     => 0,
             generate => 0,
             test     => 0,
             compare  => undef,
            );

GetOptions(
    "h|help"     => \$opts{'help'},
    "compare=s"  => \$opts{'compare'},
    "generate"   => \$opts{'generate'},
    "test"       => \$opts{'test'},
)  or usage();

# Give usage message.
usage() if $opts{'help'};

#
# Figure out number of tests that will run
#
my $ntests = 178;
if ( defined($opts{'compare'}) ) {
   $ntests += 99;
}
plan( tests=>$ntests );

# Check for unparsed arguments
if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
}
my $mode = "standard";
system( "../configure -s" );

my $bldnml = "../build-namelist -verbose";
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
eval{ system( "$bldnml > $tempfile 2>&1 " ); };
is( $@, '', "plain build-namelist" );
$cfiles->checkfilesexist( "default", $mode );
# Compare to baseline
if ( defined($opts{'compare'}) ) {
   $cfiles->doNOTdodiffonfile( "$tempfile", "default", $mode );
   $cfiles->comparefiles( "default", $mode, $opts{'compare'} );
}
$cfiles->copyfiles( "default", $mode );
&cleanup();
# Simple test -- run all the list options
foreach my $options ( "clm_demand", "rcp",      "res", 
                      "sim_year",   "use_case" ) {
   eval{ system( "$bldnml -${options} list > $tempfile 2>&1 " ); };
   my $result = `cat $tempfile`;
   my $expect;
   if ( $options !~ /use_case/ ) {
      $expect = "valid values for $options";
   } else {
      $expect = "use cases:";
   }
   $expect    = "/^build-namelist - $expect/";
   like( $result, $expect, "$options list" );
   is( (-f "lnd_in"), undef, "Check that lnd_in file does NOT exist" );
   &cleanup();
}
# Exercise a bunch of options
my $options = "-co2_ppmv 250 -glc_nec 10 -glc_grid gland5 -glc_smb .false.";
   $options .= " -res 0.9x1.25 -rtm off -rtm_tstep 10800 -rcp 2.6";
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
# drydep and megan namelists
my @mfiles = ( "lnd_in", "drv_flds_in", $tempfile );
my $mfiles = NMLTest::CompFiles->new( $cwd, @mfiles );
foreach my $options ( "-drydep", "-megan", "-drydep -megan" ) {
   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
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
# irrig, verbose, clm_demand, rcp, test, sim_year, use_case, l_ncpl
my $startfile = "clmrun.clm2\%inst_string.r.1964-05-27-00000.nc";
my $startnoin = "clmrun.clm2.r.1964-05-27-00000.nc";
my $inst      = "_0009";
foreach my $options ( "-irrig", "-verbose", "-rcp 2.6", "-test", "-sim_year 1850",
                      "-use_case 1850_control", "-l_ncpl 1", 
                      "-clm_startfile $startfile -clm_start_type startup", 
                      "-clm_startfile $startfile  -clm_start_type startup -inst_string $inst", 
                     ) {
   my $file = $startfile;
   if ( $options =~ /-inst_string/ ) {
      $file =~ s/\%inst_string/$inst/;
      system( "touch $file" );
   }
   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "options: $options" );
   $cfiles->checkfilesexist( "$options", $mode );
   system( "diff lnd_in lnd_in.default" );
   $cfiles->shownmldiff( "default", $mode );
   my $finidat = `grep finidat lnd_in`;
   if (      $options =~ /-inst_string/    ) {
      like( $finidat, "/$file/",      "$options" );
      system( "/bin/rm $file" );
   } elsif ( $options =~ /-clm_start_file/ ) {
      like( $finidat, "/$startnoin/", "$options" );
   }
   if ( $options eq "-l_ncpl 1" ) {
      my $dtime = `grep dtime lnd_in`;
      like( $dtime, "/ 86400\$/", "$options" );
   }
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$options", $mode );
      $cfiles->comparefiles( "$options", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$options", $mode );
   }
   &cleanup();
}
# Failure testing, do things that SHOULD fail
my $finidat  = "thing.nc";
system( "touch $finidat" );

my %failtest = ( 
     "coldstart but with IC file"=>{ options=>"-clm_start_type cold",
                                     namelst=>"finidat='$finidat'",
                                   },
     "coldstart but clm_startfile"=>{ options=>"-clm_start_type cold -clm_startfile file.nc",
                                     namelst=>"",
                                   },
     "l_ncpl is zero"            =>{ options=>"-l_ncpl 0",
                                     namelst=>"",
                                   },
     "l_ncpl not integer"        =>{ options=>"-l_ncpl 1.0",
                                     namelst=>"",
                                   },
     "both l_ncpl and dtime"     =>{ options=>"-l_ncpl 24",
                                     namelst=>"dtime=1800",
                                   },
     "both co2_type and on nml"  =>{ options=>"-co2_type constant",
                                     namelst=>"co2_type='prognostic'",
                                   },
     "both lnd_frac and on nml"  =>{ options=>"-lnd_frac domain.nc",
                                     namelst=>"fatmlndfrc='frac.nc'",
                                   },
     "both start_file and finidat"=>{ options=>"-clm_startfile file.nc -clm_start_type startup",
                                     namelst=>"finidat='finidat.nc'",
                                   },
     "both start_file and nrevsn"=>{ options=>"-clm_startfile file.nc -clm_start_type branch",
                                     namelst=>"nrevsn='nrevsn.nc'",
                                   },
     "branch but NO nrevsn"      =>{ options=>"-clm_start_type branch",
                                     namelst=>"",
                                   },
     "glc_nec inconsistent"      =>{ options=>"-glc_nec 10",
                                     namelst=>"maxpatch_glcmec=5",
                                   },
     "glc_smb inconsistent"      =>{ options=>"-glc_nec 10 -glc_smb .true.",
                                     namelst=>"glc_smb=.false.",
                                   },
     "glc_grid inconsistent"     =>{ options=>"-glc_nec 10 -glc_grid gland10",
                                     namelst=>"glc_grid='gland5'",
                                   },
     "rtm inconsistent"          =>{ options=>"-rtm off",
                                     namelst=>"do_rtm=.true.",
                                   },
     "rtm tstep inconsistent"    =>{ options=>"-rtm_tstep 10800",
                                     namelst=>"rtm_nsteps=1",
                                   },
               );
foreach my $key ( keys(%failtest) ) {
   my $options  = $failtest{$key}{"options"};
   my $namelist = $failtest{$key}{"namelst"};
   eval{ system( "$bldnml $options -namelist \"&clmexp $namelist /\" > $tempfile 2>&1 " ); };
   isnt( $?, 0, $key );
   system( "cat $tempfile" );
}
# Check for ALL resolutions with CN
my $mode = "CN";
system( "../configure -s -bgc cn" );
my $reslist = `../queryDefaultNamelist.pl -res list -s`;
my @resolutions = split( / /, $reslist );
my @regional;
foreach my $res ( @resolutions ) {
   $options  = "-res $res";
   if ( $res eq "512x1024" ) { 
      $options .= " -sim_year 1850"; 
   } elsif ( $res =~ /^([0-9]+x[0-9]+_[a-zA-Z]+)$/ ) {
      push( @regional, $res );
      next;
   } elsif ( $res eq "0.5x0.5"     ||
             $res eq "3x3min"      ||
             $res eq "5x5min"      ||
             $res eq "10x10min"    ||
             $res eq "0.33x0.33"  ) {
      next;
   }
   eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );
   $cfiles->checkfilesexist( "$options", $mode );
   system( "diff lnd_in lnd_in.default.standard" );
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
# Run over all use-cases...
my $list = `$bldnml -use_case list 2>&1 | grep "use case"`;
my @usecases;
if ( $list =~ /build-namelist - use cases: (.+)$/ ) {
  my @usecases  = split( / /, $list );
} else {
  die "ERROR:: Trouble getting list of use-cases\n";
}
foreach my $usecase ( @usecases ) {
   $options = "-use_case $usecase ";
   eval{ system( "$bldnml $options  > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );
   $cfiles->checkfilesexist( "$options", $mode );
   system( "diff lnd_in lnd_in.default.standard" );
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

# Run over single-point regional cases
foreach my $res ( @regional ) {
   $mode = "$res";
   system( "../configure -s -sitespf_pt $res" );
   eval{ system( "$bldnml > $tempfile 2>&1 " ); };
   is( $@, '', "$res" );
   $cfiles->checkfilesexist( "$res", $mode );
   system( "diff lnd_in lnd_in.default.standard" );
   $cfiles->shownmldiff( "default", "standard" );
   if ( defined($opts{'compare'}) ) {
      $cfiles->doNOTdodiffonfile( "$tempfile", "$res", $mode );
      $cfiles->comparefiles( "$res", $mode, $opts{'compare'} );
   }
   if ( defined($opts{'generate'}) ) {
      $cfiles->copyfiles( "$res", $mode );
   }
   &cleanup();
}

# Check for crop resolutions
my $mode = "crop";
system( "../configure -s -crop on -bgc cn" );
my @crop_res = ( "10x15", "1.9x2.5" );
foreach my $res ( @crop_res ) {
   $options = "-res $res";
   eval{ system( "$bldnml $options  > $tempfile 2>&1 " ); };
   is( $@, '', "$options" );
   $cfiles->checkfilesexist( "$options", $mode );
   system( "diff lnd_in lnd_in.default.standard" );
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

# Check for glc_mec resolutions
my $mode = "standard";
system( "../configure -s" );
#my @glc_res = ( "48x96", "0.9x1.25", "1.9x2.5" ); # T31 NOT fully functional yet
my @glc_res = ( "0.9x1.25", "1.9x2.5" );
my @use_cases = ( "1850-2100_rcp2.6_glacierMEC_transient",
                  "1850-2100_rcp4.5_glacierMEC_transient",
                  "1850-2100_rcp6_glacierMEC_transient",
                  "1850-2100_rcp8.5_glacierMEC_transient",
                  "1850_glacierMEC_control",
                  "2000_glacierMEC_control",
                  "20thC_glacierMEC_transient",
                 );
my $GLC_NEC         = 10;
foreach my $res ( @glc_res ) {
   foreach my $usecase ( @usecases ) {
      $options = "-glc_nec $GLC_NEC -res $res -use_case $usecase ";
      eval{ system( "$bldnml $options > $tempfile 2>&1 " ); };
      is( $@, '', "$options" );
      $cfiles->checkfilesexist( "$options", $mode );
      system( "diff lnd_in lnd_in.default.standard" );
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

print "Successully ran all testing for build-namelist\n\n";

system( "/bin/rm $finidat" );

&cleanup( "config" );
system( "/bin/rm lnd_in.default" );
system( "/bin/rm $tempfile" );

sub cleanup {
#
# Cleanup files created
#
  my $type = shift;

  print "Cleanup files created\n";
  if ( defined($type) ) {
     if ( $type eq "config" ) {
        system( "/bin/rm Filepath config_cache.xml CESM_cppdefs" );
     }
  } else {
     system( "/bin/rm $tempfile *_in" );
  }
}
