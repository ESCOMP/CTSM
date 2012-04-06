#!/usr/bin/env perl

# Test methods of the build-namelist script.
# Just tests for configure mode using CESM.
# Try to test that all the different options at least work.
# Test that inconsistentcies are appropriately caught.

#########################

use Test::More tests => 44;

#########################

use strict;

system( "../configure -s" );

my $bldnml = "../build-namelist";

my $tempfile = "temp_file.txt";
if ( -f $tempfile ) {
  system( "/bin/rm $tempfile" );
}

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
eval{ system( "$bldnml" ); };
is( $@, '', "plain build-namelist" );
is( (-f "lnd_in"),       1, "lnd_in file exists" );
system( "/bin/cp  lnd_in lnd_in.default" );
system( "/bin/cat lnd_in.default" );
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
eval{ system( "$bldnml $options " ); };
is( $@, '', "options: $options" );
is( (-f "lnd_in"),       1, "lnd_in file exists" );
system( "diff lnd_in lnd_in.default" );
&cleanup();
# drydep and megan namelists
foreach my $options ( "-drydep", "-megan", "-drydep -megan" ) {
   eval{ system( "$bldnml $options " ); };
   is( $@, '', "options: $options" );
   is( (-f "lnd_in"),       1, "lnd_in file exists" );
   is( (-f "drv_flds_in"),  1, "flds file exists" );
   system( "diff lnd_in lnd_in.default" );
   system( "/bin/cat drv_flds_in" );
   &cleanup();
}
# irrig, verbose, clm_demand, rcp, test, sim_year, use_case
foreach my $options ( "-irrig", "-verbose", "-rcp 2.6", " -test", "-sim_year 1850",
                      "-use_case 1850_control" ) {
   eval{ system( "$bldnml $options " ); };
   is( $@, '', "options: $options" );
   is( (-f "lnd_in"),       1, "lnd_in file exists" );
   system( "diff lnd_in lnd_in.default" );
   &cleanup();
}
# Failure testing, do things that SHOULD fail
my $finidat  = "thing.nc";
system( "touch $finidat" );

my %failtest = ( 
     "coldstart but with IC file"=>{ options=>"-clm_start_type cold",
                                     namelst=>"finidat='$finidat'",
                                     expect=>"/^build-namelist ERROR::/" },
     "branch but NO nrevsn"      =>{ options=>"-clm_start_type branch",
                                     namelst=>"",
                                     expect=>"/^build-namelist ERROR::/" },
     "glc_nec inconsistent"      =>{ options=>"-glc_nec 10",
                                     namelst=>"maxpatch_glcmec=5",
                                     expect=>"/^build-namelist:: maxpatch_glcmec/" },
     "glc_smb inconsistent"      =>{ options=>"-glc_nec 10 -glc_smb .true.",
                                     namelst=>"glc_smb=.false.",
                                     expect=>"/^build-namelist:: glc_smb/" },
     "glc_grid inconsistent"     =>{ options=>"-glc_nec 10 -glc_grid gland10",
                                     namelst=>"glc_grid='gland5'",
                                     expect=>"/^build-namelist:: glc_grid/" },
     "rtm inconsistent"          =>{ options=>"-rtm off",
                                     namelst=>"do_rtm=.true.",
                                     expect=>"/^build-namelist ERROR:: when RTM/" },
     "rtm tstep inconsistent"    =>{ options=>"-rtm_tstep 10800",
                                     namelst=>"rtm_nsteps=1",
                                     expect=>"/^build-namelist:: rtm_nstep/" },
               );
foreach my $key ( keys(%failtest) ) {
   my $options  = $failtest{$key}{"options"};
   my $namelist = $failtest{$key}{"namelst"};
   eval{ system( "$bldnml $options -namelist \"&clmexp $namelist /\" > $tempfile 2>&1 " ); };
   my $result = `cat $tempfile`;
   like( $result, $failtest{$key}{"expect"}, $key );
}

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
