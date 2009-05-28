#!/usr/bin/env perl
#=======================================================================
#
#  This is a script to run the aerosol and nitrogen deposition regrid
#  scripts for all resolutions and simulation years.
#
# Usage:
#
# runDepositionRegrid.pl
#
#=======================================================================

use Cwd;
use strict;
#use diagnostics;
use English;

#-----------------------------------------------------------------------------------------------
my $nl_definition_file = "../../bld/namelist_files/namelist_definition.xml";
(-f "$nl_definition_file")  or  die <<"EOF";
** Cannot find namelist definition file \"$nl_definition_file\" **
EOF
print "Using namelist definition file $nl_definition_file\n";

my @dirs = ( "../../../../../scripts/ccsm_utils/Tools/perl5lib" );
unshift @INC, @dirs;

require Build::NamelistDefinition;

my $definition = Build::NamelistDefinition->new($nl_definition_file);

my @resolutions = $definition->get_valid_values( "res"      );
my @all_urb = ( "'1x1_camdenNJ'","'1x1_vancouverCAN'", "'1x1_mexicocityMEX'",
                "'1x1_asphaltjungleNJ'", "'1x1_urbanc_alpha'" );

my $runit = 1;

foreach my $res ( @resolutions ) {

   my @sim_years;
   my $all_urb = 0;
   foreach my $urb_res ( @all_urb ) {
     if ( $res eq $urb_res ) {
        $all_urb    = 1;
     }
   }
   if ( $all_urb == 1 ) {
      @sim_years = ( "2000" );
   } else {
      @sim_years = ( "1850", "2000" );
   }

   foreach my $sim_yr ( @sim_years ) {
      my $syscmd = "env RES=$res SIM_YR=$sim_yr ncl aerdepregrid.ncl";
      print( "$syscmd\n" );
      if ( $runit ) { system( $syscmd ); }
      my $syscmd = "env RES=$res SIM_YR=$sim_yr ncl ndepregrid.ncl";
      print( "$syscmd\n" );
      if ( $runit ) { system( $syscmd ); }
   }
}
my @trans_res = ( "0.9x1.25", "4x5" );

my $sim_yr = "1850-2000";

foreach my $res ( @trans_res ) {
   my $syscmd = "env RES=$res SIM_YR=$sim_yr ncl aerdepregrid.ncl";
   print( "$syscmd\n" );
   if ( $runit ) { system( $syscmd ); }
   my $syscmd = "env RES=$res SIM_YR=$sim_yr ncl ndepregrid.ncl";
   print( "$syscmd\n" );
   if ( $runit ) { system( $syscmd ); }
}
