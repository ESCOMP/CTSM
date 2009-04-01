#!/usr/bin/env perl
#
# Oct/30/2008                                         Erik Kluzek
#
# mksurfdata.pl Perl script to make surface datasets for all resolutions.
#
#
use strict;
use IO::File;

   #
   # Set disk location to send files to, and list resolutions to operate over, set filenames, and short-date-name
   #
   my $CSMDATA = "/fs/cgd/csm/inputdata";
   my @hresols = ( "360x720","128x256", "0.47x0.63" );
   my @hresols = ( "64x128","48x96","32x64","8x16","0.9x1.25",
                   "1.9x2.5","2.65x3.33","4x5","10x15","5x5_amazon", "1x1_tropicAtl", "1x1_camdenNJ","1x1_vancouverCAN",
                   "1x1_mexicocityMEX", "1x1_asphaltjungleNJ", "1x1_brazil", "1x1_urbanc_alpha" );
   my $nl = "namelist";
   my $sdate = "c" . `date +%y%m%d`;
   chomp( $sdate );

   my $urban = 1;
   my @ncfiles;
   my @lfiles;
   my $cfile = "clm.input_data_files";
   if ( -f "$cfile" ) {
      `mv $cfile ${cfile}.previous`;
   }
   my $cfh = IO::File->new;
   $cfh->open( ">$cfile" ) or die "** can't open file: $cfile\n";
   system( "\rm -f $cfile" );
   system( "touch $cfile" );
   print $cfh <<"EOF";
#! /bin/csh -f
set CSMDATA = $CSMDATA
EOF
   system( "chmod +x $cfile" );
   my $svnrepo = "https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata";
   my $svnmesg = "'Update fsurdat files with mksurfdata'";
   my $surfdir = "lnd/clm2/surfdata";

   system( "/bin/rm surfdata_*.nc surfdata_*.log" );

   #
   # Loop over all resolutions listed
   #
   foreach my $res ( @hresols ) {
      #
      # Query the XML default file database to get the appropriate griddata file
      #
      my $griddata = `../../bld/queryDefaultNamelist.pl -res $res -csmdata $CSMDATA -onlyfiles -silent -justvalue -demand -var fatmgrid`;
      if ( $? != 0 ) {
         die "ERROR:: fatmgrid file NOT found\n";
      }
      chomp( $griddata );
      print "res = $res griddata = $griddata\n";
      my $desc;
      my @all_urb = ( "1x1_camdenNJ","1x1_vancouverCAN", "1x1_mexicocityMEX", 
                      "1x1_asphaltjungleNJ", "1x1_urbanc_alpha" );
      my $all_urb = "  all_urban = .false.";
      #
      # Loop over each sim_year
      #
      foreach my $sim_year ( 1850, 2000 ) {
         #
         # Create namelist file
         #
         my $fh = IO::File->new;
         $fh->open( ">$nl" ) or die "** can't open file: $nl\n";
         print $fh <<"EOF";
&clmexp
 mksrf_fgrid        = '$griddata'
 mksrf_fsoitex      = '$CSMDATA/lnd/clm2/rawdata/mksrf_soitex.10level.c010119.nc'
 mksrf_forganic     = '$CSMDATA/lnd/clm2/rawdata/mksrf_organic.10level.0.5deg.081112.nc'
 mksrf_flanwat      = '$CSMDATA/lnd/clm2/rawdata/mksrf_lanwat.050425.nc'
 mksrf_fmax         = '$CSMDATA/lnd/clm2/rawdata/mksrf_fmax.070406.nc'
 mksrf_fglacier     = '$CSMDATA/lnd/clm2/rawdata/mksrf_glacier.060929.nc'
 mksrf_ftopo        = '$CSMDATA/lnd/clm2/rawdata/mksrf_topo.10min.c080912.nc'
 mksrf_ffrac        = '$CSMDATA/lnd/clm2/griddata/fracdata_10min_USGS_071205.nc'
 mksrf_fdynuse      = ' '
 $all_urb
EOF
         my $urbdesc;
         if ( $urban ) {
            $urbdesc = "urb3den";
            print $fh <<"EOF";
 mksrf_furban       = '$CSMDATA/lnd/clm2/rawdata/mksrf_urban_3den_0.5x0.5_simyr2000.c090223_v1.nc'
EOF
         } else {
            $urbdesc = "nourb";
            print $fh <<"EOF";
 mksrf_furban       = '$CSMDATA/lnd/clm2/rawdata/mksrf_urban.060929.nc'
EOF
         }
         if ( $res =~ /[1-9]x[1-9]_[a-zA-Z0-9]/ ) {
            print $fh <<"EOF";
 mksrf_gridtype     = 'regional'
EOF
         }
         $desc = "simyr$sim_year";
         my $sdate = "c090313";
         if (      $sim_year == 2000 ) {
            $sdate = "c090320";
         } elsif ( $sim_year == 1850 ) {
            $sdate = "c090220";
         }
         print $fh <<"EOF";
 mksrf_fvegtyp      = '$CSMDATA/lnd/clm2/rawdata/mksrf_pft_0.5x0.5_$desc.$sdate.nc'
 mksrf_fsoicol      = '$CSMDATA/lnd/clm2/rawdata/mksrf_soilcol_0.5x0.5_$desc.$sdate.nc'
 mksrf_flai         = '$CSMDATA/lnd/clm2/rawdata/mksrf_lai_0.5x0.5_$desc.$sdate.nc'
/
EOF
         if ( $sim_year != 2005 && $sim_year != 2000 && $sim_year != 1990 && 
              $sim_year != 1870 && $sim_year != 1850 ) {
            die "Bad sim_year = $sim_year, expecting: 1850, 1870, 1990, 2000, or 2005\n";
         }
         $fh->close;
         print "namelist: $nl\n";
         $fh->open( "<$nl" ) or die "** can't open file: $nl\n";
         while( $_ = <$fh> ) {
           print $_;
         }
         $fh->close;
         #
         # Run mksurfdata with the namelist file
         #
         print "mksurfdata < $nl\n";
         system( "mksurfdata < $nl" );
         if ( $? ) { die "ERROR in mksurfdata: $?\n"; }

         #
         # Check that files were created and rename them to CSMDATA
         #
         @ncfiles  = glob( "surfdata_*.nc" );
         if ( $#ncfiles != 0 ) {
           die "ERROR surfdata netcdf file was NOT created!\n";
         }
         chomp( $ncfiles[0] );
         @lfiles = glob( "surfdata_*.log" );
         chomp( $lfiles[0] );
         my $lsvnmesg = "$svnmesg $urbdesc $desc";
         if ( -f "$ncfiles[0]" && -f "$lfiles[0]" ) {
            my $ofile = "surfdata_${res}_${urbdesc}_${desc}_${sdate}";
            my $mvcmd = "/bin/mv -f $ncfiles[0]  $CSMDATA/$surfdir/$ofile.nc";
            print "$mvcmd\n";
            system( "$mvcmd" );
            my $mvcmd = "/bin/mv -f $lfiles[0] $CSMDATA/$surfdir/$ofile.log";
            print "$mvcmd\n";
            system( "$mvcmd" );
            print $cfh "# FILE = \$DIN_LOC_ROOT/$surfdir/$ofile.nc\n";
            print $cfh "svn import -m $lsvnmesg \$CSMDATA/$surfdir/$ofile.nc $svnrepo/$surfdir/$ofile.nc\n";
            print $cfh "# FILE = \$DIN_LOC_ROOT/$surfdir/$ofile.log\n";
            print $cfh "svn import -m $lsvnmesg \$CSMDATA/$surfdir/$ofile.log $svnrepo/$surfdir/$ofile.log\n";

         } else {
           die "ERROR files were NOT created: nc=$ncfiles[0] log=$lfiles[0]\n";
         }
         if ( -f "$ncfiles[0]" || -f "$lfiles[0]" ) {
           die "ERROR files were NOT moved: nc=$ncfiles[0] log=$lfiles[0]\n";
         }
      }
   }
   close( $cfh );
   print "Successfully created fsurdat files\n";
