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
   my @hresols = ( "360x720","128x256","64x128","48x96","32x64","8x16","0.47x0.63","0.9x1.25",
                   "1.9x2.5","2.65x3.33","4x5","10x15","5x1_amazon", "1x1_tropicAtl", "1x1_camdenNJ","1x1_vancouverCAN",
                   "1x1_mexicocityMEX", "1x1_asphaltjungleNJ", "1x1_brazil", "1x1_urbanc_alpha" );
   my $nl = "namelist";
   my $sdate = "c" . `date +%y%m%d`;
   chomp( $sdate );

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
      #
      # Loop over each sim_year
      #
      #foreach my $sim_year ( 1992, 2000, 1890 ) {
      foreach my $sim_year ( 1992 ) {
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
 mksrf_furban       = '$CSMDATA/lnd/clm2/rawdata/mksrf_urban.060929.nc'
 mksrf_fdynuse      = ' '
EOF
         if ( $res =~ /[1-9]x[1-9]_[a-zA-Z0-9]/ ) {
            print $fh <<"EOF";
 mksrf_gridtype     = 'regional'
EOF
         }
         if ( $sim_year == 1992 ) {
            $desc = "";
            print $fh <<"EOF";
 mksrf_fvegtyp      = '$CSMDATA/lnd/clm2/rawdata/mksrf_pft.081008.nc'
 mksrf_fsoicol      = '$CSMDATA/lnd/clm2/rawdata/mksrf_soilcol.081008.nc'
 mksrf_flai         = '$CSMDATA/lnd/clm2/rawdata/mksrf_lai.081008.nc'
/
         elsif ( $sim_year == 2000 ) {
            $desc = "mcrop2000";
            print $fh <<"EOF";
 mksrf_fvegtyp      = '$CSMDATA/lnd/clm2/rawdata/mksrf_pft_mcrop2000.c081031.nc'
 mksrf_fsoicol      = '$CSMDATA/lnd/clm2/rawdata/mksrf_soilcol_mcrop2000.c081031.nc'
 mksrf_flai         = '$CSMDATA/lnd/clm2/rawdata/mksrf_lai_mcrop2000.c081031.nc'
/
EOF
         } else {
            $desc = "potveg";
            print $fh <<"EOF";
 mksrf_fvegtyp      = '$CSMDATA/lnd/clm2/rawdata/mksrf_pft_potveg.c081009.nc'
 mksrf_fsoicol      = '$CSMDATA/lnd/clm2/rawdata/mksrf_soilcol_potveg.c081009.nc'
 mksrf_flai         = '$CSMDATA/lnd/clm2/rawdata/mksrf_lai_potveg.c081009.nc'
/
EOF
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
         if ( -f "$ncfiles[0]" && -f "$lfiles[0]" ) {
            my $ofile = "surfdata_${res}_${desc}_${sdate}";
            my $mvcmd = "/bin/mv -f $ncfiles[0]  $CSMDATA/$surfdir/$ofile.nc";
            system( "$mvcmd" );
            my $mvcmd = "/bin/mv -f $lfiles[0] $CSMDATA/$surfdir/$ofile.log";
            system( "$mvcmd" );
            print $cfh "# FILE = \$DIN_LOC_ROOT/$surfdir/$ofile.nc\n";
            print $cfh "svn import -m $svnmesg \$CSMDATA/$surfdir/$ofile.nc $svnrepo/$surfdir/.\n";

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
