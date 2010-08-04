#!/usr/bin/env perl
#
# Oct/30/2008                                         Erik Kluzek
#
# mksurfdata.pl Perl script to make surface datasets for all resolutions.
#
#
use Cwd;
use strict;
use English;
use IO::File;
use Getopt::Long;


#Figure out where configure directory is and where can use the XML/Lite module from
my $ProgName;
($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!; # name of program
my $ProgDir = $1;                         # name of directory where program lives

my $cwd = getcwd();  # current working directory
my $scrdir;

if ($ProgDir) { $scrdir = $ProgDir; }
else { $scrdir = $cwd; }

#-----------------------------------------------------------------------------------------------
# Add $scrdir to the list of paths that Perl searches for modules
my @dirs = ( $scrdir, "$scrdir/../../../../../scripts/ccsm_utils/Tools/perl5lib",
           );
unshift @INC, @dirs;
my $result = eval "require XML::Lite";
if ( ! defined($result) ) {
   die <<"EOF";
** Cannot find perl module \"XML/Lite.pm\" from directories: @dirs **
EOF
}
my $result = eval "require Build::NamelistDefinition";
if ( ! defined($result) ) {
   die <<"EOF";
** Cannot find perl module \"Build/NamelistDefinition.pm\" from directories: @dirs **
EOF
}
my $nldef_file     = "$scrdir/../../bld/namelist_files/namelist_definition.xml";

my $definition = Build::NamelistDefinition->new( $nldef_file );

my $CSMDATA = "/fs/cgd/csm/inputdata";
my $PFTDATA = "/cgd/tss";

my %opts = ( 
               hgrid=>"all", 
               rcp=>"-999.9", 
               debug=>0,
               years=>"1850,2000",
               help=>0,
               nomv=>undef,
               usrname=>"",
               csmdata=>$CSMDATA,
               pftdata=>$PFTDATA,
           );

#-----------------------------------------------------------------------------------------------
sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options]
OPTIONS
     -dinlc [or -l]                Enter the directory location for inputdata 
                                   (default $opts{'csmdata'})
     -debug [or -d]                Don't actually run -- just print out what 
                                   would happen if ran.
     -years [or -y]                Simulation year(s) to run over (by default $opts{'years'}) 
                                   (can also be a simulation year range: i.e. 1850-2000)
     -help  [or -h]                Display this help.
     -pftlc [or -p]                Enter directory location for pft data
                                   (default $opts{'pftdata'})
     -nomv                         Don't move the files to inputdata after completion.
     -res   [or -r] "resolution"   Resolution(s) to use for files (by default $opts{'hgrid'} ).
     -rcp   [or -c] "rep-con-path" Representative concentration pathway(s) to use for 
                                   future scenarios 
                                   (by default $opts{'rcp'}, where -999.9 means historical ).
     -usrname "clm_usrdat_name"    CLM user data name to find grid file with.

NOTE: years, res, and rcp can be comma delimited lists.

EOF
}

#-----------------------------------------------------------------------------------------------

   my $cmdline = "@ARGV";
   GetOptions(
        "r|res=s"      => \$opts{'hgrid'},
        "c|rcp=s"      => \$opts{'rcp'},
        "l|dinlc=s"    => \$opts{'csmdata'},
        "p|pftlc=s"    => \$opts{'pftdata'},
        "nomv"         => \$opts{'nomv'},
        "d|debug"      => \$opts{'debug'},
        "y|years=s"    => \$opts{'years'},
        "h|help"       => \$opts{'help'},
        "usrname=s"    => \$opts{'usrname'},
   ) or usage();

   # Check for unparsed arguments
   if (@ARGV) {
       print "ERROR: unrecognized arguments: @ARGV\n";
       usage();
   }
   if ( $opts{'help'} ) {
       usage();
   }
   # If csmdata was changed from the default
   if ( $CSMDATA ne $opts{'csmdata'} ) {
      $CSMDATA = $opts{'csmdata'};
   }
   my $pftdata = $opts{'pftdata'};
   #
   # Set disk location to send files to, and list resolutions to operate over, 
   # set filenames, and short-date-name
   #
   my @hresols;
   my @all_hresols = $definition->get_valid_values( "res" );
   if ( $opts{'hgrid'} eq "all" ) {
      @hresols = @all_hresols;
   } else {
      @hresols = split( ",", $opts{'hgrid'} );
      # Check that resolutions are valid
      foreach my $res ( @hresols ) {
         if ( ! $definition->is_valid_value( "res", "'$res'" ) ) {
            if ( $opts{'usrname'} eq ""  || $res ne $opts{'usrname'} ) {
               print "** Invalid resolution: $res\n";
               usage();
            }
         }
      }
   }
   #
   # Set years to run over
   #
   my @years   = split( ",", $opts{'years'} );
   # Check that resolutions are valid
   foreach my $sim_year ( @years ) {
      if ( ! $definition->is_valid_value( "sim_year", $sim_year ) ) {
         if ( ! $definition->is_valid_value( "sim_year_range", "'$sim_year'" ) ) {
            print "** Invalid simulation year or simulation year range: $sim_year\n";
            usage();
         }
      }
   }
   #
   # Set rcp to use
   #
   my @rcpaths = split( ",", $opts{'rcp'} );
   # Check that rcp is valid
   foreach my $rcp ( @rcpaths  ) {
      if ( ! $definition->is_valid_value( "rcp", $rcp ) ) {
         if ( ! $definition->is_valid_value( "rcp", "$rcp" ) ) {
            print "** Invalid rcp: $rcp\n";
            usage();
         }
      }
   }

   my $nl = "namelist";
   my $sdate = "c" . `date +%y%m%d`;
   chomp( $sdate );

   my @ncfiles;
   my @lfiles;
   my @pfiles;
   my $cfile = "clm.input_data_files";
   if ( -f "$cfile" ) {
      `/bin/mv -f $cfile ${cfile}.previous`;
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
   my $svnmesg = "Update fsurdat files with mksurfdata";
   my $surfdir = "lnd/clm2/surfdata";

   system( "/bin/rm -f surfdata_*.nc surfdata_*.log" );

   #
   # Loop over all resolutions listed
   #
   foreach my $res ( @hresols ) {
      #
      # Query the XML default file database to get the appropriate griddata file
      #
      my $queryopts = "-res $res -csmdata $CSMDATA -onlyfiles -silent -justvalue";
      my $usrnam    = "";
      if ( $opts{'usrname'} ne "" && $res eq $opts{'usrname'} ) {
         $usrnam    = "-usrname ".$opts{'usrname'};
      }
      my $griddata  = `$scrdir/../../bld/queryDefaultNamelist.pl $queryopts $usrnam -var fatmgrid`;
      if ( $? != 0 ) {
         die "ERROR:: fatmgrid file NOT found in XML database\n";
      }
      chomp( $griddata );
      if ( ! -f "$griddata" ) {
         die "ERROR:: fatmgrid file NOT found: $griddata\n";
      }
      print "res = $res griddata = $griddata\n";
      my $desc;
      my $desc_yr0;
      #
      # Check if all urban single point dataset
      #
      my @all_urb = ( "1x1_camdenNJ","1x1_vancouverCAN", "1x1_mexicocityMEX", 
                      "1x1_asphaltjungleNJ", "1x1_urbanc_alpha" );
      my $all_urb = ".false.";
      my $urb_pt  = 0;
      foreach my $urb_res ( @all_urb ) {
         if ( $res eq $urb_res ) {
            $all_urb = ".true.";
            $urb_pt  = 1;
         }
      }
      #
      # Always run at double precision for output
      #
      my $double = ".true.";
      #
      # Loop over each sim_year
      #
      RCP: foreach my $rcp ( @rcpaths ) {
         #
         # Loop over each sim_year
         #
         SIM_YEAR: foreach my $sim_year ( @years ) {
            #
            # Skip if urban unless sim_year=2000
            #
            if ( $urb_pt && $sim_year != 2000 ) {
               print "For urban -- skip this simulation year = $sim_year\n";
               next SIM_YEAR;
            }
            #
            # If year is 1850-2000 actually run 1850-2005
            #
            if ( $sim_year eq "1850-2000" ) {
               my $actual = "1850-2005";
               print "For $sim_year actually run $actual\n";
               $sim_year = $actual;
            }
            #
            # Create namelist file
            #
            my $fh = IO::File->new;
            $fh->open( ">$nl" ) or die "** can't open file: $nl\n";
            print $fh <<"EOF";
&clmexp
 mksrf_gridnm       = '$res'
 mksrf_fgrid        = '$griddata'
 mksrf_fsoitex      = '$CSMDATA/lnd/clm2/rawdata/mksrf_soitex.10level.c010119.nc'
 mksrf_forganic     = '$CSMDATA/lnd/clm2/rawdata/mksrf_organic.10level.0.5deg.081112.nc'
 mksrf_flanwat      = '$CSMDATA/lnd/clm2/rawdata/mksrf_lanwat.050425.nc'
 mksrf_fmax         = '$CSMDATA/lnd/clm2/rawdata/mksrf_fmax.070406.nc'
 mksrf_fglacier     = '$CSMDATA/lnd/clm2/rawdata/mksrf_glacier.060929.nc'
 mksrf_fvocef       = '$CSMDATA/lnd/clm2/rawdata/mksrf_vocef.c060502.nc'
 mksrf_ftopo        = '$CSMDATA/lnd/clm2/rawdata/mksrf_topo.10min.c080912.nc'
 mksrf_ffrac        = '$CSMDATA/lnd/clm2/griddata/fracdata_10min_USGS_071205.nc'
 outnc_double       = $double
 all_urban          = $all_urb
EOF
            my $urbdesc = "urb3den";
            if ( ! $urb_pt ) {
               print $fh <<"EOF";
 mksrf_furban       = '$CSMDATA/lnd/clm2/rawdata/mksrf_urban_3den_0.5x0.5_simyr2000.c090223_v1.nc'
EOF
            } else {
               #
               # Query the XML default file database to get the appropriate furbinp file
               #
               my $urbdata = `$scrdir/../../bld/queryDefaultNamelist.pl $queryopts -var fsurdat -filenameonly`;
               if ( $? != 0 ) {
                  die "ERROR:: furbinp file NOT found\n";
               }
               chomp( $urbdata );
               print $fh <<"EOF";
 mksrf_furban       = '$CSMDATA/lnd/clm2/surfdata/$urbdata'
EOF
            }
            my $resol = "";
            if ( $res =~ /[1-9]x[1-9](pt|)_[a-zA-Z0-9]/ ) {
               print $fh <<"EOF";
 mksrf_gridtype     = 'regional'
EOF
            }
            if ( $res ne "1x1_tropicAtl" ) {
               $resol = "-res 0.5x0.5";
            }
            my $sim_yr0 = $sim_year;
            my $sim_yrn = $sim_year;
            if ( $sim_year =~ /([0-9]+)-([0-9]+)/ ) {
               $sim_yr0 = $1;
               $sim_yrn = $2;
            }
            my $vegtyp = `$scrdir/../../bld/queryDefaultNamelist.pl $queryopts $resol -options sim_year=$sim_yr0 -var mksrf_fvegtyp -namelist clmexp`;
            chomp( $vegtyp );
            if ( $rcp == -999.9 ) {
               $desc     = sprintf( "hist_simyr%4.4d-%4.4d", $sim_yr0, $sim_yrn );
               $desc_yr0 = sprintf( "simyr%4.4d",            $sim_yr0  );
            } else {
               $desc     = sprintf( "%s%2.1f_simyr%4.4d-%4.4d", "rcp", $rcp, $sim_yr0, $sim_yrn );
               $desc_yr0 = sprintf( "%s%2.1f_simyr%4.4d",       "rcp", $rcp, $sim_yr0  );
            }
            my $pftdyntext_file = "pftdyn_$desc.txt";
            my $fhpftdyn = IO::File->new;
            $fhpftdyn->open( ">$pftdyntext_file" ) or die "** can't open file: $pftdyntext_file\n";
            print "Writing out pftdyn text file: $pftdyntext_file\n";
            for( my $yr = $sim_yr0; $yr <= $sim_yrn; $yr++ ) {
              my $vegtypyr = `$scrdir/../../bld/queryDefaultNamelist.pl $queryopts $resol -options sim_year=$yr,rcp=$rcp -var mksrf_fvegtyp -namelist clmexp`;
              chomp( $vegtypyr );
              $vegtypyr =~ s#^$PFTDATA#$pftdata#;
              printf $fhpftdyn "%-125.125s %4.4d\n", $vegtypyr, $yr;
              if ( $yr % 100 == 0 ) {
                 print "year: $yr\n";
              }
            }
            $fhpftdyn->close;
            print "Done writing file\n";

            print $fh <<"EOF";
 mksrf_fvegtyp      = '$vegtyp'
 mksrf_fsoicol      = '$pftdata/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/mksrf_soilcol_global_c090324.nc'
 mksrf_flai         = '$pftdata/pftlandusedyn.0.5x0.5.simyr1850-2005.c090630/mksrf_lai_global_c090506.nc'
 mksrf_fdynuse      = '$pftdyntext_file'
/
EOF
            $fh->close;
            print "resolution: $res rcp=$rcp sim_year = $sim_year\n";
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
            my $filehead;
            my $pfilehead;
            if ( ! $opts{'debug'} ) {
               system( "mksurfdata < $nl" );
               if ( $? ) { die "ERROR in mksurfdata: $?\n"; }
            } else {
               $filehead  = "surfdata_$res";
               $pfilehead = "surfdata.pftdyn_testfile";
               system( "touch $filehead.nc" );
               system( "touch $pfilehead.nc" );
               system( "touch $filehead.log" );
            }
            #
            # Check that files were created
            #
            @ncfiles  = glob( "surfdata_$res.nc" );
            if ( $#ncfiles != 0 ) {
              die "ERROR surfdata netcdf file was NOT created!\n";
            }
            chomp( $ncfiles[0] );
            @lfiles = glob( "surfdata_$res.log" );
            chomp( $lfiles[0] );
            @pfiles = glob( "surfdata.pftdyn_$res.nc" );
            chomp( $pfiles[0] );
            if ( $#pfiles != 0 ) {
              die "ERROR surfdata pftdyn netcdf file was NOT created!\n";
            }
            #
            # If urban point, append grid and frac file on top of surface dataset
            #
            if ( $urb_pt ) {
               my $cmd = "ncks -A $griddata $ncfiles[0]";
               print "$cmd\n";
               if ( ! $opts{'debug'} ) { system( $cmd ); }
               my $fracdata = `$scrdir/../../bld/queryDefaultNamelist.pl $queryopts -var fatmlndfrc`;
               if ( $? != 0 ) {
                  die "ERROR:: fatmlndfrc file NOT found\n";
               }
               chomp( $fracdata );
               $cmd = "ncks -A $fracdata $ncfiles[0]";
               print "$cmd\n";
               if ( ! $opts{'debug'} ) { system( $cmd ); }
            }
            #
            # Rename files to CSMDATA
            #
            my $lsvnmesg = "'$svnmesg $urbdesc $desc'";
            if ( -f "$ncfiles[0]" && -f "$lfiles[0]" ) {
               my $outdir = "$CSMDATA/$surfdir";
               if ( defined($opts{'nomv'}) ) {
                  $outdir = ".";
               }
               my $ofile = "surfdata_${res}_${desc_yr0}_${sdate}";
               my $mvcmd = "/bin/mv -f $ncfiles[0]  $outdir/$ofile.nc";
               print "$mvcmd\n";
               if ( ! $opts{'debug'} || ! defined($opts{'nomv'}) ) {
                  system( "$mvcmd" );
                  chmod( 0444, "$outdir/$ofile.nc" );
               }
               my $mvcmd = "/bin/mv -f $lfiles[0] $outdir/$ofile.log";
               print "$mvcmd\n";
               if ( ! $opts{'debug'} || ! defined($opts{'nomv'}) ) {
                  system( "$mvcmd" );
                  chmod( 0444, "$outdir/$ofile.log" );
               }
               if ( ! defined($opts{'nomv'}) ) {
                  print $cfh "# FILE = \$DIN_LOC_ROOT/$surfdir/$ofile.nc\n";
                  print $cfh "svn import -m $lsvnmesg \$CSMDATA/$surfdir/$ofile.nc " . 
                             "$svnrepo/$surfdir/$ofile.nc\n";
                  print $cfh "# FILE = \$DIN_LOC_ROOT/$surfdir/$ofile.log\n";
                  print $cfh "svn import -m $lsvnmesg \$CSMDATA/$surfdir/$ofile.log " .
                             "$svnrepo/$surfdir/$ofile.log\n";
               }
               # If running a transient case
               if ( $sim_year ne $sim_yr0 ) {
                  $ofile = "surfdata.pftdyn_${res}_${desc}_${sdate}";
                  $mvcmd = "/bin/mv -f $pfiles[0] $outdir/$ofile.nc";
                  print "$mvcmd\n";
                  if ( ! $opts{'debug'} || ! defined($opts{'nomv'}) ) {
                     system( "$mvcmd" );
                     chmod( 0444, "$outdir/$ofile.nc" );
                  }
                  if ( ! defined($opts{'nomv'}) ) {
                     print $cfh "# FILE = \$DIN_LOC_ROOT/$surfdir/$ofile.nc\n";
                     print $cfh "svn import -m $lsvnmesg \$CSMDATA/$surfdir/$ofile.nc " .
                                "$svnrepo/$surfdir/$ofile.nc\n";
                  }
               }
   
            } else {
              die "ERROR files were NOT created: nc=$ncfiles[0] log=$lfiles[0]\n";
            }
            if ( (! $opts{'debug'}) && (-f "$ncfiles[0]" || -f "$lfiles[0]") ) {
              die "ERROR files were NOT moved: nc=$ncfiles[0] log=$lfiles[0]\n";
            }
            if ( ! $opts{'debug'} ) {
               system( "/bin/rm -f $filehead.nc $filehead.log $pfilehead.nc" );
            }
         } # End of sim_year loop
      }    # End of rcp loop
   }
   close( $cfh );
   print "Successfully created fsurdat files\n";
