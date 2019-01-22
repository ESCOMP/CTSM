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

my $debug = 0;

#-----------------------------------------------------------------------------------------------
# Add $scrdir to the list of paths that Perl searches for modules
my @dirs = ( "$scrdir/../../cime/utils/perl5lib",
             "$scrdir/../../../../cime/utils/perl5lib"
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
my $nldef_file     = "$scrdir/../../bld/namelist_files/namelist_definition_ctsm.xml";

my $definition = Build::NamelistDefinition->new( $nldef_file );

my $CSMDATA = "/glade/p/cesm/cseg/inputdata";

my %opts = ( 
               hgrid=>"all", 
               ssp_rcp=>"hist", 
               debug=>0,
               exedir=>undef,
               allownofile=>undef,
               crop=>1,
               fast_maps=>0,
               hirespft=>undef,
               years=>"1850,2000",
               glc_nec=>10,
               merge_gis=>undef,
               inlandwet=>undef,
               help=>0,
               no_surfdata=>0,
               pft_override=>undef,
               pft_frc=>undef,
               pft_idx=>undef,
               soil_override=>undef,
               soil_cly=>undef,
               soil_snd=>undef,
               soil_col=>undef,
               soil_fmx=>undef,
	       outnc_double=>undef,
	       outnc_dims=>"2",     
               usrname=>"",
               rundir=>"$cwd",
               usr_mapdir=>"../mkmapdata",
               dynpft=>undef,
               csmdata=>$CSMDATA,
               urban_skip_abort_on_invalid_data_check=>undef,
           );

my $numpft = 78;

#-----------------------------------------------------------------------------------------------
sub usage {
    die <<EOF;
SYNOPSIS 

     For supported resolutions:
     $ProgName -res <res>  [OPTIONS]
        -res [or -r] "resolution" is the supported resolution(s) to use for files (by default $opts{'hgrid'} ).

      
     For unsupported, user-specified resolutions:	
     $ProgName -res usrspec -usr_gname <user_gname> -usr_gdate <user_gdate>  [OPTIONS]
        -usr_gname "user_gname"    User resolution name to find grid file with 
                                   (only used if -res is set to 'usrspec')
        -usr_gdate "user_gdate"    User map date to find mapping files with
                                   (only used if -res is set to 'usrspec')
                                   NOTE: all mapping files are assumed to be in mkmapdata
                                    - and the user needs to have invoked mkmapdata in 
                                      that directory first
        -usr_mapdir "mapdirectory" Directory where the user-supplied mapping files are
                                   Default: $opts{'usr_mapdir'}

OPTIONS
     -allownofile                  Allow the script to run even if one of the input files
                                   does NOT exist.
     -dinlc [or -l]                Enter the directory location for inputdata 
                                   (default $opts{'csmdata'})
     -debug [or -d]                Do not actually run -- just print out what 
                                   would happen if ran.
     -dynpft "filename"            Dynamic PFT/harvesting file to use 
                                   (rather than create it on the fly) 
                                   (must be consistent with first year)
     -fast_maps                    Toggle fast mode which doesn't use the large mapping files
     -glc_nec "number"             Number of glacier elevation classes to use (by default $opts{'glc_nec'})
     -merge_gis                    If you want to use the glacier dataset that merges in
                                   the Greenland Ice Sheet data that CISM uses (typically
                                   used only if consistency with CISM is important)
     -hirespft                     If you want to use the high-resolution pft dataset rather 
                                   than the default lower resolution dataset
                                   (low resolution is at half-degree, high resolution at 3minute)
                                   (hires only available for present-day [2000])
     -exedir "directory"           Directory where mksurfdata_map program is
                                   (by default assume it is in the current directory)
     -inlandwet                    If you want to allow inland wetlands
     -no-crop                      Create datasets without the extensive list of prognostic crop types
     -no_surfdata                  Do not output a surface dataset
                                   This is useful if you only want a landuse_timeseries file
     -years [or -y] "years"        Simulation year(s) to run over (by default $opts{'years'}) 
                                   (can also be a simulation year range: i.e. 1850-2000)
     -help  [or -h]                Display this help.

     -rundir "directory"           Directory to run in
                                   (by default current directory $opts{'rundir'})

     -ssp_rcp "scenario-name"      Shared Socioeconomic Pathway and Representative Concentration Pathway Scenario name(s).
                                   "hist" for historical, otherwise in form of SSPn-m.m where n is the SSP number
                                   and m.m is the radiative forcing in W/m^2 at the peak or 2100.
     
     -usrname "clm_usrdat_name"    CLM user data name to find grid file with.

      NOTE: years, res, and ssp_rcp can be comma delimited lists.


OPTIONS to override the mapping of the input gridded data with hardcoded input

     -pft_frc "list of fractions"  Comma delimited list of percentages for veg types
     -pft_idx "list of veg index"  Comma delimited veg index for each fraction
     -soil_cly "% of clay"         % of soil that is clay
     -soil_col "soil color"        Soil color (1 [light] to 20 [dark])
     -soil_fmx "soil fmax"         Soil maximum saturated fraction (0-1)
     -soil_snd "% of sand"         % of soil that is sand

OPTIONS to work around bugs?
     -urban_skip_abort_on_invalid_data_check
                                   do not abort on an invalid data check in urban.
                                   Added 2015-01 to avoid recompiling as noted in
                                   /glade/p/cesm/cseg/inputdata/lnd/clm2/surfdata_map/README_c141219

EOF
}

sub check_soil {
#
# check that the soil options are set correctly
#
  foreach my $type ( "soil_cly", "soil_snd" ) {
     if ( ! defined($opts{$type} )  ) {
        die "ERROR: Soil variables were set, but $type was NOT set\n";
     }
  }
  #if ( $opts{'soil_col'} < 0 || $opts{'soil_col'} > 20 ) {
  #   die "ERROR: Soil color is out of range = ".$opts{'soil_col'}."\n";
  #}
  my $texsum = $opts{'soil_cly'} + $opts{'soil_snd'};
  my $loam   = 100.0 - $texsum;
  if ( $texsum < 0.0 || $texsum > 100.0 ) {
     die "ERROR: Soil textures are out of range: clay = ".$opts{'soil_cly'}.
         " sand = ".$opts{'soil_snd'}." loam = $loam\n";
  }
}

sub check_soil_col_fmx {
#
# check that the soil color or soil fmax option is set correctly
#
  if ( defined($opts{'soil_col'}) ) {
     if ( $opts{'soil_col'} < 0 || $opts{'soil_col'} > 20 ) {
        die "ERROR: Soil color is out of range = ".$opts{'soil_col'}."\n";
     }
  }
  if ( defined($opts{'soil_fmx'}) ) {
     if ( $opts{'soil_fmx'} < 0.0 || $opts{'soil_fmx'} > 1.0 ) {
        die "ERROR: Soil fmax is out of range = ".$opts{'soil_fmx'}."\n";
     }
  }
}

sub check_pft {
#
# check that the pft options are set correctly
#
  # Eliminate starting and ending square brackets
  $opts{'pft_idx'} =~ s/^\[//;
  $opts{'pft_idx'} =~ s/\]$//;
  $opts{'pft_frc'} =~ s/^\[//;
  $opts{'pft_frc'} =~ s/\]$//;
  foreach my $type ( "pft_idx", "pft_frc" ) {
     if ( ! defined($opts{$type} ) ) {
        die "ERROR: PFT variables were set, but $type was NOT set\n";
     }
  }
  my @pft_idx     = split( /,/, $opts{'pft_idx'} );
  my @pft_frc     = split( /,/, $opts{'pft_frc'} );
  if ( $#pft_idx != $#pft_frc ) {
     die "ERROR: PFT arrays are different sizes: pft_idx and pft_frc\n";
  }
  my $sumfrc = 0.0;
  for( my $i = 0; $i <= $#pft_idx; $i++ ) {
     # check index in range
     if ( $pft_idx[$i] < 0 || $pft_idx[$i] > $numpft ) {
         die "ERROR: pft_idx out of range = ".$opts{'pft_idx'}."\n";
     }
     # make sure there are no duplicates
     for( my $j = 0; $j < $i; $j++ ) {
        if ( $pft_idx[$i] == $pft_idx[$j] ) {
            die "ERROR: pft_idx has duplicates = ".$opts{'pft_idx'}."\n";
        }
     }
     # check fraction in range
     if ( $pft_frc[$i] <= 0.0 || $pft_frc[$i] > 100.0 ) {
         die "ERROR: pft_frc out of range (>0.0 and <=100.0) = ".$opts{'pft_frc'}."\n";
     }
     $sumfrc = $sumfrc + $pft_frc[$i];
  }
  # check that fraction sums up to 100%
  if ( abs( $sumfrc - 100.0) > 1.e-6 ) {
      die "ERROR: pft_frc does NOT add up to 100% = ".$opts{'pft_frc'}."\n";
  }
  
}
 
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

sub write_transient_timeseries_file {
  my ($transient, $desc, $sim_yr0, $sim_yrn, $queryfilopts, $resol, $resolhrv, $ssp_rcp, $mkcrop, $sim_yr_surfdat) = @_;

  my $strlen = 195;
  my $dynpft_format = "%-${strlen}.${strlen}s %4.4d\n";
  my $landuse_timeseries_text_file = "";
  if ( $transient ) {
    if ( ! defined($opts{'dynpft'}) && ! $opts{'pft_override'} ) {
      $landuse_timeseries_text_file = "landuse_timeseries_$desc.txt";
      my $fh_landuse_timeseries = IO::File->new;
      $fh_landuse_timeseries->open( ">$landuse_timeseries_text_file" ) or die "** can't open file: $landuse_timeseries_text_file\n";
      print "Writing out landuse_timeseries text file: $landuse_timeseries_text_file\n";
      for( my $yr = $sim_yr0; $yr <= $sim_yrn; $yr++ ) {
        my $vegtypyr = `$scrdir/../../bld/queryDefaultNamelist.pl $queryfilopts $resol -options sim_year=$yr,ssp-rcp=${ssp_rcp}${mkcrop} -var mksrf_fvegtyp -namelist clmexp`;
        chomp( $vegtypyr );
        printf $fh_landuse_timeseries $dynpft_format, $vegtypyr, $yr;
        my $hrvtypyr = `$scrdir/../../bld/queryDefaultNamelist.pl $queryfilopts $resolhrv -options sim_year=$yr,ssp-rcp=${ssp_rcp}${mkcrop} -var mksrf_fvegtyp -namelist clmexp`;
        chomp( $hrvtypyr );
        printf $fh_landuse_timeseries $dynpft_format, $hrvtypyr, $yr;
        if ( $yr % 100 == 0 ) {
          print "year: $yr\n";
        }
      }
      $fh_landuse_timeseries->close;
      print "Done writing file\n";
    } elsif ( $opts{'pft_override'} && defined($opts{'dynpft'}) ) {
      $landuse_timeseries_text_file = $opts{'dynpft'};
    } else {
      $landuse_timeseries_text_file = "landuse_timeseries_override_$desc.txt";
      my $fh_landuse_timeseries = IO::File->new;
      $fh_landuse_timeseries->open( ">$landuse_timeseries_text_file" ) or die "** can't open file: $landuse_timeseries_text_file\n";
      my $frstpft = "<pft_f>$opts{'pft_frc'}</pft_f>" . 
        "<pft_i>$opts{'pft_idx'}</pft_i>" .
        "<harv>0,0,0,0,0</harv><graz>0</graz>";
      print "Writing out landuse_timeseries text file: $landuse_timeseries_text_file\n";
      if ( (my $len = length($frstpft)) > $strlen ) {
        die "ERROR PFT line is too long ($len): $frstpft\n";
      }
      # NOTE(wjs, 2014-12-04) Using sim_yr_surfdat here rather than
      # sim_yr0. As far as I can tell, it seems somewhat arbitrary which one
      # we use, but sim_yr_surfdat seems more like what's intended.
      printf $fh_landuse_timeseries $dynpft_format, $frstpft, $sim_yr_surfdat;
      $fh_landuse_timeseries->close;
      print "Done writing file\n";
    }
  }
  return $landuse_timeseries_text_file;
}

sub write_namelist_file {
   my ($namelist_fname, $logfile_fname, $fsurdat_fname, $fdyndat_fname,
      $glc_nec, $griddata, $map, $datfil, $double,
      $all_urb, $no_inlandwet, $vegtyp, $hrvtyp, 
      $landuse_timeseries_text_file, $setnumpft) = @_;


  my $gitdescribe = `cd $scrdir; git describe; cd -`;
  chomp( $gitdescribe );
  my $fh = IO::File->new;
  $fh->open( ">$namelist_fname" ) or die "** can't open file: $namelist_fname\n";
  print $fh <<"EOF";
&clmexp
 nglcec           = $glc_nec
 mksrf_fgrid      = '$griddata'
 map_fpft         = '$map->{'veg'}'
 map_fglacier     = '$map->{'glc'}'
 map_fglacierregion = '$map->{'glcregion'}'
 map_fsoicol      = '$map->{'col'}'
 map_furban       = '$map->{'urb'}'
 map_fmax         = '$map->{'fmx'}'
 map_forganic     = '$map->{'org'}'
 map_flai         = '$map->{'lai'}'
 map_fharvest     = '$map->{'hrv'}'
 map_flakwat      = '$map->{'lak'}'
 map_fwetlnd      = '$map->{'wet'}'
 map_fvocef       = '$map->{'voc'}'
 map_fsoitex      = '$map->{'tex'}'
 map_furbtopo     = '$map->{'utp'}'
 map_fgdp         = '$map->{'gdp'}'
 map_fpeat        = '$map->{'peat'}'
 map_fsoildepth   = '$map->{'soildepth'}'
 map_fabm         = '$map->{'abm'}'
 map_fvic         = '$map->{'vic'}'
 map_fch4         = '$map->{'ch4'}'
 mksrf_fsoitex    = '$datfil->{'tex'}'
 mksrf_forganic   = '$datfil->{'org'}'
 mksrf_flakwat    = '$datfil->{'lak'}'
 mksrf_fwetlnd    = '$datfil->{'wet'}'
 mksrf_fmax       = '$datfil->{'fmx'}'
 mksrf_fglacier   = '$datfil->{'glc'}'
 mksrf_fglacierregion = '$datfil->{'glcregion'}'
 mksrf_fvocef     = '$datfil->{'voc'}'
 mksrf_furbtopo   = '$datfil->{'utp'}'
 mksrf_fgdp       = '$datfil->{'gdp'}'
 mksrf_fpeat      = '$datfil->{'peat'}'
 mksrf_fsoildepth = '$datfil->{'soildepth'}'
 mksrf_fabm       = '$datfil->{'abm'}'
 mksrf_fvic       = '$datfil->{'vic'}'
 mksrf_fch4       = '$datfil->{'ch4'}'
 outnc_double   = $double
 all_urban      = $all_urb
 no_inlandwet   = $no_inlandwet
 mksrf_furban   = '$datfil->{'urb'}'
 gitdescribe    = '$gitdescribe'
EOF
  if ( ! $opts{'fast_maps'} ) {
    print $fh <<"EOF";
 map_ftopostats   = '$map->{'topostats'}'
 mksrf_ftopostats = '$datfil->{'topostats'}'
EOF
  } else {
    print $fh <<"EOF";
 std_elev         = 371.0d00
EOF
  }
  if ( defined($opts{'soil_override'}) ) {
    print $fh <<"EOF";
 soil_clay     = $opts{'soil_cly'}
 soil_sand     = $opts{'soil_snd'}
EOF
  }
  if ( defined($opts{'pft_override'}) ) {
    print $fh <<"EOF";
 pft_frc      = $opts{'pft_frc'}
 pft_idx      = $opts{'pft_idx'}
EOF
  }

  print $fh <<"EOF";
 mksrf_fvegtyp  = '$vegtyp'
 mksrf_fhrvtyp  = '$hrvtyp'
 mksrf_fsoicol  = '$datfil->{'col'}'
 mksrf_flai     = '$datfil->{'lai'}'
EOF

  # Note that some of the file names in the following may be empty strings
  # (except for logfile_fname)
  print $fh <<"EOF";
 fsurdat        = '$fsurdat_fname'
 fsurlog        = '$logfile_fname'
 mksrf_fdynuse  = '$landuse_timeseries_text_file'
 fdyndat        = '$fdyndat_fname'
EOF

  if ( $setnumpft ) {
    print $fh <<"EOF";
 $setnumpft
EOF
  }

  if ( $opts{'urban_skip_abort_on_invalid_data_check'} ) {
    print $fh <<"EOF";
 urban_skip_abort_on_invalid_data_check = .true.
EOF
  }
  # end the namelist
  print $fh <<"EOF";
/
EOF

  $fh->close;
  # 
  # Print namelist file 
  $fh->open( "<$namelist_fname" ) or die "** can't open file: $namelist_fname\n";
  while( $_ = <$fh> ) {
    print $_;
  }
  $fh->close;
}

#-----------------------------------------------------------------------------------------------

   my $cmdline = "@ARGV";
   GetOptions(
        "allownofile"  => \$opts{'allownofile'},
        "r|res=s"      => \$opts{'hgrid'},
        "usr_gname=s"  => \$opts{'usr_gname'},
        "usr_gdate=s"  => \$opts{'usr_gdate'},
        "usr_mapdir=s" => \$opts{'usr_mapdir'},
        "crop!"        => \$opts{'crop'},
        "hirespft"     => \$opts{'hirespft'},
        "l|dinlc=s"    => \$opts{'csmdata'},
        "d|debug"      => \$opts{'debug'},
        "fast_maps"    => \$opts{'fast_maps'},
        "dynpft=s"     => \$opts{'dynpft'},
        "y|years=s"    => \$opts{'years'},
        "exedir=s"     => \$opts{'exedir'},
        "h|help"       => \$opts{'help'},
        "usrname=s"    => \$opts{'usrname'},
        "glc_nec=i"    => \$opts{'glc_nec'},
        "merge_gis"    => \$opts{'merge_gis'},
        "inlandwet"    => \$opts{'inlandwet'},
        "no_surfdata"  => \$opts{'no_surfdata'},
        "pft_frc=s"    => \$opts{'pft_frc'},
        "pft_idx=s"    => \$opts{'pft_idx'},
        "ssp_rcp=s"    => \$opts{'ssp_rcp'},
        "rundir=s"     => \$opts{'rundir'},
        "soil_col=i"   => \$opts{'soil_col'},
        "soil_fmx=f"   => \$opts{'soil_fmx'},
        "soil_cly=f"   => \$opts{'soil_cly'},
        "soil_snd=f"   => \$opts{'soil_snd'},
        "urban_skip_abort_on_invalid_data_check" => \$opts{'urban_skip_abort_on_invalid_data_check'},
     ) or usage();

   # Check for unparsed arguments
   if (@ARGV) {
       print "ERROR: unrecognized arguments: @ARGV\n";
       usage();
   }
   if ( $opts{'help'} ) {
       usage();
   }

   chdir( $opts{'rundir'} ) or die "** can't change to directory: $opts{'rundir'}\n";
   # If csmdata was changed from the default
   if ( $CSMDATA ne $opts{'csmdata'} ) {
      $CSMDATA = $opts{'csmdata'};
   }
   my $glc_nec = $opts{'glc_nec'};
   if ( $glc_nec <= 0 ) {
      print "** glc_nec must be at least 1\n";
      usage();
   }
   my $no_inlandwet = ".true.";
   if (defined($opts{'inlandwet'})) {
      $no_inlandwet = ".false.";
   }
   #
   # Set disk location to send files to, and list resolutions to operate over, 
   # set filenames, and short-date-name
   #
   my @hresols;
   my $mapdate; 
   if ( $opts{'hgrid'} eq "all" ) {
      my @all_hresols = $definition->get_valid_values( "res" );
      @hresols = @all_hresols;
   } elsif ( $opts{'hgrid'} eq "usrspec" ) {
      @hresols = $opts{'usr_gname'}; 
      $mapdate = $opts{'usr_gdate'}; 
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
     if ("-" eq substr($sim_year, 4, 1)) {
       # range of years for transient run
       if ( ! $definition->is_valid_value( "sim_year_range", "'$sim_year'" ) ) {
         print "** Invalid simulation simulation year range: $sim_year\n";
         usage();
       }
     } else {
       # single year.
       if ( ! $definition->is_valid_value( "sim_year", $sim_year ) ) {
         print "** Invalid simulation year: $sim_year\n";
         usage();
       }
     }
   }
   #
   # Set ssp-rcp to use
   #
   my @rcpaths = split( ",", $opts{'ssp_rcp'} );
   # Check that ssp-rcp is valid
   foreach my $ssp_rcp ( @rcpaths  ) {
      if ( ! $definition->is_valid_value( "ssp-rcp", $ssp_rcp ) ) {
          print "** Invalid ssp_rcp: $ssp_rcp\n";
          usage();
       }
   }

   # CMIP series input data is corresponding to
   my $cmip_series = "CMIP6";
   # Check if soil set
   if ( defined($opts{'soil_cly'}) || 
        defined($opts{'soil_snd'}) ) {
       &check_soil( );
       $opts{'soil_override'} = 1;
   }
   # Check if pft set
   if ( ! $opts{'crop'} ) { $numpft = 16; }   # First set numpft if crop is off
   if ( defined($opts{'pft_frc'}) || defined($opts{'pft_idx'}) ) {
       &check_pft( );
       $opts{'pft_override'} = 1;
   }
   # Check if dynpft set and is valid filename
   if ( defined($opts{'dynpft'}) ) {
       if ( ! -f $opts{'dynpft'} ) {
          print "** Dynamic PFT file does NOT exist: $opts{'dynpft'}\n";
          usage();
       }
   }

   my $sdate = "c" . `date +%y%m%d`;
   chomp( $sdate );

   my $cfile = "clm.input_data_list";
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
   my $surfdir = "lnd/clm2/surfdata";

   # string to add to options for crop off or on
   my $mkcrop_off = ",crop='on'";
   my $mkcrop_on  = ",crop='on'";

   #
   # Loop over all resolutions listed
   #
   foreach my $res ( @hresols ) {
      #
      # Query the XML default file database to get the appropriate files
      #
      my $queryopts, my $queryfilopts; 
      if ( $opts{'hgrid'} eq "usrspec" ) {
	  $queryopts = "-csmdata $CSMDATA -silent -justvalue";
      } else {
	  $queryopts = "-res $res -csmdata $CSMDATA -silent -justvalue";
      }
      $queryfilopts = "$queryopts -onlyfiles ";
      my $mkcrop = $mkcrop_off;
      my $setnumpft = "";
      $mkcrop    = $mkcrop_on;
      $setnumpft = "numpft = $numpft";
      my $usrnam    = "";
      if ( $opts{'usrname'} ne "" && $res eq $opts{'usrname'} ) {
         $usrnam    = "-usrname ".$opts{'usrname'};
      }
      #
      # Mapping files
      #
      my %map; my %hgrd; my %lmsk; my %datfil;
      my $hirespft = "off";
      if ( defined($opts{'hirespft'}) ) {
         $hirespft = "on";
      }
      my $merge_gis = "off";
      if ( defined($opts{'merge_gis'}) ) {
         $merge_gis = "on";
      }
      my $mopts  = "$queryopts -namelist default_settings $usrnam";
      my $mkopts = "-csmdata $CSMDATA -silent -justvalue -namelist clmexp $usrnam";
      my @typlist = ( "lak", "veg", "voc", "tex", "col", "hrv",
                        "fmx", "lai", "urb", "org", "glc", "glcregion", "utp", "wet",
		        "gdp", "peat","soildepth","abm", "vic", "ch4");
      if ( ! $opts{'fast_maps'} ) {
         push( @typlist, "topostats" );
      }
      foreach my $typ ( @typlist ) {
         my $lmask = `$scrdir/../../bld/queryDefaultNamelist.pl $mopts -options type=$typ,mergeGIS=$merge_gis,hirespft=$hirespft -var lmask`;
         $lmask = trim($lmask);
         my $hgrid_cmd = "$scrdir/../../bld/queryDefaultNamelist.pl $mopts -options type=$typ,hirespft=$hirespft -var hgrid";
         my $hgrid = `$hgrid_cmd`;
         if ($debug) {
           print "query to determine hgrid:\n    $hgrid_cmd \n\n";
         }
         $hgrid = trim($hgrid);
         my $filnm = `$scrdir/../../bld/queryDefaultNamelist.pl $mopts -options type=$typ -var mksrf_filename`;
         $filnm = trim($filnm);
         $hgrd{$typ} = $hgrid;
         $lmsk{$typ} = $lmask;
	 if ( $opts{'hgrid'} eq "usrspec" ) {
	     $map{$typ} = $opts{'usr_mapdir'}."/map_${hgrid}_${lmask}_to_${res}_nomask_aave_da_c${mapdate}\.nc";
	 } else {
	     $map{$typ} = `$scrdir/../../bld/queryDefaultNamelist.pl $queryfilopts -namelist clmexp -options frm_hgrid=$hgrid,frm_lmask=$lmask,to_hgrid=$res,to_lmask=nomask -var map`;
	 }	     
         $map{$typ} = trim($map{$typ});
         if ( $map{$typ} !~ /[^ ]+/ ) {
            die "ERROR: could NOT find a mapping file for this resolution: $res and type: $typ at $hgrid and $lmask.\n";
         }
         if ( ! defined($opts{'allownofile'}) && ! -f $map{$typ} ) {
            die "ERROR: mapping file for this resolution does NOT exist ($map{$typ}).\n";
          }
         my $typ_cmd = "$scrdir/../../bld/queryDefaultNamelist.pl $mkopts -options hgrid=$hgrid,lmask=$lmask,mergeGIS=$merge_gis$mkcrop -var $filnm";
         $datfil{$typ} = `$typ_cmd`;
         $datfil{$typ} = trim($datfil{$typ});
         if ( $datfil{$typ} !~ /[^ ]+/ ) {
            die "ERROR: could NOT find a $filnm data file for this resolution: $hgrid and type: $typ and $lmask.\n$typ_cmd\n\n";
         }
         if ( ! defined($opts{'allownofile'}) && ! -f $datfil{$typ} ) {
            die "ERROR: data file for this resolution does NOT exist ($datfil{$typ}).\n";
         }
      }
      #
      # Grid file from the pft map file or grid if not found
      #
      my $griddata    = trim($map{'veg'});
      if ( $griddata eq "" ) {
         $griddata = `$scrdir/../../bld/queryDefaultNamelist.pl $queryfilopts $usrnam -var fatmgrid`;
         if ( $griddata eq "" ) {
            die "ERROR: could NOT find a grid data file for this resolution: $res.\n";
         }
      }
      my $desc;
      my $desc_surfdat;
      #
      # Check if all urban single point dataset
      #
      my @all_urb = ( "1x1_camdenNJ","1x1_vancouverCAN", "1x1_mexicocityMEX", "1x1_urbanc_alpha" );
      my $all_urb = ".false.";
      my $urb_pt  = 0;
      foreach my $urb_res ( @all_urb ) {
         if ( $res eq $urb_res ) {
            $all_urb = ".true.";
            if ( $res ne "1x1_camdenNJ" ) { $urb_pt  = 1; }
         }
      }
      #
      # Always run at double precision for output
      #
      my $double = ".true.";
      #
      # Loop over each sim_year
      #
      RCP: foreach my $ssp_rcp ( @rcpaths ) {
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
            # If year is 1850-2000 actually run 1850-2015
            #
            if ( $sim_year eq "1850-2000" ) {
               my $actual = "1850-2015";
               print "For $sim_year actually run $actual\n";
               $sim_year = $actual;
            }
            my $urbdesc = "urb3den";
            my $resol    = "-res $hgrd{'veg'}";
            my $resolhrv = "-res $hgrd{'hrv'}";
            my $sim_yr0 = $sim_year;
            my $sim_yrn = $sim_year;
            my $transient = 0;
            if ( $sim_year =~ /([0-9]+)-([0-9]+)/ ) {
               $sim_yr0 = $1;
               $sim_yrn = $2;
               $transient = 1;
            }
            # determine simulation year to use for the surface dataset:
            my $sim_yr_surfdat = $sim_yr0;
            
            my $cmd    = "$scrdir/../../bld/queryDefaultNamelist.pl $queryfilopts $resol -options sim_year=${sim_yr_surfdat}$mkcrop -var mksrf_fvegtyp -namelist clmexp";
            my $vegtyp = `$cmd`;
            chomp( $vegtyp );
            if ( $vegtyp eq "" ) {
               die "** trouble getting vegtyp file with: $cmd\n";
            }
            my $cmd    = "$scrdir/../../bld/queryDefaultNamelist.pl $queryfilopts $resolhrv -options sim_year=${sim_yr_surfdat}$mkcrop -var mksrf_fvegtyp -namelist clmexp";
            my $hrvtyp = `$cmd`;
            chomp( $hrvtyp );
            if ( $hrvtyp eq "" ) {
               die "** trouble getting hrvtyp file with: $cmd\n";
            }
            my $options = "";
            my $crpdes  = sprintf("%2.2dpfts", $numpft);
            if ( $numpft == 16 ) {
               $crpdes .= "_Irrig";
            }
            if ( $mkcrop ne "" ) {
               $options = "-options $mkcrop";
            }
            $desc         = sprintf( "%s_%s_%s_simyr%4.4d-%4.4d", $ssp_rcp, $crpdes, $cmip_series, $sim_yr0, $sim_yrn );
            $desc_surfdat = sprintf( "%s_%s_%s_simyr%4.4d",       $ssp_rcp, $crpdes, $cmip_series, $sim_yr_surfdat  );

            my $fsurdat_fname_base = "";
            my $fsurdat_fname = "";
            if ( ! $opts{'no_surfdata'} ) {
               $fsurdat_fname_base = "surfdata_${res}_${desc_surfdat}_${sdate}";
               $fsurdat_fname = "${fsurdat_fname_base}.nc";
            }

            my $fdyndat_fname_base = "";
            my $fdyndat_fname = "";
            if ($transient) {
               $fdyndat_fname_base = "landuse.timeseries_${res}_${desc}_${sdate}";
               $fdyndat_fname = "${fdyndat_fname_base}.nc";
            }

            if (!$fsurdat_fname && !$fdyndat_fname) {
               die("ERROR: Tried to run mksurfdata_map without creating either a surface dataset or a landuse.timeseries file")
            }

            my $logfile_fname;
            my $namelist_fname;
            if ($fsurdat_fname_base) {
               $logfile_fname = "${fsurdat_fname_base}.log";
               $namelist_fname = "${fsurdat_fname_base}.namelist";
            }
            else {
               $logfile_fname = "${fdyndat_fname_base}.log";
               $namelist_fname = "${fdyndat_fname_base}.namelist";
            }

            my ($landuse_timeseries_text_file) = write_transient_timeseries_file(
                 $transient, $desc, $sim_yr0, $sim_yrn,
                 $queryfilopts, $resol, $resolhrv, $ssp_rcp, $mkcrop,
                 $sim_yr_surfdat);

            print "CSMDATA is $CSMDATA \n";
            print "resolution: $res ssp-rcp=$ssp_rcp sim_year = $sim_year\n";
            print "namelist: $namelist_fname\n";
            
            write_namelist_file(
                 $namelist_fname, $logfile_fname, $fsurdat_fname, $fdyndat_fname,
                 $glc_nec, $griddata, \%map, \%datfil, $double,
                 $all_urb, $no_inlandwet, $vegtyp, $hrvtyp, 
                 $landuse_timeseries_text_file, $setnumpft);

            #
            # Delete previous versions of files that will be created
            #
            system( "/bin/rm -f $fsurdat_fname $logfile_fname" );
            #
            # Run mksurfdata_map with the namelist file
            #
            my $exedir = $scrdir;
            if ( defined($opts{'exedir'}) ) {
               $exedir = $opts{'exedir'};
            }
            print "$exedir/mksurfdata_map < $namelist_fname\n";
            if ( ! $opts{'debug'} ) {
               system( "$exedir/mksurfdata_map < $namelist_fname" );
               if ( $? ) { die "ERROR in mksurfdata_map: $?\n"; }
            }
	    print "\n===========================================\n\n";

            #
            # If urban point, overwrite urban variables from previous surface dataset to this one
            #
            if ( $urb_pt && ! $opts{'no_surfdata'} ) {
               my $prvsurfdata = `$scrdir/../../bld/queryDefaultNamelist.pl $queryopts -var fsurdat`;
               if ( $? != 0 ) {
                  die "ERROR:: previous surface dataset file NOT found\n";
               }
               chomp( $prvsurfdata );
               my $varlist = "CANYON_HWR,EM_IMPROAD,EM_PERROAD,EM_ROOF,EM_WALL,HT_ROOF,THICK_ROOF,THICK_WALL,T_BUILDING_MIN,WIND_HGT_CANYON,WTLUNIT_ROOF,WTROAD_PERV,ALB_IMPROAD_DIR,ALB_IMPROAD_DIF,ALB_PERROAD_DIR,ALB_PERROAD_DIF,ALB_ROOF_DIR,ALB_ROOF_DIF,ALB_WALL_DIR,ALB_WALL_DIF,TK_ROOF,TK_WALL,TK_IMPROAD,CV_ROOF,CV_WALL,CV_IMPROAD,NLEV_IMPROAD,PCT_URBAN,URBAN_REGION_ID";
               print "Overwrite urban parameters with previous surface dataset values\n";
               $cmd = "ncks -A -v $varlist $prvsurfdata $fsurdat_fname";
               print "$cmd\n";
               if ( ! $opts{'debug'} ) { system( $cmd ); }
            }

         } # End of sim_year loop
      }    # End of ssp_rcp loop
   }
   close( $cfh );
   print "Successfully created fsurdat files\n";
