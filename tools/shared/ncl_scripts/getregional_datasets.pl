#!/usr/bin/env perl
#=======================================================================
#
# Extract out regional datasets from the global datasets.
#
# Usage:
#
# getregional_datasets.pl
#
#  Erik Kluzek
#  Aug/28/2009
#  $Id$
#  $HeadURL;
#
#=======================================================================

use Cwd;
use strict;
#use diagnostics;
use English;
use Getopt::Long;
use IO::File;

#-----------------------------------------------------------------------------------------------
# Set the directory that contains this scripts.  If the command was issued using a 
# relative or absolute path, that path is in $ProgDir.  Otherwise assume the
# command was issued from the current working directory.

(my $ProgName = $0) =~ s!(.*)/!!;      # name of this script
my $ProgDir = $1;                      # name of directory containing this script -- may be a
                                       # relative or absolute path, or null if the script is in
                                       # the user's PATH
my $cmdline = "@ARGV";                 # Command line arguments to script
my $cwd = getcwd();                    # current working directory
my $scrdir;                            # absolute pathname of directory that contains this script
my $nm = "ProgName::";                 # name to use if script dies
if ($ProgDir) { 
    $scrdir = absolute_path($ProgDir);
} else {
    $scrdir = $cwd;
}

# Default resolution

my $res        = "1.9x2.5";
my $rcp        = "-999.9";
my $mask       = "gx1v6";
my $sim_year   = "2000";
my $sim_yr_rng = "constant";
my $mycsmdata  = $ENV{'HOME'} . "/inputdata";

#-----------------------------------------------------------------------------------------------

sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options]	   Extracts out files for a single box region from the global
                                   grid for the region of interest. Choose a box determined by
                                   the NorthEast and SouthWest corners.
REQUIRED OPTIONS
     -mydataid "name" [or -id]     Your name for the region that will be extracted.               (REQUIRED)
                                   Recommended name: grid-size_global-resolution_location (?x?pt_f??_????)
                                   (i.e. 12x13pt_f19_alaskaUSA for 12x13 grid cells from the f19 global resolution over Alaska)
     -NE_corner "lat,lon" [or -ne] North East corner latitude and longitude                       (REQUIRED)
     -SW_corner "lat,lon" [or -sw] South West corner latitude and longitude                       (REQUIRED)
OPTIONS
     -debug [or -d]                Just debug by printing out what the script would do.
                                   This can be useful to find the size of the output area.
     -help [or -h]                 Print usage to STDOUT.
     -mask "landmask"              Type of land-mask (i.e. navy, gx3v7, gx1v6 etc.) (default $mask)
     -mycsmdata "dir"              Root directory of where to put your csmdata. 
                                   (default $mycsmdata or value of CSMDATA env variable)
     -nomv                         Do NOT move datasets to final location, just leave them in current directory
     -res "resolution"             Global horizontal resolution to extract data from (default $res).
     -rcp "pathway"                Representative concentration pathway for future scenarios 
                                   Only used when simulation year range ends in a future
                                   year, such as 2100.
                                   (default $rcp).
     -sim_year   "year"            Year to simulate for input datasets (i.e. 1850, 2000) (default $sim_year)
     -sim_yr_rng "year-range"      Range of years for transient simulations 
                                   (i.e. 1850-2000, 1850-2100,  or constant) (default $sim_yr_rng)

     -verbose [or -v]              Make output more verbose.
EOF
}

sub get_latlon {
#
# Return the latitude and longitude of the input string and validate it
#
  my $string = shift;
  my $desc   = shift;

  my $lat = undef;
  my $lon = undef;
  my $valreal1 = "[+-]?[0-9]*\.?[0-9]+[EedDqQ]?[0-9+-]*";

  if ( $string =~ /^($valreal1)\s*,\s*($valreal1)$/ ) {
     $lat = $1;
     $lon = $2;
  } else {
     die <<"EOF";
** $ProgName - Error in entering latitude/longitude for $desc **
EOF
  }
  if ( ($lat < -90.) || ($lat >  90.0) ) {
     die <<"EOF";
** $ProgName - Bad value for latitude (=$lat) for $desc **
EOF
  }
  if ( ($lon < 0.)   || ($lon > 360.0) ) {
     die <<"EOF";
** $ProgName - Bad value for longitude  (=$lat) for $desc **
EOF
  }
  return( $lat, $lon );

}

#-----------------------------------------------------------------------------------------------

# Process command-line options.

my %opts = ( 
              mask             => $mask,
              sim_year         => $sim_year,
              sim_yr_rng       => $sim_yr_rng,
              mycsmdata        => undef,
              mydataid         => undef,
              SW_corner        => undef,
              NE_corner        => undef,
              res              => $res,
              rcp              => $rcp,
              help             => 0, 
              nomv             => 0,
              verbose          => 0,
              debug            => 0,
           );
GetOptions(
    "sw|SW_corner=s"   => \$opts{'SW_corner'},
    "ne|NE_corner=s"   => \$opts{'NE_corner'},
    "mycsmdata=s"      => \$opts{'mycsmdata'},
    "id|mydataid=s"    => \$opts{'mydataid'},
    "sim_year=i"       => \$opts{'sim_year'},
    "mask=s"           => \$opts{'mask'},
    "res=s"            => \$opts{'res'},
    "rcp=f"            => \$opts{'rcp'},
    "nomv"             => \$opts{'nomv'},
    "sim_yr_rng=s"     => \$opts{'sim_yr_rng'},
    "h|help"           => \$opts{'help'},
    "d|debug"          => \$opts{'debug'},
    "v|verbose"        => \$opts{'verbose'},
)  or usage();

# Give usage message.
usage() if $opts{'help'};

# Check for unparsed arguments
if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
}

if ( ! defined($opts{'SW_corner'}) || ! defined($opts{'NE_corner'}) ) {
    print "ERROR: MUST set both SW_corner and NE_corner\n";
    usage();
}
if ( ! defined($opts{'mydataid'}) ) {
    print "ERROR: MUST set mydataid\n";
    usage();
}

my ($S_lat,$W_lon) = get_latlon( $opts{'SW_corner'}, "SW" );
my ($N_lat,$E_lon) = get_latlon( $opts{'NE_corner'}, "NE" );

if ( $N_lat <= $S_lat ) {
    print "ERROR: NE corner latitude less than or equal to SW corner latitude\n";
    usage();
}
if ( $E_lon <= $W_lon ) {
    print "ERROR: NE corner longitude less than or equal to SW corner longitude\n";
    usage();
}

my $inputdata_rootdir;

if (defined($opts{'mycsmdata'})) {
    $inputdata_rootdir = $opts{'mycsmdata'};
}
elsif (defined $ENV{'CSMDATA'}) {
    $inputdata_rootdir = $ENV{'CSMDATA'};
}
else {
    die "$ProgName - ERROR: Personal CESM inputdata root directory must be specified by either -mycsmdata argument\n" .
        " or by the CSMDATA environment variable. :";
}
(-d $inputdata_rootdir)  or  die <<"EOF";
** $ProgName - CESM inputdata root is not a directory: \"$inputdata_rootdir\" **
EOF
$ENV{'DIN_LOC_ROOT'} = $inputdata_rootdir;

print "CESM inputdata root directory: $inputdata_rootdir\n";

#-----------------------------------------------------------------------------------------------
my $debug;
if ( $opts{'debug'} ) {
  $debug = "DEBUG=TRUE";  
}
my $nomv;
if ( $opts{'nomv'} ) {
  $nomv = "NOMV=TRUE";  
}
my $print;
if ( $opts{'verbose'} ) {
  $print = "PRINT=TRUE";  
}

my $cmd = "env S_LAT=$S_lat W_LON=$W_lon N_LAT=$N_lat E_LON=$E_lon " . 
          "SIM_YR=$opts{'sim_year'} SIM_YR_RNG=$opts{'sim_yr_rng'} MASK=$opts{'mask'} " .
          "CLM_USRDAT_NAME=$opts{'mydataid'} RCP=$opts{'rcp'} RES=$opts{'res'} MYCSMDATA=$inputdata_rootdir " .
          "$debug $print $nomv ncl $scrdir/getregional_datasets.ncl";

print "Execute: $cmd\n";
system( $cmd );

#-------------------------------------------------------------------------------

sub absolute_path {
#
# Convert a pathname into an absolute pathname, expanding any . or .. characters.
# Assumes pathnames refer to a local filesystem.
# Assumes the directory separator is "/".
#
  my $path = shift;
  my $cwd = getcwd();  # current working directory
  my $abspath;         # resulting absolute pathname

# Strip off any leading or trailing whitespace.  (This pattern won't match if
# there's embedded whitespace.
  $path =~ s!^\s*(\S*)\s*$!$1!;

# Convert relative to absolute path.

  if ($path =~ m!^\.$!) {          # path is "."
      return $cwd;
  } elsif ($path =~ m!^\./!) {     # path starts with "./"
      $path =~ s!^\.!$cwd!;
  } elsif ($path =~ m!^\.\.$!) {   # path is ".."
      $path = "$cwd/..";
  } elsif ($path =~ m!^\.\./!) {   # path starts with "../"
      $path = "$cwd/$path";
  } elsif ($path =~ m!^[^/]!) {    # path starts with non-slash character
      $path = "$cwd/$path";
  }

  my ($dir, @dirs2);
  my @dirs = split "/", $path, -1;   # The -1 prevents split from stripping trailing nulls
                                     # This enables correct processing of the input "/".

  # Remove any "" that are not leading.
  for (my $i=0; $i<=$#dirs; ++$i) {
      if ($i == 0 or $dirs[$i] ne "") {
          push @dirs2, $dirs[$i];
      }
  }
  @dirs = ();

  # Remove any "."
  foreach $dir (@dirs2) {
      unless ($dir eq ".") {
          push @dirs, $dir;
      }
  }
  @dirs2 = ();

  # Remove the "subdir/.." parts.
  foreach $dir (@dirs) {
    if ( $dir !~ /^\.\.$/ ) {
        push @dirs2, $dir;
    } else {
        pop @dirs2;   # remove previous dir when current dir is ..
    }
  }
  if ($#dirs2 == 0 and $dirs2[0] eq "") { return "/"; }
  $abspath = join '/', @dirs2;
  return( $abspath );
}

#-------------------------------------------------------------------------------

