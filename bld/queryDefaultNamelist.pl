#!/usr/bin/env perl
#=======================================================================
#
#  This is a script to read the CLM namelist XML file
#
# Usage:
#
# queryDefaultNamelist.pl [options]
#
# To get help on options and usage:
#
# queryDefaultNamelist.pl -help
#
#=======================================================================

use Cwd;
use strict;
#use diagnostics;
use Getopt::Long;
use English;

#-----------------------------------------------------------------------------------------------

#Figure out where configure directory is and where can use the XML/Lite module from
my $ProgName;
($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!; # name of program
my $ProgDir = $1;                         # name of directory where program lives

my $cwd = getcwd();  # current working directory
my $cfgdir;

if ($ProgDir) { $cfgdir = $ProgDir; }
else { $cfgdir = $cwd; }

#-----------------------------------------------------------------------------------------------
# Add $cfgdir to the list of paths that Perl searches for modules
my @dirs = ( $cfgdir, "$cfgdir/perl5lib", "$cfgdir/../../../../scripts/ccsm_utils/Tools/perl5lib"
);
unshift @INC, @dirs;
my $result = eval "require XML::Lite";
if ( ! defined($result) ) {
   die <<"EOF";
** Cannot find perl module \"XML/Lite.pm\" from directories: @dirs **
EOF
}
require Build::Config;
require queryDefaultXML;

# Defaults
my $namelist = "clm_inparm";
my $file = "$cfgdir/DefaultCLM_INPARM_Namelist.xml";
my $config = "config_cache.xml";


sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options]
OPTIONS
     -config "file"                       CLM build configuration file created by configure.
     -ccsm                                CCSM mode set csmdata to \$DIN_LOC_ROOT.
     -csmdata "dir"                       Directory for head of csm inputdata.
     -demand                              Demand that something is returned.
     -file "file"                         Input xml file to read in by default ($file).
     -help  [or -h]                       Display this help.
     -justvalue                           Just display the values (NOT key = value).
     -var "varname"                       Variable name to match.
     -namelist "namelistname"             Namelist name to read in by default ($namelist).
     -onlyfiles                           Only output filenames.
     -options "item=value,item2=value2"   Set options to query for when matching.
                                          (comma delimited, with equality to set value).
     -res  "resolution"                   Resolution to use for files.
     -scpto "dir"                         Copy files to given directory.
                                          (also enables -onlyfiles, -test and -justvalue options)
     -silent [or -s]                      Don't do any extra printing.
     -test   [or -t]                      Test that files exists.
EXAMPLES

  To list all fsurdat files that match the resolution: 10x15:

  $ProgName -var "fsurdat" -res 10x15

  To only list files that match T42 resolution (or are for all resolutions)

  $ProgName  -onlyfiles  -res 64x128

  To test that all of the files exist on disk under a different default inputdata

  $ProgName  -onlyfiles -test -csmdata /spin/proj/ccsm/inputdata

  To query the CAM file for files that match particular configurations

  $ProgName  -onlyfiles \
  -file \$CAM_CFGDIR/DefaultCAM_INPARM_Namelist.xml -namelist cam_inparm \
  -options DYNAMICS=fv,PLEV=26,PHYSICS=cam1,OCEANMODEL=DOM,CHEMISTRY=

  To copy files that match to a different location

  $ProgName -onlyfiles -test -scpto bangkok:/fs/cgd/csm/inputdata

EOF
}

#-----------------------------------------------------------------------------------------------

  my %opts = ( file       => $file,
               namelist   => $namelist,
               var        => undef,
               res        => undef,
               config     => undef,
               ccsm       => undef, 
               csmdata    => undef,
               demand     => undef,
               test       => undef,
               onlyfiles  => undef,
               justvalues => undef,
               scpto      => undef,
               silent     => undef,
               help       => undef,
               options    => undef,
             );

  my $cmdline = "@ARGV";
  GetOptions(
        "f|file=s"     => \$opts{'file'},
        "n|namelist=s" => \$opts{'namelist'},
        "v|var=s"      => \$opts{'var'},
        "r|res=s"      => \$opts{'res'},
        "config=s"     => \$opts{'config'},
        "ccsm"         => \$opts{'ccsm'},
        "csmdata=s"    => \$opts{'csmdata'},
        "demand"       => \$opts{'demand'},
        "options=s"    => \$opts{'options'},
        "t|test"       => \$opts{'test'},
        "onlyfiles"    => \$opts{'onlyfiles'},
        "justvalues"   => \$opts{'justvalues'},
        "scpto=s"      => \$opts{'scpto'},
        "s|silent"     => \$opts{'silent'},
        "h|elp"        => \$opts{'help'},
  ) or usage();

  # Check for unparsed arguments
  if (@ARGV) {
      print "ERROR: unrecognized arguments: @ARGV\n";
      usage();
  }
  if ( $opts{'help'} ) {
      usage();
  }
  # Set if should do extra printing or not (if silent mode is not set)
  my $printing = 1;
  if ( defined($opts{'silent'}) ) {
      $printing = 0;
  }
  # Get list of options from command-line into the settings hash
  my %settings;
  if ( defined($opts{'options'}) ) {
     my @optionlist = split( ",", $opts{'options'} );
     foreach my $item ( @optionlist ) {
        my ($key,$value) = split( "=", $item );
        $settings{$key} = $value;
     }
  }
  my $csmdata = "";
  if ( defined($opts{'scpto'}) ) {
     if ( ! defined($opts{'justvalues'}) ) { print "When -scpto option used, -justvalues is set as well\n"; }
     if ( ! defined($opts{'onlyfiles'}) )  { print "When -scpto option used, -onlyfiles is set as well\n"; }
     if ( ! defined($opts{'test'}) )       { print "When -scpto option used, -test is set as well\n"; }
     $opts{'justvalues'} = 1;
     $opts{'onlyfiles'}  = 1;
     $opts{'test'}       = 1;
  }
  my %inputopts;
  $inputopts{file}     = $opts{file};
  $inputopts{namelist} = $opts{namelist};
  $inputopts{printing} = $printing;
  $inputopts{ProgName} = $ProgName;
  $inputopts{cmdline}  = $cmdline;
  if ( ! defined($opts{csmdata}) ) {
     $inputopts{csmdata} = "default";
  } else {
     $inputopts{csmdata} = $opts{csmdata};
  }
  if ( defined($opts{ccsm}) ) {
     $inputopts{csmdata} = '$DIN_LOC_ROOT';
  }
  if ( ! defined($opts{config}) ) {
     $inputopts{config} = "noconfig";
  } else {
     $inputopts{config} = $opts{config};
  }
  if ( ! defined($opts{res}) ) {
     $inputopts{res} = "any";
  } else {
     $inputopts{res} = $opts{res};
  }

  my $defaults_ref = &queryDefaultXML::ReadDefaultXMLFile( \%inputopts, \%settings );
  my %defaults = %$defaults_ref;
  my @keys = keys(%defaults);
  if ( defined($opts{'var'}) ) {
     my @newkeys;
     foreach my $var ( @keys ) {
        if ( $var eq $opts{'var'} ) {
            @newkeys = ( $opts{var} );
        }
     }
     @keys = @newkeys;
  }
  if ( defined($opts{'demand'}) && ($#keys == -1) ) {
     die "($ProgName $cmdline) ERROR:: demand option is set and nothing was found.\n";
  }
  my $print;
  foreach my $var ( @keys ) {
     $print = 1;
     my $value   = $defaults{$var}{value};
     my $isadir  = $defaults{$var}{isdir};
     my $isafile = $defaults{$var}{isfile};
     # If onlyfiles option set do NOT print if is NOT a file
     if ( defined($opts{'onlyfiles'}) && (! $isafile) ) {
        $print = undef;
     }
     # If is a directory
     if ( $isadir ) {
        # Test that this directory exists
        if ( defined($opts{'test'})  && defined($print) ) {
           print "Test that directory $value exists\n" if $printing;
           if ( ! -d "$value" ) {
              die "($ProgName) ERROR:: directory $value does NOT exist!\n";
           }
        }
     }
     # If is a file
     if ( $isafile ) {
        # Test that this file exists
        if ( defined($opts{'test'})  && defined($print) ) {
           print "Test that file $value exists\n" if $printing;
           if ( ! -f "$value" ) {
              die "($ProgName) ERROR:: file $value does NOT exist!\n";
           }
        }
        # Copy to remote location
        if ( defined($opts{'scpto'}) && defined($print) ) {
           my $csmdata = $inputopts{csmdata};
           my $data = $value;
           if ( $value =~ /^$csmdata\/(.+)$/ ) {
              $data = $1;
           } else {
              die "($ProgName) ERROR:: parsing $csmdata directory out of value: $value\n";
           }
           print "scp $csmdata/$data $opts{'scpto'}/$data\n";
           system( "scp $csmdata/$data $opts{'scpto'}/$data" )
        }
     }
     # If a file or directory -- add quotes around the value
     if ( (! defined($opts{'justvalues'}) ) && ($isafile || $isadir) ) {
       $value = "\'$value\'";
     }
     if ( defined($print) ) {
        if ( ! defined($opts{'justvalues'})  ) {
           print "$var = ";
        }
        print "$value\n";
     }
  }
  if ( $printing && defined($opts{'test'}) ) {
     print "\n\nTesting was successful\n\n"
  }

