#!/usr/bin/env perl
#=======================================================================
#
#  This is a script to copy missing files in the CLM default namelist XML file
#  from a remote machine.
#
# Usage:
#
# scpFromDefaultNamelist.pl [options]
#
# To get help on options and usage:
#
# scpFromDefaultNamelist.pl -help
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
my @dirs = ( $cfgdir, "$cfgdir/perl5lib",
             "$cfgdir/../../../../scripts/ccsm_utils/Tools/perl5lib",
             "$cfgdir/../../../../models/utils/perl5lib",
           );
unshift @INC, @dirs;
my $result = eval "require XML::Lite";
if ( ! defined($result) ) {
   die <<"EOF";
** Cannot find perl module \"XML/Lite.pm\" from directories: @dirs **
EOF
}
require queryDefaultXML;

# Defaults
my $file = "$cfgdir/DefaultCLM_INPARM_Namelist.xml";
my $namelist = "clm_inparm";

sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options]
OPTIONS
     -csmdata "dir"                       Directory for head of csm inputdata.
     -file "file"                         Input xml file to read in by default ($file).
     -help  [or -h]                       Display this help.
     -namelist "namelistname"             Namelist name to read in by default ($namelist).
     -res  "resolution1,resolution2,..."  List of resolution to use for files.
                                          (At least one resolution is required).
     -scpfrom "remote-directory"          Remote directory to scp copy from.
                                          (required unless use -list option)
     -scpscript "filename"                Rather than copy files -- list files that 
                                          need to be copied to a file that can be 
                                          invoked as a script (file will need to be
                                          copied to remote site that contains the files).
     -silent [or -s]                      Don't do any extra printing.
EXAMPLES

  Copy all the files needed for resolution 10x15 from bangkok:

  $ProgName -res 10x15 -scpfrom bangkok.cgd:/fs/cgd/csm/inputdata

  Copy all the files needed for resolutions 10x15,4x5, and 64x128 from jaguar.
  with silent mode on so no extra printing is done:

  $ProgName -res 10x15,4x5,64x128 -scpfrom jaguar.ccs.ornl.gov:/fs/cgd/csm/inputdata -s

  Create a script that will be run on bluevista that contains the missing files
  and will copy them over to this machine.

  $ProgName -res 64x128,0.9x1.25,10x15 -scpscript missingfiles.csh

  copy "missingfiles.csh" over to bluevista (may need to edit remote hostname to copy to)

  on bluevista 

  setenv CSMDATA /fs/cgd/csm/inputdata
  ./missingfiles.csh
  
EOF
}

#-----------------------------------------------------------------------------------------------

sub GetListofNeededFiles {
#
# Get list of files that are needed to be copied to disk from the XML file.
#
  my $inputopts_ref = shift;
  my $settings_ref  = shift;
  my $files_ref     = shift;
  
  my %inputopts = %$inputopts_ref;
  my %settings  = %$settings_ref;

  my $defaults_ref = &queryDefaultXML::ReadDefaultXMLFile( $inputopts_ref, $settings_ref );
  my %defaults = %$defaults_ref;
  my @keys     = keys(%defaults);
  my $csmdata  = $defaults{'csmdata'}{'value'};
  $$inputopts_ref{'csmdata'} = $csmdata;
  my $printing = $inputopts{'printing'};
  foreach my $var ( @keys ) {
     my $value   = $defaults{$var}{'value'};
     my $isafile = $defaults{$var}{'isfile'};
     # If is a file
     if ( $isafile ) {
        # Test that this file exists
        if ( -f "$value" ) {
           print "File $value exists\n" if $printing;
        } else {
        # If doesn't exist add it to the list of files to copy
           $value    =~ m#$csmdata/(.+?)/([^/]+)$#;
           my $dir   = $1;
           my $file  = $2;
           my $cfile = $inputopts{'scpfrom'} . "$dir/$file";
           my @dfiles;
           if ( defined($$files_ref{$dir}) ) {
              my $dir_ref = $$files_ref{$dir};
              @dfiles     = @$dir_ref;
              my $match = 0;
              foreach my $i ( @dfiles ) {
                if ( $i eq $cfile ) { $match = 1; }
              }
              if ( $match == 0 ) { push( @dfiles, $cfile ); }
           } else {
              @dfiles = ( "$cfile" );
           }
           if ( ! defined($$files_ref{$dir}) ) {
              print "           ADD $cfile to list to copy\n";
           }
           $$files_ref{$dir} = \@dfiles;
        }
     }
  }
}

#-----------------------------------------------------------------------------------------------

  my %opts = ( file       => $file,
               namelist   => $namelist,
               res        => undef,
               csmdata    => undef,
               silent     => undef,
               scpfrom    => undef,
               help       => undef,
               scpscript  => undef,
             );

  my $cmdline = "@ARGV";
  GetOptions(
        "f|file=s"     => \$opts{'file'},
        "n|namelist=s" => \$opts{'namelist'},
        "r|res=s"      => \$opts{'res'},
        "csmdata=s"    => \$opts{'csmdata'},
        "scpscript=s"  => \$opts{'scpscript'},
        "scpfrom=s"    => \$opts{'scpfrom'},
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
  #
  # Handle the scpscript case 
  if ( defined($opts{'scpscript'}) ) { 
     $opts{'scpfrom'} = '$CSMDATA/'; 
  }
  # Check for required arguments
  foreach my $req ( "scpfrom", "res" ) {
     if ( ! defined($opts{$req}) ) {
         print "ERROR: $req NOT set and it is a required argument\n";
         usage();
     }
  }
  my %inputopts;
  $inputopts{'file'}     = $opts{'file'};
  $inputopts{'namelist'} = $opts{'namelist'};
  $inputopts{'printing'} = $printing;
  $inputopts{'ProgName'} = $ProgName;
  $inputopts{'cmdline'}  = $cmdline;
  $inputopts{'scpfrom'}  = $opts{'scpfrom'};
  if ( ! defined($opts{'csmdata'}) ) {
     $inputopts{'csmdata'} = "default";
  } else {
     $inputopts{'csmdata'} = $opts{'csmdata'};
  }
  $inputopts{'config'} = "noconfig";
  my @resolutions = split( /,/, $opts{'res'} );

  my %files;
  #
  # Loop over all resolutions asked for
  #
  foreach my $res ( @resolutions ) {
     $inputopts{'res'} = $res;
     print "Resolution = $res\n" if $printing;
     my %settings;
     #
     # Loop for all possible land masks
     #
     foreach my $mask ( "USGS", "gx3v5", "gx1v3", "gx1v4", "gx1v5", "navy", "test", "tx0.1v2" ) {
        print "Mask = $mask \n" if $printing;
        $settings{'mask'} = $mask;
        #
        # Loop over all possible simulation years 
        #
        foreach my $sim_year ( "1890", "2000", "2100" ) {
           print "Mask = $mask \n" if $printing;
           $settings{'sim_year'} = $sim_year;   # 1890, 2000, 2100

           #
           # Loop over all possible BGC seetings
           #
           foreach my $bgc ( "none", "cn", "casa", "dgvm" ) {
              print "BGC = $bgc\n" if $printing;
              $settings{'BGC'} = $bgc;
              &GetListofNeededFiles( \%inputopts, \%settings, \%files );
           }
        }
     }
  }
  #
  # Loop over directories that need to have files copied into
  #
  my $hostname;
  my $csmdata = $inputopts{'csmdata'};
  if ( defined($opts{'scpscript'}) ) { 
     open( OUT, ">".$opts{'scpscript'} ) || die "ERROR: trouble opening output file:" .
                                                 $opts{'scpscript'};
     $hostname = `hostname`;
     chomp( $hostname );
     print OUT  "#!/usr/bin/csh -f\n";
     print OUT  "#\n";
     print OUT  "# Script to copy needed files to remote site: $hostname\n";
     print OUT  "# First copy this script over to a machine that has the missing files.\n";
     print OUT  "# Second make sure the remote hostname is correct.\n";
     print OUT  "# Third set the env variable CSMDATA to the root of CCSM inputdata.\n";
     print OUT  "# Lastly invoke the script and enter passwords as needed.\n";
     print OUT  "#\n";
  }
  foreach my $dir ( sort(keys(%files)) ) {
     my $files_ref = $files{$dir};
     my @files     = @$files_ref;
     if ( ! defined($opts{'scpscript'}) ) {
        print   "scp @files $csmdata/$dir/.\n";
        system( "scp @files $csmdata/$dir/." )
     } else {
        print OUT  "scp @files $hostname:$csmdata/$dir/.\n";
     }
  }
  if ( defined($opts{'scpscript'}) ) { 
     close( OUT );
     chmod( 0755, $opts{'scpscript'} ) || die "ERROR:: error changing execute permission on file: " . $opts{'scpscript'};
  }
  if ( $printing ) {
     print "\n\nSuccessful\n\n"
  }
