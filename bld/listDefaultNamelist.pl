#!/usr/bin/env perl
#=======================================================================
#
#  This is a script to svn export missing files in the CLM default 
#  namelist XML file
#
# Usage:
#
# listDefaultNamelist.pl [options]
#
# To get help on options and usage:
#
# listDefaultNamelist.pl -help
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
# The namelist defaults file contains default values for all required namelist variables.
my @nl_defaults_files = ( "$cfgdir/namelist_files/namelist_defaults_overall.xml",
                          "$cfgdir/namelist_files/namelist_defaults_clm.xml",
                          "$cfgdir/namelist_files/namelist_defaults_drv.xml",
                          "$cfgdir/namelist_files/namelist_defaults_datm.xml" );
my $list = "clm.input_data_list";

sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName [options]
OPTIONS
     -help  [or -h]                       Display this help.
     -csmdata [or -d]                     Path to CSMDATA.
     -res  "resolution1,resolution2,..."  List of resolution to use for files.
                                          (At least one resolution is required)
                                          (If res is "all" will run over all resolutions)
     -silent [or -s]                      Don't do any extra printing.
EXAMPLES

  List all the files needed for resolution 10x15 to the file $list.

  $ProgName -res 10x15

  List all the files needed for resolutions 10x15,4x5, and 64x128 from jaguar.
  with silent mode on so no extra printing is done:

  $ProgName -res 10x15,4x5,64x128 -s

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
  
  my $defaults_ref = &queryDefaultXML::ReadDefaultXMLFile( $inputopts_ref, $settings_ref );
  my @keys     = keys(%$defaults_ref);
  my $csmdata  = $$inputopts_ref{'csmdata'};
  my $printing = $$inputopts_ref{'printing'};
  foreach my $var ( @keys ) {
     my $value   = $$defaults_ref{$var}{'value'};
     my $isafile = $$defaults_ref{$var}{'isfile'};
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
           my $cfile = $$inputopts_ref{'scpfrom'} . "$dir/$file";
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

  my %opts = ( 
               res        => undef,
               silent     => undef,
               csmdata    => "default",
               list       => $list,
               help       => undef,
             );

  my $cmdline = "@ARGV";
  GetOptions(
        "d|csmdata=s"  => \$opts{'csmdata'},
        "r|res=s"      => \$opts{'res'},
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
  # Check for required arguments
  #
  foreach my $req ( "res", "list" ) {
     if ( ! defined($opts{$req}) ) {
         print "ERROR: $req NOT set and it is a required argument\n";
         usage();
     }
  }
  my %inputopts;
  $inputopts{'nldef_file'}     = "$cfgdir/namelist_files/namelist_definition.xml";
  $inputopts{'empty_cfg_file'} = "$cfgdir/config_files/config_definition.xml";

  my $definition = Build::NamelistDefinition->new( $inputopts{'nldef_file'} );
  my $cfg = Build::Config->new( $inputopts{'empty_cfg_file'} );

  # Resolutions...
  my @resolutions;
  if ( $opts{'res'} eq "all" ) {
     @resolutions = $definition->get_valid_values( "res", 'noquotes'=>1 );
  } else {
     @resolutions = split( /,/, $opts{'res'} );
  }

  # Input options
  $inputopts{'files'}     = \@nl_defaults_files;
  $inputopts{'printing'}  = $printing;
  $inputopts{'ProgName'}  = $ProgName;
  $inputopts{'cmdline'}   = $cmdline;
  $inputopts{'cfgdir'}    = $cfgdir;
  $inputopts{'csmdata'}   = $opts{'csmdata'};
  $inputopts{'config'}    = "noconfig";
  my %files;
  #
  # Loop over all resolutions asked for: 1.9x2.5, 10x15, 64x128 etc.
  #
  foreach my $res ( @resolutions ) {
     if ( ! $definition->is_valid_value( "res", "'$res'" )  ) {
        die "ERROR: Input resolution: $res is NOT a valid resolution\n";
     }
     $inputopts{'hgrid'} = $res;
     print "Resolution = $res\n" if $printing;
     my %settings;
     #
     # Loop for all possible land masks: USGS, gx1v6, gx3v5 etc.
     #
     foreach my $mask ( $definition->get_valid_values( "mask", 'noquotes'=>1 ) ) {
        print "Mask = $mask \n" if $printing;
        $settings{'mask'} = $mask;
        #
        # Loop over all possible simulation year: 1890, 2000, 2100 etc.
        #
        foreach my $sim_year ( $definition->get_valid_values( "sim_year", 'noquotes'=>1 ) ) {
           print "sim_year = $sim_year\n" if $printing;
           $settings{'sim_year'} = $sim_year;   

           #
           # Loop over all possible BGC seetings: none, cn, casa etc.
           #
           foreach my $bgc ( $cfg->get_valid_values( "bgc" ) ) {
              print "bgc = $bgc\n" if $printing;
              $settings{'bgc'} = $bgc;
              $inputopts{'namelist'} = "clm_inparm";
              &GetListofNeededFiles( \%inputopts, \%settings, \%files );
           }
        }
        #
        # Now also do the datm namelist (only for mask and resolution)
        #
        my $nml = "datm_dshr_in";
        print "nml = $nml\n" if $printing;
        $inputopts{'namelist'} = $nml;
        &GetListofNeededFiles( \%inputopts, \%settings, \%files );
     }
  }
  #
  # Loop over directories that need to have files copied into
  #
  my $hostname;
  my $csmdata = $inputopts{'csmdata'};
  open( OUT, ">$list" ) || die "ERROR: trouble opening output file: $list";
  foreach my $dir ( sort(keys(%files)) ) {
     if ( $dir eq "."     ) { next; }
     if ( $dir eq "/"     ) { next; }
     if ( $dir eq "\n"    ) { next; }
     if ( $dir eq ""      ) { next; }
     if ( ! defined($dir) ) { next; }
     my $files_ref = $files{$dir};
     my @files     = @$files_ref;
     foreach my $file ( @files ) {
        if ( $file !~ /\n$/ ) { $file = "$file\n"; }
        print OUT  "file = \$DIN_LOC_ROOT/$file";
     }
  }
  close( OUT );
  if ( $printing ) {
     print "\n\nSuccessful\n\n"
  }
