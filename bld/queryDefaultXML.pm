#=======================================================================
#
#  This is a perl module to read in a DefaultNamelist XML file.
#
#=======================================================================
use strict;
package queryDefaultXML;

#-------------------------------------------------------------------------------

sub convert_keys_touppercase
{
   my (%hash) = @_;

   my %outhash;
   foreach my $key ( keys(%hash) ) {
      $_ = $key;
      tr/a-z/A-Z/;
      $outhash{$_} = $hash{$key};
   }
   return( %outhash );
}

#-------------------------------------------------------------------------------

sub read_cfg_file
#
# Read in the configuration cache XML file on the build-time configuration
#
{
    my ($file) = @_;

    if ( ! -f "$file" ) {
      die "file $file does not exist\n";
    }
    my $xml = XML::Lite->new( $file );
    my $root = $xml->root_element();

    # Check for valid root node
    my $name = $root->get_name();
    $name eq "config_bld" or die
        "file $file is not a CLM configuration file\n";

    # Get source and build directories
    my $dirs = $xml->elements_by_name( "directories" );
    my %dirs = $dirs->get_attributes();
    my %DIRS = convert_keys_touppercase( %dirs );

    # Get cppvars
    my $cppvars = $xml->elements_by_name( "cppvars" );
    my %cppvars = $cppvars->get_attributes();
    my %CPPVARS = convert_keys_touppercase( %cppvars );

    # Get settings for Makefile (parallelism and library locations)
    my $make = $xml->elements_by_name( "makefile" );
    my %make = $make->get_attributes();
    my %MAKE = convert_keys_touppercase( %make );

    return %DIRS, %CPPVARS, %MAKE;
}

#-------------------------------------------------------------------------------

sub ReadDefaultXMLFile {
#
# Read in the default XML file for the default namelist settings
#
  my $opts_ref     = shift;
  my $settings_ref = shift;

  my %opts = %$opts_ref;
  my %settings = %$settings_ref;

  # Error check that input and opts hash has the expected variables
  my $ProgName     = $opts{'ProgName'};
  my $nm = "${ProgName}::ReadDefaultXMLFile";
  my @required_list = ( "file", "config", "namelist", "csmdata", "res", 
                        "printing", "ProgName", "cmdline" );
  foreach my $var ( @required_list ) {
     if ( ! defined($opts{$var}) ) {
        die "ERROR($nm): Required input variable $var was not found\n";
     }
  }
  my $printing = $opts{'printing'};
  my $cmdline  = $opts{'cmdline'};
  # Initialize some local variables
  my $file     = $opts{'file'};
  my $namelist = $opts{'namelist'};

  my %bld;
  if ( $opts{'config'} ne "noconfig" ) {
     %bld = read_cfg_file( $opts{'config'} );
  }

  # Open file
  print "($nm) Read: $file\n" if $printing;
  my $xml = XML::Lite->new( $file );
  if ( ! defined($xml) ) {
    die "ERROR($nm): Trouble opening or reading $file\n";
  }
  #
  # Find the namelist element for this namelist
  #
  my $elm = $xml->root_element( );
  my @list = $xml->elements_by_name( $namelist );
  if ( $#list < 0 ) {
    die "ERROR($nm): could not find the main $namelist namelist element in $file\n";
  }
  if ( $#list != 0 ) {
    die "ERROR($nm): $namelist namelist element in $file is duplicated, there should only be one\n";
  }
  $elm = $list[0];
  my @children = $elm->get_children();
  if ( $#children < 0 ) {
    die "ERROR($nm): There are no sub-elements to the $namelist element in $file\n";
  }
  #
  # Get the names for each element in the list into a hash
  #
  my %names;
  my @settings_keys = keys( %settings );
  foreach my $child ( @children ) {
    my $name = $child->get_name();
    $names{$name} = $child->get_text();
  }
  #
  # Go through the sub-elements to the namelist element
  #
  my %defaults;
  foreach my $child ( @children ) {
    #
    # Get the attributes for each namelist element
    # The attributes describe either config settings that need to match
    # or other namelist elements that need to match
    #
    my %atts = $child->get_attributes;
    # Name of element, and it's associated value
    my $name = $child->get_name();
    my $value =  $child->get_text();
    $value =~ s/\n//g;   # Get rid of extra returns
    my @keys = keys(%atts);
    my $set = 1;
    # If csmdata set to default and file contains it
    if ( ($$opts_ref{'csmdata'} eq "default")  && ($name eq "csmdata") ) {
      $$opts_ref{'csmdata'} = $value;
    }
    # If Version
    if ( $name eq "version" ) {
       if ( $value =~ 
            /(URL: https:\/\/[a-zA-Z0-9._-]+\/)([a-zA-Z0-9\/._-]+)(\/bld\/.+)/
) {
          $value = $2;
       }
    }
    if ( $#keys >= 0 ) {
      #
      # Check that all values match the appropriate settings
      #
      foreach my $key ( @keys ) {
         # Match resolution
         if ( $opts{'res'} ne "any" ) {
            if ( ($key eq "RESOLUTION") && ($atts{$key} ne $opts{'res'} ) ) {
               $set = undef;
            }
         }
         # Match any options from the build config file
         if ( $opts{'config'} ne "noconfig" ) {
            foreach my $optionkey ( keys(%bld) ) {
               if ( ($key eq "$optionkey") && ($atts{$key} ne $bld{$optionkey} ) ) {
                  $set  = undef;
               }
            }
         }
         # If ignore_ic_date option set and this is either: ic_ymd or ic_tod, continue
         if ( defined($settings{ignore_ic_date}) && (($key eq "ic_ymd") || 
              ($key eq "ic_tod")) ) {
             next;
         }

         # Match any options set from command line
         if ( $#settings_keys > -1 ) {
            foreach my $optionkey ( @settings_keys ) {
               # If ignore_ic_year option set and this is: ic_ymd continue
               if ( defined($settings{ignore_ic_year}) && ($key eq "ic_ymd") &&
                  ("ic_ymd" eq "$optionkey") ) {
                 my $ic_md     = $atts{$key};
                 my $set_ic_md = $settings{$optionkey};
                 my $ic_yr     = int($ic_md     / 10000);
                 my $set_ic_yr = int($set_ic_md / 10000);
                 $ic_md        = $ic_md     - $ic_yr    *10000;
                 $set_ic_md    = $set_ic_md - $set_ic_yr*10000;
                 if ( $ic_md ne $set_ic_md ) {
                    $set = undef;
                 }
                 next;
               }
               if ( ($key eq "$optionkey") && ($atts{$key} ne $settings{$optionkey} ) ) {
                  $set = undef;
               }
            }
         }
      }
    }
    my $isafile = 0;
    my $isadir = 0;
    #
    # If is a directory (has slashes, is csmdata or a var with dir in name)
    if ( $value =~ /\// && (($name eq "csmdata") || ($name =~ /dir/)) && 
         ($name ne "version") ) {
      if ( $name eq "csmdata" ) {
         $value = $$opts_ref{'csmdata'};
         $isadir = 1;
      } elsif ( $name eq "offline_atmdir" ) {
         $value = "$$opts_ref{'csmdata'}/$value";
         $isadir = 1;
      } else {
         $isadir = 1;
      }
    #
    # Otherwise if it is a file (has slashes in value)
    #
    } elsif ( $value =~ /\// && ($name ne "version") ) {
      $value = "$$opts_ref{'csmdata'}/$value";
      $isafile = 1;
    # For settings that are NOT files
    } else {
      $isafile = 0;
    }
    # Print out
    if ( defined($set) && (! exists($defaults{$name}{value}) ) ) { 
       $defaults{$name}{value}  = $value;
       $defaults{$name}{isfile} = $isafile;
       $defaults{$name}{isdir}  = $isadir;
    }
  }
  return( \%defaults );
}

1 # To make use or require happy
