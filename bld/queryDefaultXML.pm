#=======================================================================
#
#  This is a perl module to read in a DefaultNamelist XML file.
#
#=======================================================================
use strict;
use Build::Config;
use Build::NamelistDefinition;
use Build::NamelistDefaults;
use Build::Namelist;

package queryDefaultXML;

#-------------------------------------------------------------------------------

sub read_cfg_file
#
# Read in the configuration cache XML file on the build-time configuration
#
{
    my ($file, $empty_cfg_file, $printing, $settings_ref) = @_;

    my $cfg;
    my %config;
    if ( $file eq "noconfig" ) {
       print "No configuration cache file to read in.\n" if $printing;
       $cfg = Build::Config->new( $empty_cfg_file );
    } elsif ( -f "$file" ) {
       $cfg = Build::Config->new($file);
    } else {
       die "Bad filename entered: $file does NOT exist or can not open it.\n";
    }
    #
    # Make sure certain variables are set (in case this is a non-CLM config file)
    #
    if ( $cfg->get('spmd') eq "" ) { $config{'spmd'} = 0; }
    if ( $cfg->get('smp')  eq "" ) { $config{'smp'}  = 0; }
    if ( $cfg->is_valid_name('cam_bld') ) { $config{'mode'} = "ccsm_seq_cam"; }
    #
    # If cam config file -- set CLM specific things unknown by CAM
    #
    if ( $config{'mode'} eq "ccsm_seq_cam" ) {
       print "CAM configuration file\n" if $printing;
       my @deflist = ("dust", "voc", "ad_spinup", "exit_spinup", "supln" );
       foreach my $item ( @deflist ) {
          if ( ! defined($config{$item}) ) { $config{$item} = "off"; }
       }
       my $cppdefs = $config{'cppdefs'};
       my @cpp     = split( " ", $cppdefs );
       foreach my $cpp ( @cpp ) {
          if (      $cpp =~ /^[ ]*-D([^ ]+)[ ]*$/ ) {
            $cpp = $1;
          } elsif ( $cpp =~ /^[ ]*-U([^ ]+)[ ]*$/ ) {
            next;
          } else {
            die "$cppdefs in configuration cache file: $file is not in a readible format\n";
          }
          foreach my $item ( @deflist ) {
              if ( $cpp eq $item ) {
                $config{$item} = "on";
              }
          }
          if ( $cpp eq "cn" ) {
            $config{'bgc'} = $cpp;
          }
          if ( $cpp eq "casa" ) {
            $config{'bgc'} = $cpp;
          }
          if ( $cpp eq "DGVM" ) {
            $config{'bgc'} = $cpp;
          }
       }
    }
    foreach my $key ( keys( %config ) ) {
       if ( $cfg->is_valid_name( $key ) ) {
          $cfg->set( $key, $config{$key} );
       }
    }
    foreach my $key ( $cfg->get_names( ) ) {
       if ( defined($$settings_ref{$key}) ) {
          if ( $cfg->is_valid_name( $key ) ) {
             $cfg->set( $key, $$settings_ref{$key} );
          }
       }
    }
    return( $cfg );
}

#-------------------------------------------------------------------------------

sub ReadDefaultXMLFile {
#
# Read in the default XML file for the default namelist settings
#
  my $opts_ref     = shift;
  my $settings_ref = shift;

  # Error check that input and opts hash has the expected variables
  my $ProgName     = $$opts_ref{'ProgName'};
  my $nm = "${ProgName}::ReadDefaultXMLFile";
  my @required_list = ( "file", "nldef_file", "empty_cfg_file", "config", "namelist", 
                        "csmdata", "hgrid", "printing", "ProgName", "cmdline",      
                        "cfgdir"  );
  foreach my $var ( @required_list ) {
     if ( ! defined($$opts_ref{$var}) ) {
        die "ERROR($nm): Required input variable $var was not found\n";
     }
  }
  my $printing = $$opts_ref{'printing'};
  my $cmdline  = $$opts_ref{'cmdline'};
  # Initialize some local variables
  my $file               = $$opts_ref{'file'};
  my $nl_definition_file = $$opts_ref{'nldef_file'};
  my $empty_config_file  = $$opts_ref{'empty_cfg_file'};
  my $namelist           = $$opts_ref{'namelist'};

  my $cfg = read_cfg_file( $$opts_ref{'config'}, $$opts_ref{'empty_cfg_file'},
                           $printing, $settings_ref );

  #
  # Set up options to send to namelist defaults object
  #
  my %nlopts;
  foreach my $var ( keys( %$settings_ref) ) {
     $nlopts{$var} = $$settings_ref{$var};
  }
  if ( $$opts_ref{'hgrid'} ne "any" ) {
     $nlopts{'hgrid'} = $$opts_ref{'hgrid'};
  }
  #
  # Loop through all variables in file
  #
  print "($nm) Read: $file\n" if $printing;
  my %defaults;
  my $nldefaults = Build::NamelistDefaults->new($file, $cfg);
  my $definition = Build::NamelistDefinition->new($nl_definition_file);
  if ( $$opts_ref{'csmdata'} eq "default" ) {
    $$opts_ref{'csmdata'} = $nldefaults->get_value( "csmdata", \%nlopts );
  } 
  foreach my $name ( $nldefaults->get_variable_names() ) {
    my $value   = $nldefaults->get_value( $name, \%nlopts );
    my $isafile = 0;
    if ( $definition->is_input_pathname($name) ) {
       if ( $value && $name eq "furbinp" ) {
          $value   = $$opts_ref{'cfgdir'} . "/" . $value;
       } elsif ( $value ) {
          $value   = $$opts_ref{'csmdata'} . "/" . $value;
       }
       $isafile = 1;
    }
    my $isadir  = 0;
    my $isastr  = 0;
    if (  $definition->get_str_len($name) > 0 ) {
       $isastr = 1;
    }
    #
    # If is a directory (is a file and csmdata or a var with dir in name)
    #
    if ( $isafile  && (($name eq "csmdata") || ($name =~ /dir/)) ) {
      if ( $name eq "csmdata" ) {
         $value = $$opts_ref{'csmdata'};
         $isadir = 1;
      } else {
         $isadir = 1;
      }
    }
    # Return hash with the results
    my $group = $definition->get_group_name( $name );
    if ( $group eq $namelist && $value && (! exists($defaults{$name}{'value'}))  ) { 
       $defaults{$name}{'value'}  = $value;
       $defaults{$name}{'isfile'} = $isafile;
       $defaults{$name}{'isdir'}  = $isadir;
       $defaults{$name}{'isstr'}  = $isastr;
    }
  }
  return( \%defaults );
}

1 # To make use or require happy
