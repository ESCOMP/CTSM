package config_files::clm_phys_vers;
my $pkg_nm = 'config_files::clm_phys_vers';
#-----------------------------------------------------------------------------------------------
#
# SYNOPSIS
#
# require config_files::clm_phys_vers;
#
# my $phys = config_files::clm_phys_vers->new("clm5_0");
# print $phys->as_string();
#
# DESCRIPTION
#
# Enter the physics version as a string, with a list of valid versions, and have the ability to convert it to 
# different formats.
#
# COLLABORATORS: None
# 
#-----------------------------------------------------------------------------------------------
#
# Date        Author                  Modification
# 03/06/2014  Erik Kluzek             creation
#
#--------------------------------------------------------------------------------------------

use strict;
use bigint;
#use warnings;
#use diagnostics;

my @version_strings = ("clm4_5", "clm5_0");

#-------------------------------------------------------------------------------

sub new {
    # Constructor, enter version string as argument
    my $class       = shift;
    my $vers_string = shift;

    my $nm   = "$class\:\:new";
    my $self = {};
    bless($self, $class);
    $self->__validate_vers__( $vers_string );
    $self->{'vers_string'} = $vers_string;
    return( $self );
}

#-------------------------------------------------------------------------------

sub __validate_vers__ {
   # Make sure the version string is a valid one
   my $class       = shift;
   my $vers_string = shift;

   my $found = undef;
   foreach  my $i (0..$#version_strings) {
      if ( $vers_string eq $version_strings[$i] ) {
         $found = 1;
         last;
      }
   }
   if ( ! defined($found) ) {
      die "NOT a valid CLM version: $vers_string\n";
   }
}

#-------------------------------------------------------------------------------

sub as_string {
# Return the physics version as a string
  my $self = shift;

  my $phys = $self->{'vers_string'};
  return( $phys );
}

#-----------------------------------------------------------------------------------------------
# Unit testing of above
#-----------------------------------------------------------------------------------------------
if ( ! defined(caller) && $#ARGV == -1 ) {
   package phys_vers_unit_tester;

   require Test::More;
   Test::More->import( );

   plan( tests=>2 );

   sub testit {
      print "unit tester\n";
      my %lastv;
      my @vers_list = ( "clm4_5", "clm5_0" );
      foreach my $vers ( @vers_list ) {
         my $phys = config_files::clm_phys_vers->new($vers);
         isa_ok($phys, "config_files::clm_phys_vers", "created clm_phys_vers object");
         print "$vers: string: ".$phys->as_string()."\n";
      }
   }
}

#-----------------------------------------------------------------------------------------------
# Determine if you should run the unit test or if this is being called from a require statement
#-----------------------------------------------------------------------------------------------

if ( defined(caller) ) {
   1   # to make use or require happy
} elsif ( $#ARGV == -1 ) {
   &phys_vers_unit_tester::testit();
}
