#!/usr/bin/env perl
# -*- mode: cperl -*-
#-----------------------------------------------------------------------------------------------

require 5;

use strict;
use warnings;
use diagnostics;

use Test::More;

BEGIN {
  #
  # setup paths for cesm and local utility modules at compile time.
  #
  # Assumes that we are running from clm/bld/test....
  #
  use Cwd qw(getcwd abs_path);

  my $cwd        = getcwd();
  my $cesmroot   = abs_path("../../../../");

  my @dirs = ("../", 
	      "$cwd/perl5lib",
	      "$cesmroot/cime/utils/perl5lib",
	      "$cesmroot/cime/scripts/Tools");

  unshift @INC, @dirs;

  # check for a couple of modules from the paths we added.
  use_ok("CLMBuildNamelist");
  use_ok("Config::SetupTools");
  use_ok("Test::Class");
}


# ccsm perl modules
require Build::Config;
require Build::NamelistDefinition;
require Build::NamelistDefaults;
require Build::Namelist;
require config_files::clm_phys_vers;
require Config::SetupTools;

# local perl modules
use Test::Class;
use Test::Exception;

require CLMBuildNamelist;

use Test::Class::Load "./t";

Test::Class->runtests;

