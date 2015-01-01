package test_do_transient_pfts;

# Unit tests for function: setup_logic_do_transient_pfts

use strict;
use Data::Dumper;
use Test::More;
use Test::Exception;

use CLMBuildNamelist qw(setup_logic_do_transient_pfts);

use parent qw(Test::Class);

#-------------------------------------------------------------------------------
#
# Common test fixture for all tests:
#
#-------------------------------------------------------------------------------
sub startup : Test(startup => 4) {
  my $self = shift;
  # provide common fixture for all tests, only created once at the
  # start of the tests.
  $self->{cfg} = Build::Config->new("t/input/config_cache_clm4_5_test.xml");
  isnt($self->{cfg}, undef, (caller(0))[3] . " : config object created.");

  $self->{definition} = Build::NamelistDefinition->new("t/input/namelist_definition_clm4_5_test.xml");
  isnt($self->{definition}, undef, (caller(0))[3] . " : namelist_definition object created.");

  $self->{defaults} = Build::NamelistDefaults->new("t/input/namelist_defaults_clm4_5_test.xml");
  isnt($self->{defaults}, undef,  (caller(0))[3] . " : namelist_defaults object created.");

  $self->{physv} = config_files::clm_phys_vers->new( $self->{cfg}->get('phys') );
  isnt($self->{physv}, undef,  (caller(0))[3] . " : phys_vers object created.");

  $self->{test_files} = 0;
  $self->{nl_flags} = {};
}

sub shutdown : Test(shutdown) {
  # cleanup the single instance test fixtures
}

sub setup : Test(setup => 1) {
  my $self = shift;
  # provide common fixture for all tests, create fresh for each test

  $self->{nl} = Build::Namelist->new();
  isnt($self->{nl}, undef, (caller(0))[3] . " : empty namelist object created.");

  # Set use_cndv and use_ed so that they don't conflict with do_transient_pfts
  $self->set_value('use_cndv', '.false.');
  $self->set_value('use_ed', '.false.');
}

sub teardown : Test(teardown) {
  # clean up after test
}

#-------------------------------------------------------------------------------
#
# Other common routines
#
#-------------------------------------------------------------------------------

sub set_value {
   # Set the value of a namelist option
   my ($self, $var, $value) = @_;
   
   my $group = $self->{definition}->get_group_name($var);
   $self->{nl}->set_variable_value($group, $var, $value);
}

sub set_nontransient {
   # Set up flanduse_timeseries for a non-transient case
   my $self = shift;

   $self->set_value('flanduse_timeseries', ' ');
}

sub set_transient {
   # Set up flanduse_timeseries for a transient case
   my $self = shift;

   $self->set_value('flanduse_timeseries', 'foo.nc');
}

sub get_do_transient_pfts {
   my $self = shift;

   return $self->{nl}->get_value("do_transient_pfts");
}

sub setup_logic_do_transient_pfts {
   my $self = shift;

   CLMBuildNamelist::setup_logic_do_transient_pfts($self->{test_files}, $self->{nl_flags}, $self->{definition}, $self->{defaults}, $self->{nl}, $self->{physv});
}

#-------------------------------------------------------------------------------
#
# tests
#
#-------------------------------------------------------------------------------

sub test_do_transient_pfts__default_transient : Tests {
   my $self = shift;

   my $msg = "Test default value for do_transient_pfts in a transient case.\n";

   $self->set_transient;
   $self->setup_logic_do_transient_pfts;
   my $result = $self->get_do_transient_pfts;
   is($result, '.true.') || diag($msg);
}

sub test_do_transient_pfts__default_nontransient : Tests {
   my $self = shift;

   my $msg = "Test default value for do_transient_pfts in a non-transient case.\n";

   $self->set_nontransient;
   $self->setup_logic_do_transient_pfts;
   my $result = $self->get_do_transient_pfts;
   is($result, undef) || diag($msg);
}

sub test_do_transient_pfts__default_cndv : Tests {
   my $self = shift;

   my $msg = "Test default value for do_transient_pfts in a cndv case.\n";

   $self->set_transient;
   $self->set_value('use_cndv', '.true.');
   $self->setup_logic_do_transient_pfts;
   my $result = $self->get_do_transient_pfts;
   is($result, undef) || diag($msg);
}

sub test_do_transient_pfts__default_ed : Tests {
   my $self = shift;

   my $msg = "Test default value for do_transient_pfts in an ED case.\n";

   $self->set_transient;
   $self->set_value('use_ed', '.true.');
   $self->setup_logic_do_transient_pfts;
   my $result = $self->get_do_transient_pfts;
   is($result, undef) || diag($msg);
}

sub test_do_transient_pfts__override_default : Tests {
   my $self = shift;

   my $msg = "Test ability to set value to false when the default is true.\n";
   
   $self->set_transient;
   $self->set_value('do_transient_pfts', '.false.');
   $self->setup_logic_do_transient_pfts;
   my $result = $self->get_do_transient_pfts;
   is($result, '.false.') || diag($msg);
}

sub test_do_transient_pfts__override_default_not_allowed : Tests {
   my $self = shift;

   my $msg = "Make sure overriding the default is not allowed when running a non-transient case.\n";
   
   $self->set_nontransient;
   $self->set_value('do_transient_pfts', '.true.');
   dies_ok(sub {$self->setup_logic_do_transient_pfts}) || diag($msg);
}


1;
