package namelist_files::LogMessages;
my $pkg_nm = 'namelist_files::LogMessages';
#-----------------------------------------------------------------------------------------------
#
# SYNOPSIS
#
# require namelist_files::LogMessages;
#
# my %opts;
# my $log = namelist_files::LogMessages->new("ProgName", \%opts);
# $log->message("message to print");
# $log->verbose_message("message to print only if verbose mode is on");
# $log->warning("Warning message");
# $log->exit_message("clean exit");
# $log->fatal_error("die with fatal error");
# $log->final_exit("Final message to send (and exit");
#
#
# DESCRIPTION
#
# Handles log messages for perl. Sets up log messages according to verbose
# or silent setting. It also handles warnings printing them, but on finalization
# aborting unless ignore_warnings was set.
#
# COLLABORATORS: None
# 
#-----------------------------------------------------------------------------------------------
#
# Date        Author                  Modification
# 10/06/2017  Erik Kluzek             creation
#
#--------------------------------------------------------------------------------------------

use strict;
#use warnings;
#use diagnostics;

#-------------------------------------------------------------------------------

sub new {
    my $class    = shift;
    my $ProgName = shift;
    my %opts     = %{shift()};

    my $nm   = "$class\:\:new";
    my $self = {};
    bless($self, $class);
    $self->{'nwarns'}    = 0;
    $self->{'verbosity'} = 1;
    $self->{'NO_EXIT'}   = $opts{'NO_EXIT'};
    $self->{'ProgName'}  = $ProgName;
    $self->{'ignore_warnings'} = $opts{'ignore_warnings'};
    $self->__set_print_level( \%opts );
    return( $self );
}


#-------------------------------------------------------------------------------

sub __set_print_level {
  my $self        = shift;
  # Define print levels:
  # 0 - only issue fatal error messages
  # 1 - only informs what files are created (default)
  # 2 - verbose
  my %opts = %{shift()};

  if ( $opts{'silent'} && $opts{'verbose'} ) {
     $self->fatal_error( "Can not set both the -silent and the -verbose options -- set one or the other" );
  }
  my $verbosity = 1;
  if ($opts{'silent'})  { $verbosity = 0; }
  if ($opts{'verbose'}) { $verbosity = 2; }
  $self->{'verbosity'} = $verbosity;
  $self->{'print_verbose'} = 2;
}

#-------------------------------------------------------------------------------

sub message {
  my $self = shift;
  my ($message) = @_;
  if ($self->{'verbosity'} > 0) {
    print "$message\n";
  }
}

#-------------------------------------------------------------------------------

sub verbose_message {
  my $self = shift;

  my ($message) = @_;
  if ($self->{'verbosity'} >= $self->{'print_verbose'}) {
    print "$message\n";
  }
}
#-------------------------------------------------------------------------------

sub nwarns {
  my $self = shift;

  return( $self->{'nwarns'} );
}

#-------------------------------------------------------------------------------

sub final_exit {
  my $self = shift;
  my ($message) = @_;
  if ( $self->{'nwarns'} > 0 ) {
    $self->message( "\n\nYou ran with the -ignore_warnings options and allowed $self->{'nwarns'} to go past\n" );
  }
  $self->verbose_message( $message );
  if ( $self->{'NO_EXIT'} ) {
    die
  } else {
    exit;
  }
}

#-------------------------------------------------------------------------------
# Some simple subroutines to do a clean exit, print warning, or a fatal error

sub exit_message {
  my $self = shift;
  my ($message) = @_;
  print "$self->{ProgName} : $message\n";
  if ( $self->{'NO_EXIT'} ) {
    die
  } else {
    exit;
  }
}

#-------------------------------------------------------------------------------

sub warning {
  my $self      = shift;
  my $message   = shift;

  $self->{'nwarns'} = $self->{'nwarns'} + 1;
  my $func_name = (caller(1))[3];
  if ( $self->{'ignore_warnings'} ) {
    print "Warning : $self->{ProgName}::${func_name}() : $message\n\n";
  } else {
    die "Warning : $self->{ProgName}::${func_name}() : $message\n" . 
        " -- Add -ignore_warnings option to CLM_BLDNML_OPTS to ignore this warning\n\n";
  }
}

#-------------------------------------------------------------------------------

sub fatal_error {
  my $self        = shift;
  my ($message) = @_;
  my $func_name = (caller(1))[3];
  die "ERROR : $self->{ProgName}::${func_name}() : $message\n";
}

#-------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
# Unit testing of above
#-----------------------------------------------------------------------------------------------
if ( ! defined(caller) && $#ARGV == -1 ) {
   package LogMessage_unit_tester;

   require Test::More;
   Test::More->import( );

   plan( tests=>11 );

   sub testit {
      print "unit tester\n";
      my %opts;
      my $message;

      # Standard verbose level, test all methods
      $opts{'NO_EXIT'}         = 1;
      my $log = namelist_files::LogMessages->new("ProgName", \%opts);
      isa_ok($log, "namelist_files::LogMessages", "Created LogMessages object");
      $log->message("message to print");
      $log->verbose_message("YOU SHOULD NOT SEE THIS MESSAGE BECAUSE IT IS VERBOSE AND VERBOSE NOT ON");
      $message = "Warning message";
      is ( $log->nwarns(), 0, "Make sure have zero warnings" );
      eval{ $log->warning($message); };
      like( $@, qr/$message/, "check that a warning dies without ignore_warnings option" );
      is ( $log->nwarns(), 1, "Make sure have one warning" );
      $message = "die with fatal error";
      eval{ $log->fatal_error($message); };
      like( $@, qr/$message/, "check that a fatal_error dies" );
      $message = "exit with exit message";
      eval{ $log->exit_message($message); };
      like( $@, qr/Died/, "check that a exit_message exits" );
      $message = "Final message to send";
      eval{ $log->final_exit($message); };
      like( $@, qr/Died/, "check that a final exits" );

      # Test ignore_warnings option and verbose mode
      $opts{'ignore_warnings'} = 1;
      $opts{'verbose'}         = 1;
      $opts{'NO_EXIT'}         = 1;
      $log = namelist_files::LogMessages->new("ProgName", \%opts);
      isa_ok($log, "namelist_files::LogMessages", "Created LogMessages object");
      $log->verbose_message("message to print only if verbose mode is on");
      $log->warning("Warning message");
      $log->warning("Warning message2");
      $log->warning("Warning message3");
      $log->warning("Warning message4");
      $log->warning("Warning message5");
      is ( $log->nwarns(), 5, "Make sure have five warnings" );
      eval{ $log->final_exit($message); };
      print "content: $@\n";
      like( $@, qr/Died/, "check that a final_exit with warning exits" );
      # silent mode
      $opts{'ignore_warnings'} = 0;
      $opts{'verbose'}         = 0;
      $opts{'silent'}          = 1;
      $opts{'NO_EXIT'}         = 1;
      $log = namelist_files::LogMessages->new("ProgName", \%opts);
      $log->message("YOU SHOULD NOT SEE THIS MESSAGE BECAUSE SILENT MODE IS ON");
      $log->verbose_message("YOU SHOULD NOT SEE THIS VERBOSE MESSAGE BECAUSE SILENT MODE IS ON");
      # Should die with error if both silent and verbose mode is on
      $opts{'ignore_warnings'} = 0;
      $opts{'verbose'}         = 1;
      $opts{'silent'}          = 1;
      $opts{'NO_EXIT'}         = 1;
      eval{ $log = namelist_files::LogMessages->new("ProgName", \%opts); };
      print "content: $@\n";
      like( $@, qr/ERROR : /, "check that died if both verbose and silent mode is on" );
      print "\nSuccessfully ran all tests\n";
   }
}

#-----------------------------------------------------------------------------------------------
# Determine if you should run the unit test or if this is being called from a require statement
#-----------------------------------------------------------------------------------------------

if ( defined(caller) ) {
   1   # to make use or require happy
} elsif ( $#ARGV == -1 ) {
   &LogMessage_unit_tester::testit();
}
