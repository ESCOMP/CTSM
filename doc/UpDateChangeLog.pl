#!/usr/bin/env perl
#=======================================================================
#
#  This is a script to update the ChangeLog
#
# Usage:
#
# perl ChangeLog tag-name One-line summary
#
#
#=======================================================================

use strict;
#use warnings;
#use diagnostics;

use English;

my $ProgName;
($ProgName = $PROGRAM_NAME) =~ s!(.*)/!!; # name of program
my $ProgDir = $1;                         # name of directory where program lives

sub usage {
    die <<EOF;
SYNOPSIS
     $ProgName <tag-name> <one-line-summary>
ARGUMENTS
     <tag-name>             Tag name of tag to document
     <one-line-summary>     Short summary description of this tag
EOF
}

if ( $#ARGV != 1 ) {
  print "ERROR: wrong number of arguments: $ARGV\n";
  usage();
}

my $tag = $ARGV[0];
my $sum = $ARGV[1];

if ( $tag !~ /clm[0-9]+_(expa|[0-9]+)_[0-9]+/ ) {
  print "ERROR: bad tagname: $tag\n";
  usage();
}
my $EDITOR = $ENV{EDITOR};
if ( $EDITOR =~ "" ) {
  print "ERROR: editor NOT set -- set the env variable EDITOR to the text editor you would like to use\n";
  usage();
}


my $template      = "ChangeLog_template";
my $changelog     = "ChangeLog";
my $changelog_tmp = "ChangeLog.tmp";

my $user = $ENV{USER};
if ( $user !~ /.+/ )  {
  die "ERROR: Could not get user name: $user";
}
my @list = getpwnam( $user );
my $fullname = $list[6];
my $date = `date`;
chomp( $date );

if ( $date !~ /.+/ )  {
  die "ERROR: Could not get date: $date\n";
}

open( TL, "<$template"     )  || die "ERROR:: trouble opening file: $template";
open( CL, "<$changelog"     ) || die "ERROR:: trouble opening file: $changelog";
open( FH, ">$changelog_tmp" ) || die "ERROR:: trouble opening file: $changelog_tmp";

while( $_ = <TL> ) {
  if (      $_ =~ /Tag name:/ ) {
     chomp( $_ );
     print FH "$_ $tag\n";
  } elsif ( $_ =~ /Originator/ ) {
     chomp( $_ );
     print FH "$_ $user ($fullname)\n";
  } elsif ( $_ =~ /Date:/ ) {
     chomp( $_ );
     print FH "$_ $date\n";
  } elsif ( $_ =~ /One-line Summary:/ ) {
     chomp( $_ );
     print FH "$_ $sum\n";
  } else {
     print FH $_;
  }
}
while( $_ = <CL> ) {
  print FH $_;
}
close( TL );
close( CL );
close( FH );
system( "/bin/mv $changelog_tmp $changelog" );
system( "$EDITOR $changelog" );
