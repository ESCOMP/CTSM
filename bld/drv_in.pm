#
#	drv_in.pm			Erik Kluzek
#
#	Perl module to deal with the "virtual" drv_in namelist.
#
#------------------------------------------------------------------------
#
#	Description of methods:
#
#       new ----------------------- Constructor.
#       set_output_values --------- Move namelist values to real namelists.
#       get_descrips -------------- Get descriptions of real namelists.
#
#	$Id$
#
use strict;
#use diagnostics;
use Cwd;
use XML::Lite;

package drv_in;

use namelist;
@drv_in::ISA = 'namelist';
%drv_in::NL;
#
# Extend the namelist class to figure out which namelist the control driver input
# information goes.
#
sub new {
#
# Constructor
#
  my $class      = shift;
  my $optsref    = shift;
  my $config     = shift;
  my $LNDNLref   = shift;
  my $PROFNLref  = shift;
  my $CSMNLref   = shift;
  my $TIMNLref   = shift;
  my $ATMSHRNLref= shift;

  my %bld = &queryDefaultXML::read_cfg_file( $config );
  my $deffile;
  my $MODE = $bld{'MODE'};
  my %deffiles = ( ccsm_seq     => "SeqCCSMDrvInNamelistsDescriptions.xml",
                   ccsm_seq_cam => "CAMSeqCCSMDrvInNamelistsDescriptions.xml",
                   offline      => "CLMDrvInNamelistsDescriptions.xml",
                   ext_ccsm_con => "CCSMDrvInNamelistsDescriptions.xml"
                 );
  my $deffile = $deffiles{$MODE};
  my $interactive = 0;
  my $file        = "do_not_write_file";
  my $printlev    = $$optsref{'printlev'};
  my $cfgdir      = $$optsref{'cfgdir'};
  my $self = $class->SUPER::new( "drv_in", $file,  \%drv_in::NL, $printlev );
  $self->{'descript_name'} = $MODE;
  $self->{'descript_file'} = "$cfgdir/$deffile";
  $self->{'cfgdir'}        = $cfgdir;
  $self->{'config'}        = $config;

  if ( ref($PROFNLref) ne "prof_inparm" ) {
    die "ERROR($class): Object sent to drv_in constructor not a prof_inparm object\n";
  }
  if ( ref($LNDNLref) ne "clm_inparm" ) {
    die "ERROR($class): Object sent to drv_in constructor not a clm_inparm object\n";
  }
  $self->{'PROFNL'}       = $PROFNLref;
  $self->{'LNDNL'}        = $LNDNLref;
  if ( $MODE eq "ccsm_seq" ) {
     if ( ref($CSMNLref) ne "ccsm_inparm" ) {
       die "ERROR($class): Object sent to drv_in constructor not a ccsm_inparm object\n";
     }
     if ( ref($TIMNLref) ne "timemgr_inparm" ) {
       die "ERROR($class): Object sent to drv_in constructor not a timemgr_inparm object\n";
     }
     if ( ref($ATMSHRNLref) ne "datm_dshr_in" ) {
       die "ERROR($class): Object sent to drv_in constructor not a datm_dshr_in object\n";
     }
     $self->{'CSMNL'}        = $CSMNLref;
     $self->{'TIMNL'}        = $TIMNLref;
     $self->{'ATMSHRNL'}     = $ATMSHRNLref;
  }

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub set_output_values {
#
#  set output values
#
  my $self  = shift;

  my $class      = ref($self);
  my $NLref          = $self->{'NLREF'};
  # Read in namelist description file
  $self->get_descrips( );
  #
  # Loop through each of the namelist items
  # Move any namelist items from DRV_IN namelist to given namelist
  #
  my %namelists = ( ccsm_inparm=>"CSMNL",  timemgr_inparm=>"TIMNL", 
                    prof_inparm=>"PROFNL", clm_inparm    =>"LNDNL",
                    datm_dshr_in=>"ATMSHRNL" );
  my $move_namelists = "";
  foreach my $namel ( keys( %namelists ) ) {
     my $Listref  = $self->{$namel};
     if ( defined($Listref) )  {
        my $move_namelists .= " $namel";
        my $NLOBJref = $namelists{$namel};
        my $NAMELref = $self->{$NLOBJref};
        foreach my $key ( @$Listref ) {
           my $value = $NAMELref->Value($key);
           if ( ! defined($value) && defined($$NLref{$key}) ) {
              $NAMELref->SetValue($key,$$NLref{$key});
              $$NLref{$key} = undef;
           }
        }
     }
  }
  #
  # Move whatever is leftover to the clm_inparm namelist
  #
  my $namel    = "clm_inparm";
  my $NLOBJ    = $namelists{$namel};
  my $NLCLMref = $self->{$NLOBJ};
  foreach my $key ( keys(%$NLref) ) {
     my $value = $NLCLMref->Value($key);
     if ( ! defined($value) && defined($$NLref{$key}) ) {
        print << "EOF";

Any namelist items in the drv_in 'virtual' namelist should be moved over to a
real driver namelist used in this case. If there are any namelist items left over
in the drv_in namelist -- this indicates an error.

EOF
        if ( $move_namelists ne "" ) {
           print "drv_in namelist items moved to namelists: $move_namelists\n";
        }
        die "ERROR($class):: some drv_in namelist items are leftover: $key\n";
     }
  }

}
  
#============================================================================

sub get_descrips {
# 
#  Parse the description XML file that describes which namelist has what variables
#  
  my $self  = shift;

  my $class = ref($self);
  my $nm = "$class\:\:get_descrips";
  
  my $optsref = $self->{'optsref'};
  my $file    = $self->{'descript_file'};
  print "($nm) Read: $file\n" if ($self->{'printlev'}>2);
  my $xml = XML::Lite->new( $file );
  if ( ! defined($xml) ) {
    die "ERROR($nm): Trouble opening or reading file: $file\n";
  }
  #
  # Find the main name element
  # 
  my $elm = $xml->root_element( );
  my $main = $self->{'descript_name'};
  my @list = $xml->elements_by_name( $main );
  if ( $#list < 0 ) {
    die "ERROR($nm): could not find the main $main element in $file\n";
  }
  if ( $#list != 0 ) {
    die "ERROR($nm): $main main element in $file is duplicated, there should only be
one\n";
  }  $elm = $list[0];
  my @children = $elm->get_children();
  if ( $#children < 0 ) {
    die "ERROR($nm): There are no sub-elements to the $main element in $file\n";
  }
  #
  # Get the names for each element in the list into a hash
  #
  my %names;
  foreach my $child ( @children ) {
    my $name = $child->get_name();
    $names{$name} = $child->get_text();
  }
  #
  # Go through the sub-elements to the main name element
  #
  foreach my $child ( @children ) {
    #
    # Get the attributes for each namelist which are the namelist elements
    #
    my %atts = $child->get_attributes;
    # Name of element, and it's associated value
    my $name = $child->get_name();
    my $value =  $child->get_text();
    $value =~ s/\n//g;   # Get rid of extra returns 
    my @keys = keys(%atts);
    my $match = 1;
    if ( $#keys >= 0 ) {
       $self->{$name} = \@keys;
    }
  }
}

#============================================================================

1   # to make use or require happy
