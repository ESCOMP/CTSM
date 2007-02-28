#include <misc.h>
#include <preproc.h>

module clm_varctl

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varctl
!
! !DESCRIPTION:
! Module containing run control variables
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Run control variables
!
  character(len=256) :: caseid                  ! case id
  character(len=256) :: ctitle                  ! case title
  integer :: nsrest                             ! 0: initial run. 1: restart: 3: branch
  logical, public :: brnch_retain_casename = .false. ! true => allow case name to remain the same for branch run
                                                     ! by default this is not allowed
!
! Initial file variables
!
  character(len= 8) :: hist_crtinic             ! if set to 'MONTHLY' or 'YEARLY', write initial cond. file
!
! Long term archive variables
!
  character(len=256) :: archive_dir             ! long term archive directory (can be mass store)
  character(len=  8) :: mss_wpass               ! mass store write password for output files
  integer            :: mss_irt                 ! mass store retention period
!
! Run input files
!
  character(len=256) :: finidat                 ! initial conditions file name
  character(len=256) :: fsurdat                 ! surface data file name
  character(len=256) :: fatmgrid                ! atm grid file name
  character(len=256) :: fatmlndfrc              ! lnd frac file on atm grid
  character(len=256) :: fatmtopo                ! topography on atm grid
  character(len=256) :: flndtopo                ! topography on lnd grid
  character(len=256) :: fndepdat                ! static nitrogen deposition data file name
  character(len=256) :: fndepdyn                ! dynamic nitrogen deposition data file name
  character(len=256) :: fpftdyn                 ! dynamic landuse dataset
  character(len=256) :: fpftcon                 ! ASCII data file with PFT physiological constants
  character(len=256) :: nrevsn                  ! restart data file name for branch run
  character(len=256) :: frivinp_rtm             ! RTM input data file name
  character(len=256) :: offline_atmdir          ! directory for input offline model atm data forcing files (Mass Store ok)
!
! Landunit logic
!
  logical :: create_crop_landunit               ! true => separate crop landunit is not created by default
  logical :: allocate_all_vegpfts               ! true => allocate memory for all possible vegetated pfts on
                                                ! vegetated landunit if at least one pft has nonzero weight
!
! BGC logic
!
  character(len=16) :: co2_type                 ! values of 'prognostic','diagnostic','constant'
!
! Physics
!
  integer :: irad                               ! solar radiation frequency (iterations)
  logical :: wrtdia                             ! true => write global average diagnostics to std out
  logical :: csm_doflxave                       ! true => only communicate with flux coupler on albedo calc time steps
!
! single column control variables
!
  logical :: single_column                      ! true => single column mode
  real(r8):: scmlat			        ! single column lat
  real(r8):: scmlon			        ! single column lon
!
! Rtm control variables
!
  integer :: rtm_nsteps                         ! if > 1, average rtm over rtm_nsteps time steps
!
! Decomp control variables
!
  integer :: nsegspc                            ! number of segments per clump for decomp
!
! Derived variables (run, history and restart file)
!
  character(len=256) :: rpntdir                 ! directory name for local restart pointer file
  character(len=256) :: rpntfil                 ! file name for local restart pointer file
  character(len=256) :: version                 ! model version number
!
! Error growth perturbation limit
!
  real(r8) :: pertlim                           ! perturbation limit when doing error growth test
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein and Gordon Bonan
! 1 June 2004, Peter Thornton: added fnedpdat for nitrogen deposition data
!
!EOP
!-----------------------------------------------------------------------

end module clm_varctl
