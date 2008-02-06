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
  character(len=256) :: hostname                ! Hostname of machine running on
  character(len=256) :: username                ! username of user running program
  character(len=256) :: version                 ! version of program
  character(len=256) :: source                  ! description of this source
  character(len=256) :: conventions             ! dataset conventions

!
! Unit Numbers
!
  integer :: iulog = 6        ! "stdout" log file unit number, default is 6
!
! Initial file variables
!
  character(len= 8) :: hist_crtinic             ! if set to 'MONTHLY' or 'YEARLY', write initial cond. file
!
! Output NetCDF files
!
  logical :: outnc_large_files                  ! large file support for output NetCDF files
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
! offline atmosphere data cycling controls
!
  integer :: cycle_begyr                        ! first year of offline atm data (e.g. 1948)
  integer :: cycle_nyr                          ! number of years of offline atm data to cycle
  
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
  real(r8) :: co2_ppmv                          ! atmospheric CO2 molar ratio (by volume) (umol/mol)

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
!
! Error growth perturbation limit
!
  real(r8) :: pertlim                           ! perturbation limit when doing error growth test


!
! History File control
!
  logical :: hist_pioflag    ! turns on and off hist with pio
  logical :: ncd_lowmem2d    ! turns on low memory 2d writes in clm hist
  logical :: ncd_pio_def     ! default pio use setting
  logical :: ncd_pio_UseRearranger  ! use MCT or box
  logical :: ncd_pio_UseBoxRearr    ! use box
  logical :: ncd_pio_SerialCDF      ! write with pio serial netcdf mode
  logical :: ncd_pio_IODOF_rootonly ! write history in pio from root only
  integer :: ncd_pio_DebugLevel     ! pio debug level
  integer :: ncd_pio_num_iotasks    ! num of iotasks to use

!
! !REVISION HISTORY:
! Created by Mariana Vertenstein and Gordon Bonan
! 1 June 2004, Peter Thornton: added fnedpdat for nitrogen deposition data
!
!EOP
!-----------------------------------------------------------------------

end module clm_varctl
