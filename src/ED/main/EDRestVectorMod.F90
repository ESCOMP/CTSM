module EDRestVectorMod

#include "shr_assert.h"

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_sys_mod     , only : shr_sys_abort
  use clm_varctl      , only : iulog
  use spmdMod         , only : masterproc
  use decompMod       , only : bounds_type
  use pftconMod       , only : pftcon
  use EDTypesMod      , only : area, cohorts_per_col, numpft_ed, numWaterMem, cp_nclmax, numCohortsPerPatch
  use EDTypesMod      , only : ncwd, invalidValue, cp_nlevcan
  use EDTypesMod      , only : ed_site_type, ed_patch_type, ed_cohort_type
  use abortutils      , only : endrun

  !
  implicit none
  private
  !
  ! integer constants for storing logical data
  integer, parameter :: old_cohort = 0
  integer, parameter :: new_cohort = 1  
  !
  ! ED cohort data as a type of vectors
  !
  type, public :: EDRestartVectorClass
     !
     ! for vector start and stop, equivalent to begCohort and endCohort 
     !
     integer   :: vectorLengthStart
     integer   :: vectorLengthStop

     logical   ::  DEBUG = .false.
     !
     ! add ED vectors that need to be written for Restarts
     !

     ! required to map cohorts and patches to/fro
     ! vectors/LinkedLists
     integer,  pointer :: numPatchesPerCol(:)
     integer,  pointer :: cohortsPerPatch(:)
     !
     ! cohort data
     !
     real(r8), pointer :: balive(:)
     real(r8), pointer :: bdead(:) 
     real(r8), pointer :: bl(:) 
     real(r8), pointer :: br(:) 
     real(r8), pointer :: bstore(:) 
     real(r8), pointer :: canopy_layer(:) 
     real(r8), pointer :: canopy_trim(:) 
     real(r8), pointer :: dbh(:) 
     real(r8), pointer :: hite(:) 
     real(r8), pointer :: laimemory(:) 
     real(r8), pointer :: leaf_md(:)  ! this can probably be removed
     real(r8), pointer :: root_md(:)  ! this can probably be removed
     real(r8), pointer :: n(:) 
     real(r8), pointer :: gpp_acc(:) 
     real(r8), pointer :: npp_acc(:) 
     real(r8), pointer :: gpp(:) 
     real(r8), pointer :: npp(:) 
     real(r8), pointer :: npp_leaf(:) 
     real(r8), pointer :: npp_froot(:) 
     real(r8), pointer :: npp_bsw(:) 
     real(r8), pointer :: npp_bdead(:) 
     real(r8), pointer :: npp_bseed(:) 
     real(r8), pointer :: npp_store(:) 
     real(r8), pointer :: bmort(:) 
     real(r8), pointer :: hmort(:) 
     real(r8), pointer :: cmort(:) 
     real(r8), pointer :: imort(:) 
     real(r8), pointer :: fmort(:) 
     real(r8), pointer :: ddbhdt(:) 
     real(r8), pointer :: resp_tstep(:) 
     integer,  pointer :: pft(:) 
     integer,  pointer :: status_coh(:)
     integer,  pointer :: isnew(:)
     !
     ! patch level restart vars
     ! indexed by ncwd
     !
     real(r8), pointer :: cwd_ag(:) 
     real(r8), pointer :: cwd_bg(:)
     !
     ! indexed by pft
     !
     real(r8), pointer :: leaf_litter(:)
     real(r8), pointer :: root_litter(:)
     real(r8), pointer :: leaf_litter_in(:)
     real(r8), pointer :: root_litter_in(:)
     real(r8), pointer :: seed_bank(:)
     !
     ! indext by nclmax
     !
     real(r8), pointer :: spread(:)
     !
     ! one per patch
     !
     real(r8), pointer :: livegrass(:)   ! this can probably be removed
     real(r8), pointer :: age(:)
     real(r8), pointer :: areaRestart(:)
     real(r8), pointer :: f_sun(:)
     real(r8), pointer :: fabd_sun_z(:)
     real(r8), pointer :: fabi_sun_z(:)
     real(r8), pointer :: fabd_sha_z(:)
     real(r8), pointer :: fabi_sha_z(:)
     !
     ! site level restart vars
     !
     real(r8), pointer :: water_memory(:) 
     real(r8), pointer :: old_stock(:) 
     real(r8), pointer :: cd_status(:) 
     real(r8), pointer :: dd_status(:)
     real(r8), pointer :: ED_GDD_site(:)
     real(r8), pointer :: ncd(:)   
     real(r8), pointer :: leafondate(:)   
     real(r8), pointer :: leafoffdate(:)   
     real(r8), pointer :: dleafondate(:)   
     real(r8), pointer :: dleafoffdate(:) 
     real(r8), pointer :: acc_NI(:) 
             
     
   contains
     !
     ! implement getVector and setVector
     !
     procedure :: setVectors    
     procedure :: getVectors
     !
     ! restart calls 
     !
     procedure :: doVectorIO
     !
     ! clean up pointer arrays
     !
     procedure :: deleteEDRestartVectorClass
     !
     ! utility routines
     !
     procedure :: convertCohortListToVector
     procedure :: createPatchCohortStructure
     procedure :: convertCohortVectorToList
     procedure :: printIoInfoLL
     procedure :: printDataInfoLL
     procedure :: printDataInfoVector

  end type EDRestartVectorClass

  ! Fortran way of getting a user-defined ctor
  interface EDRestartVectorClass
     module procedure newEDRestartVectorClass
  end interface EDRestartVectorClass

  ! 
  ! non type-bound procedures
  !
  public :: EDRest

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-------------------------------------------------------------------------------!

contains

   !--------------------------------------------!
   ! Type-Bound Procedures Here:
   !--------------------------------------------!

  !-------------------------------------------------------------------------------!
  subroutine deleteEDRestartVectorClass( this )
    !
    ! !DESCRIPTION:
    ! provide clean-up routine of allocated pointer arrays
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(EDRestartVectorClass), intent(inout)      :: this
    !
    ! !LOCAL VARIABLES:
    deallocate(this%numPatchesPerCol )
    deallocate(this%cohortsPerPatch )
    deallocate(this%balive )
    deallocate(this%bdead )
    deallocate(this%bl )
    deallocate(this%br )
    deallocate(this%bstore )
    deallocate(this%canopy_layer )
    deallocate(this%canopy_trim )
    deallocate(this%dbh )
    deallocate(this%hite )
    deallocate(this%laimemory )
    deallocate(this%leaf_md )
    deallocate(this%root_md )
    deallocate(this%n )
    deallocate(this%gpp_acc )
    deallocate(this%npp_acc )
    deallocate(this%gpp )
    deallocate(this%npp )
    deallocate(this%npp_leaf )
    deallocate(this%npp_froot )
    deallocate(this%npp_bsw )
    deallocate(this%npp_bdead )
    deallocate(this%npp_bseed )
    deallocate(this%npp_store )
    deallocate(this%bmort )
    deallocate(this%hmort )
    deallocate(this%cmort )
    deallocate(this%imort )
    deallocate(this%fmort )
    deallocate(this%ddbhdt )
    deallocate(this%resp_tstep )
    deallocate(this%pft )
    deallocate(this%status_coh )
    deallocate(this%isnew )
    deallocate(this%cwd_ag )
    deallocate(this%cwd_bg )
    deallocate(this%leaf_litter )
    deallocate(this%root_litter )
    deallocate(this%leaf_litter_in )
    deallocate(this%root_litter_in )
    deallocate(this%seed_bank )
    deallocate(this%spread )
    deallocate(this%livegrass )
    deallocate(this%age )
    deallocate(this%areaRestart )
    deallocate(this%f_sun )
    deallocate(this%fabd_sun_z )
    deallocate(this%fabi_sun_z )
    deallocate(this%fabd_sha_z )
    deallocate(this%fabi_sha_z )
    deallocate(this%water_memory )
    deallocate(this%old_stock )
    deallocate(this%cd_status )
    deallocate(this%dd_status )
    deallocate(this%ED_GDD_site )
    deallocate(this%ncd )
    deallocate(this%leafondate )
    deallocate(this%leafoffdate )
    deallocate(this%dleafondate )
    deallocate(this%dleafoffdate )    
    deallocate(this%acc_NI )

  end subroutine deleteEDRestartVectorClass

  !-------------------------------------------------------------------------------!
  function newEDRestartVectorClass( bounds )
    !
    ! !DESCRIPTION:
    ! provide user-defined ctor, with array length argument
    ! allocate memory for vector to write
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)  :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:
    type(EDRestartVectorClass) :: newEDRestartVectorClass
    integer                    :: retVal = 99
    integer, parameter         :: allocOK = 0
    !-----------------------------------------------------------------------

    associate( new => newEDRestartVectorClass)

      ! set class variables
      new%vectorLengthStart = bounds%begCohort
      new%vectorLengthStop  = bounds%endCohort

      ! Column level variables
 
      allocate(new%numPatchesPerCol &
           (bounds%begc:bounds%endc), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%numPatchesPerCol(:) = invalidValue
      
      allocate(new%old_stock &
           (bounds%begc:bounds%endc), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%old_stock(:) = 0.0_r8

      allocate(new%cd_status &
           (bounds%begc:bounds%endc), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%cd_status(:) = 0_r8
      
       allocate(new%dd_status &
           (bounds%begc:bounds%endc), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%dd_status(:) = 0_r8
 
      allocate(new%ncd &
           (bounds%begc:bounds%endc), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%ncd(:) = 0_r8
 

      allocate(new%leafondate &
           (bounds%begc:bounds%endc), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%leafondate(:) = 0_r8
      
      allocate(new%leafoffdate &
           (bounds%begc:bounds%endc), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%leafoffdate(:) = 0_r8     
      
      allocate(new%dleafondate &
           (bounds%begc:bounds%endc), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%dleafondate(:) = 0_r8
      
      allocate(new%dleafoffdate &
           (bounds%begc:bounds%endc), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%dleafoffdate(:) = 0_r8        

      allocate(new%acc_NI &
           (bounds%begc:bounds%endc), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%acc_NI(:) = 0_r8        

      allocate(new%ED_GDD_site &
           (bounds%begc:bounds%endc), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%ED_GDD_site(:) = 0_r8  


      ! cohort level variables



      allocate(new%cohortsPerPatch &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%cohortsPerPatch(:) = invalidValue

      allocate(new%balive &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%balive(:) = 0.0_r8

      allocate(new%bdead &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%bdead(:) = 0.0_r8

      allocate(new%bl &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%bl(:) = 0.0_r8

      allocate(new%br &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%br(:) = 0.0_r8

      allocate(new%bstore &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%bstore(:) = 0.0_r8

      allocate(new%canopy_layer &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%canopy_layer(:) = 0.0_r8

      allocate(new%canopy_trim &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%canopy_trim(:) = 0.0_r8

      allocate(new%dbh &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%dbh(:) = 0.0_r8

      allocate(new%hite &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%hite(:) = 0.0_r8

      allocate(new%laimemory &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%laimemory(:) = 0.0_r8

      allocate(new%leaf_md &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%leaf_md(:) = 0.0_r8

      allocate(new%root_md &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%root_md(:) = 0.0_r8

      allocate(new%n &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%n(:) = 0.0_r8

      allocate(new%gpp_acc &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%gpp_acc(:) = 0.0_r8

      allocate(new%npp_acc &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%npp_acc(:) = 0.0_r8

      allocate(new%gpp &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%gpp(:) = 0.0_r8

      allocate(new%npp &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%npp(:) = 0.0_r8

      allocate(new%npp_leaf &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%npp_leaf(:) = 0.0_r8

      allocate(new%npp_froot &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%npp_froot(:) = 0.0_r8

      allocate(new%npp_bsw &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%npp_bsw(:) = 0.0_r8

      allocate(new%npp_bdead &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%npp_bdead(:) = 0.0_r8

      allocate(new%npp_bseed &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%npp_bseed(:) = 0.0_r8

      allocate(new%npp_store &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%npp_store(:) = 0.0_r8

      allocate(new%bmort &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%bmort(:) = 0.0_r8

      allocate(new%hmort &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%hmort(:) = 0.0_r8

      allocate(new%cmort &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%cmort(:) = 0.0_r8

      allocate(new%imort &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%imort(:) = 0.0_r8

      allocate(new%fmort &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%fmort(:) = 0.0_r8

      allocate(new%ddbhdt &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%ddbhdt(:) = 0.0_r8

      allocate(new%resp_tstep &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%resp_tstep(:) = 0.0_r8

      allocate(new%pft &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%pft(:) = 0

      allocate(new%status_coh &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%status_coh(:) = 0

      allocate(new%isnew &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%isnew(:) = new_cohort

      ! 
      ! some patch level variables that are required on restart
      !
      allocate(new%cwd_ag &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%cwd_ag(:) = 0.0_r8

      allocate(new%cwd_bg &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%cwd_bg(:) = 0.0_r8

      allocate(new%leaf_litter &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%leaf_litter(:) = 0.0_r8

      allocate(new%root_litter &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%root_litter(:) = 0.0_r8

      allocate(new%leaf_litter_in &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%leaf_litter_in(:) = 0.0_r8

      allocate(new%root_litter_in &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%root_litter_in(:) = 0.0_r8

      allocate(new%seed_bank &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%seed_bank(:) = 0.0_r8

      allocate(new%spread &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%spread(:) = 0.0_r8

      allocate(new%livegrass &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%livegrass(:) = 0.0_r8

      allocate(new%age &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%age(:) = 0.0_r8

      allocate(new%areaRestart &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%areaRestart(:) = 0.0_r8

      allocate(new%f_sun &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%f_sun(:) = 0.0_r8

      allocate(new%fabd_sun_z &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%fabd_sun_z(:) = 0.0_r8

      allocate(new%fabi_sun_z &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%fabi_sun_z(:) = 0.0_r8

      allocate(new%fabd_sha_z &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%fabd_sha_z(:) = 0.0_r8

      allocate(new%fabi_sha_z &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%fabi_sha_z(:) = 0.0_r8

      !
      ! Site level variable stored with cohort indexing 
      ! (to accomodate the second dimension)
      !

      allocate(new%water_memory &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      SHR_ASSERT(( retVal == allocOK ), errMsg(sourcefile, __LINE__))
      new%water_memory(:) = 0.0_r8
    

    end associate

  end function newEDRestartVectorClass

  !-------------------------------------------------------------------------------!
  subroutine setVectors( this, bounds, nsites, sites, fcolumn )
    !
    ! !DESCRIPTION:
    ! implement setVectors
    !
    ! !USES:
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(EDRestartVectorClass) , intent(inout)      :: this
    type(bounds_type)           , intent(in)         :: bounds 
    integer                     , intent(in)         :: nsites
    type(ed_site_type)          , intent(in), target :: sites(nsites)
    integer                     , intent(in)         :: fcolumn(nsites)
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    if ( masterproc ) write(iulog,*) 'edtime setVectors ',get_nstep()

    !if (this%DEBUG) then
    !  call this%printIoInfoLL ( bounds, sites, nsites )
    !  call this%printDataInfoLL ( bounds, sites, nsites )
    !end if

    call this%convertCohortListToVector ( bounds, nsites, sites, fcolumn )

    if (this%DEBUG) then
       call this%printIoInfoLL ( bounds, nsites, sites, fcolumn )
       call this%printDataInfoLL ( bounds, nsites, sites )

       ! RGK: Commenting this out because it is calling several
       ! variables over the wrong indices
!       call this%printDataInfoVector (  )
    end if

  end subroutine setVectors

  !-------------------------------------------------------------------------------!
  subroutine getVectors( this, bounds, nsites, sites, fcolumn)
    !
    ! !DESCRIPTION:
    ! implement getVectors
    !
    ! !USES:
    use clm_time_manager , only : get_nstep
    use EDCLMLinkMod     , only : ed_clm_type
    use EDMainMod        , only : ed_update_site
    !
    ! !ARGUMENTS:
    class(EDRestartVectorClass) , intent(inout)         :: this
    type(bounds_type)           , intent(in)            :: bounds 
    integer                     , intent(in)            :: nsites
    type(ed_site_type)          , intent(inout), target :: sites(nsites)
    integer                     , intent(in)            :: fcolumn(nsites)
    !
    ! !LOCAL VARIABLES:
    integer :: s
    !-----------------------------------------------------------------------

    if (this%DEBUG) then
       write(iulog,*) 'edtime getVectors ',get_nstep()
    end if

    call this%createPatchCohortStructure ( bounds, nsites, sites, fcolumn )

    call this%convertCohortVectorToList ( bounds, nsites , sites, fcolumn)

    do s = 1,nsites
       call ed_update_site( sites(s) )
    end do

    if (this%DEBUG) then
       call this%printIoInfoLL ( bounds, nsites, sites, fcolumn )
       call this%printDataInfoLL ( bounds, nsites, sites )
       call this%printDataInfoVector (  )
    end if

  end subroutine getVectors

  !-------------------------------------------------------------------------------!
  subroutine doVectorIO( this, ncid, flag  )
    !
    ! !DESCRIPTION:
    ! implement VectorIO
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_int, ncd_double
    use restUtilMod, only : restartvar
    use clm_varcon,  only : namec, nameCohort
    use spmdMod,     only : iam
    !
    ! !ARGUMENTS:
    class(EDRestartVectorClass), intent(inout) :: this
    type(file_desc_t), intent(inout)           :: ncid   ! netcdf id
    character(len=*) , intent(in)              :: flag   !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical             :: readvar
    character(len=16)   :: coh_dimName  = trim(nameCohort)
    character(len=16)   :: col_dimName  = trim(namec)
    !-----------------------------------------------------------------------


    if(this%DEBUG) then
       write(iulog,*) 'flag:',flag
       write(iulog,*) 'dimname:',col_dimName
       write(iulog,*) 'readvar:',readvar
       write(iulog,*) 'associated?',associated(this%numPatchesPerCol)
       write(iulog,*) ''
       write(iulog,*) 'col size:',size(this%numPatchesPerCol)
       write(iulog,*) 'col lbound:',lbound(this%numPatchesPerCol)
       write(iulog,*) 'col ubound:',ubound(this%numPatchesPerCol)
       
       write(iulog,*) 'coh size:',size(this%cohortsPerPatch)
       write(iulog,*) 'coh lbound:',lbound(this%cohortsPerPatch)
       write(iulog,*) 'coh ubound:',ubound(this%cohortsPerPatch)
       write(iulog,*) ''
    end if
    
    call restartvar(ncid=ncid, flag=flag, varname='ed_io_numPatchesPerCol', xtype=ncd_int,  &
         dim1name=col_dimName, &
         long_name='Num patches per column', units='unitless', &
         interpinic_flag='interp', data=this%numPatchesPerCol, &
         readvar=readvar)
    
    call restartvar(ncid=ncid, flag=flag, varname='ed_old_stock', xtype=ncd_double,  &
         dim1name=col_dimName, &
         long_name='ed cohort - old_stock', units='unitless', &
         interpinic_flag='interp', data=this%old_stock, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_cd_status', xtype=ncd_double,  &
         dim1name=col_dimName, &
         long_name='ed cold dec status', units='unitless', &
         interpinic_flag='interp', data=this%cd_status, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_dd_status', xtype=ncd_double,  &
         dim1name=col_dimName, &
         long_name='ed drought dec status', units='unitless', &
         interpinic_flag='interp', data=this%dd_status, &
         readvar=readvar)         

    call restartvar(ncid=ncid, flag=flag, varname='ed_chilling_days', xtype=ncd_double,  &
         dim1name=col_dimName, &
         long_name='ed chilling day counter', units='unitless', &
         interpinic_flag='interp', data=this%ncd, &
         readvar=readvar)       

    call restartvar(ncid=ncid, flag=flag, varname='ed_leafondate', xtype=ncd_double,  &
         dim1name=col_dimName, &
         long_name='ed leafondate', units='unitless', &
         interpinic_flag='interp', data=this%leafondate, &
         readvar=readvar)         
     
    call restartvar(ncid=ncid, flag=flag, varname='ed_leafoffdate', xtype=ncd_double,  &
         dim1name=col_dimName, &
         long_name='ed leafoffdate', units='unitless', &
         interpinic_flag='interp', data=this%leafoffdate, &
         readvar=readvar) 

    call restartvar(ncid=ncid, flag=flag, varname='ed_dleafondate', xtype=ncd_double,  &
         dim1name=col_dimName, &
         long_name='ed dleafondate', units='unitless', &
         interpinic_flag='interp', data=this%dleafondate, &
         readvar=readvar)         
                   
    call restartvar(ncid=ncid, flag=flag, varname='ed_dleafoffdate', xtype=ncd_double,  &
         dim1name=col_dimName, &
         long_name='ed dleafoffdate', units='unitless', &
         interpinic_flag='interp', data=this%dleafoffdate, &
         readvar=readvar) 
     
    call restartvar(ncid=ncid, flag=flag, varname='ed_acc_NI', xtype=ncd_double,  &
         dim1name=col_dimName, &
         long_name='ed nesterov index', units='unitless', &
         interpinic_flag='interp', data=this%acc_NI, &
         readvar=readvar) 

    call restartvar(ncid=ncid, flag=flag, varname='ed_gdd_site', xtype=ncd_double,  &
         dim1name=col_dimName, &
         long_name='ed GDD site', units='unitless', &
         interpinic_flag='interp', data=this%ED_GDD_site, &
         readvar=readvar) 

    call restartvar(ncid=ncid, flag=flag, varname='ed_balive', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort ed_balive', units='unitless', &
         interpinic_flag='interp', data=this%balive, &
         readvar=readvar)

    !
    ! cohort level vars
    !
    call restartvar(ncid=ncid, flag=flag, varname='ed_io_cohortsPerPatch', xtype=ncd_int,  &
         dim1name=coh_dimName, &
         long_name='cohorts per patch, indexed by numPatchesPerCol', units='unitless', &
         interpinic_flag='interp', data=this%cohortsPerPatch, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_bdead', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - bdead', units='unitless', &
         interpinic_flag='interp', data=this%bdead, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_bl', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - bl', units='unitless', &
         interpinic_flag='interp', data=this%bl, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_br', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - br', units='unitless', &
         interpinic_flag='interp', data=this%br, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_bstore', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - bstore', units='unitless', &
         interpinic_flag='interp', data=this%bstore, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_canopy_layer', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - canopy_layer', units='unitless', &
         interpinic_flag='interp', data=this%canopy_layer, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_canopy_trim', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - canopy_trim', units='unitless', &
         interpinic_flag='interp', data=this%canopy_trim, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_dbh', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - dbh', units='unitless', &
         interpinic_flag='interp', data=this%dbh, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_hite', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - hite', units='unitless', &
         interpinic_flag='interp', data=this%hite, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_laimemory', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - laimemory', units='unitless', &
         interpinic_flag='interp', data=this%laimemory, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_leaf_md', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - leaf_md', units='unitless', &
         interpinic_flag='interp', data=this%leaf_md, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_root_md', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - root_md', units='unitless', &
         interpinic_flag='interp', data=this%root_md, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_n', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - n', units='unitless', &
         interpinic_flag='interp', data=this%n, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_gpp_acc', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - gpp_acc', units='unitless', &
         interpinic_flag='interp', data=this%gpp_acc, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_npp_acc', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - npp_acc', units='unitless', &
         interpinic_flag='interp', data=this%npp_acc, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_gpp', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - gpp', units='unitless', &
         interpinic_flag='interp', data=this%gpp, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_npp', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - npp', units='unitless', &
         interpinic_flag='interp', data=this%npp, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_npp_leaf', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - npp_leaf', units='unitless', &
         interpinic_flag='interp', data=this%npp_leaf, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_npp_froot', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - npp_froot', units='unitless', &
         interpinic_flag='interp', data=this%npp_froot, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_npp_bsw', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - npp_bsw', units='unitless', &
         interpinic_flag='interp', data=this%npp_bsw, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_npp_bdead', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - npp_bdead', units='unitless', &
         interpinic_flag='interp', data=this%npp_bdead, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_npp_bseed', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - npp_bseed', units='unitless', &
         interpinic_flag='interp', data=this%npp_bseed, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_npp_store', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - npp_store', units='unitless', &
         interpinic_flag='interp', data=this%npp_store, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_bmort', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - bmort', units='unitless', &
         interpinic_flag='interp', data=this%bmort, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_hmort', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - hmort', units='unitless', &
         interpinic_flag='interp', data=this%hmort, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_cmort', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - cmort', units='unitless', &
         interpinic_flag='interp', data=this%cmort, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_imort', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - imort', units='unitless', &
         interpinic_flag='interp', data=this%imort, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_fmort', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - fmort', units='unitless', &
         interpinic_flag='interp', data=this%fmort, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_ddbhdt', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - ddbhdt', units='unitless', &
         interpinic_flag='interp', data=this%ddbhdt, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_resp_tstep', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - resp_tstep', units='unitless', &
         interpinic_flag='interp', data=this%resp_tstep, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_pft', xtype=ncd_int,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - pft', units='unitless', &
         interpinic_flag='interp', data=this%pft, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_status_coh', xtype=ncd_int,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - status_coh', units='unitless', &
         interpinic_flag='interp', data=this%status_coh, &
         readvar=readvar)
    
    call restartvar(ncid=ncid, flag=flag, varname='ed_isnew', xtype=ncd_int,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - isnew', units='unitless', &
         interpinic_flag='interp', data=this%isnew, &
         readvar=readvar)

    !
    ! patch level vars
    !

    call restartvar(ncid=ncid, flag=flag, varname='ed_cwd_ag', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - cwd_ag', units='unitless', &
         interpinic_flag='interp', data=this%cwd_ag, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_cwd_bg', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - cwd_bg', units='unitless', &
         interpinic_flag='interp', data=this%cwd_bg, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_leaf_litter', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - leaf_litter', units='unitless', &
         interpinic_flag='interp', data=this%leaf_litter, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_root_litter', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - root_litter', units='unitless', &
         interpinic_flag='interp', data=this%root_litter, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_leaf_litter_in', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - leaf_litter_in', units='unitless', &
         interpinic_flag='interp', data=this%leaf_litter_in, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_root_litter_in', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - root_litter_in', units='unitless', &
         interpinic_flag='interp', data=this%root_litter_in, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_seed_bank', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - seed_bank', units='unitless', &
         interpinic_flag='interp', data=this%seed_bank, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_spread', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - spread', units='unitless', &
         interpinic_flag='interp', data=this%spread, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_livegrass', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - livegrass', units='unitless', &
         interpinic_flag='interp', data=this%livegrass, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_age', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - age', units='unitless', &
         interpinic_flag='interp', data=this%age, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_area', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - area', units='unitless', &
         interpinic_flag='interp', data=this%areaRestart, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_f_sun', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - f_sun', units='unitless', &
         interpinic_flag='interp', data=this%f_sun, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_fabd_sun_z', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - fabd_sun_z', units='unitless', &
         interpinic_flag='interp', data=this%fabd_sun_z, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_fabi_sun_z', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - fabi_sun_z', units='unitless', &
         interpinic_flag='interp', data=this%fabi_sun_z, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_fabd_sha_z', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - fabd_sha_z', units='unitless', &
         interpinic_flag='interp', data=this%fabd_sha_z, &
         readvar=readvar)

    call restartvar(ncid=ncid, flag=flag, varname='ed_fabi_sha_z', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed patch - fabi_sha_z', units='unitless', &
         interpinic_flag='interp', data=this%fabi_sha_z, &
         readvar=readvar)
    !
    ! site level vars
    !

    call restartvar(ncid=ncid, flag=flag, varname='ed_water_memory', xtype=ncd_double,  &
         dim1name=coh_dimName, &
         long_name='ed cohort - water_memory', units='unitless', &
         interpinic_flag='interp', data=this%water_memory, &
         readvar=readvar)

  end subroutine doVectorIO

  !-------------------------------------------------------------------------------!
  subroutine printDataInfoVector( this )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(EDRestartVectorClass), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    character(len=32)   :: methodName = 'PDIV '
    integer :: iSta, iSto
    !-----------------------------------------------------------------------

    ! RGK: changed the vector end-point on column variables to match the start point
    !      this avoids exceeding bounds on the last column of the dataset

    iSta = this%vectorLengthStart
    iSto = iSta + 1

    write(iulog,*) trim(methodName)//' :: this%vectorLengthStart ', &
         this%vectorLengthStart
    write(iulog,*) trim(methodName)//' :: this%vectorLengthStop  ', &
         this%vectorLengthStop

    write(iulog,*) ' PDIV chk ',iSta,iSto
    write(iulog,*) trim(methodName)//' :: balive ', &
         this%balive(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: bdead ', &
         this%bdead(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: bl ', &
         this%bl(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: br ', &
         this%br(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: bstore ', &
         this%bstore(iSta:iSto)

    write(iulog,*) trim(methodName)//' :: canopy_layer ', &
         this%canopy_layer(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: canopy_trim ', &
         this%canopy_trim(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: dbh ', &
         this%dbh(iSta:iSto)

    write(iulog,*) trim(methodName)//' :: hite ', &
         this%hite(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: laimemory ', &
         this%laimemory(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: leaf_md ', &
         this%leaf_md(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: root_md ', &
         this%root_md(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: n ', &
         this%n(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: gpp_acc ', &
         this%gpp_acc(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: npp_acc ', &
         this%npp_acc(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: gpp ', &
         this%gpp(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: npp ', &
         this%npp(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: npp_leaf ', &
         this%npp_leaf(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: npp_froot ', &
         this%npp_froot(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: npp_bsw ', &
         this%npp_bsw(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: npp_bdead ', &
         this%npp_bdead(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: npp_bseed ', &
         this%npp_bseed(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: npp_store ', &
         this%npp_store(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: bmort ', &
         this%bmort(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: hmort ', &
         this%hmort(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: cmort ', &
         this%cmort(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: imort ', &
         this%imort(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: fmort ', &
         this%fmort(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: ddbhdt ', &
         this%ddbhdt(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: resp_tstep ', &
         this%resp_tstep(iSta:iSto)

    write(iulog,*) trim(methodName)//' :: pft ', &
         this%pft(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: status_coh ', &
         this%status_coh(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: isnew ', &
         this%isnew(iSta:iSto)

    write(iulog,*) trim(methodName)//' :: cwd_ag ', &
         this%cwd_ag(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: cwd_bg ', &
         this%cwd_bg(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: leaf_litter ', &
         this%leaf_litter(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: root_litter ', &
         this%root_litter(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: leaf_litter_in ', &
         this%leaf_litter_in(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: root_litter_in ', &
         this%root_litter_in(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: seed_bank ', &
         this%seed_bank(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: spread ', &
         this%spread(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: livegrass ', &
         this%livegrass(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: age ', &
         this%age(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: area ', &
         this%areaRestart(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: f_sun ', &
         this%f_sun(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: fabd_sun_z ', &
         this%fabd_sun_z(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: fabi_sun_z ', &
         this%fabi_sun_z(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: fabd_sha_z ', &
         this%fabd_sha_z(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: fabi_sha_z ', &
         this%fabi_sha_z(iSta:iSto)
    write(iulog,*) trim(methodName)//' :: water_memory ', &
         this%water_memory(iSta:iSto)
 
    write(iulog,*) trim(methodName)//' :: old_stock ', &
         this%old_stock(iSta:iSta)
    write(iulog,*) trim(methodName)//' :: cd_status', &
         this%cd_status(iSta:iSta)
    write(iulog,*) trim(methodName)//' :: dd_status', &
         this%cd_status(iSta:iSto)  
    write(iulog,*) trim(methodName)//' :: ED_GDD_site', &
         this%ED_GDD_site(iSta:iSto)  
    write(iulog,*) trim(methodName)//' :: ncd', &
         this%ncd(iSta:iSta)   
    write(iulog,*) trim(methodName)//' :: leafondate', &
         this%leafondate(iSta:iSta) 
     write(iulog,*) trim(methodName)//' :: leafoffdate', &
         this%leafoffdate(iSta:iSta) 
    write(iulog,*) trim(methodName)//' :: dleafondate', &
         this%dleafondate(iSta:iSta)  
      write(iulog,*) trim(methodName)//' :: dleafoffdate', &
         this%dleafoffdate(iSta:iSta) 
    write(iulog,*) trim(methodName)//' :: acc_NI', &
         this%acc_NI(iSta:iSta)                          
         
  end subroutine printDataInfoVector

  !-------------------------------------------------------------------------------!
  subroutine printDataInfoLL( this, bounds, nsites, sites ) 
    !
    ! !DESCRIPTION:
    ! counts the total number of cohorts over all p levels (ed_patch_type) so we
    ! can allocate vectors, copy from LL -> vector and read/write restarts.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(EDRestartVectorClass) , intent(inout)      :: this
    type(bounds_type)           , intent(in)         :: bounds 
    integer                     , intent(in)         :: nsites
    type(ed_site_type)          , intent(in), target :: sites(nsites)
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type),  pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort
    integer :: s
    integer :: totalCohorts
    integer :: numCohort
    integer :: numPatches,totPatchCount
    character(len=32) :: methodName = 'printDataInfoLL '
    !-----------------------------------------------------------------------

    totalCohorts = 0
    totPatchCount = 1

    write(iulog,*) 'vecLenStart ',this%vectorLengthStart

    do s = 1,nsites

       currentPatch => sites(s)%oldest_patch

          numPatches = 1

          do while(associated(currentPatch))
             currentCohort => currentPatch%shortest

             numCohort = 0

             do while(associated(currentCohort))  

                totalCohorts = totalCohorts + 1

                write(iulog,*) trim(methodName)//' balive '       ,totalCohorts,currentCohort%balive
                write(iulog,*) trim(methodName)//' bdead '        ,totalCohorts,currentCohort%bdead
                write(iulog,*) trim(methodName)//' bl '           ,totalCohorts,currentCohort%bl
                write(iulog,*) trim(methodName)//' br '           ,totalCohorts,currentCohort%br
                write(iulog,*) trim(methodName)//' bstore '       ,totalCohorts,currentCohort%bstore
                write(iulog,*) trim(methodName)//' canopy_layer ' ,totalCohorts,currentCohort%canopy_layer
                write(iulog,*) trim(methodName)//' canopy_trim '  ,totalCohorts,currentCohort%canopy_trim
                write(iulog,*) trim(methodName)//' dbh '          ,totalCohorts,currentCohort%dbh
                write(iulog,*) trim(methodName)//' hite '         ,totalCohorts,currentCohort%hite
                write(iulog,*) trim(methodName)//' laimemory '    ,totalCohorts,currentCohort%laimemory
                write(iulog,*) trim(methodName)//' leaf_md '      ,totalCohorts,currentCohort%leaf_md
                write(iulog,*) trim(methodName)//' root_md '      ,totalCohorts,currentCohort%root_md
                write(iulog,*) trim(methodName)//' n '            ,totalCohorts,currentCohort%n
                write(iulog,*) trim(methodName)//' gpp_acc '      ,totalCohorts,currentCohort%gpp_acc
                write(iulog,*) trim(methodName)//' npp_acc '      ,totalCohorts,currentCohort%npp_acc
                write(iulog,*) trim(methodName)//' gpp '          ,totalCohorts,currentCohort%gpp
                write(iulog,*) trim(methodName)//' npp '          ,totalCohorts,currentCohort%npp
                write(iulog,*) trim(methodName)//' npp_leaf '     ,totalCohorts,currentCohort%npp_leaf
                write(iulog,*) trim(methodName)//' npp_froot '    ,totalCohorts,currentCohort%npp_froot
                write(iulog,*) trim(methodName)//' npp_bsw '      ,totalCohorts,currentCohort%npp_bsw
                write(iulog,*) trim(methodName)//' npp_bdead '    ,totalCohorts,currentCohort%npp_bdead
                write(iulog,*) trim(methodName)//' npp_bseed '    ,totalCohorts,currentCohort%npp_bseed
                write(iulog,*) trim(methodName)//' npp_store '    ,totalCohorts,currentCohort%npp_store
                write(iulog,*) trim(methodName)//' bmort '        ,totalCohorts,currentCohort%bmort
                write(iulog,*) trim(methodName)//' hmort '        ,totalCohorts,currentCohort%hmort
                write(iulog,*) trim(methodName)//' cmort '        ,totalCohorts,currentCohort%cmort
                write(iulog,*) trim(methodName)//' imort '        ,totalCohorts,currentCohort%imort
                write(iulog,*) trim(methodName)//' fmort '        ,totalCohorts,currentCohort%fmort
                write(iulog,*) trim(methodName)//' ddbhdt '       ,totalCohorts,currentCohort%ddbhdt
                write(iulog,*) trim(methodName)//' resp_tstep '     ,totalCohorts,currentCohort%resp_tstep
                write(iulog,*) trim(methodName)//' pft '          ,totalCohorts,currentCohort%pft
                write(iulog,*) trim(methodName)//' status_coh '   ,totalCohorts,currentCohort%status_coh
                write(iulog,*) trim(methodName)//' isnew '        ,totalCohorts,currentCohort%isnew

                numCohort = numCohort + 1

                currentCohort => currentCohort%taller
             enddo ! currentCohort do while

             write(iulog,*) trim(methodName)//': numpatches for col ',&
                  numPatches

             write(iulog,*) trim(methodName)//': patches and cohorts ',&
                  totPatchCount,numCohort

             write(iulog,*) trim(methodName)//' cwd_ag '         ,currentPatch%cwd_ag
             write(iulog,*) trim(methodName)//' cwd_bg '         ,currentPatch%cwd_bg
             write(iulog,*) trim(methodName)//' leaf_litter '    ,currentPatch%leaf_litter
             write(iulog,*) trim(methodName)//' root_litter '    ,currentPatch%root_litter
             write(iulog,*) trim(methodName)//' leaf_litter_in ' ,currentPatch%leaf_litter_in
             write(iulog,*) trim(methodName)//' root_litter_in ' ,currentPatch%root_litter_in
             write(iulog,*) trim(methodName)//' seed_bank '      ,currentPatch%seed_bank
             write(iulog,*) trim(methodName)//' spread '         ,currentPatch%spread
             write(iulog,*) trim(methodName)//' livegrass '      ,currentPatch%livegrass
             write(iulog,*) trim(methodName)//' age '            ,currentPatch%age
             write(iulog,*) trim(methodName)//' area '           ,currentPatch%area
             write(iulog,*) trim(methodName)//' f_sun (sum) '     ,sum(currentPatch%f_sun)
             write(iulog,*) trim(methodName)//' fabd_sun_z (sum) '     ,sum(currentPatch%fabd_sun_z)
             write(iulog,*) trim(methodName)//' fabi_sun_z (sum) '     ,sum(currentPatch%fabi_sun_z)
             write(iulog,*) trim(methodName)//' fabd_sha_z (sum) '     ,sum(currentPatch%fabd_sha_z)
             write(iulog,*) trim(methodName)//' fabi_sha_z (sum) '     ,sum(currentPatch%fabi_sha_z)

             write(iulog,*) trim(methodName)//' old_stock '      ,sites(s)%old_stock
             write(iulog,*) trim(methodName)//' cd_status '      ,sites(s)%status
             write(iulog,*) trim(methodName)//' dd_status '      ,sites(s)%dstatus
             write(iulog,*) trim(methodName)//' ncd '            ,sites(s)%ncd
             write(iulog,*) trim(methodName)//' leafondate '     ,sites(s)%leafondate
             write(iulog,*) trim(methodName)//' leafoffdate '    ,sites(s)%leafoffdate
             write(iulog,*) trim(methodName)//' dleafondate '    ,sites(s)%dleafondate
             write(iulog,*) trim(methodName)//' dleafoffdate '   ,sites(s)%dleafoffdate
             write(iulog,*) trim(methodName)//' acc_NI'          ,sites(s)%acc_NI
             write(iulog,*) trim(methodName)//' ED_GDD_site '    ,sites(s)%ED_GDD_site

             currentPatch => currentPatch%younger

             totPatchCount = totPatchCount + 1
             numPatches = numPatches + 1
          enddo ! currentPatch do while

       write(iulog,*) trim(methodName)//' water_memory ',sites(s)%water_memory(1)

    enddo

    write(iulog,*) trim(methodName)//': total cohorts ',totalCohorts

  end subroutine printDataInfoLL

  !-------------------------------------------------------------------------------!
  subroutine printIoInfoLL( this, bounds, nsites, sites, fcolumn ) 
    !!
    ! !DESCRIPTION:
    ! for debugging.  prints some IO info regarding cohorts/patches
    ! currently prints cohort level variables
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(EDRestartVectorClass) , intent(inout)      :: this
    type(bounds_type)           , intent(in)         :: bounds 
    integer                     , intent(in)         :: nsites
    type(ed_site_type)          , intent(in), target :: sites(nsites)
    integer                     , intent(in)         :: fcolumn(nsites)
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type),  pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort
    integer s
    integer totalCohorts
    integer numCohort
    integer numPatches,totPatchCount
    character(len=32)   :: methodName = 'printIoInfoLL '
    !-----------------------------------------------------------------------

    totalCohorts = 0
    totPatchCount = 1

    write(iulog,*) 'vecLenStart ',this%vectorLengthStart

    do s = 1,nsites
       
       currentPatch => sites(s)%oldest_patch

       numPatches = 1

       do while(associated(currentPatch))
          currentCohort => currentPatch%shortest
          
          write(iulog,*) trim(methodName)//': found column with patch(s) ',fcolumn(s)
          
          numCohort = 0
          
          do while(associated(currentCohort))  
             
             totalCohorts = totalCohorts + 1
             numCohort = numCohort + 1
             
             write(iulog,*) trim(methodName)//' balive       ',numCohort,currentCohort%balive
             write(iulog,*) trim(methodName)//' bdead        ',currentCohort%bdead
             write(iulog,*) trim(methodName)//' bl           ',currentCohort%bl
             write(iulog,*) trim(methodName)//' br           ',currentCohort%br
             write(iulog,*) trim(methodName)//' bstore       ',currentCohort%bstore
             write(iulog,*) trim(methodName)//' canopy_layer ',currentCohort%canopy_layer
             write(iulog,*) trim(methodName)//' canopy_trim  ',currentCohort%canopy_trim
             write(iulog,*) trim(methodName)//' dbh          ',currentCohort%dbh
             write(iulog,*) trim(methodName)//' hite         ',currentCohort%hite
             write(iulog,*) trim(methodName)//' laimemory    ',currentCohort%laimemory
             write(iulog,*) trim(methodName)//' leaf_md      ',currentCohort%leaf_md
             write(iulog,*) trim(methodName)//' root_md      ',currentCohort%root_md
             write(iulog,*) trim(methodName)//' n            ',currentCohort%n
             write(iulog,*) trim(methodName)//' gpp_acc      ',currentCohort%gpp_acc
             write(iulog,*) trim(methodName)//' npp_acc      ',currentCohort%npp_acc
             write(iulog,*) trim(methodName)//' gpp      ',currentCohort%gpp
             write(iulog,*) trim(methodName)//' npp      ',currentCohort%npp
             write(iulog,*) trim(methodName)//' npp_leaf      ',currentCohort%npp_leaf
             write(iulog,*) trim(methodName)//' npp_froot      ',currentCohort%npp_froot
             write(iulog,*) trim(methodName)//' npp_bsw      ',currentCohort%npp_bsw
             write(iulog,*) trim(methodName)//' npp_bdead      ',currentCohort%npp_bdead
             write(iulog,*) trim(methodName)//' npp_bseed      ',currentCohort%npp_bseed
             write(iulog,*) trim(methodName)//' npp_store      ',currentCohort%npp_store
             write(iulog,*) trim(methodName)//' bmort      ',currentCohort%bmort
             write(iulog,*) trim(methodName)//' hmort      ',currentCohort%hmort
             write(iulog,*) trim(methodName)//' cmort      ',currentCohort%cmort
             write(iulog,*) trim(methodName)//' imort      ',currentCohort%imort
             write(iulog,*) trim(methodName)//' fmort      ',currentCohort%fmort
             write(iulog,*) trim(methodName)//' ddbhdt      ',currentCohort%ddbhdt
             write(iulog,*) trim(methodName)//' resp_tstep     ',currentCohort%resp_tstep
             write(iulog,*) trim(methodName)//' pft          ',currentCohort%pft
             write(iulog,*) trim(methodName)//' status_coh   ',currentCohort%status_coh
             write(iulog,*) trim(methodName)//' isnew        ',currentCohort%isnew
             
             currentCohort => currentCohort%taller
          enddo ! currentCohort do while
          
          write(iulog,*) trim(methodName)//': numpatches for column ',numPatches
          write(iulog,*) trim(methodName)//': patches and cohorts ',totPatchCount,numCohort

          currentPatch => currentPatch%younger

          totPatchCount = totPatchCount + 1
          numPatches = numPatches + 1
       enddo ! currentPatch do while
    enddo
    
    return
  end subroutine printIoInfoLL

  !-------------------------------------------------------------------------------!
  subroutine convertCohortListToVector( this, bounds, nsites, sites, fcolumn ) 
    !
    ! !DESCRIPTION:
    ! counts the total number of cohorts over all p levels (ed_patch_type) so we
    ! can allocate vectors, copy from LL -> vector and read/write restarts.
    !
    ! !USES:
    use EDTypesMod, only : cp_nclmax
    !
    ! !ARGUMENTS:
    class(EDRestartVectorClass) , intent(inout)      :: this
    type(bounds_type)           , intent(in)         :: bounds 
    integer                     , intent(in)         :: nsites
    type(ed_site_type)          , intent(in), target :: sites(nsites)
    integer                     , intent(in)         :: fcolumn(nsites)
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type), pointer  :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort
    integer ::  s, c
    integer ::  totalCohorts ! number of cohorts starting from 1
    integer ::  countCohort  ! number of cohorts starting from
    ! vectorLengthStart
    integer :: numCohort
    integer :: numPatches
    integer :: totPatchCount, offsetTotPatchCount
    integer :: countPft
    integer :: countNcwd
    integer :: countWaterMem
    integer :: countNclmax
    integer :: countSunZ
    integer :: i,j,k
    integer :: incrementOffset
    !-----------------------------------------------------------------------

    totalCohorts = 0

!    if(fcolumn(1).eq.bounds%begc .and. &
!          (fcolumn(1)-1)*cohorts_per_col+1.ne.bounds%begCohort) then
!        write(iulog,*) 'fcolumn(1) in this clump, points to the first column of the clump'
!        write(iulog,*) 'but the assumption on first cohort index does not jive'
!        call endrun(msg=errMsg(sourcefile, __LINE__))
!    end if


    do s = 1,nsites
       
       ! Calculate the offsets
       ! fcolumn is the global column index of the current site.
       ! For the first site, if that site aligns with the first column index
       ! in the clump, than the offset should be be equal to begCohort
       
       c = fcolumn(s)

!       incrementOffset     = (c-1)*cohorts_per_col + 1
!       countCohort         = (c-1)*cohorts_per_col + 1
!       countPft            = (c-1)*cohorts_per_col + 1
!       countNcwd           = (c-1)*cohorts_per_col + 1
!       countNclmax         = (c-1)*cohorts_per_col + 1
!       countWaterMem       = (c-1)*cohorts_per_col + 1
!       countSunZ           = (c-1)*cohorts_per_col + 1

       incrementOffset     = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1
       countCohort         = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1
       countPft            = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1
       countNcwd           = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1
       countNclmax         = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1
       countWaterMem       = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1
       countSunZ           = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1

       currentPatch => sites(s)%oldest_patch

       ! new column, reset num patches
       numPatches = 0

       do while(associated(currentPatch))
          
          ! found patch, increment
          numPatches = numPatches + 1
          
          currentCohort => currentPatch%shortest
          
          ! new patch, reset num cohorts
          numCohort = 0
          
          do while(associated(currentCohort))

             ! found cohort, increment
             numCohort        = numCohort    + 1
             totalCohorts     = totalCohorts + 1
             
             if (this%DEBUG) then
                write(iulog,*) 'CLTV countCohort ', countCohort
                write(iulog,*) 'CLTV vecLenStart ', this%vectorLengthStart
                write(iulog,*) 'CLTV vecLenStop  ', this%vectorLengthStop
             endif
                
             this%balive(countCohort)       = currentCohort%balive
             this%bdead(countCohort)        = currentCohort%bdead
             this%bl(countCohort)           = currentCohort%bl
             this%br(countCohort)           = currentCohort%br
             this%bstore(countCohort)       = currentCohort%bstore
             this%canopy_layer(countCohort) = currentCohort%canopy_layer
             this%canopy_trim(countCohort)  = currentCohort%canopy_trim
             this%dbh(countCohort)          = currentCohort%dbh
             this%hite(countCohort)         = currentCohort%hite
             this%laimemory(countCohort)    = currentCohort%laimemory
             this%leaf_md(countCohort)      = currentCohort%leaf_md
             this%root_md(countCohort)      = currentCohort%root_md
             this%n(countCohort)            = currentCohort%n
             this%gpp_acc(countCohort)      = currentCohort%gpp_acc
             this%npp_acc(countCohort)      = currentCohort%npp_acc
             this%gpp(countCohort)      = currentCohort%gpp
             this%npp(countCohort)      = currentCohort%npp
             this%npp_leaf(countCohort)      = currentCohort%npp_leaf
             this%npp_froot(countCohort)      = currentCohort%npp_froot
             this%npp_bsw(countCohort)      = currentCohort%npp_bsw
             this%npp_bdead(countCohort)      = currentCohort%npp_bdead
             this%npp_bseed(countCohort)      = currentCohort%npp_bseed
             this%npp_store(countCohort)      = currentCohort%npp_store
             this%bmort(countCohort)      = currentCohort%bmort
             this%hmort(countCohort)      = currentCohort%hmort
             this%cmort(countCohort)      = currentCohort%cmort
             this%imort(countCohort)      = currentCohort%imort
             this%fmort(countCohort)      = currentCohort%fmort
             this%ddbhdt(countCohort)      = currentCohort%ddbhdt
             this%resp_tstep(countCohort)     = currentCohort%resp_tstep
             this%pft(countCohort)          = currentCohort%pft
             this%status_coh(countCohort)   = currentCohort%status_coh
             if ( currentCohort%isnew ) then
                this%isnew(countCohort)        = new_cohort
             else
                this%isnew(countCohort)        = old_cohort
             endif
             
             if (this%DEBUG) then
                write(iulog,*) 'CLTV offsetNumCohorts II ',countCohort, &
                      numCohort
             endif
             
             countCohort = countCohort + 1
             
             currentCohort => currentCohort%taller
             
          enddo ! currentCohort do while

          !
          ! deal with patch level fields here
          !
          this%livegrass(incrementOffset)   = currentPatch%livegrass
          this%age(incrementOffset)         = currentPatch%age
          this%areaRestart(incrementOffset) = currentPatch%area
          
          ! set cohorts per patch for IO
          this%cohortsPerPatch( incrementOffset ) = numCohort
          
          if (this%DEBUG) then
             write(iulog,*) 'offsetNumCohorts III ' &
                   ,countCohort,cohorts_per_col, numCohort
          endif
          !
          ! deal with patch level fields of arrays here
          !
          ! these are arrays of length numpft_ed, each patch contains one
          ! vector so we increment 
          do i = 1,numpft_ed 
             this%leaf_litter(countPft)    = currentPatch%leaf_litter(i)
             this%root_litter(countPft)    = currentPatch%root_litter(i)
             this%leaf_litter_in(countPft) = currentPatch%leaf_litter_in(i)
             this%root_litter_in(countPft) = currentPatch%root_litter_in(i)
             this%seed_bank(countPft)      = currentPatch%seed_bank(i)
             countPft = countPft + 1
          end do
          
          do i = 1,ncwd ! ncwd currently 4
             this%cwd_ag(countNcwd) = currentPatch%cwd_ag(i)
             this%cwd_bg(countNcwd) = currentPatch%cwd_bg(i)
             countNcwd = countNcwd + 1
          end do
          
          do i = 1,cp_nclmax ! cp_nclmax currently 2
             this%spread(countNclmax)         = currentPatch%spread(i)
             countNclmax = countNclmax + 1
          end do
          
          if (this%DEBUG) write(iulog,*) 'CLTV countSunZ 1 ',countSunZ
          
          if (this%DEBUG) write(iulog,*) 'CLTV 1186 ',cp_nlevcan,numpft_ed,cp_nclmax
          
          do k = 1,cp_nlevcan ! cp_nlevcan currently 40
             do j = 1,numpft_ed ! numpft_ed currently 2
                do i = 1,cp_nclmax ! cp_nclmax currently 2
                   this%f_sun(countSunZ)       = currentPatch%f_sun(i,j,k)
                   this%fabd_sun_z(countSunZ)  = currentPatch%fabd_sun_z(i,j,k)
                   this%fabi_sun_z(countSunZ)  = currentPatch%fabi_sun_z(i,j,k)
                   this%fabd_sha_z(countSunZ)  = currentPatch%fabd_sha_z(i,j,k)
                   this%fabi_sha_z(countSunZ)  = currentPatch%fabi_sha_z(i,j,k)
                   countSunZ = countSunZ + 1
                end do
             end do
          end do
          
          if (this%DEBUG) write(iulog,*) 'CLTV countSunZ 2 ',countSunZ
          
          incrementOffset = incrementOffset + numCohortsPerPatch
          
          ! reset counters so that they are all advanced evenly. Currently
          ! the offset is 10, the max of numpft_ed, ncwd, cp_nclmax,
          ! countWaterMem and the number of allowed cohorts per patch
          countPft      = incrementOffset
          countNcwd     = incrementOffset
          countNclmax   = incrementOffset
          countCohort   = incrementOffset
          countSunZ     = incrementOffset
          
          if (this%DEBUG) then
             write(iulog,*) 'CLTV incrementOffset ', incrementOffset
             write(iulog,*) 'CLTV cohorts_per_col ', cohorts_per_col
             write(iulog,*) 'CLTV numCohort ', numCohort
             write(iulog,*) 'CLTV totalCohorts ', totalCohorts
          end if
          
          currentPatch => currentPatch%younger
          
       enddo ! currentPatch do while

       this%old_stock(c)    = sites(s)%old_stock
       this%cd_status(c)    = sites(s)%status
       this%dd_status(c)    = sites(s)%dstatus
       this%ncd(c)          = sites(s)%ncd 
       this%leafondate(c)   = sites(s)%leafondate
       this%leafoffdate(c)  = sites(s)%leafoffdate
       this%dleafondate(c)  = sites(s)%dleafondate
       this%dleafoffdate(c) = sites(s)%dleafoffdate
       this%acc_NI(c)       = sites(s)%acc_NI
       this%ED_GDD_site(c)  = sites(s)%ED_GDD_site
       
       ! set numpatches for this column
       this%numPatchesPerCol(c)  = numPatches
       
       do i = 1,numWaterMem ! numWaterMem currently 10
          this%water_memory( countWaterMem ) = sites(s)%water_memory(i)
          countWaterMem = countWaterMem + 1
       end do
       
    enddo
    
    if (this%DEBUG) then
       write(iulog,*) 'CLTV total cohorts ',totalCohorts
    end if
    
    return
 end subroutine convertCohortListToVector

  !-------------------------------------------------------------------------------!
  subroutine createPatchCohortStructure( this, bounds, nsites, sites, fcolumn ) 
    !
    ! !DESCRIPTION:
    ! counts the total number of cohorts over all p levels (ed_patch_type) so we
    ! can allocate vectors, copy from LL -> vector and read/write restarts.
    !
    ! !USES:
    use EDPatchDynamicsMod ,  only : zero_patch
    use EDGrowthFunctionsMod, only : Dbh
    use EDCohortDynamicsMod,  only : create_cohort
    use EDInitMod          ,  only : zero_site
    use EDParamsMod        ,  only : ED_val_maxspread
    use EDPatchDynamicsMod ,  only : create_patch
    use GridcellType       ,  only : grc
    use ColumnType         ,  only : col
    !
    ! !ARGUMENTS:
    class(EDRestartVectorClass) , intent(inout)         :: this
    type(bounds_type)           , intent(in)            :: bounds 
    integer                     , intent(in)            :: nsites
    type(ed_site_type)          , intent(inout), target :: sites(nsites)
    integer                     , intent(in)            :: fcolumn(nsites)
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type) , pointer  :: newp
    type(ed_cohort_type), allocatable :: temp_cohort
    real(r8) :: cwd_ag_local(ncwd),cwd_bg_local(ncwd),spread_local(cp_nclmax)
    real(r8) :: leaf_litter_local(numpft_ed),root_litter_local(numpft_ed)
    real(r8) :: seed_bank_local(numpft_ed)
    real(r8) :: age !notional age of this patch
    integer  :: cohortstatus
    integer  :: s  ! site index
    integer  :: c  ! column index
    integer  :: g  ! grid index
    integer  :: patchIdx,currIdx, fto, ft
    !-----------------------------------------------------------------------

   

    cwd_ag_local      = 0.0_r8 !ED_val_init_litter   !arbitrary value for litter pools. kgC m-2           ! 
    cwd_bg_local      = 0.0_r8 !ED_val_init_litter
    leaf_litter_local = 0.0_r8
    root_litter_local = 0.0_r8
    age               = 0.0_r8
    spread_local      = ED_val_maxspread

    ! 
    ! loop over model grid cells and create patch/cohort structure based on
    ! restart data
    !
    do s = 1,nsites

       c = fcolumn(s)
       if( (s-1) .ne. (c-bounds%begc) ) then
          write(iulog,*) 'NAT COLUMNS REALLY ARENT MONOTONICALLY INCREASING'
          write(iulog,*) s,c,bounds%begc,s-1,c-bounds%begc
       end if

       g = col%gridcell(c)
       
       currIdx = bounds%begCohort + (c-bounds%begc)*cohorts_per_col + 1
!       currIdx = (c-1)*cohorts_per_col + 1  ! global cohort index at the head of the column

       call zero_site( sites(s) )
       ! 
       ! set a few items that are necessary on restart for ED but not on the 
       ! restart file
       !

       sites(s)%lat = grc%latdeg(g)
       sites(s)%lon = grc%londeg(g)
       sites(s)%ncd = 0.0_r8

       if (this%numPatchesPerCol(c)<0 .or. this%numPatchesPerCol(c)>10000) then
          write(iulog,*) 'a column was expected to contain a valid number of patches'
          write(iulog,*) '0 is a valid number, but this column seems uninitialized',this%numPatchesPerCol(c)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       ! Initialize the site pointers to null
       sites(s)%youngest_patch         => null()                 
       sites(s)%oldest_patch           => null()

       do patchIdx = 1,this%numPatchesPerCol(c)

          if (this%DEBUG) then
             write(iulog,*) 'create patch ',patchIdx
             write(iulog,*) 'patchIdx 1-numCohorts : ',this%cohortsPerPatch(currIdx)
          end if

          ! create patch
          allocate(newp)    

          ! make new patch
          call create_patch(sites(s), newp, age, area, &
               spread_local, cwd_ag_local, cwd_bg_local,  &
               leaf_litter_local, root_litter_local, seed_bank_local) 

          newp%siteptr => sites(s)

          ! give this patch a unique patch number
          newp%patchno = patchIdx

          do fto = 1, this%cohortsPerPatch(currIdx)

             allocate(temp_cohort)

             temp_cohort%n = 700.0_r8
             temp_cohort%balive = 0.0_r8
             temp_cohort%bdead = 0.0_r8
             temp_cohort%bstore = 0.0_r8
             temp_cohort%laimemory = 0.0_r8
             temp_cohort%canopy_trim = 0.0_r8
             temp_cohort%canopy_layer = 1.0_r8

             ! set the pft (only 2 used in ed) based on odd/even cohort
             ! number
             ft=2
             if ((mod(fto, 2)  ==  0 )) then
                ft=1
             endif

             cohortstatus = newp%siteptr%status

             if(pftcon%stress_decid(ft) == 1)then !drought decidous, override status. 
                cohortstatus = newp%siteptr%dstatus
             endif

             temp_cohort%hite = 1.25_r8
             ! the dbh function should only take as an argument, the one
             ! item it needs, not the entire cohort...refactor
             temp_cohort%dbh = Dbh(temp_cohort) + 0.0001_r8*ft

             write(iulog,*) 'EDRestVectorMod.F90::createPatchCohortStructure call create_cohort '

             call create_cohort(newp, ft, temp_cohort%n, temp_cohort%hite, temp_cohort%dbh, &
                  temp_cohort%balive, temp_cohort%bdead, temp_cohort%bstore,  &
                  temp_cohort%laimemory, cohortstatus, temp_cohort%canopy_trim, newp%NCL_p)

             deallocate(temp_cohort)

          enddo ! ends loop over fto

          !
          ! insert this patch with cohorts into the site pointer.  At this
          ! point just insert the new patch in the youngest position
          !
          if (patchIdx == 1) then ! nothing associated yet. first patch is pointed to by youngest and oldest

             if (this%DEBUG) write(iulog,*) 'patchIdx = 1 ',patchIdx

             sites(s)%youngest_patch         => newp                   
             sites(s)%oldest_patch           => newp                        
             sites(s)%youngest_patch%younger => null()
             sites(s)%youngest_patch%older   => null()
             sites(s)%oldest_patch%younger   => null()
             sites(s)%oldest_patch%older     => null()

          else if (patchIdx == 2) then ! add second patch to list

             if (this%DEBUG) write(iulog,*) 'patchIdx = 2 ',patchIdx

             sites(s)%youngest_patch         => newp
             sites(s)%youngest_patch%younger => null()
             sites(s)%youngest_patch%older   => sites(s)%oldest_patch
             sites(s)%oldest_patch%younger   => sites(s)%youngest_patch
             sites(s)%oldest_patch%older     => null()

          else ! more than 2 patches, insert patch into youngest slot

             if (this%DEBUG) write(iulog,*) 'patchIdx > 2 ',patchIdx

             newp%older                      => sites(s)%youngest_patch
             sites(s)%youngest_patch%younger => newp
             newp%younger                    => null()
             sites(s)%youngest_patch         => newp

          endif

          currIdx = currIdx + numCohortsPerPatch

       enddo ! ends loop over patchIdx

    enddo ! ends loop over s

 end subroutine createPatchCohortStructure

  !-------------------------------------------------------------------------------!
  subroutine convertCohortVectorToList( this, bounds, nsites, sites, fcolumn ) 
    !
    ! !DESCRIPTION:
    ! counts the total number of cohorts over all p levels (ed_patch_type) so we
    ! can allocate vectors, copy from LL -> vector and read/write restarts.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(EDRestartVectorClass) , intent(inout)         :: this
    type(bounds_type)           , intent(in)            :: bounds 
    integer                     , intent(in)            :: nsites
    type(ed_site_type)          , intent(inout), target :: sites(nsites)
    integer                     , intent(in)            :: fcolumn(nsites)

    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type), pointer :: currentPatch
    type (ed_cohort_type),pointer :: currentCohort
    integer :: c, s
    integer :: totalCohorts ! number of cohorts starting from 0
    integer :: countCohort  ! number of cohorts starting from
    ! vectorLengthStart
    integer :: numCohort
    integer :: numPatches
    integer :: countPft
    integer :: countNcwd
    integer :: countWaterMem
    integer :: countNclmax
    integer :: countSunZ
    integer :: i,j,k
    integer :: incrementOffset
    !-----------------------------------------------------------------------

    totalCohorts = 0
    
    do s = 1,nsites
       
       c = fcolumn(s)

       incrementOffset = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1

       countCohort     = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1
       countPft        = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1
       countNcwd       = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1
       countNclmax     = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1
       countWaterMem   = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1
       countSunZ       = bounds%begCohort+(c-bounds%begc)*cohorts_per_col + 1

       currentPatch => sites(s)%oldest_patch
       
       ! new grid cell, reset num patches
       numPatches = 0

       do while(associated(currentPatch))

          ! found patch, increment
          numPatches = numPatches + 1
          
          currentCohort => currentPatch%shortest
          
          ! new patch, reset num cohorts
          numCohort = 0
          
          do while(associated(currentCohort))        
             
             ! found cohort, increment
             numCohort        = numCohort    + 1
             totalCohorts     = totalCohorts + 1
             
             if (this%DEBUG) then
                write(iulog,*) 'CVTL countCohort ',countCohort, this%vectorLengthStart, this%vectorLengthStop
             endif
             
             currentCohort%balive = this%balive(countCohort)
             currentCohort%bdead = this%bdead(countCohort)
             currentCohort%bl = this%bl(countCohort)
             currentCohort%br = this%br(countCohort)
             currentCohort%bstore = this%bstore(countCohort)
             currentCohort%canopy_layer = this%canopy_layer(countCohort)
             currentCohort%canopy_trim = this%canopy_trim(countCohort)
             currentCohort%dbh = this%dbh(countCohort)
             currentCohort%hite = this%hite(countCohort)
             currentCohort%laimemory = this%laimemory(countCohort)
             currentCohort%leaf_md = this%leaf_md(countCohort)
             currentCohort%root_md = this%root_md(countCohort)
             currentCohort%n = this%n(countCohort)
             currentCohort%gpp_acc = this%gpp_acc(countCohort)
             currentCohort%npp_acc = this%npp_acc(countCohort)
             currentCohort%gpp = this%gpp(countCohort)
             currentCohort%npp = this%npp(countCohort)
             currentCohort%npp_leaf = this%npp_leaf(countCohort)
             currentCohort%npp_froot = this%npp_froot(countCohort)
             currentCohort%npp_bsw = this%npp_bsw(countCohort)
             currentCohort%npp_bdead = this%npp_bdead(countCohort)
             currentCohort%npp_bseed = this%npp_bseed(countCohort)
             currentCohort%npp_store = this%npp_store(countCohort)
             currentCohort%bmort = this%bmort(countCohort)
             currentCohort%hmort = this%hmort(countCohort)
             currentCohort%cmort = this%cmort(countCohort)
             currentCohort%imort = this%imort(countCohort)
             currentCohort%fmort = this%fmort(countCohort)
             currentCohort%ddbhdt = this%ddbhdt(countCohort)
             currentCohort%resp_tstep = this%resp_tstep(countCohort)
             currentCohort%pft = this%pft(countCohort)
             currentCohort%status_coh = this%status_coh(countCohort)
             currentCohort%isnew = ( this%isnew(countCohort) .eq. new_cohort )
             
             if (this%DEBUG) then
                write(iulog,*) 'CVTL II ',countCohort, &
                      numCohort
             endif
             
             countCohort = countCohort + 1
             
             currentCohort => currentCohort%taller
             
          enddo ! current cohort do while
          
          
          ! FIX(SPM,032414) move to init if you can...or make a new init function
          currentPatch%leaf_litter(:)    = 0.0_r8
          currentPatch%root_litter(:)    = 0.0_r8
          currentPatch%leaf_litter_in(:) = 0.0_r8
          currentPatch%root_litter_in(:) = 0.0_r8
          currentPatch%seed_bank(:)      = 0.0_r8
          currentPatch%spread(:)         = 0.0_r8
          
          !
          ! deal with patch level fields here
          !
          currentPatch%livegrass  = this%livegrass(incrementOffset)
          currentPatch%age        = this%age(incrementOffset) 
          currentPatch%area       = this%areaRestart(incrementOffset) 
          
          
          
          ! set cohorts per patch for IO
          
          if (this%DEBUG) then
             write(iulog,*) 'CVTL III ' &
                   ,countCohort,cohorts_per_col, numCohort
          endif
          !
          ! deal with patch level fields of arrays here
          !
          ! these are arrays of length numpft_ed, each patch contains one
          ! vector so we increment 
          do i = 1,numpft_ed  ! numpft_ed currently 2
             currentPatch%leaf_litter(i)    = this%leaf_litter(countPft)    
             currentPatch%root_litter(i)    = this%root_litter(countPft)    
             currentPatch%leaf_litter_in(i) = this%leaf_litter_in(countPft) 
             currentPatch%root_litter_in(i) = this%root_litter_in(countPft) 
             currentPatch%seed_bank(i)      = this%seed_bank(countPft) 
             countPft = countPft + 1
          end do
          
          do i = 1,ncwd ! ncwd currently 4
             currentPatch%cwd_ag(i) = this%cwd_ag(countNcwd)
             currentPatch%cwd_bg(i) = this%cwd_bg(countNcwd)
             countNcwd = countNcwd + 1
          end do
          
          do i = 1,cp_nclmax ! cp_nclmax currently 2
             currentPatch%spread(i) = this%spread(countNclmax) 
             countNclmax  = countNclmax + 1
          end do
          
          if (this%DEBUG) write(iulog,*) 'CVTL countSunZ 1 ',countSunZ
          
          do k = 1,cp_nlevcan ! cp_nlevcan currently 40
             do j = 1,numpft_ed ! numpft_ed currently 2
                do i = 1,cp_nclmax ! cp_nclmax currently 2
                   currentPatch%f_sun(i,j,k)      = this%f_sun(countSunZ) 
                   currentPatch%fabd_sun_z(i,j,k) = this%fabd_sun_z(countSunZ) 
                   currentPatch%fabi_sun_z(i,j,k) = this%fabi_sun_z(countSunZ) 
                   currentPatch%fabd_sha_z(i,j,k) = this%fabd_sha_z(countSunZ) 
                   currentPatch%fabi_sha_z(i,j,k) = this%fabi_sha_z(countSunZ) 
                   countSunZ = countSunZ + 1
                end do
             end do
          end do
          
          if (this%DEBUG) write(iulog,*) 'CVTL countSunZ 2 ',countSunZ
          
          incrementOffset = incrementOffset + numCohortsPerPatch

          ! and the number of allowed cohorts per patch (currently 200)
          countPft      = incrementOffset
          countNcwd     = incrementOffset
          countNclmax   = incrementOffset
          countCohort   = incrementOffset
          countSunZ     = incrementOffset
          
          if (this%DEBUG) then
             write(iulog,*) 'CVTL incrementOffset ', incrementOffset
             write(iulog,*) 'CVTL cohorts_per_col ', cohorts_per_col
             write(iulog,*) 'CVTL numCohort ', numCohort
             write(iulog,*) 'CVTL totalCohorts ', totalCohorts
          end if
          
          currentPatch => currentPatch%younger
          
       enddo ! currentPatch do while
       
       do i = 1,numWaterMem
          sites(s)%water_memory(i) = this%water_memory( countWaterMem )
          countWaterMem = countWaterMem + 1
       end do
       
       sites(s)%old_stock      = this%old_stock(c)
       sites(s)%status         = this%cd_status(c)
       sites(s)%dstatus        = this%dd_status(c)
       sites(s)%ncd            = this%ncd(c)
       sites(s)%leafondate     = this%leafondate(c)
       sites(s)%leafoffdate    = this%leafoffdate(c)
       sites(s)%dleafondate    = this%dleafondate(c)
       sites(s)%dleafoffdate   = this%dleafoffdate(c)
       sites(s)%acc_NI         = this%acc_NI(c)
       sites(s)%ED_GDD_site    = this%ED_GDD_site(c)
       
    enddo

    if (this%DEBUG) then
       write(iulog,*) 'CVTL total cohorts ',totalCohorts
    end if

  end subroutine convertCohortVectorToList

  !--------------------------------------------!
  ! Non Type-Bound Procedures Here:
  !--------------------------------------------!

  !-------------------------------------------------------------------------------!
  subroutine EDRest ( bounds, nsites, sites, fcolumn, ncid, flag )
    !
    ! !DESCRIPTION:
    ! Read/write ED restart data
    ! EDRest called from restFileMod.F90
    !
    ! !USES:

    use ncdio_pio    , only : file_desc_t
    use EDCLMLinkMod , only : ed_clm_type
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)            :: bounds  ! bounds
    type(file_desc_t)       , intent(inout)         :: ncid    ! netcdf id
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)   ! The site vector
    integer                 , intent(in)            :: fcolumn(nsites)
    character(len=*)        , intent(in)            :: flag    !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    type(EDRestartVectorClass) :: ervc

    !-----------------------------------------------------------------------
    !
    ! Note: ed_allsites_inst already exists and is allocated in clm_instInit
    !

    ervc = newEDRestartVectorClass( bounds )

    if (ervc%DEBUG) then
        write(iulog,*) 'EDRestVectorMod:EDRest flag ',flag
    end if

    if ( flag == 'write' ) then
       call ervc%setVectors( bounds, nsites, sites, fcolumn )
    endif

    call ervc%doVectorIO( ncid, flag )

    if ( flag == 'read' ) then
       call ervc%getVectors( bounds, nsites, sites, fcolumn )
    endif

    call ervc%deleteEDRestartVectorClass ()

  end subroutine EDRest

end module EDRestVectorMod
