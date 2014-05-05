module EDRestVectorMod

   use shr_kind_mod        , only : r8 => shr_kind_r8
   use shr_assert_mod      , only : assert => shr_assert 
   use shr_log_mod         , only : errMsg => shr_log_errMsg
   use shr_sys_mod         , only : shr_sys_abort
   use EDtypesMod          , only : site, patch, cohort, ncwd, invalidValue, gridcell_edstate_type
   use EDtypesMod          , only : AREA, cohorts_per_gcell, numpft_ed, numWaterMem, nclmax

   implicit none
   save
   private

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
      integer,  pointer :: cellWithPatch(:)
      integer,  pointer :: numPatchesPerCell(:)
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
      integer,  pointer :: pft(:) 
      integer,  pointer :: status(:)
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
      !
      ! site level restart vars
      !
      real(r8), pointer :: water_memory(:) 
      real(r8), pointer :: old_stock(:) 
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
   end interface

   ! 
   ! non type-bound procedures
   !
   public :: EDRest

contains

   !--------------------------------------------!
   ! Type-Bound Procedures Here:
   !--------------------------------------------!

   !
   ! provide clean-up routine of allocated pointer arrays
   !
   subroutine deleteEDRestartVectorClass( this )

      class(EDRestartVectorClass), intent(inout)      :: this

      deallocate(this%cellWithPatch )
      deallocate(this%numPatchesPerCell )
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
      deallocate(this%pft )
      deallocate(this%status )
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
      deallocate(this%water_memory )
      deallocate(this%old_stock )

   end subroutine deleteEDRestartVectorClass

   !
   ! provide user-defined ctor, with array length argument
   ! allocate memory for vector to write
   !
   function newEDRestartVectorClass( bounds )

      use decompMod     , only : bounds_type

      implicit none

      type(bounds_type) , intent(in)  :: bounds ! bounds

      type(EDRestartVectorClass) :: newEDRestartVectorClass
      integer                    :: retVal = 99
      integer, parameter         :: allocOK = 0

      associate( new => newEDRestartVectorClass)

      ! set class variables
      new%vectorLengthStart = bounds%begCohort
      new%vectorLengthStop  = bounds%endCohort

      ! 
      ! cohort level variables that are required on restart
      !

      allocate(new%cellWithPatch &
           (bounds%begg:bounds%endg), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%cellWithPatch(:) = 0

      allocate(new%numPatchesPerCell &
           (bounds%begg:bounds%endg), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%numPatchesPerCell(:) = invalidValue

      allocate(new%cohortsPerPatch &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%cohortsPerPatch(:) = invalidValue

      allocate(new%balive &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%balive(:) = 0_r8

      allocate(new%bdead &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%bdead(:) = 0_r8

      allocate(new%bl &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%bl(:) = 0_r8

      allocate(new%br &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%br(:) = 0_r8

      allocate(new%bstore &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%bstore(:) = 0_r8

      allocate(new%canopy_layer &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%canopy_layer(:) = 0_r8

      allocate(new%canopy_trim &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%canopy_trim(:) = 0_r8

      allocate(new%dbh &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%dbh(:) = 0_r8

      allocate(new%hite &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%hite(:) = 0_r8

      allocate(new%laimemory &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%laimemory(:) = 0_r8

      allocate(new%leaf_md &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%leaf_md(:) = 0_r8

      allocate(new%root_md &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%root_md(:) = 0_r8

      allocate(new%n &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%n(:) = 0_r8

      allocate(new%gpp_acc &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%gpp_acc(:) = 0_r8

      allocate(new%npp_acc &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%npp_acc(:) = 0_r8

      allocate(new%pft &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%pft(:) = 0

      allocate(new%status &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%status(:) = 0

      ! 
      ! some patch level variables that are required on restart
      !
      allocate(new%cwd_ag &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%cwd_ag(:) = 0_r8

      allocate(new%cwd_bg &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%cwd_bg(:) = 0_r8

      allocate(new%leaf_litter &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%leaf_litter(:) = 0_r8

      allocate(new%root_litter &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%root_litter(:) = 0_r8

      allocate(new%leaf_litter_in &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%leaf_litter_in(:) = 0_r8

      allocate(new%root_litter_in &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%root_litter_in(:) = 0_r8

      allocate(new%seed_bank &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%seed_bank(:) = 0_r8

      allocate(new%spread &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%spread(:) = 0_r8

      allocate(new%livegrass &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%livegrass(:) = 0_r8

      allocate(new%age &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%age(:) = 0_r8

      allocate(new%areaRestart &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%areaRestart(:) = 0_r8

      !
      ! site level variable
      !

      allocate(new%water_memory &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%water_memory(:) = 0_r8

      allocate(new%old_stock &
           (new%vectorLengthStart:new%vectorLengthStop), stat=retVal)
      call assert(( retVal == allocOK ), errMsg(__FILE__, __LINE__))
      new%old_stock(:) = 0_r8

      end associate

   end function newEDRestartVectorClass

   ! 
   ! implement setVectors
   !
   subroutine setVectors( this, bounds, geds_local )

      use decompMod,         only : bounds_type
      use clm_time_manager , only : get_nstep
      use clm_varctl,        only : iulog

      implicit none

      class(EDRestartVectorClass), intent(inout)      :: this
      type(bounds_type)          , intent(in)         :: bounds 
      type(gridcell_edstate_type), target, intent(in) :: geds_local( bounds%begg: )

      write(iulog,*) 'edtime setVectors ',get_nstep()

      if (this%DEBUG) then
         call this%printIoInfoLL ( bounds, geds_local )
         call this%printDataInfoLL ( bounds, geds_local )
      end if

      call this%convertCohortListToVector ( bounds, geds_local )

      if (this%DEBUG) then
         call this%printDataInfoVector (  )
      end if

   end subroutine setVectors

   ! 
   ! implement getVectors
   !
   subroutine getVectors( this, bounds, geds_local )

      use decompMod ,        only : bounds_type
      use clm_time_manager , only : get_nstep
      use clm_varctl       , only : iulog

      use EDInitMod ,        only : ed_init_sites
      use EDCLMLinkMod,      only : clm_ed_link
      use EDMainMod,         only : ed_update_sites

      implicit none

      class(EDRestartVectorClass), intent(inout)         :: this
      type(bounds_type)          , intent(in)            :: bounds 
      type(gridcell_edstate_type), target, intent(inout) :: geds_local( bounds%begg: )

      if (this%DEBUG) then
         write(iulog,*) 'edtime getVectors ',get_nstep()
         call this%printDataInfoVector (  )
      end if

      call ed_init_sites( bounds, geds_local )

      call this%createPatchCohortStructure ( bounds, geds_local )
      call this%convertCohortVectorToList ( bounds, geds_local )

      call ed_update_sites( bounds, geds_local )
      call clm_ed_link( bounds, geds_local )

      if (this%DEBUG) then
         call this%printIoInfoLL ( bounds, geds_local )
      end if

   end subroutine getVectors

   ! 
   ! implement doVectorIO
   !
   subroutine doVectorIO( this, ncid, flag  )

      use ncdio_pio  , only : file_desc_t, ncd_int, ncd_double
      use restUtilMod, only : restartvar
      use clmType,     only : nameg, nameCohort
      use decompMod,   only : get_clmlevel_gsmap
      use spmdMod,     only : iam
      use mct_mod,     only : mct_gsMap, mct_gsmap_OP

      implicit none

      class(EDRestartVectorClass), intent(inout) :: this
      type(file_desc_t), intent(inout)           :: ncid   ! netcdf id
      character(len=*) , intent(in)              :: flag   !'read' or 'write'

      logical             :: readvar
      character(len=16)   :: dimName  = trim(nameCohort)

      type(mct_gsMap),pointer       :: gsmap   ! global seg map
      integer, pointer,dimension(:) :: gsmOP   ! gsmap ordered points

      call get_clmlevel_gsmap(clmlevel='cohort',gsmap=gsmap)
      call mct_gsmap_OP(gsmap,iam,gsmOP)

      !
      ! cohort level vars
      !

      !===========!


      call restartvar(ncid=ncid, flag=flag, varname='ed_io_cellWithPatch', xtype=ncd_int,  &
           dim1name=nameg, &
           long_name='1 if a gridcell has a patch', units='1=true,0=false', &
           interpinic_flag='interp', data=this%cellWithPatch, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_io_numPatchesPerCell', xtype=ncd_int,  &
           dim1name=nameg, &
           long_name='works with ed_cellWithPatch.  num patches per gridcell', units='unitless', &
           interpinic_flag='interp', data=this%numPatchesPerCell, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_io_cohortsPerPatch', xtype=ncd_int,  &
           dim1name=dimName, &
           long_name='list of cohorts per patch.  indexed by numPatchesPerCell', units='unitless', &
           interpinic_flag='interp', data=this%cohortsPerPatch, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_balive', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort ed_balive', units='unitless', &
           interpinic_flag='interp', data=this%balive, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_bdead', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - bdead', units='unitless', &
           interpinic_flag='interp', data=this%bdead, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_bl', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - bl', units='unitless', &
           interpinic_flag='interp', data=this%bl, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_br', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - br', units='unitless', &
           interpinic_flag='interp', data=this%br, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_bstore', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - bstore', units='unitless', &
           interpinic_flag='interp', data=this%bstore, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_canopy_layer', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - canopy_layer', units='unitless', &
           interpinic_flag='interp', data=this%canopy_layer, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_canopy_trim', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - canopy_trim', units='unitless', &
           interpinic_flag='interp', data=this%canopy_trim, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_dbh', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - dbh', units='unitless', &
           interpinic_flag='interp', data=this%dbh, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_hite', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - hite', units='unitless', &
           interpinic_flag='interp', data=this%hite, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_laimemory', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - laimemory', units='unitless', &
           interpinic_flag='interp', data=this%laimemory, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_leaf_md', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - leaf_md', units='unitless', &
           interpinic_flag='interp', data=this%leaf_md, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_root_md', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - root_md', units='unitless', &
           interpinic_flag='interp', data=this%root_md, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_n', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - n', units='unitless', &
           interpinic_flag='interp', data=this%n, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_gpp_acc', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - gpp_acc', units='unitless', &
           interpinic_flag='interp', data=this%gpp_acc, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_npp_acc', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - npp_acc', units='unitless', &
           interpinic_flag='interp', data=this%npp_acc, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_pft', xtype=ncd_int,  &
           dim1name=dimName, &
           long_name='ed cohort - pft', units='unitless', &
           interpinic_flag='interp', data=this%pft, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_status', xtype=ncd_int,  &
           dim1name=dimName, &
           long_name='ed cohort - status', units='unitless', &
           interpinic_flag='interp', data=this%status, &
           readvar=readvar)

      !
      ! patch level vars
      !

      call restartvar(ncid=ncid, flag=flag, varname='ed_cwd_ag', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - cwd_ag', units='unitless', &
           interpinic_flag='interp', data=this%cwd_ag, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_cwd_bg', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - cwd_bg', units='unitless', &
           interpinic_flag='interp', data=this%cwd_bg, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_leaf_litter', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - leaf_litter', units='unitless', &
           interpinic_flag='interp', data=this%leaf_litter, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_root_litter', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - root_litter', units='unitless', &
           interpinic_flag='interp', data=this%root_litter, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_leaf_litter_in', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - leaf_litter_in', units='unitless', &
           interpinic_flag='interp', data=this%leaf_litter_in, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_root_litter_in', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - root_litter_in', units='unitless', &
           interpinic_flag='interp', data=this%root_litter_in, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_seed_bank', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - seed_bank', units='unitless', &
           interpinic_flag='interp', data=this%seed_bank, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_spread', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - spread', units='unitless', &
           interpinic_flag='interp', data=this%spread, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_livegrass', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - livegrass', units='unitless', &
           interpinic_flag='interp', data=this%livegrass, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_age', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - age', units='unitless', &
           interpinic_flag='interp', data=this%age, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_area', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - area', units='unitless', &
           interpinic_flag='interp', data=this%areaRestart, &
           readvar=readvar)

      !
      ! site level vars
      !

      call restartvar(ncid=ncid, flag=flag, varname='ed_water_memory', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - water_memory', units='unitless', &
           interpinic_flag='interp', data=this%water_memory, &
           readvar=readvar)

      call restartvar(ncid=ncid, flag=flag, varname='ed_old_stock', xtype=ncd_double,  &
           dim1name=dimName, &
           long_name='ed cohort - old_stock', units='unitless', &
           interpinic_flag='interp', data=this%old_stock, &
           readvar=readvar)

   end subroutine doVectorIO

  ! ============================================================================
  !
  ! ============================================================================
   subroutine printDataInfoVector( this )

      use clm_varctl    , only : iulog

      implicit none

      class(EDRestartVectorClass), intent(inout) :: this

      character(len=32)   :: methodName = 'PDIV '
      integer :: iSta, iSto

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

      write(iulog,*) trim(methodName)//' :: pft ', &
           this%pft(iSta:iSto)
      write(iulog,*) trim(methodName)//' :: status ', &
           this%status(iSta:iSto)

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
      write(iulog,*) trim(methodName)//' :: water_memory ', &
           this%water_memory(iSta:iSto)
      write(iulog,*) trim(methodName)//' :: old_stock ', &
           this%old_stock(iSta:iSto)

   end subroutine printDataInfoVector

  ! ============================================================================
  !
  ! ============================================================================
  subroutine printDataInfoLL( this, bounds, geds_local ) 

      !
      ! counts the total number of cohorts over all p levels (patch) so we
      ! can allocate vectors, copy from LL -> vector and read/write restarts.
      !

      use decompMod,  only : bounds_type
      use clm_varctl, only : iulog

      implicit none

      class(EDRestartVectorClass), intent(inout)      :: this
      type(bounds_type),           intent(in)         :: bounds 
      type(gridcell_edstate_type), target, intent(in) :: geds_local( bounds%begg: )

      type (site) ,  pointer :: currentSite
      type (patch),  pointer :: currentPatch
      type (cohort), pointer :: currentCohort

      integer g
      integer totalCohorts
      integer numCohort
      integer numPatches,totPatchCount

      character(len=32)   :: methodName = 'printDataInfoLL '

      totalCohorts = 0
      totPatchCount = 1

      write(iulog,*) 'vecLenStart ',this%vectorLengthStart

      g = bounds%begg
      do while(g <= bounds%endg)
         currentSite => geds_local(g)%spnt

         if(currentSite%istheresoil == 1)then
            currentPatch => currentSite%oldest_patch

            numPatches = 1

            do while(associated(currentPatch))
               currentCohort => currentPatch%shortest

               write(iulog,*) trim(methodName)//':: found gcell with patch(s) ',g

               numCohort = 0

               do while(associated(currentCohort))  

                  totalCohorts = totalCohorts + 1

                  write(iulog,*) trim(methodName)//' balive ',totalCohorts,currentCohort%balive
                  write(iulog,*) trim(methodName)//' bdead ',totalCohorts,currentCohort%bdead
                  write(iulog,*) trim(methodName)//' bl ',totalCohorts,currentCohort%bl
                  write(iulog,*) trim(methodName)//' br ',totalCohorts,currentCohort%br
                  write(iulog,*) trim(methodName)//' bstore ',totalCohorts,currentCohort%bstore
                  write(iulog,*) trim(methodName)//' canopy_layer ',totalCohorts,currentCohort%canopy_layer
                  write(iulog,*) trim(methodName)//' canopy_trim ',totalCohorts,currentCohort%canopy_trim
                  write(iulog,*) trim(methodName)//' dbh ',totalCohorts,currentCohort%dbh
                  write(iulog,*) trim(methodName)//' hite ',totalCohorts,currentCohort%hite
                  write(iulog,*) trim(methodName)//' laimemory ',totalCohorts,currentCohort%laimemory
                  write(iulog,*) trim(methodName)//' leaf_md ',totalCohorts,currentCohort%leaf_md
                  write(iulog,*) trim(methodName)//' root_md ',totalCohorts,currentCohort%root_md
                  write(iulog,*) trim(methodName)//' n ',totalCohorts,currentCohort%n
                  write(iulog,*) trim(methodName)//' gpp_acc ', &
                       totalCohorts,currentCohort%gpp_acc
                  write(iulog,*) trim(methodName)//' npp_acc ', &
                       totalCohorts,currentCohort%npp_acc
                  write(iulog,*) trim(methodName)//' pft ',totalCohorts,currentCohort%pft
                  write(iulog,*) trim(methodName)//' status ',totalCohorts,currentCohort%status

                  numCohort = numCohort + 1

                  currentCohort => currentCohort%taller
               enddo ! currentCohort do while

               write(iulog,*) trim(methodName)//': numpatches for gcell ',currentSite%clmgcell, numPatches
               write(iulog,*) trim(methodName)//': patches and cohorts ',totPatchCount,numCohort

               write(iulog,*) trim(methodName)//' cwd_ag ',currentPatch%cwd_ag
               write(iulog,*) trim(methodName)//' cwd_bg ',currentPatch%cwd_bg
               write(iulog,*) trim(methodName)//' leaf_litter ',currentPatch%leaf_litter
               write(iulog,*) trim(methodName)//' root_litter ',currentPatch%root_litter
               write(iulog,*) trim(methodName)//' leaf_litter_in ',currentPatch%leaf_litter_in
               write(iulog,*) trim(methodName)//' root_litter_in ',currentPatch%root_litter_in
               write(iulog,*) trim(methodName)//' seed_bank ',currentPatch%seed_bank
               write(iulog,*) trim(methodName)//' spread ',currentPatch%spread
               write(iulog,*) trim(methodName)//' livegrass ',currentPatch%livegrass
               write(iulog,*) trim(methodName)//' age ',currentPatch%age
               write(iulog,*) trim(methodName)//' area ',currentPatch%area
               write(iulog,*) trim(methodName)//' old_stock ',currentSite%old_stock

               currentPatch => currentPatch%younger

               totPatchCount = totPatchCount + 1
               numPatches = numPatches + 1
            enddo ! currentPatch do while
         endif
         g = g + 1

         write(iulog,*) trim(methodName)//' water_memory ',currentSite%water_memory(1)

      enddo

      write(iulog,*) trim(methodName)//': total cohorts ',totalCohorts

   end subroutine printDataInfoLL

  ! ============================================================================
  !
  ! ============================================================================
  subroutine printIoInfoLL( this, bounds, geds_local ) 

      !
      ! counts the total number of cohorts over all p levels (patch) so we
      ! can allocate vectors, copy from LL -> vector and read/write restarts.
      !

      use decompMod,   only : bounds_type
      use clm_varctl , only : iulog

      implicit none

      class(EDRestartVectorClass), intent(inout)      :: this
      type(bounds_type),           intent(in)         :: bounds 
      type(gridcell_edstate_type), target, intent(in) :: geds_local( bounds%begg: )

      type (site) ,  pointer :: currentSite
      type (patch),  pointer :: currentPatch
      type (cohort), pointer :: currentCohort

      integer g
      integer totalCohorts
      integer numCohort
      integer numPatches,totPatchCount

      character(len=32)   :: methodName = 'printIoInfoLL '

      totalCohorts = 0
      totPatchCount = 1

      write(iulog,*) 'vecLenStart ',this%vectorLengthStart

      g = bounds%begg
      do while(g <= bounds%endg)
         currentSite => geds_local(g)%spnt

         if(currentSite%istheresoil == 1)then
            currentPatch => currentSite%oldest_patch

            numPatches = 1

            do while(associated(currentPatch))
               currentCohort => currentPatch%shortest

               write(iulog,*) trim(methodName)//': found gcell with patch(s) ',g

               numCohort = 0

               do while(associated(currentCohort))  

                  totalCohorts = totalCohorts + 1
                  numCohort = numCohort + 1

                  currentCohort => currentCohort%taller
               enddo ! currentCohort do while

               write(iulog,*) trim(methodName)//': numpatches for gcell ',currentSite%clmgcell, numPatches
               write(iulog,*) trim(methodName)//': patches and cohorts ',totPatchCount,numCohort

               currentPatch => currentPatch%younger

               totPatchCount = totPatchCount + 1
               numPatches = numPatches + 1
            enddo ! currentPatch do while
         endif
         g = g + 1
      enddo

   end subroutine printIoInfoLL

  ! ============================================================================
  !
  ! ============================================================================
  subroutine convertCohortListToVector( this, bounds, geds_local ) 

      !
      ! counts the total number of cohorts over all p levels (patch) so we
      ! can allocate vectors, copy from LL -> vector and read/write restarts.
      !

      use decompMod,  only : bounds_type
      use clm_varpar, only : nclmax
      use clm_varctl, only : iulog

      implicit none

      class(EDRestartVectorClass), intent(inout)      :: this
      type(bounds_type),           intent(in)         :: bounds 
      type(gridcell_edstate_type), target, intent(in) :: geds_local( bounds%begg: )

      type (site) , pointer  :: currentSite
      type (patch), pointer  :: currentPatch
      type (cohort), pointer :: currentCohort

      integer ::  g
      integer ::  totalCohorts ! number of cohorts starting from 1
      integer ::  offsetNumCohorts ! number of cohorts starting from
                               ! vectorLengthStart
      integer :: numCohort
      integer :: numPatches
      integer :: totPatchCount, offsetTotPatchCount
      integer :: countPft
      integer :: countNcwd
      integer :: countWaterMem
      integer :: countNclmax
      integer :: i

      totalCohorts = 1
      offsetNumCohorts = this%vectorLengthStart

      totPatchCount = 1
      offsetTotPatchCount = this%vectorLengthStart

      countPft = this%vectorLengthStart

      countNcwd = this%vectorLengthStart

      countNclmax = this%vectorLengthStart

      countWaterMem=this%vectorLengthStart

      g = bounds%begg
      do while(g <= bounds%endg)

         currentSite => geds_local(g)%spnt

         if(currentSite%istheresoil == 1)then
            currentPatch => currentSite%oldest_patch

            numPatches = 1

            do while(associated(currentPatch))
               currentCohort => currentPatch%shortest

               numCohort = 0

               do while(associated(currentCohort))
              
                  if ( offsetNumCohorts  >  this%vectorLengthStop ) then
                     call shr_sys_abort( 'error in convertCohortListToVector ::&
                     overrun of number of cohorts.  Try increasing cohorts for&
                     IO' )               
                  endif

                  this%balive(offsetNumCohorts) = currentCohort%balive
                  this%bdead(offsetNumCohorts) = currentCohort%bdead
                  this%bl(offsetNumCohorts) = currentCohort%bl
                  this%br(offsetNumCohorts) = currentCohort%br
                  this%bstore(offsetNumCohorts) = currentCohort%bstore
                  this%canopy_layer(offsetNumCohorts) = currentCohort%canopy_layer
                  this%canopy_trim(offsetNumCohorts) = currentCohort%canopy_trim
                  this%dbh(offsetNumCohorts) = currentCohort%dbh
                  this%hite(offsetNumCohorts) = currentCohort%hite
                  this%laimemory(offsetNumCohorts) = currentCohort%laimemory
                  this%leaf_md(offsetNumCohorts) = currentCohort%leaf_md
                  this%root_md(offsetNumCohorts) = currentCohort%root_md
                  this%n(offsetNumCohorts) = currentCohort%n
                  this%gpp_acc(offsetNumCohorts) = currentCohort%gpp_acc
                  this%npp_acc(offsetNumCohorts) = currentCohort%npp_acc
                  this%pft(offsetNumCohorts) = currentCohort%pft
                  this%status(offsetNumCohorts) = currentCohort%status

                  numCohort = numCohort + 1
                  totalCohorts = totalCohorts + 1
                  offsetNumCohorts = offsetNumCohorts + 1

                  currentCohort => currentCohort%taller

               enddo ! currentCohort do while

               offsetNumCohorts = offsetNumCohorts + &
                  cohorts_per_gcell - numCohort

               ! set numpatches for this gcell
               this%numPatchesPerCell( currentSite%clmgcell )  = numPatches
               ! set cohorts per patch
               this%cohortsPerPatch(offsetTotPatchCount) = numCohort

               ! deal with patch level field here...required for restart

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

               do i = 1,ncwd
                  this%cwd_ag(countNcwd) = currentPatch%cwd_ag(i)
                  this%cwd_bg(countNcwd) = currentPatch%cwd_bg(i)
                  countNcwd = countNcwd + 1
               end do

               do i = 1,nclmax
                  this%spread(countNclmax)         = currentPatch%spread(i)
                  countNclmax = countNclmax + 1
               end do

               this%livegrass(offsetTotPatchCount)   = currentPatch%livegrass
               this%age(offsetTotPatchCount)         = currentPatch%age
               this%areaRestart(offsetTotPatchCount) = currentPatch%area
               this%old_stock(offsetTotPatchCount)   = currentSite%old_stock

               currentPatch => currentPatch%younger

               totPatchCount = totPatchCount + 1
               offsetTotPatchCount = offsetTotPatchCount + cohorts_per_gcell
               numPatches = numPatches + 1
            enddo ! currentPatch do while

            ! set which gridcells have patches/cohorts
            this%cellWithPatch( currentSite%clmgcell )  = 1

            g = g + 1

            do i = 1,numWaterMem
              this%water_memory( countWaterMem ) = currentSite%water_memory(i)
              countWaterMem = countWaterMem + 1
            end do

            countWaterMem = countWaterMem + cohorts_per_gcell - numWaterMem 
            countNcwd = countNcwd + cohorts_per_gcell - ncwd 

         endif ! is there soil check

      enddo

      if (this%DEBUG) then
         write(iulog,*) 'total cohorts ',totalCohorts
      end if

   end subroutine convertCohortListToVector

  ! ============================================================================
  !
  ! ============================================================================
  subroutine createPatchCohortStructure( this, bounds, geds_local ) 

      !
      ! counts the total number of cohorts over all p levels (patch) so we
      ! can allocate vectors, copy from LL -> vector and read/write restarts.
      !

      use decompMod          ,  only : bounds_type
      use EDPatchDynamicsMod ,  only : zero_patch
      use EDGrowthFunctionsMod, only : Dbh
      use EDCohortDynamicsMod,  only : create_cohort
      use EDParamsMod        ,  only : ED_val_maxspread
      use EDPatchDynamicsMod ,  only : create_patch
      use clmtype            ,  only : pftcon
      use clm_varctl         ,  only : iulog

      implicit none

      class(EDRestartVectorClass), intent(inout)         :: this
      type(bounds_type),           intent(in)            :: bounds 
      type(gridcell_edstate_type), target, intent(inout) :: geds_local( bounds%begg: )

      type (site)  , pointer  :: currentSite
      type (patch) , pointer  :: newp

      type(cohort), allocatable :: dc

      real(r8) cwd_ag_local(ncwd),cwd_bg_local(ncwd),spread_local(nclmax)
      real(r8) leaf_litter_local(numpft_ed),root_litter_local(numpft_ed)
      real(r8) seed_bank_local(numpft_ed)
      real(r8) age !notional age of this patch

      integer :: cohortstatus
      integer :: g,patchIdx,currIdx, fto, ft

      currIdx = this%vectorLengthStart

      cwd_ag_local = 0.0_r8 !ED_val_init_litter   !arbitrary value for litter pools. kgC m-2           ! 
      cwd_bg_local = 0.0_r8 !ED_val_init_litter
      leaf_litter_local = 0.0_r8
      root_litter_local = 0.0_r8
      spread_local = ED_val_maxspread
      age = 0.0_r8
       
      ! default edinit creates one patch with 2 cohorts.  take this into account
      do g = bounds%begg, bounds%endg

         if (this%DEBUG) then
            write(iulog,*) 'cellWithPatch ',this%cellWithPatch(g),this%numPatchesPerCell(g)
         end if

         currentSite => geds_local(g)%spnt
            
           do patchIdx = 1,this%numPatchesPerCell(g)

              if (this%DEBUG) then
                 write(iulog,*) 'create patch ',patchIdx
                 write(iulog,*) 'patchIdx 1-numCohorts : ',this%cohortsPerPatch(currIdx)
              end if

              ! create patch, then create cohorts here.  1st test is with
              ! default 1 patch+2 cohorts

              allocate(newp)    
              call zero_patch(newp)

              currentSite%youngest_patch => newp 
              !
              ! make new patch
              !
              ! FIX(SPM,032414) need to handle the creation of multiple patches for restarts
              call create_patch(currentSite,newp,age,AREA,spread_local,cwd_ag_local,cwd_bg_local, &
                   leaf_litter_local,root_litter_local,seed_bank_local) 

              newp%patchno = 1

              ! FIX(SPM,032414) need to handle more than 2 cohorts per patch for restarts
              do fto = 1,this%cohortsPerPatch(currIdx)

                 allocate(dc)

                 dc%n = 700_r8
                 dc%balive = 0_r8
                 dc%bdead = 0_r8
                 dc%bstore = 0_r8
                 dc%laimemory = 0_r8
                 dc%canopy_trim = 0_r8
                 dc%canopy_layer = 1_r8

                 ! set the pft (only 2 used in ed) based on odd/even cohort
                 ! number
                 ft=2
                 if ((mod(fto,2)  ==  0 )) then
                    ft=1
                 endif

                 cohortstatus = newp%siteptr%status
  
                 if(pftcon%stress_decid(ft) == 1)then !drought decidous, override status. 
                    cohortstatus = newp%siteptr%dstatus
                 endif

                 dc%hite = 1.25_r8
                 ! the dbh function should only take as an argument, the one
                 ! item it needs, not the entire cohort...refactor
                 dc%dbh = Dbh(dc) + 0.0001_r8*ft

                 call create_cohort(ft,dc%n,dc%hite,dc%dbh,dc%balive,dc%bdead,dc%bstore, &
                      dc%laimemory,cohortstatus,dc%canopy_trim,newp%NCL_p,newp)

                  deallocate(dc)

              enddo

              currentSite%youngest_patch => newp                   
              currentSite%oldest_patch => newp                        
                   
              currIdx = currIdx + cohorts_per_gcell
            
           enddo
            
       enddo

   end subroutine createPatchCohortStructure

  ! ============================================================================
  !
  ! ============================================================================
  subroutine convertCohortVectorToList( this, bounds, geds_local ) 

      !
      ! counts the total number of cohorts over all p levels (patch) so we
      ! can allocate vectors, copy from LL -> vector and read/write restarts.
      !

      use decompMod          , only : bounds_type
      use clm_varpar         , only : nclmax
      use clm_varctl         , only : iulog

      implicit none

      class(EDRestartVectorClass), intent(inout)         :: this
      type(bounds_type),           intent(in)            :: bounds 
      type(gridcell_edstate_type), target, intent(inout) :: geds_local( bounds%begg: )

      type (site) , pointer :: currentSite
      type (patch), pointer :: currentPatch
      type (cohort),pointer :: currentCohort

      integer :: g
      integer :: offsetNumCohorts
      integer :: countPft
      integer :: countNcwd
      integer :: countWaterMem
      integer :: i
      integer :: offsetTotPatchCount
      integer :: numCohort
      integer :: countNclmax

      offsetNumCohorts = this%vectorLengthStart

      offsetTotPatchCount = this%vectorLengthStart

      countPft = this%vectorLengthStart

      countNcwd=this%vectorLengthStart

      countNclmax = this%vectorLengthStart

      countWaterMem=this%vectorLengthStart

      g = bounds%begg
      do while(g <= bounds%endg)

         currentSite => geds_local(g)%spnt

         if(currentSite%istheresoil == 1)then

            if (this%DEBUG) then
               write(iulog,*) 'CCVTL I ',g,this%cellWithPatch(g)
            end if

            currentSite%clmgcell = g
            currentPatch => currentSite%oldest_patch

            do while(associated(currentPatch))
               currentCohort => currentPatch%shortest

               numCohort = 0

               do while(associated(currentCohort))        

                  currentCohort%balive = this%balive(offsetNumCohorts)
                  currentCohort%bdead = this%bdead(offsetNumCohorts)
                  currentCohort%bl = this%bl(offsetNumCohorts)
                  currentCohort%br = this%br(offsetNumCohorts)
                  currentCohort%bstore = this%bstore(offsetNumCohorts)
                  currentCohort%canopy_layer = this%canopy_layer(offsetNumCohorts)
                  currentCohort%canopy_trim = this%canopy_trim(offsetNumCohorts)
                  currentCohort%dbh = this%dbh(offsetNumCohorts)
                  currentCohort%hite = this%hite(offsetNumCohorts)
                  currentCohort%laimemory = this%laimemory(offsetNumCohorts)
                  currentCohort%leaf_md = this%leaf_md(offsetNumCohorts)
                  currentCohort%root_md = this%root_md(offsetNumCohorts)
                  currentCohort%n = this%n(offsetNumCohorts)
                  currentCohort%gpp_acc = this%gpp_acc(offsetNumCohorts)
                  currentCohort%npp_acc = this%npp_acc(offsetNumCohorts)
                  currentCohort%pft = this%pft(offsetNumCohorts)
                  currentCohort%status = this%status(offsetNumCohorts)

                  numCohort = numCohort + 1
                  offsetNumCohorts=offsetNumCohorts+1

                  currentCohort => currentCohort%taller

               enddo ! currentPatch do while

               offsetNumCohorts = offsetNumCohorts + &
                    cohorts_per_gcell - numCohort

               ! FIX(SPM,032414) move to init if you can...or make a new init function
               currentPatch%leaf_litter(:)    = 0.0_r8
               currentPatch%root_litter(:)    = 0.0_r8
               currentPatch%leaf_litter_in(:) = 0.0_r8
               currentPatch%root_litter_in(:) = 0.0_r8
               currentPatch%seed_bank(:)      = 0.0_r8
               currentPatch%spread(:)         = 0.0_r8

               do i = 1,numpft_ed
                  currentPatch%leaf_litter(i)    = this%leaf_litter(countPft)    
                  currentPatch%root_litter(i)    = this%root_litter(countPft)    
                  currentPatch%leaf_litter_in(i) = this%leaf_litter_in(countPft) 
                  currentPatch%root_litter_in(i) = this%root_litter_in(countPft) 
                  currentPatch%seed_bank(i)      = this%seed_bank(countPft) 
                  countPft = countPft + 1
               end do

               do i = 1,ncwd
                  currentPatch%cwd_ag(i) = this%cwd_ag(countNcwd)
                  currentPatch%cwd_bg(i) = this%cwd_bg(countNcwd)
                  countNcwd = countNcwd + 1
               end do

               do i = 1,nclmax
                  currentPatch%spread(i) = this%spread(countNclmax) 
                  countNclmax  = countNclmax + 1
               end do

               currentPatch%livegrass = this%livegrass(offsetTotPatchCount)
               currentPatch%age       = this%age(offsetTotPatchCount) 
               currentPatch%area      = this%areaRestart(offsetTotPatchCount) 
               currentSite%old_stock  = this%old_stock(offsetTotPatchCount)

               currentPatch => currentPatch%younger

               offsetTotPatchCount = offsetTotPatchCount + cohorts_per_gcell

            enddo
         endif

         g = g + 1

         do i = 1,numWaterMem
            currentSite%water_memory(i) = this%water_memory( countWaterMem )
            countWaterMem = countWaterMem + 1
         end do

         countWaterMem = countWaterMem + cohorts_per_gcell - numWaterMem 
         countNcwd = countNcwd + cohorts_per_gcell - ncwd 

      enddo

   end subroutine convertCohortVectorToList

   !--------------------------------------------!
   ! Non Type-Bound Procedures Here:
   !--------------------------------------------!

   !
   ! EDRest called from restFileMod.F90
   !
   subroutine EDRest ( bounds, ncid, flag )
      !
      ! Read/write ED restart data
      !
      use ncdio_pio  , only : file_desc_t
      use decompMod  , only : bounds_type
      use EDtypesMod , only : gridCellEdState

      implicit none

      !
      type(bounds_type), intent(in)    :: bounds  ! bounds
      type(file_desc_t), intent(inout) :: ncid    ! netcdf id
      character(len=*) , intent(in)    :: flag    !'read' or 'write'
      !
      type(EDRestartVectorClass) :: ervc
      !-----------------------------------------------------------------------
 
      ervc = newEDRestartVectorClass( bounds )

      if ( flag == 'write' ) then
         !
         ! gridCellEdState already exists and is allocated in ed_init
         !
         call ervc%setVectors( bounds, gridCellEdState )
      endif

      call ervc%doVectorIO( ncid, flag )

      if ( flag == 'read' ) then
         ! 
         ! allocate gridcell ed site structure at the top level for restart runs
         !
         allocate (gridCellEdState(bounds%begg:bounds%endg))

         call ervc%getVectors( bounds, gridCellEdState )

      endif

      call ervc%deleteEDRestartVectorClass ()

   end subroutine EDRest

end module EDRestVectorMod
