module CanopyStateType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_infnan_mod  , only : nan => shr_infnan_nan, shr_infnan_isnan, assignment(=)
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun
  use decompMod       , only : bounds_type
  use landunit_varcon , only : istsoil, istcrop
  use clm_varpar      , only : nlevcan, nvegwcs
  use clm_varcon      , only : spval
  use clm_varctl      , only : iulog, use_cn, use_fates, use_fates_sp, use_hydrstress
  use LandunitType    , only : lun
  use PatchType       , only : patch
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: CanopyState_type

     integer  , pointer :: frac_veg_nosno_patch     (:)   ! patch fraction of vegetation not covered by snow (0 OR 1) [-]
     integer  , pointer :: frac_veg_nosno_alb_patch (:)   ! patch fraction of vegetation not covered by snow (0 OR 1) [-]

     real(r8) , pointer :: tlai_patch               (:)   ! patch canopy one-sided leaf area index, no burying by snow
     real(r8) , pointer :: tsai_patch               (:)   ! patch canopy one-sided stem area index, no burying by snow
     real(r8) , pointer :: elai_patch               (:)   ! patch canopy one-sided leaf area index with burying by snow
     real(r8) , pointer :: esai_patch               (:)   ! patch canopy one-sided stem area index with burying by snow

     real(r8) , pointer :: tlai_input_patch          (:)   ! patch canopy one-sided leaf area index driver data for SP mode (no burying by snow)
     real(r8) , pointer :: tsai_input_patch          (:)   ! patch canopy one-sided stem area index driver data for SP mode (no burying by snow)
     real(r8) , pointer :: htop_input_patch          (:)   ! patch canopy height driver data for SP mode
     real(r8) , pointer :: hbot_input_patch          (:)   ! patch canopy bottom driver data for SP mode
     
     real(r8) , pointer :: elai240_patch            (:)   ! patch canopy one-sided leaf area index with burying by snow average over 10days
     real(r8) , pointer :: laisun_patch             (:)   ! patch patch sunlit projected leaf area index
     real(r8) , pointer :: laisha_patch             (:)   ! patch patch shaded projected leaf area index
     real(r8) , pointer :: laisun_z_patch           (:,:) ! patch patch sunlit leaf area for canopy layer
     real(r8) , pointer :: laisha_z_patch           (:,:) ! patch patch shaded leaf area for canopy layer
     real(r8) , pointer :: mlaidiff_patch           (:)   ! patch difference between lai month one and month two (for dry deposition of chemical tracers)
     real(r8) , pointer :: annlai_patch             (:,:) ! patch 12 months of monthly lai from input data set (for dry deposition of chemical tracers)
     real(r8) , pointer :: stem_biomass_patch       (:)   ! Aboveground stem biomass (kg/m**2)
     real(r8) , pointer :: leaf_biomass_patch       (:)   ! Aboveground leaf biomass  (kg/m**2)
     real(r8) , pointer :: htop_patch               (:)   ! patch canopy top (m)
     real(r8) , pointer :: hbot_patch               (:)   ! patch canopy bottom (m)
     real(r8) , pointer :: z0m_patch                (:)   ! patch momentum roughness length (m)
     real(r8) , pointer :: displa_patch             (:)   ! patch displacement height (m)
     real(r8) , pointer :: fsun_patch               (:)   ! patch sunlit fraction of canopy
     real(r8) , pointer :: fsun24_patch             (:)   ! patch 24hr average of sunlit fraction of canopy
     real(r8) , pointer :: fsun240_patch            (:)   ! patch 240hr average of sunlit fraction of canopy

     real(r8) , pointer :: dleaf_patch              (:)   ! patch characteristic leaf width (diameter) [m]
                                                          ! for non-ED/FATES this is the same as pftcon%dleaf()
     real(r8) , pointer :: rscanopy_patch           (:)   ! patch canopy stomatal resistance (s/m) (ED specific)

     real(r8) , pointer :: vegwp_patch              (:,:) ! patch vegetation water matric potential (mm)
     real(r8) , pointer :: vegwp_ln_patch           (:,:) ! patch vegetation water matric potential at local noon (mm)
     real(r8) , pointer :: vegwp_pd_patch           (:,:) ! patch predawn vegetation water matric potential (mm)

     real(r8)           :: leaf_mr_vcm = spval            ! Scalar constant of leaf respiration with Vcmax

   contains

     procedure, public  :: Init
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, public  :: ReadNML
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: UpdateAccVars
     procedure, public  :: Restart

     procedure, public  :: SetNMLForTesting ! Set namelist for unit-testing

  end type CanopyState_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(canopystate_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

    if ( this%leaf_mr_vcm == spval ) then
       call endrun(msg="ERROR canopystate Init called before ReadNML"//errmsg(sourcefile, __LINE__))
    end if

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(canopystate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%frac_veg_nosno_patch     (begp:endp))           ; this%frac_veg_nosno_patch     (:)   = huge(1)
    allocate(this%frac_veg_nosno_alb_patch (begp:endp))           ; this%frac_veg_nosno_alb_patch (:)   = 0
    allocate(this%tlai_input_patch         (begp:endp))           ; this%tlai_input_patch         (:)   = nan
    allocate(this%tsai_input_patch         (begp:endp))           ; this%tsai_input_patch         (:)   = nan
    allocate(this%htop_input_patch         (begp:endp))           ; this%htop_input_patch         (:)   = nan
    allocate(this%hbot_input_patch         (begp:endp))           ; this%hbot_input_patch         (:)   = nan
    allocate(this%tlai_patch               (begp:endp))           ; this%tlai_patch               (:)   = nan
    allocate(this%tsai_patch               (begp:endp))           ; this%tsai_patch               (:)   = nan
    allocate(this%elai_patch               (begp:endp))           ; this%elai_patch               (:)   = nan
    allocate(this%elai240_patch            (begp:endp))           ; this%elai240_patch            (:)   = nan
    allocate(this%esai_patch               (begp:endp))           ; this%esai_patch               (:)   = nan
    allocate(this%laisun_patch             (begp:endp))           ; this%laisun_patch             (:)   = nan
    allocate(this%laisha_patch             (begp:endp))           ; this%laisha_patch             (:)   = nan
    allocate(this%laisun_z_patch           (begp:endp,1:nlevcan)) ; this%laisun_z_patch           (:,:) = nan
    allocate(this%laisha_z_patch           (begp:endp,1:nlevcan)) ; this%laisha_z_patch           (:,:) = nan
    allocate(this%mlaidiff_patch           (begp:endp))           ; this%mlaidiff_patch           (:)   = nan
    allocate(this%annlai_patch          (12,begp:endp))           ; this%annlai_patch             (:,:) = nan
    allocate(this%stem_biomass_patch       (begp:endp))           ; this%stem_biomass_patch       (:)   = nan
    allocate(this%leaf_biomass_patch       (begp:endp))           ; this%leaf_biomass_patch       (:)   = nan
    allocate(this%htop_patch               (begp:endp))           ; this%htop_patch               (:)   = nan
    allocate(this%hbot_patch               (begp:endp))           ; this%hbot_patch               (:)   = nan
    allocate(this%z0m_patch                (begp:endp))           ; this%z0m_patch                (:)   = nan
    allocate(this%displa_patch             (begp:endp))           ; this%displa_patch             (:)   = nan
    allocate(this%fsun_patch               (begp:endp))           ; this%fsun_patch               (:)   = nan
    allocate(this%fsun24_patch             (begp:endp))           ; this%fsun24_patch             (:)   = nan
    allocate(this%fsun240_patch            (begp:endp))           ; this%fsun240_patch            (:)   = nan

    allocate(this%dleaf_patch              (begp:endp))           ; this%dleaf_patch              (:)   = nan
    allocate(this%rscanopy_patch           (begp:endp))           ; this%rscanopy_patch           (:)   = nan
!    allocate(this%gccanopy_patch           (begp:endp))           ; this%gccanopy_patch           (:)   = 0.0_r8
    allocate(this%vegwp_patch              (begp:endp,1:nvegwcs)) ; this%vegwp_patch              (:,:) = nan
    allocate(this%vegwp_ln_patch           (begp:endp,1:nvegwcs)) ; this%vegwp_ln_patch           (:,:) = nan
    allocate(this%vegwp_pd_patch           (begp:endp,1:nvegwcs)) ; this%vegwp_pd_patch           (:,:) = nan
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod   , only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(canopystate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: begp, endp
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    this%elai_patch(begp:endp) = spval
    call hist_addfld1d (fname='ELAI', units='m^2/m^2', &
        avgflag='A', long_name='exposed one-sided leaf area index', &
         ptr_patch=this%elai_patch)

    this%esai_patch(begp:endp) = spval
    call hist_addfld1d (fname='ESAI', units='m^2/m^2', &
         avgflag='A', long_name='exposed one-sided stem area index', &
         ptr_patch=this%esai_patch)

    this%laisun_patch(begp:endp) = spval
    call hist_addfld1d (fname='LAISUN', units='m^2/m^2', &
         avgflag='A', long_name='sunlit projected leaf area index', &
         ptr_patch=this%laisun_patch, set_urb=0._r8)

    this%laisha_patch(begp:endp) = spval
    call hist_addfld1d (fname='LAISHA', units='m^2/m^2', &
         avgflag='A', long_name='shaded projected leaf area index', &
         ptr_patch=this%laisha_patch, set_urb=0._r8)

    this%stem_biomass_patch(begp:endp) = spval
    call hist_addfld1d (fname='AGSB', units='kg/m^2', &
         avgflag='A', long_name='Aboveground stem biomass', &
         ptr_patch=this%stem_biomass_patch, default='inactive')

    this%leaf_biomass_patch(begp:endp) = spval
    call hist_addfld1d (fname='AGLB', units='kg/m^2', &
         avgflag='A', long_name='Aboveground leaf biomass', &
         ptr_patch=this%leaf_biomass_patch, default='inactive')

    if (use_cn .or. use_fates) then
       this%fsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='FSUN', units='proportion', &
            avgflag='A', long_name='sunlit fraction of canopy', &
            ptr_patch=this%fsun_patch, default='inactive')

       this%hbot_patch(begp:endp) = spval
       call hist_addfld1d (fname='HBOT', units='m', &
            avgflag='A', long_name='canopy bottom', &
            ptr_patch=this%hbot_patch, default='inactive')

       this%displa_patch(begp:endp) = spval
       call hist_addfld1d (fname='DISPLA', units='m', &
            avgflag='A', long_name='displacement height (vegetated landunits only)', &
            ptr_patch=this%displa_patch, default='inactive', l2g_scale_type='veg')

       if(use_fates_sp)then
          this%htop_input_patch(begp:endp) = spval
          call hist_addfld1d (fname='HTOP', units='m', &
              avgflag='A', long_name='HTOP weights for SP mode', &
              ptr_patch=this%htop_input_patch)
       else
          this%htop_patch(begp:endp) = spval
          call hist_addfld1d (fname='HTOP', units='m', &
              avgflag='A', long_name='canopy top', &
              ptr_patch=this%htop_patch)
       endif
    endif

    if(use_fates_sp)then
      this%tlai_input_patch(begp:endp) = spval
      call hist_addfld1d (fname='TLAI', units='m', &
          avgflag='A', long_name='TLAI weights for SP mode', &
          ptr_patch=this%tlai_input_patch)

      this%tsai_input_patch(begp:endp) = spval
      call hist_addfld1d (fname='TSAI', units='m', &
          avgflag='A', long_name='TSAI weights for SP mode', &
          ptr_patch=this%tsai_input_patch)

    else
       this%tlai_patch(begp:endp) = spval
       call hist_addfld1d (fname='TLAI', units='m^2/m^2', &
           avgflag='A', long_name='total projected leaf area index', &
           ptr_patch=this%tlai_patch)

       this%tsai_patch(begp:endp) = spval
       call hist_addfld1d (fname='TSAI', units='m^2/m^2', &
           avgflag='A', long_name='total projected stem area index', &
           ptr_patch=this%tsai_patch)

    endif !FATES_SP

    this%z0m_patch(begp:endp) = spval
    call hist_addfld1d (fname='Z0MV_DENSE', units='m', &
         avgflag='A', long_name='roughness length over vegetation, momentum, for dense canopy', &
          ptr_patch=this%z0m_patch, default='inactive', l2g_scale_type='veg')

    ! Accumulated fields
    this%fsun24_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSUN24', units='K',  &
         avgflag='A', long_name='fraction sunlit (last 24hrs)', &
         ptr_patch=this%fsun24_patch, default='inactive')

    this%fsun240_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSUN240', units='K',  &
         avgflag='A', long_name='fraction sunlit (last 240hrs)', &
         ptr_patch=this%fsun240_patch, default='inactive')

    this%elai240_patch(begp:endp) = spval
    call hist_addfld1d (fname='LAI240', units='m^2/m^2', &
         avgflag='A', long_name='240hr average of leaf area index', &
         ptr_patch=this%elai240_patch, default='inactive')

    ! Ed specific field
    if ( use_fates ) then
       this%rscanopy_patch(begp:endp) = spval
       call hist_addfld1d (fname='RSCANOPY', units=' s m-1',  &
            avgflag='A', long_name='canopy resistance', &
            ptr_patch=this%rscanopy_patch, set_lake=0._r8, set_urb=0._r8)
    end if

!    call hist_addfld1d (fname='GCCANOPY', units='none',  &
!         avgflag='A', long_name='Canopy Conductance: mmol m-2 s-1', &
!         ptr_patch=this%GCcanopy_patch, set_lake=0._r8, set_urb=0._r8)

    if ( use_hydrstress ) then
       this%vegwp_patch(begp:endp,:) = spval
       call hist_addfld2d (fname='VEGWP',  units='mm', type2d='nvegwcs', &
            avgflag='A', long_name='vegetation water matric potential for sun/sha canopy,xyl,root segments', &
            ptr_patch=this%vegwp_patch)
       this%vegwp_ln_patch(begp:endp,:) = spval
       call hist_addfld2d (fname='VEGWPLN',  units='mm', type2d='nvegwcs', &
            avgflag='A', long_name='vegetation water matric potential for sun/sha canopy,xyl,root at local noon', &
            ptr_patch=this%vegwp_ln_patch, default='active')
       this%vegwp_pd_patch(begp:endp,:) = spval
       call hist_addfld2d (fname='VEGWPPD',  units='mm', type2d='nvegwcs', avgflag='A', &
            long_name='predawn vegetation water matric potential for sun/sha canopy,xyl,root', &
            ptr_patch=this%vegwp_pd_patch, default='active')
    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES
    use accumulMod  , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(canopystate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !---------------------------------------------------------------------

    this%fsun24_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSUN24', units='fraction',                                        &
         desc='24hr average of diffuse solar radiation',  accum_type='runmean', accum_period=-1,   &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%fsun240_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSUN240', units='fraction',                                       &
         desc='240hr average of diffuse solar radiation',  accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%elai240_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='LAI240', units='m2/m2',                                             &
         desc='240hr average of leaf area index',  accum_type='runmean', accum_period=-10,      &
         subgrid_type='pft', numlev=1, init_value=0._r8)

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(canopystate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg="extract_accum_hist allocation error for rbufslp"//&
            errMsg(sourcefile, __LINE__))
    endif

    ! Determine time step
    nstep = get_nstep()

    call extract_accum_field ('FSUN24', rbufslp, nstep)
    this%fsun24_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('FSUN240', rbufslp, nstep)
    this%fsun240_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('LAI240', rbufslp, nstep)
    this%elai240_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('FSUN24', rbufslp, nstep)
    this%fsun24_patch(begp:endp) = rbufslp(begp:endp)

    deallocate(rbufslp)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine ReadNML( this, NLFilename )
    !
    ! Read in canopy parameter namelist
    !
    ! USES:
    use shr_mpi_mod   , only : shr_mpi_bcast
    use abortutils    , only : endrun
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    use clm_varctl    , only : iulog
    use shr_log_mod   , only : errMsg => shr_log_errMsg
    !
    ! ARGUMENTS:
    implicit none
    class(canopystate_type)      :: this
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    ! LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    real(r8) :: leaf_mr_vcm         ! Scalar of leaf respiration to vcmax
    character(len=32) :: subname = 'CanopyStateType::ReadNML'  ! subroutine name
    !-----------------------------------------------------------------------
    namelist / clm_canopy_inparm / leaf_mr_vcm

    ! ----------------------------------------------------------------------
    ! Read namelist from input namelist filename
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in clm_canopy_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_canopy_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_canopy_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading clm_canopy_inparm namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR finding clm_canopy_inparm namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )

    end if

    ! Broadcast namelist variables read in
    call shr_mpi_bcast(leaf_mr_vcm, mpicom)
    this%leaf_mr_vcm = leaf_mr_vcm

  end subroutine ReadNML

 !-----------------------------------------------------------------------

  subroutine SetNMLForTesting( this )
     !
     ! Set canopy parameter namelist control settings for unit-testing
     !
     class(canopystate_type)      :: this
     ! LOCAL VARIABLES:
     !-----------------------------------------------------------------------
 
     
     this%leaf_mr_vcm = 0.015_r8
   
   end subroutine SetNMLForTesting

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! USES
    use clm_time_manager, only : get_nstep
    use accumulMod      , only : update_accum_field, extract_accum_field
    use abortutils      , only : endrun
    !
    ! !ARGUMENTS:
    class(canopystate_type)             :: this
    type(bounds_type)      , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g,p                       ! indices
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: ier                       ! error status
    integer :: begp, endp
    real(r8), pointer :: rbufslp(:)      ! temporary single level - patch level
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level patch field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    ! Accumulate and extract fsun24 & fsun240
    do p = begp,endp
       rbufslp(p) = this%fsun_patch(p)
    end do
    call update_accum_field  ('FSUN24' , rbufslp              , nstep)
    call extract_accum_field ('FSUN24' , this%fsun24_patch    , nstep)
    call update_accum_field  ('FSUN240', rbufslp              , nstep)
    call extract_accum_field ('FSUN240', this%fsun240_patch   , nstep)

    ! Accumulate and extract elai240
    do p = begp,endp
       rbufslp(p) = this%elai_patch(p)
    end do
    call update_accum_field  ('LAI240', rbufslp               , nstep)
    call extract_accum_field ('LAI240', this%elai240_patch    , nstep)

    deallocate(rbufslp)

  end subroutine UpdateAccVars

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !ARGUMENTS:
    class(canopystate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: p,l,c,g
    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       l = patch%landunit(p)

       this%tlai_patch(p)        = 0._r8
       this%tsai_patch(p)        = 0._r8
       this%elai_patch(p)        = 0._r8
       this%esai_patch(p)        = 0._r8
       this%stem_biomass_patch(p)= 0._r8
       this%leaf_biomass_patch(p)= 0._r8
       this%htop_patch(p)        = 0._r8
       this%hbot_patch(p)        = 0._r8
       this%vegwp_patch(p,:)     = -2.5e4_r8

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%laisun_patch(p) = 0._r8
          this%laisha_patch(p) = 0._r8
       end if

       this%tlai_input_patch(p)       = 0._r8
       this%tsai_input_patch(p)       = 0._r8
       this%htop_input_patch(p)       = 0._r8
       this%hbot_input_patch(p)       = 0._r8
       
       ! needs to be initialized to spval to avoid problems when averaging for the accum
       ! field
       this%fsun_patch(p) = spval
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_double, ncd_int
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(canopystate_type) :: this
    type(bounds_type) , intent(in)    :: bounds
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,p,c,iv ! indices
    logical :: readvar      ! determine if variable is on initial file
    integer :: begp, endp
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    call restartvar(ncid=ncid, flag=flag, varname='FRAC_VEG_NOSNO_ALB', xtype=ncd_int,  &
         dim1name='pft', long_name='fraction of vegetation not covered by snow (0 or 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frac_veg_nosno_alb_patch)

    call restartvar(ncid=ncid, flag=flag, varname='tlai', xtype=ncd_double,  &
         dim1name='pft', long_name='one-sided leaf area index, no burying by snow', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tlai_patch)

    call restartvar(ncid=ncid, flag=flag, varname='tsai', xtype=ncd_double,  &
         dim1name='pft', long_name='one-sided stem area index, no burying by snow', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tsai_patch)
         
     call restartvar(ncid=ncid, flag=flag, varname='tlai_input', xtype=ncd_double,  &
         dim1name='pft', long_name='sp mode driver data for one-sided leaf area index, no burying by snow', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tlai_input_patch)

    call restartvar(ncid=ncid, flag=flag, varname='tsai_input', xtype=ncd_double,  &
         dim1name='pft', long_name='sp mode driver data for one-sided stem area index, no burying by snow', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tsai_input_patch)

    call restartvar(ncid=ncid, flag=flag, varname='elai', xtype=ncd_double,  &
         dim1name='pft', long_name='one-sided leaf area index, with burying by snow', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%elai_patch)

    call restartvar(ncid=ncid, flag=flag, varname='esai', xtype=ncd_double,  &
         dim1name='pft', long_name='one-sided stem area index, with burying by snow', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%esai_patch)

    call restartvar(ncid=ncid, flag=flag, varname='stem_biomass', xtype=ncd_double,  &
         dim1name='pft', long_name='stem biomass', units='kg/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%stem_biomass_patch)

    call restartvar(ncid=ncid, flag=flag, varname='leaf_biomass', xtype=ncd_double,  &
         dim1name='pft', long_name='leaf biomass', units='kg/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%leaf_biomass_patch)

    call restartvar(ncid=ncid, flag=flag, varname='htop', xtype=ncd_double,  &
         dim1name='pft', long_name='canopy top', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%htop_patch)
         
     call restartvar(ncid=ncid, flag=flag, varname='htop_input', xtype=ncd_double,  &
         dim1name='pft', long_name='sp-mode driver data for canopy top', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%htop_input_patch)

    call restartvar(ncid=ncid, flag=flag, varname='hbot', xtype=ncd_double,  &
         dim1name='pft', long_name='canopy bottom', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%hbot_patch)
         
     call restartvar(ncid=ncid, flag=flag, varname='hbot_input', xtype=ncd_double,  &
         dim1name='pft', long_name='sp-mode driver data for canopy bottom', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%hbot_input_patch)

    call restartvar(ncid=ncid, flag=flag, varname='mlaidiff', xtype=ncd_double,  &
         dim1name='pft', long_name='difference between lai month one and month two', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%mlaidiff_patch)

    call restartvar(ncid=ncid, flag=flag, varname='fsun', xtype=ncd_double,  &
         dim1name='pft', long_name='sunlit fraction of canopy', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fsun_patch)



    if (flag=='read' )then
       do p = bounds%begp,bounds%endp
          if (shr_infnan_isnan(this%fsun_patch(p)) ) then
             this%fsun_patch(p) = spval
          end if
       end do
    end if

    if ( use_hydrstress ) then
       call restartvar(ncid=ncid, flag=flag, varname='vegwp', xtype=ncd_double,  &
            dim1name='pft', dim2name='vegwcs', switchdim=.true., &
            long_name='vegetation water matric potential', units='mm', &
            scale_by_thickness=.false., &
            interpinic_flag='interp', readvar=readvar, data=this%vegwp_patch)

       call restartvar(ncid=ncid, flag=flag, varname='VEGWPLN', xtype=ncd_double,  &
            dim1name='pft', dim2name='vegwcs', switchdim=.false., &
            long_name='vegetation water matric potential for sun/sha canopy,xyl,root at local noon', units='mm', &
            scale_by_thickness=.false., &
            interpinic_flag='skip', readvar=readvar, data=this%vegwp_ln_patch)

       call restartvar(ncid=ncid, flag=flag, varname='VEGWPPD', xtype=ncd_double,  &
            dim1name='pft', dim2name='vegwcs', switchdim=.false., &
            long_name='predawn vegetation water matric potential for sun/sha canopy,xyl,root', units='mm', &
            scale_by_thickness=.false., &
            interpinic_flag='skip', readvar=readvar, data=this%vegwp_pd_patch)

    end if

  end subroutine Restart

end module CanopyStateType
