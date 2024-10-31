module glc2lndMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle arrays used for exchanging data from glc to clm.
  !
  ! !USES:
#include "shr_assert.h"
  use decompMod      , only : bounds_type, subgrid_level_gridcell
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : maxpatch_glc
  use clm_varctl     , only : iulog, glc_do_dynglacier
  use clm_varcon     , only : spval, ispval
  use abortutils     , only : endrun
  use GridcellType   , only : grc
  use LandunitType   , only : lun
  use ColumnType     , only : col
  use landunit_varcon, only : istice
  use glcBehaviorMod , only : glc_behavior_type
  !
  ! !REVISION HISTORY:
  ! Created by William Lipscomb, Dec. 2007, based on clm_atmlnd.F90.
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! glc -> land variables structure
  type, public :: glc2lnd_type

     ! ------------------------------------------------------------------------
     ! Public data
     ! ------------------------------------------------------------------------

     ! Where we should do SMB-related runoff routing that is appropriate for having a dynamic icesheet underneath.
     real(r8), pointer :: glc_dyn_runoff_routing_grc (:) => null()

     ! ------------------------------------------------------------------------
     ! Private data
     ! ------------------------------------------------------------------------

     type(glc_behavior_type), pointer, private :: glc_behavior  ! reference to the glc_behavior instance

     real(r8), pointer, private :: frac_grc    (:,:) => null()
     real(r8), pointer, private :: topo_grc    (:,:) => null()
     real(r8), pointer, private :: hflx_grc    (:,:) => null()

     ! Area in which values from GLC are valid, received from GLC (0-1)
     real(r8), pointer, private :: icemask_grc (:)   => null()

     ! icemask_coupled_fluxes_grc is like icemask_grc, but the mask only contains icesheet
     ! points that potentially send non-zero fluxes to the coupler. i.e., it does not
     ! contain icesheets that are diagnostic only, because for those diagnostic ice sheets
     ! (which do not send calving fluxes to the coupler), we need to use the non-dynamic
     ! form of runoff routing in CLM in order to conserve water properly.
     !
     ! (However, note that this measure of "diagnostic-only" does not necessarily
     ! correspond to whether CLM is updating its glacier areas there - for example, we
     ! could theoretically have an icesheet whose areas are evolving, and CLM is updating
     ! its glacier areas to match, but where we're zeroing out the fluxes sent to the
     ! coupler, and so we're using the non-dynamic form of runoff routing in CLM.)
     real(r8), pointer, private :: icemask_coupled_fluxes_grc (:)  => null()

   contains

     ! ------------------------------------------------------------------------
     ! Public routines
     ! ------------------------------------------------------------------------

     procedure, public  :: Init
     procedure, public  :: Clean

     ! In each timestep, these routines should be called in order (though they don't need
     ! to be called all at once):
     ! - set_glc2lnd_fields
     ! - update_glc2lnd_fracs
     ! - update_glc2lnd_topo
     procedure, public  :: set_glc2lnd_fields_nuopc ! set coupling fields sent from glc to lnd
     procedure, public  :: update_glc2lnd_fracs     ! update subgrid fractions based on input from GLC
     procedure, public  :: update_glc2lnd_topo      ! update topographic heights

     ! For unit testing only:
     procedure, public  :: for_test_set_glc2lnd_fields_directly  ! set glc2lnd fields directly in a unit testing context

     ! ------------------------------------------------------------------------
     ! Private routines
     ! ------------------------------------------------------------------------

     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

     ! call various wrapup routines at the end of setting glc2lnd fields
     procedure, private :: set_glc2lnd_fields_wrapup

     ! sanity-check icemask from GLC
     procedure, private :: check_glc2lnd_icemask

     ! sanity-check icemask_coupled_fluxes from GLC
     procedure, private :: check_glc2lnd_icemask_coupled_fluxes

     ! update glc_dyn_runoff_routing field based on input from GLC
     procedure, private :: update_glc2lnd_dyn_runoff_routing

  end type glc2lnd_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, glc_behavior)

    class(glc2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds
    type(glc_behavior_type), intent(in), target :: glc_behavior

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds, glc_behavior)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize glc variables required by the land
    !
    ! !ARGUMENTS:
    class (glc2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begg,endg
    !------------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    allocate(this%frac_grc    (begg:endg,0:maxpatch_glc)) ;   this%frac_grc    (:,:) = nan
    allocate(this%topo_grc    (begg:endg,0:maxpatch_glc)) ;   this%topo_grc    (:,:) = nan
    allocate(this%hflx_grc    (begg:endg,0:maxpatch_glc)) ;   this%hflx_grc    (:,:) = nan
    allocate(this%icemask_grc (begg:endg))                   ;   this%icemask_grc (:)   = nan
    allocate(this%icemask_coupled_fluxes_grc (begg:endg))    ;   this%icemask_coupled_fluxes_grc (:)   = nan
    allocate(this%glc_dyn_runoff_routing_grc (begg:endg))    ;   this%glc_dyn_runoff_routing_grc (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(glc2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg

    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    begg = bounds%begg
    endg = bounds%endg

    this%icemask_grc(begg:endg) = spval
    call hist_addfld1d (fname='ICE_MODEL_FRACTION',  units='unitless',  &
         avgflag='A', long_name='Ice sheet model fractional coverage', &
         ptr_gcell=this%icemask_grc, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, glc_behavior)
    !
    ! !USES:
    use domainMod      , only : ldomain
    !
    ! !ARGUMENTS:
    class(glc2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds
    type(glc_behavior_type), intent(in), target :: glc_behavior
    !
    ! !LOCAL VARIABLES:
    integer :: begg, endg

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    begg = bounds%begg
    endg = bounds%endg

    this%glc_behavior => glc_behavior

    this%frac_grc(begg:endg, :) = 0.0_r8
    this%topo_grc(begg:endg, :) = 0.0_r8
    this%hflx_grc(begg:endg, :) = 0.0_r8

    ! When running with a stub glc model, it's important that icemask_grc be initialized
    ! to 0 everywhere. With an active glc model, icemask_grc will be updated in the first
    ! time step, and it isn't needed before then, so it's safe to initialize it to 0.
    ! Since icemask is 0, icemask_coupled_fluxes needs to be 0, too (and the latter is
    ! safest in case we aren't coupled to CISM, to ensure that we use the uncoupled form
    ! of runoff routing).
    this%icemask_grc(begg:endg) = 0.0_r8
    this%icemask_coupled_fluxes_grc(begg:endg) = 0.0_r8

    call this%update_glc2lnd_dyn_runoff_routing(bounds)

  end subroutine InitCold


  !-----------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !DESCRIPTION:
    ! Deallocate memory in this object
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------

    deallocate(this%frac_grc)
    deallocate(this%topo_grc)
    deallocate(this%hflx_grc)
    deallocate(this%icemask_grc)
    deallocate(this%icemask_coupled_fluxes_grc)
    deallocate(this%glc_dyn_runoff_routing_grc)

  end subroutine Clean

  !-----------------------------------------------------------------------
  subroutine set_glc2lnd_fields_nuopc(this, bounds, glc_present, &
       frac_grc, topo_grc, hflx_grc, icemask_grc, icemask_coupled_fluxes_grc)
    !
    ! !DESCRIPTION:
    ! Set coupling fields sent from glc to lnd
    !
    ! If glc_present is true, then the given fields are all assumed to be valid; if
    ! glc_present is false, then these fields are ignored.
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(inout) :: this
    type(bounds_type)  , intent(in)    :: bounds
    logical  , intent(in) :: glc_present                              ! true if running with a non-stub glc model
    real(r8) , intent(in) :: frac_grc(bounds%begg:,0:)                ! ice-covered field for each elevation class
    real(r8) , intent(in) :: topo_grc(bounds%begg:,0:)                ! topo field for each elevation class
    real(r8) , intent(in) :: hflx_grc(bounds%begg:,0:)                ! heat flux field for each elevation class
    real(r8) , intent(in) :: icemask_grc(bounds%begg:)                ! icemask field
    real(r8) , intent(in) :: icemask_coupled_fluxes_grc(bounds%begg:) ! icemask_coupled_fluxes field
    !
    ! !LOCAL VARIABLES:
    integer :: g
    integer :: ice_class

    character(len=*), parameter :: subname = 'set_glc2lnd_fields_nuopc'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(frac_grc, 1) == bounds%endg), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(topo_grc, 1) == bounds%endg), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(hflx_grc, 1) == bounds%endg), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(icemask_grc,1) == bounds%endg), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(icemask_coupled_fluxes_grc,1) == bounds%endg), sourcefile, __LINE__)

    if (glc_present) then
       do g = bounds%begg, bounds%endg
          do ice_class = 0, maxpatch_glc
             this%frac_grc(g,ice_class)  = frac_grc(g,ice_class)
             this%topo_grc(g,ice_class)  = topo_grc(g,ice_class)
             this%hflx_grc(g,ice_class)  = hflx_grc(g,ice_class)
          end do
          this%icemask_grc(g)  = icemask_grc(g)
          this%icemask_coupled_fluxes_grc(g)  = icemask_coupled_fluxes_grc(g)
       end do

       call this%set_glc2lnd_fields_wrapup(bounds)
    else
       if (glc_do_dynglacier) then
          call endrun(' ERROR: With glc_present false (e.g., a stub glc model), glc_do_dynglacier must be false '// &
               errMsg(sourcefile, __LINE__))
       end if
    end if

  end subroutine set_glc2lnd_fields_nuopc

  !-----------------------------------------------------------------------
  subroutine for_test_set_glc2lnd_fields_directly(this, bounds, &
       topo, icemask)
    !
    ! !DESCRIPTION:
    ! Set glc2lnd fields directly in a unit testing context
    !
    ! This currently only provides a mechanism to set fields that are actually needed in
    ! our unit tests. More could be added later.
    !
    ! Also: In contrast to the production version (set_glc2lnd_fields), this does NOT
    ! currently update glc2lnd_dyn_runoff_routing (because doing so would require having a
    ! sensible glc_behavior, which we may not have; and also, we currently don't need this
    ! field in a unit testing context). (Note: If we eventually want/need to update
    ! glc2lnd_dyn_runoff_routing, and thus need a fully sensible glc_behavior, then we
    ! should extract the self-calls at the end of set_glc2lnd_fields
    ! (check_glc2lnd_icemask, check_glc2lnd_icemask_coupled_fluxes,
    ! update_glc2lnd_dyn_runoff_routing) into a private routine like
    ! set_glc2lnd_fields_wrapup, which could be called by both set_glc2lnd_fields and this
    ! routine.)
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(inout) :: this
    type(bounds_type)  , intent(in)    :: bounds
    real(r8), intent(in), optional :: topo( bounds%begg: , 0: )  ! topographic height [gridcell, elevclass]
    real(r8), intent(in), optional :: icemask( bounds%begg: )
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'for_test_set_glc2lnd_fields_directly'
    !-----------------------------------------------------------------------

    if (present(topo)) then
       SHR_ASSERT_ALL_FL((ubound(topo) == (/bounds%endg, maxpatch_glc/)), sourcefile, __LINE__)
       this%topo_grc(bounds%begg:bounds%endg, 0:maxpatch_glc) = topo(bounds%begg:bounds%endg, 0:maxpatch_glc)
    end if

    if (present(icemask)) then
       SHR_ASSERT_ALL_FL((ubound(icemask) == (/bounds%endg/)), sourcefile, __LINE__)
       this%icemask_grc(bounds%begg:bounds%endg) = icemask(bounds%begg:bounds%endg)
    end if

  end subroutine for_test_set_glc2lnd_fields_directly

  !-----------------------------------------------------------------------
  subroutine set_glc2lnd_fields_wrapup(this, bounds)
    !
    ! !DESCRIPTION:
    ! Call various wrapup routines at the end of setting glc2lnd fields
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(inout) :: this
    type(bounds_type)  , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'set_glc2lnd_fields_wrapup'
    !-----------------------------------------------------------------------

    call this%check_glc2lnd_icemask(bounds)
    call this%check_glc2lnd_icemask_coupled_fluxes(bounds)
    call this%update_glc2lnd_dyn_runoff_routing(bounds)

  end subroutine set_glc2lnd_fields_wrapup

  !-----------------------------------------------------------------------
  subroutine check_glc2lnd_icemask(this, bounds)
    !
    ! !DESCRIPTION:
    ! Do a sanity check on the icemask received from CISM via coupler.
    !
    ! !USES:
    use domainMod , only : ldomain
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(in) :: this
    type(bounds_type)       , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! grid cell index

    character(len=*), parameter :: subname = 'check_glc2lnd_icemask'
    !-----------------------------------------------------------------------

    do g = bounds%begg, bounds%endg

       if (this%icemask_grc(g) > 0._r8) then

          ! Ensure that, within the icemask, there are no points that have (non-virtual
          ! and compute-SMB). This is important for two reasons:
          !
          ! (1) To ensure that, in grid cells where we're producing SMB, we have SMB for
          !     all elevation classes, so that the downscaling / vertical interpolation
          !     can be done correctly.
          !
          ! (2) To avoid conservation issues, we want to ensure that, in grid cells where
          !     we're producing SMB and are dynamically coupled to the ice sheet (if 2-way
          !     coupling is enabled), glacier areas are remaining in-sync with glc. (Note
          !     that has_virtual_columns_grc dictates where we're able to keep glacier
          !     areas in sync with glc.) (In principle, I think this one could check
          !     icemask_coupled_fluxes rather than icemask; we check icemask because we
          !     needed to check icemask for the other reason anyway; this is okay because
          !     icemask_coupled_fluxes is a subset of icemask.)
          if (this%glc_behavior%melt_replaced_by_ice_grc(g) .and. &
               .not. this%glc_behavior%has_virtual_columns_grc(g)) then
             write(iulog,'(a)') subname//' ERROR: Within the icemask, there cannot be any points that have'
             write(iulog,'(a)') '(non-virtual and compute-SMB).'
             write(iulog,'(a)') 'Ensure that GLACIER_REGION on the surface dataset and the namelist items,'
             write(iulog,'(a)') 'glacier_region_behavior and glacier_region_melt_behavior are all set correctly:'
             write(iulog,'(a)') 'Typically, the region encompassing the active GLC domain should specify'
             write(iulog,'(a)') 'glacier_region_behavior="virtual" and glacier_region_melt_behavior="replaced_by_ice".'
             write(iulog,'(a)') '(But it is also okay for part of the GLC domain to have'
             write(iulog,'(a)') 'glacier_region_melt_behavior="remains_in_place"; this part of the domain can have'
             write(iulog,'(a)') 'any setting for glacier_region_behavior.)'
             call endrun(subgrid_index=g, subgrid_level=subgrid_level_gridcell, msg=errMsg(sourcefile, __LINE__))
          end if

       end if
    end do

  end subroutine check_glc2lnd_icemask

  !-----------------------------------------------------------------------
  subroutine check_glc2lnd_icemask_coupled_fluxes(this, bounds)
    !
    ! !DESCRIPTION:
    ! Do a sanity check on the icemask_coupled_fluxes field received from CISM via coupler.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(in) :: this
    type(bounds_type)  , intent(in) :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! grid cell index

    character(len=*), parameter :: subname = 'check_glc2lnd_icemask_coupled_fluxes'
    !-----------------------------------------------------------------------

    do g = bounds%begg, bounds%endg

       ! Ensure that icemask_coupled_fluxes is a subset of icemask. This is helpful to
       ! ensure that the consistency checks that are done on glc behaviors within the
       ! icemask (in check_glc2lnd_icemask) also apply within the icemask_coupled_fluxes
       ! region. Other than that convenience, there currently is no code in CLM that
       ! depends on this relationship, but it seems helpful to ensure that this intuitive
       ! relationship holds, so that code developed in the future can rely on it.
       if (this%icemask_coupled_fluxes_grc(g) > 0._r8 .and. this%icemask_grc(g) == 0._r8) then
          write(iulog,*) subname//' ERROR: icemask_coupled_fluxes must be a subset of icemask.'
          call endrun(subgrid_index=g, subgrid_level=subgrid_level_gridcell, msg=errMsg(sourcefile, __LINE__))
       end if

    end do

  end subroutine check_glc2lnd_icemask_coupled_fluxes

  !-----------------------------------------------------------------------
  subroutine update_glc2lnd_dyn_runoff_routing(this, bounds)
    !
    ! !DESCRIPTION:
    ! Update glc_dyn_runoff_routing field based on updated icemask_coupled_fluxes field
    !
    ! !USES:
    use domainMod , only : ldomain
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(inout) :: this
    type(bounds_type)  , intent(in) :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! grid cell index

    character(len=*), parameter :: subname = 'update_glc2lnd_dyn_runoff_routing'
    !-----------------------------------------------------------------------

    ! Wherever we have an icesheet that is computing and sending fluxes to the coupler -
    ! which particularly means it is computing a calving flux - we will use the
    ! "glc_dyn_runoff_routing" scheme, with 0 < glc_dyn_runoff_routing <= 1.
    ! In these places, all or part of the snowcap flux goes to CISM rather than the runoff model.
    ! In other places - including places where CISM is not running at all, as well as places
    ! where CISM is running in diagnostic-only mode and therefore is not sending a calving flux -
    ! we have glc_dyn_runoff_routing = 0, and the snowcap flux goes to the runoff model.
    ! This is needed to conserve water correctly in the absence of a calving flux.
    !
    ! In places where we are not computing SMB, we also have glc_dyn_runoff_routing = 0.
    ! Currently glc_dyn_runoff_routing is only used where we're computing SMB, but if it
    ! were ever used elsewhere, it seems best to have it set to 0 there: this seems
    ! consistent with the fact that we zero out the SMB flux sent to GLC in that region.
    ! (However, it's possible that, once we start actually using glc_dyn_runoff_routing
    ! for some purpose outside the do_smb filter, we'll discover that this logic should
    ! be changed.)

    do g = bounds%begg, bounds%endg

       if (this%glc_behavior%melt_replaced_by_ice_grc(g)) then
          ! As noted in the comments at the top of this routine, we only set
          ! glc_dyn_runoff_routing where we are computing SMB

          ! Set glc_dyn_runoff_routing_grc(g) to a value in the range [0,1].
          !
          ! This value gives the grid cell fraction that is deemed to be coupled to the
          ! dynamic ice sheet model. For this fraction of the grid cell, snowcap fluxes are
          ! sent to the ice sheet model. The remainder of the grid cell sends snowcap fluxes
          ! to the runoff model.
          !
          ! Note: The coupler (in prep_glc_mod.F90) assumes that the fraction coupled to the
          !       dynamic ice sheet model is min(lfrac, Sg_icemask_l), where lfrac is the
          !       "frac" component of fraction_lx, and Sg_icemask_l is obtained by mapping
          !       Sg_icemask_g from the glc to the land grid. Here, ldomain%frac is
          !       equivalent to lfrac, and this%icemask_grc is equivalent to Sg_icemask_l.
          !       However, here we use icemask_coupled_fluxes_grc, so that we route all snow
          !       capping to runoff in areas where the ice sheet is not generating calving
          !       fluxes. In addition, here we need to divide by lfrac, because the coupler
          !       multiplies by it later (and, for example, if lfrac = 0.1 and
          !       icemask_coupled_fluxes = 1, we want all snow capping to go to the ice
          !       sheet model, not to the runoff model).
          !
          ! Note: In regions where CLM overlaps the CISM domain, this%icemask_grc(g) typically
          !       is nearly equal to ldomain%frac(g). So an alternative would be to simply set
          !       glc_dyn_runoff_routing_grc(g) = icemask_grc(g).
          !       The reason to cap glc_dyn_runoff_routing at lfrac is to avoid sending the
          !       ice sheet model a greater mass of water (in the form of snowcap fluxes)
          !       than is allowed to fall on a CLM grid cell that is part ocean.

          ! TODO(wjs, 2017-05-08) Ideally, we wouldn't have this duplication in logic
          ! between the coupler and CLM. The best solution would be to have the coupler
          ! itself do the partitioning of the snow capping flux between the ice sheet model
          ! and the runoff model. A next-best solution would be to have the coupler send a
          ! field to CLM telling it what fraction of snow capping should go to the runoff
          ! model in each grid cell.

          if (ldomain%frac(g) == 0._r8) then
             ! Avoid divide by 0; note that, in this case, the amount going to runoff isn't
             ! important for system-wide conservation, so we could really choose anything we
             ! want.
             this%glc_dyn_runoff_routing_grc(g) = this%icemask_coupled_fluxes_grc(g)
          else
             this%glc_dyn_runoff_routing_grc(g) = &
                  min(ldomain%frac(g), this%icemask_coupled_fluxes_grc(g)) / &
                  ldomain%frac(g)
          end if

       else  ! .not. this%glc_behavior%melt_replaced_by_ice_grc(g)
          ! As noted in the comments at the top of this routine, we set
          ! glc_dyn_runoff_routing to 0 where we are not computing SMB. (This assumes that
          ! gridcells where we compute SMB are the same as gridcells for which
          ! melt_replaced_by_ice is true.)
          this%glc_dyn_runoff_routing_grc(g) = 0._r8
       end if
          
    end do

  end subroutine update_glc2lnd_dyn_runoff_routing



  !-----------------------------------------------------------------------
  subroutine update_glc2lnd_fracs(this, bounds)
    !
    ! !DESCRIPTION:
    ! Update subgrid fractions based on input from GLC (via the coupler)
    !
    ! The weights updated here are some col%wtlunit and lun%wtgcell values
    !
    ! If glc_do_dynglacier is false, nothing is changed
    !
    ! !USES:
    use column_varcon     , only : col_itype_to_ice_class
    use subgridWeightsMod , only : set_landunit_weight
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(in) :: this
    type(bounds_type)  , intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g,c                              ! indices
    real(r8):: area_ice                         ! area of the ice landunit
    integer :: l_ice                            ! index of the ice landunit
    integer :: ice_class                        ! current ice class (1..maxpatch_glc)
    logical :: frac_assigned(1:maxpatch_glc) ! whether this%frac has been assigned for each elevation class
    logical :: error                            ! if an error was found

    character(len=*), parameter :: subname = 'update_glc2lnd_fracs'
    !-----------------------------------------------------------------------

    if (glc_do_dynglacier) then
       do g = bounds%begg, bounds%endg
          ! Values from GLC are only valid within the icemask, so we only update CLM's
          ! areas there. Also, we only update areas where the glacier region behavior is
          ! 'virtual', because that's the only region where we are guaranteed to have all
          ! of the elevation classes we need in order to remain in sync. (Note that, for
          ! conservation purposes, it's important that we update areas in all regions
          ! where we're fully-two-way-coupled to the icesheet and we're computing SMB;
          ! this requirement is checked in check_glc2lnd_icemask.) (This conditional
          ! should be kept consistent with the conditional in update_glc2lnd_topo.)
          if (this%icemask_grc(g) > 0._r8 .and. this%glc_behavior%has_virtual_columns_grc(g)) then

             ! Set total ice landunit area
             area_ice = sum(this%frac_grc(g, 1:maxpatch_glc))
             call set_landunit_weight(g, istice, area_ice)

             ! If new landunit area is greater than 0, then update column areas
             ! (If new landunit area is 0, col%wtlunit is arbitrary, so we might as well keep the existing values)
             if (area_ice > 0) then
                ! Determine index of the ice landunit
                l_ice = grc%landunit_indices(istice, g)
                if (l_ice == ispval) then
                   write(iulog,*) subname//' ERROR: no ice landunit found within the icemask, for g = ', g
                   call endrun(subgrid_index=g, subgrid_level=subgrid_level_gridcell, &
                        msg="no ice landunit found within the icemask")
                end if

                frac_assigned(1:maxpatch_glc) = .false.
                do c = lun%coli(l_ice), lun%colf(l_ice)
                   ice_class = col_itype_to_ice_class(col%itype(c))
                   col%wtlunit(c) = this%frac_grc(g, ice_class) / lun%wtgcell(l_ice)
                   frac_assigned(ice_class) = .true.
                end do

                ! Confirm that all elevation classes that have non-zero area according to
                ! this%frac have been assigned to a column in CLM's data structures
                error = .false.
                do ice_class = 1, maxpatch_glc
                   if (this%frac_grc(g, ice_class) > 0._r8 .and. &
                        .not. frac_assigned(ice_class)) then
                      error = .true.
                   end if
                end do
                if (error) then
                   write(iulog,*) subname//' ERROR: at least one glc column has non-zero area from the coupler,'
                   write(iulog,*) 'but there was no slot in memory for this column; g = ', g
                   write(iulog,*) 'this%frac_grc(g, 1:maxpatch_glc) = ', &
                        this%frac_grc(g, 1:maxpatch_glc)
                   write(iulog,*) 'frac_assigned(1:maxpatch_glc) = ', &
                        frac_assigned(1:maxpatch_glc)
                   call endrun(subgrid_index=g, subgrid_level=subgrid_level_gridcell, &
                        msg="at least one glc column has non-zero area from cpl but has no slot in memory")
                end if  ! error
             end if  ! area_ice > 0
          end if  ! this%icemask_grc(g) > 0 .and. this%glc_behavior%has_virtual_columns_grc(g)
       end do  ! g
    end if  ! glc_do_dynglacier

  end subroutine update_glc2lnd_fracs

  !-----------------------------------------------------------------------
  subroutine update_glc2lnd_topo(this, bounds, topo_col, needs_downscaling_col)
    !
    ! !DESCRIPTION:
    ! Update column-level topographic heights based on input from GLC (via the coupler).
    !
    ! Also updates the logical array, needs_downscaling_col: Sets this array to true
    ! anywhere where topo_col is updated, because these points will need downscaling.
    ! (Leaves other array elements in needs_downscaling_col untouched.)
    !
    ! If glc_do_dynglacier is false, then both topographic heights and
    ! needs_downscaling_col are left unchanged.
    !
    ! !USES:
    use landunit_varcon , only : istice
    use column_varcon   , only : col_itype_to_ice_class
    !
    ! !ARGUMENTS:
    class(glc2lnd_type) , intent(in)    :: this
    type(bounds_type)   , intent(in)    :: bounds                   ! bounds
    real(r8)            , intent(inout) :: topo_col( bounds%begc: ) ! topographic height (m)
    logical             , intent(inout) :: needs_downscaling_col( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    integer :: c, l, g      ! indices
    integer :: ice_class    ! current ice class (1..maxpatch_glc)

    character(len=*), parameter :: subname = 'update_glc2lnd_topo'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(topo_col) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(needs_downscaling_col) == (/bounds%endc/)), sourcefile, __LINE__)

    if (glc_do_dynglacier) then
       do c = bounds%begc, bounds%endc
          l = col%landunit(c)
          g = col%gridcell(c)

          ! Values from GLC are only valid within the icemask, so we only update CLM's
          ! topo values there. Also, consistently with the conditional in
          ! update_glc2lnd_fracs, we only update topo values where the glacier region
          ! behavior is 'virtual': it could be problematic to update topo values in a
          ! grid cell where we're not updating areas.
          if (this%icemask_grc(g) > 0._r8 .and. this%glc_behavior%has_virtual_columns_grc(g)) then
             if (lun%itype(l) == istice) then
                ice_class = col_itype_to_ice_class(col%itype(c))
             else
                ! If not on a glaciated column, assign topography to the bare-land value determined by GLC.
                ice_class = 0
             end if

             ! Note that we do downscaling over all column types. This is for consistency:
             ! interpretation of results would be difficult if some non-glacier column types
             ! were downscaled but others were not.
             !
             ! BUG(wjs, 2016-11-15, bugz 2377) Actually, do not downscale over urban points:
             ! this currently isn't allowed because the urban code references some
             ! non-downscaled, gridcell-level atmospheric forcings
             if (.not. lun%urbpoi(l)) then
                topo_col(c) = this%topo_grc(g, ice_class)
                needs_downscaling_col(c) = .true.
             end if
          end if
       end do
    end if

  end subroutine update_glc2lnd_topo

end module glc2lndMod
