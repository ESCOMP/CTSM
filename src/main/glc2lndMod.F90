module glc2lndMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle arrays used for exchanging data from glc to clm.
  !
  ! !USES:
#include "shr_assert.h"
  use decompMod      , only : bounds_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : maxpatch_glcmec
  use clm_varctl     , only : iulog, glc_do_dynglacier
  use clm_varcon     , only : nameg, spval, ispval
  use abortutils     , only : endrun
  use GridcellType   , only : grc 
  use LandunitType   , only : lun
  use ColumnType     , only : col
  use landunit_varcon, only : istice_mec
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

     ! Where we should do runoff routing that is appropriate for having a dynamic icesheet underneath.
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
     procedure, public  :: set_glc2lnd_fields       ! set coupling fields sent from glc to lnd
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

    allocate(this%frac_grc    (begg:endg,0:maxpatch_glcmec)) ;   this%frac_grc    (:,:) = nan
    allocate(this%topo_grc    (begg:endg,0:maxpatch_glcmec)) ;   this%topo_grc    (:,:) = nan
    allocate(this%hflx_grc    (begg:endg,0:maxpatch_glcmec)) ;   this%hflx_grc    (:,:) = nan
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
  subroutine set_glc2lnd_fields(this, bounds, glc_present, x2l, &
       index_x2l_Sg_ice_covered, index_x2l_Sg_topo, index_x2l_Flgg_hflx, &
       index_x2l_Sg_icemask, index_x2l_Sg_icemask_coupled_fluxes)
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
    logical  , intent(in) :: glc_present                         ! true if running with a non-stub glc model
    real(r8) , intent(in) :: x2l(:, bounds%begg: )               ! driver import state to land model [field, gridcell]
    integer  , intent(in) :: index_x2l_Sg_ice_covered( 0: )      ! indices of ice-covered field in x2l, for each elevation class
    integer  , intent(in) :: index_x2l_Sg_topo( 0: )             ! indices of topo field in x2l, for each elevation class
    integer  , intent(in) :: index_x2l_Flgg_hflx( 0: )           ! indices of heat flux field in x2l, for each elevation class
    integer  , intent(in) :: index_x2l_Sg_icemask                ! index of icemask field in x2l
    integer  , intent(in) :: index_x2l_Sg_icemask_coupled_fluxes ! index of icemask_coupled_fluxes field in x2l
    !
    ! !LOCAL VARIABLES:
    integer :: g
    integer :: icemec_class

    character(len=*), parameter :: subname = 'set_glc2lnd_fields'
    !-----------------------------------------------------------------------

    SHR_ASSERT((ubound(x2l, 2) == bounds%endg), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(index_x2l_Sg_ice_covered) == (/maxpatch_glcmec/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(index_x2l_Sg_topo) == (/maxpatch_glcmec/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(index_x2l_Flgg_hflx) == (/maxpatch_glcmec/)), errMsg(sourcefile, __LINE__))

    if (glc_present) then
       do g = bounds%begg, bounds%endg
          do icemec_class = 0, maxpatch_glcmec
             this%frac_grc(g,icemec_class)  = x2l(index_x2l_Sg_ice_covered(icemec_class),g)
             this%topo_grc(g,icemec_class)  = x2l(index_x2l_Sg_topo(icemec_class),g)
             this%hflx_grc(g,icemec_class)  = x2l(index_x2l_Flgg_hflx(icemec_class),g)
          end do
          this%icemask_grc(g)  = x2l(index_x2l_Sg_icemask,g)
          this%icemask_coupled_fluxes_grc(g)  = x2l(index_x2l_Sg_icemask_coupled_fluxes,g)
       end do

       call this%check_glc2lnd_icemask(bounds)
       call this%check_glc2lnd_icemask_coupled_fluxes(bounds)
       call this%update_glc2lnd_dyn_runoff_routing(bounds)
    else
       if (glc_do_dynglacier) then
          call endrun(' ERROR: With glc_present false (e.g., a stub glc model), glc_do_dynglacier must be false '// &
               errMsg(sourcefile, __LINE__))
       end if
    end if

  end subroutine set_glc2lnd_fields

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
       SHR_ASSERT_ALL((ubound(topo) == (/bounds%endg, maxpatch_glcmec/)), errMsg(sourcefile, __LINE__))
       this%topo_grc(bounds%begg:bounds%endg, 0:maxpatch_glcmec) = topo(bounds%begg:bounds%endg, 0:maxpatch_glcmec)
    end if

    if (present(icemask)) then
       SHR_ASSERT_ALL((ubound(icemask) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))
       this%icemask_grc(bounds%begg:bounds%endg) = icemask(bounds%begg:bounds%endg)
    end if

  end subroutine for_test_set_glc2lnd_fields_directly

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

          ! Ensure that icemask is a subset of has_virtual_columns. This is needed because
          ! we allocated memory based on has_virtual_columns, so it is a problem if the
          ! ice sheet tries to expand beyond the area defined by has_virtual_columns.
          if (.not. this%glc_behavior%has_virtual_columns_grc(g)) then
             write(iulog,'(a)') subname//' ERROR: icemask must be a subset of has_virtual_columns.'
             write(iulog,'(a)') 'Ensure that the glacier_region_behavior namelist item is set correctly.'
             write(iulog,'(a)') '(It should specify "virtual" for the region corresponding to the GLC domain.)'
             write(iulog,'(a)') 'If glacier_region_behavior is set correctly, then you can fix this problem'
             write(iulog,'(a)') 'by modifying GLACIER_REGION on the surface dataset.'
             write(iulog,'(a)') '(Expand the region that corresponds to the GLC domain'
             write(iulog,'(a)') '- i.e., the region specified as "virtual" in glacier_region_behavior.)'
             call endrun(decomp_index=g, clmlevel=nameg, msg=errMsg(sourcefile, __LINE__))
          end if

          ! Ensure that icemask is a subset of melt_replaced_by_ice. This is needed
          ! because we only compute SMB in the region given by melt_replaced_by_ice
          ! (according to the logic for building the do_smb filter), and we need SMB
          ! everywhere inside the icemask.
          if (.not. this%glc_behavior%melt_replaced_by_ice_grc(g)) then
             write(iulog,'(a)') subname//' ERROR: icemask must be a subset of melt_replaced_by_ice.'
             write(iulog,'(a)') 'Ensure that the glacier_region_melt_behavior namelist item is set correctly.'
             write(iulog,'(a)') '(It should specify "replaced_by_ice" for the region corresponding to the GLC domain.)'
             write(iulog,'(a)') 'If glacier_region_behavior is set correctly, then you can fix this problem'
             write(iulog,'(a)') 'by modifying GLACIER_REGION on the surface dataset.'
             write(iulog,'(a)') '(Expand the region that corresponds to the GLC domain'
             write(iulog,'(a)') '- i.e., the region specified as "replaced_by_ice" in glacier_region_melt_behavior.)'
             call endrun(decomp_index=g, clmlevel=nameg, msg=errMsg(sourcefile, __LINE__))
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

       ! Ensure that icemask_coupled_fluxes is a subset of icemask. Although there
       ! currently is no code in CLM that depends on this relationship, it seems helpful
       ! to ensure that this intuitive relationship holds, so that code developed in the
       ! future can rely on it.
       if (this%icemask_coupled_fluxes_grc(g) > 0._r8 .and. this%icemask_grc(g) == 0._r8) then
          write(iulog,*) subname//' ERROR: icemask_coupled_fluxes must be a subset of icemask.'
          call endrun(decomp_index=g, clmlevel=nameg, msg=errMsg(sourcefile, __LINE__))
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

    do g = bounds%begg, bounds%endg

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

       if (this%glc_dyn_runoff_routing_grc(g) > 0.0_r8) then

          ! Ensure that glc_dyn_runoff_routing is a subset of melt_replaced_by_ice. This
          ! is needed because glacial melt is only sent to the runoff stream in the region
          ! given by melt_replaced_by_ice (because the latter is used to create the do_smb
          ! filter, and the do_smb filter controls where glacial melt is computed).
          if (.not. this%glc_behavior%melt_replaced_by_ice_grc(g)) then
             write(iulog,'(a)') subname//' ERROR: icemask_coupled_fluxes must be a subset of melt_replaced_by_ice.'
             write(iulog,'(a)') 'Ensure that the glacier_region_melt_behavior namelist item is set correctly.'
             write(iulog,'(a)') '(It should specify "replaced_by_ice" for the region corresponding to the GLC domain.)'
             write(iulog,'(a)') 'If glacier_region_behavior is set correctly, then you can fix this problem'
             write(iulog,'(a)') 'by modifying GLACIER_REGION on the surface dataset.'
             write(iulog,'(a)') '(Expand the region that corresponds to the GLC domain'
             write(iulog,'(a)') '- i.e., the region specified as "replaced_by_ice" in glacier_region_melt_behavior.)'
             call endrun(decomp_index=g, clmlevel=nameg, msg=errMsg(sourcefile, __LINE__))
          end if
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
    use column_varcon     , only : col_itype_to_icemec_class
    use subgridWeightsMod , only : set_landunit_weight
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(in) :: this
    type(bounds_type)  , intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g,c                              ! indices
    real(r8):: area_ice_mec                     ! area of the ice_mec landunit
    integer :: l_ice_mec                        ! index of the ice_mec landunit
    integer :: icemec_class                     ! current icemec class (1..maxpatch_glcmec)
    logical :: frac_assigned(1:maxpatch_glcmec) ! whether this%frac has been assigned for each elevation class
    logical :: error                            ! if an error was found
    
    character(len=*), parameter :: subname = 'update_glc2lnd_fracs'
    !-----------------------------------------------------------------------

    if (glc_do_dynglacier) then
       do g = bounds%begg, bounds%endg
          ! Values from GLC are only valid within the icemask, so we only update CLM's areas there
          if (this%icemask_grc(g) > 0._r8) then

             ! Set total icemec landunit area
             area_ice_mec = sum(this%frac_grc(g, 1:maxpatch_glcmec))
             call set_landunit_weight(g, istice_mec, area_ice_mec)

             ! If new landunit area is greater than 0, then update column areas
             ! (If new landunit area is 0, col%wtlunit is arbitrary, so we might as well keep the existing values)
             if (area_ice_mec > 0) then
                ! Determine index of the glc_mec landunit
                l_ice_mec = grc%landunit_indices(istice_mec, g)
                if (l_ice_mec == ispval) then
                   write(iulog,*) subname//' ERROR: no ice_mec landunit found within the icemask, for g = ', g
                   call endrun()
                end if

                frac_assigned(1:maxpatch_glcmec) = .false.
                do c = lun%coli(l_ice_mec), lun%colf(l_ice_mec)
                   icemec_class = col_itype_to_icemec_class(col%itype(c))
                   col%wtlunit(c) = this%frac_grc(g, icemec_class) / lun%wtgcell(l_ice_mec)
                   frac_assigned(icemec_class) = .true.
                end do

                ! Confirm that all elevation classes that have non-zero area according to
                ! this%frac have been assigned to a column in CLM's data structures
                error = .false.
                do icemec_class = 1, maxpatch_glcmec
                   if (this%frac_grc(g, icemec_class) > 0._r8 .and. &
                        .not. frac_assigned(icemec_class)) then
                      error = .true.
                   end if
                end do
                if (error) then
                   write(iulog,*) subname//' ERROR: at least one glc_mec column has non-zero area from the coupler,'
                   write(iulog,*) 'but there was no slot in memory for this column; g = ', g
                   write(iulog,*) 'this%frac_grc(g, 1:maxpatch_glcmec) = ', &
                        this%frac_grc(g, 1:maxpatch_glcmec)
                   write(iulog,*) 'frac_assigned(1:maxpatch_glcmec) = ', &
                        frac_assigned(1:maxpatch_glcmec)
                   call endrun()
                end if  ! error
             end if  ! area_ice_mec > 0
          end if  ! this%icemask_grc(g) > 0
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
    use landunit_varcon , only : istice_mec
    use column_varcon   , only : col_itype_to_icemec_class
    !
    ! !ARGUMENTS:
    class(glc2lnd_type) , intent(in)    :: this
    type(bounds_type)   , intent(in)    :: bounds                   ! bounds
    real(r8)            , intent(inout) :: topo_col( bounds%begc: ) ! topographic height (m)
    logical             , intent(inout) :: needs_downscaling_col( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    integer :: c, l, g      ! indices
    integer :: icemec_class ! current icemec class (1..maxpatch_glcmec)

    character(len=*), parameter :: subname = 'update_glc2lnd_topo'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(topo_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(needs_downscaling_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    if (glc_do_dynglacier) then
       do c = bounds%begc, bounds%endc
          l = col%landunit(c)
          g = col%gridcell(c)

          ! Values from GLC are only valid within the icemask, so we only update CLM's topo values there
          if (this%icemask_grc(g) > 0._r8) then
             if (lun%itype(l) == istice_mec) then
                icemec_class = col_itype_to_icemec_class(col%itype(c))
             else
                ! If not on a glaciated column, assign topography to the bare-land value determined by GLC.
                icemec_class = 0
             end if

             ! Note that we do downscaling over all column types. This is for consistency:
             ! interpretation of results would be difficult if some non-glacier column types
             ! were downscaled but others were not.
             !
             ! BUG(wjs, 2016-11-15, bugz 2377) Actually, do not downscale over urban points:
             ! this currently isn't allowed because the urban code references some
             ! non-downscaled, gridcell-level atmospheric forcings
             if (.not. lun%urbpoi(l)) then
                topo_col(c) = this%topo_grc(g, icemec_class)
                needs_downscaling_col(c) = .true.
             end if
          end if
       end do
    end if

  end subroutine update_glc2lnd_topo

end module glc2lndMod

