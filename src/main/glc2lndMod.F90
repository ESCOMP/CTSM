module glc2lndMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle arrays used for exchanging data from glc to clm.
  ! For now glc datais send and received on the lnd decomposition and grid.
  !
  ! The fields sent from the lnd component to the glc component via
  !  the coupler are labeled 's2x', or sno to coupler.
  ! The fields received by the lnd component from the glc component
  !  via the coupler are labeled 'x2s', or coupler to sno.
  ! 'Sno' is a misnomer in that the exchanged data are related to
  !  the ice beneath the snow, not the snow itself.  But by CESM convention,
  ! 'ice' refers to sea ice, not land ice.
  !
  ! !USES:
#include "shr_assert.h"
  use decompMod      , only : bounds_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : maxpatch_glcmec
  use clm_varctl     , only : iulog
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

     real(r8), pointer :: frac_grc    (:,:) => null()
     real(r8), pointer :: topo_grc    (:,:) => null()
     real(r8), pointer :: hflx_grc    (:,:) => null()

     ! TODO(wjs, 2016-04-01) If the setting of icemask and icemask_coupled_fluxes were
     ! moved out of lnd_import_export into this module, then these two variables could be
     ! made private.

     ! Area in which GLC model can accept surface mass balance, received from glc (0-1)
     real(r8), pointer :: icemask_grc (:)   => null()

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
     real(r8), pointer :: icemask_coupled_fluxes_grc (:)  => null()

     ! Where we should do runoff routing that is appropriate for having a dynamic icesheet underneath.
     logical , pointer :: glc_dyn_runoff_routing_grc (:) => null()

   contains

     ! ------------------------------------------------------------------------
     ! Public routines
     ! ------------------------------------------------------------------------

     procedure, public  :: Init
     procedure, public  :: Clean
     procedure, public  :: update_glc2lnd_non_topo  ! update everything except topographic heights
     procedure, public  :: update_glc2lnd_topo      ! update topographic heights

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

     ! update subgrid fractions based on input from GLC
     procedure, private :: update_glc2lnd_fracs

  end type glc2lnd_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, glc_behavior)

    class(glc2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    type(glc_behavior_type), intent(in) :: glc_behavior

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
    allocate(this%glc_dyn_runoff_routing_grc (begg:endg))    ;   this%glc_dyn_runoff_routing_grc (:)   = .false.

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
    
    if (maxpatch_glcmec > 0) then
       this%icemask_grc(begg:endg) = spval
       call hist_addfld1d (fname='ICE_MODEL_FRACTION',  units='unitless',  &
            avgflag='I', long_name='Ice sheet model fractional coverage', &
            ptr_gcell=this%icemask_grc)
    end if

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
    type(glc_behavior_type), intent(in) :: glc_behavior
    !
    ! !LOCAL VARIABLES:
    integer :: begg, endg
    integer :: g

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    begg = bounds%begg
    endg = bounds%endg
    
    this%frac_grc(begg:endg, :) = 0.0_r8
    this%topo_grc(begg:endg, :) = 0.0_r8
    this%hflx_grc(begg:endg, :) = 0.0_r8

    ! Since we don't have GLC's icemask yet in initialization, we use has_virtual_columns
    ! as a rough initial guess. Note that has_virtual_columns is guaranteed to be a
    ! superset of the icemask.
    do g = begg, endg
       if (glc_behavior%has_virtual_columns_grc(g)) then
          this%icemask_grc(g) = 1._r8
       else
          this%icemask_grc(g) = 0._r8
       end if
    end do

    ! initialize icemask_coupled_fluxes to 0; this seems safest in case we aren't coupled
    ! to CISM (to ensure that we use the uncoupled form of runoff routing)
    this%icemask_coupled_fluxes_grc(begg:endg) = 0.0_r8
    
    call this%update_glc2lnd_dyn_runoff_routing(bounds, glc_behavior)

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
  subroutine update_glc2lnd_non_topo(this, bounds, glc_behavior)
    !
    ! !DESCRIPTION:
    ! Update values to derived-type CLM variables based on input from GLC (via the coupler)
    !
    ! This does NOT update topographic heights: those are updated in the separate
    ! update_glc2lnd_topo routine
    !
    ! icemask, icemask_coupled_fluxes, and glc_dyn_runoff_routing are always updated
    ! (although note that this routine should only be called when
    ! create_glacier_mec_landunit is true, or some similar condition; this should be
    ! controlled in a conditional around the call to this routine); fracs are updated if
    ! glc_do_dynglacier is true
    !
    ! !USES:
    use clm_varctl , only : glc_do_dynglacier
    !
    ! !ARGUMENTS:
    class(glc2lnd_type)     , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    type(glc_behavior_type) , intent(in)    :: glc_behavior
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'update_glc2lnd_non_topo'
    !-----------------------------------------------------------------------

    ! Note that nothing is needed to update icemask or icemask_coupled_fluxes here,
    ! because these values have already been set in lnd_import_export. However, we do
    ! some sanity-checking of those fields here.
    call this%check_glc2lnd_icemask(bounds, glc_behavior)
    call this%check_glc2lnd_icemask_coupled_fluxes(bounds)

    call this%update_glc2lnd_dyn_runoff_routing(bounds, glc_behavior)

    if (glc_do_dynglacier) then
       call this%update_glc2lnd_fracs(bounds)
    end if

  end subroutine update_glc2lnd_non_topo

  !-----------------------------------------------------------------------
  subroutine check_glc2lnd_icemask(this, bounds, glc_behavior)
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
    type(glc_behavior_type) , intent(in) :: glc_behavior
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
          if (.not. glc_behavior%has_virtual_columns_grc(g)) then
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
          if (.not. glc_behavior%melt_replaced_by_ice_grc(g)) then
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
  subroutine update_glc2lnd_dyn_runoff_routing(this, bounds, glc_behavior)
    !
    ! !DESCRIPTION:
    ! Update glc_dyn_runoff_routing field based on updated icemask_coupled_fluxes field
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(glc2lnd_type), intent(inout) :: this
    type(bounds_type)  , intent(in) :: bounds ! bounds
    type(glc_behavior_type) , intent(in) :: glc_behavior
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! grid cell index
    
    character(len=*), parameter :: subname = 'update_glc2lnd_dyn_runoff_routing'
    !-----------------------------------------------------------------------

    ! Wherever we have an icesheet that is computing and sending fluxes to the coupler -
    ! which particularly means it is computing a calving flux - we will use the
    ! "glc_dyn_runoff_routing" scheme. In other places - including places where CISM is
    ! not running at all, as well as places where CISM is running in diagnostic-only mode
    ! and therefore is not sending a calving flux - we use the alternative
    ! glc_dyn_runoff_routing=false scheme. This is needed to conserve water correctly in
    ! the absence of a calving flux.

    do g = bounds%begg, bounds%endg
       if (this%icemask_coupled_fluxes_grc(g) > 0._r8) then
          this%glc_dyn_runoff_routing_grc(g) = .true.
       else
          this%glc_dyn_runoff_routing_grc(g) = .false.
       end if

       if (this%glc_dyn_runoff_routing_grc(g)) then

          ! Ensure that glc_dyn_runoff_routing is a subset of melt_replaced_by_ice. This
          ! is needed because glacial melt is only sent to the runoff stream in the region
          ! given by melt_replaced_by_ice (because the latter is used to create the do_smb
          ! filter, and the do_smb filter controls where glacial melt is computed).
          if (.not. glc_behavior%melt_replaced_by_ice_grc(g)) then
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
    ! This should only be called when create_glacier_mec_landunit is true, or some
    ! similar condition (this should be controlled in a conditional around the call to
    ! this routine).
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

  end subroutine update_glc2lnd_topo

end module glc2lndMod

