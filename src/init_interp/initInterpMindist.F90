module initInterpMindist

  ! ------------------------------------------------------------------------
  ! This module contains most of the "interesting" logic of initInterp, in terms of
  ! finding the input column (or landunit, patch, etc.) to use as a template for each
  ! output column (etc.).
  !
  ! This is in a separate module to facilitate unit testing, since the full initInterp
  ! involves some awkward dependencies.
  ! ------------------------------------------------------------------------

  use shr_kind_mod   , only: r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varctl     , only: iulog
  use abortutils     , only: endrun
  use spmdMod        , only: masterproc
  use clm_varcon     , only: spval, re

  implicit none
  private
  save

  ! Public methods

  public :: set_mindist

  ! Public types

  type, public :: subgrid_special_indices_type
     integer :: ipft_not_vegetated
     integer :: icol_vegetated_or_bare_soil
     integer :: ilun_vegetated_or_bare_soil
     integer :: ilun_crop
     integer :: ilun_landice_multiple_elevation_classes
   contains
     procedure :: is_vegetated_landunit  ! returns true if the given landunit type is natural veg or crop
  end type subgrid_special_indices_type

  type, public :: subgrid_type
     character(len=16) :: name               ! pft, column, landunit, gridcell
     integer , pointer :: ptype(:) => null() ! used for patch type 
     integer , pointer :: ctype(:) => null() ! used for patch or col type
     integer , pointer :: ltype(:) => null() ! used for pft, col or lun type
     real(r8), pointer :: topoglc(:) => null()
     real(r8), pointer :: lat(:)
     real(r8), pointer :: lon(:)
     real(r8), pointer :: coslat(:)
   contains
     procedure :: print_point  ! print info about one point
  end type subgrid_type

  ! Private methods

  private :: do_fill_missing_with_natveg
  private :: is_sametype
  private :: is_baresoil

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine print_point(this, index, unit)
    !
    ! !DESCRIPTION:
    ! Print info about one point in a subgrid_type object
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(subgrid_type), intent(in) :: this
    integer            , intent(in) :: index
    integer            , intent(in) :: unit ! unit to which we should write the info
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'print_point'
    !-----------------------------------------------------------------------

    write(unit,*) 'subgrid level, index = ',&
         this%name, index
    if (associated(this%ltype)) then
       write(unit,*) 'ltype: ', this%ltype(index)
    end if
    if (associated(this%ctype)) then
       write(unit,*) 'ctype: ', this%ctype(index)
    end if
    if (associated(this%ptype)) then
       write(unit,*) 'ptype: ', this%ptype(index)
    end if
    
  end subroutine print_point


  !=======================================================================

  subroutine set_mindist(begi, endi, bego, endo, activei, activeo, subgridi, subgrido, &
       subgrid_special_indices, fill_missing_with_natveg, mindist_index)

    ! --------------------------------------------------------------------
    ! arguments
    integer            , intent(in)  :: begi, endi 
    integer            , intent(in)  :: bego, endo 
    logical            , intent(in)  :: activei(begi:endi) 
    logical            , intent(in)  :: activeo(bego:endo) 
    type(subgrid_type) , intent(in)  :: subgridi
    type(subgrid_type) , intent(in)  :: subgrido
    type(subgrid_special_indices_type), intent(in) :: subgrid_special_indices

    ! If false: if an output type cannot be found in the input, code aborts
    ! If true: if an output type cannot be found in the input, fill with closest natural
    ! veg column (using bare soil for patch-level variables)
    !
    ! NOTE: always treated as true for natural veg and crop landunits/columns/patches in
    ! the output - e.g., if we can't find the right column type to fill crop, we always
    ! use the closest natural veg column, regardless of the value of this flag.
    logical            , intent(in)  :: fill_missing_with_natveg

    integer            , intent(out) :: mindist_index(bego:endo)
    !
    ! local variables
    real(r8) :: dx,dy
    real(r8) :: distmin,dist,hgtdiffmin,hgtdiff    
    integer  :: nsizei, nsizeo
    integer  :: ni,no,nmin,ier,n,noloc
    logical  :: closest
    ! --------------------------------------------------------------------

    mindist_index(bego:endo) = 0
    distmin = spval

!$OMP PARALLEL DO PRIVATE (ni,no,n,nmin,distmin,dx,dy,dist,closest,hgtdiffmin,hgtdiff)
    do no = bego,endo

       ! Only interpolate onto active points. Otherwise, the mere act of running
       ! init_interp (e.g., of a file onto itself) would lead to changes in a bunch of
       ! inactive points - e.g., going from their cold start initial conditions to some
       ! spunup initial conditions (from the closest active point of that type). This
       ! could potentially lead to different behavior in a transient run, if those points
       ! later became active; that's undesirable.
       if (activeo(no)) then 

          nmin    = 0
          distmin = spval
          hgtdiffmin = spval
          do ni = begi,endi
             if (activei(ni)) then
                if (is_sametype(ni, no, subgridi, subgrido, subgrid_special_indices)) then
                   dy = abs(subgrido%lat(no)-subgridi%lat(ni))*re
                   dx = abs(subgrido%lon(no)-subgridi%lon(ni))*re * &
                        0.5_r8*(subgrido%coslat(no)+subgridi%coslat(ni))
                   dist = dx*dx + dy*dy
                   if (associated(subgridi%topoglc) .and. associated(subgrido%topoglc)) then
                      hgtdiff = abs(subgridi%topoglc(ni) - subgrido%topoglc(no))
                   end if
                   closest = .false.
                   if ( dist < distmin ) then
                      closest = .true.
                      distmin = dist
                      nmin = ni
                      if (associated(subgridi%topoglc) .and. associated(subgrido%topoglc)) then
                         hgtdiffmin = hgtdiff
                      end if
                   end if
                   if (.not. closest) then
                      ! For glc_mec points, we first find the closest point in lat-lon
                      ! space (above). Then, within that closest point, we find the
                      ! closest column in topographic space; this second piece is done
                      ! here. Note that this ordering means that we could choose a column
                      ! with a very different topographic height from the target column,
                      ! if it is closer in lat-lon space than any similar-height columns.
                      if (associated(subgridi%topoglc) .and. associated(subgrido%topoglc)) then
                         hgtdiff = abs(subgridi%topoglc(ni) - subgrido%topoglc(no))
                         if ((dist == distmin) .and. (hgtdiff < hgtdiffmin)) then
                            closest = .true.
                            hgtdiffmin = hgtdiff
                            distmin = dist
                            nmin = ni
                         end if
                      end if
                   end if
                end if
             end if
          end do
          
          ! If output type is not contained in input dataset, then use closest bare soil,
          ! if this point is one for which we fill missing with natveg.
          if ( distmin == spval .and. &
               do_fill_missing_with_natveg( &
               fill_missing_with_natveg, no, subgrido, subgrid_special_indices)) then
             do ni = begi, endi
                if (activei(ni)) then
                   if ( is_baresoil(ni, subgridi, subgrid_special_indices)) then
                      dy = abs(subgrido%lat(no)-subgridi%lat(ni))*re
                      dx = abs(subgrido%lon(no)-subgridi%lon(ni))*re * &
                           0.5_r8*(subgrido%coslat(no)+subgridi%coslat(ni))
                      dist = dx*dx + dy*dy
                      if ( dist < distmin )then
                         distmin = dist
                         nmin = ni
                      end if
                   end if
                end if
             end do
          end if

          ! Error conditions
          if ( distmin == spval )then
             write(iulog,*) 'ERROR initInterp set_mindist: &
                  &Cannot find any input points matching output point:'
             call subgrido%print_point(no, iulog)
             write(iulog,*) ' '
             write(iulog,*) 'Consider rerunning with the following in user_nl_clm:'
             write(iulog,*) 'init_interp_fill_missing_with_natveg = .true.'
             write(iulog,*) 'However, note that this will fill all missing types in the output'
             write(iulog,*) 'with the closest natural veg column in the input'
             write(iulog,*) '(using bare soil for patch-level variables).'
             write(iulog,*) 'So, you should consider whether that is what you want.'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if

          mindist_index(no) = nmin

       end if ! end if activeo block
    end do
!$OMP END PARALLEL DO
    
  end subroutine set_mindist

  !-----------------------------------------------------------------------
  function do_fill_missing_with_natveg(fill_missing_with_natveg, &
       no, subgrido, subgrid_special_indices)
    !
    ! !DESCRIPTION:
    ! Returns true if the given output point, if missing, should be filled with the
    ! closest natural veg point.
    !
    ! !ARGUMENTS:
    logical :: do_fill_missing_with_natveg  ! function result

    ! whether we should fill ALL missing points with natveg
    logical, intent(in) :: fill_missing_with_natveg

    integer           , intent(in)  :: no
    type(subgrid_type), intent(in)  :: subgrido
    type(subgrid_special_indices_type), intent(in) :: subgrid_special_indices
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'do_fill_missing_with_natveg'
    !-----------------------------------------------------------------------

    if (subgrido%name == 'gridcell') then
       ! It makes no sense to try to fill missing with natveg for gridcell-level values
       do_fill_missing_with_natveg = .false.
    else if (fill_missing_with_natveg) then
       ! User has asked for all missing points to be filled with natveg
       do_fill_missing_with_natveg = .true.
    else if (subgrid_special_indices%is_vegetated_landunit(subgrido%ltype(no))) then
       ! Even if user hasn't asked for it, we fill missing vegetated points (natural veg
       ! and crop) with the closest natveg point. This is mainly to support the common
       ! use case of interpolating non-crop to crop, but also supports adding a new PFT
       ! type.
       do_fill_missing_with_natveg = .true.
    else
       do_fill_missing_with_natveg = .false.
    end if

  end function do_fill_missing_with_natveg


  !=======================================================================

  logical function is_sametype (ni, no, subgridi, subgrido, subgrid_special_indices)

    ! --------------------------------------------------------------------
    ! arguments
    integer           , intent(in)  :: ni 
    integer           , intent(in)  :: no 
    type(subgrid_type), intent(in)  :: subgridi
    type(subgrid_type), intent(in)  :: subgrido
    type(subgrid_special_indices_type), intent(in) :: subgrid_special_indices
    ! --------------------------------------------------------------------

    is_sametype = .false.

    if (trim(subgridi%name) == 'pft' .and. trim(subgrido%name) == 'pft') then
       if ( subgridi%ltype(ni) == subgrid_special_indices%ilun_landice_multiple_elevation_classes .and. &
            subgrido%ltype(no) == subgrid_special_indices%ilun_landice_multiple_elevation_classes) then
          is_sametype = .true.
       else if (subgrid_special_indices%is_vegetated_landunit(subgrido%ltype(no))) then
          ! If the output type is natural veg or crop, then just look for the correct PFT,
          ! without regard for what column or landunit it's on (as long as it's on either
          ! the natural veg or crop landunit). This is needed to handle the generic crop
          ! properly when interpolating from non-crop to crop, or vice versa.
          !
          ! TODO(wjs, 2015-09-15) If we ever allow the same PFT to appear on multiple
          ! columns within a given grid cell, then this logic will need to be made
          ! somewhat more complex: e.g., preferably take something from the same column
          ! type, but if we can't find anything from the same column type, then ignore
          ! column type.

          if (subgrid_special_indices%is_vegetated_landunit(subgridi%ltype(ni)) .and. &
               subgridi%ptype(ni) == subgrido%ptype(no)) then
             is_sametype = .true.
          end if
       else if (subgridi%ptype(ni) == subgrido%ptype(no) .and. &
                subgridi%ctype(ni) == subgrido%ctype(no) .and. &
                subgridi%ltype(ni) == subgrido%ltype(no)) then
          is_sametype = .true.
       end if
    else if (trim(subgridi%name) == 'column' .and. trim(subgrido%name) == 'column') then
       if ( subgridi%ltype(ni) == subgrid_special_indices%ilun_landice_multiple_elevation_classes  .and. &
            subgrido%ltype(no) == subgrid_special_indices%ilun_landice_multiple_elevation_classes ) then
          is_sametype = .true.
       else if (subgridi%ctype(ni) == subgrido%ctype(no) .and. &
                subgridi%ltype(ni) == subgrido%ltype(no)) then
          is_sametype = .true.
       end if
    else if (trim(subgridi%name) == 'landunit' .and. trim(subgrido%name) == 'landunit') then
       if (subgridi%ltype(ni) == subgrido%ltype(no)) then
          is_sametype = .true.
       end if
    else if (trim(subgridi%name) == 'gridcell' .and. trim(subgrido%name) == 'gridcell') then
       is_sametype = .true.
    else 
       if (masterproc) then
          write(iulog,*)'ERROR interpinic: is_sametype check on input and output type not supported'
          write(iulog,*)'typei = ',trim(subgridi%name)
          write(iulog,*)'typeo = ',trim(subgrido%name)
       end if
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

  end function is_sametype

  !=======================================================================

  logical function is_baresoil (n, subgrid, subgrid_special_indices)

    ! --------------------------------------------------------------------
    ! arguments
    integer           , intent(in)  :: n 
    type(subgrid_type), intent(in)  :: subgrid
    type(subgrid_special_indices_type), intent(in) :: subgrid_special_indices
    ! --------------------------------------------------------------------

    is_baresoil = .false.

    if (subgrid%name == 'pft') then
       if (subgrid%ptype(n) == subgrid_special_indices%ipft_not_vegetated .and. &
            subgrid%ctype(n) == subgrid_special_indices%icol_vegetated_or_bare_soil .and. &
            subgrid%ltype(n) == subgrid_special_indices%ilun_vegetated_or_bare_soil) then
          is_baresoil = .true.
       end if
    else if (subgrid%name == 'column') then
       if (subgrid%ctype(n) == subgrid_special_indices%icol_vegetated_or_bare_soil .and. &
            subgrid%ltype(n) == subgrid_special_indices%ilun_vegetated_or_bare_soil) then
          is_baresoil = .true.
       end if
    else if (subgrid%name == 'landunit') then
       if (subgrid%ltype(n) == subgrid_special_indices%ilun_vegetated_or_bare_soil) then
          is_baresoil = .true.
       end if
    else 
       if (masterproc) then
          write(iulog,*)'ERROR interpinic: is_baresoil subgrid type ',subgrid%name,' not supported'
       end if
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

  end function is_baresoil

  !-----------------------------------------------------------------------
  function is_vegetated_landunit(this, ltype)
    !
    ! !DESCRIPTION:
    ! Returns true if the given landunit type is vegetated: either natural veg or crop
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    logical :: is_vegetated_landunit  ! function result
    class(subgrid_special_indices_type), intent(in) :: this
    integer, intent(in) :: ltype   ! landunit type of interest
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'is_vegetated_landunit'
    !-----------------------------------------------------------------------

    if (ltype == this%ilun_vegetated_or_bare_soil .or. &
         ltype == this%ilun_crop) then
       is_vegetated_landunit = .true.
    else
       is_vegetated_landunit = .false.
    end if

  end function is_vegetated_landunit


end module initInterpMindist
