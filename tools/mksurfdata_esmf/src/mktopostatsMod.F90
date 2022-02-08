module mktopostatsMod

  !-----------------------------------------------------------------------
  ! make various topography statistics
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use mkdomainMod , only : domain_checksame

  implicit none
  private

  public :: mktopostats            ! make topo stddev & mean slope

  type(ESMF_DynamicMask) :: dynamicMask

!===============================================================
contains
!===============================================================

  subroutine mktopostats(ldomain, mapfname, datfname, ndiag, topo_stddev_o, slope_o, std_elev)
    !
    ! make various topography statistics
    !
    use mkdiagnosticsMod, only : output_diagnostics_continuous, output_diagnostics_continuous_outonly
    use mkchecksMod, only : min_bad, max_bad
    !
    ! !ARGUMENTS:
    type(domain_type) , intent(in) :: ldomain
    character(len=*)  , intent(in) :: mapfname          ! input mapping file name
    character(len=*)  , intent(in) :: datfname          ! input data file name
    integer           , intent(in) :: ndiag             ! unit number for diag out
    real(r8)          , intent(in) :: std_elev          ! standard deviation of elevation (m) to use when not using input file
    real(r8)          , intent(out):: topo_stddev_o(:)  ! output grid: standard deviation of elevation (m)
    real(r8)          , intent(out):: slope_o(:)        ! output grid: slope (degrees)
    !
    ! !LOCAL VARIABLES:
    type(gridmap_type)    :: tgridmap
    type(domain_type)     :: tdomain            ! local domain
    real(r8), allocatable :: data_i(:)          ! data on input grid
    integer  :: ncid,varid                      ! input netCDF id's
    integer  :: ier                             ! error status
    logical  :: bypass_reading                  ! If should bypass reading dataset and just use a global value

    real(r8), parameter :: min_valid_topo_stddev = 0._r8

    real(r8), parameter :: min_valid_slope = 0._r8
    real(r8), parameter :: max_valid_slope = 90._r8

    character(len=32) :: subname = 'mktopostats'
    !-----------------------------------------------------------------------

    write (6,*) 'Attempting to make Topography statistics.....'
    if ( std_elev >= 0.0_r8 )then
       bypass_reading = .true.
       write (6,*) '    By pass the reading and just use global values'
    else
       bypass_reading = .false.
    end if

    ! Create a dynamic mask object
    ! The dynamic mask object further holds a pointer to the routine that will be called in order to
    ! handle dynamically masked elements - in this case its DynMaskProc (see below)
    call ESMF_DynamicMaskSetR8R8R8(dynamicMask, dynamicSrcMaskValue=czero, dynamicMaskRoutine=DynMaskProc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle - the normal way

    ! Do the regridding the normal way to create the mean

    ! Compute the standard deviation
    ! Use the dynamic masking could be used to do this. 
    ! dst_array(no) = dst_array(no) + wt * (src_array(ni) - weighted_means(no))**2 part,
    ! and then just do a plain sparse matrix multiply on a src Field of all 1.0 to do the
    ! weight sum part, and then do the divide and sqrt locally on each PET.

    !  One issue is how to get the weight_means() in for each
    !  dest. location. I think that you could pass them in via the
    !  dst_array and then pull them out before doing the calculation, but
    !  to do so that you to stop it from zeroing out the
    !  dst array which you can do by setting
    !  zeroregion=ESMF_REGION_EMPTY.This would also depend on the
    !  calculation happening only once for each destination location,
    !  which I would guess is true, but Gerhard can confirm.

    ! -----------------------------------------------------------------
    ! Open input file, allocate memory for input data
    ! -----------------------------------------------------------------

    write(6,*)'Open Topography file: ', trim(datfname)
    call check_ret(nf_open(datfname, 0, ncid), subname)

    allocate(data_i(tdomain%ns), stat=ier)
    if (ier/=0) call abort()

    ! -----------------------------------------------------------------
    ! Make topography standard deviation
    ! -----------------------------------------------------------------

    call check_ret(nf_inq_varid (ncid, 'ELEVATION', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
    call gridmap_areastddev(tgridmap, data_i, topo_stddev_o, nodata=0._r8)

    call output_diagnostics_continuous_outonly(topo_stddev_o, tgridmap, "Topo Std Dev", "m", ndiag)
 else
    write (6,*) '    Set std deviation of topography to ', std_elev
    topo_stddev_o = std_elev
 end if

 ! Check validity of output data
 if (min_bad(topo_stddev_o, min_valid_topo_stddev, 'topo_stddev')) then
    call abort()
 end if


 ! -----------------------------------------------------------------
 ! Regrid slope
 ! -----------------------------------------------------------------

 if ( .not. bypass_reading )then
    call check_ret(nf_inq_varid (ncid, 'SLOPE', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, data_i), subname)

    ! Subr. gridmap_areaave_no_srcmask should NOT be used in general. We have
    ! kept it to support the rare raw data files for which we have masking on
    ! the mapping file and, therefore, we do not explicitly pass the src_mask
    ! as an argument. In general, users are advised to use subroutine
    ! gridmap_areaave_srcmask.
    call gridmap_areaave_no_srcmask(tgridmap, data_i, slope_o, nodata=0._r8)

    call output_diagnostics_continuous(data_i, slope_o, tgridmap, "Slope", "degrees", ndiag, tdomain%mask, tgridmap%frac_dst)
 else
    write (6,*) '    Set slope of topography to ', 0.0_r8
    slope_o = 0.0_r8
 end if
 ! Check validity of output data
 if (min_bad(slope_o, min_valid_slope, 'slope') .or. &
      max_bad(slope_o, max_valid_slope, 'slope')) then
    call abort()
 end if


 ! -----------------------------------------------------------------
 ! Close files and deallocate dynamic memory
 ! -----------------------------------------------------------------

 if ( .not. bypass_reading )then
    call check_ret(nf_close(ncid), subname)
    call domain_clean(tdomain) 
    call gridmap_clean(tgridmap)
    deallocate (data_i)
 end if

 write (6,*) 'Successfully made Topography statistics'
 write (6,*)

end subroutine mktopostats

  !================================================================================================
  subroutine dynMaskProc(dynamicMaskList, dynamicSrcMaskValue, dynamicDstMaskValue, rc)

    use ESMF, only : ESMF_RC_ARG_BAD

    ! input/output arguments
    type(ESMF_DynamicMaskElementR8R8R8) , pointer              :: dynamicMaskList(:)
    real(ESMF_KIND_R8)                  , intent(in), optional :: dynamicSrcMaskValue
    real(ESMF_KIND_R8)                  , intent(in), optional :: dynamicDstMaskValue
    integer                             , intent(out)          :: rc

    ! local variables
    integer  :: i, j
    real(ESMF_KIND_R8)  :: renorm
    !---------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Below - ONLY if you do NOT have the source masked out then do
    ! the regridding (which is done explicitly here)

    if (associated(dynamicMaskList)) then
       do i=1, size(dynamicMaskList)
          dynamicMaskList(i)%dstElement = czero ! set to zero
          renorm = 0.d0 ! reset
          if (dynamicSrcMaskValue /= dynamicMaskList(i)%srcElement(j)) then
             do j = 1, size(dynamicMaskList(i)%factor)
                dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement + &
                     (dynamicMaskList(i)%factor(j) * dynamicMaskList(i)%srcElement(j))
                renorm = renorm + dynamicMaskList(i)%factor(j)
          enddo
          if (renorm > 0.d0) then
             dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement / renorm
          else if (present(dynamicSrcMaskValue)) then
             dynamicMaskList(i)%dstElement = dynamicSrcMaskValue
          else
             rc = ESMF_RC_ARG_BAD  ! error detected
             return
          endif
       enddo
    endif

  end subroutine DynOcnMaskProc

end module mktopostatsMod
