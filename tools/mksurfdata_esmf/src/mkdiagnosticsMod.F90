module mkdiagnosticsMod

  !-----------------------------------------------------------------------
  ! Output diagnostics to log file
  !-----------------------------------------------------------------------

  use ESMF
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_abort
  use mkutilsMod   , only : chkerr
  use mkvarctl     , only : mpicom, root_task

  implicit none
  private

#include <mpif.h>

  public :: output_diagnostics_area               ! output diagnostics for field that is % of grid area
  public :: output_diagnostics_continuous         ! output diagnostics for a continuous (real-valued) field
  public :: output_diagnostics_continuous_outonly ! output diagnostics for a continuous (real-valued) field
                                                  ! just on the output grid
  public :: output_diagnostics_index              ! output diagnostics for an index field

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine output_diagnostics_area(area_i, area_o, mask_i, data_i, data_o, name, percent, ndiag)

    ! Output diagnostics for a field that gives either fraction or percent of grid cell area

    use mkvarpar, only : re

    ! input/output variables
    real(r8)         , intent(in) :: area_i(:) 
    real(r8)         , intent(in) :: area_o(:) 
    integer          , intent(in) :: mask_i(:)
    real(r8)         , intent(in) :: data_i(:)    ! data on input grid
    real(r8)         , intent(in) :: data_o(:)    ! data on output grid
    character(len=*) , intent(in) :: name         ! name of field
    logical          , intent(in) :: percent      ! is field specified as percent? (alternative is fraction)
    integer          , intent(in) :: ndiag

    ! local variables:
    real(r8) :: loc_gdata_i ! local_global sum of input data
    real(r8) :: loc_gdata_o ! local_global sum of output data
    real(r8) :: gdata_i     ! global sum of input data
    real(r8) :: gdata_o     ! global sum of output data
    real(r8) :: loc_garea_i ! local global sum of input area
    real(r8) :: loc_garea_o ! local global sum of output area
    real(r8) :: garea_i     ! global sum of input area
    real(r8) :: garea_o     ! global sum of output area
    integer  :: ns_i, ns_o  ! sizes of input & output grids
    integer  :: ni,no,k     ! indices
    integer  :: ier         ! error code
    character(len=*), parameter :: subname = "output_diagnostics_area"
    !------------------------------------------------------------------------------

    ns_i = size(data_i)
    ns_o = size(data_o)

    ! Error check for array size consistencies
    if (size(mask_i) /= ns_i) then
       write(6,*) subname//' ERROR: incorrect size of mask_i'
       write(6,*) 'size(mask_i) = ', size(mask_i)
       write(6,*) 'ns_i = ', ns_i
       call shr_sys_abort()
    end if

    ! Sums on input grid
    loc_gdata_i = 0.
    loc_garea_i = 0.
    do ni = 1,ns_i
       loc_garea_i = loc_garea_i + area_i(ni) * re**2
       loc_gdata_i = loc_gdata_i + data_i(ni) * area_i(ni) * mask_i(ni) * re**2
    end do
    call mpi_reduce(loc_gdata_i,gdata_i,1,MPI_REAL8,MPI_SUM,0,mpicom,ier)

    ! Sums on output grid
    loc_gdata_o = 0.
    loc_garea_o = 0.
    do no = 1,ns_o
       loc_garea_o = loc_garea_o + area_o(no) * re**2
       loc_gdata_o = loc_gdata_o + data_o(no) * area_o(no) * re**2
    end do
    call mpi_reduce(loc_gdata_o,gdata_o,1,MPI_REAL8,MPI_SUM,0,mpicom,ier)

    ! Correct units
    if (percent) then
       gdata_i = gdata_i / 100._r8
       gdata_o = gdata_o / 100._r8
    end if

    ! Diagnostic output
    if (root_task) then
       write(ndiag,*)
       write(ndiag,*)
       write(ndiag,'(1x,70a1)') ('.',k=1,70)
       write(ndiag,'(a)') ' diagnostics for '//trim(name)
    write (ndiag,201)
201    format (1x,'surface type   input grid area  output grid area'/ &
               1x,'                 10**6 km**2      10**6 km**2   ')
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,*)
       write (ndiag,202) name,          gdata_i*1.e-06, gdata_o*1.e-06
       write (ndiag,202) 'all surface', garea_i*1.e-06, garea_o*1.e-06
202    format (1x,a12,           f14.3,f17.3)
    end if

  end subroutine output_diagnostics_area

  !===============================================================
  subroutine output_diagnostics_continuous(area_i, area_o, mask_i, frac_o, &
       data_i, data_o, name, units, ndiag)

    ! Output diagnostics for a continuous field (but not area, for
    ! which there is a different routine)

    use mkvarpar, only : re

    ! !ARGUMENTS:
    real(r8)         , intent(in) :: area_i(:)
    real(r8)         , intent(in) :: area_o(:)
    integer          , intent(in) :: mask_i(:)
    real(r8)         , intent(in) :: frac_o(:)
    real(r8)         , intent(in) :: data_i(:)    ! data on input grid
    real(r8)         , intent(in) :: data_o(:)    ! data on output grid
    character(len=*) , intent(in) :: name         ! name of field
    character(len=*) , intent(in) :: units        ! units of field
    integer          , intent(in) :: ndiag

    ! !LOCAL VARIABLES:
    real(r8) :: loc_gdata_i     ! local sum of input data
    real(r8) :: loc_gdata_o     ! local sum of output data
    real(r8) :: gdata_i         ! global sum of input data
    real(r8) :: gdata_o         ! global sum of output data
    real(r8) :: loc_gwt_i       ! local global sum of input weights (area * frac)
    real(r8) :: loc_gwt_o       ! local global sum of output weights (area * frac)
    real(r8) :: gwt_i           ! global sum of input weights (area * frac)
    real(r8) :: gwt_o           ! global sum of output weights (area * frac)
    integer  :: ns_i, ns_o      ! sizes of input & output grids
    integer  :: ni,no,k         ! indices
    integer  :: ier             ! error code
    character(len=*), parameter :: subname = "output_diagnostics_continuous"
    !------------------------------------------------------------------------------

    ! Error check for array size consistencies
    if (size(data_i) /= ns_i .or. size(data_o) /= ns_o) then
       write(6,*) subname//' ERROR: array size inconsistencies for ', trim(name)
       write(6,*) 'size(data_i) = ', size(data_i)
       write(6,*) 'ns_i         = ', ns_i
       write(6,*) 'size(data_o) = ', size(data_o)
       write(6,*) 'ns_o         = ', ns_o
       call shr_sys_abort()
    end if
    if (size(frac_o) /= ns_o) then
       write(6,*) subname//' ERROR: incorrect size of frac_o'
       write(6,*) 'size(frac_o) = ', size(frac_o)
       write(6,*) 'ns_o = ', ns_o
       call shr_sys_abort()
    end if
    if (size(mask_i) /= ns_i) then
       write(6,*) subname//' ERROR: incorrect size of mask_i'
       write(6,*) 'size(mask_i) = ', size(mask_i)
       write(6,*) 'ns_i = ', ns_i
       call shr_sys_abort()
    end if

    ! Sums on input grid
    loc_gdata_i = 0.
    loc_gwt_i = 0.
    do ni = 1,ns_i
       loc_gdata_i = loc_gdata_i + data_i(ni) * area_i(ni) * mask_i(ni)
       loc_gwt_i = loc_gwt_i + area_i(ni) * mask_i(ni)
    end do
    call mpi_reduce(loc_gdata_i,gdata_i,1,MPI_REAL8,MPI_SUM,0,mpicom,ier)

    ! Sums on output grid
    loc_gdata_o = 0.
    loc_gwt_o = 0.
    do no = 1,ns_o
       loc_gdata_o = loc_gdata_o + data_o(no) * area_o(no) * frac_o(no)
       loc_gwt_o = loc_gwt_o + area_o(no) * frac_o(no)
    end do
    call mpi_reduce(loc_gdata_o,gdata_o,1,MPI_REAL8,MPI_SUM,0,mpicom,ier)

    ! Correct units
    gdata_i = gdata_i / gwt_i
    gdata_o = gdata_o / gwt_o

    ! Diagnostic output
    if (root_task) then
       write (ndiag,*)
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,2001)
       write (ndiag,2002) units, units
2001   format (1x,'   parameter              input grid          output grid')
2002   format (1x,'                 ', a24, a24)
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,*)
       write (ndiag,2003) name,          gdata_i, gdata_o
2003   format (1x,a12,           f22.3,f17.3)
    end if

  end subroutine output_diagnostics_continuous

  !===============================================================
  subroutine output_diagnostics_continuous_outonly(area_o, frac_o, data_o, name, units, ndiag)
    !
    ! Output diagnostics for a continuous field, just on the output grid
    ! This is used when the average of the field on the input grid is not of interest (e.g.,
    ! when the output quantity is the standard deviation of the input field)
    !
    use mkvarpar, only : re

    ! input/output variables
    real(r8)          , intent(in) :: area_o(:) 
    real(r8)          , intent(in) :: frac_o(:)
    real(r8)          , intent(in) :: data_o(:)    ! data on output grid
    character(len=*)  , intent(in) :: name         ! name of field
    character(len=*)  , intent(in) :: units        ! units of field
    integer           , intent(in) :: ndiag        ! unit number for diagnostic output

    ! local variables:
    real(r8) :: gdata_o         ! global sum of output data
    real(r8) :: gwt_o           ! global sum of output weights (area * frac)
    integer  :: ns_o            ! size of output grid
    integer  :: no,k            ! indices
    character(len=*), parameter :: subname = "output_diagnostics_continuous_outonly"
    !------------------------------------------------------------------------------

    ns_o = size(area_o)

    ! Error check for array size consistencies
    if (size(data_o) /= ns_o) then
       write(6,*) subname//' ERROR: array size inconsistencies for ', trim(name)
       write(6,*) 'size(data_o) = ', size(data_o)
       write(6,*) 'ns_o         = ', ns_o
       call shr_sys_abort()
    end if

    ! Sums on output grid
    gdata_o = 0.
    gwt_o = 0.
    do no = 1,ns_o
       gdata_o = gdata_o + data_o(no)*area_o(no)*frac_o(no)
       gwt_o = gwt_o + area_o(no)*frac_o(no)
    end do

    ! Correct units
    gdata_o = gdata_o / gwt_o

    ! Diagnostic output
    if (root_task) then
       write (ndiag,*)
       write (ndiag,*)
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,2001)
       write (ndiag,2002) units
2001   format (1x,'   parameter              output grid')
2002   format (1x,'                 ', a24)
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,*)
       write (ndiag,2003) name,          gdata_o
2003   format (1x,a12,           f22.3)
    end if
       
  end subroutine output_diagnostics_continuous_outonly

  !===============================================================
  subroutine output_diagnostics_index(area_i, area_o, mask_i, frac_o, data_i, data_o, name, minval, maxval, ndiag)
    !
    ! Output diagnostics for an index field: area of each index in input and output
    !
    use mkvarpar, only : re

    ! input/output variables
    real(r8)           , intent(in) :: area_i(:) 
    real(r8)           , intent(in) :: area_o(:) 
    integer            , intent(in) :: mask_i(:)
    real(r8)           , intent(in) :: frac_o(:)
    integer            , intent(in) :: data_i(:) ! data on input grid
    integer            , intent(in) :: data_o(:) ! data on output grid
    character(len=*)   , intent(in) :: name      ! name of field
    integer            , intent(in) :: minval    ! minimum valid value
    integer            , intent(in) :: maxval    ! minimum valid value
    integer            , intent(in) :: ndiag     ! unit number for diagnostic output

    ! local variables:
    integer               :: ns_i, ns_o     ! sizes of input & output grids
    integer               :: ni, no, k      ! indices
    real(r8), allocatable :: loc_garea_i(:) ! input grid: global area of each index
    real(r8), allocatable :: loc_garea_o(:) ! output grid: global area of each index
    real(r8), allocatable :: garea_i(:)     ! input grid: global area of each index
    real(r8), allocatable :: garea_o(:)     ! output grid: global area of each index
    integer               :: ier            ! error status
    character(len=*), parameter :: subname = 'output_diagnostics_index'
    !-----------------------------------------------------------------------

    ! Error check for array size consistencies

    ns_i = size(area_i)
    ns_o = size(area_o)

    if (size(data_i) /= ns_i .or. size(data_o) /= ns_o) then
       write(6,*) subname//' ERROR: array size inconsistencies for ', trim(name)
       write(6,*) 'size(data_i) = ', size(data_i)
       write(6,*) 'ns_i         = ', ns_i
       write(6,*) 'size(data_o) = ', size(data_o)
       write(6,*) 'ns_o         = ', ns_o
       call shr_sys_abort()
    end if
    if (size(frac_o) /= ns_o) then
       write(6,*) subname//' ERROR: incorrect size of frac_o'
       write(6,*) 'size(frac_o) = ', size(frac_o)
       write(6,*) 'ns_o = ', ns_o
       call shr_sys_abort()
    end if
    if (size(mask_i) /= ns_i) then
       write(6,*) subname//' ERROR: incorrect size of mask_i'
       write(6,*) 'size(mask_i) = ', size(mask_i)
       write(6,*) 'ns_i = ', ns_i
       call shr_sys_abort()
    end if

    ! Sum areas on input grid
    allocate(loc_garea_i(minval:maxval), stat=ier)
    allocate(garea_i(minval:maxval), stat=ier)
    if (ier/=0) call shr_sys_abort()
    loc_garea_i(:) = 0.
    do ni = 1, ns_i
       k = data_i(ni)
       if (k >= minval .and. k <= maxval) then
          loc_garea_i(k) = loc_garea_i(k) + area_i(ni) * mask_i(ni) * re**2
       end if
    end do
    call mpi_reduce(loc_garea_i, garea_i, size(garea_i), size(garea_i), MPI_REAL8, MPI_SUM, 0, mpicom, ier)

    ! Sum areas on output grid
    allocate(loc_garea_o(minval:maxval), stat=ier)
    allocate(garea_o(minval:maxval), stat=ier)
    if (ier/=0) call shr_sys_abort()
    loc_garea_o(:) = 0.
    do no = 1, ns_o
       k = data_o(no)
       if (k >= minval .and. k <= maxval) then
          loc_garea_o(k) = garea_o(k) + area_o(no) * frac_o(no) * re**2
       end if
    end do
    call mpi_reduce(loc_garea_o, garea_o, size(garea_o), size(garea_o), MPI_REAL8, MPI_SUM, 0, mpicom, ier)

    ! Write results
    if (root_task) then
       write (ndiag,*)
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,2001)
2001   format (1x,'index      input grid area  output grid area',/ &
               1x,'               10**6 km**2       10**6 km**2')
       write (ndiag,'(1x,70a1)') ('.',k=1,70)
       write (ndiag,*)
       do k = minval, maxval
          write (ndiag,2002) k, garea_i(k)*1.e-06, garea_o(k)*1.e-06
2002      format (1x,i9,f17.3,f18.3)
       end do
    end if

  end subroutine output_diagnostics_index

end module mkdiagnosticsMod
