module mkutilsMod


  ! General-purpose utilities
  use ESMF
  use shr_kind_mod, only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod , only : shr_sys_abort

  implicit none
  private

  ! PUBLIC MEMBER FUNCTIONS:
  public :: normalize_classes_by_gcell  ! renormalize array so values are given as % of total grid cell area
  public :: slightly_below
  public :: slightly_above
  public :: get_filename  !Returns filename given full pathname
  public :: chkerr

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine normalize_classes_by_gcell(classes_pct_tot, sums, classes_pct_gcell)
    !
    ! Renormalizes an array (gcell x class) so that values are given as % of total grid cell area
    !
    ! Specifically: Given (1) an array specifying the % cover of different classes, as a % of
    ! some total ('classes_pct_tot'), and (2) a vector giving these totals ('sums'), expressed
    ! as % of grid cell area: Returns an array ('classes_pct_gcell') of the same
    ! dimensionality as classes_pct_tot, where the values now give the % cover of each class
    ! as a % of total grid cell area.
    !
    ! The size of 'sums' should match the size of the first dimension in 'classes_pct_tot' and
    ! 'classes_pct_gcell'
    !
    ! For example, if classes_pct_tot(n,i) gives the % of the urban area in grid cell n that is
    ! in urban class #i, and sums(n) gives the % of grid cell n that is urban, then
    ! classes_pct_gcell(n,i) will give the % of the total area of grid cell n that is in urban
    ! class #i.
    !
    ! input/output variables
    real(r8), intent(in) :: classes_pct_tot(:,:)   ! % cover of classes as % of total
    real(r8), intent(in) :: sums(:)                ! totals, as % of grid cell
    real(r8), intent(out):: classes_pct_gcell(:,:) ! % cover of classes as % of grid cell
    !
    ! local variables
    integer :: n, n_max
    character(len=*), parameter :: subname = "normalize_classes_by_gcell"
    !------------------------------------------------------------------------------

    ! Error-check inputs

    n_max = size(sums)
    if (size(classes_pct_tot, 1)   /= n_max .or. &
        size(classes_pct_gcell, 1) /= n_max) then
       write(6,*) subname//' ERROR: array size mismatch'
       write(6,*) 'size(sums)                 = ', n_max
       write(6,*) 'size(classes_pct_tot, 1)   = ', size(classes_pct_tot, 1)
       write(6,*) 'size(classes_pct_gcell, 1) = ', size(classes_pct_gcell, 1)
       call shr_sys_abort()
    end if

    if (size(classes_pct_tot, 2) /= size(classes_pct_gcell, 2)) then
       write(6,*) subname//' ERROR: array size mismatch'
       write(6,*) 'size(classes_pct_tot, 2)   = ', size(classes_pct_tot, 2)
       write(6,*) 'size(classes_pct_gcell, 2) = ', size(classes_pct_gcell, 2)
       call shr_sys_abort()
    end if

    ! Do the work

    do n = 1, n_max
       classes_pct_gcell(n,:) = classes_pct_tot(n,:) * (sums(n)/100._r4)
    end do
  end subroutine normalize_classes_by_gcell

  !===============================================================
  logical function slightly_below(a, b, eps)

    ! Returns true if a is slightly below b; false if a is significantly below b or if a is
    ! greater than or equal to b
    ! if provided, eps gives the relative error allowed for checking the "slightly"
    ! condition; if not provided, the tolerance defaults to the value given by eps_default

    ! !ARGUMENTS:
    real(r8), intent(in) :: a
    real(r8), intent(in) :: b
    real(r8), intent(in), optional :: eps

    ! !LOCAL VARIABLES:
    real(r8) :: l_eps
    real(r8), parameter :: eps_default = 1.e-15_r8  ! default relative error tolerance
    !------------------------------------------------------------------------------

    if (present(eps)) then
       l_eps = eps
    else
       l_eps = eps_default
    end if

    if (a < b .and. (b - a)/b < l_eps) then
       slightly_below = .true.
    else
       slightly_below = .false.
    end if

  end function slightly_below

  !===============================================================
  logical function slightly_above(a, b, eps)

    ! Returns true if a is slightly above b; false if a is significantly above b or if a is
    ! less than or equal to b
    !
    ! if provided, eps gives the relative error allowed for checking the "slightly"
    ! condition; if not provided, the tolerance defaults to the value given by eps_default

    ! input/output variables
    real(r8), intent(in) :: a
    real(r8), intent(in) :: b
    real(r8), intent(in), optional :: eps

    ! local variables:
    real(r8) :: l_eps
    real(r8), parameter :: eps_default = 1.e-15_r8  ! default relative error tolerance
    !------------------------------------------------------------------------------

    if (present(eps)) then
       l_eps = eps
    else
       l_eps = eps_default
    end if

    if (a > b .and. (a - b)/b < l_eps) then
       slightly_above = .true.
    else
       slightly_above = .false.
    end if

  end function slightly_above

  !===============================================================
  logical function chkerr(rc, line, file)
    integer          , intent(in) :: rc
    integer          , intent(in) :: line
    character(len=*) , intent(in) :: file

    ! local variables
    integer :: lrc
    chkerr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       chkerr = .true.
    endif
  end function chkerr

  !===============================================================
  character(len=256) function get_filename (fulpath)
    ! Returns filename given full pathname

    ! input/output variables
    character(len=*), intent(in)  :: fulpath !full pathname

    ! local variables:
    integer :: i    !loop index
    integer :: klen !length of fulpath character string
    !------------------------------------------------------------------------

    klen = len_trim(fulpath)
    do i = klen, 1, -1
       if (fulpath(i:i) == '/') go to 10
    end do
    i = 0
10  get_filename = fulpath(i+1:klen)

  end function get_filename


end module mkutilsMod
