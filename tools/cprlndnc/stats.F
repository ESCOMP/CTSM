module stats
  use precision

  integer :: npos      ! # of possible differences
  integer :: ndif      ! # of differences

  real(r8) :: difmx    ! max diff
  real(r8) :: rdifmx   ! maximum relative difference
  real(r8) :: xmx(2)    ! maximum field value
  real(r8) :: xmn(2)    ! minimum field value
  real(r8) :: rms      ! root mean square difference
  real(r8) :: xbar(2)   ! mean field value
  real(r8) :: rdbar    ! mean relative difference
  real(r8) :: rdlnbar  ! mean alog(difference)
  real(r8) :: diffmw   ! mass weighted diff

  real(r8) :: dmxsv(2)  ! values corresponding to difmx
  real(r8) :: rdmxsv(2) ! values corresponding to rdifmx

  integer :: imx(2)     ! i-index of field maximum
  integer :: jmx(2)     ! j-index of field maximum
  integer :: kmx(2)     ! k-index of field maximum
  integer :: imn(2)     ! i-index of field minimum
  integer :: jmn(2)     ! j-index of field minimum
  integer :: kmn(2)     ! k-index of field minimum
  integer :: idmxsv    ! i-index of max difference
  integer :: jdmxsv    ! j-index of max difference
  integer :: kdmxsv    ! k-index of max difference
  integer :: irdmxsv   ! i-index of max relative difference
  integer :: jrdmxsv   ! j-index of max relative difference
  integer :: krdmxsv   ! k-index of max relative difference

end module stats
