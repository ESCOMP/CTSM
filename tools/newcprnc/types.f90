module types
!
! Derived type containers for field stats (varinfo) and difference stats (diffinfo)
!
   use prec, only: r8

   implicit none

PUBLIC

   include 'netcdf.inc'

   type varinfo
      character(len=NF_MAX_NAME) :: name          ! variable name
      character(len=NF_MAX_NAME) :: dimnames(3)   ! spatial dimension names
 
      integer :: npossible     ! number of field values
      integer :: nfill         ! number of fillvalues
      integer :: imax          ! 1st dimension index where field max occurs
      integer :: jmax          ! 2nd dimension index where field max occurs
      integer :: kmax          ! 3rd dimension index where field max occurs
      integer :: imin          ! 1st dimension index where field min occurs
      integer :: jmin          ! 2nd dimension index where field min occurs
      integer :: kmin          ! 3rd dimension index where field min occurs

      real(r8) :: avg          ! field average value
      real(r8) :: arrmax       ! field max
      real(r8) :: arrmin       ! field min
   end type varinfo

   type diffinfo
      character(len=NF_MAX_NAME) :: name          ! variable name
      character(len=NF_MAX_NAME) :: dimnames(3)   ! spatial dimension names

      real(r8) :: dmax         ! max diff
      real(r8) :: rdmax        ! max relative diff
      real(r8) :: dvals(2)     ! field values where max diff occurs
      real(r8) :: rdvals(2)    ! field values where max relative diff occurs
      real(r8) :: rms          ! rms difference
      real(r8) :: mwrms        ! mass-weighted rms diff (unimplemented)
      real(r8) :: rdbar        ! mean relative difference
      real(r8) :: rdlnbar      ! used in computation of number of matching digits

      integer :: npossible     ! number of non-fillvalue locations
      integer :: ndiff         ! number of diffs
      integer :: imax          ! 1st dimension index where biggest diff occurs
      integer :: jmax          ! 2nd dimension index where biggest diff occurs
      integer :: kmax          ! 3rd dimension index where biggest diff occurs
      integer :: irmax         ! 1st dimension index where biggest relative diff occurs
      integer :: jrmax         ! 2nd dimension index where biggest relative diff occurs
      integer :: krmax         ! 3rd dimension index where biggest relative diff occurs
   end type diffinfo
end module types
