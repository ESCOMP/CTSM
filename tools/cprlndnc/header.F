module header
  use precision
!
! Netcdf id
!
  integer, dimension(2) :: ncid = -1
!
! Dimensions
!
  integer, dimension(2) :: londimid = -1
  integer, dimension(2) :: latdimid = -1
  integer, dimension(2) :: levdimid = -1
  integer, dimension(2) :: unlimdimid = -1
!
! Dimension variables
!
  integer, dimension(2) :: lonid = -1
  integer, dimension(2) :: latid = -1
  integer, dimension(2) :: levid = -1
  integer, dimension(2) :: timeid = -1

  real(r8), allocatable :: lon(:,:), lat(:,:), lev(:,:)
  real(r8), allocatable :: time(:), time2(:)
!
! Character variables
!
  character*16, dimension(2) :: case
  character*80, dimension(2) :: title

  integer, dimension(2) :: date_writtenid
  integer, dimension(2) :: time_writtenid

  character*8, allocatable :: date_written(:,:)
  character*8, allocatable :: time_written(:,:)
!
! Integer variables
!
  integer, dimension(2) :: ntime

  integer, dimension(2) :: nstephid = -1
  integer, dimension(2) :: ncdateid = -1
  integer, dimension(2) :: ncsecid = -1

  integer, allocatable :: nsteph(:,:)
  integer, allocatable :: ncdate(:,:)
  integer, allocatable :: ncsec(:,:)
!
! Real variables
!
  integer, dimension(2) :: areaid = -1

  real(r8), allocatable :: area(:,:,:)

end module header
