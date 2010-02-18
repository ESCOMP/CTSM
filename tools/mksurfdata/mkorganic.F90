!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkorganic
!
! !INTERFACE:
subroutine mkorganic(lsmlon, lsmlat, fname, ndiag, organic_o)
!
! !DESCRIPTION:
! make organic matter dataset
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use fileutils   , only : getfil
  use domainMod   , only : domain_type,domain_clean,domain_setptrs
  use creategridMod, only : read_domain
  use mkvarpar	
  use mkvarsur    , only : ldomain
  use mkvarctl    
  use areaMod     , only : areaini,areaave,gridmap_type,gridmap_clean  
  use ncdio
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: lsmlon, lsmlat          ! clm grid resolution
  character(len=*), intent(in) :: fname           ! input dataset file name
  integer , intent(in) :: ndiag                   ! unit number for diag out
  real(r8), intent(out):: organic_o(lsmlon,lsmlat,nlevsoi)    ! output grid:
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! 
! Author: David Lawrence
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: nlon_i                          ! input grid : lon points
  integer  :: nlat_i                          ! input grid : lat points

  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: organic_i(:,:,:)  ! input grid: total column organic matter
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  real(r8), allocatable :: omlev_i(:,:)       ! input grid: organic matter on lev
  real(r8), allocatable :: omlev_o(:,:)       ! output grid: organic matter on lev
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: gomlev_i                         ! input  grid: global organic on lev
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: gomlev_o                         ! output grid: global organic on lev
  real(r8) :: garea_o                         ! output grid: global area

  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,n,m                           ! indices
  integer  :: lev                             ! level index
  integer   :: nlay                            ! number of soil layers
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=256) locfn                    ! local dataset file name
  character(len=32) :: subname = 'mkorganic'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make organic matter dataset .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call getfil (fname, locfn, 0)

  call read_domain(tdomain,locfn)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(locfn, 0, ncid), subname)

  call check_ret(nf_inq_dimid  (ncid, 'number_of_layers', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, nlay), subname)

  allocate(organic_i(nlon_i,nlat_i,nlay), omlev_i(nlon_i,nlat_i), &
           omlev_o(lsmlon,lsmlat),stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'ORGANIC', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, organic_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Compute local fields _o

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat))

  mask_i = 1.0_r8
  mask_o = 1.0_r8

  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  mask_i = float(tdomain%mask(:,:))
  call areaave(mask_i,mask_o,tgridmap)

  call gridmap_clean(tgridmap)
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  ldomain%frac = mask_o

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  call areaave(organic_i,organic_o,tgridmap)

  do lev = 1,nlevsoi

     ! Check for conservation

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        if ((organic_o(io,jo,lev)) > 130.000001_r8) then
           write (6,*) 'MKORGANIC error: organic = ',organic_o(io,jo,lev), &
                ' greater than 130.000001 for column, row = ',io,jo
           call shr_sys_flush(6)
           stop
        end if
     enddo
     enddo

     ! Diagnostic output

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)
     write (ndiag,2002) gomlev_i*1.e-06,gomlev_o*1.e-06
     write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'organic    ',f14.3,f17.3)
2004 format (1x,'all surface ',f14.3,f17.3)

     if (lsmlat > 1) then
        k = lsmlat/2
        write (ndiag,*)
        write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
        write (ndiag,'(f10.3,a14)')ldomain%area(1,k)*1.e-06,' x 10**6 km**2'
        write (ndiag,*)
     endif
     call shr_sys_flush(ndiag)

     write (6,*) 'Successfully made organic matter, level = ', lev
     write (6,*)
     call shr_sys_flush(6)

  end do   ! lev

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (organic_i, omlev_i, omlev_o)
  deallocate (mask_i,mask_o)

end subroutine mkorganic

