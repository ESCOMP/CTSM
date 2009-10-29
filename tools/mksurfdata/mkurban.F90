!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurban
!
! !INTERFACE:
subroutine mkurban(lsmlon, lsmlat, fname, ndiag, urbn_o)
!
! !DESCRIPTION:
! make percent urban
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
  real(r8), intent(out):: urbn_o(lsmlon,lsmlat)    ! output grid: %urban
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: nlon_i                          ! input grid : lon points
  integer  :: nlat_i                          ! input grid : lat points

  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: urbn_i(:,:)        ! input grid: percent urbn
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  real(r8), allocatable :: fld_i(:,:)         ! input grid: dummy field
  real(r8), allocatable :: fld_o(:,:)         ! output grid: dummy field
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: gurbn_i                          ! input  grid: global urbn
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: gurbn_o                          ! output grid: global urbn
  real(r8) :: garea_o                         ! output grid: global area

  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,n,m                           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001_r8             ! max error: sum overlap wts ne 1
  character(len=256) locfn                    ! local dataset file name
  character(len=32) :: subname = 'mkurban'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %urban .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call getfil (fname, locfn, 0)

  call read_domain(tdomain,locfn)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(locfn, 0, ncid), subname)

  allocate(urbn_i(nlon_i,nlat_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'PCT_URBAN', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, urbn_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Compute local fields _o

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat))
  allocate( fld_i(nlon_i,nlat_i), fld_o(lsmlon,lsmlat))

  mask_i = 1.0_r8
  mask_o = 1.0_r8
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  mask_i = float(tdomain%mask(:,:))
  call areaave(mask_i,mask_o,tgridmap)

  call gridmap_clean(tgridmap)
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  call areaave(urbn_i,urbn_o,tgridmap)

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
        if (urbn_o(io,jo) < 0.1_r8) urbn_o(io,jo) = 0._r8
        if (all_urban )             urbn_o(io,jo) = 100.00_r8
  enddo
  enddo

  ! Check for conservation

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     if ((urbn_o(io,jo)) > 100.000001_r8) then
        write (6,*) 'MKURBAN error: urban = ',urbn_o(io,jo), &
                ' greater than 100.000001 for column, row = ',io,jo
        call shr_sys_flush(6)
        stop
     end if
  enddo
  enddo

  ! Global sum of output field -- must multiply by fraction of
  ! output grid that is land as determined by input grid

  sum_fldi = 0.0_r8
  do ji = 1, nlat_i
  do ii = 1, nlon_i
    fld_i(ii,ji) = ((ji-1)*nlon_i + ii) * tdomain%mask(ii,ji)
    sum_fldi = sum_fldi + tdomain%area(ii,ji) * fld_i(ii,ji)
  enddo
  enddo

  call areaave(fld_i,fld_o,tgridmap)

  sum_fldo = 0._r8
  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     fld_o(io,jo) = fld_o(io,jo)*mask_o(io,jo)
     sum_fldo = sum_fldo + ldomain%area(io,jo) * fld_o(io,jo)
  end do
  end do

  ! -----------------------------------------------------------------
  ! Error check1
  ! Compare global sum fld_o to global sum fld_i.
  ! -----------------------------------------------------------------

  if ( .not. all_urban .and. trim(mksrf_gridtype) == 'global') then
     if ( abs(sum_fldo/sum_fldi-1._r8) > relerr ) then
        write (6,*) 'MKURBAN error: input field not conserved'
        write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
        write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
        stop
     end if
  end if

  ! -----------------------------------------------------------------
  ! Error check2
  ! Compare global areas on input and output grids
  ! -----------------------------------------------------------------

  ! Input grid

  gurbn_i = 0._r8
  garea_i = 0._r8

  do ji = 1, nlat_i
  do ii = 1, nlon_i
     garea_i = garea_i + tdomain%area(ii,ji)
     gurbn_i = gurbn_i + urbn_i(ii,ji)*(tdomain%area(ii,ji)/100._r8) * &
                         tdomain%frac(ii,ji)
  end do
  end do

  ! Output grid

  gurbn_o = 0._r8
  garea_o = 0._r8

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     garea_o = garea_o + ldomain%area(io,jo)
     gurbn_o = gurbn_o + urbn_o(io,jo)*(ldomain%area(io,jo)/100._r8) * &
                         ldomain%frac(io,jo)
  end do
  end do

  ! Diagnostic output

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'Urban Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)
! write (ndiag,2002) gurbn_i*1.e-06,gurbn_o*1.e-06
  write (ndiag,2003) gurbn_i*1.e-06,gurbn_o*1.e-06
  write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'urban       ',f14.3,f17.3)
2003 format (1x,'urban       ',f14.3,f22.8)
2004 format (1x,'all surface ',f14.3,f17.3)

  if (lsmlat > 1) then
     k = lsmlat/2
     write (ndiag,*)
     write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
     write (ndiag,'(f10.3,a14)')ldomain%area(1,k)*1.e-06,' x 10**6 km**2'
     write (ndiag,*)
  endif
  call shr_sys_flush(ndiag)

  write (6,*) 'Successfully made %urban'
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (urbn_i)
  deallocate (mask_i,mask_o)
  deallocate ( fld_i, fld_o)

end subroutine mkurban

