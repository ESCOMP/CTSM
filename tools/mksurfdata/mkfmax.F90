!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkfmax
!
! !INTERFACE:
subroutine mkfmax(lsmlon, lsmlat, fname, ndiag, fmax_o)
!
! !DESCRIPTION:
! make percent fmax
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
  real(r8), intent(out):: fmax_o(lsmlon,lsmlat)    ! output grid: %fmax
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Revised: Nan Rosenbloom - used mkglacier.F90 as template.
! Original Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: nlon_i                          ! input grid : lon points
  integer  :: nlat_i                          ! input grid : lat points

  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: fmax_i(:,:)        ! input grid: percent fmax
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  real(r8), allocatable :: fld_i(:,:)         ! input grid: dummy field
  real(r8), allocatable :: fld_o(:,:)         ! output grid: dummy field
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: gfmax_i                         ! input  grid: global fmax
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: gfmax_o                         ! output grid: global fmax
  real(r8) :: garea_o                         ! output grid: global area

  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,n,m                           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=256) locfn                    ! local dataset file name
  character(len=32) :: subname = 'mkfmax'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %fmax .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call getfil (fname, locfn, 0)

  call read_domain(tdomain,locfn)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(locfn, 0, ncid), subname)

  allocate(fmax_i(nlon_i,nlat_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'FMAX', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, fmax_i), subname)

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
  ldomain%frac = mask_o

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  call areaave(fmax_i,fmax_o,tgridmap)

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
        if (fmax_o(io,jo) == 0.0) then
		fmax_o(io,jo) = .365783
		! fmax_o(io,jo) = globalAvg
	end if
        if (fmax_o(io,jo) == -999.99) then
		fmax_o(io,jo) = .365783
		! fmax_o(io,jo) = globalAvg
	end if
  enddo
  enddo

  ! Check for conservation

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     if ((fmax_o(io,jo)) > 1.000001_r8) then
        write (6,*) 'MKFMAX error: fmax = ',fmax_o(io,jo), &
                ' greater than 1.000001 for column, row = ',io,jo
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

  sum_fldo = 0.
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

  if ( trim(mksrf_gridtype) == 'global') then
     if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
        write (6,*) 'MKFMAX error: input field not conserved'
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

  gfmax_i = 0.
  garea_i = 0.

  do ji = 1, nlat_i
  do ii = 1, nlon_i
     garea_i = garea_i + tdomain%area(ii,ji)
     gfmax_i = gfmax_i + fmax_i(ii,ji)*(tdomain%area(ii,ji)/100.) * &
                                        tdomain%frac(ii,ji)
  end do
  end do

  ! Output grid

  gfmax_o = 0.
  garea_o = 0.

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     garea_o = garea_o + ldomain%area(io,jo)
     gfmax_o = gfmax_o + fmax_o(io,jo)*(ldomain%area(io,jo)/100.) * &
                                        ldomain%frac(io,jo)
     if ( (ldomain%frac(io,jo) < 0.0) .or. (ldomain%frac(io,jo) > 1.0001) )then
        write(6,*) "ERROR:: frac out of range: ", ldomain%frac(io,jo), io, jo
        stop
     end if
  end do
  end do

  ! Diagnostic output

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'Maximum Fractional Saturated Area Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)
  write (ndiag,2002) gfmax_i*1.e-06,gfmax_o*1.e-06
  write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'fmax    ',f14.3,f17.3)
2004 format (1x,'all surface ',f14.3,f17.3)

  if (lsmlat > 1) then
     k = lsmlat/2
     write (ndiag,*)
     write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
     write (ndiag,'(f10.3,a14)')ldomain%area(1,k)*1.e-06,' x 10**6 km**2'
     write (ndiag,*)
  endif
  call shr_sys_flush(ndiag)

  write (6,*) 'Successfully made %fmax'
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (fmax_i)
  deallocate (mask_i,mask_o)
  deallocate ( fld_i, fld_o)

end subroutine mkfmax

