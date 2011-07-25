module mkglacierMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkglacierMod
!
! !DESCRIPTION:
! Make percent glacier.
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!-----------------------------------------------------------------------

! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush

  implicit none

  private

! !PUBLIC MEMBER FUNCTIONS:
  public ::  mkglacier   ! read in the glacier dataset and convert to percent-glacier

!EOP
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkglacier
!
! !INTERFACE:
subroutine mkglacier(lsmlon, lsmlat, fname, ndiag, zero_out, glac_o )
!
! !DESCRIPTION:
! make percent glacier
!
! !USES:
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
  logical , intent(in) :: zero_out                ! if should zero glacier out
  integer , intent(in) :: ndiag                   ! unit number for diag out
  real(r8), intent(out):: glac_o(lsmlon,lsmlat)   ! output grid: %glacier
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

  real(r8), allocatable :: glac_i(:,:)        ! input grid: percent glac
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  real(r8), allocatable :: fld_i(:,:)         ! input grid: dummy field
  real(r8), allocatable :: fld_o(:,:)         ! output grid: dummy field
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: gglac_i                          ! input  grid: global glac
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: gglac_o                          ! output grid: global glac
  real(r8) :: garea_o                         ! output grid: global area

  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,n,m                           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=256) locfn                    ! local dataset file name
  character(len=32) :: subname = 'mkglacier'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %glacier .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call read_domain(tdomain,fname)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(locfn, 0, ncid), subname)

  allocate(glac_i(nlon_i,nlat_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'PCT_GLACIER', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, glac_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Compute local fields _o

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat))
  allocate( fld_i(nlon_i,nlat_i), fld_o(lsmlon,lsmlat))

  mask_i = 1.0_r8
  mask_o = 1.0_r8
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  if ( zero_out )then

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        glac_o(io,jo) = 0.
     enddo
     enddo

  else

     mask_i = float(tdomain%mask(:,:))
     call areaave(mask_i,mask_o,tgridmap)

     call gridmap_clean(tgridmap)
     call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

     ! Area-average percent cover on input grid to output grid 
     ! and correct according to land landmask
     ! Note that percent cover is in terms of total grid area.

     call areaave(glac_i,glac_o,tgridmap)

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        if (glac_o(io,jo) < 1.) glac_o(io,jo) = 0.
     enddo
     enddo

  end if

  ! Check for conservation

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     if ((glac_o(io,jo)) > 100.000001_r8) then
        write (6,*) 'MKGLACIER error: glacier = ',glac_o(io,jo), &
                ' greater than 100.000001 for column, row = ',io,jo
        call shr_sys_flush(6)
        stop
     end if
  enddo
  enddo

  ! Some error checking and writing of global values before and after the regrid

  if ( .not. zero_out )then

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
           write (6,*) 'MKGLACIER error: input field not conserved'
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

     gglac_i = 0.
     garea_i = 0.

     do ji = 1, nlat_i
     do ii = 1, nlon_i
        garea_i = garea_i + tdomain%area(ii,ji)
        gglac_i = gglac_i + glac_i(ii,ji)*(tdomain%area(ii,ji)/100.) * &
                            tdomain%frac(ii,ji)
     end do
     end do

     ! Output grid

     gglac_o = 0.
     garea_o = 0.

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        garea_o = garea_o + ldomain%area(io,jo)
        gglac_o = gglac_o + glac_o(io,jo)*(ldomain%area(io,jo)/100.) * &
                            ldomain%frac(io,jo)
     end do
     end do

     ! Diagnostic output

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'Glacier Output'
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)
     write (ndiag,2002) gglac_i*1.e-06,gglac_o*1.e-06
     write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'glaciers    ',f14.3,f17.3)
2004 format (1x,'all surface ',f14.3,f17.3)

     if (lsmlat > 1) then
        k = lsmlat/2
        write (ndiag,*)
        write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
        write (ndiag,'(f10.3,a14)')ldomain%area(1,k)*1.e-06,' x 10**6 km**2'
        write (ndiag,*)
     endif
     call shr_sys_flush(ndiag)

  end if

  write (6,*) 'Successfully made %glacier'
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (glac_i)
  deallocate (mask_i,mask_o)
  deallocate ( fld_i, fld_o)

end subroutine mkglacier

!-----------------------------------------------------------------------

end module mkglacierMod
