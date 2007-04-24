module mkpftMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkpft
!
! !DESCRIPTION:
! Make PFT data
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkpft
!
! !INTERFACE:
subroutine mkpft(lsmlon, lsmlat, fpft, ndiag, pctlnd_o, pctpft_o, pct_pft_i)
!
! !DESCRIPTION:
! Make PFT data
! This dataset consists of the %cover of the [numpft]+1 PFTs used by
! the model. The input %cover pertains to the "vegetated" portion of the
! grid cell and sums to 100. The real portion of each grid cell
! covered by each PFT is the PFT cover times the fraction of the
! grid cell that is land. This is the quantity preserved when
! area-averaging from the input (1/2 degree) grid to the models grid.
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
  character(len=*), intent(in) :: fpft            ! input pft dataset file name
  integer , intent(in) :: ndiag                   ! unit number for diag out
  real(r8), intent(out):: pctlnd_o(lsmlon,lsmlat) ! output grid:%land/gridcell
  real(r8), intent(out):: pctpft_o(lsmlon,lsmlat,0:numpft)  ! PFT cover 
                                                  ! (% of vegetated area)
  real(r8), optional, pointer :: pct_pft_i(:,:,:) ! Plant function type on input grid
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer  :: nlon_i                          ! input grid : lon points
  integer  :: nlat_i                          ! input grid : lat points
  integer  :: numpft_i                        ! num of plant types input data

  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: pctpft_i(:,:,:)    ! input grid: PFT percent
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  real(r8), allocatable :: fld_i(:,:)         ! input grid: dummy field
  real(r8), allocatable :: fld_o(:,:)         ! output grid: dummy field
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: wst(0:numpft)                   ! as pft_o at specific io, jo
  real(r8) :: wst_sum                         ! sum of %pft
  real(r8) :: gpft_o(0:numpft)                ! output grid: global area pfts
  real(r8) :: garea_o                         ! output grid: global area
  real(r8) :: gpft_i(0:numpft)                ! input grid: global area pfts
  real(r8) :: garea_i                         ! input grid: global area

  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,n,m                           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=256) locfn                    ! local dataset file name

  character(len=35)  veg(0:numpft)            ! vegetation types
  data veg( 0) /'not vegetated'                      /
  data veg( 1) /'needleleaf evergreen temperate tree'/
  data veg( 2) /'needleleaf evergreen boreal tree'   /
  data veg( 3) /'needleleaf deciduous boreal tree'   /
  data veg( 4) /'broadleaf evergreen tropical tree'  /
  data veg( 5) /'broadleaf evergreen temperate tree' /
  data veg( 6) /'broadleaf deciduous tropical tree'  /
  data veg( 7) /'broadleaf deciduous temperate tree' /
  data veg( 8) /'broadleaf deciduous boreal tree'    /
  data veg( 9) /'broadleaf evergreen shrub'          /
  data veg(10) /'broadleaf deciduous temperate shrub'/
  data veg(11) /'broadleaf deciduous boreal shrub'   /
  data veg(12) /'c3 arctic grass'                    /
  data veg(13) /'c3 non-arctic grass'                /
  data veg(14) /'c4 grass'                           /
  data veg(15) /'corn'                               /
  data veg(16) /'wheat'                              /
  character(len=32) :: subname = 'mkpft'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make PFTs .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input PFT file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read PCT_PFT

  call getfil (fpft, locfn, 0)

  call read_domain(tdomain,locfn)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(locfn, 0, ncid), subname)

  call check_ret(nf_inq_dimid  (ncid, 'pft', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, numpft_i), subname)

  if (numpft_i .ne. numpft+1) then
     write(6,*)'MKPFT: parameter numpft+1= ',numpft+1, &
          'does not equal input dataset numpft= ',numpft_i
     call abort()
  endif

  allocate(pctpft_i(nlon_i,nlat_i,0:numpft), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'PCT_PFT', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, pctpft_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Compute pctlnd_o, pctpft_o

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat))
  allocate( fld_i(nlon_i,nlat_i), fld_o(lsmlon,lsmlat))

  mask_i = 1.0_r8
  mask_o = 1.0_r8
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  mask_i = 100._r8*float(tdomain%mask(:,:))
  call areaave(mask_i,pctlnd_o,tgridmap)

  mask_i = mask_i   / 100._r8
  mask_o = pctlnd_o / 100._r8
  ldomain%frac = mask_o
  call gridmap_clean(tgridmap)
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  ! Area-average percent cover on input grid [pctpft_i] to output grid 
  ! [pctpft_o] and correct [pctpft_o] according to land landmask
  ! Note that percent cover is in terms of total grid area.

  do m = 0, numpft
     fld_i(:,:) = pctpft_i(:,:,m)
     call areaave(fld_i,fld_o,tgridmap)
     pctpft_o(:,:,m) = fld_o(:,:)
     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
       if (pctlnd_o(io,jo) < 1.0e-6) then
          if (m == 0) then
             pctpft_o(io,jo,m) = 100._r8
          else
             pctpft_o(io,jo,m) = 0._r8
          endif
       end if
     enddo
     enddo
  enddo

  ! Error check: percents should sum to 100 for land grid cells

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     wst_sum = 0.
     do m = 0, numpft
        wst_sum = wst_sum + pctpft_o(io,jo,m)
     enddo
     if (abs(wst_sum-100.) > 0.000001_r8) then
        write (6,*) 'MKPFT error: pft = ', &
             (pctpft_o(io,jo,m), m = 0, numpft), &
             ' do not sum to 100. at column, row = ',io,jo, &
             ' but to ', wst_sum
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

  ! Send back percent plat function types on input grid if asked for

  if ( present(pct_pft_i) )then
     allocate(pct_pft_i(nlon_i,nlat_i,0:numpft), stat=ier)
     if (ier/=0) call abort()
     pct_pft_i(:,:,:) =  pctpft_i(:,:,:)
  end if

  ! -----------------------------------------------------------------
  ! Error check1
  ! Compare global sum fld_o to global sum fld_i.
  ! -----------------------------------------------------------------

  if ( trim(mksrf_gridtype) == 'global') then
     if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
        write (6,*) 'MKPFT error: input field not conserved'
        write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
        write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
        stop
     end if
  end if

  ! -----------------------------------------------------------------
  ! Error check2
  ! Compare global areas on input and output grids
  ! -----------------------------------------------------------------

  ! input grid

  gpft_i(:) = 0.
  garea_i = 0.
  do ji = 1, nlat_i
  do ii = 1, nlon_i
     garea_i = garea_i + tdomain%area(ii,ji)
     do m = 0, numpft
        gpft_i(m) = gpft_i(m) + pctpft_i(ii,ji,m)*tdomain%area(ii,ji) * &
                                tdomain%frac(ii,ji)
     end do
  end do
  end do

  ! output grid

  gpft_o(:) = 0.
  garea_o = 0.
  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     garea_o = garea_o + ldomain%area(io,jo)
     do m = 0, numpft
        gpft_o(m) = gpft_o(m) + pctpft_o(io,jo,m)*ldomain%area(io,jo) * &
                                ldomain%frac(io,jo)
     end do
  end do
  end do

  ! comparison

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'PFTs Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,1001)
1001 format (1x,'plant type     ',20x,' input grid area',' output grid area',/ &
             1x,33x,'     10**6 km**2','      10**6 km**2')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)
  do m = 0, numpft
     write (ndiag,1002) veg(m), gpft_i(m)*1.e-06/100.,gpft_o(m)*1.e-06/100.
  end do
1002 format (1x,a35,f16.3,f17.3)

  if (lsmlat > 1) then
     k = lsmlat/2
     write (ndiag,*)
     write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
     write (ndiag,'(f10.3,a14)')ldomain%area(1,k)*1.e-06,' x 10**6 km**2'
     write (ndiag,*)
  endif
  call shr_sys_flush(ndiag)

  write (6,*) 'Successfully made PFTs'
  write (6,*)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (pctpft_i)
  deallocate (mask_i,mask_o)
  deallocate ( fld_i, fld_o)

end subroutine mkpft

end module mkpftMod
