module mklaiMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mklai
!
! !DESCRIPTION:
! Make LAI/SAI/height data
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!-----------------------------------------------------------------------
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use areaMod     , only : gridmap_type
  use mkvarctl    

  implicit none

  private

  public  :: mklai

  private :: pft_laicheck

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mklai
!
! !INTERFACE:
subroutine mklai(lsmlon, lsmlat, fname, firrig, ndiag, ncido, &
                 ni, nj, pctpft_i )
!
! !DESCRIPTION:
! Make LAI/SAI/height data
! Portions of this code could be moved out of the month loop
! for improved efficiency
!
! !USES:
  use domainMod   , only : domain_type,domain_clean,domain_setptrs
  use creategridMod, only : read_domain
  use mkvarsur    , only : ldomain
  use mkvarpar    , only : numstdpft
  use areaMod     , only : areaini,areaave, gridmap_clean
  use ncdio
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: lsmlon, lsmlat          ! clm grid resolution
  character(len=256), intent(in) :: fname         ! input dataset file name
  character(len=256), intent(in) :: firrig        ! %irrigated area filename
  integer , intent(in) :: ndiag                   ! unit number for diag out
  integer , intent(in) :: ncido                   ! output netcdf file id
  integer , intent(in) :: ni                      ! number of long dimension of pft index
  integer , intent(in) :: nj                      ! number of lat dimension of pft index
  real(r8), intent(in) :: pctpft_i(ni,nj,0:numpft)! % plant function types on input grid
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

  integer  :: numpft_i                        ! number of plant types on input
  real(r8) :: glai_o(0:numpft)                ! output grid: global area pfts
  real(r8) :: gsai_o(0:numpft)                ! output grid: global area pfts
  real(r8) :: ghgtt_o(0:numpft)               ! output grid: global area pfts
  real(r8) :: ghgtb_o(0:numpft)               ! output grid: global area pfts
  real(r8) :: glai_i(0:numpft)                ! input grid: global area pfts
  real(r8) :: gsai_i(0:numpft)                ! input grid: global area pfts
  real(r8) :: ghgtt_i(0:numpft)               ! input grid: global area pfts
  real(r8) :: ghgtb_i(0:numpft)               ! input grid: global area pfts

  real(r8), allocatable :: mlai_o(:,:,:)      ! monthly lai
  real(r8), allocatable :: msai_o(:,:,:)      ! monthly sai
  real(r8), allocatable :: mhgtt_o(:,:,:)     ! monthly height (top)
  real(r8), allocatable :: mhgtb_o(:,:,:)     ! monthly height (bottom)
  real(r8), allocatable :: mlai_max(:,:,:)    ! monthly lai
  real(r8), allocatable :: msai_max(:,:,:)    ! monthly sai
  real(r8), allocatable :: mhgtt_max(:,:,:)   ! monthly height (top)
  real(r8), allocatable :: mhgtb_max(:,:,:)   ! monthly height (bottom)
  real(r8), allocatable :: mlai_i(:,:,:)      ! monthly lai in
  real(r8), allocatable :: msai_i(:,:,:)      ! monthly sai in
  real(r8), allocatable :: mhgtt_i(:,:,:)     ! monthly height (top) in
  real(r8), allocatable :: mhgtb_i(:,:,:)     ! monthly height (bottom) in
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  real(r8), allocatable :: fld_i(:,:)         ! input grid: dummy field
  real(r8), allocatable :: fld_o(:,:)         ! output grid: dummy field
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: garea_o                         ! output grid: global area

  integer,  allocatable :: laimask(:,:,:)     ! lai+sai output mask for each plant function type
  integer  :: mwts                            ! number of weights
  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,l,n,m                         ! indices
  integer  :: ncidi,dimid,varid               ! input netCDF id's
  integer  :: beg4d(4),len4d(4)               ! netCDF variable edges
  integer  :: dim4id(4)                       ! netcdf ids
  integer  :: ntim                            ! number of input time samples
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=256) :: name                  ! name of attribute
  character(len=256) :: unit                  ! units of attribute
  character(len= 32) :: subname = 'mklai'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make LAIs/SAIs/heights .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call read_domain(tdomain,fname)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(fname, 0, ncidi), subname)

  call check_ret(nf_inq_dimid(ncidi, 'pft', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, numpft_i), subname)
  if (numpft_i /= numpft+1) then
     write(6,*)'MKLAI: parameter numpft+1= ',numpft+1, &
          'does not equal input dataset numpft= ',numpft_i
     stop
  endif

  call check_ret(nf_inq_dimid(ncidi, 'time', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, ntim), subname)
  if (ntim /= 12) then
     write(6,*)'MKLAI: must have 12 time samples on input data'
     call abort()
  endif

!  --- at bottom of routine ---
!  call check_ret(nf_close(ncidi), subname)

  ! Compute local fields _o

  allocate(mlai_i(nlon_i,nlat_i,0:numpft), &
           msai_i(nlon_i,nlat_i,0:numpft), &
           mhgtt_i(nlon_i,nlat_i,0:numpft), &
           mhgtb_i(nlon_i,nlat_i,0:numpft), stat=ier)
  if (ier /= 0) then
     write(6,*)'mklai allocation error'; call abort()
  end if

  allocate(mlai_o(lsmlon,lsmlat,0:numpft), &
           msai_o(lsmlon,lsmlat,0:numpft), &
           mhgtt_o(lsmlon,lsmlat,0:numpft), &
           mhgtb_o(lsmlon,lsmlat,0:numpft), stat=ier)
  if (ier /= 0) then
     write(6,*)'mklai allocation error'; call abort()
  end if

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat), stat=ier)
  if (ier /= 0) then
     write(6,*)'mklai allocation error'; call abort()
  end if
  allocate( fld_i(nlon_i,nlat_i), fld_o(lsmlon,lsmlat), stat=ier)
  if (ier /= 0) then
     write(6,*)'mklai allocation error'; call abort()
  end if

  mask_i = 1.0_r8
  mask_o = 1.0_r8
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  mask_i = float(tdomain%mask(:,:))
  call areaave(mask_i,mask_o,tgridmap)

  call gridmap_clean(tgridmap)
  call areaini(tdomain,ldomain,tgridmap,fracin=mask_i,fracout=mask_o)

  allocate( laimask(ni,nj,0:numpft), stat=ier )
  if (ier /= 0) then
     write(6,*)'mklai allocation error'; call abort()
  end if
  laimask(:,:,:) = 0
  do m = 1, ntim

     ! Get input data for the month

     beg4d(1) = 1 ; len4d(1) = nlon_i
     beg4d(2) = 1 ; len4d(2) = nlat_i
     beg4d(3) = 1 ; len4d(3) = numpft+1
     beg4d(4) = m ; len4d(4) = 1

     call check_ret(nf_inq_varid (ncidi, 'MONTHLY_LAI', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, beg4d, len4d, mlai_i(:,:,0:numpft)), subname)

     call check_ret(nf_inq_varid (ncidi, 'MONTHLY_SAI', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, beg4d, len4d, msai_i(:,:,0:numpft)), subname)

     call check_ret(nf_inq_varid (ncidi, 'MONTHLY_HEIGHT_TOP', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, beg4d, len4d, mhgtt_i(:,:,0:numpft)), subname)

     call check_ret(nf_inq_varid (ncidi, 'MONTHLY_HEIGHT_BOT', varid), subname)
     call check_ret(nf_get_vara_double (ncidi, varid, beg4d, len4d, mhgtb_i(:,:,0:numpft)), subname)

     mlai_o(:,:,:)  = 0.
     msai_o(:,:,:)  = 0.
     mhgtt_o(:,:,:) = 0.
     mhgtb_o(:,:,:) = 0.

     do l = 0, numpft

        fld_i(:,:) = mlai_i(:,:,l)
        call areaave(    fld_i,fld_o,tgridmap)
        mlai_o(:,:,l) = fld_o(:,:)

        fld_i(:,:) = msai_i(:,:,l)
        call areaave    (fld_i,fld_o,tgridmap)
        msai_o(:,:,l) = fld_o(:,:)

        fld_i(:,:) = mhgtt_i(:,:,l)
        call areaave    (fld_i,fld_o,tgridmap)
        mhgtt_o(:,:,l) = fld_o(:,:)

        fld_i(:,:) = mhgtb_i(:,:,l)
        call areaave    (fld_i,fld_o,tgridmap)
        mhgtb_o(:,:,l) = fld_o(:,:)

!tcx?        where (ldomain%mask(:,:) == 0)
!           mlai_o (:,:,l) = 0.
!           msai_o (:,:,l) = 0.
!           mhgtt_o(:,:,l) = 0.
!           mhgtb_o(:,:,l) = 0.
!        endwhere
     enddo

     ! if irrigation dataset present, copy LAI,SAI,Heights from PFT=15 (non-irrigated) 
     ! into PFT=16 (irrigated)
     if (firrig /= ' ') then      
        write(6,*) 'Irrigation dataset present; Copying crop (PFT=15) LAI, SAI, and heights',&
                   ' into irrigated crop (PFT=16) '
        mlai_o(:,:,16)  = mlai_o(:,:,15)
        msai_o(:,:,16)  = msai_o(:,:,15)
        mhgtt_o(:,:,16) = mhgtt_o(:,:,15)
        mhgtb_o(:,:,16) = mhgtb_o(:,:,15)
     endif

     ! -----------------------------------------------------------------
     ! Output model resolution LAI/SAI/HEIGHT data
     ! -----------------------------------------------------------------

     beg4d(1) = 1  ;  len4d(1) = lsmlon
     beg4d(2) = 1  ;  len4d(2) = lsmlat
     beg4d(3) = 1  ;  len4d(3) = numpft+1
     beg4d(4) = m  ;  len4d(4) = 1
     
     call check_ret(nf_inq_varid(ncido, 'MONTHLY_LAI', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, beg4d, len4d, mlai_o), subname)
     
     call check_ret(nf_inq_varid(ncido, 'MONTHLY_SAI', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, beg4d, len4d, msai_o), subname)
     
     call check_ret(nf_inq_varid(ncido, 'MONTHLY_HEIGHT_TOP', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, beg4d, len4d, mhgtt_o), subname)
     
     call check_ret(nf_inq_varid(ncido, 'MONTHLY_HEIGHT_BOT', varid), subname)
     call check_ret(nf_put_vara_double(ncido, varid, beg4d, len4d, mhgtb_o), subname)

     call check_ret(nf_inq_varid(ncido, 'time', varid), subname)
     call check_ret(nf_put_vara_int(ncido, varid, beg4d(4), len4d(4), m), subname)

     call check_ret(nf_sync(ncido), subname)

     ! -----------------------------------------------------------------
     ! Error check2
     ! Compare global areas on input and output grids
     ! -----------------------------------------------------------------

     ! Input grid global area

     glai_i(:)  = 0.
     gsai_i(:)  = 0.
     ghgtt_i(:) = 0.
     ghgtb_i(:) = 0.
     garea_i    = 0.

     do ji = 1, nlat_i
     do ii = 1, nlon_i
        garea_i = garea_i + tdomain%area(ii,ji)
     end do
     end do

     do l = 0, numpft
     do ji = 1, nlat_i
     do ii = 1, nlon_i
        glai_i(l)  = glai_i(l) + mlai_i(ii,ji,l)*tdomain%area(ii,ji) * &
                                 tdomain%frac(ii,ji)
        gsai_i(l)  = gsai_i(l) + msai_i(ii,ji,l)*tdomain%area(ii,ji) * &
                                 tdomain%frac(ii,ji)
        ghgtt_i(l) = ghgtt_i(l)+ mhgtt_i(ii,ji,l)*tdomain%area(ii,ji) * &
                                 tdomain%frac(ii,ji)
        ghgtb_i(l) = ghgtb_i(l)+ mhgtb_i(ii,ji,l)*tdomain%area(ii,ji) * &
                                 tdomain%frac(ii,ji)
     end do
     end do
     end do

     ! Output grid global area

     glai_o(:)  = 0.
     gsai_o(:)  = 0.
     ghgtt_o(:) = 0.
     ghgtb_o(:) = 0.
     garea_o    = 0.

     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        garea_o = garea_o + ldomain%area(io,jo)
     end do
     end do

     do l = 0, numpft
     do jo = 1, ldomain%nj
     do io = 1, ldomain%numlon(jo)
        glai_o(l)  = glai_o(l) + mlai_o(io,jo,l)*ldomain%area(io,jo) * &
                                 ldomain%frac(io,jo)
        gsai_o(l)  = gsai_o(l) + msai_o(io,jo,l)*ldomain%area(io,jo) * &
                                 ldomain%frac(io,jo)
        ghgtt_o(l) = ghgtt_o(l)+ mhgtt_o(io,jo,l)*ldomain%area(io,jo) * &
                                 ldomain%frac(io,jo)
        ghgtb_o(l) = ghgtb_o(l)+ mhgtb_o(io,jo,l)*ldomain%area(io,jo) * &
                                 ldomain%frac(io,jo)
     end do
     end do
     end do

     ! Comparison

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('=',k=1,70)
     write (ndiag,*) 'LAI Output for month ',m
     write (ndiag,'(1x,70a1)') ('=',k=1,70)

     write (ndiag,*)
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,1001)
1001 format (1x,'PFT input grid area output grid area',/ &
             1x,3x,'     10**6 km**2','      10**6 km**2')
     write (ndiag,'(1x,70a1)') ('.',k=1,70)
     write (ndiag,*)

     do l = 0, numpft
        write (ndiag,1002) l, glai_i(l)*1.e-06*1.e-02,glai_o(l)*1.e-06*1.e-02
1002    format (1x,i3,f16.3,f17.3)
     end do

     write (6,*) 'Successfully made LAIs/SAIs/heights for month ', m
     write (6,*)
     call shr_sys_flush(6)

  enddo

  call check_ret(nf_close(ncidi), subname)

  ! consistency check that PFT and LAI+SAI make sense

  !call pft_laicheck( ni, nj, pft_i, laimask )

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (mlai_o,msai_o,mhgtt_o,mhgtb_o)
  deallocate (mlai_i,msai_i,mhgtt_i,mhgtb_i)
  deallocate (mask_i,mask_o)
  deallocate ( fld_i, fld_o)

end subroutine mklai

!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
subroutine pft_laicheck( ni, nj, pctpft_i, laimask )
! !USES:
!
! !DESCRIPTION:
!
! consistency check that PFT and LAI+SAI make sense
!
! !ARGUMENTS:
  implicit none
  integer ,           intent(in) :: ni, nj                  ! input PFT grid resolution
  real(r8),           intent(in) :: pctpft_i(ni,nj,0:numpft)! % plant function types
  integer,            intent(in) :: laimask(ni,nj,0:numpft) ! mask where LAI+SAI > 0
!EOP

  character(len=*), parameter :: subName="pft_laicheck"
  integer :: ii, ji, l, n, nc      ! Indices
!-----------------------------------------------------------------------

  do l  = 0, numpft
     n  = 0
     nc = 0
     do ji = 1, nj
     do ii = 1, ni
        if ( pctpft_i(ii,ji,l) > 0.0_r8 ) nc = nc + 1
        if ( (pctpft_i(ii,ji,l) > 0.0_r8) .and. (laimask(ii,ji,l) /= 1) )then
           write (6,*) subName//' :: warning: pft and LAI+SAI mask not consistent!'
           write (6,*) 'ii,jj,l   = ', ii, ji, l
           write (6,*) 'pctpft_i  = ',pctpft_i(ii,ji,l)
           write (6,*) 'laimask   = ', laimask(ii,ji,l)
           n = n + 1
        end if
     end do
     end do
     if ( n > max(4,nc/4) ) then
        write (6,*) subName//' :: pft/LAI+SAI inconsistency over more than 25% land-cover'
        write (6,*) '# inconsistent points, total PFT pts, total LAI+SAI pts = ', &
                     n, nc, sum(laimask(:,:,l))
        stop
     end if
  end do

end subroutine pft_laicheck

!-----------------------------------------------------------------------
  
end module mklaiMod
