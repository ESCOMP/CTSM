module mkvocefMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkfmaxMod
!
! !DESCRIPTION:
! Make fmax for surface dataset
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
  public :: mkvocef  ! Get the percentage emissions for VOC for different
                     ! land cover types

!EOP
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkvocef
!
! !INTERFACE:
subroutine mkvocef (lsmlon, lsmlat, fvocef, ndiag, ef_btr_o,ef_fet_o,ef_fdt_o,ef_shr_o,ef_grs_o,ef_crp_o)
!
! !DESCRIPTION:
! make volatile organic coumpunds (VOC) emission factors.
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
  integer , intent(in) :: lsmlon, lsmlat               ! clm grid resolution
  character(len=*), intent(in) :: fvocef               ! input soicol dataset file name
  integer , intent(in) :: ndiag                        ! unit number for diagnostic output
  real(r8), intent(out):: ef_btr_o(lsmlon,lsmlat)      ! output grid: EFs for broadleaf trees
  real(r8), intent(out):: ef_fet_o(lsmlon,lsmlat)      ! output grid: EFs for fineleaf evergreen
  real(r8), intent(out):: ef_fdt_o(lsmlon,lsmlat)      ! output grid: EFs for fineleaf deciduous
  real(r8), intent(out):: ef_shr_o(lsmlon,lsmlat)      ! output grid: EFs for shrubs
  real(r8), intent(out):: ef_grs_o(lsmlon,lsmlat)      ! output grid: EFs for grasses
  real(r8), intent(out):: ef_crp_o(lsmlon,lsmlat)      ! output grid: EFs for crops
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Colette L. Heald
! 17 Jul 2007 F Vitt -- updated to pftintdat06_clm3_5_05 and corrected indexing of ef_*_i arrarys
!
!EOP
!
! !LOCAL VARIABLES:
  integer  :: nlon_i                          ! input grid : lon points
  integer  :: nlat_i                          ! input grid : lat points

  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: area_i(:,:)        ! input grid: cell area
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  real(r8), allocatable :: fld_i(:,:)         ! input grid: dummy field
  real(r8), allocatable :: fld_o(:,:)         ! output grid: dummy field
  real(r8) :: sum_fldo                        ! global sum of dummy input fld
  real(r8) :: sum_fldi                        ! global sum of dummy input fld

  integer, allocatable  :: temp_i(:,:)        ! input grid: EFs for broadleaf trees
  real(r8), allocatable :: ef_btr_i(:,:)      ! input grid: EFs for broadleaf trees
  real(r8), allocatable :: ef_fet_i(:,:)      ! input grid: EFs for fineleaf evergreen
  real(r8), allocatable :: ef_fdt_i(:,:)      ! input grid: EFs for fineleaf deciduous
  real(r8), allocatable :: ef_shr_i(:,:)      ! input grid: EFs for shrubs
  real(r8), allocatable :: ef_grs_i(:,:)      ! input grid: EFs for grasses
  real(r8), allocatable :: ef_crp_i(:,:)      ! input grid: EFs for crops

  integer :: jn,in                            ! latitute index for switch
  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,n,m                           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001_r8             ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mkvocef'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make VOC emission factors .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input Emission Factors
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call read_domain(tdomain,fvocef)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(fvocef, 0, ncid), subname)

  allocate(area_i(nlon_i,nlat_i),temp_i(nlon_i,nlat_i),&
       ef_btr_i(nlon_i,nlat_i), ef_fet_i(nlon_i,nlat_i), ef_fdt_i(nlon_i,nlat_i), &
       ef_shr_i(nlon_i,nlat_i), ef_grs_i(nlon_i,nlat_i), ef_crp_i(nlon_i,nlat_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_inq_varid (ncid, 'ef_btr', varid), subname)
  call check_ret(nf_get_var_int (ncid, varid, temp_i), subname)
  do in = 1,720
     do jn = 1, 360
        if (in < 361) then
           ef_btr_i(in,361-jn)=temp_i(in+360,jn)*1.0_r8
	end if
	if (in > 360) then
           ef_btr_i(in,361-jn)=temp_i(in-360,jn)*1.0_r8
	end if
     end do
  end do

  call check_ret(nf_inq_varid (ncid, 'ef_fet', varid), subname)
  call check_ret(nf_get_var_int (ncid, varid, temp_i), subname)
  do in = 1,720
     do jn = 1, 360
        if (in < 361) then
  	  ef_fet_i(in,361-jn)=temp_i(in+360,jn)*1.0_r8
	end if
	if (in > 360) then
	  ef_fet_i(in,361-jn)=temp_i(in-360,jn)*1.0_r8
	end if
     end do
  end do

  call check_ret(nf_inq_varid (ncid, 'ef_fdt', varid), subname)
  call check_ret(nf_get_var_int (ncid, varid, temp_i), subname)
  do in = 1,720
     do jn = 1, 360
        if (in < 361) then
  	  ef_fdt_i(in,361-jn)=temp_i(in+360,jn)*1.0_r8
	end if
	if (in > 360) then
	  ef_fdt_i(in,361-jn)=temp_i(in-360,jn)*1.0_r8
	end if
     end do
  end do

  call check_ret(nf_inq_varid (ncid, 'ef_shr', varid), subname)
  call check_ret(nf_get_var_int (ncid, varid, temp_i), subname)
  do in = 1,720
     do jn = 1, 360
        if (in < 361) then
  	  ef_shr_i(in,361-jn)=temp_i(in+360,jn)*1.0_r8
	end if
	if (in > 360) then
	  ef_shr_i(in,361-jn)=temp_i(in-360,jn)*1.0_r8
	end if
     end do
  end do

  call check_ret(nf_inq_varid (ncid, 'ef_grs', varid), subname)
  call check_ret(nf_get_var_int (ncid, varid, temp_i), subname)
  do in = 1,720
     do jn = 1, 360
        if (in < 361) then
  	  ef_grs_i(in,361-jn)=temp_i(in+360,jn)*1.0_r8
	end if
	if (in > 360) then
	  ef_grs_i(in,361-jn)=temp_i(in-360,jn)*1.0_r8
	end if
     end do
  end do

  call check_ret(nf_inq_varid (ncid, 'ef_crp', varid), subname)
  call check_ret(nf_get_var_int (ncid, varid, temp_i), subname)
  do in = 1,720
     do jn = 1, 360
        if (in < 361) then
  	  ef_crp_i(in,361-jn)=temp_i(in+360,jn)*1.0_r8
	end if
	if (in > 360) then
	  ef_crp_i(in,361-jn)=temp_i(in-360,jn)*1.0_r8
	end if
     end do
  end do


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

  call areaave(ef_btr_i,ef_btr_o,tgridmap)
  call areaave(ef_fet_i,ef_fet_o,tgridmap)
  call areaave(ef_fdt_i,ef_fdt_o,tgridmap)
  call areaave(ef_shr_i,ef_shr_o,tgridmap)
  call areaave(ef_grs_i,ef_grs_o,tgridmap)
  call areaave(ef_crp_i,ef_crp_o,tgridmap)

  ! Check for conservation

  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)

     if ( ef_btr_o(io,jo) < 0._r8 ) then
        write (6,*) 'MKVOCEF error: EF btr = ',ef_btr_o(io,jo), &
             ' is negative for col, row = ',io,jo
        call abort()
     end if
     if ( ef_fet_o(io,jo) < 0._r8 ) then
        write (6,*) 'MKVOCEF error: EF fet = ',ef_fet_o(io,jo), &
             ' is negative for col, row = ',io,jo
        call abort()
     end if
     if ( ef_fdt_o(io,jo) < 0._r8 ) then
        write (6,*) 'MKVOCEF error: EF fdt = ',ef_fdt_o(io,jo), &
             ' is negative for col, row = ',io,jo
        call abort()
     end if
     if ( ef_shr_o(io,jo) < 0._r8 ) then
        write (6,*) 'MKVOCEF error: EF shr = ',ef_shr_o(io,jo), &
             ' is negative for col, row = ',io,jo
        call abort()
     end if
     if ( ef_grs_o(io,jo) < 0._r8 ) then
        write (6,*) 'MKVOCEF error: EF grs = ',ef_grs_o(io,jo), &
             ' is negative for col, row = ',io,jo
        call abort()
     end if
     if ( ef_crp_o(io,jo) < 0._r8 ) then
        write (6,*) 'MKVOCEF error: EF crp = ',ef_crp_o(io,jo), &
             ' is negative for col, row = ',io,jo
        call abort()
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

  if ( trim(mksrf_gridtype) == 'global') then
     if ( abs(sum_fldo/sum_fldi-1._r8) > relerr ) then
        write (6,*) 'MKVOCEF error: input field not conserved'
        write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
        write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
        stop
     end if
  end if

  write (6,*) 'Successfully made VOC Emission Factors'
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (area_i)
  deallocate (mask_i,mask_o)
  deallocate ( fld_i, fld_o)
  deallocate ( temp_i,&
       ef_btr_i, ef_fet_i, ef_fdt_i, &
       ef_shr_i, ef_grs_i, ef_crp_i )
end subroutine mkvocef

!-----------------------------------------------------------------------

end module mkvocefMod
