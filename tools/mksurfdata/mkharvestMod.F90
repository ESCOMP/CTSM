module mkharvestMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkharvest
!
! !DESCRIPTION:
! Make harvest and grazing data to add to the dynamic PFT file.
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!-----------------------------------------------------------------------
! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_sys_mod  , only : shr_sys_flush

  implicit none

  private

! !PUBLIC MEMBER FUNCTIONS:
  public mkharvest_init       ! Initialization
  public mkharvest            ! Calculate the harvest values on output grid
  public mkharvest_fieldname  ! Field name
  public mkharvest_longname   ! Long name
  public mkharvest_numtypes   ! Number of harvest types

! !PRIVATE DATA MEMBERS:

  integer, parameter :: numharv = 6   ! number of harvest and grazing fields
  integer, parameter :: harlen  = 12  ! length of strings for harvest fieldnames
  character(len=harlen), parameter  :: harvest_fieldnames(numharv) = (/ &
                                                        'HARVEST_VH1',  &
                                                        'HARVEST_VH2',  &
                                                        'HARVEST_SH1',  &
                                                        'HARVEST_SH2',  &
                                                        'HARVEST_SH3',  &
                                                        'GRAZING    '   &
                                                      /)
  character(len=CL), parameter :: string_undef = 'STRING_UNDEFINED'
  character(len=CL), save :: harvest_longnames(numharv) = string_undef


!EOP
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_init
!
! !INTERFACE:
  subroutine mkharvest_init( lsmlon, lsmlat, init_val, harvest, fharvest )
!
! !DESCRIPTION:
!              Initialization of mkharvest module.
!
! !USES:
    use fileutils    , only : getfil
    use ncdio        , only : nf_open, check_ret, nf_inq_varid, &
                              nf_close, nf_get_att_text
    implicit none
!
! !ARGUMENTS:
    integer,   intent(in)  :: lsmlon, lsmlat       ! clm output grid resolution
    real(r8),  intent(in)  :: init_val             ! initial value to set to
    real(r8),  pointer     :: harvest(:,:,:)       ! output grid: normalized harvesting
    character(len=*), intent(in)  :: fharvest      ! input harvest dataset file name
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'mkharvest_init'
    integer  :: ncid,varid                      ! input netCDF id's
    character(len=256) locfn                    ! local dataset file name
    integer  :: ifld                            ! indices
!EOP
!-----------------------------------------------------------------------

    allocate( harvest(lsmlon,lsmlat,numharv) )

    harvest(:,:,:) = init_val

    call getfil (fharvest, locfn, 0)

    call check_ret(nf_open(locfn, 0, ncid), subname)

    do ifld = 1, numharv
       call check_ret(nf_inq_varid (   ncid, mkharvest_fieldname(ifld), varid),            subname)
       call check_ret(nf_get_att_text( ncid, varid, 'long_name', harvest_longnames(ifld)), subname)
    end do

    call check_ret(nf_close(ncid), subname)

  end subroutine mkharvest_init

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_fieldname
!
! !INTERFACE:
  character(len=harlen) function mkharvest_fieldname( nfield )
!
! !DESCRIPTION:
!              Return harvest fieldname of input field number.
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    integer, intent(in) :: nfield
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'mkharvest_fieldname'
!EOP
!-----------------------------------------------------------------------

      if (      nfield < 1       )then
         write(6,*) subname, ' error nfield < 1'
         call abort()
      else if ( nfield > numharv )then
         write(6,*) subname, ' error nfield > max fields'
         call abort()
      else
         mkharvest_fieldname = harvest_fieldnames(nfield)
      end if

  end function mkharvest_fieldname

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_longname
!
! !INTERFACE:
  character(len=CL) function mkharvest_longname( nfield )
!
! !DESCRIPTION:
!              Return longname description of given input field number.
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    integer, intent(in) :: nfield
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'mkharvest_longname'
!EOP
!-----------------------------------------------------------------------

      if (      nfield < 1       )then
         write(6,*) subname, ' error nfield < 1'
         call abort()
      else if ( nfield > numharv )then
         write(6,*) subname, ' error nfield > max fields'
         call abort()
      else
         if ( trim(harvest_longnames(nfield)) .eq. trim(string_undef) )then
            write(6,*) subname, ' error harvest_longnames not set yet'
            call abort()
         end if
         mkharvest_longname = harvest_longnames(nfield)
      end if

  end function mkharvest_longname

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest_numtypes
!
! !INTERFACE:
  integer function mkharvest_numtypes( )
!
! !DESCRIPTION:
!              Return number of different harvest field types.
!
! !USES:
    implicit none
!
! !ARGUMENTS:
    character(len=*), parameter :: subname = 'mkharvest_numtypes'
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
      mkharvest_numtypes = numharv

  end function mkharvest_numtypes

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkharvest
!
! !INTERFACE:
subroutine mkharvest(lsmlon, lsmlat, fharvest, ndiag, harv_o)
!
! !DESCRIPTION:
! Make harvest data for the dynamic PFT dataset.
! This dataset consists of the normalized harvest or grazing fraction (0-1) of
! the model.
!
! !USES:
  use fileutils    , only : getfil
  use domainMod    , only : domain_type, domain_clean, domain_setptrs
  use creategridMod, only : read_domain
  use mkvarsur     , only : ldomain
  use areaMod      , only : areaini, areaave, gridmap_type, gridmap_clean  
  use ncdio        , only : nf_open, check_ret, nf_inq_varid, nf_get_var_double, &
                            nf_close, nf_get_att_text
!
! !ARGUMENTS:
  implicit none
  integer ,         intent(in)  :: lsmlon, lsmlat                ! clm output grid resolution
  character(len=*), intent(in)  :: fharvest                      ! input harvest dataset file name
  integer ,         intent(in)  :: ndiag                         ! unit number for diag out
  real(r8),         intent(out) :: harv_o(lsmlon,lsmlat,numharv) ! output grid: normalized harvesting
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: nlon_i                          ! input grid : lon points
  integer  :: nlat_i                          ! input grid : lat points

  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: harv_i(:,:,:)      ! input grid: harvest/grazing percent
  real(r8), allocatable :: mask_i(:,:)        ! input grid: mask (0, 1)
  real(r8), allocatable :: mask_o(:,:)        ! output grid: mask (0, 1)
  real(r8), allocatable :: pctlnd_o(:,:)      ! output grid: percent land 
  real(r8) :: gharv_o(numharv)                ! output grid: global area harvesting
  real(r8) :: garea_o                         ! output grid: global area
  real(r8) :: gharv_i(numharv)                ! input grid: global area harvesting
  real(r8) :: garea_i                         ! input grid: global area

  integer  :: ifld                            ! indices
  integer  :: ii,ji                           ! indices
  integer  :: io,jo                           ! indices
  integer  :: k,n,m                           ! indices
  integer  :: ncid,varid                      ! input netCDF id's
  integer  :: ier                             ! error status
  character(len=256) locfn                    ! local dataset file name

  character(len=*), parameter :: unit = '10**6 km**2' ! Output units
  real(r8), parameter :: fac = 1.e-06_r8              ! Output factor
  real(r8), parameter :: rat = fac/100._r8            ! Output factor divided by 100%
  character(len=*), parameter :: subname = 'mkharvest'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make harvest fields .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input harvesting file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read HARVEST_VH1, HARVEST_VH2, ... GRAZING etc.

  call getfil (fharvest, locfn, 0)

  call read_domain(tdomain,locfn)
  call domain_setptrs(tdomain,ni=nlon_i,nj=nlat_i)

  call check_ret(nf_open(locfn, 0, ncid), subname)

  allocate(harv_i(nlon_i,nlat_i,1:numharv), stat=ier)
  if (ier/=0) call abort()

  do ifld = 1, numharv
     call check_ret(nf_inq_varid (     ncid, mkharvest_fieldname(ifld), varid),             subname)
     call check_ret(nf_get_var_double (ncid, varid, harv_i(:,:,ifld)),                      subname)
  end do

  call check_ret(nf_close(ncid), subname)

  ! Compute pctlnd_o, harv_o

  allocate(mask_i(nlon_i,nlat_i),mask_o(lsmlon,lsmlat),pctlnd_o(lsmlon,lsmlat))

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


  ! Area-average normalized harvest on input grid [harv_i] to output grid [harv_o]

  call areaave(harv_i,harv_o,tgridmap)

  ! -----------------------------------------------------------------
  ! Error check
  ! Compare global areas on input and output grids
  ! -----------------------------------------------------------------

  ! input grid

  gharv_i(:) = 0.
  garea_i = 0.
  do ji = 1, nlat_i
  do ii = 1, nlon_i
     garea_i = garea_i + tdomain%area(ii,ji)
     do m = 1, numharv
        gharv_i(m) = gharv_i(m) + harv_i(ii,ji,m)*tdomain%area(ii,ji) * &
                                tdomain%frac(ii,ji)
     end do
  end do
  end do

  ! output grid

  gharv_o(:) = 0.
  garea_o = 0.
  do jo = 1, ldomain%nj
  do io = 1, ldomain%numlon(jo)
     garea_o = garea_o + ldomain%area(io,jo)
     do m = 1, numharv
        gharv_o(m) = gharv_o(m) + harv_o(io,jo,m)*ldomain%area(io,jo) * &
                                ldomain%frac(io,jo)
     end do
  end do
  end do

  !
  ! Write out to diagnostic output file
  !

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'Harvesting Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,1001) unit, unit
1001 format (1x,'harvest type   ',20x,' input grid area',' output grid area',/ &
             1x,33x,'     ',A,'      ',A)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)
  do m = 1, numharv
     write (ndiag,1002) mkharvest_fieldname(m), gharv_i(m)*rat,gharv_o(m)*rat
  end do
1002 format (1x,a35,f16.3,f17.3)

  if (lsmlat > 1) then
     k = lsmlat/2
     write (ndiag,*)
     write (ndiag,*) 'For reference the area on the output grid of a cell near the equator is: '
     write (ndiag,'(f10.3,2a)')ldomain%area(1,k)*fac,' x ', unit
     write (ndiag,*)
  endif
  call shr_sys_flush(ndiag)

  write (6,*) 'Successfully made harvest and grazing'
  write (6,*)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (harv_i)
  deallocate (mask_i,mask_o)
  deallocate( pctlnd_o )
!-----------------------------------------------------------------------

end subroutine mkharvest

!-----------------------------------------------------------------------

end module mkharvestMod
