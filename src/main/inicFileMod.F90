#include <misc.h>
#include <preproc.h>

module inicFileMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: inicFileMod
!
! !DESCRIPTION:
! Read CLM initial data netCDF files (used only in CAM perpetual mode)
!
! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use spmdMod        , only : masterproc
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  use ncdio

! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: inicfile_perp
!
! !REVISION HISTORY:
! Jan/2004 Mariana Vertenstein: Creation
! Jan/2010 Li Xu: Modified to correct ncd_iolocal and snow_fraction
!
!EOP
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: inicfile_perp
!
! !INTERFACE:
  subroutine inicfile_perp
!
! !DESCRIPTION: 
! Read perpetual initial data fields
!
! !USES:
    use clmtype
    use clm_varctl  , only : finidat
    use clm_varpar  , only : nlevsno, nlevgrnd, nlevsoi, nlevlak
    use clm_varcon  , only : denice, denh2o, zlnd, isturb
    use fileutils   , only : getfil
    use decompMod   , only : get_proc_bounds, get_proc_global
    use fileutils   , only : getfil
    use shr_sys_mod , only : shr_sys_flush
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: ncid                           ! netCDF dataset id
    integer :: j,c,l                          ! indices
    integer :: begp, endp                     ! per-proc beginning and ending pft indices
    integer :: begc, endc                     ! per-proc beginning and ending column indices
    integer :: begl, endl                     ! per-proc beginning and ending landunit indices
    integer :: begg, endg                     ! per-proc gridcell ending gridcell indices
    integer :: numg                           ! total number of gridcells across all processors
    integer :: numl                           ! total number of landunits across all processors
    integer :: numc                           ! total number of column across all processors
    integer :: nump                           ! total number of pfts across all processors
    logical :: readvar                        ! determine if variable is on initial file
    type(gridcell_type), pointer :: gptr      ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr      ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr      ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr      ! pointer to pft derived subtype
    integer , pointer :: clandunit(:)         ! landunit index associated with each column
    integer , pointer :: itypelun(:)          ! landunit type
    character(len=256), save :: loc_fni       ! local finidat name
    logical, save :: opened_finidat = .false. ! true => finidat was opened for read
    character(len= 32) :: subname='inicfile_perp'  ! subroutine name
!------------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (landunit-level)

    itypelun            => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    clandunit           => clm3%g%l%c%landunit

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Determine necessary processor subgrid bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    ! Read initial dataset

    if (masterproc) then
       if (.not. opened_finidat) then
          call getfil(finidat, loc_fni, 0)
          call check_ret(nf_open(loc_fni, nf_nowrite, ncid), subname)
	  write(iulog,*)trim(subname),': opened netcdf file ',loc_fni
	  call shr_sys_flush(iulog)

          call check_dim(ncid, 'gridcell', numg)
          call check_dim(ncid, 'landunit', numl)
          call check_dim(ncid, 'column'  , numc)
          call check_dim(ncid, 'pft'     , nump)
          call check_dim(ncid, 'levsno'  , nlevsno)
          call check_dim(ncid, 'levgrnd' , nlevgrnd)
          call check_dim(ncid, 'levlak'  , nlevlak) 
          opened_finidat = .true.
       else
          call check_ret(nf_open(loc_fni, nf_nowrite, ncid), subname)
       end if
    end if

    call ncd_iolocal(varname='ZSNO', data=cptr%cps%z, dim1name=namec, dim2name='levsno', &
         lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='DZSNO', data=cptr%cps%dz, dim1name=namec, dim2name='levsno', &
         lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='ZISNO', data=cptr%cps%zi, dim1name=namec, dim2name='levsno', &
         lowerb2=-nlevsno, upperb2=-1, ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='H2OSNO', data=cptr%cws%h2osno, dim1name=namec, &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='SNOWDP', data=cptr%cps%snowdp, dim1name=namec, &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='SNLSNO', data=cptr%cps%snl, dim1name=namec, &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='H2OSOI_LIQ', data=cptr%cws%h2osoi_liq, dim1name=namec, dim2name='levtot', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='H2OSOI_ICE', data=cptr%cws%h2osoi_ice, dim1name=namec, dim2name='levtot', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    if (masterproc) then
       call check_ret(nf_close(ncid), subname)
    end if

    ! Determine volumetric soil water

    do j = 1,nlevsoi
       do c = begc,endc
          l = cptr%landunit(c)
          if (.not. lptr%lakpoi(l)) then
             cptr%cws%h2osoi_vol(c,j) = cptr%cws%h2osoi_liq(c,j)/(cptr%cps%dz(c,j)*denh2o) &
                                       +cptr%cws%h2osoi_ice(c,j)/(cptr%cps%dz(c,j)*denice)
          end if
       end do
    end do

    ! ============================================================================
    ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
    ! ============================================================================

    do c = begc,endc
       l = clandunit(c)
       if (itypelun(l) == isturb) then
          ! Urban landunit use Bonan 1996 (LSM Technical Note)
          cptr%cps%frac_sno(c) = min( cptr%cps%snowdp(c)/0.05_r8, 1._r8)
       else
          ! snow cover fraction in Niu et al. 2007
          cptr%cps%frac_sno(c) = 0.0_r8
          if ( cptr%cps%snowdp(c) > 0.0_r8 ) then
            cptr%cps%frac_sno(c) = tanh(cptr%cps%snowdp(c)/(2.5_r8*zlnd* &
              (min(800._r8,cptr%cws%h2osno(c)/cptr%cps%snowdp(c))/100._r8)**1._r8) )
          endif
       end if
    end do

  end subroutine inicfile_perp

end module inicFileMod
