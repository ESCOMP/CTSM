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
  use ncdio

! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: inicfile_perp
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 1/04
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
    use clm_varpar  , only : nlevsno, nlevsoi, nlevlak
    use clm_varcon  , only : denice, denh2o, zlnd
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
!EOP
!
! !LOCAL VARIABLES:
    integer :: ncid                           ! netCDF dataset id
    integer :: j,c,l                          ! indices
    integer :: begp, endp                     ! per-proc beginning and ending pft indices
    integer :: begc, endc                     ! per-proc beginning and ending column indices
    integer :: begl, endl                     ! per-proc beginning and ending landunit indices
    integer :: begg, endg                     ! per-proc gridcell ending gridcell indices
    integer :: numg                           ! total number of gridcells across all processors
    integer :: numl                           ! total number of landunits across all processors
    integer :: numc                           ! total number of columns across all processors
    integer :: nump                           ! total number of pfts across all processors
    logical :: readvar                        ! determine if variable is on initial file
    type(gridcell_type), pointer :: gptr      ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr      ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr      ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr      ! pointer to pft derived subtype
    character(len=256), save :: loc_fni       ! local finidat name
    logical, save :: opened_finidat = .false. ! true => finidat was opened for read
    character(len= 32) :: subname='inicfile_perp'  ! subroutine name
!------------------------------------------------------------------------

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
	  write(6,*)trim(subname),': opened netcdf file ',loc_fni
#ifndef UNICOSMP
	  call shr_sys_flush(6)
#endif
          call check_dim(ncid, 'gridcell', numg)
          call check_dim(ncid, 'landunit', numl)
          call check_dim(ncid, 'column'  , numc)
          call check_dim(ncid, 'pft'     , nump)
          call check_dim(ncid, 'levsno'  , nlevsno)
          call check_dim(ncid, 'levsoi'  , nlevsoi)
          call check_dim(ncid, 'levlak'  , nlevlak) 
          opened_finidat = .true.
       else
          call check_ret(nf_open(loc_fni, nf_nowrite, ncid), subname)
       end if
    end if

    call ncd_iolocal(varname='ZSNO', data=cptr%cps%z, dim1name='column', dim2name='levsno', &
         lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='DZSNO', data=cptr%cps%dz, dim1name='column', dim2name='levsno', &
         lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='ZISNO', data=cptr%cps%zi, dim1name='column', dim2name='levsno', &
         lowerb2=-nlevsno, upperb2=-1, ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='H2OSNO', data=cptr%cws%h2osno, dim1name='column', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='SNOWDP', data=cptr%cps%snowdp, dim1name='column', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='SNOWAGE', data=cptr%cps%snowage, dim1name='column', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='SNLSNO', data=cptr%cps%snl, dim1name='column', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='H2OSOI_LIQ', data=cptr%cws%h2osoi_liq, dim1name='column', dim2name='levtot', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    call ncd_iolocal(varname='H2OSOI_ICE', data=cptr%cws%h2osoi_ice, dim1name='column', dim2name='levtot', &
         ncid=ncid, flag='read', readvar=readvar)
    if (.not. readvar) call endrun()

    if (masterproc) then
       call check_ret(nf_close(ncid), subname)
    end if

    ! Determine volumetric soil water

    do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          l = cptr%landunit(c)
          if (.not. lptr%lakpoi(l)) then
             cptr%cws%h2osoi_vol(c,j) = cptr%cws%h2osoi_liq(c,j)/(cptr%cps%dz(c,j)*denh2o) &
                                       +cptr%cws%h2osoi_ice(c,j)/(cptr%cps%dz(c,j)*denice)
          end if
       end do
    end do

!dir$ concurrent
!cdir nodep
     do c = begc,endc
        cptr%cps%frac_sno(c) = cptr%cps%snowdp(c) / (10._r8*zlnd + cptr%cps%snowdp(c))
     end do

  end subroutine inicfile_perp

end module inicFileMod
