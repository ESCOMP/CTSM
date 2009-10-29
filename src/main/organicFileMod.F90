#include <misc.h>
#include <preproc.h>

module organicFileMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: organicFileMod
!
! !DESCRIPTION:
! Contains methods for reading in organic matter data file which has 
! organic matter density for each grid point and soil level 
!
! !USES
  use abortutils   , only : endrun
  use clm_varctl   , only : iulog
  use shr_kind_mod , only : r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: organicrd  ! Read organic matter dataset
!
! !REVISION HISTORY:
! Created by David Lawrence, 4 May 2006
! Revised by David Lawrence, 21 September 2007
! Revised by David Lawrence, 14 October 2008
!
!EOP
!
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: organicrd
!
! !INTERFACE:
  subroutine organicrd(organic)
!
! !DESCRIPTION: 
! Read the organic matter dataset.
!
! !USES:
    use clm_varctl  , only : fsurdat, single_column, scmlon, scmlat
    use clm_varpar  , only : lsmlon, lsmlat, nlevsoi
    use ncdio       , only : ncd_iolocal, check_dim, check_ret
    use fileutils   , only : getfil
    use decompMod   , only : get_proc_bounds
    use spmdMod     , only : masterproc
    use clmtype     , only : grlnd
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    real(r8), pointer :: organic(:,:)         ! organic matter density (kg/m3)
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by David Lawrence, 4 May 2006
! Revised by David Lawrence, 21 September 2007
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=256) :: locfn                          ! local file name
    real(r8),pointer :: arrayl(:)                        ! generic global array
    integer  :: ncid,dimid,varid                         ! netCDF id's
    integer  :: begg,endg                                ! start/stop gridcells
    integer  :: start(3),count(3)                        ! netcdf start/count arrays
    integer  :: ier                                      ! error status 
    integer  :: n                                        ! indices
    integer  :: closelatidx,closelonidx
    real(r8) :: closelat,closelon
    character(len=32) :: subname = 'organicrd'          ! subroutine name
!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg)

    ! Initialize data to zero - no organic matter dataset

    organic(:,:)   = 0._r8
       
    ! read data if file was specified in namelist
       
    if (fsurdat /= ' ') then

       count(1) = lsmlon
       count(2) = lsmlat
       start(1) = 1
       start(2) = 1
       start(3) = 1
       count(3) = 1

       ! Obtain netcdf file and read surface data

       if (masterproc) then

          write(iulog,*) 'Attempting to read organic matter data .....'

          call getfil (fsurdat, locfn, 0)
          call check_ret(nf_open(locfn, 0, ncid), subname)
	  write(iulog,*) subname,trim(fsurdat)
          write(iulog,*) "  Expected dimensions: lsmlon=",lsmlon," lsmlat=",lsmlat
          if (.not.single_column) then
             call check_dim(ncid, 'lsmlon' , lsmlon)
             call check_dim(ncid, 'lsmlat' , lsmlat)
	  endif 
       endif 

       allocate(arrayl(begg:endg))
       do n = 1,nlevsoi 
          start(3) = n
          call ncd_iolocal(ncid,'ORGANIC','read',arrayl,grlnd,start,count)
          organic(begg:endg,n) = arrayl(begg:endg)
       enddo
       deallocate(arrayl)

    endif

    if ( masterproc )then
       write(iulog,*) 'Successfully read organic matter data'
       write(iulog,*)
    end if

  end subroutine organicrd

end module organicFileMod
