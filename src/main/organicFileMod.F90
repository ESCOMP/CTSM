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
    use clm_varctl  , only : fsurdat, single_column
    use fileutils   , only : getfil
    use spmdMod     , only : masterproc
    use clmtype     , only : grlnd
    use domainMod   , only : llatlon
    use ncdio_pio
!
! !ARGUMENTS:
    implicit none
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
    character(len=256) :: locfn                 ! local file name
    type(file_desc_t)  :: ncid                  ! netcdf id
    integer            :: dimid                 ! dimension id
    integer            :: ni,nj,ns              ! dimension sizes  
    integer            :: ier                   ! error status
    logical            :: isgrid2d              ! true => file is 2d
    logical            :: readvar               ! true => variable is on dataset
    character(len=32)  :: subname = 'organicrd' ! subroutine name
!-----------------------------------------------------------------------

    ! Initialize data to zero - no organic matter dataset

    organic(:,:)   = 0._r8
       
    ! Read data if file was specified in namelist
       
    if (fsurdat /= ' ') then
       if (masterproc) then
          write(iulog,*) 'Attempting to read organic matter data .....'
	  write(iulog,*) subname,trim(fsurdat)
       end if

       call getfil (fsurdat, locfn, 0)
       call ncd_pio_openfile (ncid, locfn, 0)

       call ncd_inqfdims (ncid, isgrid2d, ni, nj, ns)
       if (llatlon%ns /= ns .or. llatlon%ni /= ni .or. llatlon%nj /= nj) then
          write(iulog,*)trim(subname), 'llatlon and input file do not match dims '
          write(iulog,*)trim(subname), 'llatlon%ni,ni,= ',llatlon%ni,ni
          write(iulog,*)trim(subname), 'llatlon%nj,nj,= ',llatlon%nj,nj
          write(iulog,*)trim(subname), 'llatlon%ns,ns,= ',llatlon%ns,ns
          call endrun()
       end if
       
       call ncd_io(ncid=ncid, varname='ORGANIC', flag='read', data=organic, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun('organicrd: errror reading ORGANIC')

       if ( masterproc )then
          write(iulog,*) 'Successfully read organic matter data'
          write(iulog,*)
       end if
    endif

  end subroutine organicrd

end module organicFileMod
