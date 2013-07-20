
module CNSharedConstsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNSharedConstsMod
!
! !DESCRIPTION:
! 
! 
!
! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  implicit none
  save
! !PUBLIC MEMBER FUNCTIONS:
!
! !REVISION HISTORY:
!   Dec 3 2012 : Created by S. Muszala
!
!EOP
!-----------------------------------------------------------------------

   ! CNConstShareInst.  PGI wants the type decl. public but the instance
   ! is indeed protected.  A generic private statement at the start of the module
   ! overrides the protected functionality with PGI

   type, public  :: CNConstShareType
      real(r8) :: Q10   ! temperature dependence
      real(r8) :: froz_q10    ! separate q10 for frozen soil respiration rates
      real(r8) :: decomp_depth_efolding ! e-folding depth for reduction in decomposition (m) 
   end type CNConstShareType

   type(CNConstShareType),protected :: CNConstShareInst
  
contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNConstReadShared
!
! !INTERFACE:
subroutine CNConstReadShared(ncid)
!
! !DESCRIPTION:
! 
! 
!
! !USES:
   use shr_kind_mod , only: r8 => shr_kind_r8
   use ncdio_pio    , only : file_desc_t,ncd_io
   use abortutils   , only: endrun
!
! !ARGUMENTS:
   implicit none
   type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
!
! !CALLED FROM:clm_initializeMod.F90::initialize1 
! 
!
! !REVISION HISTORY:
!   Dec 3 2012 : Created by S. Muszala
! !LOCAL VARIABLES:
!
!
! !OTHER LOCAL VARIABLES:
   character(len=32)  :: subname = 'CNConstReadShared'
   character(len=100) :: errCode = 'Error reading in CN shared const file '
   logical            :: readv ! has variable been read in or not
   real(r8)           :: tempr ! temporary to read in constant
   character(len=100) :: tString ! temp. var for reading

!EOP
!-----------------------------------------------------------------------

   !
   ! netcdf read here
   !
   tString='q10_mr'
   call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
   if ( .not. readv ) call endrun( trim(subname)//trim(errCode)//trim(tString))
   CNConstShareInst%Q10=tempr

   tString='froz_q10'
   call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
   if ( .not. readv ) call endrun( trim(subname)//trim(errCode)//trim(tString))
   CNConstShareInst%froz_q10=tempr

   tString='decomp_depth_efolding'
   call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
   if ( .not. readv ) call endrun( trim(subname)//trim(errCode)//trim(tString))
   CNConstShareInst%decomp_depth_efolding=tempr

end subroutine CNConstReadShared

end module CNSharedConstsMod
