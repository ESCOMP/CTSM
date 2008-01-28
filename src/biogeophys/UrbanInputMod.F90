#include <misc.h>
#include <preproc.h>

module UrbanInputMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: UrbanInputMod
! 
! !DESCRIPTION: 
! Read in input urban data - fill in data structure urbinp
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun  
  use shr_sys_mod , only : shr_sys_flush 
!
! !PUBLIC TYPES:
  implicit none
  save

  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanInput         ! Read in urban input data

  type urbinp_t
     real(r8), pointer :: canyon_hwr(:)  
     real(r8), pointer :: wtlunit_roof(:)  
     real(r8), pointer :: wtroad_perv(:)  
     real(r8), pointer :: em_roof(:)   
     real(r8), pointer :: em_improad(:)  
     real(r8), pointer :: em_perroad(:)  
     real(r8), pointer :: em_wall(:)  
     real(r8), pointer :: alb_roof_dir(:,:)  
     real(r8), pointer :: alb_roof_dif(:,:)  
     real(r8), pointer :: alb_improad_dir(:,:)  
     real(r8), pointer :: alb_improad_dif(:,:)  
     real(r8), pointer :: alb_perroad_dir(:,:)  
     real(r8), pointer :: alb_perroad_dif(:,:)  
     real(r8), pointer :: alb_wall_dir(:,:)  
     real(r8), pointer :: alb_wall_dif(:,:)  
     real(r8), pointer :: ht_roof(:)
     real(r8), pointer :: wind_hgt_canyon(:)
     real(r8), pointer :: tk_wall(:,:)
     real(r8), pointer :: tk_roof(:,:)
     real(r8), pointer :: tk_improad(:,:)
     real(r8), pointer :: cv_wall(:,:)
     real(r8), pointer :: cv_roof(:,:)
     real(r8), pointer :: cv_improad(:,:)
     real(r8), pointer :: sandfrac_road(:,:)
     real(r8), pointer :: clayfrac_road(:,:)
     real(r8), pointer :: scalez_wall(:)
     real(r8), pointer :: scalez_roof(:)
     real(r8), pointer :: thick_wall(:)
     real(r8), pointer :: thick_roof(:)
  end type urbinp_t
  public urbinp_t

  type (urbinp_t)   , public :: urbinp  ! urban input derived type
  character(len=256), public :: furbinp ! urban input data
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanInput
!
! !INTERFACE:
  subroutine UrbanInput(mode)
!
! !DESCRIPTION: 
! Allocate memory and read in urban input data
!
! !USES:
    use clm_varpar, only : lsmlon, lsmlat, numrad, nlevsoi
    use clm_varctl, only : iulog
    use fileutils , only : getavu, relavu, getfil, opnfil
    use spmdMod   , only : masterproc
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: mode
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein July 2004
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: locfn       ! local file name
    character(len=32)  :: desc
    integer :: nw,n,k,i,j
    integer :: ier  
    real(r8) :: soilfrac
!-----------------------------------------------------------------------

    if (furbinp == ' ') then
       return
    end if

    if (mode == 'initialize') then

       ! Allocate dynamic memory

       allocate(urbinp%canyon_hwr(lsmlon*lsmlat), &  
                urbinp%wtlunit_roof(lsmlon*lsmlat), &  
                urbinp%wtroad_perv(lsmlon*lsmlat), &
                urbinp%em_roof(lsmlon*lsmlat), &     
                urbinp%em_improad(lsmlon*lsmlat), &    
                urbinp%em_perroad(lsmlon*lsmlat), &    
                urbinp%em_wall(lsmlon*lsmlat), &    
                urbinp%alb_roof_dir(lsmlon*lsmlat,numrad), &    
                urbinp%alb_roof_dif(lsmlon*lsmlat,numrad), &    
                urbinp%alb_improad_dir(lsmlon*lsmlat,numrad), &    
                urbinp%alb_perroad_dir(lsmlon*lsmlat,numrad), &    
                urbinp%alb_improad_dif(lsmlon*lsmlat,numrad), &    
                urbinp%alb_perroad_dif(lsmlon*lsmlat,numrad), &    
                urbinp%alb_wall_dir(lsmlon*lsmlat,numrad), &    
                urbinp%alb_wall_dif(lsmlon*lsmlat,numrad), &
                urbinp%ht_roof(lsmlon*lsmlat), &
                urbinp%wind_hgt_canyon(lsmlon*lsmlat), &
                urbinp%tk_wall(lsmlon*lsmlat,nlevsoi), &
                urbinp%tk_roof(lsmlon*lsmlat,nlevsoi), &
                urbinp%tk_improad(lsmlon*lsmlat,5), &
                urbinp%cv_wall(lsmlon*lsmlat,nlevsoi), &
                urbinp%cv_roof(lsmlon*lsmlat,nlevsoi), &
                urbinp%cv_improad(lsmlon*lsmlat,5), &
                urbinp%sandfrac_road(lsmlon*lsmlat,nlevsoi), &
                urbinp%clayfrac_road(lsmlon*lsmlat,nlevsoi), &
                urbinp%scalez_wall(lsmlon*lsmlat), &
                urbinp%scalez_roof(lsmlon*lsmlat), &
                urbinp%thick_wall(lsmlon*lsmlat), &
                urbinp%thick_roof(lsmlon*lsmlat), &
                stat=ier)
       if (ier /= 0) then
          write(iulog,*)'initUrbanInput: allocation error '; call endrun()
       endif

       ! Read urban data (currently assume ascii format)
       
       if (masterproc) write(iulog,*)' Reading in urban input data ...'
       n = getavu()
       call getfil (furbinp, locfn, 0)
       call opnfil (locfn, n, 'f')
       nw = 0
       do j = 1,lsmlat
          do i = 1,lsmlon
             nw = nw + 1
             read (n,*,iostat=ier) desc, urbinp%canyon_hwr(nw)
             if ( index( desc, "canyon_hwr" ) /= 1 )then
                call endrun( 'ERROR:: reading of ASCII urban input file got off track' )
             end if
             if ( (urbinp%canyon_hwr(nw) <= 0.0_r8) .or. (urbinp%canyon_hwr(nw) >= 1000.0_r8) )then
                call endrun( 'ERROR:: canyon_hwr out of range' )
             end if
             read (n,*,iostat=ier) desc, urbinp%wtlunit_roof(nw)
             if ( (urbinp%wtlunit_roof(nw) <= 0.0_r8) .or. (urbinp%wtlunit_roof(nw) >= 1.0_r8) )then
                call endrun( 'ERROR:: wtlunit_roof out of range' )
             end if
             read (n,*,iostat=ier) desc, urbinp%wtroad_perv(nw)
             if ( (urbinp%wtroad_perv(nw) < 0.0_r8) .or. (urbinp%wtroad_perv(nw) > 1.0_r8) )then
                call endrun( 'ERROR:: wtlunit_roof out of range' )
             end if
             read (n,*,iostat=ier) desc, urbinp%em_roof(nw)
             read (n,*,iostat=ier) desc, urbinp%em_improad(nw)
             read (n,*,iostat=ier) desc, urbinp%em_perroad(nw)
             read (n,*,iostat=ier) desc, urbinp%em_wall(nw)
             read (n,*,iostat=ier) desc, urbinp%alb_roof_dir(nw,1), urbinp%alb_roof_dir(nw,2)
             read (n,*,iostat=ier) desc, urbinp%alb_roof_dif(nw,1), urbinp%alb_roof_dif(nw,2)
             read (n,*,iostat=ier) desc, urbinp%alb_improad_dir(nw,1), urbinp%alb_improad_dir(nw,2)
             if ( index( desc, "alb_improad_dir" ) /= 1 )then
                call endrun( 'ERROR:: reading of ASCII urban input file got off track' )
             end if
             read (n,*,iostat=ier) desc, urbinp%alb_improad_dif(nw,1), urbinp%alb_improad_dif(nw,2)
             read (n,*,iostat=ier) desc, urbinp%alb_perroad_dir(nw,1), urbinp%alb_perroad_dir(nw,2)
             read (n,*,iostat=ier) desc, urbinp%alb_perroad_dif(nw,1), urbinp%alb_perroad_dif(nw,2)
             read (n,*,iostat=ier) desc, urbinp%alb_wall_dir(nw,1), urbinp%alb_wall_dir(nw,2)
             read (n,*,iostat=ier) desc, urbinp%alb_wall_dif(nw,1), urbinp%alb_wall_dif(nw,2)
             read (n,*,iostat=ier) desc, urbinp%ht_roof(nw)
             if ( (urbinp%ht_roof(nw) <= 0.0_r8) .or. (urbinp%ht_roof(nw) > 20.0_r8) )then
                call endrun( 'ERROR:: roof height out of range' )
             end if
             if ( index( desc, "ht_roof" ) /= 1 )then
                call endrun( 'ERROR:: reading of ASCII urban input file got off track' )
             end if
             read (n,*,iostat=ier) desc, urbinp%wind_hgt_canyon(nw)
             read (n,*,iostat=ier) desc, (urbinp%tk_wall(nw,k),k=1,nlevsoi)
             read (n,*,iostat=ier) desc, (urbinp%tk_roof(nw,k),k=1,nlevsoi)
             read (n,*,iostat=ier) desc, (urbinp%tk_improad(nw,k),k=1,5)
             read (n,*,iostat=ier) desc, (urbinp%cv_wall(nw,k),k=1,nlevsoi)
             read (n,*,iostat=ier) desc, (urbinp%cv_roof(nw,k),k=1,nlevsoi)
             read (n,*,iostat=ier) desc, (urbinp%cv_improad(nw,k),k=1,5)
             if ( index( desc, "cv_improad" ) /= 1 )then
                call endrun( 'ERROR:: reading of ASCII urban input file got off track' )
             end if
             read (n,*,iostat=ier) desc, (urbinp%sandfrac_road(nw,k),k=1,nlevsoi)
             read (n,*,iostat=ier) desc, (urbinp%clayfrac_road(nw,k),k=1,nlevsoi)
             do k = 1, nlevsoi
                soilfrac = urbinp%sandfrac_road(nw,k) + urbinp%clayfrac_road(nw,k)
                if ( (soilfrac < 0.0_r8) .or. (soilfrac >= 100.0_r8) )then
                   call endrun( 'ERROR:: sandfrac_road and/or clayfrac_road out of range' )
                end if
             end do
             read (n,*,iostat=ier) desc, urbinp%scalez_wall(nw)
             read (n,*,iostat=ier) desc, urbinp%scalez_roof(nw)
             read (n,*,iostat=ier) desc, urbinp%thick_wall(nw)
             read (n,*,iostat=ier) desc, urbinp%thick_roof(nw)
             if ( index( desc, "thick_roof" ) /= 1 )then
                call endrun( 'ERROR:: reading of ASCII urban input file got off track' )
             end if
          end do
       end do
       
       call relavu (n)
       if (masterproc) write(iulog,*)' Sucessfully read urban input data' 

       ! Write out input info

       if (masterproc) then
          nw = 0
          do j = 1,lsmlat
             do i = 1,lsmlon
                nw = nw + 1
                write(iulog,*)'nw= ', nw
                write(iulog,*)'   canyon_hwr     = ', urbinp%canyon_hwr(nw)
                write(iulog,*)'   wtlunit_roof   = ', urbinp%wtlunit_roof(nw)
                write(iulog,*)'   wtroad_perv    = ', urbinp%wtroad_perv(nw)
                write(iulog,*)'   em_roof        = ', urbinp%em_roof(nw)
                write(iulog,*)'   em_improad     = ', urbinp%em_improad(nw)
                write(iulog,*)'   em_perroad     = ', urbinp%em_perroad(nw)
                write(iulog,*)'   em_wall        = ', urbinp%em_wall(nw)
                write(iulog,*)'   alb_roof_dir   = ', urbinp%alb_roof_dir(nw,1), urbinp%alb_roof_dir(nw,2)
                write(iulog,*)'   alb_roof_dif   = ', urbinp%alb_roof_dif(nw,1), urbinp%alb_roof_dif(nw,2)
                write(iulog,*)'   alb_improad_dir= ', urbinp%alb_improad_dir(nw,1), urbinp%alb_improad_dir(nw,2)
                write(iulog,*)'   alb_improad_dif= ', urbinp%alb_improad_dif(nw,1), urbinp%alb_improad_dif(nw,2)
                write(iulog,*)'   alb_perroad_dir= ', urbinp%alb_perroad_dir(nw,1), urbinp%alb_perroad_dir(nw,2)
                write(iulog,*)'   alb_perroad_dif= ', urbinp%alb_perroad_dif(nw,1), urbinp%alb_perroad_dif(nw,2)
                write(iulog,*)'   alb_wall_dir   = ', urbinp%alb_wall_dir(nw,1), urbinp%alb_wall_dir(nw,2)
                write(iulog,*)'   alb_wall_dif   = ', urbinp%alb_wall_dif(nw,1), urbinp%alb_wall_dif(nw,2)
                write(iulog,*)'   ht_roof        = ', urbinp%ht_roof(nw)
                write(iulog,*)'   wind_hgt_canyon= ', urbinp%wind_hgt_canyon(nw)
                write(iulog,*)'   tk_wall        = ', urbinp%tk_wall(nw,:)
                write(iulog,*)'   tk_roof        = ', urbinp%tk_roof(nw,:)
                write(iulog,*)'   tk_improad     = ', urbinp%tk_improad(nw,:)
                write(iulog,*)'   cv_wall        = ', urbinp%cv_wall(nw,:)
                write(iulog,*)'   cv_roof        = ', urbinp%cv_roof(nw,:)
                write(iulog,*)'   cv_improad     = ', urbinp%cv_improad(nw,:)
                write(iulog,*)'   sandfrac_road  = ', urbinp%sandfrac_road(nw,:)
                write(iulog,*)'   clayfrac_road  = ', urbinp%clayfrac_road(nw,:)
                write(iulog,*)'   scalez_wall    = ', urbinp%scalez_wall(nw)
                write(iulog,*)'   scalez_roof    = ', urbinp%scalez_roof(nw)
                write(iulog,*)'   thick_wall    = ', urbinp%thick_wall(nw)
                write(iulog,*)'   thick_roof    = ', urbinp%thick_roof(nw)
             end do
          end do
       end if

    else if (mode == 'finalize') then

       deallocate(urbinp%canyon_hwr, &
                  urbinp%wtlunit_roof, &
                  urbinp%wtroad_perv, &
                  urbinp%em_roof, &
                  urbinp%em_improad, &
                  urbinp%em_perroad, &
                  urbinp%em_wall, &
                  urbinp%alb_roof_dir, &
                  urbinp%alb_roof_dif, &
                  urbinp%alb_improad_dir, &
                  urbinp%alb_perroad_dir, &
                  urbinp%alb_improad_dif, &
                  urbinp%alb_perroad_dif, &
                  urbinp%alb_wall_dir, &
                  urbinp%alb_wall_dif, &
                  urbinp%ht_roof, &
                  urbinp%wind_hgt_canyon, &
                  urbinp%tk_wall, &
                  urbinp%tk_roof, &
                  urbinp%tk_improad, &
                  urbinp%cv_wall, &
                  urbinp%cv_roof, &
                  urbinp%cv_improad, &
                  urbinp%sandfrac_road, &
                  urbinp%clayfrac_road, &
                  urbinp%scalez_wall, &
                  urbinp%scalez_roof, &
                  urbinp%thick_wall, &
                  urbinp%thick_roof, &
                  stat=ier)
       if (ier /= 0) then
          write(iulog,*)'initUrbanInput: deallocation error '; call endrun()
       endif

    else
       write(iulog,*)'initUrbanInput error: mode ',trim(mode),' not supported '
       call endrun()
    end if

  end subroutine UrbanInput

end module UrbanInputMod

