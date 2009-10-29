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
     real(r8), pointer :: thick_wall(:)
     real(r8), pointer :: thick_roof(:)
     integer,  pointer :: nlev_improad(:)
     real(r8), pointer :: t_building_min(:)
     real(r8), pointer :: t_building_max(:)
  end type urbinp_t
  public urbinp_t

  type (urbinp_t)   , public :: urbinp        ! urban input derived type
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
    use clm_varpar, only : lsmlon, lsmlat, numrad, nlevurb, numsolar
    use clm_varctl, only : iulog, fsurdat, single_column
    use fileutils , only : getavu, relavu, getfil, opnfil
    use spmdMod   , only : masterproc
    use ncdio     , only : ncd_iolocal, check_dim, check_ret
    use clmtype   , only : grlnd
    use decompMod , only : get_proc_bounds
    use spmdGathScatMod
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(in) :: mode
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein July 2004
! Revised by Keith Oleson for netcdf input Jan 2008
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=256) :: locfn      ! local file name
    character(len=32)  :: desc
    integer :: ncid,dimid,varid      ! netCDF id's
    integer :: begg,endg             ! start/stop gridcells
    integer :: start4d(4),count4d(4) ! netcdf start/count arrays
    integer :: start3d(3),count3d(3) ! netcdf start/count arrays
    integer :: nw,n,k,i,j,nn,mm
    integer :: ier  
    integer :: nlevurb_i             ! input grid: number of urban vertical levels
    integer :: numsolar_i            ! input grid: number of solar type (DIR/DIF)
    integer :: numrad_i              ! input grid: number of solar bands (VIS/NIR)
    integer :: ret
    real(r8),pointer :: arrayl(:)    ! generic global array
    character(len=32) :: subname = 'UrbanInput'          ! subroutine name
!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg)

    if (mode == 'initialize') then

       ! Allocate dynamic memory

       allocate(urbinp%canyon_hwr(begg:endg), &  
                urbinp%wtlunit_roof(begg:endg), &  
                urbinp%wtroad_perv(begg:endg), &
                urbinp%em_roof(begg:endg), &     
                urbinp%em_improad(begg:endg), &    
                urbinp%em_perroad(begg:endg), &    
                urbinp%em_wall(begg:endg), &    
                urbinp%alb_roof_dir(begg:endg,numrad), &    
                urbinp%alb_roof_dif(begg:endg,numrad), &    
                urbinp%alb_improad_dir(begg:endg,numrad), &    
                urbinp%alb_perroad_dir(begg:endg,numrad), &    
                urbinp%alb_improad_dif(begg:endg,numrad), &    
                urbinp%alb_perroad_dif(begg:endg,numrad), &    
                urbinp%alb_wall_dir(begg:endg,numrad), &    
                urbinp%alb_wall_dif(begg:endg,numrad), &
                urbinp%ht_roof(begg:endg), &
                urbinp%wind_hgt_canyon(begg:endg), &
                urbinp%tk_wall(begg:endg,nlevurb), &
                urbinp%tk_roof(begg:endg,nlevurb), &
                urbinp%tk_improad(begg:endg,nlevurb), &
                urbinp%cv_wall(begg:endg,nlevurb), &
                urbinp%cv_roof(begg:endg,nlevurb), &
                urbinp%cv_improad(begg:endg,nlevurb), &
                urbinp%thick_wall(begg:endg), &
                urbinp%thick_roof(begg:endg), &
                urbinp%nlev_improad(begg:endg), &
                urbinp%t_building_min(begg:endg), &
                urbinp%t_building_max(begg:endg), &
                stat=ier)
       if (ier /= 0) then
          write(iulog,*)'initUrbanInput: allocation error '; call endrun()
       endif

       ! Read urban data
       
       if (masterproc) then
          write(iulog,*)' Reading in urban input data from fsurdat file ...'
          call getfil (fsurdat, locfn, 0)
          call check_ret(nf_open(locfn, 0, ncid), subname)
          write(iulog,*) subname,trim(fsurdat)
          write(iulog,*) " Expected dimensions: lsmlon=",lsmlon," lsmlat=",lsmlat
          if (.not. single_column) then
             call check_dim(ncid, 'lsmlon', lsmlon)
             call check_dim(ncid, 'lsmlon', lsmlon)
          end if 

          call check_ret(nf_inq_dimid(ncid, 'nlevurb', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, nlevurb_i), subname)
          if (nlevurb_i /= nlevurb) then
             write(iulog,*)trim(subname)// ': parameter nlevurb= ',nlevurb, &
                           'does not equal input dataset nlevurb= ',nlevurb_i
             call endrun
          endif
          call check_ret(nf_inq_dimid(ncid, 'numsolar', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, numsolar_i), subname)
          if (numsolar_i /= numsolar) then
             write(iulog,*)trim(subname)// ': parameter numsolar= ',numsolar, &
                           'does not equal input dataset numsolar= ',numsolar_i
             call endrun
          endif
          call check_ret(nf_inq_dimid(ncid, 'numrad', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, numrad_i), subname)
          if (numrad_i /= numrad) then
             write(iulog,*)trim(subname)// ': parameter numrad= ',numrad, &
                           'does not equal input dataset numrad= ',numrad_i
             call endrun
          endif

       end if

       allocate(arrayl(begg:endg))

       call ncd_iolocal(ncid,'CANYON_HWR','read',urbinp%canyon_hwr,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: CANYON_HWR NOT on fsurdat file' )

       call ncd_iolocal(ncid,'WTLUNIT_ROOF','read',urbinp%wtlunit_roof,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: WTLUNIT_ROOF NOT on fsurdat file' )

       call ncd_iolocal(ncid,'WTROAD_PERV','read',urbinp%wtroad_perv,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: WTROAD_PERV NOT on fsurdat file' )

       call ncd_iolocal(ncid,'EM_ROOF','read',urbinp%em_roof,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: EM_ROOF NOT on fsurdat file' )

       call ncd_iolocal(ncid,'EM_IMPROAD','read',urbinp%em_improad,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: EM_IMPROAD NOT on fsurdat file' )

       call ncd_iolocal(ncid,'EM_PERROAD','read',urbinp%em_perroad,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: EM_PERROAD NOT on fsurdat file' )

       call ncd_iolocal(ncid,'EM_WALL','read',urbinp%em_wall,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: EM_WALL NOT on fsurdat file' )

       call ncd_iolocal(ncid,'HT_ROOF','read',urbinp%ht_roof,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: HT_ROOF NOT on fsurdat file' )

       call ncd_iolocal(ncid,'WIND_HGT_CANYON','read',urbinp%wind_hgt_canyon,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: WIND_HGT_CANYON NOT on fsurdat file' )

       call ncd_iolocal(ncid,'THICK_WALL','read',urbinp%thick_wall,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: THICK_WALL NOT on fsurdat file' )

       call ncd_iolocal(ncid,'THICK_ROOF','read',urbinp%thick_roof,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: THICK_ROOF NOT on fsurdat file' )

       call ncd_iolocal(ncid,'NLEV_IMPROAD','read',urbinp%nlev_improad,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: NLEV_IMPROAD NOT on fsurdat file' )

       call ncd_iolocal(ncid,'T_BUILDING_MIN','read',urbinp%t_building_min,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: T_BUILDING_MIN NOT on fsurdat file' )

       call ncd_iolocal(ncid,'T_BUILDING_MAX','read',urbinp%t_building_max,grlnd,status=ret)
       if (ret /= 0) call endrun( trim(subname)//' ERROR: T_BUILDING_MAX NOT on fsurdat file' )

       start4d(1) = 1
       count4d(1) = lsmlon
       start4d(2) = 1
       count4d(2) = lsmlat
       start4d(3) = 1
       count4d(3) = 1
       start4d(4) = 1
       count4d(4) = 1

       do mm = 1,numsolar
         do nn = 1,numrad
            start4d(3) = nn
            start4d(4) = mm
            call ncd_iolocal(ncid,'ALB_IMPROAD','read',arrayl,grlnd,start4d,count4d,status=ret)
            if (ret /= 0) call endrun( trim(subname)//' ERROR: ALB_IMPROAD NOT on fsurdat file' )
            if (mm .eq. 1) then
              urbinp%alb_improad_dir(begg:endg,nn) = arrayl(begg:endg)
            else
              urbinp%alb_improad_dif(begg:endg,nn) = arrayl(begg:endg)
            end if
            call ncd_iolocal(ncid,'ALB_PERROAD','read',arrayl,grlnd,start4d,count4d,status=ret)
            if (ret /= 0) call endrun( trim(subname)//' ERROR: ALB_PERROAD NOT on fsurdat file' )
            if (mm .eq. 1) then
              urbinp%alb_perroad_dir(begg:endg,nn) = arrayl(begg:endg)
            else
              urbinp%alb_perroad_dif(begg:endg,nn) = arrayl(begg:endg)
            end if
            call ncd_iolocal(ncid,'ALB_ROOF','read',arrayl,grlnd,start4d,count4d,status=ret)
            if (ret /= 0) call endrun( trim(subname)//' ERROR: ALB_ROOF NOT on fsurdat file' )
            if (mm .eq. 1) then
              urbinp%alb_roof_dir(begg:endg,nn) = arrayl(begg:endg)
            else
              urbinp%alb_roof_dif(begg:endg,nn) = arrayl(begg:endg)
            end if
            call ncd_iolocal(ncid,'ALB_WALL','read',arrayl,grlnd,start4d,count4d,status=ret)
            if (ret /= 0) call endrun( trim(subname)//' ERROR: ALB_WALL NOT on fsurdat file' )
            if (mm .eq. 1) then
              urbinp%alb_wall_dir(begg:endg,nn) = arrayl(begg:endg)
            else
              urbinp%alb_wall_dif(begg:endg,nn) = arrayl(begg:endg)
            end if
         end do
       end do

       start3d(1) = 1
       count3d(1) = lsmlon
       start3d(2) = 1
       count3d(2) = lsmlat
       start3d(3) = 1
       count3d(3) = 1
       do nn = 1,nlevurb
          start3d(3) = nn
          call ncd_iolocal(ncid,'TK_IMPROAD','read',arrayl,grlnd,start3d,count3d,status=ret)
          if (ret /= 0) call endrun( trim(subname)//' ERROR: TK_IMPROAD NOT on fsurdat file' )
          urbinp%tk_improad(begg:endg,nn) = arrayl(begg:endg)
          call ncd_iolocal(ncid,'TK_ROOF','read',arrayl,grlnd,start3d,count3d,status=ret)
          if (ret /= 0) call endrun( trim(subname)//' ERROR: TK_ROOF NOT on fsurdat file' )
          urbinp%tk_roof(begg:endg,nn) = arrayl(begg:endg)
          call ncd_iolocal(ncid,'TK_WALL','read',arrayl,grlnd,start3d,count3d,status=ret)
          if (ret /= 0) call endrun( trim(subname)//' ERROR: TK_WALL NOT on fsurdat file' )
          urbinp%tk_wall(begg:endg,nn) = arrayl(begg:endg)
          call ncd_iolocal(ncid,'CV_IMPROAD','read',arrayl,grlnd,start3d,count3d,status=ret)
          if (ret /= 0) call endrun( trim(subname)//' ERROR: CV_IMPROAD NOT on fsurdat file' )
          urbinp%cv_improad(begg:endg,nn) = arrayl(begg:endg)
          call ncd_iolocal(ncid,'CV_ROOF','read',arrayl,grlnd,start3d,count3d,status=ret)
          if (ret /= 0) call endrun( trim(subname)//' ERROR: CV_ROOF NOT on fsurdat file' )
          urbinp%cv_roof(begg:endg,nn) = arrayl(begg:endg)
          call ncd_iolocal(ncid,'CV_WALL','read',arrayl,grlnd,start3d,count3d,status=ret)
          if (ret /= 0) call endrun( trim(subname)//' ERROR: CV_WALL NOT on fsurdat file' )
          urbinp%cv_wall(begg:endg,nn) = arrayl(begg:endg)
       end do

       deallocate(arrayl)
       
       if (masterproc) then  
          call check_ret(nf_close(ncid), subname)
          write(iulog,*)' Sucessfully read urban input data' 
          write(iulog,*)
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
                  urbinp%thick_wall, &
                  urbinp%thick_roof, &
                  urbinp%nlev_improad, &
                  urbinp%t_building_min, &
                  urbinp%t_building_max, &
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

