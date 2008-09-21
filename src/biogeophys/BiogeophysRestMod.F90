#include <misc.h>
#include <preproc.h>

module BiogeophysRestMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: BiogeophysRestMod
!
! !DESCRIPTION:
! Reads from or biogeophysics restart/initial data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
!
! !PUBLIC TYPES:
  implicit none
! save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: BiogeophysRest
!
! !REVISION HISTORY:
! 2005-06-12: Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BiogeophysRest
!
! !INTERFACE:
  subroutine BiogeophysRest( ncid, flag )
!
! !DESCRIPTION:
! Read/Write biogeophysics information to/from restart file.
!
! !USES:
    use clmtype
    use ncdio
    use decompMod       , only : get_proc_bounds
    use clm_varpar      , only : nlevsoi, nlevsno
    use clm_varcon      , only : denice, denh2o
    use clm_varctl      , only : allocate_all_vegpfts, nsrest
    use initSurfAlbMod  , only : do_initsurfalb
    use clm_time_manager, only : is_first_step
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid           ! netcdf id
    character(len=*), intent(in) :: flag  ! 'read' or 'write'
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
! 12/11/2003, Peter Thornton: Added cps%coszen, pps%gdir, and pps%omega
!   for new sunlit/shaded canopy algorithm (in SUNSHA ifdef block)
! 4/25/2005, Peter Thornton: Removed the SUNSHA ifdefs, since this is now the
!   default code behavior.
! 6/12/2005, Moved to netcdf format and renamed file
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: p,c,l,g,j    ! indices
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    logical :: readvar      ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Note - for the snow interfaces, are only examing the snow interfaces
    ! above zi=0 which is why zisno and zsno have the same level dimension below
    ! (Note - for zisno, zi(0) is set to 0 in routine iniTimeConst)
    
    ! pft weight wrt gridcell 

    if (allocate_all_vegpfts) then
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='PFT_WTGCELL', xtype=nf_double,  &
               dim1name='pft', &
               long_name='pft weight relative to corresponding gridcell', units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_iolocal(varname='PFT_WTGCELL', data=pptr%wtgcell, &
               dim1name=namep, &
               ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun()
          end if
       end if
    end if

    ! pft weight wrt landunit

    if (allocate_all_vegpfts) then
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='PFT_WTLUNIT', xtype=nf_double,  &
               dim1name='pft', &
               long_name='pft weight relative to corresponding landunit', units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_iolocal(varname='PFT_WTLUNIT', data=pptr%wtlunit, &
               dim1name=namep, &
               ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun()
          end if
       end if
    end if

    ! pft weight wrt column

    if (allocate_all_vegpfts) then
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='PFT_WTCOL', xtype=nf_double,  &
               dim1name='pft', &
               long_name='pft weight relative to corresponding column', units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_iolocal(varname='PFT_WTCOL', data=pptr%wtcol, &
               dim1name=namep, &
               ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) call endrun()
          end if
       end if
    end if

    ! pft energy flux - eflx_lwrad_out

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='EFLX_LWRAD_OUT', xtype=nf_double,  &
            dim1name='pft', &
            long_name='emitted infrared (longwave) radiation', units='watt/m^2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='EFLX_LWRAD_OUT', data=pptr%pef%eflx_lwrad_out, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - snow levels

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='SNLSNO', xtype=nf_int,  &
            dim1name='column', &
            long_name='number of snow layers', units='unitless')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='SNLSNO', data=cptr%cps%snl, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - snowdp

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='SNOWDP', xtype=nf_double, &
            dim1name='column', &
            long_name='snow depth', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='SNOWDP', data=cptr%cps%snowdp, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - snowage

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='SNOWAGE', xtype=nf_double,  &
            dim1name='column', &
            long_name='snow age', units='unitless')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='SNOWAGE', data=cptr%cps%snowage, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - wa

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='WA', xtype=nf_double,  &
            dim1name='column', &
            long_name='water in the unconfined aquifer', units='mm')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='WA', data=cptr%cws%wa, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - wt

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='WT', xtype=nf_double,  &
            dim1name='column', &
            long_name='total water storage', units='mm')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='WT', data=cptr%cws%wt, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - zwt

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ZWT', xtype=nf_double,  &
            dim1name='column', &
            long_name='water table depth', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ZWT', data=cptr%cws%zwt, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type physical state variable - frac_sno

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frac_sno', xtype=nf_double,  &
            dim1name='column',&
            long_name='fraction of ground covered by snow (0 to 1)',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frac_sno', data=cptr%cps%frac_sno, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type physical state variable - dzsno

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='DZSNO', xtype=nf_double,  &
            dim1name='column', dim2name='levsno', &
            long_name='snow layer thickness', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='DZSNO', data=cptr%cps%dz, &
            dim1name=namec, dim2name='levsno', &
            lowerb2=-nlevsno+1, upperb2=0, ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type physical state variable - zsno

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ZSNO', xtype=nf_double,  &
            dim1name='column', dim2name='levsno',&
            long_name='snow layer depth', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ZSNO', data=cptr%cps%z, &
            dim1name=namec, dim2name='levsno', lowerb2=-nlevsno+1, upperb2=0, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type physical state variable - zisno

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ZISNO', xtype=nf_double,  &
            dim1name='column', dim2name='levsno', &
            long_name='snow interface depth', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ZISNO', data=cptr%cps%zi, &
            dim1name=namec, dim2name='levsno', lowerb2=-nlevsno, upperb2=-1, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type physical state variable - coszen

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='coszen', xtype=nf_double,  &
            dim1name='column', &
            long_name='cosine of solar zenith angle', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='coszen', data=cptr%cps%coszen, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - gdir

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gdir', xtype=nf_double,  &
            dim1name='pft', &
            long_name='leaf projection in solar direction (0 to 1)', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gdir', data=pptr%pps%gdir, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - omega

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='omega', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='fraction of intercepted radiation that is scattered (0 to 1)', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='omega', data=pptr%pps%omega, &
            dim1name=namep, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_roof_dir

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_roof_dir', xtype=nf_double,  &
            dim1name='landunit', dim2name='numrad', &
            long_name='direct solar absorbed by roof per unit ground area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sabs_roof_dir', data=lptr%lps%sabs_roof_dir, &
            dim1name=namel, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_roof_dif

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_roof_dif', xtype=nf_double,  &
            dim1name='landunit', dim2name='numrad', &
            long_name='diffuse solar absorbed by roof per unit ground area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sabs_roof_dif', data=lptr%lps%sabs_roof_dif, &
            dim1name=namel, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_sunwall_dir

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_sunwall_dir', xtype=nf_double,  &
            dim1name='landunit', dim2name='numrad', &
            long_name='direct solar absorbed by sunwall per unit wall area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sabs_sunwall_dir', data=lptr%lps%sabs_sunwall_dir, &
            dim1name=namel, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_sunwall_dif

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_sunwall_dif', xtype=nf_double,  &
            dim1name='landunit', dim2name='numrad', &
            long_name='diffuse solar absorbed by sunwall per unit wall area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sabs_sunwall_dif', data=lptr%lps%sabs_sunwall_dif, &
            dim1name=namel, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_shadewall_dir

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_shadewall_dir', xtype=nf_double,  &
            dim1name='landunit', dim2name='numrad', &
            long_name='direct solar absorbed by shadewall per unit wall area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sabs_shadewall_dir', data=lptr%lps%sabs_shadewall_dir, &
            dim1name=namel, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_shadewall_dif

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_shadewall_dif', xtype=nf_double,  &
            dim1name='landunit', dim2name='numrad', &
            long_name='diffuse solar absorbed by shadewall per unit wall area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sabs_shadewall_dif', data=lptr%lps%sabs_shadewall_dif, &
            dim1name=namel, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_improad_dir

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_improad_dir', xtype=nf_double,  &
            dim1name='landunit', dim2name='numrad', &
            long_name='direct solar absorbed by impervious road per unit ground area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sabs_improad_dir', data=lptr%lps%sabs_improad_dir, &
            dim1name=namel, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_improad_dif

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_improad_dif', xtype=nf_double,  &
            dim1name='landunit', dim2name='numrad', &
            long_name='diffuse solar absorbed by impervious road per unit ground area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sabs_improad_dif', data=lptr%lps%sabs_improad_dif, &
            dim1name=namel, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_perroad_dir

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_perroad_dir', xtype=nf_double,  &
            dim1name='landunit', dim2name='numrad', &
            long_name='direct solar absorbed by pervious road per unit ground area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sabs_perroad_dir', data=lptr%lps%sabs_perroad_dir, &
            dim1name=namel, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - sabs_perroad_dif

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sabs_perroad_dif', xtype=nf_double,  &
            dim1name='landunit', dim2name='numrad', &
            long_name='diffuse solar absorbed by pervious road per unit ground area per unit incident flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sabs_perroad_dif', data=lptr%lps%sabs_perroad_dif, &
            dim1name=namel, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - vf_sr

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vf_sr', xtype=nf_double,  &
            dim1name='landunit', &
            long_name='view factor of sky for road',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='vf_sr', data=lptr%lps%vf_sr, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - vf_wr

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vf_wr', xtype=nf_double,  &
            dim1name='landunit', &
            long_name='view factor of one wall for road',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='vf_wr', data=lptr%lps%vf_wr, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - vf_sw

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vf_sw', xtype=nf_double,  &
            dim1name='landunit', &
            long_name='view factor of sky for one wall',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='vf_sw', data=lptr%lps%vf_sw, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - vf_rw

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vf_rw', xtype=nf_double,  &
            dim1name='landunit', &
            long_name='view factor of road for one wall',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='vf_rw', data=lptr%lps%vf_rw, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - vf_ww

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vf_ww', xtype=nf_double,  &
            dim1name='landunit', &
            long_name='view factor of opposing wall for one wall',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='vf_ww', data=lptr%lps%vf_ww, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - taf

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='taf', xtype=nf_double,  &
            dim1name='landunit', &
            long_name='urban canopy air temperature',units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='taf', data=lptr%lps%taf, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! landunit type physical state variable - qaf

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='qaf', xtype=nf_double,  &
            dim1name='landunit', &
            long_name='urban canopy specific humidity',units='kg/kg')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='qaf', data=lptr%lps%qaf, &
            dim1name=namel, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - albd

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='albd', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='surface albedo (direct) (0 to 1)',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='albd', data=pptr%pps%albd, &
            dim1name=namep, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             if (nsrest == 0) do_initsurfalb = .true.
          end if
       end if
    end if

    ! pft type physical state variable - albi

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='albi', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='surface albedo (diffuse) (0 to 1)',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='albi', data=pptr%pps%albi, &
            dim1name=namep, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          else
             if (nsrest == 0) do_initsurfalb = .true.
          end if
       end if
    end if

    ! column type physical state variable - albgrd

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='albgrd', xtype=nf_double,  &
            dim1name='column', dim2name='numrad', &
            long_name='ground albedo (direct) (0 to 1)',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='albgrd', data=cptr%cps%albgrd, &
            dim1name=namec, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type physical state variable - albgri

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='albgri', xtype=nf_double,  &
            dim1name='column', dim2name='numrad', &
            long_name='ground albedo (indirect) (0 to 1)',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='albgri', data=cptr%cps%albgri, &
            dim1name=namec, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! column water state variable - h2osno

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='H2OSNO', xtype=nf_double,  &
            dim1name='column', &
            long_name='snow water', units='kg/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='H2OSNO', data=cptr%cws%h2osno, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - h2osoi_liq

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='H2OSOI_LIQ', xtype=nf_double,  &
            dim1name='column', dim2name='levtot', &
            long_name='liquid water', units='kg/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='H2OSOI_LIQ', data=cptr%cws%h2osoi_liq, &
            dim1name=namec, dim2name='levtot', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column water state variable - h2osoi_ice

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='H2OSOI_ICE', xtype=nf_double,   &
            dim1name='column', dim2name='levtot', long_name='ice lens', units='kg/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='H2OSOI_ICE', data=cptr%cws%h2osoi_ice, &
            dim1name=namec, dim2name='levtot', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! column energy state variable - t_grnd

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_GRND', xtype=nf_double,  &
            dim1name='column', &
            long_name='ground temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_GRND', data=cptr%ces%t_grnd, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! column urban energy state variable - eflx_urban_ac

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='URBAN_AC', xtype=nf_double,  &
            dim1name='column', &
            long_name='urban air conditioning flux', units='watt/m^2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='URBAN_AC', data=cptr%cef%eflx_urban_ac, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! column urban energy state variable - eflx_urban_heat

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='URBAN_HEAT', xtype=nf_double,  &
            dim1name='column', &
            long_name='urban heating flux', units='watt/m^2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='URBAN_HEAT', data=cptr%cef%eflx_urban_heat, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! pft energy state variable - t_ref2m_min

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MIN', xtype=nf_double,  &
            dim1name='pft', &
            long_name='daily minimum of average 2 m height surface air temperature (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_REF2M_MIN', data=pptr%pes%t_ref2m_min, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_max

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MAX', xtype=nf_double,  &
            dim1name='pft', &
            long_name='daily maximum of average 2 m height surface air temperature (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_REF2M_MAX', data=pptr%pes%t_ref2m_max, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_min_inst

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MIN_INST', xtype=nf_double,  &
            dim1name='pft', &
            long_name='instantaneous daily min of average 2 m height surface air temp (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_REF2M_MIN_INST', data=pptr%pes%t_ref2m_min_inst, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_max_inst

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MAX_INST', xtype=nf_double,  &
            dim1name='pft', &
            long_name='instantaneous daily max of average 2 m height surface air temp (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_REF2M_MAX_INST', data=pptr%pes%t_ref2m_max_inst, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

    ! pft energy state variable - t_ref2m_u

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname="T_REF2M_U", xtype=nf_double,  &
            dim1name='pft', &
            long_name='Urban 2m height surface air temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname="T_REF2M_U", data=pptr%pes%t_ref2m_u, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! column energy state variable - t_grnd_u

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_GRND_U', xtype=nf_double,  &
            dim1name='column', &
            long_name='urban ground temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_GRND_U', data=cptr%ces%t_grnd_u, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! pft energy state variable - t_ref2m_min_u

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MIN_U', xtype=nf_double,  &
            dim1name='pft', &
            long_name='urban daily minimum of average 2 m height surface air temperature (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_REF2M_MIN_U', data=pptr%pes%t_ref2m_min_u, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_max_u

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MAX_U', xtype=nf_double,  &
            dim1name='pft', &
            long_name='urban daily maximum of average 2 m height surface air temperature (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_REF2M_MAX_U', data=pptr%pes%t_ref2m_max_u, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_min_inst_u

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MIN_INST_U', xtype=nf_double,  &
            dim1name='pft', &
            long_name='urban instantaneous daily min of average 2 m height surface air temp (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_REF2M_MIN_INST_U', data=pptr%pes%t_ref2m_min_inst_u, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_max_inst_u

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MAX_INST_U', xtype=nf_double,  &
            dim1name='pft', &
            long_name='urban instantaneous daily max of average 2 m height surface air temp (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_REF2M_MAX_INST_U', data=pptr%pes%t_ref2m_max_inst_u, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

    ! pft energy state variable - t_ref2m_r

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname="T_REF2M_R", xtype=nf_double,  &
            dim1name='pft', &
            long_name='Rural 2m height surface air temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname="T_REF2M_R", data=pptr%pes%t_ref2m_r, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! column energy state variable - t_grnd_r

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_GRND_R', xtype=nf_double,  &
            dim1name='column', &
            long_name='rural ground temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_GRND_R', data=cptr%ces%t_grnd_r, &
            dim1name=namec, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

   ! pft energy state variable - t_ref2m_min_r

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MIN_R', xtype=nf_double,  &
            dim1name='pft', &
            long_name='rural daily minimum of average 2 m height surface air temperature (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_REF2M_MIN_R', data=pptr%pes%t_ref2m_min_r, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_max_r

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MAX_R', xtype=nf_double,  &
            dim1name='pft', &
            long_name='rural daily maximum of average 2 m height surface air temperature (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_REF2M_MAX_R', data=pptr%pes%t_ref2m_max_r, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_min_inst_r

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MIN_INST_R', xtype=nf_double,  &
            dim1name='pft', &
            long_name='rural instantaneous daily min of average 2 m height surface air temp (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_REF2M_MIN_INST_R', data=pptr%pes%t_ref2m_min_inst_r, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

   ! pft energy state variable - t_ref2m_max_inst_r

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_REF2M_MAX_INST_R', xtype=nf_double,  &
            dim1name='pft', &
            long_name='rural instantaneous daily max of average 2 m height surface air temp (K)', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_REF2M_MAX_INST_R', data=pptr%pes%t_ref2m_max_inst_r, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    endif

    ! column energy state variable - t_soisno

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_SOISNO', xtype=nf_double,   &
            dim1name='column', dim2name='levtot', &
            long_name='soil-snow temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_SOISNO', data=cptr%ces%t_soisno, &
            dim1name=namec, dim2name='levtot', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column type energy state variable - t_lake

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_LAKE', xtype=nf_double,  &
            dim1name='column', dim2name='levlak', &
            long_name='lake temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_LAKE', data=cptr%ces%t_lake, &
            dim1name=namec, dim2name='levlak', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft physical state variable - frac_veg_nosno_alb

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='FRAC_VEG_NOSNO_ALB', xtype=nf_int,  &
            dim1name='pft',&
            long_name='fraction of vegetation not covered by snow (0 or 1)',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='FRAC_VEG_NOSNO_ALB', data=pptr%pps%frac_veg_nosno_alb, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - fwet

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='FWET', xtype=nf_double,  &
            dim1name='pft', &
            long_name='fraction of canopy that is wet (0 to 1)', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='FWET', data=pptr%pps%fwet, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - tlai

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tlai', xtype=nf_double,  &
            dim1name='pft', &
            long_name='one-sided leaf area index, no burying by snow', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='tlai', data=pptr%pps%tlai, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - tsai

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tsai', xtype=nf_double,  &
            dim1name='pft', &
            long_name='one-sided stem area index, no burying by snow', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='tsai', data=pptr%pps%tsai, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - elai

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='elai', xtype=nf_double,  &
            dim1name='pft', &
            long_name='one-sided leaf area index, with burying by snow', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='elai', data=pptr%pps%elai, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - esai

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='esai', xtype=nf_double,  &
            dim1name='pft', &
            long_name='one-sided stem area index, with burying by snow', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='esai', data=pptr%pps%esai, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - fsun

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fsun', xtype=nf_double,  &
            dim1name='pft', &
            long_name='sunlit fraction of canopy', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='fsun', data=pptr%pps%fsun, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - htop

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='htop', xtype=nf_double,  &
            dim1name='pft', &
            long_name='canopy top', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='htop', data=pptr%pps%htop, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - hbot

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='hbot', xtype=nf_double,  &
            dim1name='pft', &
            long_name='canopy botton', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='hbot', data=pptr%pps%hbot, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - fabd

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fabd', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='flux absorbed by veg per unit direct flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='fabd', data=pptr%pps%fabd, &
            dim1name=namep, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - fabi

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fabi', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='flux absorbed by veg per unit diffuse flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='fabi', data=pptr%pps%fabi, &
            dim1name=namep, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - ftdd

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ftdd', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='down direct flux below veg per unit direct flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ftdd', data=pptr%pps%ftdd, &
            dim1name=namep, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - ftid

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ftid', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='down diffuse flux below veg per unit direct flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ftid', data=pptr%pps%ftid, &
            dim1name=namep, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft type physical state variable - ftii

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ftii', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='down diffuse flux below veg per unit diffuse flux',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ftii', data=pptr%pps%ftii, &
            dim1name=namep, dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft energy state variable - t_veg

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_VEG', xtype=nf_double,  &
            dim1name='pft', &
            long_name='vegetation temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_VEG', data=pptr%pes%t_veg, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft energy state variable - t_ref2m

#if (defined CN)
    varname = 't_ref2m'
#else
    varname = 'T_REF2M'
#endif
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname=varname, xtype=nf_double,  &
            dim1name='pft', &
            long_name='2m height surface air temperature', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname=varname, data=pptr%pes%t_ref2m, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
#if (defined CN)
          call endrun()
#else
          if (allocate_all_vegpfts) then
             call endrun()
          else
             if (is_restart()) call endrun()
          end if
#endif
       end if
    end if

    ! pft type water state variable - h2ocan

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='H2OCAN', xtype=nf_double,  &
            dim1name='pft', &
            long_name='canopy water', units='kg/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='H2OCAN', data=pptr%pws%h2ocan, &
            dim1name=namep, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! ------------------------------------------------------------
    ! Determine volumetric soil water (for read only)
    ! ------------------------------------------------------------

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    if (flag == 'read' ) then
       do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
          do c = begc,endc
             cptr%cws%h2osoi_vol(c,j) = cptr%cws%h2osoi_liq(c,j)/(cptr%cps%dz(c,j)*denh2o) &
                                      + cptr%cws%h2osoi_ice(c,j)/(cptr%cps%dz(c,j)*denice)
          end do
       end do

       ! ------------------------------------------------------------
       ! If initial run -- ensure that water is properly bounded
       ! ------------------------------------------------------------

       if ( is_first_step() )then
          do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
             do c = begc,endc
                ! soil water should be ... less than saturation
                if (      cptr%cws%h2osoi_vol(c,j) > cptr%cps%watsat(c,j) )then
                   cptr%cws%h2osoi_liq(c,j) = (cptr%cps%watsat(c,j) &
                                - cptr%cws%h2osoi_ice(c,j)/(cptr%cps%dz(c,j)*denice)) * &
                                  (cptr%cps%dz(c,j)*denh2o)
                   cptr%cws%h2osoi_vol(c,j) = cptr%cps%watsat(c,j)
                ! ........................ and greater than zero
                else if ( cptr%cws%h2osoi_vol(c,j) < 0.0_r8 )then
                      cptr%cws%h2osoi_vol(c,j) = 0.0_r8
                end if
             end do
          end do
       end if
    endif

  end subroutine BiogeophysRest

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_restart
!
! !INTERFACE:
  logical function is_restart( )
!
! !DESCRIPTION:
! Determine if restart run
!
! !USES:
    use clm_varctl, only : nsrest
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    if (nsrest == 1) then
       is_restart = .true.
    else
       is_restart = .false.
    end if

  end function is_restart

end module BiogeophysRestMod
