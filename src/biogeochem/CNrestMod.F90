#include <misc.h>
#include <preproc.h>

module CNrestMod

#if (defined CN)
!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: CNrestMod
! 
! !DESCRIPTION: 
! Read/Write to/from CN info to CLM restart file. 
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNRest
!
! !REVISION HISTORY:
! 11/05/03: Module created by Peter Thornton
!
!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNRest
!
! !INTERFACE:
  subroutine CNRest ( ncid, flag )
!
! !DESCRIPTION: 
! Read/write CN restart data
!
! !USES:
    use clmtype
    use clm_atmlnd, only : clm_a2l
    use ncdio
    use clm_varpar, only : numrad
    use decompMod , only : get_proc_bounds
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid            ! restart unit 
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: Peter Thornton
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: c,p,j        ! indices 
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices 
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    real(r8):: m            ! multiplier for the exit_spinup code
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

    ! Determine necessary subgrid bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    !--------------------------------
    ! pft physical state variables 
    !--------------------------------
    
    ! laisun
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='laisun', xtype=nf_double,  &
            dim1name='pft',long_name='sunlit projected leaf area index',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='laisun', data=pptr%pps%laisun, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! laisha
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='laisha', xtype=nf_double,  &
            dim1name='pft',long_name='sunlit projected leaf area index',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='laisha', data=pptr%pps%laisha, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! slasun
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='slasun', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='slasun', data=pptr%pps%slasun, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! slasha
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='slasha', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='slasha', data=pptr%pps%slasha, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! lncsun
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='lncsun', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='lncsun', data=pptr%pps%lncsun, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! lncsha
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='lncsha', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='lncsha', data=pptr%pps%lncsha, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! vcmxsun
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vcmxsun', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='vcmxsun', data=pptr%pps%vcmxsun, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! vcmxsha
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vcmxsha', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='vcmxsha', data=pptr%pps%vcmxsha, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cisun
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cisun', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cisun', data=pptr%pps%cisun, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cisha
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cisha', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cisha', data=pptr%pps%cisha, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! alphapsnsun
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='alphapsnsun', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='alphapsnsun', data=pptr%pps%alphapsnsun, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! alphapsnsha
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='alphapsnsha', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='alphapsnsha', data=pptr%pps%alphapsnsha, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! 2-d fields (all numrad)
    
    ! eff_kid (2-d numrad)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='eff_kid', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='effective extinction coefficient for indirect from direct',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='eff_kid', data=pptr%pps%eff_kid, &
            dim1name='pft', dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! eff_kii (2-d numrad)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='eff_kii', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='effective extinction coefficient for indirect from indirect',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='eff_kii', data=pptr%pps%eff_kii, &
            dim1name='pft', dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! sun_faid (2-d numrad)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sun_faid', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='fraction sun canopy absorbed indirect from direct',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sun_faid', data=pptr%pps%sun_faid, &
            dim1name='pft', dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if
    
    ! sun_faii (2-d numrad)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sun_faii', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='fraction sun canopy absorbed indirect from indirect',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sun_faii', data=pptr%pps%sun_faii, &
            dim1name='pft', dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! sha_faid (2-d numrad)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sha_faid', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='fraction shade canopy absorbed indirect from direct',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sha_faid', data=pptr%pps%sha_faid, &
            dim1name='pft', dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! sha_faii (2-d numrad)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sha_faii', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='fraction shade canopy absorbed indirect from direct',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sha_faii', data=pptr%pps%sha_faii, &
            dim1name='pft', dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    !--------------------------------
    ! pft ecophysiological variables 
    !--------------------------------
    
    ! dormant_flag
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='dormant_flag', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='dormant_flag', data=pptr%pepv%dormant_flag, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! days_active
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='days_active', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='days_active', data=pptr%pepv%days_active, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! onset_flag
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='onset_flag', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='onset_flag', data=pptr%pepv%onset_flag, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! onset_counter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='onset_counter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='onset_counter', data=pptr%pepv%onset_counter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! onset_gddflag
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='onset_gddflag', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='onset_gddflag', data=pptr%pepv%onset_gddflag, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! onset_fdd
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='onset_fdd', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='onset_fdd', data=pptr%pepv%onset_fdd, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! onset_gdd
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='onset_gdd', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='onset_gdd', data=pptr%pepv%onset_gdd, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! onset_swi
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='onset_swi', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='onset_swi', data=pptr%pepv%onset_swi, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! offset_flag
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='offset_flag', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='offset_flag', data=pptr%pepv%offset_flag, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! offset_counter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='offset_counter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='offset_counter', data=pptr%pepv%offset_counter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! offset_fdd
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='offset_fdd', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='offset_fdd', data=pptr%pepv%offset_fdd, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! offset_swi
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='offset_swi', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='offset_swi', data=pptr%pepv%offset_swi, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! lgsf
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='lgsf', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='lgsf', data=pptr%pepv%lgsf, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! bglfr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='bglfr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='bglfr', data=pptr%pepv%bglfr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! bgtr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='bgtr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='bgtr', data=pptr%pepv%bgtr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! dayl
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='dayl', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='dayl', data=pptr%pepv%dayl, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! prev_dayl
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prev_dayl', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='prev_dayl', data=pptr%pepv%prev_dayl, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! annavg_t2m
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annavg_t2m', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='annavg_t2m', data=pptr%pepv%annavg_t2m, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! tempavg_t2m
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempavg_t2m', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='tempavg_t2m', data=pptr%pepv%tempavg_t2m, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! gpp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gpp_pepv', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gpp_pepv', data=pptr%pepv%gpp, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! availc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='availc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='availc', data=pptr%pepv%availc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! xsmrpool_recover
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='xsmrpool_recover', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='xsmrpool_recover', data=pptr%pepv%xsmrpool_recover, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! xsmrpool_c13ratio
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='xsmrpool_c13ratio', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='xsmrpool_c13ratio', data=pptr%pepv%xsmrpool_c13ratio, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! alloc_pnow
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='alloc_pnow', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='alloc_pnow', data=pptr%pepv%alloc_pnow, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! c_allometry
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='c_allometry', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='c_allometry', data=pptr%pepv%c_allometry, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! n_allometry
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='n_allometry', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='n_allometry', data=pptr%pepv%n_allometry, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! plant_ndemand
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='plant_ndemand', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='plant_ndemand', data=pptr%pepv%plant_ndemand, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! tempsum_plant_ndemand
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempsum_plant_ndemand', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='tempsum_plant_ndemand', data=pptr%pepv%tempsum_plant_ndemand, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !annsum_plant_ndemand 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_plant_ndemand', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='annsum_plant_ndemand', data=pptr%pepv%annsum_plant_ndemand, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! tempsum_retransn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempsum_retransn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='tempsum_retransn', data=pptr%pepv%tempsum_retransn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! annsum_retransn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_retransn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='annsum_retransn', data=pptr%pepv%annsum_retransn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! avail_retransn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='avail_retransn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='avail_retransn', data=pptr%pepv%avail_retransn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! plant_nalloc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='plant_nalloc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='plant_nalloc', data=pptr%pepv%plant_nalloc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! plant_calloc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='plant_calloc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='plant_calloc', data=pptr%pepv%plant_calloc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! excess_cflux
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='excess_cflux', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='excess_cflux', data=pptr%pepv%excess_cflux, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! downreg
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='downreg', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='downreg', data=pptr%pepv%downreg, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! prev_leafc_to_litter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prev_leafc_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='prev_leafc_to_litter', data=pptr%pepv%prev_leafc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! prev_frootc_to_litter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prev_frootc_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='prev_frootc_to_litter', data=pptr%pepv%prev_frootc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! tempsum_npp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempsum_npp', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='tempsum_npp', data=pptr%pepv%tempsum_npp, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! annsum_npp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_npp', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='annsum_npp', data=pptr%pepv%annsum_npp, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! rc13_canair
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='rc13_canair', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='rc13_canair', data=pptr%pepv%rc13_canair, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! rc13_psnsun
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='rc13_psnsun', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='rc13_psnsun', data=pptr%pepv%rc13_psnsun, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! rc13_psnsha
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='rc13_psnsha', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='rc13_psnsha', data=pptr%pepv%rc13_psnsha, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !--------------------------------
    ! pft carbon state variables 
    !--------------------------------
    
    ! leafc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc', data=pptr%pcs%leafc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_storage', data=pptr%pcs%leafc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_xfer', data=pptr%pcs%leafc_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc', data=pptr%pcs%frootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_storage', data=pptr%pcs%frootc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !frootc_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_xfer', data=pptr%pcs%frootc_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemc', data=pptr%pcs%livestemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemc_storage', data=pptr%pcs%livestemc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemc_xfer', data=pptr%pcs%livestemc_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemc', data=pptr%pcs%deadstemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemc_storage', data=pptr%pcs%deadstemc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemc_xfer', data=pptr%pcs%deadstemc_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootc', data=pptr%pcs%livecrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootc_storage', data=pptr%pcs%livecrootc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootc_xfer', data=pptr%pcs%livecrootc_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootc', data=pptr%pcs%deadcrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootc_storage', data=pptr%pcs%deadcrootc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootc_xfer', data=pptr%pcs%deadcrootc_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! gresp_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gresp_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gresp_storage', data=pptr%pcs%gresp_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! gresp_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gresp_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gresp_xfer', data=pptr%pcs%gresp_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool', data=pptr%pcs%cpool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! xsmrpool
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='xsmrpool', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='xsmrpool', data=pptr%pcs%xsmrpool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft_ctrunc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_ctrunc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='pft_ctrunc', data=pptr%pcs%pft_ctrunc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! dispvegc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='dispvegc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='dispvegc', data=pptr%pcs%dispvegc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! storvegc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='storvegc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='storvegc', data=pptr%pcs%storvegc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! totvegc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totvegc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totvegc', data=pptr%pcs%totvegc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! totpftc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totpftc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totpftc', data=pptr%pcs%totpftc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !--------------------------------
    ! C13 pft carbon state variables 
    !--------------------------------
    
    ! leafc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_13', data=pptr%pc13s%leafc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_storage_13', data=pptr%pc13s%leafc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_xfer_13', data=pptr%pc13s%leafc_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_13', data=pptr%pc13s%frootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_storage_13', data=pptr%pc13s%frootc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !frootc_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_xfer_13', data=pptr%pc13s%frootc_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemc_13', data=pptr%pc13s%livestemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemc_storage_13', data=pptr%pc13s%livestemc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemc_xfer_13', data=pptr%pc13s%livestemc_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemc_13', data=pptr%pc13s%deadstemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemc_storage_13', data=pptr%pc13s%deadstemc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemc_xfer_13', data=pptr%pc13s%deadstemc_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootc_13', data=pptr%pc13s%livecrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootc_storage_13', data=pptr%pc13s%livecrootc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootc_xfer_13', data=pptr%pc13s%livecrootc_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootc_13', data=pptr%pc13s%deadcrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootc_storage_13', data=pptr%pc13s%deadcrootc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootc_xfer_13', data=pptr%pc13s%deadcrootc_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! gresp_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gresp_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gresp_storage_13', data=pptr%pc13s%gresp_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! gresp_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gresp_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gresp_xfer_13', data=pptr%pc13s%gresp_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_13', data=pptr%pc13s%cpool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! xsmrpool
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='xsmrpool_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='xsmrpool_13', data=pptr%pc13s%xsmrpool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft_ctrunc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_ctrunc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='pft_ctrunc_13', data=pptr%pc13s%pft_ctrunc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! dispvegc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='dispvegc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='dispvegc_13', data=pptr%pc13s%dispvegc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! storvegc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='storvegc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='storvegc_13', data=pptr%pc13s%storvegc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! totvegc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totvegc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totvegc_13', data=pptr%pc13s%totvegc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! totpftc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totpftc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totpftc_13', data=pptr%pc13s%totpftc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !--------------------------------
    ! pft nitrogen state variables
    !--------------------------------
    
    ! leafn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafn', data=pptr%pns%leafn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafn_storage', data=pptr%pns%leafn_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafn_xfer', data=pptr%pns%leafn_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootn', data=pptr%pns%frootn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootn_storage', data=pptr%pns%frootn_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootn_xfer', data=pptr%pns%frootn_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemn', data=pptr%pns%livestemn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemn_storage', data=pptr%pns%livestemn_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemn_xfer', data=pptr%pns%livestemn_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemn', data=pptr%pns%deadstemn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !deadstemn_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemn_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemn_storage', data=pptr%pns%deadstemn_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !deadstemn_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemn_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemn_xfer', data=pptr%pns%deadstemn_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootn', data=pptr%pns%livecrootn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootn_storage', data=pptr%pns%livecrootn_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !livecrootn_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootn_xfer', data=pptr%pns%livecrootn_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootn', data=pptr%pns%deadcrootn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootn_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootn_storage', data=pptr%pns%deadcrootn_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootn_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootn_xfer', data=pptr%pns%deadcrootn_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !retransn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='retransn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='retransn', data=pptr%pns%retransn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool', data=pptr%pns%npool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! pft_ntrunc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_ntrunc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='pft_ntrunc', data=pptr%pns%pft_ntrunc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! dispvegn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='dispvegn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='dispvegn', data=pptr%pns%dispvegn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! storvegn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='storvegn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='storvegn', data=pptr%pns%storvegn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! totvegn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totvegn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totvegn', data=pptr%pns%totvegn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! totpftn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totpftn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totpftn', data=pptr%pns%totpftn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !--------------------------------
    ! pft carbon flux variables
    !--------------------------------

    ! psnsun
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='psnsun', xtype=nf_double,  &
            dim1name='pft',long_name='sunlit leaf photosynthesis',units='umol CO2/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='psnsun', data=pptr%pcf%psnsun, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! psnsha
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='psnsha', xtype=nf_double,  &
            dim1name='pft',long_name='shaded leaf photosynthesis',units='umol CO2/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='psnsha', data=pptr%pcf%psnsha, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! Gap mortality fluxes
    
    ! m_leafc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_to_litter', data=pptr%pcf%m_leafc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafc_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_storage_to_litter', data=pptr%pcf%m_leafc_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafc_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_xfer_to_litter', data=pptr%pcf%m_leafc_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_to_litter', data=pptr%pcf%m_frootc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootc_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_storage_to_litter', data=pptr%pcf%m_frootc_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootc_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_xfer_to_litter', data=pptr%pcf%m_frootc_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_to_litter', data=pptr%pcf%m_livestemc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemc_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_storage_to_litter', data=pptr%pcf%m_livestemc_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemc_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_xfer_to_litter', data=pptr%pcf%m_livestemc_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_to_litter', data=pptr%pcf%m_deadstemc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_storage_to_litter', data=pptr%pcf%m_deadstemc_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_xfer_to_litter', data=pptr%pcf%m_deadstemc_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_to_litter', data=pptr%pcf%m_livecrootc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootc_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_storage_to_litter', data=pptr%pcf%m_livecrootc_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootc_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_xfer_to_litter', data=pptr%pcf%m_livecrootc_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_to_litter', data=pptr%pcf%m_deadcrootc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_storage_to_litter', data=pptr%pcf%m_deadcrootc_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_xfer_to_litter', data=pptr%pcf%m_deadcrootc_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_gresp_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_gresp_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_gresp_storage_to_litter', data=pptr%pcf%m_gresp_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_gresp_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_gresp_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_gresp_xfer_to_litter', data=pptr%pcf%m_gresp_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! Fire fluxes 
    
    ! m_leafc_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_to_fire', data=pptr%pcf%m_leafc_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafc_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_storage_to_fire', data=pptr%pcf%m_leafc_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafc_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_xfer_to_fire', data=pptr%pcf%m_leafc_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootc_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_to_fire', data=pptr%pcf%m_frootc_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootc_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_storage_to_fire', data=pptr%pcf%m_frootc_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootc_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_xfer_to_fire', data=pptr%pcf%m_frootc_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemc_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_to_fire', data=pptr%pcf%m_livestemc_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemc_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_storage_to_fire', data=pptr%pcf%m_livestemc_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemc_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_xfer_to_fire', data=pptr%pcf%m_livestemc_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_to_fire', data=pptr%pcf%m_deadstemc_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_to_litter_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_to_litter_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_to_litter_fire', data=pptr%pcf%m_deadstemc_to_litter_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_storage_to_fire', data=pptr%pcf%m_deadstemc_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_xfer_to_fire', data=pptr%pcf%m_deadstemc_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootc_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_to_fire', data=pptr%pcf%m_livecrootc_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootc_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_storage_to_fire', data=pptr%pcf%m_livecrootc_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootc_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_xfer_to_fire', data=pptr%pcf%m_livecrootc_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_to_fire', data=pptr%pcf%m_deadcrootc_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_to_litter_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_to_litter_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_to_litter_fire', data=pptr%pcf%m_deadcrootc_to_litter_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_storage_to_fire', data=pptr%pcf%m_deadcrootc_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_xfer_to_fire', data=pptr%pcf%m_deadcrootc_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_gresp_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_gresp_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_gresp_storage_to_fire', data=pptr%pcf%m_gresp_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_gresp_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_gresp_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_gresp_xfer_to_fire', data=pptr%pcf%m_gresp_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! phenology fluxes from transfer pool
    
    ! leafc_xfer_to_leafc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_xfer_to_leafc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_xfer_to_leafc', data=pptr%pcf%leafc_xfer_to_leafc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootc_xfer_to_frootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_xfer_to_frootc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_xfer_to_frootc', data=pptr%pcf%frootc_xfer_to_frootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc_xfer_to_livestemc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_xfer_to_livestemc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemc_xfer_to_livestemc', data=pptr%pcf%livestemc_xfer_to_livestemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc_xfer_to_deadstemc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc_xfer_to_deadstemc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemc_xfer_to_deadstemc', data=pptr%pcf%deadstemc_xfer_to_deadstemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc_xfer_to_livecrootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_xfer_to_livecrootc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootc_xfer_to_livecrootc', data=pptr%pcf%livecrootc_xfer_to_livecrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc_xfer_to_deadcrootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc_xfer_to_deadcrootc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootc_xfer_to_deadcrootc', data=pptr%pcf%deadcrootc_xfer_to_deadcrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! litterfall fluxes 
    
    ! leafc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_to_litter', data=pptr%pcf%leafc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_to_litter', data=pptr%pcf%frootc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! maintenance respiration fluxes 
    
    ! leaf_mr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leaf_mr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leaf_mr', data=pptr%pcf%leaf_mr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! froot_mr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='froot_mr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='froot_mr', data=pptr%pcf%froot_mr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestem_mr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestem_mr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestem_mr', data=pptr%pcf%livestem_mr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecroot_mr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecroot_mr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecroot_mr', data=pptr%pcf%livecroot_mr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leaf_curmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leaf_curmr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leaf_curmr', data=pptr%pcf%leaf_curmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! froot_curmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='froot_curmr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='froot_curmr', data=pptr%pcf%froot_curmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestem_curmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestem_curmr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestem_curmr', data=pptr%pcf%livestem_curmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecroot_curmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecroot_curmr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecroot_curmr', data=pptr%pcf%livecroot_curmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leaf_xsmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leaf_xsmr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leaf_xsmr', data=pptr%pcf%leaf_xsmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! froot_xsmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='froot_xsmr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='froot_xsmr', data=pptr%pcf%froot_xsmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestem_xsmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestem_xsmr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestem_xsmr', data=pptr%pcf%livestem_xsmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecroot_xsmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecroot_xsmr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecroot_xsmr', data=pptr%pcf%livecroot_xsmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! photosynthesis fluxes
    
    ! psnsun_to_cpool 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='psnsun_to_cpool', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='psnsun_to_cpool', data=pptr%pcf%psnsun_to_cpool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! psnshade_to_cpool 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='psnshade_to_cpool', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='psnshade_to_cpool', data=pptr%pcf%psnshade_to_cpool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! allocation fluxes, from current GPP
    
    ! cpool_to_xsmrpool 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_xsmrpool', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_xsmrpool', data=pptr%pcf%cpool_to_xsmrpool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_leafc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_leafc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_leafc', data=pptr%pcf%cpool_to_leafc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_leafc_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_leafc_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_leafc_storage', data=pptr%pcf%cpool_to_leafc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_frootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_frootc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_frootc', data=pptr%pcf%cpool_to_frootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_frootc_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_frootc_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_frootc_storage', data=pptr%pcf%cpool_to_frootc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_livestemc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_livestemc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_livestemc', data=pptr%pcf%cpool_to_livestemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_livestemc_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_livestemc_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_livestemc_storage', data=pptr%pcf%cpool_to_livestemc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_deadstemc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_deadstemc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_deadstemc', data=pptr%pcf%cpool_to_deadstemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_deadstemc_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_deadstemc_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_deadstemc_storage', data=pptr%pcf%cpool_to_deadstemc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_livecrootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_livecrootc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_livecrootc', data=pptr%pcf%cpool_to_livecrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_livecrootc_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_livecrootc_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_livecrootc_storage', data=pptr%pcf%cpool_to_livecrootc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_deadcrootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_deadcrootc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_deadcrootc', data=pptr%pcf%cpool_to_deadcrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_deadcrootc_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_deadcrootc_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_deadcrootc_storage', data=pptr%pcf%cpool_to_deadcrootc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_gresp_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_gresp_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_gresp_storage', data=pptr%pcf%cpool_to_gresp_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! growth respiration fluxes
    
    ! cpool_leaf_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_leaf_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_leaf_gr', data=pptr%pcf%cpool_leaf_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_leaf_storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_leaf_storage_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_leaf_storage_gr', data=pptr%pcf%cpool_leaf_storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_leaf_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_leaf_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_leaf_gr', data=pptr%pcf%transfer_leaf_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_froot_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_froot_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_froot_gr', data=pptr%pcf%cpool_froot_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_froot_storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_froot_storage_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_froot_storage_gr', data=pptr%pcf%cpool_froot_storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_froot_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_froot_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_froot_gr', data=pptr%pcf%transfer_froot_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_livestem_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_livestem_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_livestem_gr', data=pptr%pcf%cpool_livestem_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_livestem_storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_livestem_storage_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_livestem_storage_gr', data=pptr%pcf%cpool_livestem_storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_livestem_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_livestem_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_livestem_gr', data=pptr%pcf%transfer_livestem_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_deadstem_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_deadstem_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_deadstem_gr', data=pptr%pcf%cpool_deadstem_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_deadstem_storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_deadstem_storage_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_deadstem_storage_gr', data=pptr%pcf%cpool_deadstem_storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_deadstem_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_deadstem_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_deadstem_gr', data=pptr%pcf%transfer_deadstem_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_livecroot_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_livecroot_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_livecroot_gr', data=pptr%pcf%cpool_livecroot_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_livecroot_storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_livecroot_storage_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_livecroot_storage_gr', data=pptr%pcf%cpool_livecroot_storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_livecroot_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_livecroot_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_livecroot_gr', data=pptr%pcf%transfer_livecroot_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_deadcroot_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_deadcroot_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_deadcroot_gr', data=pptr%pcf%cpool_deadcroot_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_deadcroot_storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_deadcroot_storage_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_deadcroot_storage_gr', data=pptr%pcf%cpool_deadcroot_storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_deadcroot_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_deadcroot_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_deadcroot_gr', data=pptr%pcf%transfer_deadcroot_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! turnover of storage to transfer pools
    
    ! leafc_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_storage_to_xfer', data=pptr%pcf%leafc_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootc_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_storage_to_xfer', data=pptr%pcf%frootc_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemc_storage_to_xfer', data=pptr%pcf%livestemc_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemc_storage_to_xfer', data=pptr%pcf%deadstemc_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootc_storage_to_xfer', data=pptr%pcf%livecrootc_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootc_storage_to_xfer', data=pptr%pcf%deadcrootc_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! gresp_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gresp_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gresp_storage_to_xfer', data=pptr%pcf%gresp_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! turnover of livewood to deadwood
    
    ! livestemc_to_deadstemc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_to_deadstemc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemc_to_deadstemc', data=pptr%pcf%livestemc_to_deadstemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc_to_deadcrootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_to_deadcrootc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootc_to_deadcrootc', data=pptr%pcf%livecrootc_to_deadcrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! summary (diagnostic) flux variables
    
    ! gpp 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gpp_pcf', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gpp_pcf', data=pptr%pcf%gpp, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! mr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='mr', data=pptr%pcf%mr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! current_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='current_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='current_gr', data=pptr%pcf%current_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_gr', data=pptr%pcf%transfer_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='storage_gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='storage_gr', data=pptr%pcf%storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gr', data=pptr%pcf%gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! ar 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ar', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ar', data=pptr%pcf%ar, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! rr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='rr', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='rr', data=pptr%pcf%rr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npp 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npp', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npp', data=pptr%pcf%npp, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! agnpp 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='agnpp', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='agnpp', data=pptr%pcf%agnpp, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! bgnpp 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='bgnpp', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='bgnpp', data=pptr%pcf%bgnpp, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! litfall 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litfall', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litfall', data=pptr%pcf%litfall, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! vegfire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vegfire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='vegfire', data=pptr%pcf%vegfire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! pft_cinputs 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_cinputs', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='pft_cinputs', data=pptr%pcf%pft_cinputs, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! pft_coutputs 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_coutputs', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='pft_coutputs', data=pptr%pcf%pft_coutputs, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! pft_fire_closs 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_fire_closs', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='pft_fire_closs', data=pptr%pcf%pft_fire_closs, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !--------------------------------
    ! C13 pft carbon flux variables
    !--------------------------------
    ! Gap mortality fluxes
    
    ! m_leafc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_to_litter_13', data=pptr%pc13f%m_leafc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafc_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_storage_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_storage_to_litter_13', data=pptr%pc13f%m_leafc_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafc_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_xfer_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_xfer_to_litter_13', data=pptr%pc13f%m_leafc_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_to_litter_13', data=pptr%pc13f%m_frootc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootc_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_storage_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_storage_to_litter_13', data=pptr%pc13f%m_frootc_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootc_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_xfer_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_xfer_to_litter_13', data=pptr%pc13f%m_frootc_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_to_litter_13', data=pptr%pc13f%m_livestemc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemc_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_storage_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_storage_to_litter_13', data=pptr%pc13f%m_livestemc_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemc_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_xfer_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_xfer_to_litter_13', data=pptr%pc13f%m_livestemc_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_to_litter_13', data=pptr%pc13f%m_deadstemc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_storage_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_storage_to_litter_13', data=pptr%pc13f%m_deadstemc_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_xfer_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_xfer_to_litter_13', data=pptr%pc13f%m_deadstemc_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_to_litter_13', data=pptr%pc13f%m_livecrootc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootc_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_storage_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_storage_to_litter_13', data=pptr%pc13f%m_livecrootc_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootc_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_xfer_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_xfer_to_litter_13', data=pptr%pc13f%m_livecrootc_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_to_litter_13', data=pptr%pc13f%m_deadcrootc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_storage_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_storage_to_litter_13', data=pptr%pc13f%m_deadcrootc_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_xfer_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_xfer_to_litter_13', data=pptr%pc13f%m_deadcrootc_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_gresp_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_gresp_storage_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_gresp_storage_to_litter_13', data=pptr%pc13f%m_gresp_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_gresp_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_gresp_xfer_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_gresp_xfer_to_litter_13', data=pptr%pc13f%m_gresp_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! Fire fluxes 
    
    ! m_leafc_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_to_fire_13', data=pptr%pc13f%m_leafc_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafc_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_storage_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_storage_to_fire_13', data=pptr%pc13f%m_leafc_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafc_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_xfer_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_xfer_to_fire_13', data=pptr%pc13f%m_leafc_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootc_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_to_fire_13', data=pptr%pc13f%m_frootc_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootc_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_storage_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_storage_to_fire_13', data=pptr%pc13f%m_frootc_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootc_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_xfer_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_xfer_to_fire_13', data=pptr%pc13f%m_frootc_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemc_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_to_fire_13', data=pptr%pc13f%m_livestemc_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemc_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_storage_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_storage_to_fire_13', data=pptr%pc13f%m_livestemc_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemc_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_xfer_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_xfer_to_fire_13', data=pptr%pc13f%m_livestemc_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_to_fire_13', data=pptr%pc13f%m_deadstemc_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_to_litter_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_to_litter_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_to_litter_fire_13', data=pptr%pc13f%m_deadstemc_to_litter_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_storage_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_storage_to_fire_13', data=pptr%pc13f%m_deadstemc_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemc_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_xfer_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_xfer_to_fire_13', data=pptr%pc13f%m_deadstemc_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootc_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_to_fire_13', data=pptr%pc13f%m_livecrootc_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootc_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_storage_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_storage_to_fire_13', data=pptr%pc13f%m_livecrootc_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootc_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_xfer_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_xfer_to_fire_13', data=pptr%pc13f%m_livecrootc_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_to_fire_13', data=pptr%pc13f%m_deadcrootc_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_to_litter_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_to_litter_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_to_litter_fire_13', data=pptr%pc13f%m_deadcrootc_to_litter_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_storage_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_storage_to_fire_13', data=pptr%pc13f%m_deadcrootc_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootc_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_xfer_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_xfer_to_fire_13', data=pptr%pc13f%m_deadcrootc_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_gresp_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_gresp_storage_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_gresp_storage_to_fire_13', data=pptr%pc13f%m_gresp_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_gresp_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_gresp_xfer_to_fire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_gresp_xfer_to_fire_13', data=pptr%pc13f%m_gresp_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! phenology fluxes from transfer pool
    
    ! leafc_xfer_to_leafc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_xfer_to_leafc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_xfer_to_leafc_13', data=pptr%pc13f%leafc_xfer_to_leafc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootc_xfer_to_frootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_xfer_to_frootc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_xfer_to_frootc_13', data=pptr%pc13f%frootc_xfer_to_frootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc_xfer_to_livestemc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_xfer_to_livestemc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemc_xfer_to_livestemc_13', data=pptr%pc13f%livestemc_xfer_to_livestemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc_xfer_to_deadstemc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc_xfer_to_deadstemc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemc_xfer_to_deadstemc_13', data=pptr%pc13f%deadstemc_xfer_to_deadstemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc_xfer_to_livecrootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_xfer_to_livecrootc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootc_xfer_to_livecrootc_13', data=pptr%pc13f%livecrootc_xfer_to_livecrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc_xfer_to_deadcrootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc_xfer_to_deadcrootc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootc_xfer_to_deadcrootc_13', data=pptr%pc13f%deadcrootc_xfer_to_deadcrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! litterfall fluxes 
    
    ! leafc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_to_litter_13', data=pptr%pc13f%leafc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootc_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_to_litter_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_to_litter_13', data=pptr%pc13f%frootc_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! maintenance respiration fluxes 
    
    ! leaf_mr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leaf_mr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leaf_mr_13', data=pptr%pc13f%leaf_mr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! froot_mr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='froot_mr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='froot_mr_13', data=pptr%pc13f%froot_mr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestem_mr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestem_mr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestem_mr_13', data=pptr%pc13f%livestem_mr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecroot_mr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecroot_mr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecroot_mr_13', data=pptr%pc13f%livecroot_mr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leaf_curmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leaf_curmr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leaf_curmr_13', data=pptr%pc13f%leaf_curmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! froot_curmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='froot_curmr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='froot_curmr_13', data=pptr%pc13f%froot_curmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestem_curmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestem_curmr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestem_curmr_13', data=pptr%pc13f%livestem_curmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecroot_curmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecroot_curmr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecroot_curmr_13', data=pptr%pc13f%livecroot_curmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leaf_xsmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leaf_xsmr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leaf_xsmr_13', data=pptr%pc13f%leaf_xsmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! froot_xsmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='froot_xsmr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='froot_xsmr_13', data=pptr%pc13f%froot_xsmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestem_xsmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestem_xsmr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestem_xsmr_13', data=pptr%pc13f%livestem_xsmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecroot_xsmr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecroot_xsmr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecroot_xsmr_13', data=pptr%pc13f%livecroot_xsmr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! photosynthesis fluxes
    
    ! psnsun_to_cpool 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='psnsun_to_cpool_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='psnsun_to_cpool_13', data=pptr%pc13f%psnsun_to_cpool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! psnshade_to_cpool 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='psnshade_to_cpool_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='psnshade_to_cpool_13', data=pptr%pc13f%psnshade_to_cpool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! allocation fluxes, from current GPP
    
    ! cpool_to_xsmrpool 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_xsmrpool_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_xsmrpool_13', data=pptr%pc13f%cpool_to_xsmrpool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_leafc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_leafc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_leafc_13', data=pptr%pc13f%cpool_to_leafc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_leafc_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_leafc_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_leafc_storage_13', data=pptr%pc13f%cpool_to_leafc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_frootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_frootc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_frootc_13', data=pptr%pc13f%cpool_to_frootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_frootc_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_frootc_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_frootc_storage_13', data=pptr%pc13f%cpool_to_frootc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_livestemc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_livestemc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_livestemc_13', data=pptr%pc13f%cpool_to_livestemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_livestemc_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_livestemc_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_livestemc_storage_13', data=pptr%pc13f%cpool_to_livestemc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_deadstemc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_deadstemc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_deadstemc_13', data=pptr%pc13f%cpool_to_deadstemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_deadstemc_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_deadstemc_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_deadstemc_storage_13', data=pptr%pc13f%cpool_to_deadstemc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_livecrootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_livecrootc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_livecrootc_13', data=pptr%pc13f%cpool_to_livecrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_livecrootc_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_livecrootc_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_livecrootc_storage_13', data=pptr%pc13f%cpool_to_livecrootc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_deadcrootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_deadcrootc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_deadcrootc_13', data=pptr%pc13f%cpool_to_deadcrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_deadcrootc_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_deadcrootc_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_deadcrootc_storage_13', data=pptr%pc13f%cpool_to_deadcrootc_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_to_gresp_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_gresp_storage_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_to_gresp_storage_13', data=pptr%pc13f%cpool_to_gresp_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! growth respiration fluxes
    
    ! cpool_leaf_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_leaf_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_leaf_gr_13', data=pptr%pc13f%cpool_leaf_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_leaf_storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_leaf_storage_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_leaf_storage_gr_13', data=pptr%pc13f%cpool_leaf_storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_leaf_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_leaf_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_leaf_gr_13', data=pptr%pc13f%transfer_leaf_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_froot_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_froot_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_froot_gr_13', data=pptr%pc13f%cpool_froot_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_froot_storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_froot_storage_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_froot_storage_gr_13', data=pptr%pc13f%cpool_froot_storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_froot_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_froot_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_froot_gr_13', data=pptr%pc13f%transfer_froot_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_livestem_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_livestem_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_livestem_gr_13', data=pptr%pc13f%cpool_livestem_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_livestem_storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_livestem_storage_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_livestem_storage_gr_13', data=pptr%pc13f%cpool_livestem_storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_livestem_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_livestem_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_livestem_gr_13', data=pptr%pc13f%transfer_livestem_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_deadstem_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_deadstem_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_deadstem_gr_13', data=pptr%pc13f%cpool_deadstem_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_deadstem_storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_deadstem_storage_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_deadstem_storage_gr_13', data=pptr%pc13f%cpool_deadstem_storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_deadstem_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_deadstem_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_deadstem_gr_13', data=pptr%pc13f%transfer_deadstem_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_livecroot_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_livecroot_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_livecroot_gr_13', data=pptr%pc13f%cpool_livecroot_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_livecroot_storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_livecroot_storage_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_livecroot_storage_gr_13', data=pptr%pc13f%cpool_livecroot_storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_livecroot_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_livecroot_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_livecroot_gr_13', data=pptr%pc13f%transfer_livecroot_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_deadcroot_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_deadcroot_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_deadcroot_gr_13', data=pptr%pc13f%cpool_deadcroot_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! cpool_deadcroot_storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_deadcroot_storage_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cpool_deadcroot_storage_gr_13', data=pptr%pc13f%cpool_deadcroot_storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_deadcroot_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_deadcroot_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_deadcroot_gr_13', data=pptr%pc13f%transfer_deadcroot_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! turnover of storage to transfer pools
    
    ! leafc_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_storage_to_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_storage_to_xfer_13', data=pptr%pc13f%leafc_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootc_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_storage_to_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_storage_to_xfer_13', data=pptr%pc13f%frootc_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemc_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_storage_to_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemc_storage_to_xfer_13', data=pptr%pc13f%livestemc_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemc_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemc_storage_to_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemc_storage_to_xfer_13', data=pptr%pc13f%deadstemc_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_storage_to_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootc_storage_to_xfer_13', data=pptr%pc13f%livecrootc_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootc_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootc_storage_to_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootc_storage_to_xfer_13', data=pptr%pc13f%deadcrootc_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! gresp_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gresp_storage_to_xfer_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gresp_storage_to_xfer_13', data=pptr%pc13f%gresp_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! turnover of livewood to deadwood
    
    ! livestemc_to_deadstemc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_to_deadstemc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemc_to_deadstemc_13', data=pptr%pc13f%livestemc_to_deadstemc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootc_to_deadcrootc 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootc_to_deadcrootc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootc_to_deadcrootc_13', data=pptr%pc13f%livecrootc_to_deadcrootc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! summary (diagnostic) flux variables
    
    ! gpp 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gpp_pc13f_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gpp_pc13f_13', data=pptr%pc13f%gpp, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! mr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='mr_13', data=pptr%pc13f%mr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! current_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='current_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='current_gr_13', data=pptr%pc13f%current_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! transfer_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='transfer_gr_13', data=pptr%pc13f%transfer_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! storage_gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='storage_gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='storage_gr_13', data=pptr%pc13f%storage_gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! gr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gr_13', data=pptr%pc13f%gr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! ar 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ar_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ar_13', data=pptr%pc13f%ar, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! rr 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='rr_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='rr_13', data=pptr%pc13f%rr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npp 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npp_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npp_13', data=pptr%pc13f%npp, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! agnpp 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='agnpp_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='agnpp_13', data=pptr%pc13f%agnpp, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! bgnpp 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='bgnpp_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='bgnpp_13', data=pptr%pc13f%bgnpp, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! litfall 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litfall_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litfall_13', data=pptr%pc13f%litfall, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! vegfire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vegfire_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='vegfire_13', data=pptr%pc13f%vegfire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! pft_cinputs 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_cinputs_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='pft_cinputs_13', data=pptr%pc13f%pft_cinputs, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! pft_coutputs 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_coutputs_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='pft_coutputs_13', data=pptr%pc13f%pft_coutputs, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! pft_fire_closs 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_fire_closs_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='pft_fire_closs_13', data=pptr%pc13f%pft_fire_closs, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !--------------------------------
    ! pft nitrogen flux variables
    !--------------------------------
    
    ! m_leafn_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafn_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafn_to_litter', data=pptr%pnf%m_leafn_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootn_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootn_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootn_to_litter', data=pptr%pnf%m_frootn_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafn_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafn_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafn_storage_to_litter', data=pptr%pnf%m_leafn_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootn_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootn_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootn_storage_to_litter', data=pptr%pnf%m_frootn_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemn_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemn_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemn_storage_to_litter', data=pptr%pnf%m_livestemn_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemn_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemn_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemn_storage_to_litter', data=pptr%pnf%m_deadstemn_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootn_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootn_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootn_storage_to_litter', data=pptr%pnf%m_livecrootn_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootn_storage_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootn_storage_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootn_storage_to_litter', data=pptr%pnf%m_deadcrootn_storage_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafn_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafn_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafn_xfer_to_litter', data=pptr%pnf%m_leafn_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootn_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootn_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootn_xfer_to_litter', data=pptr%pnf%m_frootn_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemn_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemn_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemn_xfer_to_litter', data=pptr%pnf%m_livestemn_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemn_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemn_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemn_xfer_to_litter', data=pptr%pnf%m_deadstemn_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootn_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootn_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootn_xfer_to_litter', data=pptr%pnf%m_livecrootn_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootn_xfer_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootn_xfer_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootn_xfer_to_litter', data=pptr%pnf%m_deadcrootn_xfer_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemn_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemn_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemn_to_litter', data=pptr%pnf%m_livestemn_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemn_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemn_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemn_to_litter', data=pptr%pnf%m_deadstemn_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootn_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootn_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootn_to_litter', data=pptr%pnf%m_livecrootn_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootn_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootn_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootn_to_litter', data=pptr%pnf%m_deadcrootn_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_retransn_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_retransn_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_retransn_to_litter', data=pptr%pnf%m_retransn_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafn_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafn_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafn_to_fire', data=pptr%pnf%m_leafn_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootn_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootn_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootn_to_fire', data=pptr%pnf%m_frootn_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafn_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafn_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafn_storage_to_fire', data=pptr%pnf%m_leafn_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootn_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootn_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootn_storage_to_fire', data=pptr%pnf%m_frootn_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemn_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemn_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemn_storage_to_fire', data=pptr%pnf%m_livestemn_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemn_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemn_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemn_storage_to_fire', data=pptr%pnf%m_deadstemn_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootn_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootn_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootn_storage_to_fire', data=pptr%pnf%m_livecrootn_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootn_storage_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootn_storage_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootn_storage_to_fire', data=pptr%pnf%m_deadcrootn_storage_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_leafn_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafn_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafn_xfer_to_fire', data=pptr%pnf%m_leafn_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_frootn_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootn_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootn_xfer_to_fire', data=pptr%pnf%m_frootn_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemn_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemn_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemn_xfer_to_fire', data=pptr%pnf%m_livestemn_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemn_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemn_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemn_xfer_to_fire', data=pptr%pnf%m_deadstemn_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootn_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootn_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootn_xfer_to_fire', data=pptr%pnf%m_livecrootn_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootn_xfer_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootn_xfer_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootn_xfer_to_fire', data=pptr%pnf%m_deadcrootn_xfer_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livestemn_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemn_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemn_to_fire', data=pptr%pnf%m_livestemn_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemn_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemn_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemn_to_fire', data=pptr%pnf%m_deadstemn_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadstemn_to_litter_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemn_to_litter_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemn_to_litter_fire', data=pptr%pnf%m_deadstemn_to_litter_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_livecrootn_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootn_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootn_to_fire', data=pptr%pnf%m_livecrootn_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootn_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootn_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootn_to_fire', data=pptr%pnf%m_deadcrootn_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_deadcrootn_to_litter_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootn_to_litter_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootn_to_litter_fire', data=pptr%pnf%m_deadcrootn_to_litter_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! m_retransn_to_fire 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_retransn_to_fire', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_retransn_to_fire', data=pptr%pnf%m_retransn_to_fire, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafn_xfer_to_leafn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_xfer_to_leafn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafn_xfer_to_leafn', data=pptr%pnf%leafn_xfer_to_leafn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootn_xfer_to_frootn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn_xfer_to_frootn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootn_xfer_to_frootn', data=pptr%pnf%frootn_xfer_to_frootn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemn_xfer_to_livestemn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn_xfer_to_livestemn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemn_xfer_to_livestemn', data=pptr%pnf%livestemn_xfer_to_livestemn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemn_xfer_to_deadstemn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemn_xfer_to_deadstemn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemn_xfer_to_deadstemn', data=pptr%pnf%deadstemn_xfer_to_deadstemn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootn_xfer_to_livecrootn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn_xfer_to_livecrootn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootn_xfer_to_livecrootn', data=pptr%pnf%livecrootn_xfer_to_livecrootn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootn_xfer_to_deadcrootn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootn_xfer_to_deadcrootn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootn_xfer_to_deadcrootn', data=pptr%pnf%deadcrootn_xfer_to_deadcrootn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafn_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafn_to_litter', data=pptr%pnf%leafn_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootn_to_litter 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn_to_litter', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootn_to_litter', data=pptr%pnf%frootn_to_litter, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafn_to_retransn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_to_retransn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafn_to_retransn', data=pptr%pnf%leafn_to_retransn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! retransn_to_npool 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='retransn_to_npool', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='retransn_to_npool', data=pptr%pnf%retransn_to_npool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! sminn_to_npool 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_npool', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_npool', data=pptr%pnf%sminn_to_npool, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool_to_leafn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_leafn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool_to_leafn', data=pptr%pnf%npool_to_leafn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool_to_leafn_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_leafn_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool_to_leafn_storage', data=pptr%pnf%npool_to_leafn_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool_to_frootn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_frootn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool_to_frootn', data=pptr%pnf%npool_to_frootn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool_to_frootn_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_frootn_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool_to_frootn_storage', data=pptr%pnf%npool_to_frootn_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool_to_livestemn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_livestemn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool_to_livestemn', data=pptr%pnf%npool_to_livestemn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool_to_livestemn_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_livestemn_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool_to_livestemn_storage', data=pptr%pnf%npool_to_livestemn_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool_to_deadstemn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_deadstemn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool_to_deadstemn', data=pptr%pnf%npool_to_deadstemn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool_to_deadstemn_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_deadstemn_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool_to_deadstemn_storage', data=pptr%pnf%npool_to_deadstemn_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool_to_livecrootn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_livecrootn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool_to_livecrootn', data=pptr%pnf%npool_to_livecrootn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool_to_livecrootn_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_livecrootn_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool_to_livecrootn_storage', data=pptr%pnf%npool_to_livecrootn_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool_to_deadcrootn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_deadcrootn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool_to_deadcrootn', data=pptr%pnf%npool_to_deadcrootn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! npool_to_deadcrootn_storage 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_deadcrootn_storage', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='npool_to_deadcrootn_storage', data=pptr%pnf%npool_to_deadcrootn_storage, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! leafn_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafn_storage_to_xfer', data=pptr%pnf%leafn_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! frootn_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootn_storage_to_xfer', data=pptr%pnf%frootn_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livestemn_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemn_storage_to_xfer', data=pptr%pnf%livestemn_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadstemn_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemn_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadstemn_storage_to_xfer', data=pptr%pnf%deadstemn_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! livecrootn_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootn_storage_to_xfer', data=pptr%pnf%livecrootn_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! deadcrootn_storage_to_xfer 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootn_storage_to_xfer', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='deadcrootn_storage_to_xfer', data=pptr%pnf%deadcrootn_storage_to_xfer, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! livestemn_to_deadstemn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn_to_deadstemn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemn_to_deadstemn', data=pptr%pnf%livestemn_to_deadstemn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! livestemn_to_retransn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn_to_retransn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livestemn_to_retransn', data=pptr%pnf%livestemn_to_retransn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! livecrootn_to_deadcrootn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn_to_deadcrootn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootn_to_deadcrootn', data=pptr%pnf%livecrootn_to_deadcrootn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! livecrootn_to_retransn 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn_to_retransn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='livecrootn_to_retransn', data=pptr%pnf%livecrootn_to_retransn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! ndeploy 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ndeploy', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ndeploy', data=pptr%pnf%ndeploy, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft_ninputs 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_ninputs', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='pft_ninputs', data=pptr%pnf%pft_ninputs, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft_noutputs 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_noutputs', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='pft_noutputs', data=pptr%pnf%pft_noutputs, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft_fire_nloss 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_fire_nloss', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='pft_fire_nloss', data=pptr%pnf%pft_fire_nloss, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !--------------------------------
    ! column physical state variables
    !--------------------------------
    
    ! bsw2 (2d-numrad)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='bsw2', xtype=nf_double,   &
            dim1name='column', dim2name='levsoi', &
            long_name='Clapp and Hornberger "b" for CN code', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='bsw2', data=cptr%cps%bsw2, &
            dim1name='column', dim2name='levsoi', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! psisat (2d-numrad)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='psisat', xtype=nf_double,   &
            dim1name='column', dim2name='levsoi', &
            long_name='soil water potential at saturation for CN code', units='MPa')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='psisat', data=cptr%cps%psisat, &
            dim1name='column', dim2name='levsoi', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! vwcsat (2d-numrad)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vwcsat', xtype=nf_double,   &
            dim1name='column', dim2name='levsoi', &
            long_name='volumetric water content at saturation for CN code ', units='m3/m3')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='vwcsat', data=cptr%cps%vwcsat, &
            dim1name='column', dim2name='levsoi', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! soilpsi (2d-numrad)
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soilpsi', xtype=nf_double,  &
            dim1name='column', dim2name='levsoi', &
            long_name='soil water potential',units='MPa')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soilpsi', data=cptr%cps%soilpsi, &
            dim1name='column', dim2name='levsoi', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! decl
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='decl', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='decl', data=cptr%cps%decl, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! fpi
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fpi', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='fpi', data=cptr%cps%fpi, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! fpg
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fpg', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='fpg', data=cptr%cps%fpg, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! annsum_counter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_counter', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='annsum_counter', data=cptr%cps%annsum_counter, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! cannsum_npp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cannsum_npp', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cannsum_npp', data=cptr%cps%cannsum_npp, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! wf
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='wf', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='wf', data=cptr%cps%wf, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! me
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='me', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='me', data=cptr%cps%me, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! fire_prob
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fire_prob', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='fire_prob', data=cptr%cps%fire_prob, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! mean_fire_prob
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mean_fire_prob', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='mean_fire_prob', data=cptr%cps%mean_fire_prob, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! fireseasonl
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fireseasonl', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='fireseasonl', data=cptr%cps%fireseasonl, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! farea_burned
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='farea_burned', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='farea_burned', data=cptr%cps%farea_burned, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! ann_farea_burned
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ann_farea_burned', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ann_farea_burned', data=cptr%cps%ann_farea_burned, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !--------------------------------
    ! column carbon state variables
    !--------------------------------

    ! cwdc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdc', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cwdc', data=cptr%ccs%cwdc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr1c', data=cptr%ccs%litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !litr2c 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr2c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr2c', data=cptr%ccs%litr2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr3c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr3c', data=cptr%ccs%litr3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !soil1c 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil1c', data=cptr%ccs%soil1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil2c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil2c', data=cptr%ccs%soil2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil3c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil3c', data=cptr%ccs%soil3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil4c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil4c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil4c', data=cptr%ccs%soil4c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! col_ctrunc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_ctrunc', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='col_ctrunc', data=cptr%ccs%col_ctrunc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! totlitc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totlitc', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totlitc', data=cptr%ccs%totlitc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! totsomc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totsomc', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totsomc', data=cptr%ccs%totsomc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! totecosysc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totecosysc', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totecosysc', data=cptr%ccs%totecosysc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! totcolc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totcolc', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totcolc', data=cptr%ccs%totcolc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !--------------------------------
    ! C13 column carbon state variables
    !--------------------------------

    ! cwdc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdc_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cwdc_13', data=cptr%cc13s%cwdc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr1c_13', data=cptr%cc13s%litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !litr2c 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr2c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr2c_13', data=cptr%cc13s%litr2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr3c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr3c_13', data=cptr%cc13s%litr3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !soil1c 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil1c_13', data=cptr%cc13s%soil1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil2c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil2c_13', data=cptr%cc13s%soil2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil3c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil3c_13', data=cptr%cc13s%soil3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! soil4c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil4c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil4c_13', data=cptr%cc13s%soil4c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! col_ctrunc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_ctrunc_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='col_ctrunc_13', data=cptr%cc13s%col_ctrunc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! totlitc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totlitc_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totlitc_13', data=cptr%cc13s%totlitc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! totsomc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totsomc_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totsomc_13', data=cptr%cc13s%totsomc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! totecosysc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totecosysc_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totecosysc_13', data=cptr%cc13s%totecosysc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! totcolc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totcolc_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totcolc_13', data=cptr%cc13s%totcolc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !--------------------------------
    ! column nitrogen state variables
    !--------------------------------
    
    ! cwdn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cwdn', data=cptr%cns%cwdn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !litr1n 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr1n', data=cptr%cns%litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! litr2n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr2n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr2n', data=cptr%cns%litr2n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! litr3n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr3n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr3n', data=cptr%cns%litr3n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! soil1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil1n', data=cptr%cns%soil1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! soil2n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil2n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil2n', data=cptr%cns%soil2n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! soil3n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil3n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil3n', data=cptr%cns%soil3n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! soil4n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil4n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil4n', data=cptr%cns%soil4n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! sminn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn', data=cptr%cns%sminn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! col_ntrunc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_ntrunc', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='col_ntrunc', data=cptr%cns%col_ntrunc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! totlitn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totlitn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totlitn', data=cptr%cns%totlitn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! totsomn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totsomn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totsomn', data=cptr%cns%totsomn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! totecosysn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totecosysn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totecosysn', data=cptr%cns%totecosysn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! totcoln
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totcoln', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totcoln', data=cptr%cns%totcoln, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !--------------------------------
    ! column carbon flux variables
    !--------------------------------

    ! m_leafc_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_to_litr1c', data=cptr%ccf%m_leafc_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafc_to_litr2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_to_litr2c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_to_litr2c', data=cptr%ccf%m_leafc_to_litr2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafc_to_litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_to_litr3c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_to_litr3c', data=cptr%ccf%m_leafc_to_litr3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootc_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_to_litr1c', data=cptr%ccf%m_frootc_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootc_to_litr2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_to_litr2c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_to_litr2c', data=cptr%ccf%m_frootc_to_litr2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootc_to_litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_to_litr3c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_to_litr3c', data=cptr%ccf%m_frootc_to_litr3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafc_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_storage_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_storage_to_litr1c', data=cptr%ccf%m_leafc_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootc_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_storage_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_storage_to_litr1c', data=cptr%ccf%m_frootc_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livestemc_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_storage_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_storage_to_litr1c', data=cptr%ccf%m_livestemc_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadstemc_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_storage_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_storage_to_litr1c', data=cptr%ccf%m_deadstemc_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livecrootc_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_storage_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_storage_to_litr1c', data=cptr%ccf%m_livecrootc_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadcrootc_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_storage_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_storage_to_litr1c', data=cptr%ccf%m_deadcrootc_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafc_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_xfer_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_xfer_to_litr1c', data=cptr%ccf%m_leafc_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootc_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_xfer_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_xfer_to_litr1c', data=cptr%ccf%m_frootc_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livestemc_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_xfer_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_xfer_to_litr1c', data=cptr%ccf%m_livestemc_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadstemc_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_xfer_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_xfer_to_litr1c', data=cptr%ccf%m_deadstemc_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livecrootc_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_xfer_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_xfer_to_litr1c', data=cptr%ccf%m_livecrootc_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadcrootc_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_xfer_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_xfer_to_litr1c', data=cptr%ccf%m_deadcrootc_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livestemc_to_cwdc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_to_cwdc', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_to_cwdc', data=cptr%ccf%m_livestemc_to_cwdc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadstemc_to_cwdc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_to_cwdc', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_to_cwdc', data=cptr%ccf%m_deadstemc_to_cwdc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livecrootc_to_cwdc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_to_cwdc', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_to_cwdc', data=cptr%ccf%m_livecrootc_to_cwdc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadcrootc_to_cwdc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_to_cwdc', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_to_cwdc', data=cptr%ccf%m_deadcrootc_to_cwdc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_gresp_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_gresp_storage_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_gresp_storage_to_litr1c', data=cptr%ccf%m_gresp_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_gresp_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_gresp_xfer_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_gresp_xfer_to_litr1c', data=cptr%ccf%m_gresp_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadstemc_to_cwdc_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_to_cwdc_fire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_to_cwdc_fire', data=cptr%ccf%m_deadstemc_to_cwdc_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadcrootc_to_cwdc_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_to_cwdc_fire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_to_cwdc_fire', data=cptr%ccf%m_deadcrootc_to_cwdc_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_litr1c_to_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_litr1c_to_fire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_litr1c_to_fire', data=cptr%ccf%m_litr1c_to_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_litr2c_to_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_litr2c_to_fire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_litr2c_to_fire', data=cptr%ccf%m_litr2c_to_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_litr3c_to_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_litr3c_to_fire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_litr3c_to_fire', data=cptr%ccf%m_litr3c_to_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_cwdc_to_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_cwdc_to_fire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_cwdc_to_fire', data=cptr%ccf%m_cwdc_to_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! leafc_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_to_litr1c', data=cptr%ccf%leafc_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! leafc_to_litr2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_to_litr2c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_to_litr2c', data=cptr%ccf%leafc_to_litr2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! leafc_to_litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_to_litr3c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_to_litr3c', data=cptr%ccf%leafc_to_litr3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! frootc_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_to_litr1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_to_litr1c', data=cptr%ccf%frootc_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! frootc_to_litr2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_to_litr2c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_to_litr2c', data=cptr%ccf%frootc_to_litr2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! frootc_to_litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_to_litr3c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_to_litr3c', data=cptr%ccf%frootc_to_litr3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! cwdc_to_litr2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdc_to_litr2c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cwdc_to_litr2c', data=cptr%ccf%cwdc_to_litr2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! cwdc_to_litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdc_to_litr3c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cwdc_to_litr3c', data=cptr%ccf%cwdc_to_litr3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr1_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr1_hr', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr1_hr', data=cptr%ccf%litr1_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr1c_to_soil1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr1c_to_soil1c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr1c_to_soil1c', data=cptr%ccf%litr1c_to_soil1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr2_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr2_hr', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr2_hr', data=cptr%ccf%litr2_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr2c_to_soil2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr2c_to_soil2c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr2c_to_soil2c', data=cptr%ccf%litr2c_to_soil2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr3_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr3_hr', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr3_hr', data=cptr%ccf%litr3_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr3c_to_soil3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr3c_to_soil3c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr3c_to_soil3c', data=cptr%ccf%litr3c_to_soil3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil1_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil1_hr', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil1_hr', data=cptr%ccf%soil1_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil1c_to_soil2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil1c_to_soil2c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil1c_to_soil2c', data=cptr%ccf%soil1c_to_soil2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil2_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil2_hr', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil2_hr', data=cptr%ccf%soil2_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil2c_to_soil3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil2c_to_soil3c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil2c_to_soil3c', data=cptr%ccf%soil2c_to_soil3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil3_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil3_hr', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil3_hr', data=cptr%ccf%soil3_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil3c_to_soil4c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil3c_to_soil4c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil3c_to_soil4c', data=cptr%ccf%soil3c_to_soil4c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil4_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil4_hr', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil4_hr', data=cptr%ccf%soil4_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! lithr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='lithr', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='lithr', data=cptr%ccf%lithr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! somhr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='somhr', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='somhr', data=cptr%ccf%somhr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='hr', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='hr', data=cptr%ccf%hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sr', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sr', data=cptr%ccf%sr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! er
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='er', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='er', data=cptr%ccf%er, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litfire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litfire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litfire', data=cptr%ccf%litfire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! somfire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='somfire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='somfire', data=cptr%ccf%somfire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! totfire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totfire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totfire', data=cptr%ccf%totfire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! nep
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='nep', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='nep', data=cptr%ccf%nep, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! nee
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='nee', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='nee', data=cptr%ccf%nee, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! col_cinputs
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_cinputs', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='col_cinputs', data=cptr%ccf%col_cinputs, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! col_coutputs
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_coutputs', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='col_coutputs', data=cptr%ccf%col_coutputs, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! col_fire_closs
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_fire_closs', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='col_fire_closs', data=cptr%ccf%col_fire_closs, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !--------------------------------
    ! C13 column carbon flux variables
    !--------------------------------

    ! m_leafc_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_to_litr1c_13', data=cptr%cc13f%m_leafc_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafc_to_litr2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_to_litr2c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_to_litr2c_13', data=cptr%cc13f%m_leafc_to_litr2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafc_to_litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_to_litr3c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_to_litr3c_13', data=cptr%cc13f%m_leafc_to_litr3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootc_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_to_litr1c_13', data=cptr%cc13f%m_frootc_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootc_to_litr2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_to_litr2c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_to_litr2c_13', data=cptr%cc13f%m_frootc_to_litr2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootc_to_litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_to_litr3c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_to_litr3c_13', data=cptr%cc13f%m_frootc_to_litr3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafc_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_storage_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_storage_to_litr1c_13', data=cptr%cc13f%m_leafc_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootc_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_storage_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_storage_to_litr1c_13', data=cptr%cc13f%m_frootc_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livestemc_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_storage_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_storage_to_litr1c_13', data=cptr%cc13f%m_livestemc_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadstemc_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_storage_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_storage_to_litr1c_13', data=cptr%cc13f%m_deadstemc_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livecrootc_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_storage_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_storage_to_litr1c_13', data=cptr%cc13f%m_livecrootc_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadcrootc_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_storage_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_storage_to_litr1c_13', data=cptr%cc13f%m_deadcrootc_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafc_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafc_xfer_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafc_xfer_to_litr1c_13', data=cptr%cc13f%m_leafc_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootc_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootc_xfer_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootc_xfer_to_litr1c_13', data=cptr%cc13f%m_frootc_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livestemc_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_xfer_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_xfer_to_litr1c_13', data=cptr%cc13f%m_livestemc_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadstemc_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_xfer_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_xfer_to_litr1c_13', data=cptr%cc13f%m_deadstemc_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livecrootc_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_xfer_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_xfer_to_litr1c_13', data=cptr%cc13f%m_livecrootc_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadcrootc_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_xfer_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_xfer_to_litr1c_13', data=cptr%cc13f%m_deadcrootc_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livestemc_to_cwdc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemc_to_cwdc_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemc_to_cwdc_13', data=cptr%cc13f%m_livestemc_to_cwdc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadstemc_to_cwdc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_to_cwdc_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_to_cwdc_13', data=cptr%cc13f%m_deadstemc_to_cwdc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livecrootc_to_cwdc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootc_to_cwdc_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootc_to_cwdc_13', data=cptr%cc13f%m_livecrootc_to_cwdc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadcrootc_to_cwdc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_to_cwdc_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_to_cwdc_13', data=cptr%cc13f%m_deadcrootc_to_cwdc, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_gresp_storage_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_gresp_storage_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_gresp_storage_to_litr1c_13', data=cptr%cc13f%m_gresp_storage_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_gresp_xfer_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_gresp_xfer_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_gresp_xfer_to_litr1c_13', data=cptr%cc13f%m_gresp_xfer_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadstemc_to_cwdc_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemc_to_cwdc_fire_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemc_to_cwdc_fire_13', data=cptr%cc13f%m_deadstemc_to_cwdc_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadcrootc_to_cwdc_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootc_to_cwdc_fire_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootc_to_cwdc_fire_13', data=cptr%cc13f%m_deadcrootc_to_cwdc_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_litr1c_to_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_litr1c_to_fire_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_litr1c_to_fire_13', data=cptr%cc13f%m_litr1c_to_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_litr2c_to_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_litr2c_to_fire_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_litr2c_to_fire_13', data=cptr%cc13f%m_litr2c_to_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_litr3c_to_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_litr3c_to_fire_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_litr3c_to_fire_13', data=cptr%cc13f%m_litr3c_to_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_cwdc_to_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_cwdc_to_fire_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_cwdc_to_fire_13', data=cptr%cc13f%m_cwdc_to_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! leafc_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_to_litr1c_13', data=cptr%cc13f%leafc_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! leafc_to_litr2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_to_litr2c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_to_litr2c_13', data=cptr%cc13f%leafc_to_litr2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! leafc_to_litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_to_litr3c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_to_litr3c_13', data=cptr%cc13f%leafc_to_litr3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! frootc_to_litr1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_to_litr1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_to_litr1c_13', data=cptr%cc13f%frootc_to_litr1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! frootc_to_litr2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_to_litr2c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_to_litr2c_13', data=cptr%cc13f%frootc_to_litr2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! frootc_to_litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootc_to_litr3c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootc_to_litr3c_13', data=cptr%cc13f%frootc_to_litr3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! cwdc_to_litr2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdc_to_litr2c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cwdc_to_litr2c_13', data=cptr%cc13f%cwdc_to_litr2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! cwdc_to_litr3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdc_to_litr3c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cwdc_to_litr3c_13', data=cptr%cc13f%cwdc_to_litr3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr1_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr1_hr_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr1_hr_13', data=cptr%cc13f%litr1_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr1c_to_soil1c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr1c_to_soil1c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr1c_to_soil1c_13', data=cptr%cc13f%litr1c_to_soil1c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr2_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr2_hr_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr2_hr_13', data=cptr%cc13f%litr2_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr2c_to_soil2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr2c_to_soil2c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr2c_to_soil2c_13', data=cptr%cc13f%litr2c_to_soil2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr3_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr3_hr_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr3_hr_13', data=cptr%cc13f%litr3_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr3c_to_soil3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr3c_to_soil3c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr3c_to_soil3c_13', data=cptr%cc13f%litr3c_to_soil3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil1_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil1_hr_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil1_hr_13', data=cptr%cc13f%soil1_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil1c_to_soil2c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil1c_to_soil2c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil1c_to_soil2c_13', data=cptr%cc13f%soil1c_to_soil2c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil2_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil2_hr_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil2_hr_13', data=cptr%cc13f%soil2_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil2c_to_soil3c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil2c_to_soil3c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil2c_to_soil3c_13', data=cptr%cc13f%soil2c_to_soil3c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil3_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil3_hr_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil3_hr_13', data=cptr%cc13f%soil3_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil3c_to_soil4c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil3c_to_soil4c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil3c_to_soil4c_13', data=cptr%cc13f%soil3c_to_soil4c, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil4_hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil4_hr_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil4_hr_13', data=cptr%cc13f%soil4_hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! lithr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='lithr_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='lithr_13', data=cptr%cc13f%lithr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! somhr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='somhr_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='somhr_13', data=cptr%cc13f%somhr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! hr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='hr_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='hr_13', data=cptr%cc13f%hr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sr_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sr_13', data=cptr%cc13f%sr, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! er
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='er_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='er_13', data=cptr%cc13f%er, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litfire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litfire_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litfire_13', data=cptr%cc13f%litfire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! somfire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='somfire_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='somfire_13', data=cptr%cc13f%somfire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! totfire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totfire_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='totfire_13', data=cptr%cc13f%totfire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! nep
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='nep_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='nep_13', data=cptr%cc13f%nep, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! nee
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='nee_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='nee_13', data=cptr%cc13f%nee, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! col_cinputs
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_cinputs_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='col_cinputs_13', data=cptr%cc13f%col_cinputs, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! col_coutputs
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_coutputs_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='col_coutputs_13', data=cptr%cc13f%col_coutputs, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! col_fire_closs
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_fire_closs_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='col_fire_closs_13', data=cptr%cc13f%col_fire_closs, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !--------------------------------
    ! column nitrogen flux variables
    !--------------------------------

    ! ndep_to_sminn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ndep_to_sminn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ndep_to_sminn', data=cptr%cnf%ndep_to_sminn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! nfix_to_sminn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='nfix_to_sminn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='nfix_to_sminn', data=cptr%cnf%nfix_to_sminn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafn_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafn_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafn_to_litr1n', data=cptr%cnf%m_leafn_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafn_to_litr2n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafn_to_litr2n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafn_to_litr2n', data=cptr%cnf%m_leafn_to_litr2n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafn_to_litr3n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafn_to_litr3n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafn_to_litr3n', data=cptr%cnf%m_leafn_to_litr3n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootn_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootn_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootn_to_litr1n', data=cptr%cnf%m_frootn_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootn_to_litr2n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootn_to_litr2n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootn_to_litr2n', data=cptr%cnf%m_frootn_to_litr2n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootn_to_litr3n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootn_to_litr3n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootn_to_litr3n', data=cptr%cnf%m_frootn_to_litr3n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafn_storage_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafn_storage_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafn_storage_to_litr1n', data=cptr%cnf%m_leafn_storage_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootn_storage_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootn_storage_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootn_storage_to_litr1n', data=cptr%cnf%m_frootn_storage_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livestemn_storage_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemn_storage_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemn_storage_to_litr1n', data=cptr%cnf%m_livestemn_storage_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadstemn_storage_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemn_storage_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemn_storage_to_litr1n', data=cptr%cnf%m_deadstemn_storage_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livecrootn_storage_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootn_storage_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootn_storage_to_litr1n', data=cptr%cnf%m_livecrootn_storage_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadcrootn_storage_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootn_storage_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootn_storage_to_litr1n', data=cptr%cnf%m_deadcrootn_storage_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_leafn_xfer_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_leafn_xfer_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_leafn_xfer_to_litr1n', data=cptr%cnf%m_leafn_xfer_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_frootn_xfer_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_frootn_xfer_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_frootn_xfer_to_litr1n', data=cptr%cnf%m_frootn_xfer_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livestemn_xfer_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemn_xfer_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemn_xfer_to_litr1n', data=cptr%cnf%m_livestemn_xfer_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadstemn_xfer_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemn_xfer_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemn_xfer_to_litr1n', data=cptr%cnf%m_deadstemn_xfer_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livecrootn_xfer_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootn_xfer_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootn_xfer_to_litr1n', data=cptr%cnf%m_livecrootn_xfer_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadcrootn_xfer_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootn_xfer_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootn_xfer_to_litr1n', data=cptr%cnf%m_deadcrootn_xfer_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livestemn_to_cwdn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livestemn_to_cwdn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livestemn_to_cwdn', data=cptr%cnf%m_livestemn_to_cwdn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadstemn_to_cwdn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemn_to_cwdn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemn_to_cwdn', data=cptr%cnf%m_deadstemn_to_cwdn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_livecrootn_to_cwdn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_livecrootn_to_cwdn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_livecrootn_to_cwdn', data=cptr%cnf%m_livecrootn_to_cwdn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadcrootn_to_cwdn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootn_to_cwdn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootn_to_cwdn', data=cptr%cnf%m_deadcrootn_to_cwdn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_retransn_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_retransn_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_retransn_to_litr1n', data=cptr%cnf%m_retransn_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadstemn_to_cwdn_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadstemn_to_cwdn_fire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadstemn_to_cwdn_fire', data=cptr%cnf%m_deadstemn_to_cwdn_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_deadcrootn_to_cwdn_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_deadcrootn_to_cwdn_fire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_deadcrootn_to_cwdn_fire', data=cptr%cnf%m_deadcrootn_to_cwdn_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_litr1n_to_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_litr1n_to_fire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_litr1n_to_fire', data=cptr%cnf%m_litr1n_to_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_litr2n_to_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_litr2n_to_fire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_litr2n_to_fire', data=cptr%cnf%m_litr2n_to_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_litr3n_to_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_litr3n_to_fire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_litr3n_to_fire', data=cptr%cnf%m_litr3n_to_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! m_cwdn_to_fire
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='m_cwdn_to_fire', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='m_cwdn_to_fire', data=cptr%cnf%m_cwdn_to_fire, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! leafn_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafn_to_litr1n', data=cptr%cnf%leafn_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! leafn_to_litr2n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_to_litr2n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafn_to_litr2n', data=cptr%cnf%leafn_to_litr2n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! leafn_to_litr3n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_to_litr3n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafn_to_litr3n', data=cptr%cnf%leafn_to_litr3n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! frootn_to_litr1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn_to_litr1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootn_to_litr1n', data=cptr%cnf%frootn_to_litr1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! frootn_to_litr2n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn_to_litr2n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootn_to_litr2n', data=cptr%cnf%frootn_to_litr2n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! frootn_to_litr3n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn_to_litr3n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='frootn_to_litr3n', data=cptr%cnf%frootn_to_litr3n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! cwdn_to_litr2n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdn_to_litr2n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cwdn_to_litr2n', data=cptr%cnf%cwdn_to_litr2n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! cwdn_to_litr3n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdn_to_litr3n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cwdn_to_litr3n', data=cptr%cnf%cwdn_to_litr3n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr1n_to_soil1n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr1n_to_soil1n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr1n_to_soil1n', data=cptr%cnf%litr1n_to_soil1n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_soil1n_l1
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_soil1n_l1', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_soil1n_l1', data=cptr%cnf%sminn_to_soil1n_l1, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr2n_to_soil2n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr2n_to_soil2n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr2n_to_soil2n', data=cptr%cnf%litr2n_to_soil2n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_soil2n_l2
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_soil2n_l2', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_soil2n_l2', data=cptr%cnf%sminn_to_soil2n_l2, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! litr3n_to_soil3n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='litr3n_to_soil3n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='litr3n_to_soil3n', data=cptr%cnf%litr3n_to_soil3n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_soil3n_l3
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_soil3n_l3', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_soil3n_l3', data=cptr%cnf%sminn_to_soil3n_l3, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil1n_to_soil2n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil1n_to_soil2n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil1n_to_soil2n', data=cptr%cnf%soil1n_to_soil2n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_soil2n_s1
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_soil2n_s1', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_soil2n_s1', data=cptr%cnf%sminn_to_soil2n_s1, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil2n_to_soil3n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil2n_to_soil3n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil2n_to_soil3n', data=cptr%cnf%soil2n_to_soil3n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_soil3n_s2
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_soil3n_s2', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_soil3n_s2', data=cptr%cnf%sminn_to_soil3n_s2, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil3n_to_soil4n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil3n_to_soil4n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil3n_to_soil4n', data=cptr%cnf%soil3n_to_soil4n, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_soil4n_s3
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_soil4n_s3', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_soil4n_s3', data=cptr%cnf%sminn_to_soil4n_s3, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! soil4n_to_sminn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='soil4n_to_sminn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='soil4n_to_sminn', data=cptr%cnf%soil4n_to_sminn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_denit_l1s1
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_denit_l1s1', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_denit_l1s1', data=cptr%cnf%sminn_to_denit_l1s1, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_denit_l2s2
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_denit_l2s2', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_denit_l2s2', data=cptr%cnf%sminn_to_denit_l2s2, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_denit_l3s3
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_denit_l3s3', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_denit_l3s3', data=cptr%cnf%sminn_to_denit_l3s3, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_denit_s1s2
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_denit_s1s2', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_denit_s1s2', data=cptr%cnf%sminn_to_denit_s1s2, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_denit_s2s3
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_denit_s2s3', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_denit_s2s3', data=cptr%cnf%sminn_to_denit_s2s3, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_denit_s3s4
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_denit_s3s4', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_denit_s3s4', data=cptr%cnf%sminn_to_denit_s3s4, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_denit_s4
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_denit_s4', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_denit_s4', data=cptr%cnf%sminn_to_denit_s4, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_denit_excess
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_denit_excess', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_denit_excess', data=cptr%cnf%sminn_to_denit_excess, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_leached
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_leached', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_leached', data=cptr%cnf%sminn_leached, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! potential_immob
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='potential_immob', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='potential_immob', data=cptr%cnf%potential_immob, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! actual_immob
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='actual_immob', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='actual_immob', data=cptr%cnf%actual_immob, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! sminn_to_plant
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='sminn_to_plant', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='sminn_to_plant', data=cptr%cnf%sminn_to_plant, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! supplement_to_sminn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='supplement_to_sminn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='supplement_to_sminn', data=cptr%cnf%supplement_to_sminn, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! gross_nmin
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gross_nmin', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='gross_nmin', data=cptr%cnf%gross_nmin, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! net_nmin
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='net_nmin', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='net_nmin', data=cptr%cnf%net_nmin, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! denit
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='denit', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='denit', data=cptr%cnf%denit, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! col_ninputs
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_ninputs', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='col_ninputs', data=cptr%cnf%col_ninputs, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! col_noutputs
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_noutputs', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='col_noutputs', data=cptr%cnf%col_noutputs, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! col_fire_nloss
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_fire_nloss', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='col_fire_nloss', data=cptr%cnf%col_fire_nloss, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !--------------------------------
    ! gridcell a2l variables
    !--------------------------------

    ! gricell forc_hgt_u
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='forc_hgt_u', xtype=nf_double,  &
            dim1name='gridcell',long_name='observational height of wind',units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='forc_hgt_u', data=clm_a2l%forc_hgt_u, &
            dim1name='gridcell', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

#if (defined EXIT_SPINUP)
    if (flag == 'read') then
       m = 10.0_r8
       do c = begc, endc
          clm3%g%l%c%ccs%cwdc(c)   = clm3%g%l%c%ccs%cwdc(c)   * m
          clm3%g%l%c%ccs%litr1c(c) = clm3%g%l%c%ccs%litr1c(c) * m
          clm3%g%l%c%ccs%litr2c(c) = clm3%g%l%c%ccs%litr2c(c) * m
          clm3%g%l%c%ccs%litr3c(c) = clm3%g%l%c%ccs%litr3c(c) * m
          clm3%g%l%c%ccs%soil1c(c) = clm3%g%l%c%ccs%soil1c(c) * m
          clm3%g%l%c%ccs%soil2c(c) = clm3%g%l%c%ccs%soil2c(c) * m
          clm3%g%l%c%ccs%soil3c(c) = clm3%g%l%c%ccs%soil3c(c) * m
          clm3%g%l%c%ccs%soil4c(c) = clm3%g%l%c%ccs%soil4c(c) * m
          
          ! adding code for 13C, 12/25/05, PET 
          clm3%g%l%c%cc13s%cwdc(c)   = clm3%g%l%c%cc13s%cwdc(c)   * m
          clm3%g%l%c%cc13s%litr1c(c) = clm3%g%l%c%cc13s%litr1c(c) * m
          clm3%g%l%c%cc13s%litr2c(c) = clm3%g%l%c%cc13s%litr2c(c) * m
          clm3%g%l%c%cc13s%litr3c(c) = clm3%g%l%c%cc13s%litr3c(c) * m
          clm3%g%l%c%cc13s%soil1c(c) = clm3%g%l%c%cc13s%soil1c(c) * m
          clm3%g%l%c%cc13s%soil2c(c) = clm3%g%l%c%cc13s%soil2c(c) * m
          clm3%g%l%c%cc13s%soil3c(c) = clm3%g%l%c%cc13s%soil3c(c) * m
          clm3%g%l%c%cc13s%soil4c(c) = clm3%g%l%c%cc13s%soil4c(c) * m
          
          clm3%g%l%c%cns%cwdn(c)   = clm3%g%l%c%cns%cwdn(c)   * m
          clm3%g%l%c%cns%litr1n(c) = clm3%g%l%c%cns%litr1n(c) * m
          clm3%g%l%c%cns%litr2n(c) = clm3%g%l%c%cns%litr2n(c) * m
          clm3%g%l%c%cns%litr3n(c) = clm3%g%l%c%cns%litr3n(c) * m
          clm3%g%l%c%cns%soil1n(c) = clm3%g%l%c%cns%soil1n(c) * m
          clm3%g%l%c%cns%soil2n(c) = clm3%g%l%c%cns%soil2n(c) * m
          clm3%g%l%c%cns%soil3n(c) = clm3%g%l%c%cns%soil3n(c) * m
          clm3%g%l%c%cns%soil4n(c) = clm3%g%l%c%cns%soil4n(c) * m
       end do
    end if
#endif

  end subroutine CNRest

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

#endif

end module CNrestMod

