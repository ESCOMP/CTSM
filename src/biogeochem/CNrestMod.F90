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
!
! !LOCAL VARIABLES:
!EOP
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
    integer , pointer :: iptemp(:) ! pointer to memory to be allocated
    integer :: ier                 ! error status
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Determine necessary subgrid bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    !--------------------------------
    ! pft ecophysiological variables 
    !--------------------------------
    
    ! dormant_flag
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='dormant_flag', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='dormant_flag', data=pptr%pepv%dormant_flag, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

#if (defined C13)
    ! xsmrpool_c13ratio
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='xsmrpool_c13ratio', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='xsmrpool_c13ratio', data=pptr%pepv%xsmrpool_c13ratio, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
#endif

    ! alloc_pnow
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='alloc_pnow', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='alloc_pnow', data=pptr%pepv%alloc_pnow, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! tempsum_potential_gpp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempsum_potential_gpp', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='tempsum_potential_gpp', data=pptr%pepv%tempsum_potential_gpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    !annsum_potential_gpp 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_potential_gpp', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='annsum_potential_gpp', data=pptr%pepv%annsum_potential_gpp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! tempmax_retransn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempmax_retransn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='tempmax_retransn', data=pptr%pepv%tempmax_retransn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! annmax_retransn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annmax_retransn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='annmax_retransn', data=pptr%pepv%annmax_retransn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

#if (defined C13)
    ! rc13_canair
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='rc13_canair', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='rc13_canair', data=pptr%pepv%rc13_canair, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
#endif

    !--------------------------------
    ! pft carbon state variables 
    !--------------------------------

    ! leafc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc', data=pptr%pcs%leafc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

#if (defined C13)
    !--------------------------------
    ! C13 pft carbon state variables 
    !--------------------------------
    
    ! leafc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafc_13', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafc_13', data=pptr%pc13s%leafc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
#endif

    !--------------------------------
    ! pft nitrogen state variables
    !--------------------------------
    
    ! leafn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafn', data=pptr%pns%leafn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    !--------------------------------
    ! column physical state variables
    !--------------------------------
    
    ! decl
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='decl', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='decl', data=cptr%cps%decl, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! cannavg_t2m
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cannavg_t2m', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cannavg_t2m', data=cptr%cps%cannavg_t2m, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! seedc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='seedc', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='seedc', data=cptr%ccs%seedc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! prod10c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod10c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='prod10c', data=cptr%ccs%prod10c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! prod100c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod100c', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='prod100c', data=cptr%ccs%prod100c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
#if (defined C13)
    !--------------------------------
    ! C13 column carbon state variables
    !--------------------------------

    ! cwdc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdc_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cwdc_13', data=cptr%cc13s%cwdc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! seedc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='seedc_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='seedc_13', data=cptr%cc13s%seedc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! prod10c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod10c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='prod10c_13', data=cptr%cc13s%prod10c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! prod100c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod100c_13', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='prod100c_13', data=cptr%cc13s%prod100c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
#endif
    
    !--------------------------------
    ! column nitrogen state variables
    !--------------------------------
    
    ! cwdn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cwdn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='cwdn', data=cptr%cns%cwdn, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
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
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! seedn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='seedn', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='seedn', data=cptr%cns%seedn, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! prod10n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod10n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='prod10n', data=cptr%cns%prod10n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
    ! prod100n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod100n', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='prod100n', data=cptr%cns%prod100n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if
    
#if (defined EXIT_SPINUP)
    if (flag == 'read') then
       m = 20._r8
       do c = begc, endc
          clm3%g%l%c%ccs%soil1c(c) = clm3%g%l%c%ccs%soil1c(c) * m
          clm3%g%l%c%ccs%soil2c(c) = clm3%g%l%c%ccs%soil2c(c) * m
          clm3%g%l%c%ccs%soil3c(c) = clm3%g%l%c%ccs%soil3c(c) * m
          clm3%g%l%c%ccs%soil4c(c) = clm3%g%l%c%ccs%soil4c(c) * m
          
#if (defined C13)
          ! adding code for 13C, 12/25/05, PET 
          clm3%g%l%c%cc13s%soil1c(c) = clm3%g%l%c%cc13s%soil1c(c) * m
          clm3%g%l%c%cc13s%soil2c(c) = clm3%g%l%c%cc13s%soil2c(c) * m
          clm3%g%l%c%cc13s%soil3c(c) = clm3%g%l%c%cc13s%soil3c(c) * m
          clm3%g%l%c%cc13s%soil4c(c) = clm3%g%l%c%cc13s%soil4c(c) * m
#endif
          
          clm3%g%l%c%cns%soil1n(c) = clm3%g%l%c%cns%soil1n(c) * m
          clm3%g%l%c%cns%soil2n(c) = clm3%g%l%c%cns%soil2n(c) * m
          clm3%g%l%c%cns%soil3n(c) = clm3%g%l%c%cns%soil3n(c) * m
          clm3%g%l%c%cns%soil4n(c) = clm3%g%l%c%cns%soil4n(c) * m
       end do
    end if
#endif

#if (defined CNDV)
    ! pft type dgvm physical state - crownarea
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CROWNAREA', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='CROWNAREA', data=pptr%pdgvs%crownarea, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! tempsum_litfall
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempsum_litfall', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='tempsum_litfall', data=pptr%pepv%tempsum_litfall, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! annsum_litfall
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_litfall', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='annsum_litfall', data=pptr%pepv%annsum_litfall, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! nind
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='nind', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='nind', data=pptr%pdgvs%nind, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! fpcgrid
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fpcgrid', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='fpcgrid', data=pptr%pdgvs%fpcgrid, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! fpcgridold
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fpcgridold', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='fpcgridold', data=pptr%pdgvs%fpcgridold, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! gridcell type dgvm physical state - tmomin20
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='TMOMIN20', xtype=nf_double,  &
            dim1name='gridcell',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='TMOMIN20', data=gptr%gdgvs%tmomin20, &
            dim1name=nameg, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! gridcell type dgvm physical state - agdd20
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='AGDD20', xtype=nf_double,  &
            dim1name='gridcell',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='AGDD20', data=gptr%gdgvs%agdd20, &
            dim1name=nameg, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! pft type dgvm physical state - t_mo_min
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_MO_MIN', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_MO_MIN', data=pptr%pdgvs%t_mo_min, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! present
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='present', xtype=nf_int,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       allocate (iptemp(begp:endp), stat=ier)
       if (ier /= 0) then
          write(6,*) 'CNrest: allocation error '; call endrun()
       end if
       if (flag == 'write') then
          do p = begp,endp
             iptemp(p) = 0
             if (pptr%pdgvs%present(p)) iptemp(p) = 1
          end do
       end if
       call ncd_iolocal(varname='present', data=iptemp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read') then
          if (.not. readvar) then
             if (is_restart()) call endrun
          else
             do p = begp,endp
                pptr%pdgvs%present(p) = .false.
                if (iptemp(p) == 1) pptr%pdgvs%present(p) = .true.
             end do
          end if
       end if
       deallocate (iptemp)
    end if

    ! leafcmax
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafcmax', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='leafcmax', data=pptr%pcs%leafcmax, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! heatstress
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='heatstress', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='heatstress', data=pptr%pdgvs%heatstress, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
    end if

    ! greffic
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='greffic', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='greffic', data=pptr%pdgvs%greffic, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun
       end if
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

