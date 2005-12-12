#include <misc.h>
#include <preproc.h>

module DGVMRestMod

#if (defined DGVM)
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: DGVMRestMod
!
! !DESCRIPTION:
! Read/Write to/from DGVM info to CLM restart file.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use abortutils  , only: endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: DGVMRest
!
! !REVISION HISTORY:
! Module created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: DGVMRest
!
! !INTERFACE:
  subroutine DGVMRest (ncid, flag)
!
! !DESCRIPTION:
! Read/write DGVM restart data
!
! !USES:
    use clmtype
    use ncdio
    use decompMod, only : get_proc_bounds, get_proc_global
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid            ! netcdf id
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: c,p          ! indices
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
    integer :: ier          ! error status
    logical :: readvar      ! determine if variable is on initial file
    integer , pointer :: iptemp(:) ! pointer to memory to be allocated
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
    call get_proc_global(numg, numl, numc, nump)

    ! column type physical state variable - wf
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='WF', xtype=nf_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='WF', data=cptr%cps%wf, &
            dim1name='column', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type physical state - htop
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='HTOP', xtype=nf_double,  &
            dim1name='pft',long_name='canopy top',units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='HTOP', data=pptr%pps%htop, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()
    end if

    ! pft type dgvm physical state - t_mo_min
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_MO_MIN', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T_MO_MIN', data=pptr%pdgvs%t_mo_min, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - annpsn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ANNPSN', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ANNPSN', data=pptr%pdgvs%annpsn, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - annpsnpot
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ANNPSNPOT', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ANNPSNPOT', data=pptr%pdgvs%annpsnpot, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft cflux tye - fmicr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='FMICR', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='FMICR', data=pptr%pcf%fmicr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - bm_inc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='BM_INC', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='BM_INC', data=pptr%pdgvs%bm_inc, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - afmicr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='AFMICR', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='AFMICR', data=pptr%pdgvs%afmicr, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - t10min
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T10MIN', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='T10MIN', data=pptr%pdgvs%t10min, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - tmomin20
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='TMOMIN20', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='TMOMIN20', data=pptr%pdgvs%tmomin20, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - agdd20
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='AGDD20', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='AGDD20', data=pptr%pdgvs%agdd20, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - itypveg
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ITYPVEG', xtype=nf_int,  &
            dim1name='pft', long_name='plant functional type')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='ITYPVEG', data=pptr%itype, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) call endrun()
    end if

    ! pft type dgvm physical state - fpcgrid
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='FPCGRID', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='FPCGRID', data=pptr%pdgvs%fpcgrid, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - lai_ind
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='LAI_IND', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='LAI_IND', data=pptr%pdgvs%lai_ind, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - crownarea
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CROWNAREA', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='CROWNAREA', data=pptr%pdgvs%crownarea, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - dphen
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='DPHEN', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='DPHEN', data=pptr%pdgvs%dphen, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - leafon
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='LEAFON', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='LEAFON', data=pptr%pdgvs%leafon, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - leafof
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='LEAFOF', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='LEAFOF', data=pptr%pdgvs%leafof, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - firelength
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='FIRELENGTH', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='FIRELENGTH', data=pptr%pdgvs%firelength, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - litterag
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='LITTERAG', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='LITTERAG', data=pptr%pdgvs%litterag, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - litterbg
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='LITTERBG', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='LITTERBG', data=pptr%pdgvs%litterbg, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - cpool_fast
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CPOOL_FAST', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='CPOOL_FAST', data=pptr%pdgvs%cpool_fast, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - cpool_slow
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='CPOOL_SLOW', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='CPOOL_SLOW', data=pptr%pdgvs%cpool_slow, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - k_fast_ave
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='k_fast_ave', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='k_fast_ave', data=pptr%pdgvs%k_fast_ave, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - k_slow_ave
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='K_SLOW_AVE', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='K_SLOW_AVE', data=pptr%pdgvs%k_slow_ave, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - litter_decom_ave
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='LITTER_DECOM_AVE', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='LITTER_DECOM_AVE', data=pptr%pdgvs%litter_decom_ave, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - nind
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='NIND', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='NIND', data=pptr%pdgvs%nind, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - lm_ind
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='LM_IND', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='LM_IND', data=pptr%pdgvs%lm_ind, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - sm_ind
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='SM_IND', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='SM_IND', data=pptr%pdgvs%sm_ind, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - hm_ind
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='HM_IND', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='HM_IND', data=pptr%pdgvs%hm_ind, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - rm_ind
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='RM_IND', xtype=nf_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='RM_IND', data=pptr%pdgvs%rm_ind, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
  	  if (is_restart()) call endrun
       end if	
    end if

    ! pft type dgvm physical state - present
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='PRESENT', xtype=nf_int,  &
            dim1name='pft', long_name='if pft is present is patch, 1=>yes,=>no)')
    else if (flag == 'read' .or. flag == 'write') then
       allocate (iptemp(nump), stat=ier)
       if (ier /= 0) then
          write(6,*) 'DGVMRest: allocation error '; call endrun()
       end if
       if (flag == 'write') then
!dir$ concurrent
!cdir nodep
          do p = begp,endp
             iptemp(p) = 0
             if (pptr%pdgvs%present(p)) iptemp(p) = 1
          end do
       end if
       call ncd_iolocal(varname='PRESENT', data=iptemp, &
            dim1name='pft', ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read') then
          if (.not. readvar) then
             call endrun()
          else
!dir$ concurrent
!cdir nodep
             do p = begp,endp
                pptr%pdgvs%present(p) = .false.
                if (iptemp(p) == 1) pptr%pdgvs%present(p) = .true.
             end do
          end if
       end if
       deallocate (iptemp)
    end if

  end subroutine DGVMRest

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

end module DGVMRestMod
