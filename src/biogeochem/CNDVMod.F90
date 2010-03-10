#include <misc.h>
#include <preproc.h>

module CNDVMod

#if (defined CNDV)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNDVMod
!
! !DESCRIPTION:
! Module containing routines to drive the annual dynamic vegetation
! that works with CN, reset related variables,
! and initialize/reset time invariant variables
!
! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clm_varpar          , only : maxpatch_pft, lsmlon, lsmlat, nlevsoi
  use abortutils          , only : endrun
  use CNVegStructUpdateMod, only : CNVegStructUpdate
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public dv                 ! Drives the annual dynamic vegetation that
                            ! works with CN
  public histCNDV           ! Output CNDV history file
!
! !REVISION HISTORY:
! Module modified by Sam Levis from similar module DGVMMod
! created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: dv
!
! !INTERFACE:
  subroutine dv(lbg, ubg, lbp, ubp, num_natvegp, filter_natvegp, kyr)
!
! !DESCRIPTION:
! Drives the annual dynamic vegetation that works with CN
!
! !USES:
    use clmtype
    use CNDVLightMod        , only : Light
    use CNDVEstablishmentMod, only : Establishment
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbg, ubg       ! gridcell bounds
    integer, intent(in) :: lbp, ubp       ! pft bounds
    integer, intent(inout) :: num_natvegp ! number of naturally-vegetated
                                          ! pfts in filter
    integer, intent(inout) :: filter_natvegp(ubp-lbp+1) ! filter for
                                          ! naturally-vegetated pfts
    integer, intent(in) :: kyr            ! used in routine climate20 below
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
   integer , pointer :: mxy(:)         ! pft m index (for laixy(i,j,m),etc.)
   integer , pointer :: pgridcell(:)   ! gridcell of corresponding pft
   real(r8), pointer :: fpcgrid(:)     ! foliar projective cover on gridcell (fraction)
   real(r8), pointer :: agdd(:)        ! accumulated growing degree days above 5
   real(r8), pointer :: t_mo_min(:)    ! annual min of t_mo (Kelvin)
!
! local pointers to implicit inout arguments
!
   real(r8), pointer :: tmomin20(:)         ! 20-yr running mean of tmomin
   real(r8), pointer :: agdd20(:)           ! 20-yr running mean of agdd
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: g,p                    ! indices
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    agdd20    => clm3%g%gdgvs%agdd20
    tmomin20  => clm3%g%gdgvs%tmomin20

    ! Assign local pointers to derived type members (pft-level)

    mxy       => clm3%g%l%c%p%mxy
    pgridcell => clm3%g%l%c%p%gridcell
    fpcgrid   => clm3%g%l%c%p%pdgvs%fpcgrid
    t_mo_min  => clm3%g%l%c%p%pdgvs%t_mo_min
    agdd      => clm3%g%l%c%p%pdgvs%agdd

    ! *************************************************************************
    ! S. Levis version of LPJ's routine climate20: 'Returns' tmomin20 & agdd20
    ! for use in routine bioclim, which I have placed in routine Establishment
    ! Instead of 20-yr running mean of coldest monthly temperature,
    ! use 20-yr running mean of minimum 10-day running mean
    ! *************************************************************************

    do p = lbp,ubp
       g = pgridcell(p)
       if (kyr == 2) then ! slevis: add ".and. start_type==arb_ic" here?
          tmomin20(g) = t_mo_min(p) ! NO, b/c want to be able to start dgvm
          agdd20(g) = agdd(p)       ! w/ clmi file from non-dgvm simulation
       end if
       tmomin20(g) = (19._r8 * tmomin20(g) + t_mo_min(p)) / 20._r8
       agdd20(g)   = (19._r8 * agdd20(g)   + agdd(p)    ) / 20._r8
    end do

    ! Rebuild filter of present natually-vegetated pfts after Kill()

    call BuildNatVegFilter(lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns fpcgrid and nind

    call Light(lbg, ubg, lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns updated fpcgrid, nind, crownarea, and present. Due to updated
    ! present, we do not use the natveg filter in this subroutine.

    call Establishment(lbg, ubg, lbp, ubp)

    ! Reset dgvm variables needed in next yr (too few to keep subr. dvreset)

    do p = lbp,ubp
       clm3%g%l%c%p%pcs%leafcmax(p) = 0._r8
       clm3%g%l%c%p%pdgvs%t_mo_min(p) = 1.0e+36_r8
    end do
  end subroutine dv

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: histCNDV
!
! !INTERFACE:
  subroutine histCNDV()
!
! !DESCRIPTION:
! Create CNDV history file
!
! !USES:
    use clmtype
    use ncdio
    use decompMod       , only : get_proc_bounds, get_proc_global, ldecomp
    use clm_varpar      , only : lsmlon, lsmlat, maxpatch_pft
    use domainMod       , only : ldomain,llatlon
    use clm_varctl      , only : caseid, ctitle, finidat, fsurdat, fpftcon, &
                                 frivinp_rtm
    use clm_varcon      , only : spval
    use clm_time_manager, only : get_ref_date, get_nstep, get_curr_date, &
                                 get_curr_time
    use fileutils       , only : set_filename, putfil, get_filename
    use shr_sys_mod     , only : shr_sys_getenv
    use spmdMod         , only : masterproc
    use shr_const_mod   , only : SHR_CONST_CDAY
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
   logical , pointer :: ifspecial(:)        ! true=>landunit is not vegetated (landunit-level)
   integer , pointer :: pgridcell(:)        ! gridcell index of corresponding pft (pft-level)
   integer , pointer :: plandunit(:)        ! landunit index of corresponding pft (pft-level)
   integer , pointer :: mxy(:)              ! pft m index (for laixy(i,j,m),etc.)
   real(r8), pointer :: fpcgrid(:)          ! foliar projective cover on gridcell (fraction)
   real(r8), pointer :: nind(:)             ! number of individuals (#/m**2)
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: dgvm_fn      ! dgvm history filename
    integer :: ncid                    ! dgvm netcdf file id
    integer :: omode                   ! returned mode from netCDF call
    integer :: ncprec                  ! output precision
    integer :: g,p,l                   ! indices
    integer :: begp, endp              ! per-proc beginning and ending pft indices
    integer :: begc, endc              ! per-proc beginning and ending column indices
    integer :: begl, endl              ! per-proc beginning and ending landunit indices
    integer :: begg, endg              ! per-proc gridcell ending gridcell indices
    integer :: numg,numl,numc,nump     ! global glcp cells
    integer :: ier                     ! error status
    integer :: mdcur, mscur, mcdate    ! outputs from get_curr_time
    integer :: yr,mon,day,mcsec        ! outputs from get_curr_date
    integer :: hours,minutes,secs      ! hours,minutes,seconds of hh:mm:ss
    integer :: nstep                   ! time step
    integer :: nbsec                   ! seconds components of a date
    integer :: dimid                   ! dimension, variable id
    real(r8):: time                    ! current time
    character(len=256) :: str          ! temporary string
    character(len=  8) :: curdate      ! current date
    character(len=  8) :: curtime      ! current time
    character(len= 10) :: basedate     ! base date (yyyymmdd)
    character(len=  8) :: basesec      ! base seconds
    character(len=256) :: rem_dir      ! remote (archive) directory
    character(len=256) :: rem_fn       ! remote (archive) filename
    real(r8), pointer :: rbuf2dg(:,:)  ! temporary
    integer , pointer :: ibuf2dg(:,:)  ! temporary
    character(len=32) :: subname='histCNDV'
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    ! NONE

    ! Assign local pointers to derived type members (landunit-level)

    ifspecial  => clm3%g%l%ifspecial

    ! Assign local pointers to derived subtypes components (pft-level)

    mxy       => clm3%g%l%c%p%mxy
    pgridcell => clm3%g%l%c%p%gridcell
    plandunit => clm3%g%l%c%p%landunit
    fpcgrid   => clm3%g%l%c%p%pdgvs%fpcgrid
    nind      => clm3%g%l%c%p%pdgvs%nind

    ! Determine subgrid bounds for this processor and allocate dynamic memory

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    allocate(rbuf2dg(begg:endg,maxpatch_pft), ibuf2dg(begg:endg,maxpatch_pft), stat=ier)
    if (ier /= 0) call endrun('histCNDV: allocation error for rbuf2dg, ibuf2dg')

    ! Set output precision

    ncprec = nf_double

    if (masterproc) then

       ! -----------------------------------------------------------------------
       ! Create new netCDF file. File will be in define mode
       ! Set fill mode to "no fill" to optimize performance
       ! -----------------------------------------------------------------------

       dgvm_fn = set_dgvm_filename()
       call check_ret(nf_create(dgvm_fn, nf_clobber, ncid), subname)
       call check_ret(nf_set_fill(ncid, nf_nofill, omode), subname)

       ! -----------------------------------------------------------------------
       ! Create global attributes.
       ! -----------------------------------------------------------------------

       str = 'CF1.0'
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'conventions', len_trim(str), trim(str)), subname)

       call getdatetime(curdate, curtime)
       str = 'created on ' // curdate // ' ' // curtime
       call check_ret(&
            nf_put_att_text(ncid, NF_GLOBAL,'history', len_trim(str), trim(str)), subname)

       call shr_sys_getenv('LOGNAME', str, ier)
       if (ier /= 0) call endrun('error: LOGNAME environment variable not defined')
       
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'logname',len_trim(str), trim(str)), subname)
       
       call shr_sys_getenv('HOST', str, ier)
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'host', len_trim(str), trim(str)), subname)
       
       str = 'Community Land Model: CLM3'
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'source', len_trim(str), trim(str)), subname)
       
       str = '$Name$'
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'version', len_trim(str), trim(str)), subname)
       
       str = '$Id$'
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'revision_id', len_trim(str), trim(str)), subname)

       str = ctitle
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'case_title', len_trim(str), trim(str)), subname)

       str = caseid
       call check_ret(&
            nf_put_att_text (ncid, NF_GLOBAL, 'case_id', len_trim(str), trim(str)), subname)
       
       str = get_filename(fsurdat)
       call check_ret(&
            nf_put_att_text(ncid, NF_GLOBAL, 'Surface_dataset', len_trim(str), trim(str)), subname)
       
       str = 'arbitrary initialization'
       if (finidat /= ' ') str = get_filename(finidat)
       call check_ret(&
            nf_put_att_text(ncid, NF_GLOBAL, 'Initial_conditions_dataset', len_trim(str), trim(str)), subname)
       
       str = get_filename(fpftcon)
       call check_ret(&
            nf_put_att_text(ncid, NF_GLOBAL, 'PFT_physiological_constants_dataset', len_trim(str), trim(str)), subname)
       
       if (frivinp_rtm /= ' ') then
          str = get_filename(frivinp_rtm)
          call check_ret(&
               nf_put_att_text(ncid, NF_GLOBAL, 'RTM_input_datset', len_trim(str), trim(str)), subname)
       end if

       ! -----------------------------------------------------------------------
       ! Define dimensions.
       ! -----------------------------------------------------------------------
       
       call check_ret(nf_def_dim (ncid, 'lon', lsmlon, dimid), subname)
       call check_ret(nf_def_dim (ncid, 'lat', lsmlat, dimid), subname)
       call check_ret(nf_def_dim (ncid, 'pft', maxpatch_pft, dimid), subname)
       call check_ret(nf_def_dim (ncid, 'time', nf_unlimited, dimid), subname)
       call check_ret(nf_def_dim (ncid, 'string_length', 80, dimid), subname)
       
       ! -----------------------------------------------------------------------
       ! Define variables
       ! -----------------------------------------------------------------------
       
       ! Define coordinate variables (including time)
       
       call ncd_defvar(ncid=ncid, varname='lon', xtype=ncprec, dim1name='lon', &
            long_name='coordinate longitude', units='degrees_east')
          
       call ncd_defvar(ncid=ncid, varname='lat', xtype=ncprec, dim1name='lat', &
            long_name='coordinate latitude', units='degrees_north')

       call get_curr_time(mdcur, mscur)
       call get_ref_date(yr, mon, day, nbsec)
       hours   = nbsec / 3600
       minutes = (nbsec - hours*3600) / 60
       secs    = (nbsec - hours*3600 - minutes*60)
       write(basedate,80) yr,mon,day
80     format(i4.4,'-',i2.2,'-',i2.2)
       write(basesec ,90) hours, minutes, secs
90     format(i2.2,':',i2.2,':',i2.2)
       str = 'days since ' // basedate // " " // basesec
       time = mdcur + mscur/SHR_CONST_CDAY
       
       call ncd_defvar(ncid=ncid, varname='time', xtype=nf_double, dim1name='time', &
            long_name='time', units=str)
       
       call ncd_defvar(ncid=ncid, varname='edgen', xtype=ncprec, &
            long_name='northern edge of surface grid', units='degrees_north')
          
       call ncd_defvar(ncid=ncid, varname='edgee', xtype=ncprec, &
            long_name='eastern edge of surface grid', units='degrees_east')
          
       call ncd_defvar(ncid=ncid, varname='edges', xtype=ncprec, &
            long_name='southern edge of surface grid', units='degrees_north')
          
       call ncd_defvar(ncid=ncid, varname='edgew', xtype=ncprec, &
            long_name='western edge of surface grid', units='degrees_east')

       ! Define surface grid (coordinate variables, latitude, longitude, surface type).
       
       call ncd_defvar(ncid=ncid, varname='longxy', xtype=ncprec, dim1name='lon', dim2name='lat', &
            long_name='longitude', units='degrees_east')
       
       call ncd_defvar(ncid=ncid, varname='latixy', xtype=ncprec, dim1name='lon', dim2name='lat', &
            long_name='latitude', units='degrees_north')
       
       call ncd_defvar(ncid=ncid, varname='landmask', xtype=nf_int, dim1name='lon', dim2name='lat', &
            long_name='land/ocean mask (0.=ocean and 1.=land)')
       
       ! Define time information
       
       call ncd_defvar(ncid=ncid, varname='mcdate', xtype=nf_int, dim1name='time',&
            long_name='current date (YYYYMMDD)')
       
       call ncd_defvar(ncid=ncid, varname='mcsec', xtype=nf_int, dim1name='time',&
            long_name='current seconds of current date', units='s')
       
       call ncd_defvar(ncid=ncid, varname='mdcur', xtype=nf_int, dim1name='time',&
            long_name='current day (from base day)')
       
       call ncd_defvar(ncid=ncid, varname='mscur', xtype=nf_int, dim1name='time',&
            long_name='current seconds of current day', units='s')

       call ncd_defvar(ncid=ncid, varname='nstep', xtype=nf_int, dim1name='time',&
            long_name='time step', units='s')

       ! Define time dependent variables
       
       call ncd_defvar(ncid=ncid, varname='FPCGRID', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='plant functional type cover', units='fraction of vegetated area', &
            missing_value=spval, fill_value=spval)

       call ncd_defvar(ncid=ncid, varname='NIND', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='number of individuals', units='individuals/m2 vegetated land', &
            missing_value=spval, fill_value=spval)
       
       call check_ret(nf_enddef(ncid), subname)

    end if   ! End of if-masterproc block

    ! -----------------------------------------------------------------------
    ! Write variables
    ! -----------------------------------------------------------------------

    call ncd_ioglobal(varname='edgen', data=llatlon%edges(1), ncid=ncid, flag='write')
    call ncd_ioglobal(varname='edgee', data=llatlon%edges(2), ncid=ncid, flag='write')
    call ncd_ioglobal(varname='edges', data=llatlon%edges(3), ncid=ncid, flag='write')
    call ncd_ioglobal(varname='edgew', data=llatlon%edges(4), ncid=ncid, flag='write')

    ! Write surface grid (coordinate variables, latitude, longitude, surface type).

    if (masterproc) then
       call ncd_ioglobal(varname='lon', data=llatlon%lonc, ncid=ncid, flag='write')
       call ncd_ioglobal(varname='lat', data=llatlon%latc, ncid=ncid, flag='write')
    end if
    call ncd_iolocal(varname='longxy'  , data=ldomain%lonc, ncid=ncid, &
                      flag='write', dim1name=grlnd)
    call ncd_iolocal(varname='latixy'  , data=ldomain%latc, ncid=ncid, &
                      flag='write', dim1name=grlnd)
    call ncd_iolocal(varname='landmask', data=ldomain%mask, ncid=ncid, &
                      flag='write', dim1name=grlnd)
    ! Write current date, current seconds, current day, current nstep

    call get_curr_date(yr, mon, day, mcsec)
    mcdate = yr*10000 + mon*100 + day
    nstep = get_nstep()

    call ncd_ioglobal(varname='mcdate', data=mcdate, nt=1, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='mcsec' , data=mcsec , nt=1, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='mdcur' , data=mdcur , nt=1, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='mscur' , data=mcsec , nt=1, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='nstep' , data=nstep , nt=1, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='time'  , data=time  , nt=1, ncid=ncid, flag='write')

    ! Write time dependent variables to CNDV history file

    ! The if .not. ifspecial statment below guarantees that the m index will
    ! always lie between 1 and maxpatch_pft

    rbuf2dg(:,:) = 0._r8
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(g,mxy(p)) = fpcgrid(p)*100._r8
    end do
    call ncd_iolocal(varname='FPCGRID', dim1name=grlnd, dim2name='pft', data=rbuf2dg, &
         nt=1, ncid=ncid, flag='write')

    rbuf2dg(:,:) = 0._r8
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(g,mxy(p)) = nind(p)
    end do
    call ncd_iolocal(varname='NIND',dim1name=grlnd, dim2name='pft', data=rbuf2dg, &
         nt=1, ncid=ncid, flag='write')

    ! Deallocate dynamic memory

    deallocate(rbuf2dg, ibuf2dg)

    !------------------------------------------------------------------
    ! Close and archive netcdf CNDV history file
    !------------------------------------------------------------------

    if (masterproc) then
       call check_ret(nf_close(ncid), subname)
       write(6,*)'(histCNDV): Finished writing CNDV history dataset ',&
            trim(dgvm_fn), 'at nstep = ',get_nstep()
    end if

  end subroutine histCNDV

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_dgvm_filename
!
! !INTERFACE:
  character(len=256) function set_dgvm_filename ()
!
! !DESCRIPTION:
! Determine initial dataset filenames
!
! !USES:
    use clm_varctl      , only : caseid
    use clm_time_manager, only : get_curr_date
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
!-----------------------------------------------------------------------

    call get_curr_date (yr, mon, day, sec)
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    set_dgvm_filename = "./"//trim(caseid)//".clm2.hv."//trim(cdate)//".nc"

  end function set_dgvm_filename

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BuildNatVegFilter
!
! !INTERFACE:
  subroutine BuildNatVegFilter(lbp, ubp, num_natvegp, filter_natvegp)
!
! !DESCRIPTION:
! Reconstruct a filter of naturally-vegetated PFTs for use in DGVM
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: lbp, ubp                   ! pft bounds
    integer, intent(out) :: num_natvegp                ! number of pfts in naturally-vegetated filter
    integer, intent(out) :: filter_natvegp(ubp-lbp+1)  ! pft filter for naturally-vegetated points
!
! !CALLED FROM:
! subroutine lpj in this module
!
! !REVISION HISTORY:
! Author: Forrest Hoffman
!
! !LOCAL VARIABLES:
! local pointers to implicit in arguments
    logical , pointer :: present(:)     ! whether this pft present in patch
!EOP
!
! !LOCAL VARIABLES:
    integer :: p
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (pft-level)
    present   => clm3%g%l%c%p%pdgvs%present

    num_natvegp = 0
    do p = lbp,ubp
       if (present(p)) then
          num_natvegp = num_natvegp + 1
          filter_natvegp(num_natvegp) = p
       end if
    end do

  end subroutine BuildNatVegFilter

#endif

end module CNDVMod
