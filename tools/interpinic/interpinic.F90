module interpinic

  !----------------------------------------------------------------------- 
  ! Interpolate initial conditions file from one resolution and/or landmask
  ! to anoterh resolution and/or landmask
  !----------------------------------------------------------------------- 

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

  ! Variables read in from input and output initial files

  integer :: nlevsno        ! maximum number of snow levels
  integer :: nlevsoi        ! number of soil levels
  integer :: nlevlak        ! number of lake levels
  integer :: nlevtot        ! number of soil and snow levels

  integer :: numcols        ! input file number of columns 
  integer :: numcolso       ! output file number of columns 

  integer :: numpfts        ! input file number of pfts 
  integer :: numpftso       ! output file number of pfts 

  ! Variable needed in various subroutines of this module

  logical :: do_dgvmo       !true => output initial file has dgvm data

  ! RTM river routing model

  integer, parameter :: rtmlon = 720  ! # of rtm longitudes
  integer, parameter :: rtmlat = 360  ! # of rtm latitudes
  real(r8) :: volr(rtmlon,rtmlat)     ! water volume in cell (m^3)

  ! Parameters

  real(r8), parameter :: spval = 1.e36 ! special value for missing data (ocean)

  ! Public methods

  public :: interp_filei

  ! Private methods

  private :: interp_ml_real
  private :: interp_sl_real
  private :: interp_sl_int

  SAVE

contains

  !=======================================================================

  subroutine interp_filei (fin, fout, cmdline)

    !----------------------------------------------------------------------- 
    ! Read initial data from netCDF instantaneous initial data history file 
    !-----------------------------------------------------------------------

    implicit none
    include 'netcdf.inc'

    ! ------------------------ arguments------ -----------------------------
    character(len=*), intent(in) :: fin   !input  initial dataset
    character(len=*), intent(in) :: fout   !output initial dataset
    character(len=*), intent(in) :: cmdline    !command line arguments
    ! --------------------------------------------------------------------

    ! ------------------------ local variables -----------------------------
    integer :: i,j,k,l,m,n         !loop indices    
    integer :: ncidi               !netCDF dataset id
    integer :: ncido               !output net
    integer :: dimid               !netCDF dimension id 
    integer :: varid               !netCDF variable id
    integer :: dimlen              !input dimension length       
    integer :: ret                 !netcdf return code
    logical :: do_rtmi             !true => input initial file has rtm data
    logical :: do_rtmo             !true => output initial file has rtm data
    logical :: do_dgvmi            !true => input initial file has dgvm data
    ! --------------------------------------------------------------------

    write (6,*) 'Mapping clm initial data from input to output initial files'

    ! Open input and output initial conditions files

    call wrap_open (fin , NF_NOWRITE, ncidi)
    call wrap_open (fout, NF_WRITE,   ncido)

    call addglobal (ncido, cmdline)

    ! Determine if input initial file has RTM or DGVM data

    ret = nf_inq_varid (ncidi, 'RTMVOLR', varid)
    if (ret == NF_NOERR) then
       do_rtmi = .true.
    else
       do_rtmi = .false.
    endif

    ret = nf_inq_varid (ncidi, 'ITYPVEG', varid)
    if (ret == NF_NOERR) then
       do_dgvmi = .true.
    else
       do_dgvmi = .false.
    endif

    ! Determine if output initial file has RTM or DGVM data

    ret = nf_inq_varid (ncido, 'RTMVOLR', varid)
    if (ret == NF_NOERR) then
       if (.not. do_rtmi) then 
          write(*,*) 'input initial file has no rtm data for the output file'
          stop
       end if
       do_rtmo = .true.
    else
       do_rtmo = .false.
    endif

    ret = nf_inq_varid (ncido, 'ITYPVEG', varid)
    if (ret == NF_NOERR) then
       if (.not. do_dgvmi) then 
          write(*,*) 'input initial file has no dgvm data for the output file'
          stop
       end if
       do_dgvmo = .true.
    else
       do_dgvmo = .false.
    endif

    call wrap_inq_dimid (ncidi, 'column', dimid)
    call wrap_inq_dimlen(ncidi, dimid, numcols)
    call wrap_inq_dimid (ncido, 'column', dimid)
    call wrap_inq_dimlen(ncido, dimid, numcolso)
    write (6,*) 'input numcols = ',numcols,' output numcols = ',numcolso

    call wrap_inq_dimid (ncidi, 'pft', dimid)
    call wrap_inq_dimlen(ncidi, dimid, numpfts)
    call wrap_inq_dimid (ncido, 'pft', dimid)
    call wrap_inq_dimlen(ncido, dimid, numpftso)
    write (6,*) 'input numpfts = ',numpfts,' output numpfts = ',numpftso

    call wrap_inq_dimid (ncidi, 'levsno', dimid)
    call wrap_inq_dimlen(ncidi, dimid, nlevsno)
    call wrap_inq_dimid (ncido, 'levsno', dimid)
    call wrap_inq_dimlen(ncido, dimid, dimlen)
    if (dimlen/=nlevsno) then
       write (6,*) 'error: input and output nlevsno values disagree'
       write (6,*) 'input nlevsno = ',nlevsno,' output nlevsno = ',dimlen
       stop
    end if

    call wrap_inq_dimid (ncidi, 'levsoi', dimid)
    call wrap_inq_dimlen(ncidi, dimid, nlevsoi)
    call wrap_inq_dimid (ncido, 'levsoi', dimid)
    call wrap_inq_dimlen(ncido, dimid, dimlen)
    if (dimlen/=nlevsoi) then
       write (6,*) 'error: input and output nlevsoi values disagree'
       write (6,*) 'input nlevsoi = ',nlevsoi,' output nlevsoi = ',dimlen
       stop
    end if

    call wrap_inq_dimid (ncidi, 'levlak', dimid)
    call wrap_inq_dimlen(ncidi, dimid, nlevlak)
    call wrap_inq_dimid (ncido, 'levlak', dimid)
    call wrap_inq_dimlen(ncido, dimid, dimlen)
    if (dimlen/=nlevlak) then
       write (6,*) 'error: input and output nlevlak values disagree'
       write (6,*) 'input nlevlak = ',nlevlak,' output nlevlak = ',dimlen
       stop
    end if

    call wrap_inq_dimid (ncidi, 'levtot', dimid)
    call wrap_inq_dimlen(ncidi, dimid, nlevtot)
    call wrap_inq_dimid (ncido, 'levtot', dimid)
    call wrap_inq_dimlen(ncido, dimid, dimlen)
    if (dimlen/=nlevtot) then
       write (6,*) 'error: input and output nlevtot values disagree'
       write (6,*) 'input nlevtot = ',nlevtot,' output nlevtot = ',dimlen
       stop
    end if

    if (do_rtmi) then
       call wrap_inq_dimid (ncidi, 'rtmlon', dimid)
       call wrap_inq_dimlen(ncidi, dimid, dimlen)
       if (dimlen/=rtmlon) then
          write (6,*) 'error: input rtmlon does not equal ',rtmlon; stop
       end if
       call wrap_inq_dimid (ncidi, 'rtmlat', dimid)
       call wrap_inq_dimlen (ncidi, dimid, dimlen)
       if (dimlen/=rtmlat) then
          write (6,*) 'error: input rtmlat does not equal ',rtmlat; stop
       end if
    end if
    if (do_rtmo) then
       call wrap_inq_dimid (ncido, 'rtmlon', dimid)
       call wrap_inq_dimlen(ncido, dimid, dimlen)
       if (dimlen/=rtmlon) then
          write (6,*) 'error: output rtmlon does not equal ',rtmlon; stop
       end if
       call wrap_inq_dimid (ncido, 'rtmlat', dimid)
       call wrap_inq_dimlen (ncido, dimid, dimlen)
       if (dimlen/=rtmlat) then
          write (6,*) 'error: output rtmlat does not equal ',rtmlat; stop
       end if
    endif

    ! Read input initial data and write output initial data
    ! Only examing the snow interfaces above zi=0 => zisno and zsno have
    ! the same level dimension below

    call interp_ml_real('ZISNO', ncidi, 'ZISNO', ncido, nlev=nlevsno, nvec=numcols, nveco=numcolso)

    call interp_ml_real('ZSNO', ncidi, 'ZSNO', ncido, nlev=nlevsno, nvec=numcols, nveco=numcolso)

    call interp_ml_real('DZSNO', ncidi, 'DZSNO', ncido, nlev=nlevsno, nvec=numcols, nveco=numcolso)

    call interp_ml_real('H2OSOI_LIQ', ncidi, 'H2OSOI_LIQ', ncido, nlev=nlevtot, nvec=numcols, nveco=numcolso)
                
    call interp_ml_real('H2OSOI_ICE', ncidi, 'H2OSOI_ICE', ncido, nlev=nlevtot, nvec=numcols, nveco=numcolso)
                
    call interp_ml_real('T_SOISNO', ncidi, 'T_SOISNO', ncido, nlev=nlevtot, nvec=numcols, nveco=numcolso)
                
    call interp_ml_real('T_LAKE', ncidi, 'T_LAKE', ncido, nlev=nlevlak, nvec=numcols, nveco=numcolso)
                
    call interp_sl_real('WA', ncidi, 'WA', ncido, nvec=numcols, nveco=numcolso)
    call interp_sl_real('WT', ncidi, 'WT', ncido, nvec=numcols, nveco=numcolso)
    call interp_sl_real('ZWT', ncidi, 'ZWT', ncido, nvec=numcols, nveco=numcolso)
                
    call interp_sl_real('T_GRND', ncidi, 'T_GRND', ncido, nvec=numcols, nveco=numcolso)
                
    call interp_sl_real('T_VEG', ncidi, 'T_VEG', ncido, nvec=numpfts, nveco=numpftso)
                
    call interp_sl_real('H2OCAN', ncidi, 'H2OCAN', ncido, nvec=numpfts, nveco=numpftso)
                
    call interp_sl_real('H2OSNO', ncidi, 'H2OSNO', ncido, nvec=numcols, nveco=numcolso)

    call interp_sl_real('SNOWDP', ncidi, 'SNOWDP', ncido, nvec=numcols, nveco=numcolso)
                
    call interp_sl_real('SNOWAGE', ncidi, 'SNOWAGE', ncido, nvec=numcols, nveco=numcolso)
                
    call interp_sl_int('SNLSNO', ncidi, 'SNLSNO', ncido, nvec=numcols, nveco=numcolso)

    if (do_rtmi) then
       call wrap_inq_varid (ncidi, 'RTMVOLR', varid)
       call wrap_get_var_double (ncidi, varid, volr)
       if (do_rtmo) then
          call wrap_inq_varid (ncido, 'RTMVOLR', varid)
          call wrap_put_var_double (ncido, varid, volr)
       endif
    endif

    if (do_dgvmi .and. do_dgvmo) then
       call interp_sl_int('ITYPVEG', ncidi, 'ITYPVEG', ncido, nvec=numpfts, nveco=numpftso)

       call interp_sl_real('FPCGRID', ncidi, 'FPCGRID', ncido, nvec=numpfts, nveco=numpftso)

       call interp_sl_real ('LAI_IND', ncidi, 'LAI_IND', ncido, nvec=numpfts, nveco=numpftso)
                             
       call interp_sl_real ('CROWNAREA', ncidi, 'CROWNAREA', ncido, nvec=numpfts, nveco=numpftso)
                             
       call interp_sl_real ('LITTERAG', ncidi, 'LITTERAG', ncido, nvec=numpfts, nveco=numpftso)
                             
       call interp_sl_real ('LITTERBG', ncidi, 'LITTERBG', ncido, nvec=numpfts, nveco=numpftso)
                             
       call interp_sl_real ('CPOOL_FAST', ncidi, 'CPOOL_FAST', ncido, nvec=numpfts, nveco=numpftso)
                             
       call interp_sl_real ('CPOOL_SLOW', ncidi, 'CPOOL_SLOW', ncido, nvec=numpfts, nveco=numpftso)
                             
       call interp_sl_int ('PRESENT', ncidi, 'PRESENT', ncido, nvec=numpfts, nveco=numpftso)
                             
       call interp_sl_real ('NIND', ncidi, 'NIND', ncido, nvec=numpfts, nveco=numpftso)
                             
       call interp_sl_real ('LM_IND', ncidi, 'LM_IND', ncido, nvec=numpfts, nveco=numpftso)
                             
       call interp_sl_real ('SM_IND', ncidi, 'SM_IND', ncido, nvec=numpfts, nveco=numpftso)
                             
       call interp_sl_real ('HM_IND', ncidi, 'HM_IND', ncido, nvec=numpfts, nveco=numpftso)
                             
       call interp_sl_real ('RM_IND', ncidi, 'RM_IND', ncido, nvec=numpfts, nveco=numpftso)
                             
       call interp_sl_real ('HTOP', ncidi, 'HTOP', ncido, nvec=numpfts, nveco=numpftso)
    end if

    ! Close input and output files

    call wrap_close(ncidi)
    call wrap_close(ncido)

    write (6,*) ' Successfully mapped from input to output initial files'

  end subroutine interp_filei

  !=======================================================================

  subroutine interp_ml_real (varnamei, ncidi, varnameo, ncido, nlev, nvec, nveco)

    implicit none

    ! ------------------------ arguments ---------------------------------
    character(len=*), intent(in) :: varnamei    ! input variable name 
    integer , intent(in)  :: ncidi              ! input netCdf id
    character(len=*), intent(in) :: varnameo    ! output variable name
    integer , intent(in)  :: ncido              ! output netCDF id  
    integer , intent(in)  :: nlev               ! number of levels
    integer , intent(in)  :: nvec               ! number of points
    integer , intent(in)  :: nveco              ! number of points
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer :: n,no                             !indices
    integer :: varid                            !variable id
    integer , allocatable :: typei(:)
    integer , allocatable :: typeo(:)
    real(r8), allocatable :: lati(:)
    real(r8), allocatable :: lato(:)
    real(r8), allocatable :: loni(:)
    real(r8), allocatable :: lono(:)
    real(r8), allocatable :: wti(:)
    real(r8), allocatable :: wto(:)
    real(r8), allocatable :: rbufmli (:,:)      !input array
    real(r8), allocatable :: rbufmlo (:,:)      !output array
    real(r8), allocatable :: dist(:)
    real(r8) :: dx,dy,distmin
    integer :: count
    ! --------------------------------------------------------------------

    allocate (rbufmli(nlev,nvec))
    allocate (rbufmlo(nlev,nveco))
    allocate (typei(nvec))
    allocate (typeo(nveco))
    allocate (lati(nvec))
    allocate (lato(nveco))
    allocate (loni(nvec))
    allocate (lono(nveco))
    allocate (wti(nvec))
    allocate (wto(nveco))
    allocate (dist(nvec))

    if (nvec == numcols) then
       call wrap_inq_varid (ncidi, 'cols1d_ityplun', varid)
       call wrap_get_var_int(ncidi, varid, typei)
       call wrap_inq_varid (ncidi, 'cols1d_lon', varid)
       call wrap_get_var_double(ncidi, varid, loni)
       call wrap_inq_varid (ncidi, 'cols1d_lat', varid)
       call wrap_get_var_double(ncidi, varid, lati)
       call wrap_inq_varid (ncidi, 'cols1d_wtxy', varid)
       call wrap_get_var_double(ncidi, varid, wti)
    end if

    if (nveco == numcolso) then
       call wrap_inq_varid (ncido, 'cols1d_ityplun', varid)
       call wrap_get_var_int(ncido, varid, typeo)
       call wrap_inq_varid (ncido, 'cols1d_lon', varid)
       call wrap_get_var_double(ncido, varid, lono)
       call wrap_inq_varid (ncido, 'cols1d_lat', varid)
       call wrap_get_var_double(ncido, varid, lato)
       call wrap_inq_varid (ncido, 'cols1d_wtxy', varid)
       call wrap_get_var_double(ncido, varid, wto)
    end if

    call wrap_inq_varid (ncidi, trim(varnamei), varid)
    call wrap_get_var_double(ncidi, varid, rbufmli)

    ! initialization will be overwritten in subsequent loops

    rbufmlo(:,:) = 0.

    do no = 1, nveco
       if (wto(no)>0.) then
          dist(:) = spval !initialize before each use
          do n = 1, nvec
             if (wti(n)>0. .and. typei(n) == typeo(no)) then
                dy = 3.14/180.* abs(lato(no)-lati(n))*6.37e6
                dx = 3.14/180.* abs(lono(no)-loni(n))*6.37e6 * &
                     0.5 * (cos(lato(no)*3.14/180.)+cos(lati(n)*3.14/180.))
                dist(n) = sqrt(dx*dx + dy*dy)
             end if
          end do          !output data land loop
          distmin = minval(dist)
          if (distmin==spval) then
             write(*,*) 'distmin=',spval
             stop
          end if
          count = 0
          do n = 1, nvec
             if (dist(n)==distmin) then
                count = count + 1
                rbufmlo(:,no) = rbufmli(:,n)
             end if       !found nearest neighbor
             if (count == 1) exit
          end do          !input data land loop
          if (count < 1) then
             write(*,*) 'no data was written: subroutine interp_ml_real'
             stop
          end if
       end if             !output data with positive weight
    end do                !output data land loop

    call wrap_inq_varid (ncido, trim(varnameo), varid)
    call wrap_put_var_double(ncido, varid, rbufmlo)

    write(*,*) 'written variable ',varnameo,' to output initial file'

    deallocate(rbufmli)
    deallocate(rbufmlo)
    deallocate(typei)
    deallocate(typeo)
    deallocate(lati)
    deallocate(lato)
    deallocate(loni)
    deallocate(lono)
    deallocate(wti)
    deallocate(wto)
    deallocate(dist)

  end subroutine interp_ml_real

  !=======================================================================

  subroutine interp_sl_real (varnamei, ncidi, varnameo, ncido, nvec, nveco)

    implicit none

    ! ------------------------ arguments ---------------------------------
    character(len=*), intent(in) :: varnamei    ! input variable name 
    integer , intent(in)  :: ncidi              ! input netCdf id
    character(len=*), intent(in) :: varnameo    ! output variable name
    integer , intent(in)  :: ncido              ! output netCDF id  
    integer , intent(in)  :: nvec               ! number of points
    integer , intent(in)  :: nveco              ! number of points
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer :: i,j,k,n,no,ni,noo                !indices
    integer :: varid                       !variable id
    integer keepi1(10), keepo1(10)
    integer , allocatable :: keepo2(:)
!   integer , allocatable :: ixyo(:) !slevis: for diagnostic work
!   integer , allocatable :: jxyo(:) !slevis: for diagnostic work
    integer , allocatable :: vtypei(:)
    integer , allocatable :: vtypeo(:)
    integer , allocatable :: typei(:)
    integer , allocatable :: typeo(:)
    real(r8), allocatable :: lati(:)
    real(r8), allocatable :: lato(:)
    real(r8), allocatable :: loni(:)
    real(r8), allocatable :: lono(:)
    real(r8), allocatable :: wti(:)
    real(r8), allocatable :: wto(:)
    real(r8), allocatable :: rbufsli (:)   !input array
    real(r8), allocatable :: rbufslo (:)   !output array
    real(r8), allocatable :: dist(:)
    real(r8) :: dx,dy,distmin
    integer :: count
    ! --------------------------------------------------------------------

    allocate (rbufsli(nvec))
    allocate (rbufslo(nveco))
    allocate (keepo2(nveco))
!   allocate (ixyo(nveco)) !slevis: for diagnostic work
!   allocate (jxyo(nveco)) !slevis: for diagnostic work
    allocate (vtypei(nvec))
    allocate (vtypeo(nveco))
    allocate (typei(nvec))
    allocate (typeo(nveco))
    allocate (lati(nvec))
    allocate (lato(nveco))
    allocate (loni(nvec))
    allocate (lono(nveco))
    allocate (wti(nvec))
    allocate (wto(nveco))
    allocate (dist(nvec))

   if (nvec == numcols) then
       call wrap_inq_varid (ncidi, 'cols1d_ityplun', varid)
       call wrap_get_var_int(ncidi, varid, typei)
       call wrap_inq_varid (ncidi, 'cols1d_lon', varid)
       call wrap_get_var_double(ncidi, varid, loni)
       call wrap_inq_varid (ncidi, 'cols1d_lat', varid)
       call wrap_get_var_double(ncidi, varid, lati)
       call wrap_inq_varid (ncidi, 'cols1d_wtxy', varid)
       call wrap_get_var_double(ncidi, varid, wti)
    else if (nvec == numpfts) then
       call wrap_inq_varid (ncidi, 'pfts1d_itypveg', varid)
       call wrap_get_var_int(ncidi, varid, vtypei)
       call wrap_inq_varid (ncidi, 'pfts1d_ityplun', varid)
       call wrap_get_var_int(ncidi, varid, typei)
       call wrap_inq_varid (ncidi, 'pfts1d_lon', varid)
       call wrap_get_var_double(ncidi, varid, loni)
       call wrap_inq_varid (ncidi, 'pfts1d_lat', varid)
       call wrap_get_var_double(ncidi, varid, lati)
       call wrap_inq_varid (ncidi, 'pfts1d_wtxy', varid)
       call wrap_get_var_double(ncidi, varid, wti)
    end if

    if (nveco == numcolso) then
       call wrap_inq_varid (ncido, 'cols1d_ityplun', varid)
       call wrap_get_var_int(ncido, varid, typeo)
       call wrap_inq_varid (ncido, 'cols1d_lon', varid)
       call wrap_get_var_double(ncido, varid, lono)
       call wrap_inq_varid (ncido, 'cols1d_lat', varid)
       call wrap_get_var_double(ncido, varid, lato)
       call wrap_inq_varid (ncido, 'cols1d_wtxy', varid)
       call wrap_get_var_double(ncido, varid, wto)
    else if (nveco == numpftso) then
!      call wrap_inq_varid (ncido, 'pfts1d_ixy', varid) !slevis start
!      call wrap_get_var_int(ncido, varid, ixyo)
!      call wrap_inq_varid (ncido, 'pfts1d_jxy', varid)
!      call wrap_get_var_int(ncido, varid, jxyo)        !slevis end
       call wrap_inq_varid (ncido, 'pfts1d_itypveg', varid)
       call wrap_get_var_int(ncido, varid, vtypeo)
       call wrap_inq_varid (ncido, 'pfts1d_ityplun', varid)
       call wrap_get_var_int(ncido, varid, typeo)
       call wrap_inq_varid (ncido, 'pfts1d_lon', varid)
       call wrap_get_var_double(ncido, varid, lono)
       call wrap_inq_varid (ncido, 'pfts1d_lat', varid)
       call wrap_get_var_double(ncido, varid, lato)
       call wrap_inq_varid (ncido, 'pfts1d_wtxy', varid)
       call wrap_get_var_double(ncido, varid, wto)
    end if

    call wrap_inq_varid (ncidi, trim(varnamei), varid)
    call wrap_get_var_double(ncidi, varid, rbufsli)

    ! initialization will be overwritten in subsequent loops

    rbufslo(:) = 0.
    if (varnameo=='T_VEG') rbufslo(:) = 283.
    if (varnameo=='T_GRND') rbufslo(:) = 283.

    keepo2(:) = 0
    do no = 1, nveco
    if (keepo2(no) == 0) then
    if (.not. do_dgvmo                          .or. &
        (do_dgvmo                              .and. &
         ( nveco==numcolso                      .or. & !column variable
          (nveco==numpftso .and. vtypeo(no)>14) .or. & !pft var + crop
          (nveco==numpftso .and. typeo(no)>1)  ))) then!pft var + not veg lunit
       if (wto(no)>0.) then
          dist(:) = spval !initialize before each use
          do n = 1, nvec
             if (wti(n)>0. .and. typei(n) == typeo(no)) then
                dy = 3.14/180.* abs(lato(no)-lati(n))*6.37e6
                dx = 3.14/180.* abs(lono(no)-loni(n))*6.37e6 * &
                     0.5*(cos(lato(no)*3.14/180.)+cos(lati(n)*3.14/180.))
                dist(n) = sqrt(dx*dx + dy*dy)
             end if
          end do          !output data land loop
          distmin = minval(dist)
          if (distmin==spval) then
             write(*,*) 'distmin, no, varnameo=',spval, no, varnameo
             stop
          end if
          count = 0
          do n = 1, nvec
             if (dist(n)==distmin) then
                count = count + 1
                if (varnameo=='HTOP' .or. varnameo=='FPCGRID') then
                   rbufslo(no) = rbufslo(no)
                else
                   rbufslo(no) = rbufsli(n)
                end if
             end if       ! found nearest neighbor
             if (count == 1) exit
          end do          !output data land loop
          if (count < 1) then
             write(*,*) 'no data was written: subroutine interp_sl_real'
             stop
          end if
       end if             !output data with positive weight
    else if (nveco==numpftso .and. typeo(no)==1 .and. vtypeo(no)<15) then
       i = 0
       do noo = 1, nveco
          if (typeo(noo)==1 .and. vtypeo(noo)<15 .and. &
              lato(noo)==lato(no) .and. lono(noo)==lono(no)) then
             i = i + 1
             if (i>10) then
                write(*,*) 'i>10 WHY?'
                stop
             end if
             keepo1(i) = noo
             keepo2(noo) = 1 !use keepo2 to skip the pfts already done
          end if
       end do
       dist(:) = spval !initialize before each use
       do n = 1, nvec
          if (wti(n)>0. .and. typei(n) == 1 .and. vtypei(n)<15) then
             dy = 3.14/180.* abs(lato(no)-lati(n))*6.37e6
             dx = 3.14/180.* abs(lono(no)-loni(n))*6.37e6 * &
                  0.5*(cos(lato(no)*3.14/180.)+cos(lati(n)*3.14/180.))
             dist(n) = sqrt(dx*dx + dy*dy)
          end if
       end do          !input data land loop
       distmin = minval(dist)
       if (distmin==spval) then
          write(*,*) 'distmin=',spval
          stop
       end if
       count = 0
       do n = 1, nvec
          if (dist(n)==distmin) then
             count = count + 1
             j = 0
             do ni = 1, nvec
                if (typei(ni)==1 .and. vtypei(ni)<15 .and. lati(ni)==lati(n) .and. loni(ni)==loni(n)) then
                   j = j + 1
                   if (j>10) then
                      write(*,*) 'j>10 WHY?'
                      stop
                   end if
                   keepi1(j) = ni
                end if
             end do
             if (j/=i) then
                write(6,*) 'j/=i, where j,i=', j,i
                stop
             end if
          end if       !found nearest neighbor
          if (count == 1) exit
       end do          !input data loop
       do k = 1, i
          rbufslo(keepo1(k)) = rbufsli(keepi1(k))
       end do          !loop to write input to output
       if (count < 1) then
          write(*,*) 'subroutine interp_sl_real: no data was written'
          stop
       end if
    else
       write(*,*) 'subroutine interp_sl_real: no data written to typeo,vtypeo,no=', typeo(no),vtypeo(no),no
       stop
    end if
    end if                !check whether to skip this no value 
    end do                !output data land loop

    call wrap_inq_varid (ncido, trim(varnameo), varid)
    call wrap_put_var_double(ncido, varid, rbufslo)

    write(*,*) 'written variable ',varnameo,' to output initial file'

    deallocate(rbufsli)
    deallocate(rbufslo)
    deallocate(keepo2)
!   deallocate(ixyo) !slevis
!   deallocate(jxyo) !slevis
    deallocate(vtypei)
    deallocate(vtypeo)
    deallocate(typei)
    deallocate(typeo)
    deallocate(lati)
    deallocate(lato)
    deallocate(loni)
    deallocate(lono)
    deallocate(wti)
    deallocate(wto)
    deallocate(dist)

  end subroutine interp_sl_real

  !=======================================================================

  subroutine interp_sl_int (varnamei, ncidi, varnameo, ncido, nvec, nveco)

    implicit none
    include 'netcdf.inc'

    ! ------------------------ arguments ---------------------------------
    character(len=*), intent(in) :: varnamei
    integer , intent(in)  :: ncidi
    character(len=*), intent(in) :: varnameo
    integer , intent(in)  :: ncido
    integer , intent(in)  :: nvec               ! number of points
    integer , intent(in)  :: nveco              ! number of points
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer :: i,j,k,n,no,ni,noo                !indices
    integer :: varid                            !variable id
    integer keepi1(10), keepo1(10)
    integer , allocatable :: keepo2(:)
    integer , allocatable :: vtypei(:)
    integer , allocatable :: vtypeo(:)
    integer , allocatable :: typei(:)
    integer , allocatable :: typeo(:)
    real(r8), allocatable :: lati(:)
    real(r8), allocatable :: lato(:)
    real(r8), allocatable :: loni(:)
    real(r8), allocatable :: lono(:)
    real(r8), allocatable :: wti(:)
    real(r8), allocatable :: wto(:)
    integer , allocatable :: ibufsli (:)        !input array
    integer , allocatable :: ibufslo (:)        !output array
    real(r8), allocatable :: dist(:)
    real(r8) :: dx,dy,distmin
    integer :: count
    ! --------------------------------------------------------------------

    allocate (ibufsli(nvec))
    allocate (ibufslo(nveco))
    allocate (keepo2(nveco))
    allocate (vtypei(nvec))
    allocate (vtypeo(nveco))
    allocate (typei(nvec))
    allocate (typeo(nveco))
    allocate (lati(nvec))
    allocate (lato(nveco))
    allocate (loni(nvec))
    allocate (lono(nveco))
    allocate (wti(nvec))
    allocate (wto(nveco))
    allocate (dist(nvec))

    if (nvec == numcols) then
       call wrap_inq_varid (ncidi, 'cols1d_ityplun', varid)
       call wrap_get_var_int(ncidi, varid, typei)
       call wrap_inq_varid (ncidi, 'cols1d_lon', varid)
       call wrap_get_var_double(ncidi, varid, loni)
       call wrap_inq_varid (ncidi, 'cols1d_lat', varid)
       call wrap_get_var_double(ncidi, varid, lati)
       call wrap_inq_varid (ncidi, 'cols1d_wtxy', varid)
       call wrap_get_var_double(ncidi, varid, wti)
    else if (nvec == numpfts) then
       call wrap_inq_varid (ncidi, 'pfts1d_itypveg', varid)
       call wrap_get_var_int(ncidi, varid, vtypei)
       call wrap_inq_varid (ncidi, 'pfts1d_ityplun', varid)
       call wrap_get_var_int(ncidi, varid, typei)
       call wrap_inq_varid (ncidi, 'pfts1d_lon', varid)
       call wrap_get_var_double(ncidi, varid, loni)
       call wrap_inq_varid (ncidi, 'pfts1d_lat', varid)
       call wrap_get_var_double(ncidi, varid, lati)
       call wrap_inq_varid (ncidi, 'pfts1d_wtxy', varid)
       call wrap_get_var_double(ncidi, varid, wti)
    end if

    if (nveco == numcolso) then
       call wrap_inq_varid (ncido, 'cols1d_ityplun', varid)
       call wrap_get_var_int(ncido, varid, typeo)
       call wrap_inq_varid (ncido, 'cols1d_lon', varid)
       call wrap_get_var_double(ncido, varid, lono)
       call wrap_inq_varid (ncido, 'cols1d_lat', varid)
       call wrap_get_var_double(ncido, varid, lato)
       call wrap_inq_varid (ncido, 'cols1d_wtxy', varid)
       call wrap_get_var_double(ncido, varid, wto)
    else if (nveco == numpftso) then
       call wrap_inq_varid (ncido, 'pfts1d_itypveg', varid)
       call wrap_get_var_int(ncido, varid, vtypeo)
       call wrap_inq_varid (ncido, 'pfts1d_ityplun', varid)
       call wrap_get_var_int(ncido, varid, typeo)
       call wrap_inq_varid (ncido, 'pfts1d_lon', varid)
       call wrap_get_var_double(ncido, varid, lono)
       call wrap_inq_varid (ncido, 'pfts1d_lat', varid)
       call wrap_get_var_double(ncido, varid, lato)
       call wrap_inq_varid (ncido, 'pfts1d_wtxy', varid)
       call wrap_get_var_double(ncido, varid, wto)
    end if

    call wrap_inq_varid (ncidi, trim(varnamei), varid)
    call wrap_get_var_int (ncidi, varid, ibufsli)

    ! initialization will be overwritten in subsequent loops

    ibufslo(:) = 0

    keepo2(:) = 0
    do no = 1, nveco
    if (keepo2(no) == 0) then
    if (.not. do_dgvmo                          .or. &
        (do_dgvmo                              .and. &
         ( nveco==numcolso                      .or. & !column variable
          (nveco==numpftso .and. vtypeo(no)>14) .or. & !pft var + crop
          (nveco==numpftso .and. typeo(no)>1)  ))) then!pft var + not veg lunit
       if (wto(no)>0.) then
          dist(:) = spval !initialize before each use
          do n = 1, nvec
             if (wti(n)>0. .and. typei(n) == typeo(no)) then
                dy = 3.14/180.* abs(lato(no)-lati(n))*6.37e6
                dx = 3.14/180.* abs(lono(no)-loni(n))*6.37e6 * &
                     0.5*(cos(lato(no)*3.14/180.)+cos(lati(n)*3.14/180.))
                dist(n) = sqrt(dx*dx + dy*dy)
             end if
          end do          !input data land loop
          distmin = minval(dist)
          if (distmin==spval) then
             write(*,*) 'distmin=',spval
             stop
          end if
          count = 0
          do n = 1, nvec
             if (dist(n)==distmin) then
                count = count + 1
                if (varnameo=='PRESENT') then
                   ibufslo(no) = 0
                else if (varnameo=='ITYPVEG') then
                   ibufslo(no) = vtypeo(no)
                else
                   ibufslo(no) = ibufsli(n)
                end if
             end if       !found nearest neighbor
             if (count == 1) exit
          end do          !input data land loop
          if (count < 1) then
             write(*,*) 'no data was written: subroutine '
             stop
          end if
       end if             !input data with positive weight
    else if (nveco==numpftso .and. typeo(no)==1 .and. vtypeo(no)<15) then
       i = 0
       do noo = 1, nveco
          if (typeo(noo)==1 .and. vtypeo(noo)<15 .and. &
              lato(noo)==lato(no) .and. lono(noo)==lono(no)) then
             i = i + 1
             if (i>10) then
                write(*,*) 'i>10 WHY?'
                stop
             end if
             keepo1(i) = noo
             keepo2(noo) = 1 !use keepo2 to skip the pfts already done
          end if
       end do
       dist(:) = spval !initialize before each use
       do n = 1, nvec
          if (wti(n)>0. .and. typei(n) == 1 .and. vtypei(n)<15) then
             dy = 3.14/180.* abs(lato(no)-lati(n))*6.37e6
             dx = 3.14/180.* abs(lono(no)-loni(n))*6.37e6 * &
                  0.5*(cos(lato(no)*3.14/180.)+cos(lati(n)*3.14/180.))
             dist(n) = sqrt(dx*dx + dy*dy)
          end if
       end do          !input data land loop
       distmin = minval(dist)
       if (distmin==spval) then
          write(*,*) 'distmin=',spval
          stop
       end if
       count = 0
       do n = 1, nvec
          if (dist(n)==distmin) then
             count = count + 1
             j = 0
             do ni = 1, nvec
                if (typei(ni)==1 .and. vtypei(ni)<15 .and. lati(ni)==lati(n) .and. loni(ni)==loni(n)) then
                   j = j + 1
                   if (j>10) then
                      write(*,*) 'j>10 WHY?'
                      stop
                   end if
                   keepi1(j) = ni
                end if
             end do
             if (j/=i) then
                write(6,*) 'j/=i, where j,i=', j,i
                stop
             end if
          end if       !found nearest neighbor
          if (count == 1) exit
       end do          !input data oop
       do k = 1, i
          ibufslo(keepo1(k)) = ibufsli(keepi1(k))
       end do          !loop to write input to output
       if (count < 1) then
          write(*,*) 'subroutine interp_sl_int: no data was written'
          stop
       end if
    else
       write(*,*) 'subroutine interp_sl_int: no data written to typeo,vtypeo,no=', typeo(no),vtypeo(no),no
       stop
    end if
    end if                !check whether to skip this no value
    end do                !input data land loop

    call wrap_inq_varid (ncido, trim(varnameo), varid)
    call wrap_put_var_int (ncido, varid, ibufslo)

    write(*,*) 'written variable ',varnameo,' to output initial file'

    deallocate (ibufsli)
    deallocate (ibufslo)
    deallocate (keepo2)
    deallocate (vtypei)
    deallocate (vtypeo)
    deallocate (typei)
    deallocate (typeo)
    deallocate (lati)
    deallocate (lato)
    deallocate (loni)
    deallocate (lono)
    deallocate (wti)
    deallocate (wto)
    deallocate (dist)

  end subroutine interp_sl_int

end module interpinic
