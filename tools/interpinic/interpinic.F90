module interpinic

  !----------------------------------------------------------------------- 
  ! Interpolate initial conditions file from one resolution and/or landmask
  ! to another resolution and/or landmask
  !----------------------------------------------------------------------- 

  use shr_kind_mod   , only: r8 => shr_kind_r8
  use shr_const_mod  , only: SHR_CONST_PI, SHR_CONST_REARTH
  implicit none

  private

  ! Variables read in from input and output initial files

  integer :: nlevsno        ! maximum number of snow levels
  integer :: nlevsoi        ! number of soil levels
  integer :: nlevlak        ! number of lake levels
  integer :: nlevtot        ! number of soil and snow levels

  integer :: numcols        ! input file number of columns 
  integer :: numcolso       ! output file number of columns 

  integer :: numpfts        ! input file number of pfts 
  integer :: numpftso       ! output file number of pfts 

  ! RTM river routing model

  integer, parameter :: rtmlon = 720  ! # of rtm longitudes
  integer, parameter :: rtmlat = 360  ! # of rtm latitudes
  real(r8) :: volr(rtmlon,rtmlat)     ! water volume in cell (m^3)

  ! Other parameter sizes
  integer, parameter :: nlive  = 3    ! # of live pools
  integer, parameter :: npools = 12   ! # of C pools
  integer, parameter :: numrad = 2    ! # of radiation bands

  ! Parameters

  real(r8), parameter :: spval = 1.e36_r8 ! special value for missing data (ocean)
  real(r8), parameter :: converttorad = SHR_CONST_PI/180._r8
  real(r8), parameter :: radius_earth = SHR_CONST_REARTH
  integer,  parameter :: croptype     = 15
  integer,  parameter :: nonveg       = 2

  ! Public methods

  public :: interp_filei

  ! Private methods

  private :: interp_ml_real
  private :: interp_sl_real
  private :: interp_sl_int
  private :: initDistCols
  private :: findMinDistPFTs
  private :: findMinDistCols

  ! Private data
 
  real(r8), allocatable, save :: distCols(:,:)               ! Distances between old and new columns
  integer , allocatable, save :: PFT2colindxi(:)             ! PFT to column index input grid
  integer , allocatable, save :: PFT2colindxo(:)             ! PFT to column index output grid
  logical ,              save :: allPFTSfromSameGC = .false. ! Get all PFTS from the same gridcells
  logical ,              save :: noAbortIfDNE      = .false. ! Do NOT abort if some input data does not exist

  SAVE

contains

  !=======================================================================

  subroutine interp_filei (fin, fout, cmdline)

    !----------------------------------------------------------------------- 
    ! Read initial data from netCDF instantaneous initial data history file 
    !-----------------------------------------------------------------------

    use netcdf
#ifdef AIX
    use IEEE_ARITHMETIC, only: IEEE_IS_NAN
#endif

    implicit none
    include 'netcdf.inc'

    ! ------------------------ arguments------ -----------------------------
    character(len=256), intent(in) :: fin   !input  initial dataset
    character(len=256), intent(in) :: fout   !output initial dataset
    character(len=256), intent(in) :: cmdline    !command line arguments
    ! --------------------------------------------------------------------

    ! ------------------------ local variables -----------------------------
    integer :: i,j,k,l,m,n         !loop indices    
    integer :: ncidi               !netCDF dataset id
    integer :: nvecin              !input vector length
    integer :: nvars               !number of variables
    integer :: nvecout             !output vector length
    integer :: ncido               !output net
    integer :: dimid               !netCDF dimension id 
    integer :: dimidpft            !netCDF dimension id PFT
    integer :: dimidcols           !netCDF dimension id columns
    integer :: dimidlak            !netCDF dimension id lake
    integer :: dimidsoi            !netCDF dimension id soil depth
    integer :: dimidsno            !netCDF dimension id snow depth
    integer :: dimidtot            !netCDF dimension id total
    integer :: dimidrad            !netCDF dimension id numrad
    integer :: dimidrtmlat         !netCDF dimension id rtmlat
    integer :: dimidrtmlon         !netCDF dimension id rtmlon
    integer :: dimidnlive          !netCDF dimension id nlive
    integer :: dimidnpools         !netCDF dimension id npools
    integer :: varid               !netCDF variable id
    integer :: varido              !netCDF variable id
    integer :: xtype               !netCDF variable type
    integer :: ndims               !netCDF number of dimensions
    integer :: dimids(3)           !netCDF dimension ids
    integer :: dimlen              !input dimension length       
    integer :: ret                 !netcdf return code
    character(len=256) :: varname  !variable name
    real(r8), allocatable :: rbufmlo (:,:) !output array
    !--------------------------------------------------------------------

    write (6,*) 'Mapping clm initial data from input to output initial files'

    ! Open input and output initial conditions files

    ret = nf90_open(fin,  NF90_NOWRITE, ncidi )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_open(fout, NF90_WRITE,   ncido )
    if (ret/=NF90_NOERR) call handle_error (ret)

    call addglobal (ncido, cmdline)

    ret = nf90_inq_dimid(ncidi, "column", dimidcols )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncidi, dimidcols, len=numcols)
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inq_dimid(ncido, "column", dimid )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncido, dimid, len=numcolso)
    if (ret/=NF90_NOERR) call handle_error (ret)
    write (6,*) 'input numcols = ',numcols,' output numcols = ',numcolso

    ret = nf90_inq_dimid(ncidi, "pft", dimidpft )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncidi, dimidpft, len=numpfts)
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inq_dimid(ncido, "pft", dimid )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncido, dimid, len=numpftso)
    if (ret/=NF90_NOERR) call handle_error (ret)
    write (6,*) 'input numpfts = ',numpfts,' output numpfts = ',numpftso

    ret = nf90_inq_dimid(ncidi, "levsno", dimidsno )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncidi, dimidsno, len=nlevsno)
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inq_dimid(ncido, "levsno", dimid )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncido, dimid, len=dimlen)
    if (ret/=NF90_NOERR) call handle_error (ret)
    if (dimlen/=nlevsno) then
       write (6,*) 'error: input and output nlevsno values disagree'
       write (6,*) 'input nlevsno = ',nlevsno,' output nlevsno = ',dimlen
       stop
    end if

    ret = nf90_inq_dimid(ncidi, "levsoi", dimidsoi )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncidi, dimidsoi, len=nlevsoi)
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inq_dimid(ncido, "levsoi", dimid )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncido, dimid, len=dimlen)
    if (ret/=NF90_NOERR) call handle_error (ret)
    if (dimlen/=nlevsoi) then
       write (6,*) 'error: input and output nlevsoi values disagree'
       write (6,*) 'input nlevsoi = ',nlevsoi,' output nlevsoi = ',dimlen
       stop
    end if

    ret = nf90_inq_dimid(ncidi, "levlak", dimidlak )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncidi, dimidlak, len=nlevlak)
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inq_dimid(ncido, "levlak", dimid )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncido, dimid, len=dimlen)
    if (ret/=NF90_NOERR) call handle_error (ret)
    if (dimlen/=nlevlak) then
       write (6,*) 'error: input and output nlevlak values disagree'
       write (6,*) 'input nlevlak = ',nlevlak,' output nlevlak = ',dimlen
       stop
    end if

    ret = nf90_inq_dimid(ncidi, "levtot", dimidtot )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncidi, dimidtot, len=nlevtot)
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inq_dimid(ncido, "levtot", dimid )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncido, dimid, len=dimlen)
    if (ret/=NF90_NOERR) call handle_error (ret)
    if (dimlen/=nlevtot) then
       write (6,*) 'error: input and output nlevtot values disagree'
       write (6,*) 'input nlevtot = ',nlevtot,' output nlevtot = ',dimlen
       stop
    end if

    ret = nf90_inq_dimid(ncidi, "numrad", dimidrad)
    if (ret/=NF90_NOERR) call handle_error (ret)

    ! If CASA data exist on input file
    ret = nf_inq_varid (ncidi, 'livefr', varid)
    if (ret == NF_NOERR) then
       ret = nf90_inq_dimid(ncidi, "nlive", dimidnlive)
       if (ret/=NF90_NOERR) call handle_error (ret)
       ret = nf90_inquire_dimension(ncidi, dimidnlive, len=dimlen)
       if (ret/=NF90_NOERR) call handle_error (ret)
       if (dimlen/=nlive) then
          write (6,*) 'error: input nlive does not equal ',nlive; stop
       end if
       ret = nf90_inq_dimid(ncidi, "npools", dimidnpools)
       if (ret/=NF90_NOERR) call handle_error (ret)
       ret = nf90_inquire_dimension(ncidi, dimidnpools, len=dimlen)
       if (ret/=NF90_NOERR) call handle_error (ret)
       if (dimlen/=npools) then
          write (6,*) 'error: input npools does not equal ',npools; stop
       end if
    else
       dimidnpools = -1
       dimidnlive  = -1
    end if

    ! numrad dimension

    ret = nf90_inquire_dimension(ncido, dimidrad, len=dimlen)
    if (ret/=NF90_NOERR) call handle_error (ret)
    if (dimlen/=numrad) then
       write (6,*) 'error: output numrad dimension size does not equal ',numrad; stop
    end if
#ifdef AIX
    allocate( rbufmlo(numrad,numpftso) )
#endif

    ! If RTM data exists on input file
    ret = nf_inq_varid (ncidi, 'RTMVOLR', varid)
    if (ret == NF_NOERR) then
       ret = nf90_inq_dimid(ncidi, "rtmlon", dimidrtmlon)
       if (ret/=NF90_NOERR) call handle_error (ret)
       ret = nf90_inquire_dimension(ncidi, dimidrtmlon, len=dimlen)
       if (ret/=NF90_NOERR) call handle_error (ret)
       if (dimlen/=rtmlon) then
          write (6,*) 'error: input rtmlon does not equal ',rtmlon; stop
       end if
       ret = nf90_inq_dimid(ncidi, "rtmlat", dimidrtmlat)
       if (ret/=NF90_NOERR) call handle_error (ret)
       ret = nf90_inquire_dimension(ncidi, dimidrtmlat, len=dimlen)
       if (ret/=NF90_NOERR) call handle_error (ret)
       if (dimlen/=rtmlat) then
          write (6,*) 'error: input rtmlat does not equal ',rtmlat; stop
       end if
    else
       dimidrtmlat = -1
       dimidrtmlon = -1
    end if
    ! If RTM data exists on output file
    ret = nf_inq_varid (ncido, 'RTMVOLR', varid)
    if (ret == NF_NOERR) then
       ret = nf90_inq_dimid(ncido, "rtmlon", dimid)
       if (ret/=NF90_NOERR) call handle_error (ret)
       ret = nf90_inquire_dimension(ncido, dimid, len=dimlen)
       if (ret/=NF90_NOERR) call handle_error (ret)
       if (dimlen/=rtmlon) then
          write (6,*) 'error: output rtmlon does not equal ',rtmlon; stop
       end if
       ret = nf90_inq_dimid(ncido, "rtmlat", dimid)
       if (ret/=NF90_NOERR) call handle_error (ret)
       ret = nf90_inquire_dimension(ncido, dimid, len=dimlen)
       if (ret/=NF90_NOERR) call handle_error (ret)
       if (dimlen/=rtmlat) then
          write (6,*) 'error: output rtmlat does not equal ',rtmlat; stop
       end if
    endif
    !
    ! Check if DGVM data exists on the dataset and if so -- use all PFT data from same grid cell
    ! Otherwise use data for the same veg type from potentially different grid-cells
    !
    ret = nf90_inq_varid(ncidi, 'ITYPVEG', varid)
    if (ret==NF90_ENOTVAR)then
       allPFTSfromSameGC = .false.
    else if (ret/=NF90_NOERR)then
       call handle_error (ret)
    else
       allPFTSfromSameGC = .true.
    end if

    ! Initialize the column distances
    call initDistCols( ncidi, ncido )
    ! Prepare find distances to function
    call findMinDistCols( ncidi, ncido, 0, nmin=i, allocate=.true. )
    call findMinDistPFTs( ncidi, ncido, 0, nmin=i, allocate=.true. )

    ! Read input initial data and write output initial data
    ! Only examing the snow interfaces above zi=0 => zisno and zsno have
    ! the same level dimension below

    ret = nf90_inquire(ncidi, nVariables=nvars )
    if (ret/=NF90_NOERR) call handle_error (ret)
    do i = 1, nvars
       varid = i
       ret = nf_inq_varname(ncidi, varid, varname )
       if (ret/=NF_NOERR) call handle_error (ret)

       if ( index(varname,"timemgr_"           ) /= 0 ) cycle
       if ( index(varname,"PFT_"               ) /= 0 ) cycle
       if ( index(varname,"1d"                 ) /= 0 ) cycle
       if ( index(varname,"EFLX_LWRAD_OUT"     ) /= 0 ) cycle
       if ( index(varname,"FRAC_VEG_NOSNO_ALB" ) /= 0 ) cycle
       if ( index(varname,"tlai"               ) /= 0 ) cycle
       if ( index(varname,"tsai"               ) /= 0 ) cycle
       if ( index(varname,"elai"               ) /= 0 ) cycle
       if ( index(varname,"esai"               ) /= 0 ) cycle
       if ( index(varname,"T_REF"              ) /= 0 ) cycle
       if ( index(varname,"TREF"               ) /= 0 ) cycle
       if ( index(varname,"t_ref2m"            ) /= 0 ) cycle
       ret = nf90_inq_varid(ncido, varname, varido)
       if (ret==NF90_ENOTVAR)then
          cycle
       else if (ret/=NF90_NOERR)then
          call handle_error (ret)
       end if
       ret = nf90_inquire_variable( ncidi, varid, xtype=xtype, ndims=ndims, &
                                    dimids=dimids )
       if (ret/=NF90_NOERR) call handle_error (ret)
       ! For 1D variables
       if ( ndims == 1 ) then
          if ( dimids(1) == dimidcols )then
             nvecin  = numcols
             nvecout = numcolso
          else if ( dimids(1) == dimidpft )then
             nvecin  = numpfts
             nvecout = numpftso
          else
             write (6,*) 'Skip 1D variable with unknown dimension: ', trim(varname)
             cycle
          end if
          if ( xtype == NF90_INT )then
             call interp_sl_int( varname, ncidi, ncido, nvec=nvecin, nveco=nvecout )
          else if ( xtype == NF90_DOUBLE )then
             call interp_sl_real( varname, ncidi, ncido, nvec=nvecin, nveco=nvecout )
          else
             write (6,*) 'error: variable is not of type double or integer'; stop
          end if
       ! For RTM variables
       else if ( (dimids(2) == dimidrtmlat) .and. (dimids(1) == dimidrtmlon) )then
          if ( index(varname,"RTMVOLR"  ) /= 1 ) cycle  ! If anything BUT RTMVOLR -- go to next variable
          ret = nf90_inq_varid(ncidi, varname, varid)
          if (ret/=NF90_NOERR) call handle_error (ret)
          ret = nf90_get_var(ncidi, varid, volr)
          if (ret/=NF90_NOERR) call handle_error (ret)
          ret = nf90_inq_varid(ncido, varname, varid)
          if (ret==NF90_ENOTVAR)then
             cycle
          else if (ret/=NF90_NOERR)then
             call handle_error (ret)
          end if
#ifdef AIX
!$OMP PARALLEL DO PRIVATE (i,j)
          do j  = 1, rtmlat
             do l  = 1, rtmlon
                if ( IEEE_IS_NAN(volr(l,j)) ) volr(l,j) = spval
             end do
          end do
!$OMP END PARALLEL DO
#endif
          write (6,*) 'RTM variable copied over: ', trim(varname)
          ret = nf90_put_var(ncido, varid, volr)
          if (ret/=NF90_NOERR) call handle_error (ret)
       ! For 2D variables
       else if ( ndims == 2 )then
          if ( xtype /= NF90_DOUBLE )then
             write (6,*) 'error: variable is not of double type'; stop
          end if
          if ( dimids(1) == dimidrad )then
#ifdef AIX
             if ( dimids(2) == dimidpft )then
                ret = nf90_inq_varid(ncido, varname, varid)
                if (ret==NF90_ENOTVAR)then
                   cycle
                else if (ret/=NF90_NOERR)then
                   call handle_error (ret)
                end if
                ret = nf90_get_var(ncido, varid, rbufmlo)
                if (ret/=NF90_NOERR) call handle_error (ret)
!$OMP PARALLEL DO PRIVATE (n,k)
                do n  = 1, numpftso
                   do k = 1, numrad
                      if ( IEEE_IS_NAN(rbufmlo(k,n)) ) rbufmlo(k,n) = spval
                   end do
                end do
!$OMP END PARALLEL DO
                ret = nf90_put_var(ncido, varid, rbufmlo)
                if (ret/=NF90_NOERR) call handle_error (ret)
                write (6,*) 'copied and cleaned variable with numrad dimension: ', trim(varname)
             else
                write (6,*) 'skipping variable with numrad dimension: ', trim(varname)
             end if
#else
             write (6,*) 'skipping variable with numrad dimension: ', trim(varname)
#endif
             cycle
          end if
          if ( dimids(2) /= dimidcols .and. dimids(2) /= dimidpft )then
             write (6,*) 'error: variable = ', varname
             write (6,*) 'error: variables second dimension is not recognized'; stop
          end if
          if ( dimids(1) == dimidlak )then
             call interp_ml_real(varname, ncidi, ncido, &
                                 nlev=nlevlak, nvec=numcols, nveco=numcolso)
          else if ( dimids(1) == dimidtot )then
             call interp_ml_real(varname, ncidi, ncido, &
                                 nlev=nlevtot, nvec=numcols, nveco=numcolso)
          else if ( dimids(1) == dimidsno )then
             call interp_ml_real(varname, ncidi, ncido, &
                                 nlev=nlevsno, nvec=numcols, nveco=numcolso)
          else if ( dimids(1) == dimidsoi )then
             call interp_ml_real(varname, ncidi, ncido, &
                                 nlev=nlevsoi, nvec=numcols, nveco=numcolso)
          else if ( dimids(1) == dimidnlive)then
             call interp_ml_real(varname, ncidi, ncido, &
                                 nlev=nlive, nvec=numpfts, nveco=numpftso)
          else if ( dimids(1) == dimidnpools)then
             call interp_ml_real(varname, ncidi, ncido, &
                                 nlev=npools, nvec=numpfts, nveco=numpftso)
          else
             write (6,*) 'error: variable = ', varname
             write (6,*) 'error: variables first dimension is not recognized'; stop
          end if
       else
          write (6,*) 'skipping variable NOT 1 or 2D: ', trim(varname)
       end if
    end do

    ! Close input and output files

    ret = nf90_close( ncidi )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_close( ncido )
    if (ret/=NF90_NOERR) call handle_error (ret)

    write (6,*) ' Successfully mapped from input to output initial files'

  end subroutine interp_filei

  !=======================================================================

  subroutine initDistCols( ncidi, ncido )
    ! Initialize the column distances
    use netCDF
    implicit none

    ! ------------------------ arguments ---------------------------------
    integer , intent(in)  :: ncidi              ! input netCdf id
    integer , intent(in)  :: ncido              ! output netCDF id  
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    real(r8), allocatable :: lati(:)
    real(r8), allocatable :: lato(:)
    real(r8), allocatable :: loni(:)
    real(r8), allocatable :: lono(:)
    real(r8), allocatable :: latpfti(:)
    real(r8), allocatable :: latpfto(:)
    real(r8), allocatable :: lonpfti(:)
    real(r8), allocatable :: lonpfto(:)
    real(r8) :: dx,dy
    integer  :: i, n, no
    integer  :: ret                !NetCDF return code
    integer  :: varid              !netCDF variable id
    ! --------------------------------------------------------------------

    !
    ! Distances for columns
    !
    write(6,*) 'Calculating distances for initialization'

    allocate (lati(numcols))
    allocate (lato(numcolso))
    allocate (loni(numcols))
    allocate (lono(numcolso))
    call wrap_inq_varid (ncidi, 'cols1d_lon', varid)
    call wrap_get_var_double(ncidi, varid, loni)
    call wrap_inq_varid (ncidi, 'cols1d_lat', varid)
    call wrap_get_var_double(ncidi, varid, lati)

    call wrap_inq_varid (ncido, 'cols1d_lon', varid)
    call wrap_get_var_double(ncido, varid, lono)
    call wrap_inq_varid (ncido, 'cols1d_lat', varid)
    call wrap_get_var_double(ncido, varid, lato)

    allocate (distCols(numcols,numcolso))
    ! Convert to radians
    do n = 1, numcols
       lati(n) = lati(n)*converttorad
       loni(n) = loni(n)*converttorad
    end do
    do no = 1, numcolso
       lato(no) = lato(no)*converttorad
       lono(no) = lono(no)*converttorad
    end do

!$OMP PARALLEL DO PRIVATE (no,n,dx,dy)
    do no = 1, numcolso
       do n = 1, numcols
          dy = abs(lato(no)-lati(n))*radius_earth
          dx = abs(lono(no)-loni(n))*radius_earth * 0.5_r8 * (cos(lato(no))+cos(lati(n)))
          distCols(n,no) = sqrt(dx*dx + dy*dy)
       end do          !output data land loop
    end do
!$OMP END PARALLEL DO
    write(6,*) 'Distances done!'
    !
    ! Now find which PFT's match which columns
    !
    write(6,*) 'Finding which PFTs match columns!'
    allocate (PFT2colindxi(numpfts))
    allocate (PFT2colindxo(numpftso))
    allocate (latpfti(numpfts))
    allocate (latpfto(numpftso))
    allocate (lonpfti(numpfts))
    allocate (lonpfto(numpftso))
    call wrap_inq_varid (ncidi, 'pfts1d_lon', varid)
    !if (ret/=NF90_NOERR) call handle_error (ret)
    call wrap_get_var_double(ncidi, varid, lonpfti)
    call wrap_inq_varid (ncidi, 'pfts1d_lat', varid)
    call wrap_get_var_double(ncidi, varid, latpfti)
    call wrap_inq_varid (ncido, 'pfts1d_lon', varid)
    call wrap_get_var_double(ncido, varid, lonpfto)
    call wrap_inq_varid (ncido, 'pfts1d_lat', varid)
    call wrap_get_var_double(ncido, varid, latpfto)
    ! Convert to radians
    do n = 1, numpfts
       latpfti(n) = latpfti(n)*converttorad
       lonpfti(n) = lonpfti(n)*converttorad
    end do
    do no = 1, numpftso
       latpfto(no) = latpfto(no)*converttorad
       lonpfto(no) = lonpfto(no)*converttorad
    end do
    PFT2colindxi(:) = -1
    do i = 1, numpfts
       do n = 1, numcols
          if ( (lati(n) == latpfti(i)) .and. (loni(n) == lonpfti(i)) )then
            PFT2colindxi(i) = n
            exit
          end if
       end do
       if ( PFT2colindxi(i) == -1 )then
          write(6,*) 'PFT input position does NOT match a column position: PFT # ', i
          stop
       end if
    end do
    PFT2colindxo(:) = -1
    do i = 1, numpftso
       do n = 1, numcolso
          if ( (lato(n) == latpfto(i)) .and. (lono(n) == lonpfto(i)) )then
            PFT2colindxo(i) = n
            exit
          end if
       end do
       if ( PFT2colindxo(i) == -1 )then
          write(6,*) 'PFT output position does NOT match a column position: PFT # ', i
          stop
       end if
    end do
    write(6,*) 'Done finding PFT indices!'

    deallocate(latpfti)
    deallocate(latpfto)
    deallocate(lonpfti)
    deallocate(lonpfto)

    deallocate(lati)
    deallocate(lato)
    deallocate(loni)
    deallocate(lono)
  end subroutine initDistCols

  !=======================================================================

  subroutine findMinDistCols( ncidi, ncido, no, nmin, allocate )
    ! Find the minimun column distances excluding columns of different type

    use netcdf

    implicit none

    ! ------------------------ arguments ---------------------------------
    integer , intent(in)  :: ncidi              ! input netCdf id
    integer , intent(in)  :: ncido              ! output netCDF id  
    integer , intent(in)  :: no                 ! vector number
    integer,  intent(out) :: nmin               ! index of minimum distance
    logical,  intent(in), optional :: allocate  ! if just allocating
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer , allocatable, save :: typei(:)
    integer , allocatable, save :: typeo(:)
    real(r8), allocatable, save :: wti(:)
    real(r8) :: distmin                          ! Minimum distance
    integer  :: n
    integer  :: varid                            ! netCDF variable id
    ! --------------------------------------------------------------------

    if ( present(allocate) ) then
       allocate (typei(numcols))
       allocate (typeo(numcolso))
       allocate (wti(numcols))
       call wrap_inq_varid (ncidi, 'cols1d_ityplun', varid)
       call wrap_get_var_int(ncidi, varid, typei)
       call wrap_inq_varid (ncidi, 'cols1d_wtxy', varid)
       call wrap_get_var_double(ncidi, varid, wti)
       call wrap_inq_varid (ncido, 'cols1d_ityplun', varid)
       call wrap_get_var_int(ncido, varid, typeo)
    else

       distmin = spval
       nmin    = 0
       do n = 1, numcols
          if ( (wti(n) > 0.0_r8) .and. (typei(n) == typeo(no)) ) then
             if ( distCols(n,no) < distmin )then
                distmin = distCols(n,no)
                nmin = n
             end if
          end if
       end do
       if ( distmin == spval )then
          write(*,*) 'Can not find the closest column: no = ', no
          stop
       end if
    end if

  end subroutine findMinDistCols

  !=======================================================================

  subroutine findMinDistPFTs( ncidi, ncido, no, nmin, allocate )
    ! Find the PFT distances based on the column distances already calculated

    use netcdf

    implicit none

    ! ------------------------ arguments ---------------------------------
    integer , intent(in)  :: ncidi              ! input netCdf id
    integer , intent(in)  :: ncido              ! output netCDF id  
    integer , intent(in)  :: no                 ! vector number
    integer,  intent(out) :: nmin               ! index of minimum distance
    logical,  intent(in), optional :: allocate  ! if just allocating
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer , allocatable, save :: typei(:)
    integer , allocatable, save :: typeo(:)
    integer , allocatable, save :: vtypei(:)
    integer , allocatable, save :: vtypeo(:)
    real(r8), allocatable, save :: wti(:)
    real(r8) :: distmin                          ! Minimum distance
    integer  :: n
    integer  :: ret                              ! NetCDF return code
    integer  :: varid                            ! netCDF variable id
    ! --------------------------------------------------------------------

    !
    ! Distances for PFT's to output index no
    !
    if ( present(allocate) ) then
       allocate (typei (numpfts))
       allocate (vtypei(numpfts))
       allocate (wti   (numpfts))

       allocate (typeo (numpftso))
       allocate (vtypeo(numpftso))

       call wrap_inq_varid (ncidi, 'pfts1d_ityplun', varid)
       call wrap_get_var_int(ncidi, varid, typei)

       call wrap_inq_varid (ncidi, 'pfts1d_wtxy', varid)
       call wrap_get_var_double(ncidi, varid, wti)

       ret = nf90_inq_varid( ncidi, 'pfts1d_itypveg', varid )
       if (ret/=NF90_NOERR) call handle_error (ret)
       ret = nf90_get_var( ncidi, varid, vtypei)
       if (ret/=NF90_NOERR) call handle_error (ret)

       call wrap_inq_varid (ncido, 'pfts1d_ityplun', varid)
       call wrap_get_var_int(ncido, varid, typeo)

       ret = nf90_inq_varid( ncido, 'pfts1d_itypveg', varid )
       if (ret/=NF90_NOERR) call handle_error (ret)
       ret = nf90_get_var( ncido, varid, vtypeo)
       if (ret/=NF90_NOERR) call handle_error (ret)

       if ( size(PFT2colindxi) /= numpfts  .or. (size(PFT2colindxo) /= numpftso) )then
          write(6,*) 'PFT2colindxi or PFT2colindxo out of bounds'
          stop
       end if
    else

       nmin    = 0
       distmin = spval
       do n = 1, numpfts
          if (wti(n)>0. .and. (typei(n) == typeo(no)) ) then
             if ( allPFTSfromSameGC .or. (typeo(no) >= nonveg) .or. (vtypei(n) == vtypeo(no)) )then
                if ( distCols(PFT2colindxi(n),PFT2colindxo(no)) < distmin )then
                   distmin = distCols(PFT2colindxi(n),PFT2colindxo(no))
                   nmin    = n
                end if 
             end if
          end if
       end do
       if ( distmin == spval )then
          write(*,*) 'Can not find the closest column: no = ', no
          stop
       end if
    end if

  end subroutine findMinDistPFTs

  !=======================================================================

  subroutine interp_ml_real (varname, ncidi, ncido, nlev, nvec, nveco)

    use netcdf

    implicit none
    include 'netcdf.inc'

    ! ------------------------ arguments ---------------------------------
    character(len=*), intent(in) :: varname     ! input variable name 
    integer , intent(in)  :: ncidi              ! input netCdf id
    integer , intent(in)  :: ncido              ! output netCDF id  
    integer , intent(in)  :: nlev               ! number of levels
    integer , intent(in)  :: nvec               ! number of points
    integer , intent(in)  :: nveco              ! number of points
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer :: n,no,l                           !indices
    integer :: varid                            !variable id
    real(r8), allocatable :: rbufmli (:,:)      !input array
    real(r8), allocatable :: rbufmlo (:,:)      !output array
    real(r8), allocatable :: origValues(:)      !temporary array of original values
    real(r8), allocatable :: wto(:)
    integer :: count
    integer :: ret                              !netCDF return code
    ! --------------------------------------------------------------------

    allocate (rbufmli(nlev,nvec))
    allocate (rbufmlo(nlev,nveco))
    allocate (wto(nveco))

    if (nveco == numcolso) then
       call wrap_inq_varid (ncido, 'cols1d_wtxy', varid)
       call wrap_get_var_double(ncido, varid, wto)
    else if (nveco == numpftso) then
       ret = nf90_inq_varid( ncido, 'pfts1d_wtxy', varid )
       if (ret/=NF90_NOERR) call handle_error (ret)
       ret = nf90_get_var( ncido, varid, wto)
       if (ret/=NF90_NOERR) call handle_error (ret)
    end if

    call wrap_inq_varid (ncidi, trim(varname), varid)
    call wrap_get_var_double(ncidi, varid, rbufmli)
    ret = nf90_inq_varid (ncido, trim(varname), varid)
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_get_var( ncido, varid, rbufmlo )
    if (ret/=NF90_NOERR) call handle_error (ret)

    if (nvec == numcols) then
!$OMP PARALLEL DO PRIVATE (no,n)
       do no = 1, nveco
          if (wto(no)>0._r8) then
             call findMinDistCols( ncidi, ncido, no, nmin=n )
             rbufmlo(:,no) = rbufmli(:,n)
          else
             rbufmlo(:,no) = spval
          end if             !output data with positive weight
       end do                !output data land loop
!$OMP END PARALLEL DO
    else if (nvec == numpfts) then
!$OMP PARALLEL DO PRIVATE (no,n)
       do no = 1, nveco
          if (wto(no)>0._r8) then
             call findMinDistPFTs( ncidi, ncido, no, nmin=n )
             rbufmlo(:,no) = rbufmli(:,n)
          else
             rbufmlo(:,no) = spval
          end if             !output data with positive weight
       end do                !output data land loop
!$OMP END PARALLEL DO
    else
       write(*,*) 'no data was written: subroutine interp_ml_real'
       stop
    end if

    call wrap_inq_varid (ncido, varname, varid)
    call wrap_put_var_double(ncido, varid, rbufmlo)

    write(*,*) 'written variable ',trim(varname),' to output initial file'

    deallocate(rbufmli)
    deallocate(rbufmlo)
    deallocate(wto)

  end subroutine interp_ml_real

  !=======================================================================

  subroutine interp_sl_real (varname, ncidi, ncido, nvec, nveco)

    use netcdf
#ifdef AIX
    use IEEE_ARITHMETIC, only: IEEE_IS_NAN
#endif

    implicit none
    include 'netcdf.inc'

    ! ------------------------ arguments ---------------------------------
    character(len=256), intent(in) :: varname   ! input variable name 
    integer , intent(in)  :: ncidi              ! input netCdf id
    integer , intent(in)  :: ncido              ! output netCDF id  
    integer , intent(in)  :: nvec               ! number of points
    integer , intent(in)  :: nveco              ! number of points
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer :: i,j,k,n,no,ni,noo,l,g,gg    !indices
    integer :: varid                       !variable id
    integer :: dimid                       !dimension id
    integer , allocatable :: vtypei(:)
    integer , allocatable :: vtypeo(:)
    integer , allocatable :: typeo(:)
    real(r8), allocatable :: wto(:)
    real(r8), allocatable :: wti(:)
    real(r8), allocatable :: rbufsli (:)   !input array
    real(r8), allocatable :: rbufslo (:)   !output array
    integer , allocatable :: pft_lio(:)    !PFT to land-unit index
    integer , allocatable :: pft_gio(:)    !PFT to gridcell index
    real(r8), allocatable :: lnd_wto (:)   !Land-unit weights
    real(r8), allocatable :: grd_wto (:)   !Grid-cell weights
    real(r8), allocatable :: sumwto  (:)   !sum of PFT weights
    real(r8), allocatable :: sumlwto (:)   !sum of land PFT weights
    real(r8), allocatable :: normaliz(:)   !normalization weigth
    integer :: count
    integer :: ret                         ! NetCDF return code
    integer :: num                         ! number of gridcells NOT normalized
    integer :: numgrdso                    ! number of gridcells on output grid
    integer :: numlduso                    ! number of landcell units on output grid
    logical, save         :: initialize = .false.
    logical :: htop_var    = .false.       !If variable name is == htop
    logical :: fpcgrid_var = .false.       !If variable name is == fpcgrid
    ! --------------------------------------------------------------------

    allocate (rbufsli(nvec))
    allocate (rbufslo(nveco))

    allocate (wti(nvec))
    allocate (wto(nveco))
    if (nvec == numpfts) then
       if ( trim(varname) == 'HTOP'    ) htop_var    = .true.
       if ( trim(varname) == 'FPCGRID' ) fpcgrid_var = .true.
    end if

    if (     nvec == numcols) then
       call wrap_inq_varid (ncidi, 'cols1d_wtxy', varid)
       call wrap_get_var_double(ncidi, varid, wti)
    else if (nvec == numpfts) then
       if ( htop_var .or. fpcgrid_var )then
          allocate (vtypei(nvec))
          call wrap_inq_varid (ncidi, 'pfts1d_itypveg', varid)
          call wrap_get_var_int(ncidi, varid, vtypei)
          call wrap_inq_varid (ncidi, 'pfts1d_wtxy', varid)
          call wrap_get_var_double(ncidi, varid, wti)
       end if
    end if

    if (nveco == numcolso) then
       call wrap_inq_varid (ncido, 'cols1d_wtxy', varid)
       call wrap_get_var_double(ncido, varid, wto)
    else if (nveco == numpftso) then
       ret = nf90_inq_varid( ncido, 'pfts1d_wtxy', varid )
       if (ret/=NF90_NOERR) call handle_error (ret)
       ret = nf90_get_var( ncido, varid, wto)
       if (ret/=NF90_NOERR) call handle_error (ret)
       if ( htop_var .or. fpcgrid_var )then
          allocate (vtypeo(nveco))
          allocate (typeo(nveco))
          call wrap_inq_varid (ncido, 'pfts1d_itypveg', varid)
          call wrap_get_var_int(ncido, varid, vtypeo)
          ret = nf90_inq_varid( ncido, 'pfts1d_ityplun', varid )
          if (ret/=NF90_NOERR) call handle_error (ret)
          ret = nf90_get_var( ncido, varid, typeo)
          if (ret/=NF90_NOERR) call handle_error (ret)
       end if
    end if

    call wrap_inq_varid (ncidi, varname, varid)
    call wrap_get_var_double(ncidi, varid, rbufsli)
    ret = nf90_inq_varid( ncido, varname, varid )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_get_var( ncido, varid, rbufslo )
    if (ret/=NF90_NOERR) call handle_error (ret)

    if (      nvec == numcols )then
!$OMP PARALLEL DO PRIVATE (no,n)
       do no = 1, nveco
          if (wto(no)>0._r8) then
             call findMinDistCols( ncidi, ncido, no, nmin=n )
             rbufslo(no) = rbufsli(n)
          else
             rbufslo(no) = spval
          end if             !output data with positive weight
       end do
!$OMP END PARALLEL DO
    else if ( nvec == numpfts )then
!$OMP PARALLEL DO PRIVATE (no,n)
       do no = 1, nveco
          if (wto(no)>0._r8) then
             !
             ! If variable-name is htop or fpcgrid
             !
             if ( (htop_var .or. fpcgrid_var) )then
                !
                ! AND this is crop or non-vegetated land-unit -- set to zero
                !
                if( (vtypeo(no) >= croptype) .or. (typeo(no) >= nonveg) )then
                   rbufslo(no) = 0.0_r8
                else
                   !
                   ! Otherwise calculate it from the nearest neighbor
                   !
                   call findMinDistPFTs( ncidi, ncido, no, nmin=n )
                   rbufslo(no) = rbufsli(n)
                end if
             !
             ! Otherwise calculate it from the nearest neighbor
             !
             else
                call findMinDistPFTs( ncidi, ncido, no, nmin=n )
#ifdef AIX
                if ( IEEE_IS_NAN(rbufsli(n)) ) then
                   rbufsli(n) = spval
                end if
#endif
                rbufslo(no) = rbufsli(n)
             end if          !data type
          else
             rbufslo(no) = spval
          end if             !output data with positive weight
       end do                !output data land loop
!$OMP END PARALLEL DO
       !
       ! If variable-name is fpcgrid -- renormalize so columns in a grid cell sum to 1.0
       !
       if ( fpcgrid_var )then
          write(*,*) '  normalizing variable=',trim(varname)
          ret = nf90_inq_dimid(ncido, "landunit", dimid )
          if (ret/=NF90_NOERR) call handle_error (ret)
          ret = nf90_inquire_dimension(ncido, dimid, len=numlduso)
          if (ret/=NF90_NOERR) call handle_error (ret)

          ret = nf90_inq_dimid(ncido, "gridcell", dimid )
          if (ret/=NF90_NOERR) call handle_error (ret)
          ret = nf90_inquire_dimension(ncido, dimid, len=numgrdso)
          if (ret/=NF90_NOERR) call handle_error (ret)
   
          allocate ( lnd_wto(numlduso) )
          allocate ( pft_lio(nveco) )
          allocate ( pft_gio(nveco) )
   
          if (ret/=NF90_NOERR) call handle_error (ret)
          ret = nf90_inq_varid( ncido, "pfts1d_li", varid )
          if (ret/=NF90_NOERR) call handle_error (ret)
          ret = nf90_get_var( ncido, varid, pft_lio)
          if (ret/=NF90_NOERR) call handle_error (ret)
          ret = nf90_inq_varid( ncido, "land1d_wtxy", varid )
          if (ret/=NF90_NOERR) call handle_error (ret)
          ret = nf90_get_var( ncido, varid, lnd_wto)
          if (ret/=NF90_NOERR) call handle_error (ret)
          ret = nf90_inq_varid( ncido, "pfts1d_gi", varid )
          if (ret/=NF90_NOERR) call handle_error (ret)
          ret = nf90_get_var( ncido, varid, pft_gio)
          if (ret/=NF90_NOERR) call handle_error (ret)
   
          ! Find sum of weights for each grid cell
          allocate ( sumwto  (numgrdso) )
          allocate ( sumlwto (numgrdso) )
          allocate ( normaliz(numgrdso) )
          sumwto(:numgrdso)  = 0.0_r8
          sumlwto(:numgrdso) = 0.0_r8
          do no = 1, nveco
             l = pft_lio(no)
             g = pft_gio(no)
             if ( typeo(no) == 1 ) then
                sumwto(g)  = sumwto(g)  + rbufslo(no)*lnd_wto(l)
                sumlwto(g) = sumlwto(g) + rbufslo(no)*lnd_wto(l)
             else
                sumwto(g)  = sumwto(g) + wto(no)
             end if
          end do
          num = 0
          do g = 1, numgrdso
             if ( abs(sumwto(g) - 1.0_r8) > 1.0e-6_r8 ) num = num + 1
             if ( sumwto(g) == 0.0_r8 )then
                write(*,*) 'ERROR: sumwto == 0.0: subroutine interp_sl_real'
                stop
             end if
             if ( sumlwto(g) /= 0.0_r8 )then
                normaliz(g) = ((1.0_r8 - sumlwto(g)) - sumwto(g) ) / sumlwto(g)
             else
                normaliz(g) = 1.0_r8
             end if
          end do
          write(*,*) '  % g-cells NOT norm  =', real(100*num,r8)/real(numgrdso,r8)
   
          ! Normalize each output weight by the sum of weights for each grid cell
          do no = 1, nveco
             g = pft_gio(no)
             if ( sumlwto(g) /= 0.0_r8 )then
                rbufslo(no) = rbufslo(no)*normaliz(g)
             end if
          end do
          ! Re-check
          sumwto(:numgrdso) = 0.0_r8
          do no = 1, nveco
             l = pft_lio(no)
             g = pft_gio(no)
             if ( typeo(no) == 1 ) then
                sumwto(g) = sumwto(g) + rbufslo(no)*lnd_wto(l)
             else
                sumwto(g) = sumwto(g) + wto(no)
             end if
          end do
          num = 0
          do g = 1, numgrdso
             if ( abs(sumwto(g) - 1.0_r8) > 1.0e-6_r8 )then
                write(*,*) 'ERROR: sumwto != 1.0: subroutine interp_sl_real'
                write(*,*) 'g, sum, suml, normal = ', g, sumwto(g), sumlwto(g), normaliz(g)
                do no = 1, nveco
                   l  = pft_lio(no)
                   gg = pft_gio(no)
                   if ( gg == g )then
                      write(*,*) 'no, rbufslo, lnd_wto, wto, type: ', no, rbufslo(no), lnd_wto(l), wto(no), typeo(no)
                   end if
                end do
                stop
             end if
          end do
          deallocate(lnd_wto)
          deallocate(pft_lio)
          deallocate(pft_gio)
          deallocate(sumwto)
          deallocate(sumlwto)
          deallocate(normaliz)
       end if

    else
       write(*,*) 'subroutine interp_sl_real: no data written to typeo,vtypeo,no=', &
                   typeo(no),vtypeo(no),no
       stop
    end if

    call wrap_inq_varid (ncido, varname, varid)
    call wrap_put_var_double(ncido, varid, rbufslo)

    write(*,*) 'written variable ', trim(varname),' to output initial file'

    deallocate(rbufsli)
    deallocate(rbufslo)

    if ( allocated(vtypei) ) deallocate(vtypei)
    if ( allocated(vtypeo) ) deallocate(vtypeo)
    if ( allocated(typeo) )  deallocate(typeo)
    deallocate(wto)

  end subroutine interp_sl_real

  !=======================================================================

  subroutine interp_sl_int (varname, ncidi, ncido, nvec, nveco)

    use netcdf

    implicit none
    include 'netcdf.inc'

    ! ------------------------ arguments ---------------------------------
    character(len=256), intent(in) :: varname
    integer , intent(in)  :: ncidi
    integer , intent(in)  :: ncido
    integer , intent(in)  :: nvec               ! number of points
    integer , intent(in)  :: nveco              ! number of points
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer :: i,j,k,n,no,ni,noo                !indices
    integer :: varid                            !variable id
    integer , allocatable :: vtypei(:)
    integer , allocatable :: vtypeo(:)
    integer , allocatable :: typeo(:)
    real(r8), allocatable :: wti(:)
    real(r8), allocatable :: wto(:)
    integer , allocatable :: ibufsli (:)        !input array
    integer , allocatable :: ibufslo (:)        !output array
    integer :: count
    integer :: ret                              !NetCDF return code
    logical :: present_var = .false.            !If variable name is == present
    logical :: itypveg_var = .false.            !If variable name is == itypveg
    ! --------------------------------------------------------------------

    allocate (ibufsli(nvec))
    allocate (ibufslo(nveco))
    allocate (wto(nveco))

    if (nvec == numpfts) then
       if ( trim(varname) == 'PRESENT' ) present_var = .true.
       if ( trim(varname) == 'ITYPVEG' ) itypveg_var = .true.
    end if

    if (nvec == numcols) then
       allocate (wti(nvec))
       call wrap_inq_varid (ncidi, 'cols1d_wtxy', varid)
       call wrap_get_var_double(ncidi, varid, wti)
    else if (nvec == numpfts) then
       if ( present_var .or. itypveg_var )then
          allocate (vtypei(nvec))
          call wrap_inq_varid (ncidi, 'pfts1d_itypveg', varid)
          call wrap_get_var_int(ncidi, varid, vtypei)
       end if
    end if

    if (nveco == numcolso) then
       allocate (typeo(nveco))
       call wrap_inq_varid (ncido, 'cols1d_ityplun', varid)
       call wrap_get_var_int(ncido, varid, typeo)
       call wrap_inq_varid (ncido, 'cols1d_wtxy', varid)
       call wrap_get_var_double(ncido, varid, wto)
    else if (nveco == numpftso) then
       if ( present_var .or. itypveg_var )then
          allocate (typeo(nveco))
          allocate (vtypeo(nveco))
          call wrap_inq_varid (ncido, 'pfts1d_itypveg', varid)
          call wrap_get_var_int(ncido, varid, vtypeo)
          call wrap_inq_varid (ncido, 'pfts1d_ityplun', varid)
          call wrap_get_var_int(ncido, varid, typeo)
          call wrap_inq_varid (ncido, 'pfts1d_wtxy', varid)
          call wrap_get_var_double(ncido, varid, wto)
       end if
    end if

    call wrap_inq_varid (ncidi, varname, varid)
    call wrap_get_var_int (ncidi, varid, ibufsli)

    ret = nf90_inq_varid( ncido, varname, varid )
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_get_var( ncido, varid, ibufslo )
    if (ret/=NF90_NOERR) call handle_error (ret)

    if (      nvec == numcols )then
!$OMP PARALLEL DO PRIVATE (no,n)
       do no = 1, nveco
          if (wto(no)>0._r8) then
             call findMinDistCols( ncidi, ncido, no, nmin=n )
             ibufslo(no) = ibufsli(n)
          end if             !output data with positive weight
       end do
!$OMP END PARALLEL DO
    else if ( nvec == numpfts )then
!$OMP PARALLEL DO PRIVATE (no,n) SHARED(present_var,itypveg_var)
       do no = 1, nveco
          if (wto(no)>0._r8) then
             !
             ! If variable-name is present or itypveg 
             ! AND this is crop or non-vegetated land-unit
             !
             if ( (present_var .or. itypveg_var) .and. &
                  ((vtypeo(no) >= croptype) .or. (typeo(no) >= nonveg)) )then
                if ( present_var ) ibufslo(no) = 0
                if ( itypveg_var ) ibufslo(no) = vtypeo(no)
             !
             ! Otherwise calculate it from the nearest neighbor
             !
             else
                call findMinDistPFTs( ncidi, ncido, no, nmin=n )
                ibufslo(no) = ibufsli(n)
             end if          !variable type
          end if             !output data with positive weight
       end do                !output data land loop
!$OMP END PARALLEL DO
    else
       write(*,*) 'subroutine interp_sl_int: no data written to typeo,vtypeo,no=', &
                   typeo(no),vtypeo(no),no
       stop
    end if

    call wrap_inq_varid (ncido, varname, varid)
    call wrap_put_var_int (ncido, varid, ibufslo)

    write(*,*) 'written variable ', trim(varname),' to output initial file'

    deallocate (ibufsli)
    deallocate (ibufslo)
    deallocate (wto)
    if ( allocated(vtypei) ) deallocate (vtypei)
    if ( allocated(vtypeo) ) deallocate (vtypeo)
    if ( allocated(typeo) )  deallocate (typeo)
    if ( allocated(wti) )    deallocate (wti)

  end subroutine interp_sl_int

end module interpinic
