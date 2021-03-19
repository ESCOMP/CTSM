module lnd_set_decomp_and_domain

  use shr_kind_mod , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use spmdMod      , only : masterproc
  use clm_varctl   , only : iulog
  use perf_mod     , only : t_startf, t_stopf, t_barrierf

  implicit none
  private ! except

  ! public member routines
  public :: lnd_set_decomp_and_domain_from_surfrd

  ! private member routines
  private :: surfrd_get_globmask  ! Reads global land mask (needed for setting domain decomp)
  private :: surfrd_get_grid      ! Read grid/ladnfrac data into domain (after domain decomp)

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine lnd_set_decomp_and_domain_from_surfrd(noland, ni, nj)

    ! Initialize ldecomp and ldomain data types

    use clm_varpar    , only: nlevsoi
    use clm_varctl    , only: fatmlndfrc, use_soil_moisture_streams
    use decompInitMod , only: decompInit_lnd, decompInit_lnd3D
    use decompMod     , only: bounds_type, get_proc_bounds
    use domainMod     , only: ldomain, domain_init, domain_check

    ! input/output variables
    logical, intent(out) :: noland
    integer, intent(out) :: ni, nj ! global grid sizes

    ! local variables
    integer ,pointer  :: amask(:)   ! global land mask
    integer           :: begg, endg ! processor bounds
    type(bounds_type) :: bounds     ! bounds
    character(len=32) :: subname = 'lnd_set_decomp_and_domain_from_surfrd'
    !-----------------------------------------------------------------------

    ! Read in global land grid and land mask (amask)- needed to set decomposition
    ! global memory for amask is allocate in surfrd_get_glomask - must be deallocated below
    if (masterproc) then
       write(iulog,*) 'Attempting to read global land mask from ',trim(fatmlndfrc)
    endif

    ! Get global mask, ni and nj
    call surfrd_get_globmask(filename=fatmlndfrc, mask=amask, ni=ni, nj=nj)

    ! Exit early if no valid land points
    if ( all(amask == 0) )then
       if (masterproc) write(iulog,*) trim(subname)//': no valid land points do NOT run clm'
       noland = .true.
       return
    else
       noland = .false.
    end if

    ! Initialize ldecomp data type
    ! Determine ctsm gridcell decomposition and processor bounds for gridcells
    call decompInit_lnd(ni, nj, amask)
    deallocate(amask)
    if (use_soil_moisture_streams) call decompInit_lnd3D(ni, nj, nlevsoi)

    ! Initialize bounds for just gridcells
    ! Remaining bounds (landunits, columns, patches) will be determined
    ! after the call to decompInit_glcp - so get_proc_bounds is called
    ! twice and the gridcell information is just filled in twice
    call get_proc_bounds(bounds)

    ! Get grid cell bounds values
    begg = bounds%begg
    endg = bounds%endg

    ! Initialize ldomain data type
    if (masterproc) then
       write(iulog,*) 'Attempting to read ldomain from ',trim(fatmlndfrc)
    endif
    call surfrd_get_grid(begg, endg, ldomain, fatmlndfrc)
    if (masterproc) then
       call domain_check(ldomain)
    endif
    ldomain%mask = 1  !!! TODO - is this needed?

  end subroutine lnd_set_decomp_and_domain_from_surfrd

  !-----------------------------------------------------------------------
  subroutine surfrd_get_globmask(filename, mask, ni, nj)

    ! Read the surface dataset grid related information
    ! This is used to set the domain decomposition - so global data is read here

    use fileutils  , only : getfil
    use ncdio_pio  , only : ncd_io, ncd_pio_openfile, ncd_pio_closefile, ncd_inqfdims, file_desc_t
    use abortutils , only : endrun
    use shr_log_mod, only : errMsg => shr_log_errMsg

    ! input/output variables
    character(len=*), intent(in)    :: filename  ! grid filename
    integer         , pointer       :: mask(:)   ! grid mask
    integer         , intent(out)   :: ni, nj    ! global grid sizes

    ! local variables
    logical               :: isgrid2d
    integer               :: dimid,varid ! netCDF id's
    integer               :: ns          ! size of grid on file
    integer               :: n,i,j       ! index
    integer               :: ier         ! error status
    type(file_desc_t)     :: ncid        ! netcdf id
    character(len=256)    :: varname     ! variable name
    character(len=256)    :: locfn       ! local file name
    logical               :: readvar     ! read variable in or not
    integer , allocatable :: idata2d(:,:)
    character(len=32) :: subname = 'surfrd_get_globmask' ! subroutine name
    !-----------------------------------------------------------------------

    if (filename == ' ') then
       mask(:) = 1
    else
       ! Check if file exists
       if (masterproc) then
          if (filename == ' ') then
             write(iulog,*) trim(subname),' ERROR: filename must be specified '
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif
       end if

       ! Open file
       call getfil( filename, locfn, 0 )
       call ncd_pio_openfile (ncid, trim(locfn), 0)

       ! Determine dimensions and if grid file is 2d or 1d
       call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)
       if (masterproc) then
          write(iulog,*)'lat/lon grid flag (isgrid2d) is ',isgrid2d
       end if
       allocate(mask(ns))
       mask(:) = 1
       if (isgrid2d) then
          ! Grid is 2d
          allocate(idata2d(ni,nj))
          idata2d(:,:) = 1
          call ncd_io(ncid=ncid, varname='LANDMASK', data=idata2d, flag='read', readvar=readvar)
          if (.not. readvar) then
             call ncd_io(ncid=ncid, varname='mask', data=idata2d, flag='read', readvar=readvar)
          end if
          if (readvar) then
             do j = 1,nj
                do i = 1,ni
                   n = (j-1)*ni + i
                   mask(n) = idata2d(i,j)
                enddo
             enddo
          end if
          deallocate(idata2d)
       else
          ! Grid is not 2d
          call ncd_io(ncid=ncid, varname='LANDMASK', data=mask, flag='read', readvar=readvar)
          if (.not. readvar) then
             call ncd_io(ncid=ncid, varname='mask', data=mask, flag='read', readvar=readvar)
          end if
       end if
       if (.not. readvar) call endrun( msg=' ERROR: landmask not on fatmlndfrc file'//errMsg(sourcefile, __LINE__))

       ! Close file
       call ncd_pio_closefile(ncid)
    end if

  end subroutine surfrd_get_globmask

  !-----------------------------------------------------------------------
  subroutine surfrd_get_grid(begg, endg, ldomain, filename, glcfilename)

    ! Read the surface dataset grid related information:
    ! This is called after the domain decomposition has been created
    ! - real latitude  of grid cell (degrees)
    ! - real longitude of grid cell (degrees)

    use clm_varcon   , only : spval, re, grlnd
    use domainMod    , only : domain_type, domain_init, domain_clean, lon1d, lat1d
    use fileutils    , only : getfil
    use abortutils   , only : endrun
    use shr_log_mod  , only : errMsg => shr_log_errMsg
    use ncdio_pio    , only : file_desc_t, var_desc_t, ncd_pio_openfile, ncd_pio_closefile
    use ncdio_pio    , only : ncd_io, check_var, ncd_inqfdims, check_dim_size, ncd_inqdid, ncd_inqdlen
    use pio

    ! input/output variables
    integer                    , intent(in)    :: begg, endg
    type(domain_type)          , intent(inout) :: ldomain     ! domain to init
    character(len=*)           , intent(in)    :: filename    ! grid filename
    character(len=*) ,optional , intent(in)    :: glcfilename ! glc mask filename

    ! local variables
    type(file_desc_t)     :: ncid               ! netcdf id
    integer               :: beg                ! local beg index
    integer               :: end                ! local end index
    integer               :: ni,nj,ns           ! size of grid on file
    integer               :: dimid,varid        ! netCDF id's
    integer               :: start(1), count(1) ! 1d lat/lon array sections
    integer               :: ier,ret            ! error status
    logical               :: readvar            ! true => variable is on input file
    logical               :: isgrid2d           ! true => file is 2d lat/lon
    logical               :: istype_domain      ! true => input file is of type domain
    real(r8), allocatable :: rdata2d(:,:)       ! temporary
    character(len=16)     :: vname              ! temporary
    character(len=256)    :: locfn              ! local file name
    integer               :: n                  ! indices
    character(len=32) :: subname = 'surfrd_get_grid'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       if (filename == ' ') then
          write(iulog,*) trim(subname),' ERROR: filename must be specified '
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    end if

    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Determine dimensions
    call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)

    ! Determine isgrid2d flag for domain
    call domain_init(ldomain, isgrid2d=isgrid2d, ni=ni, nj=nj, nbeg=begg, nend=endg)

    ! Determine type of file - old style grid file or new style domain file
    call check_var(ncid=ncid, varname='xc', readvar=readvar)
    if (readvar)then
        istype_domain = .true.
    else
        istype_domain = .false.
    end if

    ! Read in area, lon, lat
    if (istype_domain) then
       call ncd_io(ncid=ncid, varname= 'area', flag='read', data=ldomain%area, &
            dim1name=grlnd, readvar=readvar)
       ! convert from radians**2 to km**2
       ldomain%area = ldomain%area * (re**2)
       if (.not. readvar) call endrun( msg=' ERROR: area NOT on file'//errMsg(sourcefile, __LINE__))
       call ncd_io(ncid=ncid, varname= 'xc', flag='read', data=ldomain%lonc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: xc NOT on file'//errMsg(sourcefile, __LINE__))
       call ncd_io(ncid=ncid, varname= 'yc', flag='read', data=ldomain%latc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: yc NOT on file'//errMsg(sourcefile, __LINE__))
    else
       call endrun( msg=" ERROR: can no longer read non domain files" )
    end if

    if (isgrid2d) then
       allocate(rdata2d(ni,nj), lon1d(ni), lat1d(nj))
       if (istype_domain) vname = 'xc'
       call ncd_io(ncid=ncid, varname=trim(vname), data=rdata2d, flag='read', readvar=readvar)
       lon1d(:) = rdata2d(:,1)
       if (istype_domain) vname = 'yc'
       call ncd_io(ncid=ncid, varname=trim(vname), data=rdata2d, flag='read', readvar=readvar)
       lat1d(:) = rdata2d(1,:)
       deallocate(rdata2d)
    end if

    ! Check lat limited to -90,90
    if (minval(ldomain%latc) < -90.0_r8 .or. &
        maxval(ldomain%latc) >  90.0_r8) then
       write(iulog,*) trim(subname),' WARNING: lat/lon min/max is ', &
            minval(ldomain%latc),maxval(ldomain%latc)
    endif
    if ( any(ldomain%lonc < 0.0_r8) )then
       call endrun( msg=' ERROR: lonc is negative (see https://github.com/ESCOMP/ctsm/issues/507)' &
            //errMsg(sourcefile, __LINE__))
    endif
    call ncd_io(ncid=ncid, varname='mask', flag='read', data=ldomain%mask, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun( msg=' ERROR: LANDMASK NOT on fracdata file'//errMsg(sourcefile, __LINE__))
    end if
    call ncd_io(ncid=ncid, varname='frac', flag='read', data=ldomain%frac, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun( msg=' ERROR: LANDFRAC NOT on fracdata file'//errMsg(sourcefile, __LINE__))
    end if

    call ncd_pio_closefile(ncid)

  end subroutine surfrd_get_grid

end module lnd_set_decomp_and_domain
