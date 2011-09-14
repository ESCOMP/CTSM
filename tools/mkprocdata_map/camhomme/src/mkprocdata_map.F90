program mkprocdata_map

  use netcdf
  use shr_kind_mod, only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use fileutils
  use gridmapMod
  implicit none

  integer :: nDimensions    ! number of dimensions defined for this netCDF dataset.
  integer :: nVariables     ! number of variables defined for this netCDF dataset.
  integer :: nAttributes    ! number of global attributes defined for this netCDF dataset.
  integer :: unlimitedDimID ! ID of the unlimited dimension, if there is one
                            ! If no unlimited length dimension has been defined, 
                            ! -1 is returned.
  integer :: ier            ! returned error code
  integer :: nlon,nlat
  integer :: dimid
  integer :: dimid_grid
  integer :: dimid_lon, dimid_lat
  integer :: ncidi, ncido
  integer :: varidi, varido
  integer :: dimlen
  integer :: xtype
  integer :: ndimsi,ndimso
  integer :: dimidsi(4)     ! dimension id array
  integer :: dimidso(4)     ! dimension id array
  integer :: nAtts          ! number of variable attributes
  integer :: n,nv,nt,nd,na  ! indices
  integer :: attlen
  character(len=128):: dim1name
  character(len=128):: dimname
  character(len=128):: varname
  character(len=128):: attname 
  character(len=256):: locfn
  character(len=256):: filei
  character(len=256):: fileo
  character(len=256):: fmap
  character(len=128):: cattvalue
  real(r8)          :: dattvalue
  real(r4)          :: fattvalue
  integer           :: iattvalue
  type(gridmap_type):: tgridmap
  logical :: mapvar 
  real(r8), allocatable :: array2d(:,:)
  real(r8), allocatable :: lon(:),lat(:)

  !--------------------------------------------------------
  ! Read in mapping file - will have frac in it - and input and output domain
  !--------------------------------------------------------

  namelist /mkprocdata_map_in/ &
       filei,fileo,fmap

  read(5, mkprocdata_map_in, iostat=ier)
  if (ier /= 0) then
      write(6,*)'error: namelist input resulted in error code ',ier
      call abort()
  endif

  dim1name = 'ncol'

  call getfil (fmap, locfn, 0)

  call handle_ncerr(nf90_open(locfn, NF90_NOWRITE, ncidi))
  call handle_ncerr(nf90_inq_dimid(ncidi, 'ni_b', dimid))
  call handle_ncerr(nf90_inquire_dimension(ncidi, dimid, len=nlon))
  call handle_ncerr(nf90_inq_dimid(ncidi, 'nj_b', dimid))
  call handle_ncerr(nf90_inquire_dimension(ncidi, dimid, len=nlat))
  call handle_ncerr(nf90_close(ncidi))

  call gridmap_mapread(tgridmap, locfn)
  allocate(lon(nlon),lat(nlat))

  allocate(array2d(nlon,nlat))
  array2d(:,:) = reshape(tgridmap%xc_dst,(/nlon,nlat/))
  lon(:) = array2d(:,1)
  array2d(:,:) = reshape(tgridmap%yc_dst,(/nlon,nlat/))
  lat(:) = array2d(1,:)
  deallocate(array2d)

  call getfil (filei, locfn, 0)

  !--------------------------------------------------------
  ! Create output file (put it in define mode)
  !--------------------------------------------------------

  call handle_ncerr(nf90_create(fileo, NF90_64BIT_OFFSET, ncido))

  !--------------------------------------------------------
  ! Define output dimensions - creating file puts it in define mode 
  !--------------------------------------------------------

  call handle_ncerr(nf90_open(locfn, NF90_NOWRITE, ncidi))
  call handle_ncerr(nf90_inquire(ncidi, nDimensions, nVariables, &
       nAttributes, unlimitedDimId))

  write(6,*)'ndims= ',nDimensions,' nVariables= ',nVariables 
  do nd = 1,nDimensions
     ! Determine input dimensions
     call handle_ncerr(nf90_inquire_dimension(ncidi, dimid=nd, name=dimname, len=dimlen))
 
     ! Define output variables
     if (dimname == 'time') then
        call handle_ncerr(nf90_def_dim(ncido, name=dimname, len=nf90_unlimited, dimid=dimid))
     else
        if (trim(dimname) == dim1name) then
           dimid_grid= nd
           call handle_ncerr(nf90_def_dim(ncido, name='lon', len=nlon, dimid=dimid_lon)) 
           call handle_ncerr(nf90_def_dim(ncido, name='lat', len=nlat, dimid=dimid_lat)) 
        else
           call handle_ncerr(nf90_def_dim(ncido, name=dimname, len=dimlen, dimid=dimid))
        end if
     end if
     write(6,*)'n = ',nd,' dimname= ',trim(dimname)
  end do

  !--------------------------------------------------------
  ! Define output variables
  !--------------------------------------------------------

  ! Loop over input variables
  do nv = 1,nVariables

     ! Determine input variable
     call handle_ncerr(nf90_Inquire_Variable(ncid=ncidi, varid=nv, natts=natts, &
          name=varname, ndims=ndimsi, dimids=dimidsi, xtype=xtype)) 
     
     ! Determine output dimension ids
     if (dimidsi(1) == dimid_grid) then
        ndimso = ndimsi + 1
        dimidso(1) = dimid_lon
        dimidso(2) = dimid_lat
        do n = 2,ndimsi
           call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(n), dimname))
           call handle_ncerr(nf90_inq_dimid(ncido, dimname, dimidso(n+1)))
        end do
     else
        ndimso = ndimsi
        do n = 1,ndimsi
           call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(n), dimname))
           call handle_ncerr(nf90_inq_dimid(ncido, dimname, dimidso(n))) 
        end do
     end if

     ! Define output variable and attributes 
     if (trim(varname) == 'lon') then
        call handle_ncerr(nf90_def_var(ncido, name=varname, xtype=xtype, &
             dimids=dimid_lon, varid=varido))
     else if (trim(varname) == 'lat') then
        call handle_ncerr(nf90_def_var(ncido, name=varname, xtype=xtype, &
             dimids=dimid_lat, varid=varido))
     else
        call handle_ncerr(nf90_def_var(ncido, name=varname, xtype=xtype, &
             dimids=dimidso(1:ndimso), varid=varido))
     end if

     do na = 1,natts
        call handle_ncerr(nf90_inq_attname(ncidi, varid=nv, attnum=na, name=attname))
        call handle_ncerr(nf90_inquire_attribute(ncidi, varid=nv, name=attname, &
             len=attlen, xtype=xtype))
        if (xtype == nf90_char) then
           call handle_ncerr(nf90_get_att(ncidi, varid=nv, name=attname, values=cattvalue))
           call handle_ncerr(nf90_put_att(ncido, varid=varido, name=attname, &
                values=cattvalue(1:attlen)))
        else if (xtype == nf90_double) then
           call handle_ncerr(nf90_get_att(ncidi, varid=nv, name=attname, values=dattvalue))
           call handle_ncerr(nf90_put_att(ncido, varid=varido, name=attname, values=dattvalue))
        else if (xtype == nf90_float) then
           call handle_ncerr(nf90_get_att(ncidi, varid=nv, name=attname, values=fattvalue))
           call handle_ncerr(nf90_put_att(ncido, varid=varido, name=attname, values=fattvalue))
        else if (xtype == nf90_int) then
           call handle_ncerr(nf90_get_att(ncidi, varid=nv, name=attname, values=iattvalue))
           call handle_ncerr(nf90_put_att(ncido, varid=varido, name=attname, values=iattvalue))
        end if
     end do
  end do

  ! End define mode
  call handle_ncerr(nf90_enddef(ncido))  

  !--------------------------------------------------------
  ! Read in variables and write them out
  !--------------------------------------------------------

  do nv = 1,nVariables
     call outputvar(nv, nlon, nlat, ncidi, ncido, tgridmap)
  end do

  call handle_ncerr(nf90_close(ncidi))
  call handle_ncerr(nf90_close(ncido))

contains

  subroutine outputvar(nv, nlon, nlat, ncidi, ncido, tgridmap)

    integer, intent(in) :: nv
    integer, intent(in) :: nlon
    integer, intent(in) :: nlat
    integer, intent(in) :: ncidi
    integer, intent(in) :: ncido
    type(gridmap_type), intent(inout) :: tgridmap

    integer :: len1,len2,len3,len4
    integer :: n1,n2,n3,n4
    integer :: dimidsi(6)     
    integer :: dimidso(6)     
    real(r8), allocatable:: rarrayi(:),rarrayo(:)
    character(len=128):: dimname
    character(len=128):: varname
    character(len=128):: attname 
    logical :: first_time = .true.
    integer :: vid
    real(r8):: spval = 1.e36 

    varidi = nv
    call handle_ncerr(nf90_Inquire_Variable(ncid=ncidi, varid=varidi, &
         name=varname, ndims=ndimsi, dimids=dimidsi, xtype=xtype))

    write(6,*)'n = ',nv,' varname= ',trim(varname)
    if (ndimsi == 0) then
       return
    end if
    
    ! special variables not currently handled
    if (trim(varname) == 'date_written' .or. &
        trim(varname) == 'time_written' ) then
       write(6,*)'skipping writing out variable ',trim(varname)         
       return
    end if
    
    call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(1), &
         len=len1, name=dimname))
    allocate(rarrayi(len1))

    if (trim(dimname)==dim1name) then
       mapvar = .true.
    else
       mapvar = .false.
    end if
    
    varido = nv
    call handle_ncerr(nf90_Inquire_Variable(ncid=ncido, varid=varido, &
         name=varname, ndims=ndimso, dimids=dimidso, xtype=xtype))

    if (trim(varname) == 'lon') then
       write(6,*)' outputing lon'
       call handle_ncerr(nf90_put_var(ncido, varido, lon))
       return
    end if
    if (trim(varname) == 'lat') then
       write(6,*)' outputing lat'
       call handle_ncerr(nf90_put_var(ncido, varido, lat))
       return
    end if

    if (mapvar) then
       call handle_ncerr(nf90_inquire_dimension(ncido, dimidso(1), &
            len=len1, name=dimname))
       call handle_ncerr(nf90_inquire_dimension(ncido, dimidso(2), &
            len=len2, name=dimname))
       allocate(rarrayo(len1*len2))
    else
       len1 = size(rarrayi)
       allocate(rarrayo(len1))
    end if
    
    if (ndimsi == 1) then

       call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(1), len=len1))
       call handle_ncerr(nf90_get_var(ncidi, varidi, rarrayi))
       if (mapvar) then
          call gridmap_map(tgridmap, rarrayi, rarrayo)
          call handle_ncerr(nf90_put_var(ncido, varido, &
               reshape(rarrayo,(/nlon,nlat/))))
       else
          rarrayo(:)= rarrayi(:)
          call handle_ncerr(nf90_put_var(ncido, varido, rarrayo))
       end if
    
    else if (ndimsi == 2) then

       call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(1), len=len1))
       call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(2), len=len2))
       do n2 = 1,len2
          call handle_ncerr(nf90_get_var(ncidi, varidi, rarrayi, &
               start=(/1,n2/), count=(/len1,1/)))
          if (mapvar) then
             call gridmap_map(tgridmap, rarrayi, rarrayo)
             call handle_ncerr(nf90_put_var(ncido, varido, &
                  reshape(rarrayo,(/nlon,nlat/)), &
                  start=(/1,1,n2/), count=(/nlon,nlat/)))
          else
             rarrayo(:)= rarrayi(:)
             call handle_ncerr(nf90_put_var(ncido, varido, rarrayo, &
                  start=(/1,n2/), count=(/len1,1/)))
           end if
        end do

     else if (ndimsi == 3) then

        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(1), len=len1))
        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(2), len=len2))
        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(3), len=len3))
        do n2 = 1,len2
        do n3 = 1,len3           
           call handle_ncerr(nf90_get_var(ncidi, varidi, rarrayi, &
                start=(/1,n2,n3/), count=(/len1,1,1/)))
           if (mapvar) then
              call gridmap_map(tgridmap, rarrayi, rarrayo)
              call handle_ncerr(nf90_put_var(ncido, varido, &
                   reshape(rarrayo,(/nlon,nlat/)), &
                   start=(/1,1,n2,n3/), count=(/nlon,nlat/)))
           else
              rarrayo(:)= rarrayi(:)
              call handle_ncerr(nf90_put_var(ncido, varido, rarrayo, &
                   start=(/1,n2,n3/), count=(/len1,1,1/)))
           end if
        end do
        end do

     else if (ndimsi == 4)then

        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(1), len=len1))
        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(2), len=len2))
        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(3), len=len3))
        call handle_ncerr(nf90_inquire_dimension(ncidi, dimidsi(4), len=len4))
        do n2 = 1,len2
        do n3 = 1,len3           
        do n4 = 1,len4           
           call handle_ncerr(nf90_get_var(ncidi, varidi, rarrayi, &
                start=(/1,n2,n3,n4/), count=(/len1,1,1,1/)))
           if (mapvar) then
              call gridmap_map(tgridmap, rarrayi, rarrayo)
              call handle_ncerr(nf90_put_var(ncido, varido, &
                   reshape(rarrayo,(/nlon,nlat/)), &
                   start=(/ 1,1,n2,n3,n4 /), count=(/nlon,nlat,1,1,1/)))
           else
              rarrayo(:)= rarrayi(:)
              call handle_ncerr(nf90_put_var(ncido, varido, rarrayo, &
                   start=(/1,n2,n3,n4/), count=(/len1,1,1,1/)))
           end if
        end do
        end do
        end do

     else

        write(6,*)'ndimsi= ',ndimsi,' is not supported'
        stop

     end if

     deallocate(rarrayi)
     deallocate(rarrayo)

  end subroutine outputvar

  subroutine handle_err(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
       write(6,*)trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine handle_err

  subroutine handle_ncerr(ret)
     integer         , intent(in) :: ret  ! return code from netCDF library routine
     if ( ret .ne. NF90_NOERR ) then
        write(6,*) nf90_strerror( ret )
        stop
     endif
   end subroutine handle_ncerr

end program mkprocdata_map
