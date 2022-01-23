module mksoilcolMod

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort
  use mkpioMod     , only : mkpio_get_rawdata, mkpio_get_dim_lengths
  use mkpioMod     , only : mkpio_iodesc_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod    , only : regrid_rawdata
  use mkutilsMod   , only : chkerr
  use mkvarctl     , only : root_task, ndiag
  use mkvarpar     , only : nlevsoi

  implicit none
  private

  public mksoilcol      ! Set soil colors

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mkorganic(file_mesh_i, file_data_i, field_o, mesh_o, soil_color_o, nsoilcol)

    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i     ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i     ! input data file name
    type(ESMF_Field)  , intent(inout) :: field_o
    type(ESMF_Mesh)   , intent(in)    :: mesh_o
    integer           , intent(out)   :: soil_color_o(:) ! soil color classes
    integer           , intent(out)   :: nsoicol         ! number of soil colors 
    integer           , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle)         :: routehandle
    type(ESMF_Mesh)                :: mesh_i
    type(ESMF_Field)               :: field_i
    type(file_desc_t)              :: pioid
    type(var_desc_t)               :: pio_varid
    type(io_desc_t)                :: pio_iodesc
    integer                        :: pio_vartype
    real(r8), allocatable          :: organic_i(:,:)
    integer                        :: ni,no
    integer                        :: ns_i, ns_o
    integer                        :: k,l,m, ! indices
    integer                        :: rcode, ier             ! error status
    integer                        :: srcMaskValue = 0
    integer                        :: dstMaskValue = -987987 ! spval for RH mask values
    integer                        :: srcTermProcessing_Value = 0
    integer                        :: ndims
    integer, allocatable           :: dim_lengths(:)
    real(r4), allocatable          :: data_real(:,:)
    real(r8), allocatable          :: data_double(:,:)
    integer , allocatable          :: soil_color_i(:)       ! input grid: BATS soil color
    character(len=35), allocatable :: col(:)                ! name of each color
    integer                        :: ier                   ! error status
    character(len=*), parameter :: subname = 'mksoilcol'
    !-----------------------------------------------------------------------

    write (6,*) 'Attempting to make soil color classes .....'

    ! -----------------------------------------------------------------
    ! Read input file
    ! -----------------------------------------------------------------

    ns_o = ldomain%ns

    ! Obtain input grid info, read local fields

    call domain_read(tdomain,datfname)
    ns_i = tdomain%ns
    allocate(soil_color_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(frac_dst(ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort()

    write (6,*) 'Open soil color file: ', trim(datfname)
    call check_ret(nf_open(datfname, 0, ncid), subname)
    call check_ret(nf_inq_varid (ncid, 'SOIL_COLOR', varid), subname)
    call check_ret(nf_get_var_int (ncid, varid, soil_color_i), subname)
    call check_ret(nf_close(ncid), subname)

    nsoicol = maxval(soil_color_i)
    write(6,*)'nsoicol = ',nsoicol

    ! -----------------------------------------------------------------
    ! Define the model color classes: 0 to nsoicol
    ! -----------------------------------------------------------------

    if (nsoicol == 20) then
       col(0)  = 'no soil                            '
       col(1)  = 'class 1: light                     '
       col(2)  = 'class 2:                           '
       col(3)  = 'class 3:                           '
       col(4)  = 'class 4:                           '
       col(5)  = 'class 5:                           '
       col(6)  = 'class 6:                           '
       col(7)  = 'class 7:                           '
       col(8)  = 'class 8:                           '
       col(9)  = 'class 9:                           '
       col(10) = 'class 10:                          '
       col(11) = 'class 11:                          '
       col(12) = 'class 12:                          '
       col(13) = 'class 13:                          '
       col(14) = 'class 14:                          '
       col(15) = 'class 15:                          '
       col(16) = 'class 16:                          '
       col(17) = 'class 17:                          '
       col(18) = 'class 18:                          '
       col(19) = 'class 19:                          '
       col(20) = 'class 20: dark                     '
    else if (nsoicol == 8) then
       col(0) = 'no soil                            '
       col(1) = 'class 1: light                     '
       col(2) = 'class 2:                           '
       col(3) = 'class 3:                           '
       col(4) = 'class 4:                           '
       col(5) = 'class 5:                           '
       col(6) = 'class 6:                           '
       col(7) = 'class 7:                           '
       col(8) = 'class 8: dark                      '
    else
       write(6,*)'nsoicol value of ',nsoicol,' is not currently supported'
       call shr_sys_abort()
    end if

    ! Error check soil_color if it is set
    if ( soil_color /= unsetcol )then

       if ( soil_color > nsoicol )then
          write(6,*)'soil_color is out of range = ', soil_color
          call shr_sys_abort()
       end if

       do no = 1,ns_o
          soil_color_o(no) = soil_color
       end do

    else

       call gridmap_mapread(tgridmap, mapfname)

       ! Error checks for domain and map consistencies

       call domain_checksame( tdomain, ldomain, tgridmap )

       ! Obtain frac_dst
       call gridmap_calc_frac_dst(tgridmap, tdomain%mask, frac_dst)

       ! Determine dominant soil color for each output cell

       call dominant_soil_color( &
            tgridmap = tgridmap, &
            mask_i = tdomain%mask, &
            soil_color_i = soil_color_i, &
            nsoicol = nsoicol, &
            soil_color_o = soil_color_o)

       ! Global sum of output field 

       allocate(mask_r8(ns_i), stat=ier)
       if (ier/=0) call shr_sys_abort()
       mask_r8 = tdomain%mask
       call gridmap_check( tgridmap, mask_r8, frac_dst, subname )

    end if

    ! Deallocate dynamic memory

    call domain_clean(tdomain)
    if ( soil_color == unsetcol )then
       call gridmap_clean(tgridmap)
    end if
    deallocate (soil_color_i,gast_i,gast_o,col, frac_dst, mask_r8)

    write (6,*) 'Successfully made soil color classes'
    write (6,*)

  end subroutine mksoilcol

  !-----------------------------------------------------------------------
  subroutine dominant_soil_color(tgridmap, mask_i, soil_color_i, nsoicol, soil_color_o)
    !
    ! Determine the dominant soil color in each output cell
    !
    ! input/output variables
    type(gridmap_type) , intent(in)  :: tgridmap
    integer            , intent(in)  :: mask_i(:)       ! input grid: land mask (1 = land, 0 = ocean)
    integer            , intent(in)  :: soil_color_i(:) ! input grid: BATS soil color
    integer            , intent(in)  :: nsoicol         ! number of soil colors
    integer            , intent(out) :: soil_color_o(:) ! output grid: soil color classes
    !
    ! local variables:
    integer, parameter :: num = 2             ! set soil mapunit number
    integer  :: wsti(num)                     ! index to 1st and 2nd largest wst
    integer  :: k, n, ni, no, ns_i, ns_o
    real(r8) :: wt                            ! map overlap weight
    real(r8), allocatable :: wst(:,:)         ! overlap weights, by surface type
    logical :: has_color                      ! whether this grid cell has non-zero color
    integer, parameter :: miss = 99999        ! missing data indicator

    character(len=*), parameter :: subname = 'dominant_soil_color'
    !-----------------------------------------------------------------------

    ns_i = size(mask_i)
    if (size(soil_color_i) /= ns_i) then
       write(6,*) subname, ' ERROR: size of soil_color_i should match size of mask_i'
       write(6,*) 'size(mask_i), size(soil_color_i) = ', &
            size(mask_i), size(soil_color_i)
       call shr_sys_abort()
    end if

    ! find area of overlap for each soil color for each no

    ns_o = size(soil_color_o)
    allocate(wst(0:nsoicol,ns_o))
    wst(0:nsoicol,:) = 0

    ! TODO: need to do a loop to determine
    ! the maximum number of over lap cells throughout the grid 
    ! first get an array that is novr(ns_o) and fill this in - then set
    ! maxovr - to max(novr) - then allocate the array wst to be size of
    ! maxovr,ns_o or 0:nsoilcol,ns_o

    do n = 1,tgridmap%ns
       ni = tgridmap%src_indx(n)
       no = tgridmap%dst_indx(n)
       wt = tgridmap%wovr(n) * mask_i(ni)
       k  = soil_color_i(ni) * mask_i(ni)
       wst(k,no) = wst(k,no) + wt
    enddo

    soil_color_o(:) = 0
    do no = 1,ns_o

       ! If the output cell has any non-zero-colored inputs, then set the weight of
       ! zero-colored inputs to 0, to ensure that the zero-color is NOT dominant.
       if (any(wst(1:nsoicol,no) > 0.)) then
          has_color = .true.
          wst(0,no) = 0.0
       else
          has_color = .false.
       end if

       ! Rank non-zero weights by color type. wsti(1) is the most extensive
       ! color type. 

       if (has_color) then
          call mkrank (nsoicol, wst(0:nsoicol,no), miss, wsti, num)
          soil_color_o(no) = wsti(1)
       end if

       ! If land but no color, set color to 15 (in older dataset generic 
       ! soil color 4)

       if (nsoicol == 8) then
          if (soil_color_o(no)==0) then
             soil_color_o(no) = 4
          end if
       else if (nsoicol == 20) then
          if (soil_color_o(no)==0) then
             soil_color_o(no) = 15
          end if
       else
          write(6,*) 'MKSOILCOL error: unhandled nsoicol: ', nsoicol
          call shr_sys_abort()
       end if

       ! Error checks

       if (soil_color_o(no) < 0 .or. soil_color_o(no) > nsoicol) then
          write (6,*) 'MKSOILCOL error: land model soil color = ', &
               soil_color_o(no),' is not valid for lon,lat = ',no
          call shr_sys_abort()
       end if

    end do

    deallocate (wst)

  end subroutine dominant_soil_color

end module mksoilcolMod
