module mksoiltexMod

  !-----------------------------------------------------------------------
  ! Make soil data (texture)
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata, mkpio_get_dimlengths
  use mkpioMod       , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod      , only : regrid_rawdata, create_routehandle_r8, get_meshareas
  use mkutilsMod     , only : chkerr
  use mkvarctl
  use mkvarpar

  implicit none
  private           ! By default make data private

  public :: mksoiltex      ! Set soil texture

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mksoiltex(file_mesh_i, file_data_i, mesh_o, sand_o, clay_o, mapunit_o, rc)
    !
    ! make %sand and %clay from IGBP soil data, which includes
    ! igbp soil 'mapunits' and their corresponding textures
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o      ! output mesh
    real(r8)          , intent(inout) :: sand_o(:,:) ! % sand (output grid)
    real(r8)          , intent(inout) :: clay_o(:,:) ! % clay (output grid)
    integer           , intent(inout) :: mapunit_o(:)
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid
    type(var_desc_t)       :: pio_varid
    integer                :: pio_vartype
    integer                :: dimid
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer                :: n,l,m
    integer                :: gindex, lindex
    character(len=38)      :: typ                     ! soil texture based on ...
    integer                :: nlay                    ! number of soil layers
    integer                :: mapunittemp             ! temporary igbp soil mapunit
    integer                :: maxovr
    real(r8), allocatable  :: frac_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r8), allocatable  :: sand_i(:,:)             ! input grid: percent sand
    real(r8), allocatable  :: clay_i(:,:)             ! input grid: percent clay
    integer, parameter     :: num=2                   ! set soil mapunit number
    integer, parameter     :: nlsm=4                  ! number of soil textures
    character(len=38)      :: soil(0:nlsm)            ! name of each soil texture
    real(r8)               :: gast_i(0:nlsm)          ! global area, by texture type
    real(r8)               :: gast_o(0:nlsm)          ! global area, by texture type
    real(r8)               :: wt                      ! map overlap weight
    real(r8)               :: sum_fldi                ! global sum of dummy input fld
    real(r8)               :: sum_fldo                ! global sum of dummy output fld
    integer                :: rcode, ier              ! error status
    real(r8)               :: sumtex
    integer                :: mapunit_value_max
    integer                :: mapunit_value_min
    integer                :: mapunit_value
    integer                :: nmax
    integer                :: loop, nloops
    real(r4)               :: max_value 
    integer                :: max_index(1)
    real(r8), allocatable  :: mapunit_i(:)            ! input grid: igbp soil mapunits
    real(r8), allocatable  :: data_i(:,:)
    real(r8), allocatable  :: data_o(:,:)
    real(r8), allocatable  :: global_max_value(:)
    real(r8), allocatable  :: global_max_index(:)
    character(len=*), parameter :: subname = 'mksoiltex'
    !-----------------------------------------------------------------------

    if (root_task) then
       write (ndiag,'(a)') 'Attempting to make %sand and %clay .....'
    end if

    if ( soil_clay_override /= unsetsoil )then
       write(6,*) 'Replace soil clay % for all points with: ', soil_clay_override
       if ( soil_sand_override == unsetsoil )then
          write (6,*) subname//':error: soil_clay set, but NOT soil_sand'
          call shr_sys_abort()
       end if
    end if
    if ( soil_sand_override /= unsetsoil )then
       write(6,*) 'Replace soil sand % for all points with: ', soil_sand_override
       if ( soil_clay_override == unsetsoil )then
          write (6,*) subname//':error: soil_sand set, but NOT soil_clay'
          call shr_sys_abort()
       end if
       sumtex = soil_sand_override + soil_clay_override
       if ( sumtex < 0.0_r8 .or. sumtex > 100.0_r8 )then
          write (6,*) subname//':error: soil_sand and soil_clay out of bounds: sand, clay = ', &
               soil_sand_override, soil_clay_override
          call shr_sys_abort()
       end if
    end if

    if (soil_sand_override /= unsetsoil .and. soil_clay_override /= unsetsoil) then
       if (root_task) then
          write(ndiag,'(a,i8)') ' Overriding soil color for all points with: ', soil_color_override
       end if
       sand_o(:,:) = soil_sand_override
       clay_o(:,:) = soil_clay_override
       RETURN
    end if

    ! Define the model surface types: 0:4
    soil(0) = 'no soil: ocean, glacier, lake, no data'
    soil(1) = 'clays                                 '
    soil(2) = 'sands                                 '
    soil(3) = 'loams                                 '
    soil(4) = 'silts                                 '

    ! Open input data file
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Create a route handle between the input and output mesh
    call create_routehandle_r8(mesh_i, mesh_o, routehandle, rc)
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Determine ns_i and allocate data_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o and allocate data_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read in mapunit data
    allocate(mapunit_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'MAPUNITS', mesh_i, mapunit_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! TODO: Determine minimum and maximum mapunit values across all processors- for now hardwire

    ! Now determine data_i as a real 2d array - for every possible soil color create a global
    ! field with gridcells equal to 1 for that soil color and zero elsewhere
    rcode = pio_inq_dimid  (pioid, 'max_value_mapunit', dimid)
    rcode = pio_inq_dimlen (pioid, dimid, mapunit_value_max)

    mapunit_value_min = 0
    nmax = 100
    nloops = (mapunit_value_max - mapunit_value_min + nmax)/nmax
    !write(6,*)'nloops = ',nloops

    allocate(global_max_value(ns_o)) ; global_max_value(:) = -999.
    if (ier/=0) call shr_sys_abort()
    allocate(global_max_index(ns_o)) ; global_max_index(:) = 0
    if (ier/=0) call shr_sys_abort()
    allocate(data_i(nmax,ns_i))
    if (ier/=0) call shr_sys_abort()
    allocate(data_o(nmax,ns_o))
    if (ier/=0) call shr_sys_abort()

    mapunit_o(:) = 0
    do loop = 1,nloops
       data_i(:,:) = 0._r8
       do lindex = 0,nmax-1
          mapunit_value = lindex + (nmax*(loop-1))
          if (mapunit_value <= mapunit_value_max) then
             do ni = 1,ns_i
                if (int(mapunit_i(ni)) == mapunit_value) then
                   data_i(lindex+1,ni) = 1._r4
                end if
             end do
          end if
       end do
       
       call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 1, nmax, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       do no = 1,ns_o
          max_index = maxloc(data_o(:,no))
          max_value = data_o(max_index(1),no)
          if ( max_value > global_max_value(no)) then
             global_max_value(no) = max_value
             mapunit_o(no) = max_index(1) + (nmax*(loop-1)) - 1
             if (mapunit_o(no) < 0 .or. mapunit_o(no) > mapunit_value_max) then
                call shr_sys_abort('mapunit_o is invalid')
             end if
          end if
       end do
    end do

    deallocate(data_i)
    deallocate(data_o)

    ! Get dimensions from input file and allocate memory for sand_i and clay_i
    rcode = pio_inq_dimid  (pioid, 'number_of_layers', dimid)
    rcode = pio_inq_dimlen (pioid, dimid, nlay)
    allocate(sand_i(mapunit_value_max,nlay), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(clay_i(mapunit_value_max,nlay), stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! read in sand_i and clay_i (they will read in total on all processors)
    rcode = pio_inq_varid(pioid, 'PCT_SAND', pio_varid)
    rcode = pio_get_var(pioid, pio_varid, sand_i)
    rcode = pio_inq_varid(pioid, 'PCT_CLAY', pio_varid)
    rcode = pio_get_var(pioid, pio_varid, clay_i)

    ! Set soil texture as follows:
    ! a. Use dominant igbp soil mapunit based on area of overlap unless 'no data' is dominant
    ! b. If this has no data, use loam for soil texture

    do no = 1,ns_o
       if (mapunit_o(no) > 0) then
          ! valid value is obtained
          if (mapunit_o(no) > mapunit_value_max) then
             call shr_sys_abort("mapunit_o is out of bounds")
          end if
          do l = 1, nlay
             sand_o(no,l) = sand_i(mapunit_o(no),l)
             clay_o(no,l) = clay_i(mapunit_o(no),l)
          end do
       else
          ! use loam
          do l = 1, nlay
             sand_o(no,l) = 43.
             clay_o(no,l) = 18.
          end do
       end if
    end do

    ! -----------------------------------------------------------------
    ! Error check2
    ! Compare global area of each soil type on input and output grids
    ! -----------------------------------------------------------------

    ! input grid: global areas by texture class
    
!     gast_i(:) = 0.
!     do l = 1, nlay
!        do ni = 1,ns_i
!           mapunittemp = nint(mapunit_i(ni))
!           if (mapunittemp==0) then
!              typ = 'no soil: ocean, glacier, lake, no data'
!           else if (clay_i(mapunittemp,l) >= 40.) then
!              typ = 'clays'
!           else if (sand_i(mapunittemp,l) >= 50.) then
!              typ = 'sands'
!           else if (clay_i(mapunittemp,l)+sand_i(mapunittemp,l) < 50.) then
!              if (mask(ni) /= 0.) then
!                 typ = 'silts'
!              else            !if (mask(ni) == 0.) then no data
!                 typ = 'no soil: ocean, glacier, lake, no data'
!              end if
!           else
!              typ = 'loams'
!           end if
!           do m = 0, nlsm
!              if (typ == soil(m)) go to 101
!           end do
!           write (6,*) 'MKSOILTEX error: sand = ',sand_i(mapunittemp,l), &
!                ' clay = ',clay_i(mapunittemp,l), &
!                ' not assigned to soil type for input grid lon,lat,layer = ',ni,l
!           call shr_sys_abort()
! 101       continue
!           gast_i(m) = gast_i(m) + area_src(ni)*mask(ni)*re**2
!        end do
!     end do

!     ! output grid: global areas by texture class

!     gast_o(:) = 0.
!     do l = 1, nlay
!        do no = 1,ns_o
!           if (clay_o(no,l)==0. .and. sand_o(no,l)==0.) then
!              typ = 'no soil: ocean, glacier, lake, no data'
!           else if (clay_o(no,l) >= 40.) then
!              typ = 'clays'
!           else if (sand_o(no,l) >= 50.) then
!              typ = 'sands'
!           else if (clay_o(no,l)+sand_o(no,l) < 50.) then
!              typ = 'silts'
!           else
!              typ = 'loams'
!           end if
!           do m = 0, nlsm
!              if (typ == soil(m)) go to 102
!           end do
!           write (6,*) 'MKSOILTEX error: sand = ',sand_o(no,l), &
!                ' clay = ',clay_o(no,l), &
!                ' not assigned to soil type for output grid lon,lat,layer = ',no,l
!           call shr_sys_abort()
! 102       continue
!           gast_o(m) = gast_o(m) + area_dst(no)*frac_o(no)*re**2
!        end do
!     end do

!     ! Diagnostic output

!     write (ndiag,*)
!     write (ndiag,'(1x,70a1)') ('=',l=1,70)
!     write (ndiag,*) 'Soil Texture Output'
!     write (ndiag,'(1x,70a1)') ('=',l=1,70)
!     write (ndiag,*)

!     write (ndiag,*) 'The following table of soil texture classes is for comparison only.'
!     write (ndiag,*) 'The actual data is continuous %sand, %silt and %clay not textural classes'
!     write (ndiag,*)

!     write (ndiag,*)
!     write (ndiag,'(1x,70a1)') ('.',l=1,70)
!     write (ndiag,1001)
! 1001 format (1x,'soil texture class',17x,' input grid area output grid area',/ &
!          1x,33x,'     10**6 km**2','      10**6 km**2')
!     write (ndiag,'(1x,70a1)') ('.',l=1,70)
!     write (ndiag,*)

!     do l = 0, nlsm
!        write (ndiag,1002) soil(l),gast_i(l)*1.e-6,gast_o(l)*1.e-6
! 1002   format (1x,a38,f16.3,f17.3)
!     end do

    ! Deallocate dynamic memory
    deallocate (sand_i,clay_i,mapunit_i)
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made %sand and %clay'
       write (ndiag,*)
    end if

  end subroutine mksoiltex

end module mksoiltexMod
