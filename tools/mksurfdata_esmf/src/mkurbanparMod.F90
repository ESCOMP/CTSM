module mkurbanparMod

  !-----------------------------------------------------------------------
  ! Make Urban Parameter data
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4, cs => shr_kind_cs
  use shr_sys_mod  , only : shr_sys_abort
  use mkpioMod     , only : mkpio_get_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkpioMod     , only : mkpio_iodesc_output, mkpio_def_spatial_var, mkpio_wopen
  use mkpioMod     , only : mkpio_get_dimlengths, mkpio_get_rawdata
  use mkpioMod     , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod    , only : regrid_rawdata, create_routehandle_r8
  use mkutilsMod   , only : chkerr
  use mkvarctl     , only : root_task, ndiag, ispval, fsurdat, outnc_double
  use mkvarctl     , only : root_task, ndiag, mpicom, MPI_INTEGER, MPI_MAX

  implicit none
  private

  public :: mkurbanInit
  public :: mkurban
  public :: mkurbanpar
  public :: normalize_urbn_by_tot

  ! Note: normalize_urbn_by_tot could be private, but because there
  ! are associated test routines in a separate module, it needs to be public

  ! public data members
  integer, public :: numurbl           ! number of urban classes
  integer, public :: nlevurb = ispval  ! number of urban layers
  integer, public :: nregions

  ! private data members:
  ! flag to indicate nodata for index variables in output file:
  integer         , parameter :: index_nodata = 0
  character(len=*), parameter :: modname = 'mkurbanparMod'

  private :: index_nodata
  private :: modname

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkurbanInit(datfname)
    !
    ! Initialize variables needed for urban
    !
    ! input variables
    character(len=*), intent(in) :: datfname  ! input data file name (same as file used in mkurban)
    !
    ! local variables:
    type(file_desc_t) :: pioid
    integer           :: dimid
    integer           :: rcode
    character(len=*), parameter :: subname = 'mkurbanInit'
    !-----------------------------------------------------------------------

    ! Set numurbl and nlevurb
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(datfname), pio_nowrite)
    rcode = pio_inq_dimid(pioid, 'density_class', dimid)
    rcode = pio_inq_dimlen(pioid, dimid, numurbl)
    rcode = pio_inq_dimid(pioid, 'nlevurb', dimid)
    rcode = pio_inq_dimlen(pioid, dimid, dimid)
    rcode = pio_inq_dimid(pioid, 'region', dimid)
    rcode = pio_inq_dimlen(pioid, dimid, nregions)
    call pio_closefile(pioid)

  end subroutine mkurbanInit

  !===============================================================
  subroutine mkurban(file_mesh_i, file_data_i, mesh_o, &
       zero_out, urbn_o, urban_classes_o, region_o, rc)
    !
    ! make total percent urban, breakdown into urban classes, and region ID on the output grid
    !
    ! urban_classes_o(n, i) gives the percent of the urban area in grid cell n that is in class #i.
    ! This is normalized so that sum(urban_classes_o(n,:)) = 100 for all n, even for grid
    ! cells where urbn_o(n) = 0 (in the case where urbn_o(n) = 0, we come up with an
    ! arbitrary assignment of urban into the different classes).
    !
    ! See comments under the normalize_urbn_by_tot subroutine for how urban_classes_o is
    ! determined when the total % urban is 0, according to the input data. Note that this
    ! also applies when all_urban=.true., for points that have 0 urban according to the input
    ! data.
    !
    ! TODO (WJS 6-12-14): I think this could be rewritten slightly to take advantage of the
    ! new mkpctPftTypeMod (which should then be renamed to something more general; or maybe
    ! better, in terms of maintaining helpful abstractions, there could be a new type to
    ! handle urban, and both that and pct_pft_type could be build on a single set of shared
    ! code - either as a single base class or through a "has-a" mechanism). This would allow
    ! us to combine urbn_o and urban_classes_o into a single derived type variable. I think
    ! this would also replace the use of normalize_classes_by_gcell, and maybe some other
    ! urban-specific code.
    !
    use mkindexmapMod       , only : get_dominant_indices
#ifdef TODO
    use mkurbanparCommonMod , only : mkurban_pct_diagnostics
#endif
    use mkurbanparCommonMod , only : MIN_DENS
    use mkutilsMod          , only : normalize_classes_by_gcell
    use mkvarctl            , only : all_urban
    use mkvarpar
#ifdef TODO
    use mkdiagnosticsMod    , only : output_diagnostics_index
#endif

    ! input/output variables
    character(len=*) , intent(in)    :: file_mesh_i          ! input mesh file name
    character(len=*) , intent(in)    :: file_data_i          ! input data file name
    type(ESMF_Mesh)  , intent(in)    :: mesh_o               ! model mesh
    logical          , intent(in)    :: zero_out             ! if should zero urban out
    real(r8)         , intent(out)   :: urbn_o(:)            ! output grid: total % urban
    real(r8)         , intent(out)   :: urban_classes_o(:,:) ! output grid: breakdown of total urban into each class
    integer          , intent(inout) :: region_o(:)          ! output grid: region ID
    integer          , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid
    real(r8), allocatable  :: urban_classes_gcell_i(:,:) ! % input urban in each density class (% of total grid cell area)
    real(r8), allocatable  :: urban_classes_gcell_o(:,:) ! % putput urban in each density class (% of total grid cell area)
    real(r8), allocatable  :: data_i(:,:)
    real(r8), allocatable  :: data_o(:,:)
    integer , allocatable  :: region_i(:)               ! input grid: region ID
    integer                :: ni,no                     ! indices
    integer                :: n,k,l                     ! indices   
    integer                :: ns_i, ns_o                ! array sizes
    integer                :: dimlen                    ! netCDF dimension length
    integer, allocatable   :: dimlengths(:)
    integer                :: max_region                ! maximum region index
    real(r4), allocatable  :: frac_i(:)
    real(r4), allocatable  :: frac_o(:)
    integer                :: max_regions_local
    integer                :: max_regions
    integer                :: max_index(1)
    integer                :: rcode, ier                ! error status
    character(len=*), parameter :: subname = 'mkurban'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write (ndiag,'(a)') 'Attempting to make %urban .....'
    end if

    ! Open raw data file - need to do this first to obtain ungridded dimension size
    if (root_task) then
       write (ndiag,'(a)') 'Opening urban file: '//trim(file_data_i)
    end if
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Create a route handle between the input and output mesh
    call create_routehandle_r8(mesh_i, mesh_o, routehandle, rc)
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_i and allocate data_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(data_i(numurbl,ns_i),stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! Determine ns_o and allocate data_o
    ns_o = size(urbn_o, dim=1)
    if (ier/=0) call shr_sys_abort()
    allocate(data_o(numurbl,ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! Read in input data
    ! - levels are the innermost dimension for esmf fields
    ! - levels are the outermost dimension in pio reads
    ! Input data is read into (ns_i,numurbl) array and then transferred to urban_classes_gcell_i(numurbl,ns_i)
    call mkpio_get_rawdata(pioid, 'PCT_URBAN', mesh_i, urban_classes_gcell_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Regrid raw data to model resolution
    ! make percent urban on output grid, given percent urban on input grid
    ! This assumes that we're neither using all_urban or zero_out
    ! Determine urbn_o on ouput grid:
    ! Area-average percent cover on input grid to output grid
    ! and correct according to land landmask
    ! Note that percent cover is in terms of total grid area.
    call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 1, numurbl, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regrid_data for in "//trim(subname))

    ! Now Determine total % urban
    call ESMF_LogWrite("Before allocate of urban_classes_gcell_o", ESMF_LOGMSG_INFO)
    allocate(urban_classes_gcell_o(ns_o, numurbl), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call ESMF_LogWrite("Before loop", ESMF_LOGMSG_INFO)
    do l = 1,numurbl
       do n = 1,ns_o
          urban_classes_gcell_o(n,l) = data_o(l,n)
       end do
    end do

    ! Determine frac_o (regrid frac_i to frac_o)
    allocate(frac_i(ns_i))
    allocate(frac_o(ns_o))
    call mkpio_get_rawdata(pioid, 'LANDMASK', mesh_i, frac_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call regrid_rawdata(mesh_i, mesh_o, routehandle, frac_i, frac_o, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regrid_data for landmask in "//trim(subname))

    do n = 1,ns_o
       if (frac_o(n) > 0._r8) then
          urban_classes_gcell_o(n,:) = urban_classes_gcell_o(n,:) / frac_o(n)
       else
          urban_classes_gcell_o(n,:) = 0._r8
       end if
    end do
    do n = 1, ns_o
       urbn_o(n) = sum(urban_classes_gcell_o(n,:))
    end do
    call ESMF_LogWrite("After loop", ESMF_LOGMSG_INFO)

    ! Check for conservation
#ifdef TODO
    do no = 1, ns_o
       if ((urbn_o(no)) > 100.000001_r8) then
          write (6,'(a,i8,a,i8)') 'MKURBAN error: percent urban = ',urbn_o(no), &
               ' greater than 100.000001 for no = ',no
          !call shr_sys_abort()
       end if
    enddo
    !call shr_sys_abort()
#endif

    ! determine urban_classes_o
    call ESMF_LogWrite("Before normalize_urbn", ESMF_LOGMSG_INFO)
    call normalize_urbn_by_tot(urban_classes_gcell_o, urbn_o, urban_classes_o)
    call ESMF_LogWrite("After normalize_urbn", ESMF_LOGMSG_INFO)

    ! Handle special cases
    ! Note that, for all these adjustments of total urban %, we do not change anything
    ! about the breakdown into the different urban classes. In particular: when urbn_o is
    ! set to 0 for a point, the breakdown into the different urban classes is maintained
    ! as it was before.
    if (all_urban) then
       urbn_o(:) = 100._r8
    else if (zero_out) then
       urbn_o(:) = 0._r8
    else
       ! Set points to 0% if they fall below a given threshold
       do no = 1, ns_o
          if (urbn_o(no) < MIN_DENS) then
             urbn_o(no) = 0._r8
          end if
       end do
    end if

    deallocate(data_i)
    deallocate(data_o)

#ifdef TODO
    ! Print diagnostics
    ! First, recompute urban_classes_gcell_o, based on any changes we have made to urbn_o
    call normalize_classes_by_gcell(urban_classes_o, urbn_o, urban_classes_gcell_o)
    do k = 1, numurbl
       call mkurban_pct_diagnostics(ldomain, tdomain, tgridmap, &
            urban_classes_gcell_i(:,k), urban_classes_gcell_o(:,k), &
            ndiag, dens_class=k, frac_dst=frac_dst)
    end do
#endif
    if (root_task) then
       write (ndiag,'(a)') 'Successfully made %urban'
       write (ndiag,'(a)') 'Attempting to make urban region .....'
    end if

    ! ------------------------------------------------------
    ! Read in region field
    ! ------------------------------------------------------

    allocate(region_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()

    call ESMF_LogWrite("Before mkpio_get_rawdata", ESMF_LOGMSG_INFO)
    call mkpio_get_rawdata(pioid, 'REGION_ID', mesh_i, region_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite("After mkpio_get_rawdata", ESMF_LOGMSG_INFO)

    max_regions = nregions

    ! Now determine data_i as a real 2d array
    allocate(data_i(max_regions, ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort('error allocating data_i(max_regions, ns_i)')
    allocate(data_o(max_regions, ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort('error allocating data_i(max_regions, ns_o)')

    data_i(:,:) = 0._r4
    do l = 1,max_regions
       do n = 1,ns_i
          if (region_i(n) == l) then
             !data_i(l,n) = 1._r4 * mask_i(n) 
             data_i(l,n) = 1. ! DEBUG
          end if
       end do
    end do

    ! Regrid data_i to data_o
    call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 1, nregions, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Now find dominant region in each output gridcell - this is identical to the maximum index
    do n = 1,ns_o
       max_index = maxloc(data_o(:,n))
    end do
    ! do n = 1,ns_o
    !    maxindex = 1
    !    maxval = data_o(1,n)
    !    do l = 2, max_regions
    !       if (data_o(l,n) > maxval) then
    !          maxindex = l
    !          maxval = data_o(l,n)
    !       end if
    !    end do
    ! end do

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made urban region'
       write (ndiag,*)
    end if

    ! Output diagnostics
    ! call output_diagnostics_index(factorlist, factorindex, region_i, region_o, 'Urban Region ID', &
    !      1, max_region, ndiag)

    ! Close the file
    call pio_closefile(pioid)

    ! Deallocate dynamic memory & other clean up
    ! TODO: determine what to deallocate
    ! deallocate (urban_classes_gcell_i, urban_classes_gcell_o, region_i)

    call ESMF_VMLogMemInfo("Before destroy operation in "//trim(subname))
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("after destroy operation in "//trim(subname))

  end subroutine mkurban

  !===============================================================
  subroutine normalize_urbn_by_tot(classes_pct_gcell, sums, classes_pct_tot)
    !
    ! Normalizes urban class areas to produce % cover of each class, as % of total urban area
    !
    ! Specifically: Given (1) an array specifying the % cover of each urban class, as a % of
    ! the total grid cell area ('classes_pct_gcell'), and (2) a vector giving the total urban
    ! area in each grid cell, expressed as % of the grid cell area: Returns an array
    ! ('classes_pct_tot') of the same dimensionality as classes_pct_gcell, where the values
    ! now give % cover of each class as a % of the total urban area.
    !
    ! Assumes that sums(n) = sum(classes_pct_gcell(n,:))
    !
    ! When sums(n) = 0, the creation of classes_pct_tot(n,:) is ambiguous. Here we use the
    ! rule that all area is assigned to the medium-density class, defined by parameter MD.
    !
    ! The returned array satisfies sum(classes_pct_tot(n,:))==100 for all n (within rounding error)

    ! input/output variables
    real(r8), intent(in) :: classes_pct_gcell(:,:) ! % cover of classes as % of grid cell
    real(r8), intent(in) :: sums(:)                ! totals, as % of grid cell
    real(r8), intent(out):: classes_pct_tot(:,:)   ! % cover of classes as % of total

    ! local variables:
    integer  :: n         ! index
    integer  :: n_max     ! number of points
    integer  :: nclasses  ! number of classes
    real(r8) :: suma      ! sum for error check

    ! index of medium-density class, which is where we assign urban areas when the total
    ! urban area is 0
    integer, parameter :: MD = 3

    ! relative error tolerance for error check
    real(r8), parameter :: relerr = 1.e-10_r8

    character(len=*), parameter :: subname = 'normalize_urbn_by_tot'
    !-----------------------------------------------------------------------

    ! Error-check inputs

    n_max = size(sums)
    if (size(classes_pct_tot, 1)   /= n_max .or. &
        size(classes_pct_gcell, 1) /= n_max) then
       write(6,*) subname//' ERROR: array size mismatch'
       write(6,*) 'size(sums)                 = ', n_max
       write(6,*) 'size(classes_pct_tot, 1)   = ', size(classes_pct_tot, 1)
       write(6,*) 'size(classes_pct_gcell, 1) = ', size(classes_pct_gcell, 1)
       call shr_sys_abort()
    end if

    if (size(classes_pct_tot, 2) /= size(classes_pct_gcell, 2)) then
       write(6,*) subname//' ERROR: array size mismatch'
       write(6,*) 'size(classes_pct_tot, 2)   = ', size(classes_pct_tot, 2)
       write(6,*) 'size(classes_pct_gcell, 2) = ', size(classes_pct_gcell, 2)
       call shr_sys_abort()
    end if

    nclasses = size(classes_pct_gcell, 2)
    if (MD > nclasses) then
       write(6,*) subname//' ERROR: MD exceeds nclasses'
       write(6,*) 'MD       = ', MD
       write(6,*) 'nclasses = ', nclasses
       call shr_sys_abort()
    end if

    ! Do the work

    do n = 1, n_max
       if (sums(n) > 0._r8) then
          classes_pct_tot(n,:) = classes_pct_gcell(n,:)/sums(n) * 100._r8
       else
          ! Creation of classes_pct_tot is ambiguous. Apply the rule that all area is
          ! assigned to the medium-density class.
          classes_pct_tot(n,:)  =   0._r8
          classes_pct_tot(n,MD) = 100._r8
       end if
    end do

    ! Error-check output: Make sure sum(classes_pct_tot(n,:)) = 100 for all n

    do n = 1, n_max
       suma = sum(classes_pct_tot(n,:))
       if (abs(suma/100._r8 - 1._r8) > relerr) then
          write(6,*) subname//' ERROR: sum does not equal 100 at point ', n
          write(6,*) 'suma = ', suma
          call shr_sys_abort()
       end if
    end do

  end subroutine normalize_urbn_by_tot

  !===============================================================
  subroutine mkurbanpar(datfname_i, mesh_o, &
       region_o, urban_classes_gcell_o, urban_skip_abort_on_invalid_data_check)
    !
    ! Make Urban Parameter data
    !
    ! Note that, in a grid cell with region_o==r, parameter values are filled from region r
    ! for ALL density classes. Thus, the parameter variables have a numurbl dimension along
    ! with their other dimensions.
    !
    ! Note that we will have a 'nodata' value (given by the fill_val value associated with
    ! each parameter) wherever (1) we have a nodata value for region_o, or (2) the parameter
    ! has nodata for the given region/density combination in the input lookup table.
    !
    use mkindexmapMod, only : dim_slice_type, lookup_2d_netcdf
    use mkvarpar
    !
    ! input/output variables
    character(len=*)  , intent(in) :: datfname_i                ! input data file name
    type(ESMF_Mesh)   , intent(in) :: mesh_o                    ! model mesh
    integer           , intent(in) :: region_o(:)               ! output grid: region ID (length: ns_o)
    real(r8)          , intent(in) :: urban_classes_gcell_o(:,:) ! output grid: percent urban in each density class
                                                                ! (% of total grid cell area) (dimensions: ns_o, numurbl)
    logical           , intent(in) :: urban_skip_abort_on_invalid_data_check

    ! local variables
    ! Type to store information about each urban parameter
    type param
       character(len=32) :: name          ! name in input & output files
       real(r8)          :: fill_val      ! value to put where we have no data in output
       logical           :: check_invalid ! should we check whether there are any invalid data in the output?
    end type param
    real(r8), allocatable      :: data_scalar_o(:,:)   ! output array for parameters with no extra dimensions
    real(r8), allocatable      :: data_rad_o(:,:,:,:)  ! output array for parameters dimensioned by numrad & numsolar
    real(r8), allocatable      :: data_levurb_o(:,:,:) ! output array for parameters dimensioned by nlevurb
    integer , allocatable      :: unity_dens_o(:,:)    ! artificial density indices
    integer                    :: nlevurb_i            ! input  grid: number of urban vertical levels
    integer                    :: numsolar_i           ! input  grid: number of solar type (DIR/DIF)
    integer                    :: numrad_i             ! input  grid: number of solar bands (VIS/NIR)
    integer                    :: m,n,no,ns_o,p,k      ! indices
    type(file_desc_t)          :: pioid_i
    type(file_desc_t)          :: pioid_o
    type(var_desc_t)           :: pio_varid
    integer                    :: pio_vartype
    type(io_desc_t)            :: pio_iodesc
    integer                    :: dimid 
    integer                    :: ier, rcode, rc       ! error status
    character(len=cs)          :: varname              ! variable name
    integer                    :: xtype                ! external type
    ! information on extra dimensions for lookup tables greater than 2-d:
    type(dim_slice_type), allocatable :: extra_dims(:)

    ! value to put where we have no data in output variables, for real-valued parameters
    real(r8), parameter :: fill_val_real = 0._r8

    ! To add a new urban parameter, simply add an element to one of the below lists
    ! (params_scalar, params_rad or params_levurb)

    ! Urban parameters with no extra dimensions
    type(param), parameter :: params_scalar(13) = &
         (/ param('CANYON_HWR'      , fill_val_real, .true.), &  ! 1
            param('EM_IMPROAD'      , fill_val_real, .true.), &  ! 2
            param('EM_PERROAD'      , fill_val_real, .true.), &  ! 3
            param('EM_ROOF'         , fill_val_real, .true.), &  ! 4
            param('EM_WALL'         , fill_val_real, .true.), &  ! 5
            param('HT_ROOF'         , fill_val_real, .true.), &  ! 6
            param('THICK_ROOF'      , fill_val_real, .true.), &  ! 7
            param('THICK_WALL'      , fill_val_real, .true.), &  ! 8
            param('T_BUILDING_MIN'  , fill_val_real, .true.), &  ! 9
            param('WIND_HGT_CANYON' , fill_val_real, .true.), &  ! 10
            param('WTLUNIT_ROOF'    , fill_val_real, .true.), &  ! 11
            param('WTROAD_PERV'     , fill_val_real, .true.), &  ! 12

            ! Note that NLEV_IMPROAD is written as an integer, meaning that type conversion occurs
            ! by truncation. Thus we expect the values in the NLEV_IMPROAD lookup table to be exact;
            ! e.g., if a value were 1.99999 rather than 2.0000, it would be written as 1 instead of 2
            ! Also note: we use fill_val=-1 rather than 0, because 0 appears in the lookup table
            param('NLEV_IMPROAD'    , -1            , .true.) /)  ! 13

    ! Urban parameters dimensioned by numrad & numsolar
    type(param), parameter :: params_rad(4) = &
         (/ param('ALB_IMPROAD' , fill_val_real, .true.), &  ! 1
            param('ALB_PERROAD' , fill_val_real, .true.), &  ! 2
            param('ALB_ROOF'    , fill_val_real, .true.), &  ! 3
            param('ALB_WALL'    , fill_val_real, .true.) /)  ! 4

    ! suffix for variables dimensioned by numsolar, for each value of numsolar:
    character(len=8), parameter :: solar_suffix(numsolar) = (/'_DIR', '_DIF'/)

    ! Urban parameters dimensioned by nlevurb
    type(param), parameter :: params_levurb(6) = &
         (/ param('TK_ROOF', fill_val_real, .true.), &    ! 1
            param('TK_WALL', fill_val_real, .true.), &    ! 2
            param('CV_ROOF', fill_val_real, .true.), &    ! 3
            param('CV_WALL', fill_val_real, .true.), &    ! 4

            ! Impervious road thermal conductivity and heat capacity have varying levels of
            ! data. Thus, we expect to find some missing values in the lookup table -- we
            ! do not want to treat that as an error -- thus, we set check_invalid=.false.
            param('CV_IMPROAD', fill_val_real, .false.), & ! 5
            param('TK_IMPROAD', fill_val_real, .false.) /) ! 6


    character(len=*), parameter :: subname = 'mkurbanpar'
    !-----------------------------------------------------------------------

    write (6,*) 'Attempting to make Urban Parameters .....'

    ! Determine & error-check array sizes
    ns_o = size(region_o)
    if (size(urban_classes_gcell_o, 1) /= ns_o) then
       write(6,*) modname//':'//subname//' ERROR: array size mismatch'
       write(6,*) 'size(region_o) = ', size(region_o)
       write(6,*) 'size(urban_classes_gcell_o, 1) = ', size(urban_classes_gcell_o, 1)
       call shr_sys_abort()
    end if
    if (size(urban_classes_gcell_o, 2) /= numurbl) then
       write(6,*) modname//':'//subname//' ERROR: array size mismatch'
       write(6,*) 'size(urban_classes_gcell_o, 2) = ', size(urban_classes_gcell_o, 2)
       write(6,*) 'numurbl = ', numurbl
    end if


    ! Read dimensions from input file
    if (root_task) then
       write (ndiag,'(a)') 'Opening input urban parameter file: '//trim(datfname_i)
    end if
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(datfname_i), pio_nowrite)
    rcode = pio_inq_dimid(pioid_i, 'nlevurb', dimid)
    rcode = pio_inq_dimlen(pioid_i, dimid, nlevurb_i)
    rcode = pio_inq_dimid(pioid_i, 'numsolar', dimid)
    rcode = pio_inq_dimlen(pioid_i, dimid, numsolar_i)
    rcode = pio_inq_dimid(pioid_i, 'numrad', dimid)
    rcode = pio_inq_dimlen(pioid_i, dimid, numrad_i)

    if (nlevurb_i /= nlevurb) then
       write(6,*)'MKURBANPAR: parameter nlevurb= ',nlevurb, &
            'does not equal input dataset nlevurb= ',nlevurb_i
       call shr_sys_abort()
    endif
    if (numsolar_i /= numsolar) then
       write(6,*)'MKURBANPAR: parameter numsolar= ',numsolar, &
            'does not equal input dataset numsolar= ',numsolar_i
       call shr_sys_abort()
    endif
    if (numrad_i /= numrad) then
       write(6,*)'MKURBANPAR: parameter numrad= ',numrad, &
            'does not equal input dataset numrad= ',numrad_i
       call shr_sys_abort()
    endif

    ! Create an array that will hold the density indices
    ! In a given grid cell, we output parameter values for all density classes, for the
    ! region of that grid cell. In order to do this while still using the lookup_2d
    ! routine, we create a dummy unity_dens_o array that contains the density values
    ! passed to the lookup routine.

    allocate(unity_dens_o(ns_o, numurbl))
    do k = 1, numurbl
       unity_dens_o(:,k) = k
    end do

    ! Open output file (this will go into define mode)
    call mkpio_wopen(trim(fsurdat), clobber=.false., pioid=pioid_o)

    if ( outnc_double ) then
       xtype = PIO_DOUBLE
    else
       xtype = PIO_REAL
    end if

    ! Define urbanoutput variables
    call mkpio_def_spatial_var(pioid_o, varname='CANYON_HWR', xtype=xtype, &
         lev1name='numurbl', &
         long_name='canyon height to width ratio', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='EM_IMPROAD', xtype=xtype, &
         lev1name='numurbl', &
         long_name='emissivity of impervious road', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='EM_PERROAD', xtype=xtype, &
         lev1name='numurbl', &
         long_name='emissivity of pervious road', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='EM_ROOF', xtype=xtype, &
         lev1name='numurbl', &
         long_name='emissivity of roof', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='EM_WALL', xtype=xtype, &
         lev1name='numurbl', &
         long_name='emissivity of wall', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='HT_ROOF', xtype=xtype, &
         lev1name='numurbl', &
         long_name='height of roof', units='meters')

    call mkpio_def_spatial_var(pioid_o, varname='THICK_ROOF', xtype=xtype, &
         lev1name='numurbl', &
         long_name='thickness of roof', units='meters')

    call mkpio_def_spatial_var(pioid_o, varname='THICK_WALL', xtype=xtype, &
         lev1name='numurbl', &
         long_name='thickness of wall', units='meters')

    call mkpio_def_spatial_var(pioid_o, varname='T_BUILDING_MIN', xtype=xtype, &
         lev1name='numurbl', &
         long_name='minimum interior building temperature', units='K')

    call mkpio_def_spatial_var(pioid_o, varname='WIND_HGT_CANYON', xtype=xtype, &
         lev1name='numurbl', &
         long_name='height of wind in canyon', units='meters')

    call mkpio_def_spatial_var(pioid_o, varname='WTLUNIT_ROOF', xtype=xtype, &
         lev1name='numurbl', &
         long_name='fraction of roof', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='WTROAD_PERV', xtype=xtype, &
         lev1name='numurbl', &
         long_name='fraction of pervious road', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='ALB_IMPROAD_DIR', xtype=xtype, &
         lev1name='numurbl', lev2name='numrad', &
         long_name='direct albedo of impervious road', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='ALB_IMPROAD_DIF', xtype=xtype, &
         lev1name='numurbl', lev2name='numrad', &
         long_name='diffuse albedo of impervious road', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='ALB_PERROAD_DIR', xtype=xtype, &
         lev1name='numurbl', lev2name='numrad', &
         long_name='direct albedo of pervious road', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='ALB_PERROAD_DIF', xtype=xtype, &
         lev1name='numurbl', lev2name='numrad', &
         long_name='diffuse albedo of pervious road', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='ALB_ROOF_DIR', xtype=xtype, &
         lev1name='numurbl', lev2name='numrad', &
         long_name='direct albedo of roof', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='ALB_ROOF_DIF', xtype=xtype, &
         lev1name='numurbl', lev2name='numrad', &
         long_name='diffuse albedo of roof', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='ALB_WALL_DIR', xtype=xtype, &
         lev1name='numurbl', lev2name='numrad', &
         long_name='direct albedo of wall', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='ALB_WALL_DIF', xtype=xtype, &
         lev1name='numurbl', lev2name='numrad', &
         long_name='diffuse albedo of wall', units='unitless')

    call mkpio_def_spatial_var(pioid_o, varname='TK_ROOF', xtype=xtype, &
         lev1name='numurbl', lev2name='nlevurb', &
         long_name='thermal conductivity of roof', units='W/m*K')

    call mkpio_def_spatial_var(pioid_o, varname='TK_WALL', xtype=xtype, &
         lev1name='numurbl', lev2name='nlevurb', &
         long_name='thermal conductivity of wall', units='W/m*K')

    call mkpio_def_spatial_var(pioid_o, varname='TK_IMPROAD', xtype=xtype, &
         lev1name='numurbl', lev2name='nlevurb', &
         long_name='thermal conductivity of impervious road', units='W/m*K')

    call mkpio_def_spatial_var(pioid_o, varname='CV_ROOF', xtype=xtype, &
         lev1name='numurbl', lev2name='nlevurb', &
         long_name='volumetric heat capacity of roof', units='J/m^3*K')

    call mkpio_def_spatial_var(pioid_o, varname='CV_WALL', xtype=xtype, &
         lev1name='numurbl', lev2name='nlevurb', &
         long_name='volumetric heat capacity of wall', units='J/m^3*K')

    call mkpio_def_spatial_var(pioid_o, varname='CV_IMPROAD', xtype=xtype, &
         lev1name='numurbl', lev2name='nlevurb', &
         long_name='volumetric heat capacity of impervious road', units='J/m^3*K')

    call mkpio_def_spatial_var(pioid_o, varname='NLEV_IMPROAD', xtype=PIO_INT, &
         lev1name='numurbl', &
         long_name='number of impervious road layers', units='unitless')

    ! End define model
    rcode = pio_enddef(pioid_o)
    
    ! ------------------------------------------------
    ! Handle urban parameters with no extra dimensions
    ! ------------------------------------------------

   allocate(data_scalar_o(ns_o, numurbl), stat=ier)
   if (ier /= 0) then
      write(6,*)'mkurbanpar allocation error'; call shr_sys_abort()
   end if
   
   do p = 1, size(params_scalar)
      ! get variable output (data_scalar_o)
      call lookup_and_check_err(pioid_i, params_scalar(p)%name, params_scalar(p)%fill_val, &
           params_scalar(p)%check_invalid, urban_skip_abort_on_invalid_data_check, &
           data_scalar_o, 0)

      ! get io descriptor for variable output 
      call mkpio_iodesc_output(pioid_o, mesh_o, params_scalar(p)%name, pio_iodesc, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//&
           trim(params_scalar(p)%name))

      ! write out variable
      rcode = pio_inq_varid(pioid_o, params_scalar(p)%name, pio_varid)
      call pio_write_darray(pioid_o, pio_varid, pio_iodesc, data_scalar_o(:,:), rcode)

      ! Free memory for io descriptor
      call pio_freedecomp(pioid_o, pio_iodesc)
   end do

   deallocate(data_scalar_o)

   ! ------------------------------------------------
   ! Handle urban parameters with no extra dimensions
   ! ------------------------------------------------

   ! Handle urban parameters dimensioned by numrad & numsolar

   allocate(data_rad_o(ns_o, numurbl, numrad, numsolar), stat=ier)
   if (ier /= 0) then
      write(6,*)'mkurbanpar allocation error'; call shr_sys_abort()
   end if

   allocate(extra_dims(2))
   extra_dims(1)%name = 'numrad'
   extra_dims(2)%name = 'numsolar'

   do p = 1, size(params_rad)

      ! Get variable output (data_rad_o)
      do m = 1,numsolar
         extra_dims(2)%val = m
         do n = 1,numrad
            extra_dims(1)%val = n
            call lookup_and_check_err(pioid_i, params_rad(p)%name, params_rad(p)%fill_val, &
                 params_rad(p)%check_invalid, urban_skip_abort_on_invalid_data_check, &
                 data_rad_o(:,:,n,m), 2, extra_dims)
         end do
      end do

      ! Special handling of numsolar: rather than outputting variables with a numsolar
      ! dimension, we output separate variables for each value of numsolar
      do m = 1,numsolar
         if (len_trim(params_rad(p)%name) + len_trim(solar_suffix(m)) > len(varname)) then
            write(6,*) 'variable name exceeds length of varname'
            write(6,*) trim(params_rad(p)%name)//trim(solar_suffix(m))
            call shr_sys_abort()
         end if

         ! Determine variable name
         varname = trim(params_rad(p)%name)//trim(solar_suffix(m))

         ! get io descriptor for variable output 
         call mkpio_iodesc_output(pioid_o, mesh_o, varname, pio_iodesc, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//&
              trim(params_scalar(p)%name))

         ! write out variable
         rcode = pio_inq_varid(pioid_o, varname, pio_varid)
         call pio_write_darray(pioid_o, pio_varid, pio_iodesc, data_rad_o(:,:,:,m), rcode)

         ! Free memory for io descriptor
         call pio_freedecomp(pioid_o, pio_iodesc)
      end do

   end do

   deallocate(data_rad_o)
   deallocate(extra_dims)

    ! ------------------------------------------------
    ! Handle urban parameters dimensioned by nlevurb
    ! ------------------------------------------------

    allocate(data_levurb_o(ns_o, numurbl, nlevurb), stat=ier)
    if (ier /= 0) then
       write(6,*)'mkurbanpar allocation error'; call shr_sys_abort()
    end if

    allocate(extra_dims(1))
    extra_dims(1)%name = 'nlevurb'

    do p = 1, size(params_levurb)
       do n = 1,nlevurb
          extra_dims(1)%val = n
          call lookup_and_check_err(pioid_i, params_levurb(p)%name, params_levurb(p)%fill_val, &
               params_levurb(p)%check_invalid, &
               urban_skip_abort_on_invalid_data_check, data_levurb_o(:,:,n), 1, extra_dims)
       end do

       ! get io descriptor for variable output 
       call mkpio_iodesc_output(pioid_o, mesh_o, params_levurb(p)%name, pio_iodesc, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//&
            trim(params_levurb(p)%name))

       ! write out variable
       rcode = pio_inq_varid(pioid_o, params_levurb(p)%name, pio_varid)
       call pio_write_darray(pioid_o, pio_varid, pio_iodesc, data_levurb_o(:,:,:), rcode)

       ! Free memory for io descriptor
       call pio_freedecomp(pioid_o, pio_iodesc)
    end do

    deallocate(data_levurb_o)
    deallocate(extra_dims)

    ! Close input data file and output surface dataset
    call pio_closefile(pioid_i)
    call pio_closefile(pioid_o)

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made Urban Parameters'
       write (ndiag,'(a)')
    end if

    deallocate(unity_dens_o)

  contains
    !------------------------------------------------------------------------------
    subroutine lookup_and_check_err(pioid, varname, fill_val, check_invalid, &
         urban_skip_abort_on_invalid_data_check, data, n_extra_dims, extra_dims)

      ! Wrapper to lookup_2d_netcdf: Loops over each density class, calling lookup_2d_netcdf
      ! with that density class and filling the appropriate slice of the data array. Also
      ! checks for any errors, aborting if there were any.
      !
      ! Note that the lookup_2d_netcdf routine is designed to work with a single value of
      ! each of the indices. However, we want to fill parameter values for ALL density
      ! classes. This is why we loop over density class in this routine.
      !
      ! Note: inherits a number of variables from the parent routine

      use mkindexmapMod, only : lookup_2d_netcdf

      ! input/output variables
      type(file_desc_t) , intent(inout) :: pioid
      character(len=*)  , intent(in)    :: varname       ! name of lookup table
      real(r8)          , intent(in)    :: fill_val      ! value to put where we have no data in output variables
      logical           , intent(in)    :: check_invalid ! should we check whether there are any invalid data in the output?
      logical           , intent(in)    :: urban_skip_abort_on_invalid_data_check
      real(r8)          , intent(out)   :: data(:,:)     ! output from lookup_2d_netcdf
      integer           , intent(in)    :: n_extra_dims  ! number of extra dimensions in the lookup table

      ! slice to use if lookup table variable has more than 2 dimensions:
      type(dim_slice_type), intent(in), optional :: extra_dims(:)

      ! Local variables:
      integer :: k,n   ! indices
      integer :: ierr  ! error return code
      !-----------------------------------------------------------------------

      do k = 1, numurbl
         ! In the following, note that unity_dens_o(:,k) has been constructed so that
         ! unity_dens_o(:,k)==k everywhere. Thus, we fill data(:,k) with the parameter
         ! values corresponding to density class k.
         ! Also note: We use invalid_okay=.true. because we fill all density classes,
         ! some of which may have invalid entries. Because doing so disables some error
         ! checking, we do our own error checking after the call.

         call lookup_2d_netcdf(pioid, varname, .true., &
              'density_class', 'region', n_extra_dims, &
              unity_dens_o(:,k), region_o, fill_val, data(:,k), ierr, &
              extra_dims=extra_dims, nodata=index_nodata, &
              invalid_okay=.true.)

         if (ierr /= 0) then
            write(6,*) modname//':'//subname//' ERROR in lookup_2d_netcdf for ', &
                 trim(varname), ' class', k, ': err=', ierr
            call shr_sys_abort()
         end if

         if (check_invalid) then
            ! Make sure we have valid parameter values wherever we have non-zero urban cover
            do n = 1, ns_o
               ! This check assumes that fill_val doesn't appear in any of the valid entries
               ! of the lookup table
               if (urban_classes_gcell_o(n,k) > 0. .and. data(n,k) == fill_val) then
                  write(6,*) modname//':'//subname//' ERROR: fill value found in output where urban cover > 0'
                  write(6,*) 'var: ', trim(varname)
                  write(6,*) 'class: ', k
                  write(6,*) 'n: ', n
                  write(6,*) 'region: ', region_o(n)
                  write(6,*) 'urban_classes_gcell_o(n,k): ', urban_classes_gcell_o(n,k)
                  if (.not. urban_skip_abort_on_invalid_data_check) then
                     ! NOTE(bja, 2015-01) added to work around a ?bug? noted in
                     ! /glade/p/cesm/cseg/inputdata/lnd/clm2/surfdata_map/README_c141219
                     call shr_sys_abort()
                  end if
               end if
            end do
         end if
      end do
    end subroutine lookup_and_check_err

  end subroutine mkurbanpar

end module mkurbanparMod
