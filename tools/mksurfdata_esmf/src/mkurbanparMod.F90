module mkurbanparMod

  !-----------------------------------------------------------------------
  ! Make Urban Parameter data
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod     , only : r8 => shr_kind_r8, r4 => shr_kind_r4, cs => shr_kind_cs
  use shr_sys_mod      , only : shr_sys_abort
  use mkpioMod         , only : mkpio_get_rawdata, pio_iotype, pio_iosystem
  use mkpioMod         , only : mkpio_iodesc_output, mkpio_get_rawdata
  use mkesmfMod        , only : regrid_rawdata, create_routehandle_r8
  use mkutilsMod       , only : chkerr
  use mkvarctl         , only : root_task, ndiag, ispval, outnc_double
  use mkutilsMod       , only : normalize_classes_by_gcell
  use mkdiagnosticsMod , only : output_diagnostics_index, output_diagnostics_continuous
  use mkindexmapMod    , only : dim_slice_type, lookup_2d_netcdf
  use mkvarpar         , only : numrad, numsolar, re

  implicit none
  private

  public :: mkurbanInit
  public :: mkurban
  public :: mkurbanpar
  public :: update_max_array_urban
  public :: normalize_urbn_by_tot
  public :: mkurban_topo                   ! Get elevation to reduce urban for high elevation areas
  public :: mkurban_pct_diagnostics        ! print diagnostics related to pct urban

  ! Note: normalize_urbn_by_tot could be private, but because there
  ! are associated test routines in a separate module, it needs to be public

  ! public data members
  integer, public :: numurbl           ! number of urban classes
  integer, public :: nlevurb = ispval  ! number of urban layers
  integer, public :: nregions

  real(r8), parameter :: MIN_DENS = 0.1_r8 ! minimum urban density (% of grid cell) - below this value, urban % is set to 0

  ! private data members:
  ! flag to indicate nodata for index variables in output file:
  integer         , parameter :: index_nodata = 0
  real(r8)        , allocatable :: frac_o_mkurban_nonorm(:)
  type(ESMF_RouteHandle) :: routehandle_mkurban_nonorm
  character(len=*), parameter :: modname = 'mkurbanparMod'

  private :: index_nodata
  private :: modname

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkurbanInit(datafname)
    !
    ! Initialize variables needed for urban
    !
    ! input variables
    character(len=*), intent(in) :: datafname  ! input data file name (same as file used in mkurban)
    !
    ! local variables:
    type(file_desc_t) :: pioid
    integer           :: dimid
    integer           :: rcode
    character(len=*), parameter :: subname = 'mkurbanInit'
    !-----------------------------------------------------------------------

    ! Set numurbl, nlevurb and nregions
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(datafname), pio_nowrite)
    rcode = pio_inq_dimid(pioid, 'density_class', dimid)
    rcode = pio_inq_dimlen(pioid, dimid, numurbl)
    rcode = pio_inq_dimid(pioid, 'nlevurb', dimid)
    rcode = pio_inq_dimlen(pioid, dimid, nlevurb)
    rcode = pio_inq_dimid(pioid, 'region', dimid)
    rcode = pio_inq_dimlen(pioid, dimid, nregions)
    call pio_closefile(pioid)

  end subroutine mkurbanInit

  !===============================================================
  subroutine mkurban(file_mesh_i, file_data_i, mesh_o, pcturb_o, &
                     urban_classes_o, region_o, rc)
    !
    ! make total percent urban, breakdown into urban classes, and region ID on the output grid
    !
    ! urban_classes_o(n, i) gives the percent of the urban area in grid cell n that is in class #i.
    ! This is normalized so that sum(urban_classes_o(n,:)) = 100 for all n, even for grid
    ! cells where pcturb_o(n) = 0 (in the case where pcturb_o(n) = 0, we come up with an
    ! arbitrary assignment of urban into the different classes).
    !
    ! See comments under the normalize_urbn_by_tot subroutine for how urban_classes_o is
    ! determined when the total % urban is 0, according to the input data.
    !
    ! TODO (WJS 6-12-14): I think this could be rewritten slightly to take advantage of the
    ! new mkpctPftTypeMod (which should then be renamed to something more general; or maybe
    ! better, in terms of maintaining helpful abstractions, there could be a new type to
    ! handle urban, and both that and pct_pft_type could be build on a single set of shared
    ! code - either as a single base class or through a "has-a" mechanism). This would allow
    ! us to combine pcturb_o and urban_classes_o into a single derived type variable. I think
    ! this would also replace the use of normalize_classes_by_gcell, and maybe some other
    ! urban-specific code.
    !
    ! uses
    use mkinputMod, only: mksrf_fdynuse
    !
    ! input/output variables
    character(len=*) , intent(in)    :: file_mesh_i          ! input mesh file name
    character(len=*) , intent(in)    :: file_data_i          ! input data file name
    type(ESMF_Mesh)  , intent(in)    :: mesh_o               ! model mesh
    real(r8)         , intent(inout) :: pcturb_o(:)            ! output grid: total % urban
    real(r8)         , intent(inout) :: urban_classes_o(:,:) ! output grid: breakdown of total urban into each class
    integer          , intent(inout) :: region_o(:)          ! output grid: region ID
    integer          , intent(out)   :: rc

    ! local variables:
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: frac_i(:)
    real(r8), allocatable  :: data_i(:,:)
    real(r8), allocatable  :: data_o(:,:)
    real(r8), allocatable  :: urban_classes_gcell_o(:,:) ! % putput urban in each density class (% of total grid cell area)
    integer , allocatable  :: region_i(:)                ! input grid: region ID
    integer                :: ni,no                      ! indices
    integer                :: ns_i, ns_o                 ! array sizes
    integer                :: n,k,l                      ! indices
    integer                :: max_regions                ! maximum region index
    integer                :: max_index(1)
    integer                :: rcode, ier                 ! error status
    character(len=*), parameter :: subname = 'mkurban'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make %urban .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
    end if

    ! Open input data file
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Read in input mesh
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the landmask from the input data file and reset the mesh mask based on that
    allocate(frac_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'LANDMASK', mesh_i, frac_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (frac_i(ni) > 0._r8) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle between the input and output mesh
    if (.not. ESMF_RouteHandleIsCreated(routehandle_mkurban_nonorm)) then
       allocate(frac_o_mkurban_nonorm(ns_o))
       call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.false., &
            routehandle=routehandle_mkurban_nonorm, frac_o=frac_o_mkurban_nonorm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))
    end if

    ! Read in input data
    ! - levels are the outermost dimension in pio reads
    ! - levels are the innermost dimension for esmf fields
    ! Input data is read into (ns_i,numurbl) array and then transferred to data_i(numurbl,ns_i)
    allocate(data_i(numurbl,ns_i),stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(data_o(numurbl,ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'PCT_URBAN', mesh_i, data_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Regrid input data to model resolution
    !
    ! Use a nonorm mapper because we're mapping a field expressed as % of the grid cell area
    call regrid_rawdata(mesh_i, mesh_o, routehandle_mkurban_nonorm, data_i, data_o, 1, numurbl, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regrid_data for in "//trim(subname))

    ! Now Determine total % urban
    ! urbn_classes_gcell_o is % urban  of total grid cell area for each density class
    allocate(urban_classes_gcell_o(ns_o, numurbl), stat=ier)
    if (ier/=0) call shr_sys_abort()
    do l = 1,numurbl
       do no = 1,ns_o
          urban_classes_gcell_o(no,l) = data_o(l,no)
       end do
    end do
    do no = 1, ns_o
       pcturb_o(no) = sum(urban_classes_gcell_o(no,:))
    end do
    ! Check for conservation
    do no = 1, ns_o
       if ((pcturb_o(no)) > 100.000001_r8) then
          write (6,'(a,d13.5,a,i8)') trim(subname)//' error: percent urban = ',pcturb_o(no), &
               ' greater than 100.000001 for no = ',no
          call shr_sys_abort()
       end if
    enddo

    ! Determine urban_classes_o
    ! Make percent urban on output grid, given percent urban on input grid
    ! Determine pcturb_o on ouput grid:
    call normalize_urbn_by_tot(urban_classes_gcell_o, pcturb_o, urban_classes_o)
    call ESMF_LogWrite("After normalize_urbn", ESMF_LOGMSG_INFO)

    ! Handle special cases
    ! Note that, for all these adjustments of total urban %, we do not change the breakdown
    ! into the different urban classes. In particular: when pcturb_o is set to 0 for a point,
    ! the breakdown into the different urban classes is maintained as it was before.
    ! Set points to 0% if they fall below a given threshold
    do no = 1, ns_o
       if (pcturb_o(no) < MIN_DENS) then
          pcturb_o(no) = 0._r8
       end if
    end do

    ! Print diagnostics
    ! TODO: call to mkurban_pct_diagnostics has to be rewritten
    ! First, recompute urban_classes_gcell_o, based on any changes we have made to pcturb_o
    ! call normalize_classes_by_gcell(urban_classes_o, pcturb_o, urban_classes_gcell_o)
    ! do k = 1, numurbl
    !    call mkurban_pct_diagnostics(ldomain, tdomain, tgridmap, &
    !         urban_classes_gcell_i(:,k), urban_classes_gcell_o(:,k), &
    !         ndiag, dens_class=k, frac_dst=frac_dst)
    ! end do

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made %urban'
    end if

    ! ------------------------------------------------------
    ! Read in region field
    ! ------------------------------------------------------

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make urban region .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
    end if

    ! Read in region_i
    allocate(region_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'REGION_ID', mesh_i, region_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite("After reading in region_id in "//trim(subname), ESMF_LOGMSG_INFO)
    max_regions = nregions
    if (root_task) then
       write(ndiag,'(a,i8)')" max urban regions = ",max_regions
    end if

    ! Create a multi-dimensional array (data_i) where each ungridded dimension corresponds to a 2d field
    ! where there is a 1 for every gridcell that has region_i equal to a given region
    if (allocated(data_i)) deallocate(data_i)
    allocate(data_i(max_regions, ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort('error allocating data_i(max_regions, ns_i)')
    data_i(:,:) = 0._r8
    do l = 1,max_regions
       do ni = 1,ns_i
          if (region_i(ni) == l) then
             data_i(l,ni) = 1._r8 * frac_i(ni)
          end if
       end do
    end do
    deallocate(frac_i)

    ! Regrid data_i to data_o
    if (allocated(data_o)) deallocate(data_o)
    allocate(data_o(max_regions, ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort('error allocating data_i(max_regions, ns_o)')
    ! This regridding could be done either with or without fracarea normalization,
    ! because we just use it to find a dominant value. We use nonorm because we already
    ! have a nonorm mapper for the sake of PCTURB and this way we don't need to make a
    ! separate mapper with fracarea normalization.
    call regrid_rawdata(mesh_i, mesh_o, routehandle_mkurban_nonorm, data_i, data_o, 1, nregions, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Now find dominant region in each output gridcell - this is identical to the maximum index
    region_o(:) = 0
    do no = 1,ns_o
       max_index = maxloc(data_o(:,no))
       if (data_o(max_index(1),no) > 0._r8) then
          region_o(no) = max_index(1)
       end if
    end do

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made urban region'
       write (ndiag,*)
    end if

    ! Output diagnostics
    call output_diagnostics_index(mesh_i, mesh_o, mask_i, frac_o_mkurban_nonorm, &
         1, max_regions, region_i, region_o, 'urban region', ndiag, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()


    ! Close the file
    call pio_closefile(pioid)

    ! Deallocate dynamic memory & other clean up
    ! TODO: determine what to deallocate
    ! deallocate (urban_classes_gcell_i, urban_classes_gcell_o, region_i)
    if (mksrf_fdynuse == ' ') then  ! ...else we will reuse it
       deallocate(frac_o_mkurban_nonorm)
       call ESMF_VMLogMemInfo("Before destroy operation in "//trim(subname))
       call ESMF_RouteHandleDestroy(routehandle_mkurban_nonorm, nogarbage = .true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    end if
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
    real(r8), intent(in)   :: classes_pct_gcell(:,:) ! % cover of classes as % of grid cell
    real(r8), intent(in)   :: sums(:)                ! totals, as % of grid cell
    real(r8), intent(inout):: classes_pct_tot(:,:)   ! % cover of classes as % of total

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
  subroutine mkurbanpar(datfname_i, pioid_o, mesh_o, &
       region_o, urban_classes_gcell_o, urban_skip_abort_on_invalid_data_check)
    !
    ! Make Urban Parameter data using the information from region_o input
    !
    ! Note that, in a grid cell with region_o==r, parameter values are filled from region r
    ! for ALL density classes. Thus, the parameter variables have a numurbl dimension along
    ! with their other dimensions.
    !
    ! Note that we will have a 'nodata' value (given by the fill_val value associated with
    ! each parameter) wherever (1) we have a nodata value for region_o, or (2) the parameter
    ! has nodata for the given region/density combination in the input lookup table.
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: datfname_i                 ! input data file name
    type(file_desc_t) , intent(inout) :: pioid_o                    ! output file pio id
    type(ESMF_Mesh)   , intent(in)    :: mesh_o                     ! model mesh
    integer           , intent(in)    :: region_o(:)                ! output grid: region ID (length: ns_o)
    real(r8)          , intent(in)    :: urban_classes_gcell_o(:,:) ! output grid: percent urban in each density class
                                                                    ! (% of total grid cell area) (dimensions: ns_o, numurbl)
    logical           , intent(in)    :: urban_skip_abort_on_invalid_data_check

    ! local variables
    ! Type to store information about each urban parameter
    type param
       character(len=32) :: name          ! name in input & output files
       real(r8)          :: fill_val      ! value to put where we have no data in output
       logical           :: check_invalid ! should we check whether there are any invalid data in the output?
    end type param
    integer , allocatable :: idata_scalar_o(:,:)  ! output array for parameters with no extra dimensions
    real(r8), allocatable :: data_scalar_o(:,:)   ! output array for parameters with no extra dimensions
    real(r8), allocatable :: data_rad_o(:,:,:,:)  ! output array for parameters dimensioned by numrad & numsolar
    real(r8), allocatable :: data_levurb_o(:,:,:) ! output array for parameters dimensioned by nlevurb
    integer , allocatable :: unity_dens_o(:,:)    ! artificial density indices
    integer               :: nlevurb_i            ! input  grid: number of urban vertical levels
    integer               :: numsolar_i           ! input  grid: number of solar type (DIR/DIF)
    integer               :: numrad_i             ! input  grid: number of solar bands (VIS/NIR)
    integer               :: m,n,no,ns_o,p,k      ! indices
    type(file_desc_t)     :: pioid_i
    type(var_desc_t)      :: pio_varid
    type(io_desc_t)       :: pio_iodesc
    integer               :: pio_vartype
    integer               :: dimid
    integer               :: ier, rcode, rc       ! error status
    character(len=cs)     :: varname              ! variable name
    integer               :: xtype                ! external type

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

    if (root_task) then
       write (ndiag,'(a)') 'Attempting to make Urban Parameters .....'
    end if

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

    ! ------------------------------------------------
    ! Handle urban parameters with no extra dimensions
    ! ------------------------------------------------

    allocate(data_scalar_o(ns_o, numurbl), stat=ier)
    data_scalar_o(:,:) = 0._r8
    if (ier /= 0) call shr_sys_abort('mkurbanpar allocation error')

    do p = 1, size(params_scalar)
       ! get variable output (data_scalar_o)
       call lookup_and_check_err(pioid_i, params_scalar(p)%name, params_scalar(p)%fill_val, &
            params_scalar(p)%check_invalid, urban_skip_abort_on_invalid_data_check, &
            data_scalar_o, n_extra_dims = 0)

       ! get io descriptor for variable output, write out variable and free memory for io descriptor
       call mkpio_iodesc_output(pioid_o, mesh_o, params_scalar(p)%name, pio_iodesc, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//&
            trim(params_scalar(p)%name))
       rcode = pio_inq_varid(pioid_o, params_scalar(p)%name, pio_varid)
       rcode = pio_inq_vartype(pioid_o, pio_varid, pio_vartype)
       if (pio_vartype == PIO_INT) then
          allocate(idata_scalar_o(ns_o, numurbl))
          idata_scalar_o(:,:) = int(data_scalar_o)
          call pio_write_darray(pioid_o, pio_varid, pio_iodesc, idata_scalar_o(:,:), rcode)
          deallocate(idata_scalar_o)
       else
          call pio_write_darray(pioid_o, pio_varid, pio_iodesc, data_scalar_o(:,:), rcode)
       end if
       call pio_freedecomp(pioid_o, pio_iodesc)
    end do

    deallocate(data_scalar_o)

    ! ------------------------------------------------
    ! Handle urban parameters dimensioned by numrad & numsolar
    ! ------------------------------------------------

    allocate(data_rad_o(ns_o, numurbl, numrad, numsolar), stat=ier)
    if (ier /= 0) call shr_sys_abort('mkurbanpar allocation error for data_rad_o')

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
                  data_rad_o(:,:,n,m), n_extra_dims=2, extra_dims=extra_dims)
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

          ! get io descriptor for variable output, write out variable and free memory for io descriptor
          call mkpio_iodesc_output(pioid_o, mesh_o, varname, pio_iodesc, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//&
               trim(params_scalar(p)%name))
          rcode = pio_inq_varid(pioid_o, varname, pio_varid)
          call pio_write_darray(pioid_o, pio_varid, pio_iodesc, data_rad_o(:,:,:,m), rcode)
          call pio_freedecomp(pioid_o, pio_iodesc)
       end do

    end do

    deallocate(data_rad_o)
    deallocate(extra_dims)

    ! ------------------------------------------------
    ! Handle urban parameters dimensioned by nlevurb
    ! ------------------------------------------------

    allocate(data_levurb_o(ns_o, numurbl, nlevurb), stat=ier)
    if (ier /= 0) call shr_sys_abort('mkurbanpar allocation error for data_levurb_o')

    allocate(extra_dims(1))
    extra_dims(1)%name = 'nlevurb'

    do p = 1, size(params_levurb)
       do n = 1,nlevurb
          extra_dims(1)%val = n
          call lookup_and_check_err(pioid_i, params_levurb(p)%name, params_levurb(p)%fill_val, &
               params_levurb(p)%check_invalid, &
               urban_skip_abort_on_invalid_data_check, data_levurb_o(:,:,n), &
               n_extra_dims=1, extra_dims=extra_dims)
       end do

       ! get io descriptor for variable output, write out variable and free memory for io descriptor
       call mkpio_iodesc_output(pioid_o, mesh_o, params_levurb(p)%name, pio_iodesc, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in generating an iodesc for '//&
            trim(params_levurb(p)%name))
       rcode = pio_inq_varid(pioid_o, params_levurb(p)%name, pio_varid)
       call pio_write_darray(pioid_o, pio_varid, pio_iodesc, data_levurb_o(:,:,:), rcode)
       call pio_freedecomp(pioid_o, pio_iodesc)
    end do

    deallocate(data_levurb_o)
    deallocate(extra_dims)

    ! Close input data file
    call pio_closefile(pioid_i)

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
      !
      ! input/output variables
      type(file_desc_t) , intent(inout) :: pioid
      character(len=*)  , intent(in)    :: varname       ! name of lookup table
      real(r8)          , intent(in)    :: fill_val      ! value to put where we have no data in output variables
      logical           , intent(in)    :: check_invalid ! should we check whether there are any invalid data in the output?
      logical           , intent(in)    :: urban_skip_abort_on_invalid_data_check
      real(r8)          , intent(inout) :: data(:,:)     ! output from lookup_2d_netcdf
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

         call lookup_2d_netcdf(pioid = pioid, &
                               tablename = varname, &
                               lookup_has_invalid = .true., &
                               dimname1 = 'density_class', &
                               dimname2 = 'region', &
                               n_extra_dims = n_extra_dims, &
                               index1 = unity_dens_o(:,k), &
                               index2 = region_o, &
                               fill_val = fill_val, &
                               data = data(:,k), &
                               ierr = ierr, &
                               extra_dims= extra_dims, &
                               nodata = index_nodata, &
                               invalid_okay = .true.)
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
                     ! /glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/surfdata_map/README_c141219
                     call shr_sys_abort()
                  end if
               end if
            end do
         end if
      end do
    end subroutine lookup_and_check_err

  end subroutine mkurbanpar

  !===============================================================
  subroutine mkurban_pct_diagnostics(area_i, area_o, mask_i, frac_o, urbn_i, urbn_o, dens_class)
    !
    ! print diagnostics related to pct urban
    ! Compare global areas on input and output grids
    !
    ! This is intended to be called after mkurban_pct, but is split out into a separate
    ! routine so that modifications to urbn_o can be made in between the two calls (e.g.,
    ! setting urbn_o to 0 wherever it is less than a certain threshold; the rules for doing
    ! this can't always be applied inline in mkurban_pct).
    !
    ! input/output variables
    real(r8)          , intent(in) :: area_i(:)
    real(r8)          , intent(in) :: area_o(:)
    integer           , intent(in) :: mask_i(:)
    real(r8)          , intent(in) :: frac_o(:)
    real(r8)          , intent(in) :: urbn_i(:)  ! input grid: percent urban
    real(r8)          , intent(in) :: urbn_o(:)  ! output grid: percent urban
    integer , intent(in), optional :: dens_class ! density class

    ! local variables:
    real(r8) :: gurbn_i ! input  grid: global urbn
    real(r8) :: garea_i ! input  grid: global area
    real(r8) :: gurbn_o ! output grid: global urbn
    real(r8) :: garea_o ! output grid: global area
    integer  :: ni,no,k ! indices
    character(len=*), parameter :: subname = 'mkurban_pct_diagnostics'
    !-----------------------------------------------------------------------

    ! Input grid
    gurbn_i = 0._r8
    garea_i = 0._r8
    do ni = 1, size(area_i)
       garea_i = garea_i + area_i(ni)*re**2
       gurbn_i = gurbn_i + urbn_i(ni)*(area_i(ni)/100._r8)* mask_i(ni)*re**2
    end do

    ! Output grid
    gurbn_o = 0._r8
    garea_o = 0._r8
    do no = 1, size(area_o)
       garea_o = garea_o + area_o(no)*re**2
       gurbn_o = gurbn_o + urbn_o(no)* (area_o(no)/100._r8)*frac_o(no)*re**2
    end do

    ! Diagnostic output
    write (ndiag,*)
    write (ndiag,'(1x,70a1)') ('=',k=1,70)
    if (present(dens_class)) then
       write (ndiag,'(1x,a,i0)') 'Urban Output -- class ', dens_class
    else
       write (ndiag,'(1x,a)') 'Urban Output'
    end if
    write (ndiag,'(1x,70a1)') ('=',k=1,70)
    write (ndiag,*)
    write (ndiag,'(1x,70a1)') ('.',k=1,70)
    write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/&
             1x,'                 10**6 km**2      10**6 km**2   ')
    write (ndiag,'(1x,70a1)') ('.',k=1,70)
    write (ndiag,*)
    write (ndiag,2003) gurbn_i*1.e-06,gurbn_o*1.e-06
    write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'urban       ',f14.3,f17.3)
2003 format (1x,'urban       ',f14.3,f22.8)
2004 format (1x,'all surface ',f14.3,f17.3)

  end subroutine mkurban_pct_diagnostics

  !===============================================================
  subroutine mkurban_topo(file_mesh_i, file_data_i, mesh_o, varname, elev_o, rc)
    !
    ! Make elevation data
    !
    ! input/output variables
    character(len=*) , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*) , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)  , intent(in)    :: mesh_o      ! model mesh
    character(len=*) , intent(in)    :: varname     ! topo variable name
    real(r8)         , intent(inout) :: elev_o(:)   ! output elevation data
    integer          , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle)      :: routehandle
    type(ESMF_Mesh)             :: mesh_i
    type(file_desc_t)           :: pioid
    real(r8), allocatable       :: frac_o(:)
    real(r8), allocatable       :: data_i(:,:)
    real(r8), allocatable       :: data_o(:,:)
    real(r8), allocatable       :: elev_i(:)  ! canyon_height to width ratio in
    integer                     :: ns_i,ns_o  ! bounds
    integer                     :: ni, no     ! indices
    integer                     :: k,l,n,m    ! indices
    character(len=CS)           :: name       ! name of attribute
    character(len=CS)           :: unit       ! units of attribute
    integer                     :: ier,rcode        ! error status
    character(len=*), parameter :: subname = 'mkelev'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make urban topo elevation .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
    end if

    ! Open input data file
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Read in input mesh
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read topo elev dataset with unit mask everywhere
    allocate(elev_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, trim(varname), mesh_i, elev_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call pio_closefile(pioid)

    ! Create a route handle between the input and output mesh
    allocate(frac_o(ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.true., &
         routehandle=routehandle, frac_o=frac_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Regrid input data to model resolution - determine elev_o on output grid
    elev_o(:) = 0.
    if (ier/=0) call shr_sys_abort()
    call regrid_rawdata(mesh_i, mesh_o, routehandle, elev_i, elev_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call output_diagnostics_continuous(mesh_i, mesh_o, elev_i, elev_o, &
       "Urban elev variable", "m", ndiag=ndiag, rc=rc, nomask=.true.)

    ! Deallocate dynamic memory
    deallocate (elev_i)
    deallocate (frac_o)
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made elevation'
       write (ndiag,'(a)')
    end if

  end subroutine mkurban_topo

  !===============================================================
  subroutine update_max_array_urban(pct_urbmax_arr,pct_urban_arr)
    !
    ! !DESCRIPTION:
    ! Update the maximum percent cover of each urban class for landuse.timeseries file
    !
    ! !ARGUMENTS:
    real(r8)         , intent(inout):: pct_urbmax_arr(:,:)           ! max percent cover of each urban class
    real(r8)         , intent(in):: pct_urban_arr(:,:)           ! percent cover of each urban class that is used to update the old pct_urbmax_arr
    !
    ! !LOCAL VARIABLES:
    integer :: n,k,ns              ! indices

    character(len=*), parameter :: subname = 'update_max_array_urban'
    !-----------------------------------------------------------------------
    ns = size(pct_urban_arr,1)
    do n = 1, ns
       do k =1, numurbl
          if (pct_urban_arr(n,k) > pct_urbmax_arr(n,k)) then
             pct_urbmax_arr(n,k) = pct_urban_arr(n,k)
          end if
       end do
    end do

  end subroutine update_max_array_urban

end module mkurbanparMod
