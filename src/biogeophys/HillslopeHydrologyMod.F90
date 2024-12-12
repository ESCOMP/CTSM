module HillslopeHydrologyMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read geomorphological parameters for hillslope columns
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use spmdMod        , only : masterproc, iam
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  use clm_varctl     , only : use_hillslope_routing
  use decompMod      , only : bounds_type
  use clm_varcon     , only : rpi
  use HillslopeHydrologyUtilsMod, only : HillslopeSoilThicknessProfile_linear

  ! !PUBLIC TYPES:
  implicit none

  private
  save

  ! !PUBLIC MEMBER FUNCTIONS:
  public hillslope_properties_init
  public InitHillslope
  public SetHillslopeSoilThickness
  public HillslopeSoilThicknessProfile
  public HillslopeSetLowlandUplandPfts
  public HillslopeDominantLowlandPft
  public HillslopePftFromFile
  public HillslopeStreamOutflow
  public HillslopeUpdateStreamWater

  integer, public :: pft_distribution_method ! Method for distributing pfts across hillslope columns
  integer, public :: soil_profile_method     ! Method for varying soil thickness across hillslope columns

  ! Streamflow methods
  integer, public, parameter :: streamflow_manning = 0
  ! Pft distribution methods
  integer, public, parameter :: pft_standard  = 0
  integer, public, parameter :: pft_from_file = 1
  integer, public, parameter :: pft_uniform_dominant_pft = 2
  integer, public, parameter :: pft_lowland_dominant_pft = 3
  integer, public, parameter :: pft_lowland_upland = 4

  ! PRIVATE
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  integer, private, parameter :: soil_profile_uniform               = 0
  integer, private, parameter :: soil_profile_from_file             = 1
  integer, private, parameter :: soil_profile_set_lowland_upland    = 2
  integer, private, parameter :: soil_profile_linear                = 3

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine hillslope_properties_init(NLFilename)
    !
    ! DESCRIPTION
    ! read in hillslope hydrology veg/soil properties namelist variables
    !
    ! !USES:
    use abortutils      , only : endrun
    use fileutils       , only : getavu, relavu
    use spmdMod         , only : mpicom, masterproc
    use shr_mpi_mod     , only : shr_mpi_bcast
    use clm_varctl      , only : iulog
    use clm_nlUtilsMod  , only : find_nlgroup_name

    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !----------------------------------------------------------------------
    integer            :: nu_nml                     ! unit for namelist file
    integer            :: nml_error                  ! namelist i/o error flag
    character(len=*), parameter :: nmlname = 'hillslope_properties_inparm'
    character(*), parameter    :: subName = "('read_hillslope_properties_namelist')"
    ! Default values for namelist
    character(len=50) :: hillslope_pft_distribution_method = 'Standard'      ! pft distribution method string
    character(len=50) :: hillslope_soil_profile_method     = 'Uniform'       ! soil thickness distribution method string
    !-----------------------------------------------------------------------

! MUST agree with name in namelist and read statement
    namelist /hillslope_properties_inparm/  &
         hillslope_pft_distribution_method, &
         hillslope_soil_profile_method

    ! Read hillslope hydrology namelist
    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'hillslope_properties_inparm', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=hillslope_properties_inparm,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(subname // ':: ERROR reading hillslope properties namelist')
          end if
       else
          call endrun(subname // ':: ERROR reading hillslope properties namelist')
       end if
       close(nu_nml)
       call relavu( nu_nml )

       if (      trim(hillslope_pft_distribution_method) == 'Standard' ) then
          pft_distribution_method = pft_standard
       else if ( trim(hillslope_pft_distribution_method) == 'FromFile' ) then
          pft_distribution_method = pft_from_file
       else if ( trim(hillslope_pft_distribution_method) == 'DominantPftUniform') then
          pft_distribution_method = pft_uniform_dominant_pft
       else if ( trim(hillslope_pft_distribution_method) == 'DominantPftLowland') then
          pft_distribution_method = pft_lowland_dominant_pft
       else if ( trim(hillslope_pft_distribution_method) == 'PftLowlandUpland') then
          pft_distribution_method = pft_lowland_upland
       else
          call endrun(msg="ERROR bad value for hillslope_pft_distribution_method in "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if

       if (      trim(hillslope_soil_profile_method) == 'Uniform' ) then
          soil_profile_method = soil_profile_uniform
       else if ( trim(hillslope_soil_profile_method) == 'FromFile' ) then
          soil_profile_method = soil_profile_from_file
       else if ( trim(hillslope_soil_profile_method) == 'SetLowlandUpland' ) then
          soil_profile_method = soil_profile_set_lowland_upland
       else if ( trim(hillslope_soil_profile_method) == 'Linear') then
          soil_profile_method = soil_profile_linear
       else
          call endrun(msg="ERROR bad value for hillslope_soil_profile_method in "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if

    end if

    call shr_mpi_bcast(pft_distribution_method, mpicom)
    call shr_mpi_bcast(soil_profile_method, mpicom)

    if (masterproc) then

       write(iulog,*) ' '
       write(iulog,*) 'hillslope_properties settings:'
       write(iulog,*) '  hillslope_pft_distribution_method = ',hillslope_pft_distribution_method
       write(iulog,*) '  hillslope_soil_profile_method     = ',hillslope_soil_profile_method

    end if

  end subroutine hillslope_properties_init

  !-----------------------------------------------------------------------
  subroutine check_aquifer_layer()
    !
    ! !DESCRIPTION:
    ! Check whether use_hillslope and use_aquifer_layer are both set
    ! The use of use_hillslope is implied by the call to this function
    ! in InitHillslope, but explicitly compare here for clarity.
    !
    ! !USES:
    use clm_varctl              , only : use_hillslope
    use SoilWaterMovementMod    , only : use_aquifer_layer
    if (use_hillslope .and. use_aquifer_layer()) then
       write(iulog,*) ' ERROR: use_hillslope and use_aquifer_layer may not be used simultaneously'
       call endrun(msg=' ERROR: use_hillslope and use_aquifer_layer cannot both be set to true' // &
            errMsg(sourcefile, __LINE__))
    end if

  end subroutine check_aquifer_layer

  !-----------------------------------------------------------------------

  subroutine InitHillslope(bounds, hillslope_file)
    !
    ! !DESCRIPTION:
    ! Initialize hillslope geomorphology from input dataset
    !
    ! !USES:
    use LandunitType    , only : lun
    use GridcellType    , only : grc
    use ColumnType      , only : col
    use clm_varctl      , only : nhillslope, max_columns_hillslope
    use spmdMod         , only : masterproc
    use fileutils       , only : getfil
    use clm_varcon      , only : spval, ispval, grlnd
    use landunit_varcon , only : istsoil
    use subgridWeightsMod , only : compute_higher_order_weights
    use ncdio_pio
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    character(len=*) , intent(in) :: hillslope_file    ! hillslope data file name
    integer,  pointer     :: ihillslope_in(:,:) ! read in - integer
    integer,  pointer     :: ncolumns_hillslope_in(:) ! read in number of columns
    integer,  allocatable :: ncolumns_hillslope(:)    ! number of hillslope columns
    integer,  allocatable :: hill_ndx(:,:)      ! hillslope index
    integer,  allocatable :: col_ndx(:,:)       ! column index
    integer,  allocatable :: col_dndx(:,:)      ! downhill column index
    integer,  allocatable :: hill_pftndx(:,:)   ! hillslope pft index []
    integer,  allocatable :: col_pftndx(:)      ! hillslope column pft index []
    real(r8), pointer     :: fhillslope_in(:,:) ! read in - float
    real(r8), allocatable :: pct_hillslope(:,:) ! percent of landunit occupied by hillslope
    real(r8), allocatable :: hill_slope(:,:)    ! hillslope slope  [m/m]
    real(r8), allocatable :: hill_aspect(:,:)   ! hillslope azimuth [radians]
    real(r8), allocatable :: hill_area(:,:)     ! hillslope area   [m2]
    real(r8), allocatable :: hill_dist(:,:)     ! hillslope length [m]
    real(r8), allocatable :: hill_width(:,:)    ! hillslope width  [m]
    real(r8), allocatable :: hill_elev(:,:)     ! hillslope height [m]
    real(r8), allocatable :: hill_bedrock(:,:)  ! hillslope bedrock depth [m]
    real(r8), pointer     :: fstream_in(:)      ! read in - 1D - float

    type(file_desc_t)     :: ncid                 ! netcdf id
    logical               :: readvar              ! check whether variable on file
    character(len=256)    :: locfn                ! local filename
    integer               :: ierr                 ! error code
    integer               :: c, l, g, i, j, ci, nh       ! indices

    real(r8)              :: ncol_per_hillslope(nhillslope) ! number of columns per hillslope
    real(r8)              :: hillslope_area(nhillslope)     ! area of hillslope
    real(r8)              :: nhill_per_landunit(nhillslope) ! total number of each representative hillslope per landunit

    character(len=*), parameter :: subname = 'InitHillslope'

    !-----------------------------------------------------------------------

    ! consistency check
    call check_aquifer_layer()

    ! Open hillslope dataset to read in data below

    call getfil (hillslope_file, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    allocate( &
         ncolumns_hillslope(bounds%begl:bounds%endl),  &
         pct_hillslope(bounds%begl:bounds%endl,nhillslope),  &
         hill_ndx     (bounds%begl:bounds%endl,max_columns_hillslope), &
         col_ndx      (bounds%begl:bounds%endl,max_columns_hillslope), &
         col_dndx     (bounds%begl:bounds%endl,max_columns_hillslope), &
         hill_slope   (bounds%begl:bounds%endl,max_columns_hillslope), &
         hill_aspect  (bounds%begl:bounds%endl,max_columns_hillslope), &
         hill_area    (bounds%begl:bounds%endl,max_columns_hillslope), &
         hill_dist    (bounds%begl:bounds%endl,max_columns_hillslope), &
         hill_width   (bounds%begl:bounds%endl,max_columns_hillslope), &
         hill_elev    (bounds%begl:bounds%endl,max_columns_hillslope), &
         col_pftndx   (bounds%begc:bounds%endc), &
         stat=ierr)

    allocate(ncolumns_hillslope_in(bounds%begg:bounds%endg))

    call ncd_io(ncid=ncid, varname='nhillcolumns', flag='read', data=ncolumns_hillslope_in, dim1name=grlnd, readvar=readvar)
    if (masterproc .and. .not. readvar) then
       call endrun( 'ERROR:: nhillcolumns not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       ncolumns_hillslope(l) = ncolumns_hillslope_in(g)
       ! vegetated landunits having nonzero hillslope columns and nonzero weight
       if (lun%wtgcell(l) > 0._r8 .and. lun%itype(l) == istsoil .and. ncolumns_hillslope_in(g) > 0) then
          do c = lun%coli(l), lun%colf(l)
             col%is_hillslope_column(c) = .true.
          enddo
       end if
    enddo
    deallocate(ncolumns_hillslope_in)

    allocate(fhillslope_in(bounds%begg:bounds%endg,nhillslope))

    call ncd_io(ncid=ncid, varname='pct_hillslope', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (masterproc .and. .not. readvar) then
       call endrun( 'ERROR:: pct_hillslope not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       pct_hillslope(l,:) = fhillslope_in(g,:)
    enddo
    deallocate(fhillslope_in)

    allocate(ihillslope_in(bounds%begg:bounds%endg,max_columns_hillslope))

    call ncd_io(ncid=ncid, varname='hillslope_index', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
    if (masterproc .and. .not. readvar) then
       call endrun( 'ERROR:: hillslope_index not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_ndx(l,:) = ihillslope_in(g,:)
    enddo

    call ncd_io(ncid=ncid, varname='column_index', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
    if (masterproc .and. .not. readvar) then
       call endrun( 'ERROR:: column_index not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       col_ndx(l,:) = ihillslope_in(g,:)
    enddo

    call ncd_io(ncid=ncid, varname='downhill_column_index', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
    if (masterproc .and. .not. readvar) then
       call endrun( 'ERROR:: downhill_column_index not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       col_dndx(l,:) = ihillslope_in(g,:)
    enddo
    deallocate(ihillslope_in)

    allocate(fhillslope_in(bounds%begg:bounds%endg,max_columns_hillslope))
    call ncd_io(ncid=ncid, varname='hillslope_slope', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (masterproc .and. .not. readvar) then
       call endrun( 'ERROR:: hillslope_slope not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
    end if

    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_slope(l,:) = fhillslope_in(g,:)
    enddo

    call ncd_io(ncid=ncid, varname='hillslope_aspect', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (masterproc .and. .not. readvar) then
       call endrun( 'ERROR:: hillslope_aspect not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
    end if

    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_aspect(l,:) = fhillslope_in(g,:)
    enddo

    call ncd_io(ncid=ncid, varname='hillslope_area', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (masterproc .and. .not. readvar) then
       call endrun( 'ERROR:: hillslope_area not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_area(l,:) = fhillslope_in(g,:)
    enddo
    call ncd_io(ncid=ncid, varname='hillslope_distance', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (masterproc .and. .not. readvar) then
       call endrun( 'ERROR:: hillslope_distance not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
    end if

    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_dist(l,:) = fhillslope_in(g,:)
    enddo

    call ncd_io(ncid=ncid, varname='hillslope_width', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (masterproc .and. .not. readvar) then
       call endrun( 'ERROR:: hillslope_width not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_width(l,:) = fhillslope_in(g,:)
    enddo

    call ncd_io(ncid=ncid, varname='hillslope_elevation', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (masterproc .and. .not. readvar) then
       call endrun( 'ERROR:: hillslope_elevation not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_elev(l,:) = fhillslope_in(g,:)
    enddo

    deallocate(fhillslope_in)

    allocate(ihillslope_in(bounds%begg:bounds%endg,max_columns_hillslope))
    call ncd_io(ncid=ncid, varname='hillslope_pftndx', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
    if (readvar) then
       allocate(hill_pftndx (bounds%begl:bounds%endl,max_columns_hillslope), stat=ierr)
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          hill_pftndx(l,:) = ihillslope_in(g,:)
       enddo
    end if

    deallocate(ihillslope_in)

    if (use_hillslope_routing) then
       allocate(fstream_in(bounds%begg:bounds%endg))

       call ncd_io(ncid=ncid, varname='hillslope_stream_depth', flag='read', data=fstream_in, dim1name=grlnd, readvar=readvar)
       if (masterproc .and. .not. readvar) then
          call endrun( 'ERROR:: hillslope_stream_depth not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          lun%stream_channel_depth(l) = fstream_in(g)
       enddo

       call ncd_io(ncid=ncid, varname='hillslope_stream_width', flag='read', data=fstream_in, dim1name=grlnd, readvar=readvar)
       if (masterproc .and. .not. readvar) then
          call endrun( 'ERROR:: hillslope_stream_width not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          lun%stream_channel_width(l) = fstream_in(g)
       enddo

       call ncd_io(ncid=ncid, varname='hillslope_stream_slope', flag='read', data=fstream_in, dim1name=grlnd, readvar=readvar)
       if (masterproc .and. .not. readvar) then
          call endrun( 'ERROR:: hillslope_stream_slope not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          lun%stream_channel_slope(l) = fstream_in(g)
       enddo

       deallocate(fstream_in)
    end if

    !  Set hillslope hydrology column level variables
    !  This needs to match how columns set up in subgridMod
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       if (lun%itype(l) == istsoil) then

          ! map external column index to internal column index
          do c = lun%coli(l), lun%colf(l)
             ! ci should span [1:nhillcolumns(l)]
             ci = c-lun%coli(l)+1

             if (col_dndx(l,ci) <= -999) then
                ! lowermost column of hillslope has no downstream neighbor
                col%cold(c) = ispval
             else
                ! relative separation should be the same
                col%cold(c) = c + (col_dndx(l,ci) - col_ndx(l,ci))
             end if
          enddo

          do c = lun%coli(l), lun%colf(l)
             ci = c-lun%coli(l)+1

             col%hillslope_ndx(c) = hill_ndx(l,ci)

             ! Find uphill neighbors (this may not actually be useful...)
             col%colu(c) = ispval
             do i = lun%coli(l), lun%colf(l)
                if (c == col%cold(i)) then
                   col%colu(c) = i
                end if
             enddo

             ! distance of lower edge of column from hillslope bottom
             col%hill_distance(c) = hill_dist(l,ci)
             ! width of lower edge of column
             col%hill_width(c)    = hill_width(l,ci)
             ! mean elevation of column relative to gridcell mean elevation
             col%hill_elev(c)     = hill_elev(l,ci)
             ! mean along-hill slope of column
             col%hill_slope(c)    = hill_slope(l,ci)
             ! area of column
             col%hill_area(c)     = hill_area(l,ci)
             ! azimuth of column
             col%hill_aspect(c)   = hill_aspect(l,ci)
             ! pft index of column
             if ( allocated(hill_pftndx) ) then
                col_pftndx(c) = hill_pftndx(l,ci)
             end if

          enddo

          ! Calculate total hillslope area on landunit and
          ! number of columns in each hillslope
          ncol_per_hillslope(:)= 0._r8
          hillslope_area(:)    = 0._r8
          do c = lun%coli(l), lun%colf(l)
             nh = col%hillslope_ndx(c)
             if (nh > 0) then
                ncol_per_hillslope(nh) = ncol_per_hillslope(nh) + 1
                hillslope_area(nh) = hillslope_area(nh) + col%hill_area(c)
             end if
          enddo

          if (use_hillslope_routing) then

             ! Total area occupied by each hillslope (m2) is
             ! grc%area(g)*1.e6*lun%wtgcell(l)*pct_hillslope(l,nh)*0.01
             ! Number of representative hillslopes per landunit
             ! is the total area divided by individual area
             ! include factor of 0.5 because a channel is shared by ~2 hillslopes

             lun%stream_channel_number(l) = 0._r8
             do nh = 1, nhillslope
                if (hillslope_area(nh) > 0._r8) then
                   nhill_per_landunit(nh) = grc%area(g)*1.e6_r8*lun%wtgcell(l) &
                        *pct_hillslope(l,nh)*0.01/hillslope_area(nh)

                   lun%stream_channel_number(l) = lun%stream_channel_number(l) &
                        + 0.5_r8 * nhill_per_landunit(nh)
                end if
             enddo

             ! Calculate steam channel length
             ! Total length of stream banks is individual widths
             ! times number of hillslopes per landunit
             ! include factor of 0.5 because a channel is shared by ~2 hillslopes
             lun%stream_channel_length(l) = 0._r8
             do c = lun%coli(l), lun%colf(l)
                if (col%cold(c) == ispval) then
                   lun%stream_channel_length(l) = lun%stream_channel_length(l) &
                        + col%hill_width(c) * 0.5_r8 * nhill_per_landunit(col%hillslope_ndx(c))
                end if
             enddo
          end if

          ! if missing hillslope information on dataset,
          ! call endrun
          if (ncolumns_hillslope(l) > 0 .and. sum(hillslope_area) == 0._r8 .and. masterproc) then
             write(iulog,*) 'Problem with input data: nhillcolumns is non-zero, but hillslope area is zero'
             write(iulog,*) 'Check hillslope data for gridcell at (lon/lat): ', grc%londeg(g),grc%latdeg(g)
             call endrun( 'ERROR:: sum of hillslope areas is zero.'//errmsg(sourcefile, __LINE__) )
          end if

          ! Recalculate column weights using input areas
          ! The higher order weights will be updated in a subsequent reweight_wrapup call
          do c = lun%coli(l), lun%colf(l)
             nh = col%hillslope_ndx(c)
             if (col%is_hillslope_column(c)) then
                col%wtlunit(c) = (col%hill_area(c)/hillslope_area(nh)) &
                     * (pct_hillslope(l,nh)*0.01_r8)
             end if
          enddo
       end if
    enddo ! end of landunit loop

    deallocate(ncolumns_hillslope,pct_hillslope,hill_ndx,col_ndx,col_dndx, &
         hill_slope,hill_area,hill_dist, &
         hill_width,hill_elev,hill_aspect)

    ! Modify pft distributions
    ! this may require modifying subgridMod/natveg_patch_exists
    ! to ensure patch exists in every gridcell
    if (pft_distribution_method == pft_from_file) then
       call HillslopePftFromFile(bounds,col_pftndx)
    else if (pft_distribution_method == pft_lowland_dominant_pft) then
       ! Specify different pfts for uplands / lowlands
       call HillslopeDominantLowlandPft(bounds)
    else if (pft_distribution_method == pft_lowland_upland) then
       ! example usage:
       ! upland_ivt  = 13 ! c3 non-arctic grass
       ! lowland_ivt = 7  ! broadleaf deciduous tree
       call HillslopeSetLowlandUplandPfts(bounds,lowland_ivt=7,upland_ivt=13)
    else if (masterproc .and. .not. (pft_distribution_method == pft_standard .or. pft_distribution_method ==pft_uniform_dominant_pft)) then
      call endrun( 'ERROR:: unrecognized hillslope_pft_distribution_method'//errmsg(sourcefile, __LINE__) )
    end if

    if ( allocated(hill_pftndx) ) then
       deallocate(hill_pftndx)
       deallocate(col_pftndx)
    end if

    ! Update higher order weights and check that weights sum to 1
    call compute_higher_order_weights(bounds)

    call ncd_pio_closefile(ncid)

  end subroutine InitHillslope

  !-----------------------------------------------------------------------

  subroutine SetHillslopeSoilThickness(bounds, hillslope_file, soil_depth_lowland_in, soil_depth_upland_in)
    !
    ! !DESCRIPTION:
    ! Set hillslope column nbedrock values
    !
    ! !USES:
    use LandunitType    , only : lun
    use GridcellType    , only : grc
    use ColumnType      , only : col
    use clm_varctl      , only : nhillslope, max_columns_hillslope
    use clm_varcon      , only : zmin_bedrock, zisoi
    use clm_varpar      , only : nlevsoi
    use spmdMod         , only : masterproc
    use fileutils       , only : getfil
    use clm_varcon      , only : spval, ispval, grlnd
    use ncdio_pio
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    character(len=*) , intent(in) :: hillslope_file    ! hillslope data file name
    real(r8), intent(in), optional  :: soil_depth_lowland_in
    real(r8), intent(in), optional  :: soil_depth_upland_in
    real(r8), pointer     :: fhillslope_in(:,:) ! read in - float

    type(file_desc_t)     :: ncid                 ! netcdf id
    logical               :: readvar              ! check whether variable on file
    character(len=256)    :: locfn                ! local filename
    integer               :: ierr                 ! error code
    integer               :: c, l, g, j, ci       ! indices

    real(r8)              :: soil_depth_lowland
    real(r8)              :: soil_depth_upland
    real(r8), parameter   :: soil_depth_lowland_default = 8.0
    real(r8), parameter   :: soil_depth_upland_default  = 8.0
    character(len=*), parameter :: subname = 'SetHillslopeSoilThickness'

    !-----------------------------------------------------------------------

    if (soil_profile_method==soil_profile_from_file) then

       ! Open hillslope dataset to read in data below
       call getfil (hillslope_file, locfn, 0)
       call ncd_pio_openfile (ncid, locfn, 0)

       allocate(fhillslope_in(bounds%begg:bounds%endg,max_columns_hillslope))
       call ncd_io(ncid=ncid, varname='hillslope_bedrock_depth', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
       if (masterproc .and. .not. readvar) then
          call endrun( 'ERROR:: soil_profile_method = "FromFile", but hillslope_bedrock not found on hillslope data set.'//errmsg(sourcefile, __LINE__) )
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          do c = lun%coli(l), lun%colf(l)
             if (col%is_hillslope_column(c) .and. col%active(c)) then
                ci = c-lun%coli(l)+1
                do j = 1,nlevsoi
                   if (zisoi(j-1) > zmin_bedrock) then
                      if (zisoi(j-1) < fhillslope_in(g,ci) &
                           .and. zisoi(j) >= fhillslope_in(g,ci)) then
                         col%nbedrock(c) = j
                      end if
                   end if
                enddo
             end if
          enddo
       enddo
       deallocate(fhillslope_in)
       call ncd_pio_closefile(ncid)

    else if (soil_profile_method==soil_profile_set_lowland_upland &
         .or. soil_profile_method==soil_profile_linear) then

       if (present(soil_depth_lowland_in)) then
          soil_depth_lowland = soil_depth_lowland_in
       else
          soil_depth_lowland = soil_depth_lowland_default
       end if

       if (present(soil_depth_upland_in)) then
          soil_depth_upland = soil_depth_upland_in
       else
          soil_depth_upland = soil_depth_upland_default
       end if

       ! Modify hillslope soil thickness profile
       call HillslopeSoilThicknessProfile(bounds,&
            soil_profile_method=soil_profile_method,&
            soil_depth_lowland_in=soil_depth_lowland,&
            soil_depth_upland_in=soil_depth_upland)

    else if (soil_profile_method /= soil_profile_uniform .and. masterproc) then
       call endrun( msg=' ERROR: unrecognized hillslope_soil_profile_method'//errMsg(sourcefile, __LINE__))

    end if

  end subroutine SetHillslopeSoilThickness

  !-----------------------------------------------------------------------
  subroutine HillslopeSoilThicknessProfile(bounds,&
       soil_profile_method,soil_depth_lowland_in,soil_depth_upland_in)
    !
    ! !DESCRIPTION:
    ! Modify soil thickness across hillslope by changing
    ! col%nbedrock
    !
    ! !USES:
    use LandunitType    , only : lun
    use GridcellType    , only : grc
    use ColumnType      , only : col
    use clm_varcon      , only : zmin_bedrock, zisoi
    use clm_varpar      , only : nlevsoi
    use spmdMod         , only : masterproc
    use fileutils       , only : getfil
    use clm_varcon      , only : spval, ispval, grlnd
    use ncdio_pio
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer,  intent(in)  :: soil_profile_method
    real(r8), intent(in), optional  :: soil_depth_lowland_in
    real(r8), intent(in), optional  :: soil_depth_upland_in

    integer               :: c, l, g, i, j
    real(r8)              :: min_hill_dist, max_hill_dist
    real(r8)              :: m, b           ! linear soil thickness slope/intercept
    real(r8)              :: soil_depth_col
    real(r8)              :: soil_depth_lowland
    real(r8)              :: soil_depth_upland
    real(r8), parameter   :: soil_depth_lowland_default = 8.0
    real(r8), parameter   :: soil_depth_upland_default  = 8.0

    character(len=*), parameter :: subname = 'HillslopeSoilThicknessProfile'

    !-----------------------------------------------------------------------

    if (present(soil_depth_lowland_in)) then
       soil_depth_lowland = soil_depth_lowland_in
    else
       soil_depth_lowland = soil_depth_lowland_default
    end if

    if (present(soil_depth_upland_in)) then
       soil_depth_upland = soil_depth_upland_in
    else
       soil_depth_upland = soil_depth_upland_default
    end if

    ! Specify lowland/upland soil thicknesses separately
    if (soil_profile_method == soil_profile_set_lowland_upland) then
       do c = bounds%begc,bounds%endc
          if (col%is_hillslope_column(c) .and. col%active(c)) then
             if (col%cold(c) /= ispval) then
                do j = 1,nlevsoi
                   if (zisoi(j-1) > zmin_bedrock) then
                      if (zisoi(j-1) < soil_depth_upland .and. zisoi(j) >= soil_depth_upland) then
                         col%nbedrock(c) = j
                      end if
                   end if
                enddo
             else
                do j = 1,nlevsoi
                   if (zisoi(j-1) > zmin_bedrock) then
                      if (zisoi(j-1) < soil_depth_lowland .and. zisoi(j) >= soil_depth_lowland) then
                         col%nbedrock(c) = j
                      end if
                   end if
                enddo
             end if
          end if
       end do
    ! Linear soil thickness profile
    else if (soil_profile_method == soil_profile_linear) then
       call HillslopeSoilThicknessProfile_linear(bounds, soil_depth_lowland, soil_depth_upland)
    else if (masterproc) then
       call endrun( 'ERROR:: invalid soil_profile_method.'//errmsg(sourcefile, __LINE__) )
    end if

  end subroutine HillslopeSoilThicknessProfile

  !------------------------------------------------------------------------
  subroutine HillslopeSetLowlandUplandPfts(bounds,lowland_ivt,upland_ivt)
    !
    ! !DESCRIPTION:
    ! Reassign patch type of each column based on whether a column
    ! is identified as a lowland or an upland.
    ! Assumes each column has a single pft.
    ! In preparation for this reassignment of patch type, only the
    ! first patch was given a non-zero weight in surfrd_hillslope
    !
    ! !USES
    use LandunitType    , only : lun
    use ColumnType      , only : col
    use clm_varcon      , only : ispval
    use clm_varpar      , only : natpft_lb
    use PatchType       , only : patch
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: upland_ivt
    integer, intent(in) :: lowland_ivt
    !
    ! !LOCAL VARIABLES:
    integer :: p,c    ! indices
    integer :: npatches_per_column

    !------------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       if (col%is_hillslope_column(c)) then
          npatches_per_column = 0
          do p = col%patchi(c), col%patchf(c)
             if (col%cold(c) == ispval) then
                ! lowland
                patch%itype(p) = lowland_ivt
             else
                ! upland
                patch%itype(p) = upland_ivt
             end if
             ! update mxy as is done in initSubgridMod.add_patch
             patch%mxy(p) = patch%itype(p) + (1 - natpft_lb)

             npatches_per_column = npatches_per_column + 1
          enddo
          if ((npatches_per_column /= 1) .and. masterproc) then
             call endrun( 'ERROR:: number of patches per hillslope column not equal to 1'//errmsg(sourcefile, __LINE__) )
          end if
       end if
    enddo

  end subroutine HillslopeSetLowlandUplandPfts

  !------------------------------------------------------------------------
  subroutine HillslopeDominantLowlandPft(bounds)
    !
    ! !DESCRIPTION:
    ! Reassign patch weights of each column based on each gridcell's
    ! two most dominant pfts on the input dataset.
    ! HillslopeTwoLargestPftIndices is called in surfrd_hillslope to
    ! prepare the patch weights for this routine.
    ! Assumes each column has a single pft.
    ! Use largest weight for lowland, 2nd largest weight for uplands
    !
    ! !USES
    use LandunitType    , only : lun
    use ColumnType      , only : col
    use decompMod       , only : get_clump_bounds, get_proc_clumps
    use clm_varcon      , only : ispval
    use PatchType       , only : patch
    use pftconMod       , only : pftcon, ndllf_evr_tmp_tree, nc3_nonarctic_grass, nc4_grass
    use array_utils     , only : find_k_max_indices
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,c    ! indices
    integer :: plow, phigh
    integer :: max_index(1)
    integer, allocatable  :: max_indices(:)    ! largest weight pft indices
    real(r8) :: sum_wtcol, sum_wtlun, sum_wtgrc

    !------------------------------------------------------------------------

    allocate(max_indices(2))
    do c = bounds%begc,bounds%endc
       if (col%is_hillslope_column(c)) then

          ! if only one pft exists, find dominant pft index and set 2nd index to the same value

          if (size(patch%wtcol(col%patchi(c):col%patchf(c))) == 1) then
             call find_k_max_indices(patch%wtcol(col%patchi(c):col%patchf(c)),1,1,max_index)
             max_indices(1) = max_index(1) + (col%patchi(c) - 1)
             max_indices(2) = max_indices(1)
          else
             call find_k_max_indices(patch%wtcol(col%patchi(c):col%patchf(c)),1,2,max_indices)
             max_indices = max_indices + (col%patchi(c) - 1)
          end if

          sum_wtcol = sum(patch%wtcol(col%patchi(c):col%patchf(c)))
          sum_wtlun = sum(patch%wtlunit(col%patchi(c):col%patchf(c)))
          sum_wtgrc = sum(patch%wtgcell(col%patchi(c):col%patchf(c)))

          patch%wtcol(col%patchi(c):col%patchf(c)) = 0._r8
          patch%wtlunit(col%patchi(c):col%patchf(c)) = 0._r8
          patch%wtgcell(col%patchi(c):col%patchf(c)) = 0._r8

          ! Put the highest stature vegetation on the lowland column
          ! non-tree and tree       ; place tree on lowland
          ! grass and shrub         ; place shrub on lowland
          ! bare soil and vegetation; place vegetation on lowland
          if ((.not. pftcon%is_tree(patch%itype(max_indices(1))) .and. pftcon%is_tree(patch%itype(max_indices(2)))) &
               .or. (pftcon%is_grass(patch%itype(max_indices(1))) .and. pftcon%is_shrub(patch%itype(max_indices(2)))) &
               .or. (patch%itype(max_indices(1)) == 0)) then
             plow = max_indices(2)
             phigh = max_indices(1)
          else
             plow = max_indices(1)
             phigh = max_indices(2)
          end if

          ! Special cases (subjective)

          ! if NET/BDT assign BDT to lowland
          if ((patch%itype(max_indices(1)) == ndllf_evr_tmp_tree) .and. pftcon%is_tree(patch%itype(max_indices(2)))) then
             plow = max_indices(2)
             phigh = max_indices(1)
          end if
          ! if C3/C4 assign C4 to lowland
          if ((patch%itype(max_indices(1)) == nc4_grass) .and. (patch%itype(max_indices(2)) == nc3_nonarctic_grass)) then
             plow = max_indices(1)
             phigh = max_indices(2)
          end if
          if ((patch%itype(max_indices(1)) == nc3_nonarctic_grass) .and. (patch%itype(max_indices(2)) == nc4_grass)) then
             plow = max_indices(2)
             phigh = max_indices(1)
          end if

          if (col%cold(c) == ispval) then
             ! lowland column
             patch%wtcol(plow)   = sum_wtcol
             patch%wtlunit(plow) = sum_wtlun
             patch%wtgcell(plow) = sum_wtgrc
          else
             ! upland columns
             patch%wtcol(phigh)   = sum_wtcol
             patch%wtlunit(phigh) = sum_wtlun
             patch%wtgcell(phigh) = sum_wtgrc
          end if
       end if
    enddo    ! end loop c
    deallocate(max_indices)

  end subroutine HillslopeDominantLowlandPft

  !------------------------------------------------------------------------
  subroutine HillslopePftFromFile(bounds,col_pftndx)
    !
    ! !DESCRIPTION:
    ! Reassign patch type using indices from data file
    ! Assumes one patch per hillslope column
    ! In preparation for this reassignment of patch type, only the
    ! first patch was given a non-zero weight in surfrd_hillslope.
    !
    ! !USES
    use ColumnType      , only : col
    use PatchType       , only : patch
    use clm_varpar      , only : natpft_lb
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in)           :: col_pftndx(:)
    !
    ! !LOCAL VARIABLES:
    integer :: p,c    ! indices
    integer :: npatches_per_column

    !------------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       if (col%is_hillslope_column(c)) then
          ! In preparation for this re-weighting of patch type
          ! only first patch was given a non-zero weight in surfrd_hillslope
          npatches_per_column = 0
          do p = col%patchi(c), col%patchf(c)
             patch%itype(p) = col_pftndx(c)
             ! update mxy as is done in initSubgridMod.add_patch
             patch%mxy(p) = patch%itype(p) + (1 - natpft_lb)
             npatches_per_column = npatches_per_column + 1
          enddo
          if ((npatches_per_column /= 1) .and. masterproc) then
             call endrun( 'ERROR:: number of patches per hillslope column not equal to 1'//errmsg(sourcefile, __LINE__) )
          end if
       end if
    enddo

  end subroutine HillslopePftFromFile

  !-----------------------------------------------------------------------
  subroutine HillslopeStreamOutflow(bounds, &
       waterstatebulk_inst, waterfluxbulk_inst,streamflow_method)
    !
    ! !DESCRIPTION:
    ! Calculate discharge from stream channel
    !
    ! !USES:
    use LandunitType    , only : lun
    use GridcellType    , only : grc
    use ColumnType      , only : col
    use WaterFluxBulkType   , only : waterfluxbulk_type
    use WaterStateBulkType  , only : waterstatebulk_type
    use spmdMod         , only : masterproc
    use clm_varcon      , only : spval, ispval, grlnd
    use landunit_varcon , only : istsoil
    use ncdio_pio
    use clm_time_manager , only : get_step_size_real
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer,  intent(in)  :: streamflow_method
    type(waterstatebulk_type), intent(inout) :: waterstatebulk_inst
    type(waterfluxbulk_type),  intent(inout) :: waterfluxbulk_inst

    integer               :: c, l, g, i, j
    integer               :: nstep
    real(r8) :: dtime                                     ! land model time step (sec)
    real(r8)              :: cross_sectional_area         ! cross sectional area of stream water (m2)
    real(r8)              :: stream_depth                 ! depth of stream water (m)
    real(r8)              :: hydraulic_radius             ! cross sectional area divided by wetted perimeter (m)
    real(r8)              :: flow_velocity                ! flow velocity (m/s)
    real(r8)              :: overbank_area                ! area of water above bankfull (m2)
    real(r8), parameter   :: manning_roughness = 0.03_r8  ! manning roughness
    real(r8), parameter   :: manning_exponent  = 0.667_r8 ! manning exponent

    integer, parameter    :: overbank_method = 1          ! method to treat overbank stream storage; 1 = increase dynamic slope, 2 = increase flow area cross section, 3 = remove instantaneously
    logical               :: active_stream
    character(len=*), parameter :: subname = 'HillslopeStreamOutflow'

    !-----------------------------------------------------------------------
    associate( &
         stream_water_volume     =>    waterstatebulk_inst%stream_water_volume_lun , &  ! Input:  [real(r8) (:)   ] stream water volume (m3)
         volumetric_streamflow   =>    waterfluxbulk_inst%volumetric_streamflow_lun  &  ! Input:  [real(r8) (:)   ] stream water discharge (m3/s)
         )

      ! Get time step
      dtime = get_step_size_real()

      do l = bounds%begl,bounds%endl
         volumetric_streamflow(l) = 0._r8

         ! Check for vegetated landunits having initialized stream channel properties
         active_stream = .false.
         if (lun%itype(l) == istsoil .and. &
              lun%stream_channel_length(l) > 0._r8 .and. &
              lun%stream_channel_width(l) > 0._r8) then
            active_stream = .true.
         end if

         if (lun%active(l) .and. active_stream) then
            ! Streamflow calculated from Manning equation
            if (streamflow_method == streamflow_manning) then
               cross_sectional_area = stream_water_volume(l) &
                    /lun%stream_channel_length(l)
               stream_depth =  cross_sectional_area &
                    /lun%stream_channel_width(l)
               hydraulic_radius = cross_sectional_area &
                    /(lun%stream_channel_width(l) + 2*stream_depth)

               if (hydraulic_radius <= 0._r8) then
                  volumetric_streamflow(l) = 0._r8
               else
                  flow_velocity = (hydraulic_radius)**manning_exponent &
                       * sqrt(lun%stream_channel_slope(l)) &
                       / manning_roughness
                  ! overbank flow
                  if (stream_depth > lun%stream_channel_depth(l)) then
                     if (overbank_method  == 1) then
                        ! try increasing dynamic slope
                        volumetric_streamflow(l) = cross_sectional_area * flow_velocity &
                             *(stream_depth/lun%stream_channel_depth(l))
                     else if (overbank_method  == 2) then
                        ! try increasing flow area cross section
                        overbank_area = (stream_depth -lun%stream_channel_depth(l)) * 30._r8 * lun%stream_channel_width(l)
                        volumetric_streamflow(l) = (cross_sectional_area + overbank_area) * flow_velocity
                     else if (overbank_method  == 3) then
                        ! try removing all overbank flow instantly
                        volumetric_streamflow(l) = cross_sectional_area * flow_velocity &
                             + (stream_depth-lun%stream_channel_depth(l)) &
                             *lun%stream_channel_width(l)*lun%stream_channel_length(l)/dtime
                     else
                        call endrun( 'ERROR:: invalid overbank_method.'//errmsg(sourcefile, __LINE__) )
                     end if

                  else
                     volumetric_streamflow(l) = cross_sectional_area * flow_velocity
                  end if

                  ! scale streamflow by number of channel reaches
                  volumetric_streamflow(l) = volumetric_streamflow(l) * lun%stream_channel_number(l)

                  volumetric_streamflow(l) = max(0._r8,min(volumetric_streamflow(l),stream_water_volume(l)/dtime))
               end if
            else
               call endrun( 'ERROR:: invalid streamflow_method'//errmsg(sourcefile, __LINE__) )
            end if
         end if ! end of istsoil
      enddo    ! end of loop over landunits

  end associate

  end subroutine HillslopeStreamOutflow

  !-----------------------------------------------------------------------
  subroutine HillslopeUpdateStreamWater(bounds, waterstatebulk_inst, &
       waterfluxbulk_inst,waterdiagnosticbulk_inst)
    !
    ! !DESCRIPTION:
    ! Calculate discharge from stream channel
    !
    ! !USES:
    use LandunitType    , only : lun
    use GridcellType    , only : grc
    use ColumnType      , only : col
    use WaterFluxBulkType   , only : waterfluxbulk_type
    use WaterStateBulkType  , only : waterstatebulk_type
    use WaterDiagnosticBulkType  , only : waterdiagnosticbulk_type
    use spmdMod         , only : masterproc
    use clm_varcon      , only : spval, ispval, grlnd
    use landunit_varcon , only : istsoil
    use clm_time_manager, only : get_step_size_real
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    type(waterstatebulk_type), intent(inout) :: waterstatebulk_inst
    type(waterfluxbulk_type),  intent(inout) :: waterfluxbulk_inst
    type(waterdiagnosticbulk_type), intent(inout) :: waterdiagnosticbulk_inst

    integer  :: c, l, g, i, j
    real(r8) :: qflx_surf_vol               ! volumetric surface runoff (m3/s)
    real(r8) :: qflx_drain_perched_vol      ! volumetric perched saturated drainage (m3/s)
    real(r8) :: qflx_drain_vol              ! volumetric saturated drainage (m3/s)
    real(r8) :: dtime                       ! land model time step (sec)
    logical  :: active_stream

    character(len=*), parameter :: subname = 'HillslopeUpdateStreamWater'

    !-----------------------------------------------------------------------
    associate( &
         stream_water_volume     =>    waterstatebulk_inst%stream_water_volume_lun, & ! Input/Output:  [real(r8) (:)   ] stream water volume (m3)
         volumetric_streamflow   =>    waterfluxbulk_inst%volumetric_streamflow_lun,& ! Input:  [real(r8) (:)   ] stream water discharge (m3/s)
         qflx_drain              =>    waterfluxbulk_inst%qflx_drain_col,           & ! Input:  [real(r8) (:)   ]  column level sub-surface runoff (mm H2O /s)
         qflx_drain_perched      =>    waterfluxbulk_inst%qflx_drain_perched_col,   & ! Input:  [real(r8) (:)   ]  column level sub-surface runoff (mm H2O /s)
         qflx_surf               =>    waterfluxbulk_inst%qflx_surf_col,            & ! Input: [real(r8) (:)   ]  total surface runoff (mm H2O /s)
         stream_water_depth      =>    waterdiagnosticbulk_inst%stream_water_depth_lun   & ! Output:  [real(r8) (:)   ] stream water depth (m)
         )

       ! Get time step
       dtime = get_step_size_real()

       do l = bounds%begl,bounds%endl

          ! Check for vegetated landunits having initialized stream channel properties
          active_stream = .false.
          if (lun%itype(l) == istsoil .and. &
               lun%stream_channel_length(l) > 0._r8 .and. &
               lun%stream_channel_width(l) > 0._r8) then
             active_stream = .true.
          end if

          if (lun%active(l) .and. active_stream) then
             g = lun%gridcell(l)
             ! the drainage terms are 'net' quantities, so summing over
             ! all columns in a hillslope is equivalent to the outflow
             ! from the lowland column
             do c = lun%coli(l), lun%colf(l)
                if (col%is_hillslope_column(c) .and. col%active(c)) then
                   qflx_surf_vol = qflx_surf(c)*1.e-3_r8 &
                        *(grc%area(g)*1.e6_r8*col%wtgcell(c))
                   qflx_drain_perched_vol = qflx_drain_perched(c)*1.e-3_r8 &
                        *(grc%area(g)*1.e6_r8*col%wtgcell(c))
                   qflx_drain_vol = qflx_drain(c)*1.e-3_r8 &
                        *(grc%area(g)*1.e6_r8*col%wtgcell(c))

                   stream_water_volume(l) = stream_water_volume(l) &
                        + (qflx_drain_perched_vol &
                         + qflx_drain_vol + qflx_surf_vol) * dtime
                end if
             enddo
             stream_water_volume(l) = stream_water_volume(l) &
                  - volumetric_streamflow(l) * dtime

             ! account for negative drainage (via searchforwater in soilhydrology)
             if (stream_water_volume(l) < 0._r8) then
                volumetric_streamflow(l) = volumetric_streamflow(l) + stream_water_volume(l)/dtime
                stream_water_volume(l) = 0._r8
             end if

             stream_water_depth(l) = stream_water_volume(l) &
                  /lun%stream_channel_length(l) &
                  /lun%stream_channel_width(l)

          end if
       enddo

    end associate

  end subroutine HillslopeUpdateStreamWater

end module HillslopeHydrologyMod
