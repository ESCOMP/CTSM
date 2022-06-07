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

  ! !PUBLIC TYPES:
  implicit none

  private   
  save

  ! !PUBLIC MEMBER FUNCTIONS:
  public hillslope_properties_init
  public InitHillslope
  public HillslopeSoilThicknessProfile
  public HillslopeSetLowlandUplandPfts
  public HillslopeDominantPft
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
    
!    pft_distribution_method = pft_standard
!    soil_profile_method     = soil_profile_uniform
    
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

    endif

    call shr_mpi_bcast(pft_distribution_method, mpicom)
    call shr_mpi_bcast(soil_profile_method, mpicom)

    if (masterproc) then

       write(iulog,*) ' '
       write(iulog,*) 'hillslope_properties settings:'
       write(iulog,*) '  hillslope_pft_distribution_method = ',hillslope_pft_distribution_method
       write(iulog,*) '  hillslope_soil_profile_method     = ',hillslope_soil_profile_method

    endif

  end subroutine hillslope_properties_init

  !-----------------------------------------------------------------------

  subroutine InitHillslope(bounds,fsurdat,glc_behavior)
    !
    ! !DESCRIPTION:
    ! Initialize hillslope geomorphology from input dataset
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
    use landunit_varcon , only : istsoil
    use ncdio_pio
    use initVerticalMod , only : setSoilLayerClass
    use reweightMod     , only : reweight_wrapup
    use subgridWeightsMod , only : compute_higher_order_weights
    use glcBehaviorMod    , only : glc_behavior_type
    use subgridWeightsMod , only : check_weights
    
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    character(len=*) , intent(in) :: fsurdat    ! surface data file name
    type(glc_behavior_type), intent(in) :: glc_behavior
    integer,  pointer     :: ihillslope_in(:,:) ! read in - integer
    integer,  pointer     :: ncolumns_hillslope_in(:) ! read in number of columns
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
    real(r8), allocatable :: hill_length(:,:)   ! hillslope length [m]
    real(r8), allocatable :: hill_width(:,:)    ! hillslope width  [m]
    real(r8), allocatable :: hill_height(:,:)   ! hillslope height [m]
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

    ! Open surface dataset to read in data below 

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)
    
    allocate( &
         pct_hillslope(bounds%begl:bounds%endl,nhillslope),  &
         hill_ndx     (bounds%begl:bounds%endl,max_columns_hillslope), &
         col_ndx      (bounds%begl:bounds%endl,max_columns_hillslope), &
         col_dndx     (bounds%begl:bounds%endl,max_columns_hillslope), &
         hill_slope   (bounds%begl:bounds%endl,max_columns_hillslope), &
         hill_aspect  (bounds%begl:bounds%endl,max_columns_hillslope), &
         hill_area    (bounds%begl:bounds%endl,max_columns_hillslope), &
         hill_length  (bounds%begl:bounds%endl,max_columns_hillslope), &
         hill_width   (bounds%begl:bounds%endl,max_columns_hillslope), &
         hill_height  (bounds%begl:bounds%endl,max_columns_hillslope), &
         col_pftndx   (bounds%begc:bounds%endc), &
         stat=ierr)

    allocate(ncolumns_hillslope_in(bounds%begg:bounds%endg))
    
    call ncd_io(ncid=ncid, varname='nhillcolumns', flag='read', data=ncolumns_hillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: nhillcolumns not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       ! vegetated landunits having nonzero hillslope columns
       if(lun%itype(l) == istsoil .and. ncolumns_hillslope_in(g) > 0) then
          do c = lun%coli(l), lun%colf(l)
             col%is_hillslope_column(c) = .true.
          enddo
       endif
    enddo
    deallocate(ncolumns_hillslope_in)

    allocate(fhillslope_in(bounds%begg:bounds%endg,nhillslope))
    
    call ncd_io(ncid=ncid, varname='pct_hillslope', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: pct_hillslope not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       pct_hillslope(l,:) = fhillslope_in(g,:)
    enddo
    deallocate(fhillslope_in)
    
    allocate(ihillslope_in(bounds%begg:bounds%endg,max_columns_hillslope))
    
    call ncd_io(ncid=ncid, varname='hillslope_index', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: hillslope_index not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_ndx(l,:) = ihillslope_in(g,:)
    enddo
    
    call ncd_io(ncid=ncid, varname='column_index', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: column_index not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       col_ndx(l,:) = ihillslope_in(g,:)
    enddo
    
    call ncd_io(ncid=ncid, varname='downhill_column_index', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: downhill_column_index not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       col_dndx(l,:) = ihillslope_in(g,:)
    enddo
    deallocate(ihillslope_in)
    
    allocate(fhillslope_in(bounds%begg:bounds%endg,max_columns_hillslope))
    call ncd_io(ncid=ncid, varname='h_slope', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: h_slope not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_slope(l,:) = fhillslope_in(g,:)
    enddo
    
    call ncd_io(ncid=ncid, varname='h_aspect', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: h_aspect not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_aspect(l,:) = fhillslope_in(g,:)
    enddo
    
    call ncd_io(ncid=ncid, varname='h_area', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: h_area not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_area(l,:) = fhillslope_in(g,:)
    enddo
    call ncd_io(ncid=ncid, varname='h_length', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: h_length not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_length(l,:) = fhillslope_in(g,:)
    enddo
    
    call ncd_io(ncid=ncid, varname='h_width', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: h_width not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_width(l,:) = fhillslope_in(g,:)
    enddo
    
    call ncd_io(ncid=ncid, varname='h_height', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       if (masterproc) then
          call endrun( 'ERROR:: h_height not found on surface data set.'//errmsg(sourcefile, __LINE__) )
       end if
    end if
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       hill_height(l,:) = fhillslope_in(g,:)
    enddo

    call ncd_io(ncid=ncid, varname='h_bedrock', flag='read', data=fhillslope_in, dim1name=grlnd, readvar=readvar)
    if (readvar) then
       allocate(hill_bedrock (bounds%begl:bounds%endl,max_columns_hillslope), stat=ierr)
       if (masterproc) then
          write(iulog,*) 'h_bedrock found on surface data set'
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          hill_bedrock(l,:) = fhillslope_in(g,:)
       enddo
    end if

    deallocate(fhillslope_in)
    
    allocate(ihillslope_in(bounds%begg:bounds%endg,max_columns_hillslope))
    call ncd_io(ncid=ncid, varname='h_pftndx', flag='read', data=ihillslope_in, dim1name=grlnd, readvar=readvar)
    if (readvar) then
       allocate(hill_pftndx (bounds%begl:bounds%endl,max_columns_hillslope), stat=ierr)
       if (masterproc) then
          write(iulog,*) 'h_pftndx found on surface data set'
       end if
       do l = bounds%begl,bounds%endl
          g = lun%gridcell(l)
          hill_pftndx(l,:) = ihillslope_in(g,:)
       enddo
    end if

    deallocate(ihillslope_in)

    if (use_hillslope_routing) then 
       allocate(fstream_in(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncid, varname='h_stream_depth', flag='read', data=fstream_in, dim1name=grlnd, readvar=readvar)
       if (readvar) then
          if (masterproc) then
             write(iulog,*) 'h_stream_depth found on surface data set'
          end if
          do l = bounds%begl,bounds%endl
             g = lun%gridcell(l)
             lun%stream_channel_depth(l) = fstream_in(g)
          enddo
       else
          if (masterproc) then
             call endrun( 'ERROR:: h_stream_depth not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       endif
       call ncd_io(ncid=ncid, varname='h_stream_width', flag='read', data=fstream_in, dim1name=grlnd, readvar=readvar)
       if (readvar) then
          if (masterproc) then
             write(iulog,*) 'h_stream_width found on surface data set'
          end if
          do l = bounds%begl,bounds%endl
             g = lun%gridcell(l)
             lun%stream_channel_width(l) = fstream_in(g)
          enddo
       else
          if (masterproc) then
             call endrun( 'ERROR:: h_stream_width not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if
       call ncd_io(ncid=ncid, varname='h_stream_slope', flag='read', data=fstream_in, dim1name=grlnd, readvar=readvar)
       if (readvar) then
          if (masterproc) then
             write(iulog,*) 'h_stream_slope found on surface data set'
          end if
          do l = bounds%begl,bounds%endl
             g = lun%gridcell(l)
             lun%stream_channel_slope(l) = fstream_in(g)
          enddo
       else
          if (masterproc) then
             call endrun( 'ERROR:: h_stream_slope not found on surface data set.'//errmsg(sourcefile, __LINE__) )
          end if
       end if

       deallocate(fstream_in)
    endif
       
    !  Set hillslope hydrology column level variables
    !  This needs to match how columns set up in subgridMod
    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       if(lun%itype(l) == istsoil) then

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
             endif
          enddo
          
          do c = lun%coli(l), lun%colf(l)
             ci = c-lun%coli(l)+1
             
             col%hillslope_ndx(c) = hill_ndx(l,ci)

             ! Find uphill neighbors (this may not actually be useful...)
             col%colu(c) = ispval
             do i = lun%coli(l), lun%colf(l)
                if(c == col%cold(i)) then
                   col%colu(c) = i
                endif
             enddo
             
             ! distance of lower edge of column from hillslope bottom
             col%hill_distance(c) = hill_length(l,ci)
             ! width of lower edge of column 
             col%hill_width(c) = hill_width(l,ci)
             ! mean elevation of column relative to gridcell mean elevation
             col%hill_elev(c) = hill_height(l,ci)
             ! mean along-hill slope of column
             col%hill_slope(c) = hill_slope(l,ci)
             ! area of column
             col%hill_area(c) = hill_area(l,ci)
             ! azimuth of column
             col%hill_aspect(c) = hill_aspect(l,ci)
             ! soil thickness of column
             if (soil_profile_method==soil_profile_from_file) then
                if ( allocated(hill_bedrock) ) then
                   do j = 1,nlevsoi
                      if(zisoi(j-1) > zmin_bedrock) then
                         if (zisoi(j-1) < hill_bedrock(l,ci) .and. zisoi(j) >= hill_bedrock(l,ci)) then
                            col%nbedrock(c) = j
                         end if
                      endif
                   enddo
                else
                   if (masterproc) then
                      call endrun( 'ERROR:: soil_profile_method = soil_profile_from_file, but h_bedrock not found on surface data set.'//errmsg(sourcefile, __LINE__) )
                   end if
                endif
             endif
             ! pft index of column
             if ( allocated(hill_pftndx) ) then
                col_pftndx(c) = hill_pftndx(l,ci)
             endif
             
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
             endif
          enddo

          if (use_hillslope_routing) then 
                       
             ! Total area occupied by each hillslope (m2) is
             ! grc%area(g)*1.e6*lun%wtgcell(l)*pct_hillslope(l,nh)*0.01
             ! Number of representative hillslopes per landunit
             ! is the total area divided by individual area

             lun%stream_channel_number(l) = 0._r8
             do nh = 1, nhillslope
                if(hillslope_area(nh) > 0._r8) then
                   nhill_per_landunit(nh) = grc%area(g)*1.e6_r8*lun%wtgcell(l) &
                        *pct_hillslope(l,nh)*0.01/hillslope_area(nh)

                   lun%stream_channel_number(l) = lun%stream_channel_number(l) &
                        + nhill_per_landunit(nh)
                endif
             enddo
             
             ! Calculate steam channel length
             ! Total length of stream banks is individual widths
             ! times number of hillslopes per landunit divided
             ! by 2 to convert from bank length to channel length

             lun%stream_channel_length(l) = 0._r8
             do c = lun%coli(l), lun%colf(l)
                if(col%cold(c) == ispval) then
                   lun%stream_channel_length(l) = lun%stream_channel_length(l) &
                        + col%hill_width(c) * nhill_per_landunit(col%hillslope_ndx(c))
                endif
             enddo
             lun%stream_channel_length(l) = 0.5_r8 * lun%stream_channel_length(l)
          endif
                       
          ! if missing hillslope information on surface dataset, fill data
          ! and recalculate hillslope_area
          if (sum(hillslope_area) == 0._r8) then
             do c = lun%coli(l), lun%colf(l)
                nh = col%hillslope_ndx(c)
                col%hill_area(c) = (grc%area(g)/real(lun%ncolumns(l),r8))*1.e6_r8 ! km2 to m2
                col%hill_width(c)    = sqrt(col%hill_area(c))
                col%hill_slope(c)    = tan((rpi/180.)*col%topo_slope(c))
                col%hill_aspect(c)   = (rpi/2.) ! east (arbitrarily chosen)
                if (nh > 0) then
                   col%hill_elev(c)     = col%topo_std(c) &
                        *((c-lun%coli(l))/ncol_per_hillslope(nh))
                   col%hill_distance(c) = sqrt(col%hill_area(c)) &
                        *((c-lun%coli(l))/ncol_per_hillslope(nh))
                   pct_hillslope(l,nh)  = 100/nhillslope
                else
                   col%hill_elev(c)     = col%topo_std(c)
                   col%hill_distance(c) = sqrt(col%hill_area(c))
                endif
             enddo
             
          endif
          
          ! Recalculate column weights using input areas
          do c = lun%coli(l), lun%colf(l)
             nh = col%hillslope_ndx(c)
             if (col%is_hillslope_column(c)) then
                col%wtlunit(c) = (col%hill_area(c)/hillslope_area(nh)) &
                     * (pct_hillslope(l,nh)*0.01_r8)
             endif
          enddo
       endif
    enddo ! end of landunit loop
    
    deallocate(pct_hillslope,hill_ndx,col_ndx,col_dndx, &
         hill_slope,hill_area,hill_length, &
         hill_width,hill_height,hill_aspect)

    if(soil_profile_method==soil_profile_from_file) then
       if ( allocated(hill_bedrock) ) then
          deallocate(hill_bedrock)
       endif
    else if (soil_profile_method==soil_profile_set_lowland_upland &
         .or. soil_profile_method==soil_profile_linear) then 
       ! Modify hillslope soil thickness profile
       call HillslopeSoilThicknessProfile(bounds,&
            soil_profile_method=soil_profile_method,&
            soil_depth_lowland_in=8.0_r8,soil_depth_upland_in=8.0_r8)
    endif

    ! Update layer classes if nbedrock has been modified
    call setSoilLayerClass(bounds)
       
    ! Modify pft distributions
    ! this may require modifying subgridMod/natveg_patch_exists
    ! to ensure patch exists in every gridcell
    if (pft_distribution_method == pft_from_file) then
       call HillslopePftFromFile(bounds,col_pftndx)
    else if (pft_distribution_method == pft_uniform_dominant_pft) then
       ! Specify single dominant pft per gridcell
       call HillslopeDominantPft(bounds)
    else if (pft_distribution_method == pft_lowland_dominant_pft) then
       ! Specify different pfts for uplands / lowlands
       call HillslopeDominantLowlandPft(bounds)
    else if (pft_distribution_method == pft_lowland_upland) then
       ! example usage: 
       ! upland_ivt  = 13 ! c3 non-arctic grass
       ! lowland_ivt = 7  ! broadleaf deciduous tree 
       call HillslopeSetLowlandUplandPfts(bounds,lowland_ivt=7,upland_ivt=13)
    endif
    
    if ( allocated(hill_pftndx) ) then
       deallocate(hill_pftndx)
       deallocate(col_pftndx)
    endif

    ! Update higher order weights and check that weights sum to 1
    call compute_higher_order_weights(bounds)
    if (masterproc) then
       write(iulog,*) 'Checking modified hillslope weights via reweight_wrapup'
    end if

    ! filters have not been allocated yet!
    ! can check_weights be called directly here?
!    call reweight_wrapup(bounds, glc_behavior)
!    call check_weights(bounds)
    call check_weights(bounds, active_only=.true.)
    
    call ncd_pio_closefile(ncid)
    
  end subroutine InitHillslope

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

    if(present(soil_depth_lowland_in)) then
       soil_depth_lowland = soil_depth_lowland_in
    else
       soil_depth_lowland = soil_depth_lowland_default
    endif
    
    if(present(soil_depth_upland_in)) then
       soil_depth_upland = soil_depth_upland_in
    else
       soil_depth_upland = soil_depth_upland_default
    endif

    ! Specify lowland/upland soil thicknesses separately
    if(soil_profile_method == soil_profile_set_lowland_upland) then
       do c = bounds%begc,bounds%endc
          if (col%is_hillslope_column(c) .and. col%active(c)) then
             if(col%cold(c) /= ispval) then 
                do j = 1,nlevsoi
                   if(zisoi(j-1) > zmin_bedrock) then
                      if (zisoi(j-1) < soil_depth_upland .and. zisoi(j) >= soil_depth_upland) then
                         col%nbedrock(c) = j
                      end if
                   end if
                enddo
             else 
                do j = 1,nlevsoi 
                   if(zisoi(j-1) > zmin_bedrock) then
                      if (zisoi(j-1) < soil_depth_lowland .and. zisoi(j) >= soil_depth_lowland) then
                         col%nbedrock(c) = j
                      end if
                   end if
                enddo
             endif
          endif
       end do
    ! Linear soil thickness profile
    else if(soil_profile_method == soil_profile_linear) then
       do l = bounds%begl,bounds%endl
          min_hill_dist = minval(col%hill_distance(lun%coli(l):lun%colf(l)))
          max_hill_dist = maxval(col%hill_distance(lun%coli(l):lun%colf(l)))
          m = (soil_depth_lowland - soil_depth_upland)/ &
               (max_hill_dist - min_hill_dist)
          b = soil_depth_upland
          
          do c =  lun%coli(l), lun%colf(l)
             if (col%is_hillslope_column(c) .and. col%active(c)) then
                soil_depth_col = m*(max_hill_dist - col%hill_distance(c)) + b

                do j = 1,nlevsoi
                   if ((zisoi(j-1) <  soil_depth_col) .and. (zisoi(j) >= soil_depth_col)) then
                      col%nbedrock(c) = j
                   end if
                enddo
             endif
          enddo
       enddo
    else
       if (masterproc) then
          call endrun( 'ERROR:: invalid soil_profile_method.'//errmsg(sourcefile, __LINE__) )
       end if
    endif
       
  end subroutine HillslopeSoilThicknessProfile

  !------------------------------------------------------------------------
  subroutine HillslopeSetLowlandUplandPfts(bounds,lowland_ivt,upland_ivt)
    !
    ! !DESCRIPTION: 
    ! Reassign patch weights such that each column has a single pft.

    !
    ! !USES
    use LandunitType    , only : lun                
    use ColumnType      , only : col                
    use clm_varcon      , only : ispval
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
    logical :: check_npatches = .true.

    !------------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       if (col%is_hillslope_column(c) .and. col%active(c)) then
          ! In preparation for this re-weighting of patch type
          ! only first patch was given a non-zero weight in surfrd_hillslope
          npatches_per_column = 0
          do p = col%patchi(c), col%patchf(c)
             ! lowland
             if(col%cold(c) == ispval) then
                patch%itype(p) = lowland_ivt
             else ! upland
                patch%itype(p) = upland_ivt
             endif
             npatches_per_column = npatches_per_column + 1
          enddo
          if (check_npatches) then
             if ((npatches_per_column /= 1) .and. masterproc) then
                call endrun( 'ERROR:: number of patches per hillslope column not equal to 1'//errmsg(sourcefile, __LINE__) )
             end if
          endif
       endif
    enddo

  end subroutine HillslopeSetLowlandUplandPfts

  !------------------------------------------------------------------------
  subroutine HillslopeDominantPft(bounds)
    !
    ! !DESCRIPTION: 
    ! Reassign patch weights such that each column has a single pft  
    ! determined by each column's most dominant pft on input dataset.  
    ! Best performance when used with n_dom_pfts = 1 (Actually, this
    ! is probably redundant to behavior with n_dom_pts = 1 and pft_distribution_method = pft_standard)

    !
    ! !USES
    use LandunitType    , only : lun                
    use ColumnType      , only : col                
    use decompMod       , only : get_clump_bounds, get_proc_clumps
    use clm_varcon      , only : ispval
    use PatchType       , only : patch
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: nc,p,pu,pl,l,c    ! indices
    integer :: pdom(1)
    real(r8) :: sum_wtcol, sum_wtlun, sum_wtgrc

    !------------------------------------------------------------------------

    do c = bounds%begc,bounds%endc
       if (col%is_hillslope_column(c) .and. col%active(c)) then
          pdom = maxloc(patch%wtcol(col%patchi(c):col%patchf(c)))
          pdom = pdom + (col%patchi(c) - 1)

          sum_wtcol = sum(patch%wtcol(col%patchi(c):col%patchf(c)))
          sum_wtlun = sum(patch%wtlunit(col%patchi(c):col%patchf(c)))
          sum_wtgrc = sum(patch%wtgcell(col%patchi(c):col%patchf(c)))

          patch%wtcol(col%patchi(c):col%patchf(c)) = 0._r8
          patch%wtlunit(col%patchi(c):col%patchf(c)) = 0._r8
          patch%wtgcell(col%patchi(c):col%patchf(c)) = 0._r8

          patch%wtcol(pdom(1)) = sum_wtcol
          patch%wtlunit(pdom(1)) = sum_wtlun
          patch%wtgcell(pdom(1)) = sum_wtgrc
       endif
    enddo    ! end loop c

  end subroutine HillslopeDominantPft

  !------------------------------------------------------------------------
  subroutine HillslopeDominantLowlandPft(bounds)
    !
    ! !DESCRIPTION: 
    ! Reassign patch weights such that each column has a single, 
    ! dominant pft.  Use largest weight for lowland, 2nd largest
    ! weight for uplands
    ! Best performance when used with n_dom_pfts = 2

    !
    ! !USES
    use LandunitType    , only : lun                
    use ColumnType      , only : col                
    use decompMod       , only : get_clump_bounds, get_proc_clumps
    use clm_varcon      , only : ispval
    use PatchType       , only : patch
    use pftconMod       , only : pftcon, ndllf_evr_tmp_tree, nc3_nonarctic_grass, nc4_grass
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: nc,p,pu,pl,l,c    ! indices
    integer :: pdom(1),psubdom(1)
    integer :: plow, phigh
    real(r8) :: sum_wtcol, sum_wtlun, sum_wtgrc
    real(r8),allocatable :: mask(:)

    !------------------------------------------------------------------------

    do c = bounds%begc,bounds%endc
       if (col%is_hillslope_column(c) .and. col%active(c)) then
          pdom = maxloc(patch%wtcol(col%patchi(c):col%patchf(c)))
          ! create mask to exclude pdom
          allocate(mask(col%npatches(c)))
          mask(:) = 1.
          mask(pdom(1)) = 0.

          pdom = pdom + (col%patchi(c) - 1)

          psubdom = maxloc(mask*patch%wtcol(col%patchi(c):col%patchf(c)))
          psubdom = psubdom + (col%patchi(c) - 1)
          deallocate(mask)

          sum_wtcol = sum(patch%wtcol(col%patchi(c):col%patchf(c)))
          sum_wtlun = sum(patch%wtlunit(col%patchi(c):col%patchf(c)))
          sum_wtgrc = sum(patch%wtgcell(col%patchi(c):col%patchf(c)))

          patch%wtcol(col%patchi(c):col%patchf(c)) = 0._r8
          patch%wtlunit(col%patchi(c):col%patchf(c)) = 0._r8
          patch%wtgcell(col%patchi(c):col%patchf(c)) = 0._r8

          ! assumes trees are 1-8, shrubs 9-11, and grasses 12-14
          ! and puts the lowest ivt on the lowland column
          if ((patch%itype(pdom(1)) > patch%itype(psubdom(1)) &
               .and. patch%itype(psubdom(1)) > 0) .or. &
               (patch%itype(pdom(1)) == 0)) then
             plow = psubdom(1)
             phigh = pdom(1)
          else
             plow = pdom(1)
             phigh = psubdom(1)
          endif

          ! Special cases (subjective)

          ! if NET/BDT assign BDT to lowland
          if ((patch%itype(pdom(1)) == ndllf_evr_tmp_tree) .and. pftcon%is_tree(patch%itype(psubdom(1)))) then
             plow = psubdom(1)        
             phigh = pdom(1)
          endif
          ! if C3/C4 assign C4 to lowland
          if ((patch%itype(pdom(1)) == nc4_grass) .and. (patch%itype(psubdom(1)) == nc3_nonarctic_grass)) then
             plow = pdom(1)        
             phigh = psubdom(1)
          endif
          if ((patch%itype(pdom(1)) == nc3_nonarctic_grass) .and. (patch%itype(psubdom(1)) == nc4_grass)) then
             plow = psubdom(1)        
             phigh = pdom(1)
          endif

          if(col%cold(c) == ispval) then
             ! lowland column
             patch%wtcol(plow)   = sum_wtcol
             patch%wtlunit(plow) = sum_wtlun
             patch%wtgcell(plow) = sum_wtgrc
          else
             ! upland columns
             patch%wtcol(phigh)   = sum_wtcol
             patch%wtlunit(phigh) = sum_wtlun
             patch%wtgcell(phigh) = sum_wtgrc
          endif
       endif
    enddo    ! end loop c

  end subroutine HillslopeDominantLowlandPft

  !------------------------------------------------------------------------
  subroutine HillslopePftFromFile(bounds,col_pftndx)
    !
    ! !DESCRIPTION: 
    ! Reassign patch weights using indices from surface data file
    ! Assumes one patch per hillslope column
    !
    ! !USES
    use ColumnType      , only : col                
    use PatchType       , only : patch

    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(in)           :: col_pftndx(:)
    !
    ! !LOCAL VARIABLES:
    integer :: p,c    ! indices
    integer :: npatches_per_column
    logical :: check_npatches = .true.
    
    !------------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       if (col%is_hillslope_column(c) .and. col%active(c)) then
          ! In preparation for this re-weighting of patch type
          ! only first patch was given a non-zero weight in surfrd_hillslope
          npatches_per_column = 0
          do p = col%patchi(c), col%patchf(c)
             patch%itype(p) = col_pftndx(c)
             npatches_per_column = npatches_per_column + 1
          enddo
          if (check_npatches) then
             if ((npatches_per_column /= 1) .and. masterproc) then
                call endrun( 'ERROR:: number of patches per hillslope column not equal to 1'//errmsg(sourcefile, __LINE__) )
             end if
          endif
       endif
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
    character(len=*), parameter :: subname = 'HillslopeStreamOutflow'

    !-----------------------------------------------------------------------
    associate(                                                            & 
         stream_water_volume     =>    waterstatebulk_inst%stream_water_volume_lun            , & ! Input:  [real(r8) (:)   ] stream water volume (m3)
         qstreamflow             =>    waterfluxbulk_inst%qstreamflow_lun               &  ! Input:  [real(r8) (:)   ] stream water discharge (m3/s)
         )

      ! Get time step
      dtime = get_step_size_real()

      do l = bounds%begl,bounds%endl
         qstreamflow(l) = 0._r8
         if(lun%itype(l) == istsoil .and. lun%active(l)) then
            ! Streamflow calculated from Manning equation
            if(streamflow_method == streamflow_manning) then
               cross_sectional_area = stream_water_volume(l) &
                    /lun%stream_channel_length(l)
               stream_depth =  cross_sectional_area &
                    /lun%stream_channel_width(l)
               hydraulic_radius = cross_sectional_area &
                    /(lun%stream_channel_width(l) + 2*stream_depth)

               if(hydraulic_radius <= 0._r8) then
                  qstreamflow(l) = 0._r8
               else
                  flow_velocity = (hydraulic_radius)**manning_exponent &
                       * sqrt(lun%stream_channel_slope(l)) &
                       / manning_roughness
                  ! overbank flow
                  if (stream_depth > lun%stream_channel_depth(l)) then
                     if (overbank_method  == 1) then
                        ! try increasing dynamic slope
                        qstreamflow(l) = cross_sectional_area * flow_velocity &
                             *(stream_depth/lun%stream_channel_depth(l))
                     else if (overbank_method  == 2) then
                        ! try increasing flow area cross section
                        overbank_area = (stream_depth -lun%stream_channel_depth(l)) * 30._r8 * lun%stream_channel_width(l)
                        qstreamflow(l) = (cross_sectional_area + overbank_area) * flow_velocity
                     else if (overbank_method  == 3) then
                        ! try removing all overbank flow instantly
                        qstreamflow(l) = cross_sectional_area * flow_velocity &
                             + (stream_depth-lun%stream_channel_depth(l)) &
                             *lun%stream_channel_width(l)*lun%stream_channel_length(l)/dtime
                     else
                        if (masterproc) then
                           call endrun( 'ERROR:: invalid overbank_method.'//errmsg(sourcefile, __LINE__) )
                        end if
                     endif
                    
                  else
                     qstreamflow(l) = cross_sectional_area * flow_velocity
                  endif

                  ! scale streamflow by number of channel reaches
                  qstreamflow(l) = qstreamflow(l) * lun%stream_channel_number(l)

                  qstreamflow(l) = max(0._r8,min(qstreamflow(l),stream_water_volume(l)/dtime))
               endif
            else
               if (masterproc) then
                  call endrun( 'ERROR:: invalid streamflow_method'//errmsg(sourcefile, __LINE__) )
               end if
            endif
         endif ! end of istsoil
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
    type(waterdiagnosticbulk_type), intent(out) :: waterdiagnosticbulk_inst
    
    integer               :: c, l, g, i, j
    real(r8) :: qflx_surf_vol                           ! volumetric surface runoff (m3/s)
    real(r8) :: qflx_drain_perched_vol                  ! volumetric perched saturated drainage (m3/s)
    real(r8) :: qflx_drain_vol                          ! volumetric saturated drainage (m3/s)
    real(r8) :: dtime                                   ! land model time step (sec)

    character(len=*), parameter :: subname = 'HillslopeUpdateStreamWater'

    !-----------------------------------------------------------------------
    associate( & 
         stream_water_volume     =>    waterstatebulk_inst%stream_water_volume_lun, & ! Input/Output:  [real(r8) (:)   ] stream water volume (m3)
         qstreamflow             =>    waterfluxbulk_inst%qstreamflow_lun      ,    & ! Input:  [real(r8) (:)   ] stream water discharge (m3/s)
         qflx_drain              =>    waterfluxbulk_inst%qflx_drain_col,           & ! Input:  [real(r8) (:)   ]  column level sub-surface runoff (mm H2O /s)
         qflx_drain_perched      =>    waterfluxbulk_inst%qflx_drain_perched_col,   & ! Input:  [real(r8) (:)   ]  column level sub-surface runoff (mm H2O /s)
         qflx_surf               =>    waterfluxbulk_inst%qflx_surf_col        ,    & ! Input: [real(r8) (:)   ]  total surface runoff (mm H2O /s)
         stream_water_depth      =>    waterdiagnosticbulk_inst%stream_water_depth_lun   & ! Output:  [real(r8) (:)   ] stream water depth (m)
         )

       ! Get time step
       dtime = get_step_size_real()

       do l = bounds%begl,bounds%endl
          if(lun%itype(l) == istsoil) then
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
                endif
             enddo
             stream_water_volume(l) = stream_water_volume(l) &
                  - qstreamflow(l) * dtime
             
             ! account for negative drainage (via searchforwater in soilhydrology)
             if(stream_water_volume(l) < 0._r8) then
                qstreamflow(l) = qstreamflow(l) + stream_water_volume(l)/dtime
                stream_water_volume(l) = 0._r8
             endif

             stream_water_depth(l) = stream_water_volume(l) &
                  /lun%stream_channel_length(l) &
                  /lun%stream_channel_width(l)

          endif
       enddo
      
    end associate

  end subroutine HillslopeUpdateStreamWater
  
end module HillslopeHydrologyMod
